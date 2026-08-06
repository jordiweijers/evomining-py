#!/usr/bin/env python3
"""
evomining.analysis

Python replacement for EvoMining's heatplot.pl and downstream analysis.
Parses the BLAST results from `evomining start` and produces:

  1. Expansion heatplot (HTML interactive)
  2. Per-enzyme-family copy count table
  3. BBH (conserved metabolism) classification
  4. antiSMASH cross-reference (cosmetic: cyan tree-leaf class only)
  5. MIBiG BLAST for the protein copies
  6. Summary tables for downstream transcript cross-referencing

Composite-ID version. Genome proteins are identified by `<genome_stem>__<locus_tag>`
(written by `evomining generate-genome-db`); the RAST `666666.<id>.peg.<N>` scheme,
`Corason_Rast.IDs`, and `org_lookup` are gone. Genome identity is the part before
the first `__`; organism names come from genome_names.tsv and products from
genome_functions.tsv.

Usage:
  evomining analyze \\
      --blast-dir         run_evomining/blast \\
      --central-db        run_evomining/Cyanos.Central \\
      --genome-names      evomining_db/genome_names.tsv \\
      --genome-functions  evomining_db/genome_functions.tsv \\
      --genomes-fasta     evomining_db/GENOMES.fasta \\
      --antismash         run_evomining/cyanoSMASH.tsv \\
      --mibig             run_evomining/MiBIG_DB.faa \\
      -o evomining_results
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from .external import Tools


# Must match evomining.genomedb.ID_SEP. Kept local (like antismashdb.py) so this
# module doesn't pull in the GenBank parser just to split a string.
ID_SEP = "__"


# ---------------------------------------------------------------------------
#  Metadata loaders (composite-ID artifacts from generate-genome-db)
# ---------------------------------------------------------------------------

def load_genome_names(names_path):
    """Read genome_names.tsv (protein_id<TAB>genome_name) and collapse it to
    {genome_stem: genome_name}.

    The file lists one row per protein; every protein of a genome carries the
    same organism name, so keying on the stem (the part before the first `__`)
    is exact. This is the composite-ID replacement for `parse_rast_ids`' first
    return value; the whole `org_lookup` mechanism is no longer needed.
    """
    names = {}
    with open(names_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) < 2:
                continue
            composite = parts[0].strip()
            if not composite or composite == "protein_id":  # header row
                continue
            stem = composite.split(ID_SEP, 1)[0]
            names[stem] = parts[1].strip()
    return names


def load_functions(functions_path, wanted=None):
    """Read genome_functions.tsv (protein_id<TAB>function) -> {composite_id: product}.

    Reads genome_functions.tsv produced by generate-genome-db.
    Because BLAST outfmt 6 truncates subject IDs at whitespace and GENOMES.fasta
    headers are bare composite IDs, this table is the authoritative source of the
    product string. `wanted` (a set of composite IDs) bounds the map to proteins
    that actually appear in the analysis.
    """
    functions = {}
    with open(functions_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) < 2:
                continue
            composite = parts[0].strip()
            if not composite or composite == "protein_id":  # header row
                continue
            if wanted is not None and composite not in wanted:
                continue
            functions[composite] = parts[1]
    return functions


def parse_central_db(central_file):
    """
    Parse Central DB FASTA headers to get enzyme family info.
    Headers: >enzyme_family|family_num|copy_name|genome_id
    Returns:
      families: dict of "pathway|family_num" -> list of seed sequences
      family_names: dict of "pathway|family_num" -> readable enzyme name

    (Unchanged by the migration: the Central DB is the seed side, not the genome
    side, so it keeps its pipe-delimited header format.)
    """
    families = defaultdict(list)
    family_names = {}
    with open(central_file) as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
                parts = header.split("|")
                if len(parts) >= 3:
                    family_key = f"{parts[0]}|{parts[1]}"
                    families[family_key].append(header)
                    if family_key not in family_names:
                        name = re.sub(r"_\d+$", "", parts[2])  # strip trailing _1, _2, _3
                        family_names[family_key] = name
    return dict(families), family_names


def parse_antismash_mapping(antismash_file):
    """
    Parse the composite-ID antiSMASH NP mapping written by
    `evomining generate-antismash-db`. Format (no header):
        <genome_stem>__<locus_tag> <TAB> cf_putative <TAB> <cluster_name>

    Returns: (set of composite protein IDs in antiSMASH clusters,
              dict composite_id -> set of cluster names).

    COSMETIC-ONLY: this membership is used solely to assign the cyan
    "Secondary metabolism (antiSMASH)" tree-leaf class. It must never influence
    copy counts, expansion calls, BBH, MIBiG, or the green-clade prediction.
    Absent file -> empty set (see the invariant note in run()).
    """
    antismash_proteins = set()
    antismash_clusters = defaultdict(set)

    if not antismash_file or not Path(antismash_file).exists():
        return antismash_proteins, antismash_clusters

    with open(antismash_file) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                composite = parts[0]
                cluster = parts[2]
                antismash_proteins.add(composite)
                antismash_clusters[composite].add(cluster)

    return antismash_proteins, antismash_clusters


def parse_subject(subject):
    """A genome-protein BLAST subject/query is the composite ID `<stem>__<locus_tag>`.

    Returns (genome_stem, composite). BLAST outfmt 6 emits the full ID (composite
    IDs contain no whitespace), so composite == subject. Returns (None, None) for
    anything without an ID separator, which is skipped upstream.

    Replaces `parse_blast_header` + `org_lookup`. It intentionally does no map
    lookups: callers resolve organism name / product from genome_names.tsv /
    genome_functions.tsv only where they need them, which keeps this a pure,
    allocation-free ID split on the BLAST hot path.
    """
    composite = subject.split()[0] if subject else ""
    if ID_SEP not in composite:
        return None, None
    genome_stem = composite.split(ID_SEP, 1)[0]
    return genome_stem, composite


def parse_central_header(header_str):
    """
    Parse a Central DB BLAST subject/query header.
    Format: enzyme_family|family_num|copy_name|genome_id
    Returns: (full_family_key, pathway, family_num, copy_name)

    full_family_key = "pathway|family_num" e.g. "3PGA_AMINOACIDS|4"
    This gives individual enzyme families rather than pathway groups.
    """
    parts = header_str.split("|")
    if len(parts) >= 3:
        pathway = parts[0]
        family_num = parts[1]
        copy_name = parts[2]
        full_key = f"{pathway}|{family_num}"
        return full_key, pathway, family_num, copy_name
    return None, None, None, None


# ---------------------------------------------------------------------------
#  BLAST parsing
# ---------------------------------------------------------------------------

def parse_blast_results(blast_file):
    """
    Parse BLAST tabular output (outfmt 6).
    Returns list of dicts with query, subject, pident, evalue, bitscore, etc.
    """
    results = []
    with open(blast_file) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 12:
                results.append({
                    "query": parts[0],
                    "subject": parts[1],
                    "pident": float(parts[2]),
                    "length": int(parts[3]),
                    "mismatch": int(parts[4]),
                    "gapopen": int(parts[5]),
                    "qstart": int(parts[6]),
                    "qend": int(parts[7]),
                    "sstart": int(parts[8]),
                    "send": int(parts[9]),
                    "evalue": float(parts[10]),
                    "bitscore": float(parts[11]),
                })
    return results


def get_best_hits(blast_results, by="query"):
    """Get best hit per query (or subject), by bitscore."""
    best = {}
    for hit in blast_results:
        key = hit[by]
        if key not in best or hit["bitscore"] > best[key]["bitscore"]:
            best[key] = hit
    return best


# ---------------------------------------------------------------------------
#  Core EvoMining analysis
# ---------------------------------------------------------------------------

def compute_copy_counts(forward_blast, genome_names, method="multi"):
    """
    Count copies of each enzyme family per genome.

    forward_blast: Central DB vs Genome DB (thereand*.blast)
    Each query is a Central DB seed, each subject is a genome protein (composite ID).

    method:
      'multi' (default): original EvoMining heatplot.pl behaviour. Each protein
          is counted in ALL families it hits (deduplicated within each family).
          Correct for curated Central DBs with well-defined pathway families.
      'best': each protein assigned to its single best-matching family only.
          Better for custom single-seed enzyme DBs (e.g. transcript-based)
          where broad domain families inflate counts under multi-family counting.

    Both methods apply the original heatplot.pl filters:
      - bitscore >= 100
      - send*100/qend > 50

    Returns:
      copy_counts:    dict[family][genome_stem] = count
      family_members: dict[family][genome_stem] = list of member dicts
    """
    copy_counts = defaultdict(lambda: defaultdict(int))
    family_members = defaultdict(lambda: defaultdict(list))

    if method == "multi":
        # Original heatplot.pl: count each protein in ALL families it hits
        seen = defaultdict(set)  # (family, genome_stem) -> set of composite IDs

        for hit in forward_blast:
            family, _, _, _ = parse_central_header(hit["query"])
            if family is None:
                continue
            genome_stem, composite = parse_subject(hit["subject"])
            if genome_stem is None:
                continue
            if genome_stem not in genome_names:
                continue
            if hit["bitscore"] < 100:
                continue
            if hit["qend"] == 0 or (hit["send"] * 100 / hit["qend"]) <= 50:
                continue

            key = (family, genome_stem)
            if composite not in seen[key]:
                seen[key].add(composite)
                copy_counts[family][genome_stem] += 1
                family_members[family][genome_stem].append({
                    "function": "unknown",
                    "evalue": hit["evalue"],
                    "bitscore": hit["bitscore"],
                    "protein_id": composite,
                })

    elif method == "best":
        # Best-family-only: each protein assigned to its highest-scoring family
        protein_best_family = {}

        for hit in forward_blast:
            family, _, _, _ = parse_central_header(hit["query"])
            if family is None:
                continue
            genome_stem, composite = parse_subject(hit["subject"])
            if genome_stem is None:
                continue
            if genome_stem not in genome_names:
                continue
            if hit["bitscore"] < 100:
                continue
            if hit["qend"] == 0 or (hit["send"] * 100 / hit["qend"]) <= 50:
                continue

            if (composite not in protein_best_family
                    or hit["bitscore"] > protein_best_family[composite][1]):
                protein_best_family[composite] = (family, hit["bitscore"], genome_stem, hit["evalue"])

        for composite, (family, bitscore, genome_stem, evalue) in protein_best_family.items():
            copy_counts[family][genome_stem] += 1
            family_members[family][genome_stem].append({
                "function": "unknown",
                "evalue": evalue,
                "bitscore": bitscore,
                "protein_id": composite,
            })

    else:
        raise ValueError(f"Unknown copy count method: {method!r}. Use 'multi' or 'best'.")

    return dict(copy_counts), dict(family_members)


def compute_bbh(forward_blast, reverse_blast):
    """
    Compute Best Bidirectional Hits (BBH).
    Forward: Central DB seeds -> Genome DB
    Reverse: Genome DB -> Central DB seeds

    A BBH exists when seed S's best hit in genome G is protein P, and P's best
    hit back in the Central DB is a seed of the same family as S.

    Returns: set of composite IDs that are BBH (conserved metabolism).
    """
    # Best genome protein per (family, genome_stem)
    forward_best = defaultdict(dict)
    for hit in forward_blast:
        family, _, _, _ = parse_central_header(hit["query"])
        if family is None:
            continue
        genome_stem, composite = parse_subject(hit["subject"])
        if genome_stem is None:
            continue

        key = (family, genome_stem)
        if key not in forward_best or hit["bitscore"] > forward_best[key]["bitscore"]:
            forward_best[key] = {**hit, "composite": composite, "genome_stem": genome_stem, "family": family}

    # Best Central family per genome protein
    reverse_best = {}
    for hit in reverse_blast:
        genome_stem, composite = parse_subject(hit["query"])
        if genome_stem is None:
            continue

        family, _, _, _ = parse_central_header(hit["subject"])
        if family is None:
            continue

        if composite not in reverse_best or hit["bitscore"] > reverse_best[composite]["bitscore"]:
            reverse_best[composite] = {"family": family, "bitscore": hit["bitscore"]}

    # Find BBH
    bbh_proteins = set()
    for (family, genome_stem), fwd in forward_best.items():
        composite = fwd["composite"]
        rev = reverse_best.get(composite)
        if rev is not None and rev["family"] == family:
            bbh_proteins.add(composite)

    return bbh_proteins


def compute_expansions(copy_counts, genome_names, sd=1.0):
    """
    Identify significantly expanded enzyme families per genome.
    Expansion = copy count > mean + `sd` * std across all genomes.

    `sd` defaults to 1.0 (evominingAnalysis_1SD.py behaviour). Pass sd=2.0 for the
    mean + 2SD variant described in the paper.
    """
    expansions = {}
    all_genomes = sorted(genome_names.keys())

    for family, counts in sorted(copy_counts.items()):
        values = [counts.get(stem, 0) for stem in all_genomes]
        mean_val = np.mean(values)          # correct mean (divide by n)
        std_val = np.std(values)            # population SD
        threshold = mean_val + sd * std_val
        expanded = [
            (stem, counts.get(stem, 0))
            for stem in all_genomes
            if counts.get(stem, 0) > threshold  # strict > not >=
        ]
        expansions[family] = {
            "mean": mean_val,
            "std": std_val,
            "threshold": threshold,
            "max": max(values),
            "min": min(values),
            "expanded_genomes": expanded,
            "n_expanded": len(expanded),
        }

    return expansions


def classify_proteins(family_members, bbh_set, antismash_set, expansions, genome_names):
    """
    Emits per copy:
      - is_bbh        : conserved central metabolism (Best Bidirectional Hit)
      - in_antismash  : falls in an antiSMASH cluster -- COSMETIC OVERLAY ONLY
      - classification: the tree-leaf colour label, retained for trees:
                          conserved (red) / antismash (cyan) /
                          transition (purple) / normal (grey)

    The `classification` label folds in antiSMASH (that is what colours cyan
    leaves), so it is the ONE antiSMASH-dependent output. Everything the pipeline
    predicts keys on `is_bbh` + `is_expanded_genome` instead, both of which are
    antiSMASH-independent -- so the prediction set is identical with or without
    --antismash. antiSMASH membership is tested against the composite-ID set.
    """
    classified = []

    for family, genome_members in family_members.items():
        exp_info = expansions.get(family, {})
        threshold = exp_info.get("threshold", float("inf"))

        for genome_stem, members in genome_members.items():
            n_copies = len(members)
            is_expanded_genome = n_copies > threshold

            for member in members:
                protein_id = member["protein_id"]           # composite id
                is_bbh = protein_id in bbh_set
                is_antismash = protein_id in antismash_set

                if is_bbh and is_antismash:
                    classification = "transition"
                elif is_bbh:
                    classification = "conserved"
                elif is_antismash:
                    classification = "antismash"
                else:
                    classification = "normal"

                classified.append({
                    "family": family,
                    "genome_stem": genome_stem,
                    "genome_name": genome_names.get(genome_stem, "unknown"),
                    "protein_id": protein_id,
                    "function": member["function"],
                    "evalue": member["evalue"],
                    "bitscore": member["bitscore"],
                    "is_bbh": is_bbh,               # conserved central metabolism (antiSMASH-independent)
                    "in_antismash": is_antismash,   # cosmetic overlay only; never gates a prediction
                    "classification": classification,
                    "n_copies_in_genome": n_copies,
                    "is_expanded_genome": is_expanded_genome,
                })

    return classified


# ---------------------------------------------------------------------------
#  MIBiG BLAST for the protein copies
# ---------------------------------------------------------------------------

def blast_vs_mibig(query_proteins, genomes_fasta, mibig_fasta, output_dir, threads=8, tools=None):
    """
    BLAST protein copies against MIBiG.
    Returns dict: composite_id -> best MIBiG hit info.
    """
    if not mibig_fasta or not Path(mibig_fasta).exists():
        print("  No MIBiG database provided, skipping NP BLAST")
        return {}

    tools = tools or Tools()
    tools.require("makeblastdb", "blastp")

    query_fasta = output_dir / "expanded_queries.faa"
    wanted_ids = {p["protein_id"] for p in query_proteins}  # composite IDs

    # Pull query sequences straight from GENOMES.fasta, whose headers are the
    # composite protein IDs (`>{stem}__{locus_tag}`) from generate-genome-db.
    # No per-genome .faa re-reading and no org_lookup any more.
    genomes_fasta = Path(genomes_fasta)
    if not genomes_fasta.exists():
        raise SystemExit(f"ERROR: GENOMES.fasta not found: {genomes_fasta}")

    sequences = {}
    current_id = None
    current_seq = []
    with open(genomes_fasta) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id and current_id in wanted_ids and current_seq:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]  # composite id
                current_seq = []
            else:
                current_seq.append(line)
        if current_id and current_id in wanted_ids and current_seq:
            sequences[current_id] = "".join(current_seq)

    if not sequences:
        print("  No protein sequences found for MIBiG BLAST")
        return {}

    with open(query_fasta, "w") as fh:
        for composite, seq in sequences.items():
            fh.write(f">{composite}\n{seq}\n")

    print(f"  BLASTing {len(sequences)} proteins vs MIBiG...")

    mibig_db = output_dir / "mibig_blastdb"
    tools.run("makeblastdb",
              ["-in", mibig_fasta, "-dbtype", "prot", "-out", mibig_db],
              check=True)

    blast_out = output_dir / "expanded_vs_mibig.blast"
    tools.run("blastp",
              ["-query", query_fasta, "-db", mibig_db, "-out", blast_out,
               "-outfmt", "6", "-evalue", "0.001",
               "-max_target_seqs", "500", "-num_threads", threads],
              check=True)

    mibig_hits = {}
    if blast_out.exists():
        results = parse_blast_results(str(blast_out))
        best = get_best_hits(results, by="query")  # best-hit-only per query protein
        for composite, hit in best.items():
            mibig_hits[composite] = {
                "mibig_subject": hit["subject"],
                "mibig_pident": hit["pident"],
                "mibig_evalue": hit["evalue"],
                "mibig_bitscore": hit["bitscore"],
            }

    print(f"  Found {len(mibig_hits)} proteins with MIBiG hits")
    return mibig_hits


# ---------------------------------------------------------------------------
#  Output: tables
# ---------------------------------------------------------------------------

def write_copy_count_table(copy_counts, genome_names, expansions, output_file):
    """Write per-genome copy count matrix as TSV."""
    families = sorted(copy_counts.keys())
    all_genomes = sorted(genome_names.keys())

    with open(output_file, "w") as fh:
        fh.write("genome_stem\tgenome_name\t" + "\t".join(families) + "\n")

        for stem in all_genomes:
            genome = genome_names.get(stem, "unknown")
            counts = [str(copy_counts.get(f, {}).get(stem, 0)) for f in families]
            fh.write(f"{stem}\t{genome}\t" + "\t".join(counts) + "\n")


def write_expansion_summary(expansions, genome_names, output_file):
    """Write expansion summary TSV."""
    with open(output_file, "w") as fh:
        fh.write("enzyme_family\tmean_copies\tstd\tthreshold\tmax_copies\t"
                 "n_expanded_genomes\texpanded_genomes\n")

        for family in sorted(expansions.keys()):
            info = expansions[family]
            expanded_str = "; ".join([
                f"{genome_names.get(stem, stem)}({cnt})"
                for stem, cnt in info["expanded_genomes"]
            ])
            fh.write(
                f"{family}\t{info['mean']:.2f}\t{info['std']:.2f}\t"
                f"{info['threshold']:.2f}\t{info['max']}\t"
                f"{info['n_expanded']}\t{expanded_str}\n"
            )


def write_classified_proteins(classified, mibig_hits, output_file):
    """Write full classified protein table."""
    with open(output_file, "w") as fh:
        fh.write(
            "family\tgenome_stem\tgenome_name\tprotein_id\tfunction\t"
            "evalue\tbitscore\tis_bbh\tin_antismash\tclassification\t"
            "n_copies\tis_expanded\t"
            "mibig_hit\tmibig_pident\tmibig_evalue\n"
        )

        for p in sorted(classified, key=lambda x: (x["family"], x["genome_stem"])):
            mibig = mibig_hits.get(p["protein_id"], {})
            fh.write(
                f"{p['family']}\t{p['genome_stem']}\t{p['genome_name']}\t"
                f"{p['protein_id']}\t{p['function']}\t"
                f"{p['evalue']}\t{p['bitscore']}\t"
                f"{p['is_bbh']}\t{p['in_antismash']}\t{p['classification']}\t"
                f"{p['n_copies_in_genome']}\t{p['is_expanded_genome']}\t"
                f"{mibig.get('mibig_subject', 'NA')}\t"
                f"{mibig.get('mibig_pident', 'NA')}\t"
                f"{mibig.get('mibig_evalue', 'NA')}\n"
            )


# ---------------------------------------------------------------------------
#  Output: interactive heatplot HTML
# ---------------------------------------------------------------------------

def write_heatplot_html(copy_counts, genome_names, expansions, output_file, family_names=None):
    """Generate an interactive heatplot as a standalone HTML file."""
    families = sorted(copy_counts.keys())
    all_genomes = sorted(genome_names.keys())
    if family_names is None:
        family_names = {}

    # Build data matrix
    data_rows = []
    for stem in all_genomes:
        genome = genome_names.get(stem, "unknown")
        counts = [copy_counts.get(f, {}).get(stem, 0) for f in families]
        data_rows.append({
            "genome_stem": stem,
            "genome": genome,
            "counts": counts,
        })

    # Compute thresholds for coloring
    thresholds = {}
    for f in families:
        info = expansions.get(f, {})
        thresholds[f] = info.get("threshold", 0)

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>EvoMining Expansion Heatplot</title>
<style>
  * {{ margin: 0; padding: 0; box-sizing: border-box; }}
  body {{ font-family: 'Courier New', monospace; background: #1a1a2e; color: #e0e0e0; padding: 20px; }}
  h1 {{ color: #e94560; margin-bottom: 5px; font-size: 1.4em; }}
  .subtitle {{ color: #888; margin-bottom: 20px; font-size: 0.85em; }}
  .controls {{ margin-bottom: 15px; display: flex; gap: 15px; align-items: center; flex-wrap: wrap; }}
  .controls label {{ font-size: 0.8em; color: #aaa; }}
  .controls select, .controls input {{ background: #16213e; color: #e0e0e0; border: 1px solid #333;
    padding: 4px 8px; font-family: inherit; font-size: 0.8em; }}
  .container {{ overflow-x: auto; overflow-y: auto; max-height: 85vh; }}
  table {{ border-collapse: collapse; font-size: 0.7em; }}
  th {{ position: sticky; top: 0; background: #16213e; color: #e94560; padding: 4px 6px;
    border: 1px solid #333; z-index: 2; white-space: nowrap; }}
  th.genome-col {{ position: sticky; left: 0; z-index: 3; background: #16213e; }}
  td {{ padding: 3px 6px; border: 1px solid #222; text-align: center; white-space: nowrap; }}
  td.genome-cell {{ position: sticky; left: 0; background: #16213e; text-align: left;
    color: #ccc; z-index: 1; max-width: 250px; overflow: hidden; text-overflow: ellipsis; }}
  .cell-0 {{ background: #1a1a2e; color: #444; }}
  .cell-normal {{ background: #0f3460; color: #7ec8e3; }}
  .cell-expanded {{ background: #e94560; color: #fff; font-weight: bold; }}
  .legend {{ margin-top: 15px; display: flex; gap: 20px; font-size: 0.8em; }}
  .legend-item {{ display: flex; align-items: center; gap: 5px; }}
  .legend-box {{ width: 16px; height: 16px; border: 1px solid #444; }}
  .stats {{ margin-top: 10px; font-size: 0.8em; color: #888; }}
  tr:hover td {{ opacity: 0.85; }}
  .tooltip {{ position: fixed; background: #16213e; border: 1px solid #e94560;
    padding: 8px 12px; font-size: 0.8em; pointer-events: none; z-index: 100;
    display: none; max-width: 350px; }}
</style>
</head>
<body>
<h1>EvoMining Expansion Heatplot</h1>
<div class="subtitle">{len(all_genomes)} genomes &times; {len(families)} enzyme families</div>

<div class="controls">
  <label>Show only expanded:
    <input type="checkbox" id="expandedOnly" onchange="filterTable()">
  </label>
  <label>Search genome:
    <input type="text" id="genomeSearch" placeholder="type to filter..." oninput="filterTable()">
  </label>
</div>

<div class="container">
<table id="heatTable">
<thead>
<tr>
  <th class="genome-col">Genome</th>
"""

    for f in families:
        name = family_names.get(f, f)
        short = name[:25]
        html += f'  <th title="{f}: {name}">{short}</th>\n'

    html += "</tr>\n</thead>\n<tbody>\n"

    for row in data_rows:
        has_expansion = any(
            row["counts"][i] > thresholds.get(families[i], float("inf"))
            for i in range(len(families))
        )
        html += f'<tr data-expanded="{1 if has_expansion else 0}" '
        html += f'data-genome="{row["genome"].lower()}">\n'
        html += f'  <td class="genome-cell" title="{row["genome"]}">{row["genome"]}</td>\n'

        for i, cnt in enumerate(row["counts"]):
            fam = families[i]
            thresh = thresholds.get(fam, float("inf"))
            if cnt == 0:
                cls = "cell-0"
            elif cnt > thresh:
                cls = "cell-expanded"
            else:
                cls = "cell-normal"
            html += f'  <td class="{cls}" title="{fam}: {cnt} copies (threshold: {thresh:.1f})">{cnt}</td>\n'

        html += "</tr>\n"

    html += """</tbody>
</table>
</div>

<div class="legend">
  <div class="legend-item"><div class="legend-box cell-0"></div> 0 copies</div>
  <div class="legend-item"><div class="legend-box cell-normal"></div> Normal copies</div>
  <div class="legend-item"><div class="legend-box cell-expanded"></div> Expanded (&gt; mean + n&sigma;)</div>
</div>

<div class="stats" id="stats"></div>

<script>
function filterTable() {
  const expandedOnly = document.getElementById('expandedOnly').checked;
  const search = document.getElementById('genomeSearch').value.toLowerCase();
  const rows = document.querySelectorAll('#heatTable tbody tr');
  let shown = 0;
  rows.forEach(r => {
    let show = true;
    if (expandedOnly && r.dataset.expanded === '0') show = false;
    if (search && !r.dataset.genome.includes(search)) show = false;
    r.style.display = show ? '' : 'none';
    if (show) shown++;
  });
  document.getElementById('stats').textContent = `Showing ${shown} of """ + str(len(all_genomes)) + """ genomes`;
}
filterTable();
</script>
</body>
</html>"""

    with open(output_file, "w") as fh:
        fh.write(html)


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------

def run(args, tools=None):
    """EvoMining analysis: copy counts, expansions, classification, MIBiG BLAST."""
    tools = tools or Tools()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    blast_dir = Path(args.blast_dir)

    # ---- Load metadata ----
    print("Loading metadata...")
    genome_names = load_genome_names(args.genome_names)   # {genome_stem: genome_name}
    print(f"  {len(genome_names)} genomes in database")

    families, family_names = parse_central_db(args.central_db)
    print(f"  {len(families)} enzyme families in Central DB")

    antismash_proteins, antismash_clusters = parse_antismash_mapping(args.antismash)
    print(f"  {len(antismash_proteins)} proteins in antiSMASH clusters (cosmetic: cyan leaves only)")

    # ---- Parse BLAST results ----
    print("\nParsing BLAST results...")
    blast_files = list(blast_dir.glob("thereand*.blast"))
    if not blast_files:
        raise SystemExit(f"ERROR: no thereand*.blast found in {blast_dir}")
    genomes_name = blast_files[0].stem[len("thereand"):]  # strip "thereand" prefix
    forward_file = blast_dir / f"thereand{genomes_name}.blast"
    reverse_file = blast_dir / f"backagain{genomes_name}.blast"

    if not forward_file.exists() or not reverse_file.exists():
        raise SystemExit(f"ERROR: BLAST files not found in {blast_dir}")

    forward_blast = parse_blast_results(str(forward_file))
    print(f"  Forward BLAST (Central vs Genomes): {len(forward_blast)} hits")

    reverse_blast = parse_blast_results(str(reverse_file))
    print(f"  Reverse BLAST (Genomes vs Central): {len(reverse_blast)} hits")

    # ---- Compute copy counts ----
    print("\nComputing copy counts per enzyme family per genome...")
    copy_counts, family_members = compute_copy_counts(forward_blast, genome_names, method=args.ccm)
    print(f"  {len(copy_counts)} families with hits")

    # ---- Compute BBH ----
    print("\nComputing Best Bidirectional Hits (conserved metabolism)...")
    bbh_set = compute_bbh(forward_blast, reverse_blast)
    print(f"  {len(bbh_set)} BBH proteins (conserved metabolism)")

    # ---- Compute expansions ----
    print(f"\nIdentifying expanded enzyme families (threshold = mean + {args.sd}*SD)...")
    expansions = compute_expansions(copy_counts, genome_names, sd=args.sd)
    n_exp = sum(1 for e in expansions.values() if e["n_expanded"])
    print(f"  {n_exp} families expanded in at least one genome")

    # ---- Classify proteins ----
    print("\nClassifying all protein copies...")
    classified = classify_proteins(family_members, bbh_set, antismash_proteins, expansions, genome_names)

    # ---- Backfill protein functions from genome_functions.tsv ----
    # BLAST outfmt 6 truncates subject IDs at whitespace and GENOMES.fasta headers
    # are bare composite IDs, so genome_functions.tsv is the only source of the
    # product string. Without it every downstream label (tree leaves included)
    # reads "unknown".
    print("\nLoading protein functions from genome_functions.tsv...")
    wanted = {p["protein_id"] for p in classified}
    functions = load_functions(args.genome_functions, wanted=wanted)
    n_filled = 0
    for prot in classified:
        func = functions.get(prot["protein_id"])
        if func:
            prot["function"] = func
            n_filled += 1
    print(f"  {n_filled}/{len(classified)} protein copies annotated")
    if classified and n_filled == 0:
        print("  WARNING: no functions recovered — is --genome-functions pointing at the "
              "genome_functions.tsv written by `evomining generate-genome-db`?", file=sys.stderr)

    class_counts = defaultdict(int)
    for prot in classified:
        class_counts[prot["classification"]] += 1
    for cls, cnt in sorted(class_counts.items()):
        print(f"  {cls}: {cnt}")

    # ---- BLAST all copies vs MIBiG ----
    # Mirrors original heatplot.pl, which BLASTs ALL proteins (tabla2 = all hits).
    query_proteins = list(classified)
    print(f"\nProteins for MIBiG BLAST: {len(query_proteins)}")

    mibig_hits = {}
    if args.mibig:
        mibig_hits = blast_vs_mibig(
            query_proteins, args.genomes_fasta, args.mibig, outdir,
            threads=args.threads, tools=tools,
        )
    # MIBiG hits select reference sequences for trees only; genome protein
    # classifications are NOT updated — matches the original pipeline.

    # ---- Write outputs ----
    print("\nWriting output files...")

    cc_file = outdir / "copy_count_matrix.tsv"
    write_copy_count_table(copy_counts, genome_names, expansions, cc_file)
    print(f"  {cc_file}")

    exp_file = outdir / "expansion_summary.tsv"
    write_expansion_summary(expansions, genome_names, exp_file)
    print(f"  {exp_file}")

    class_file = outdir / "classified_proteins.tsv"
    write_classified_proteins(classified, mibig_hits, class_file)
    print(f"  {class_file}")

    heat_file = outdir / "heatplot.html"
    write_heatplot_html(copy_counts, genome_names, expansions, heat_file, family_names)
    print(f"  {heat_file}")

    # No predictions file is written here: the EvoMining recruitment call is a
    # phylogenetic-tree judgment (green clades near MIBiG leaves, made in trees),
    # not something the copy-count/BLAST stage can determine. classified_proteins.tsv
    # already carries the per-copy candidate signals (is_bbh, in_antismash,
    # is_expanded, MIBiG hit) for anyone who wants to filter the shortlist.

    print(f"\n{'=' * 60}")
    print(f"  EvoMining Analysis Complete")
    print(f"  Genomes:              {len(genome_names)}")
    print(f"  Enzyme families:      {len(families)}")
    print(f"  Total protein copies: {len(classified)}")
    print(f"  Conserved (BBH):      {class_counts.get('conserved', 0)}")
    print(f"  In antiSMASH cluster: {class_counts.get('antismash', 0) + class_counts.get('transition', 0)}  (cosmetic overlay)")
    print(f"  Expansion threshold:  mean + {args.sd}*SD")
    print(f"{'=' * 60}")
    return outdir
