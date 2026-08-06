#!/usr/bin/env python3
"""
evomining.antismashdb
=====================

Build the antiSMASH NP mapping: which genome proteins fall inside antiSMASH-
predicted BGC regions. This tells `analyze` which enzyme-family copies sit in a
secondary-metabolite cluster (the cyan / "Secondary metabolism (antiSMASH)"
class in the trees).

Composite-ID version. antiSMASH region GBKs identify member genes by locus_tag,
and EvoMining's protein IDs are now <genome_stem>__<locus_tag>, so the mapping is
a direct prefixing -- no peg-number reconstruction, no Bakta .faa re-reading, no
RAST ids. (The old builder existed almost entirely to translate locus_tag into a
positional peg number, which composite IDs make unnecessary.)

To guarantee the emitted IDs match what `generate-genome-db` produced, the set of
valid genome stems is read from genome_names.tsv. An antiSMASH directory whose
name is not among those stems is warned about and skipped, rather than silently
emitting IDs that match nothing in the copy-count matrix.

Output (tab-separated):
    <genome_stem>__<locus_tag>    cf_putative    <cluster_name>
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


def load_genome_stems(genome_names_path):
    """Genome stems and full protein IDs the genome step produced.

    genome_names.tsv is: protein_id <TAB> genome_name, where protein_id is
    <stem>__<locus_tag>. Returns (stems, protein_ids). Stems drive antiSMASH
    directory matching; protein_ids filter emitted mappings so cyanoSMASH only
    contains locus_tags that are real EvoMining proteins (a region GBK also
    lists pseudogenes / RNA genes that were dropped from the protein DB).
    """
    stems = set()
    proteins = set()
    with open(genome_names_path) as fh:
        header = fh.readline()  # protein_id\tgenome_name
        for line in fh:
            pid = line.split("\t", 1)[0].strip()
            if not pid:
                continue
            proteins.add(pid)
            stem = pid.split("__", 1)[0]
            if stem:
                stems.add(stem)
    return stems, proteins


def find_antismash_dir(antismash_base, genome_stem):
    """Find the antiSMASH output directory for a genome.

    Searches the top level and one level of clade subdirectories for a directory
    named exactly `genome_stem`.
    """
    base = Path(antismash_base)
    if (base / genome_stem).is_dir():
        return base / genome_stem
    for subdir in base.iterdir():
        if not subdir.is_dir():
            continue
        candidate = subdir / genome_stem
        if candidate.is_dir():
            return candidate
    return None


def extract_region_locus_tags(region_gbk_path):
    """All unique locus_tags in a region GBK's CDS features."""
    locus_tags = set()
    with open(region_gbk_path) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.startswith('/locus_tag="'):
                locus_tags.add(stripped.split('"')[1])
    return locus_tags


def extract_region_name(gbk_filename):
    """A cluster name from the region GBK filename, e.g. contig_1.region001."""
    name = gbk_filename.replace(".gbk", "")
    m = re.search(r"(contig_\S+\.region\d+)", name)
    if m:
        return m.group(1)
    m = re.search(r"((?:scaffold|chromosome|chr)\S*\.region\d+)", name, re.IGNORECASE)
    if m:
        return m.group(1)
    return name


def run(args):
    """Build the antiSMASH NP mapping keyed on composite IDs."""
    stems, proteins = load_genome_stems(args.genome_names)
    print(f"Loaded {len(stems)} genome stems / {len(proteins)} proteins "
          f"from {Path(args.genome_names).name}\n")

    total_mappings = 0
    total_regions = 0
    genomes_with_bgcs = 0
    genomes_missing_antismash = 0
    skipped_non_protein = 0

    with open(args.output, "w") as out_fh:
        for stem in sorted(stems):
            as_dir = find_antismash_dir(args.antismash_dir, stem)
            if as_dir is None:
                genomes_missing_antismash += 1
                continue

            region_gbks = [g for g in sorted(as_dir.glob("*region*.gbk"))
                           if ".region" in g.name]
            if not region_gbks:
                continue

            genomes_with_bgcs += 1
            genome_mappings = 0

            for region_gbk in region_gbks:
                cluster_name = extract_region_name(region_gbk.name)
                total_regions += 1
                for tag in sorted(extract_region_locus_tags(region_gbk)):
                    composite = f"{stem}__{tag}"
                    # A region GBK also lists locus_tags for pseudogenes / RNA
                    # genes that were dropped from the protein DB. Only emit
                    # mappings for real EvoMining proteins, so every cyanoSMASH
                    # row corresponds to something in the copy-count matrix.
                    if composite not in proteins:
                        skipped_non_protein += 1
                        continue
                    out_fh.write(f"{composite}\tcf_putative\t{cluster_name}\n")
                    genome_mappings += 1

            total_mappings += genome_mappings
            if genome_mappings > 0:
                print(f"  {stem:50s}  {len(region_gbks):>3} regions  "
                      f"{genome_mappings:>5} proteins mapped")

    # Warn about antiSMASH directories that match no known genome stem -- these
    # would have produced IDs absent from the copy-count matrix.
    unmatched = _unmatched_antismash_dirs(args.antismash_dir, stems)
    if unmatched:
        print(f"\n  WARNING: {len(unmatched)} antiSMASH director(ies) match no genome "
              f"stem in {Path(args.genome_names).name} and were ignored:", file=sys.stderr)
        for name in sorted(unmatched)[:10]:
            print(f"    {name}", file=sys.stderr)
        if len(unmatched) > 10:
            print(f"    ... and {len(unmatched) - 10} more", file=sys.stderr)
        print("  (ensure antiSMASH directory names equal your cleaned .gbff stems)",
              file=sys.stderr)

    print(f"\n{'=' * 60}")
    print(f"  Genomes with BGCs:       {genomes_with_bgcs}")
    print(f"  Total BGC regions:       {total_regions}")
    print(f"  Total protein mappings:  {total_mappings}")
    print(f"  Non-protein locus_tags skipped (pseudo/RNA): {skipped_non_protein}")
    print(f"  Genomes without antiSMASH dir: {genomes_missing_antismash}")
    print(f"  Output:                  {args.output}")
    print(f"{'=' * 60}")
    return Path(args.output)


def _unmatched_antismash_dirs(antismash_base, stems):
    """antiSMASH leaf directories whose name is not a known genome stem."""
    base = Path(antismash_base)
    seen = set()
    unmatched = set()

    def consider(d):
        if d.name in stems:
            seen.add(d.name)
        else:
            # only flag dirs that actually contain region GBKs
            if any(d.glob("*region*.gbk")):
                unmatched.add(d.name)

    for entry in base.iterdir():
        if not entry.is_dir():
            continue
        # a clade subdir contains per-genome dirs; a genome dir contains region gbks
        if any(entry.glob("*region*.gbk")):
            consider(entry)
        else:
            for sub in entry.iterdir():
                if sub.is_dir():
                    consider(sub)
    return unmatched
