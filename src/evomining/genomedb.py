#!/usr/bin/env python3
"""
evomining.genomedb
==================

Build EvoMining's genome database from GenBank (.gbff/.gbk) input, given either as
a flat folder of GenBank files (stem = filename) or as a folder of per-genome
sub-folders each holding that genome's annotation (stem = folder name, e.g. Bakta
output). Both layouts are accepted transparently; see resolve_genome_inputs.

This replaces the former RAST-table builder. There is no more 13-column .txt, no
Corason_Rast.IDs, and no 666666.<id>.peg.<N> scheme. Genomes are parsed straight
from GenBank by the vendored corason-py parser (evomining.io), and EvoMining's
downstream steps consume three flat artifacts written here:

    GENOMES.fasta          protein FASTA, one record per CDS, header = composite ID
    genome_functions.tsv   composite_id <TAB> product
    genome_names.tsv       composite_id <TAB> organism_display_name

The composite protein ID is:

    <genome_stem>__<locus_tag>

i.e. the GenBank filename stem, a double underscore, and the gene's stable
identifier (locus_tag / protein_id / gene / synthesized, per the parser). The
stem namespaces the locus so that genomes which reuse locus tags (Bakta/Prokka
restart numbering per genome) do not collide when flattened into one FASTA. This
matches corason-py's own <genome>__<locus_tag> convention for its .gbk/cluster
names. The stem is used (not the organism name) so IDs stay stable regardless of
what else is in the batch, since organism names may be disambiguated at load time.

Optional --names TSV overrides organism display names for genomes with poor
GenBank metadata; it uses the same format as corason-py's --names/--rast-ids.

Clade / iTOL metadata is intentionally NOT produced here: clade assignment is a
study-specific overlay, not part of EvoMining. Build any clade colorstrips with a
separate tool from the genome lists.
"""

from __future__ import annotations

import sys
from pathlib import Path

from .io.loader import load_genomes, resolve_genome_inputs


ID_SEP = "__"


def composite_id(genome_stem: str, gene_id: str) -> str:
    return f"{genome_stem}{ID_SEP}{gene_id}"




def run(args):
    """Parse GenBank genomes and write EvoMining's flat genome artifacts."""
    input_dir = Path(args.input_dir)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Resolve the genome file list. Accepts either a flat directory of GenBank
    # files (stem = filename) or a directory of per-genome sub-folders holding
    # each genome's annotation (stem = folder name); --lists restricts stems.
    try:
        pairs = resolve_genome_inputs(input_dir)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}")

    wanted = _wanted_stems(args)
    if wanted is not None:
        pairs = [(p, s) for (p, s) in pairs if s in wanted]
        missing = wanted - {s for _, s in pairs}
        for stem in sorted(missing):
            print(f"  WARNING: no GenBank file for: {stem}", file=sys.stderr)

    if not pairs:
        raise SystemExit("ERROR: no genomes selected")

    paths = [p for p, _ in pairs]
    stems = {p: s for p, s in pairs}
    names_file = Path(args.names) if args.names else None

    print(f"Loading {len(paths)} GenBank genome(s)...")
    genomes = load_genomes(paths, names_file=names_file,
                           keep_pseudo=getattr(args, 'keep_pseudogenes', False),
                           stems=stems)

    fasta_path = outdir / "GENOMES.fasta"
    func_path = outdir / "genome_functions.tsv"
    names_path = outdir / "genome_names.tsv"

    total_genes = 0
    with open(fasta_path, "w") as fa, \
         open(func_path, "w") as fn, \
         open(names_path, "w") as nm:
        fn.write("protein_id\tfunction\n")
        nm.write("protein_id\tgenome_name\n")
        for g in genomes:
            gene_count = 0
            for gene in g.genes():
                cid = composite_id(g.id, gene.id)
                fa.write(f">{cid}\n")
                for i in range(0, len(gene.translation), 60):
                    fa.write(gene.translation[i:i + 60] + "\n")
                fn.write(f"{cid}\t{gene.product}\n")
                nm.write(f"{cid}\t{g.name}\n")
                gene_count += 1
            total_genes += gene_count
            print(f"  [{g.id}]  {g.name:40s}  {gene_count:>5} proteins  "
                  f"({len(g.contigs)} contigs)")

    print(f"""
{'=' * 60}
  Done!
  Genomes:           {len(genomes)}
  Total proteins:    {total_genes}

  Protein FASTA:     {fasta_path}
  Function map:      {func_path}
  Name map:          {names_path}
{'=' * 60}
""")
    return outdir


def _wanted_stems(args):
    """Set of genome stems to include from --lists / --list-dir, or None for all."""
    stems = set()
    have_lists = False

    if args.lists:
        have_lists = True
        for lst in args.lists:
            stems |= _read_list(Path(lst))
    if args.list_dir:
        have_lists = True
        for lst in sorted(Path(args.list_dir).glob("*.txt")):
            stems |= _read_list(lst)

    return stems if have_lists else None


def _read_list(path):
    out = set()
    with open(path) as fh:
        for line in fh:
            name = line.strip()
            if name:
                out.add(name)
    return out
