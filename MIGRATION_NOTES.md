# RAST -> GenBank migration status

## Done (tested, on the new composite-ID scheme)

- **generate-genome-db** — GenBank (.gbff) folder in; writes GENOMES.fasta,
  genome_functions.tsv, genome_names.tsv. Protein ID = `<genome_stem>__<locus_tag>`.
  Parser vendored from corason-py into evomining/io/. `--keep-pseudogenes` flag.
  Validated on the real 344-genome set.
- **start** — consumes GENOMES.fasta directly. Output: thereand*.blast (forward,
  Central vs genomes) and backagain*.blast (reverse). Subjects are composite IDs.
- **generate-antismash-db** — composite-ID mapping. Reads genome_names.tsv for the
  authoritative stem set; antiSMASH dirs not matching a stem are warned + skipped
  (option b). Output rows: `<stem>__<locus_tag>  cf_putative  <cluster_name>`.

## TODO: analyze retarget (the big one)

- link to CORASON


## Deferred for now

- GATOR-GC handoff (evominingtoGatorGC_perEFDB.py), classify_gator_windows.py,
  This selects all genomes fromt he evomining tree, then selects queries for GATOR-GC
  using nw-utilities. To do this it scans the tree branches for green leaves, then one branch backwards
  until it finds a non-green leaf. That collection should correspond to a diverging clade of expanded enzymes.
  THen, it selects a representative sequences based on the highest bitscore in the BBH blast. Not sure if this
  is biologically relevant (is this the "least" divergent sequence from the clade?). But this works.
  It then ran gator-gc using those genomes and that clade. It works, but if there are multiple expansion clades
  you sometimes get the same gator-gc output for the different clades.
