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

The ID layer must switch from RAST (`fig|666666.rast.peg.N` + org_lookup) to the
composite scheme. Concrete changes in analysis.py:

- DELETE `parse_rast_ids` / `org_lookup`. Genome identity = the part before `__`.
- Replace with reading genome_names.tsv (composite_id -> genome_name) and
  genome_functions.tsv (composite_id -> product).
- `parse_blast_header(subject)` -> `genome_stem, composite = subject.split("__", 1)`;
  return (genome_stem, composite, genome_name_from_map, function_from_map).
- `parse_antismash_mapping` -> new format: col0 = composite id, col1 = cf_putative,
  col2 = cluster_name. No org_lookup. Returns set of composite IDs + clusters.
  ANTISMASH IS COSMETIC-ONLY: it assigns the cyan tree-leaf class and NOTHING
  else. It must never affect copy counts, expansion calls, BBH, MIBiG, or the
  green-clade prediction. Empty set when --antismash is absent; pipeline output
  (predictions) identical with or without it.
- `compute_copy_counts` -> key on genome_stem (was rast_id); protein key = composite
  (was `{rast_id}.{peg}`); the `fig_id` field emitted into family_members becomes the
  composite id (was `666666.{rast_id}.{peg}`). Keep filters bitscore>=100, cov>50.
- `compute_bbh` -> same rast_id->genome_stem swap; BBH set holds composite IDs.
- `classify_proteins` -> keys on composite IDs; antiSMASH membership tested against
  the composite-ID set. Logic otherwise unchanged.
- `load_functions` -> from genome_functions.tsv (already built as io-side concept).
- MIBiG blast query FASTA -> composite-ID headers; downstream parsing follows.
- CLI: drop `--rast-ids` from analyze; add `--genome-names` / `--genome-functions`
  (conventional defaults evomining_db/genome_names.tsv, genome_functions.tsv).

## TODO after analyze: trees retarget

- Leaf label format `<stem>__<locus_tag>|<genome_name>|<function>` (was
  `666666.rast.peg|genome|function`). Read genome_names.tsv + genome_functions.tsv.
- Drop rast_to_orig indirection. classified_proteins.tsv already carries composite
  fig_id + function after the analyze retarget, so trees reads those columns.

## Deferred (explicitly out of scope for now)

- GATOR-GC handoff (evominingtoGatorGC_perEFDB.py), classify_gator_windows.py,
  patch scripts. All assume fig| IDs; migrate when GATOR-GC work resumes.
