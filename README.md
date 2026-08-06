# EvoMining

> **Migration in progress (GenBank input).** `generate-genome-db` now reads a
> folder of GenBank (`.gbff`/`.gbk`) files directly and writes composite-ID
> artifacts (`GENOMES.fasta`, `genome_functions.tsv`, `genome_names.tsv`) — no
> RAST. The GenBank parser is vendored/adapted from corason-py. `start` now consumes GENOMES.fasta directly. The
> `analyze` / `trees` steps are still being retargeted from the old
> RAST/`fig|` scheme to the new composite IDs; run `generate-genome-db` and
> inspect its artifacts, but do not rely on the full chain until this note is
> removed.

> **Pseudogenes.** `generate-genome-db` drops pseudogene CDS by default (they are
> not functional enzymes and would inflate enzyme-family copy counts). Pass
> `--keep-pseudogenes` to retain them, e.g. for continuity with older RAST-based
> runs, which kept them.



Enzyme family expansion analysis for biosynthetic gene cluster discovery — a
Python reimplementation of EvoMining (Selem-Mojica et al. 2019), packaged as a
single conda environment and a single `evomining` command.

GATOR-GC, CORASON and antiSMASH stay in their own environments. They pin
conflicting Perl/Python stacks, and they are separate steps in the workflow
rather than part of EvoMining proper.

## Install

```bash
micromamba env create -f environment.yml
micromamba activate evomining
pip install -e .          # editable; use `pip install .` for a fixed install
evomining check           # verify every external tool resolves
```

The environment provides `blast`, `muscle` (v5), `mafft`, `trimal`, `fasttree`
and `newick_utils`. Nothing is hardcoded to `/vol/local/conda_envs/...` any
more — tools are taken from `$PATH`.

If a tool has to come from a different environment (an old BLAST to reproduce a
previous run, say), point at it explicitly:

```bash
evomining check --tool-env blastp=/vol/local/conda_envs/blast+
```

which makes that one call run as `micromamba run -p <prefix> blastp ...`.

## The working directory is the workspace

There is no project file and no `init` step. The directory you run commands from
is the workspace, and every step reads and writes at fixed, predictable
locations beneath it:

```
evomining_db/GENOMES/          genome DB (13-col .txt + bare .faa)
evomining_db/Corason_Rast.IDs
enzymes_db.faa                 Central DB
antismash_db.tsv               antiSMASH NP mapping (optional)
evomining_blast/<run>/blast/   forward + reverse BLAST
evomining_analysis/            copy counts, expansions, classified_proteins.tsv
evomining_trees/               per-family trees
```

So `analyze` looks for `GENOMES/` and the BLAST output where `generate-genome-db`
and `start` put them, checks they are present, and only needs a flag when you
deviate from the layout. If a prerequisite is missing, the error names the file
and the command that produces it. Run the steps in order from one folder and no
paths need repeating.

## Workflow

```bash
# 1. genome DB (see "Genome input layout" below for the directory structure) — Bakta -> 13-column RAST format (EvoMining + CORASON)
evomining generate-genome-db \
    -i /vol/databases/genomes_cyanobacteria_Nostocales/genomes_annotated \
    --list-dir ../0.clades/symb_clades \
    --clade-map                      # optional: metadata for iTOL colorstrips

# 2. Central (enzyme) DB — from cleaned transcript headers
evomining generate-enzyme-db --fasta core_symbiont_proteins_clean.faa

# 3. antiSMASH NP mapping (optional)
evomining generate-antismash-db \
    --bakta-dir /vol/databases/.../genomes_annotated \
    --antismash-dir ../6.antismash

# 4. run — supply MIBiG to analyze (and to trees, which reads analyze's blast)
evomining start   --np-db MiBIG_DB.faa   # np name only labels the output dir
evomining analyze --mibig /path/to/MiBIG_DB.faa
evomining trees   --mibig /path/to/MiBIG_DB.faa
```

Every step after `generate-genome-db` can be run with no path arguments — the
only thing you pass repeatedly is `--mibig`, since MIBiG lives wherever you keep
it rather than in the workspace.

## Notes on behaviour

### Genome input layout

`generate-genome-db` expects one subdirectory per genome under `--input-dir`,
and the subdirectory name is the genome name (it must match the name in your
`--lists` / `--list-dir` file exactly — case-sensitive, no normalisation):

```
input/
  Nostoc_sp_PCC_7107/
    <anything>.gbff   (or .gbk)      required — exactly one
    <anything>.faa                   required — the primary protein FASTA
    <anything>.ffn                   optional — nucleotide source
```

Filenames are not inspected, so this works with Bakta, Prokka, PGAP, or any
annotator that produces a GenBank file and a protein FASTA — provided you put
one genome per directory.

**Multiple `.faa` files.** Some annotators write more than one protein FASTA.
The tool keeps the primary one (whose headers key on the same `locus_tag` used
in the GenBank CDS features) and ignores known secondary files by name:

| Annotator | kept | ignored |
|---|---|---|
| Bakta | `<name>.faa` | `<name>.hypotheticals.faa`, `<name>_proteins.faa` |
| PGAP | `annot.faa` | `annot_translated_cds.faa` |
| Prokka | `<name>.faa` | (none — single `.faa`) |

The ignore list matches the markers `hypotheticals`, `_proteins` and
`translated_cds` as a substring of the filename. If, after ignoring those, a genome directory still
contains more than one `.faa` (or more than one `.gbff`/`.gbk`), the tool stops
with an error naming the files rather than guessing — picking the wrong protein
FASTA would produce a coordinate-wrong genome DB that only fails later in
CORASON. If you use an annotator not in the table above and it emits a
secondary `.faa`, remove or rename the extra file so exactly one remains.

**`.faa` headers must key on `locus_tag`.** The primary `.faa` and the GenBank
`.gbff` are joined by `locus_tag` (the first token of each `.faa` header must
equal a `/locus_tag` in the GenBank CDS features). Bakta, Prokka and PGAP's
`annot.faa` all satisfy this. After building one genome, a quick sanity check:

```bash
awk -F'\t' 'NR>1 && ($5==0 || $12=="")' evomining_db/GENOMES/100001.txt | head
```

Rows printed mean the join or nucleotide extraction failed for that annotator's
output, and the layout or file choice needs adjusting.

**Expansion threshold.** `analyze --sd` sets it; the default is `mean + 1*SD`,
matching `evominingAnalysis_1SD.py`. Pass `--sd 2` for the mean + 2SD variant
described in the paper. Note that re-running `analyze` with a different `--sd` overwrites
`evomining_analysis/`; use `-o` to keep both.

**MIBiG.** Only `analyze` BLASTs against the MIBiG FASTA. `trees` reads it too,
but only to place reference sequences keyed off `analyze`'s blast output — so
pass the *same* MIBiG file to both. `start` never touches MIBiG; its `--np-db`
argument only labels the output directory.

**The clade map is optional.** `generate-genome-db` only writes
`genome_clade_mapping.tsv` when `--clade-map` is passed. It is metadata used for
iTOL colorstrips and heatmap annotation, not a pipeline input. Without it,
`trees` skips `itol_clades.txt` and carries on.

**Protein functions come from the 13-column `.txt`.** Since `generate-genome-db`
started writing bare `.faa` headers (`>fig|666666.100021.peg.1`), the product
description only exists in column 7 of `GENOMES/<rast_id>.txt`. BLAST tabular
output truncates subject IDs at whitespace, so it cannot carry the function
either. `analyze` now reads column 7 directly, so `classified_proteins.tsv` and
the tree leaf labels carry real annotations rather than `unknown`.

**Coordinates encode strand by order.** Minus-strand genes are written
`start > stop`, matching the RAST convention CORASON's `3_Draw.pl` expects.
Anything measuring a gene must normalise first — `evomining.rastio.span()` does
this. Never use `stop - start`.

## What changed from the standalone scripts

| Was | Now |
|---|---|
| `generateGenomeDB_corason.py` | `evomining generate-genome-db` |
| `generateEnzymeDB.py` | `evomining generate-enzyme-db` |
| `generateAntismashDB.py` | `evomining generate-antismash-db` |
| `startEvoMining.py` | `evomining start` |
| `evominingAnalysis.py` / `evominingAnalysis_1SD.py` | `evomining analyze --sd 1\|2` |
| `evominingTrees.py` | `evomining trees` |

Pipeline logic is unchanged: copy counting (bitscore ≥ 100, coverage > 50%,
each protein assigned to its best-matching family), BBH, the `color_tree.pl`
green-clade expansion, the preserved `"Secondary metabolsim (antiSMASH)"` label,
and the `CLADE_COLORS` palette all carry over verbatim.

Four things did change:

1. Hardcoded conda prefixes are gone; tools resolve from `$PATH`.
2. Paths resolve by convention from the working directory instead of being
   passed to every command.
3. The 1SD/2SD fork is a flag rather than two files kept in sync by hand.
4. Protein functions are read from the RAST `.txt` (see above).

`classify_gator_windows.py`, `evominingtoGatorGC_perEFDB.py` and the heatmap
scripts are deliberately outside this package — they run after GATOR-GC and
belong with it.
