# EvoMining

> **Pseudogenes.** `generate-genome-db` drops pseudogene CDS by default (they are
> not functional enzymes and would inflate enzyme-family copy counts). Pass
> `--keep-pseudogenes` to retain them, e.g. for continuity with older RAST-based
> runs, which kept them.

Enzyme family expansion analysis for biosynthetic gene cluster discovery — a
Python reimplementation of EvoMining (Selem-Mojica et al. 2019), packaged as a
single conda environment and a single `evomining` command.

Genomes are read directly from GenBank (`.gbff`/`.gbk`); the parser is
vendored/adapted from corason-py, and every artifact is keyed on a composite
`<genome_stem>__<locus_tag>` ID (no RAST, no `fig|` scheme).

GATOR-GC, CORASON and antiSMASH stay in their own environments. They pin
conflicting Python stacks, and they are separate steps in the workflow rather
than part of EvoMining proper.

## Install

Clone the repository and build the environment from `environment.yml` — this
installs the external tools (`blast`, `muscle` v5, `mafft`, `trimal`, `fasttree`,
`newick_utils`) *and* the `evomining` package in one step:

```bash
git clone https://github.com/bscheep/evomining-py.git
cd evomining-py
micromamba env create -f environment.yml
micromamba activate evomining
evomining check           # verify every external tool resolves
```

If you already have the external tools on your `$PATH` and only want the Python
package, install it into your existing environment:

```bash
pip install "git+https://github.com/bscheep/evomining-py.git@main"
```

Nothing is hardcoded to `/vol/local/conda_envs/...` — tools are taken from
`$PATH`. If a tool has to come from a different environment (an old BLAST to
reproduce a previous run, say), point at it explicitly:

```bash
evomining check --tool-env blastp=/vol/local/conda_envs/blast+
```

which makes that one call run as `micromamba run -p <prefix> blastp ...`.

## The working directory is the workspace

There is no project file and no `init` step. The directory you run commands from
is the workspace, and every step reads and writes at fixed, predictable
locations beneath it:

```
evomining_db/GENOMES.fasta         genome protein DB (composite <stem>__<locus_tag> IDs)
evomining_db/genome_names.tsv      composite_id -> organism display name
evomining_db/genome_functions.tsv  composite_id -> product
evomining_db/enzymes_db.faa        Central (enzyme) DB
evomining_db/antismash_db.tsv      antiSMASH NP mapping (optional)
evomining_out/<run>/blast/         forward (thereand.blast) + reverse (backagain.blast)
evomining_analysis/                copy counts, expansions, classified_proteins.tsv
evomining_trees/                   per-family trees
```

So `analyze` looks for `GENOMES.fasta`, `genome_names.tsv`, `genome_functions.tsv`
and the BLAST output where `generate-genome-db` and `start` put them, checks they
are present, and only needs a flag when you deviate from the layout. If a
prerequisite is missing, the error names the file and the command that produces
it. Run the steps in order from one folder and no paths need repeating.

## Workflow

```bash
# 1. genome DB — GenBank (.gbff/.gbk) folders -> composite-ID protein DB
evomining generate-genome-db \
    -i /vol/databases/genomes_cyanobacteria_Nostocales/genomes_annotated \
    --list-dir ../0.clades/symb_clades   # .txt lists; basename = clade label

# 2. Central (enzyme) DB — from a FASTA already in canonical EvoMining header format
#    (SUBSYSTEM|family_number|Function_querynumber|Organism)
evomining generate-enzyme-db --fasta central_enzymes.faa
#    other input modes: --fasta-dir (one file per family), --tsv, --transcripts

# 3. antiSMASH NP mapping (optional; cosmetic — labels leaves, never gates predictions)
evomining generate-antismash-db --antismash-dir ../6.antismash

# 4. run — MIBiG is supplied to analyze (and to trees, which reads analyze's blast)
evomining start
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
`translated_cds` as a substring of the filename. If, after ignoring those, a
genome directory still contains more than one `.faa` (or more than one
`.gbff`/`.gbk`), the tool stops with an error naming the files rather than
guessing — picking the wrong protein FASTA would produce a coordinate-wrong
genome DB that only fails later. If you use an annotator not in the table above
and it emits a secondary `.faa`, remove or rename the extra file so exactly one
remains.

**`.faa` headers must key on `locus_tag`.** The primary `.faa` and the GenBank
`.gbff` are joined by `locus_tag` (the first token of each `.faa` header must
equal a `/locus_tag` in the GenBank CDS features). Bakta, Prokka and PGAP's
`annot.faa` all satisfy this.

**Clade labels come from `--list-dir`.** Each `.txt` file under `--list-dir`
lists genome stems, one per line, and the file's basename is used as the clade
label for those genomes (metadata for iTOL colorstrips and heatmap annotation).
It is not a pipeline gate — without it, per-clade annotation is simply skipped.

**Expansion threshold.** `analyze --sd` sets it; the default is `mean + 1*SD`.
Pass `--sd 2` for the mean + 2SD variant described in the paper. Re-running
`analyze` with a different `--sd` overwrites `evomining_analysis/`; use `-o` to
keep both.

**MIBiG.** Only `analyze` BLASTs against the MIBiG FASTA. `trees` reads it too,
but only to place reference sequences keyed off `analyze`'s blast output — so
pass the *same* MIBiG file to both. `start` never touches MIBiG.

**antiSMASH is cosmetic.** The optional antiSMASH mapping only sets a leaf label
and the `in_antismash` flag; it never gates predictions. The candidate set is
identical with and without it.

**Protein functions come from `genome_functions.tsv`.** `GENOMES.fasta` headers
are bare composite IDs, and BLAST tabular output truncates subject IDs at
whitespace, so neither carries the product description. `analyze` reads
`genome_functions.tsv` directly, so `classified_proteins.tsv` and the tree leaf
labels carry real annotations rather than `unknown`.

## What changed from the standalone scripts

| Was | Now |
|---|---|
| `generateGenomeDB_corason.py` | `evomining generate-genome-db` |
| `generateEnzymeDB.py` | `evomining generate-enzyme-db` |
| `generateAntismashDB.py` | `evomining generate-antismash-db` |
| `startEvoMining.py` | `evomining start` |
| `evominingAnalysis.py` / `evominingAnalysis_1SD.py` | `evomining analyze --sd 1\|2` |
| `evominingTrees.py` | `evomining trees` |

Pipeline logic is unchanged: copy counting (bitscore ≥ 100, coverage > 50%),
BBH, the `color_tree.pl` green-clade expansion, the preserved antiSMASH label,
and the `CLADE_COLORS` palette all carry over.

What changed:

1. Genome input is GenBank with composite `<stem>__<locus_tag>` IDs, not RAST
   `fig|` IDs; `Corason_Rast.IDs` and the 13-column `.txt` are gone.
2. Hardcoded conda prefixes are gone; tools resolve from `$PATH`.
3. Paths resolve by convention from the working directory instead of being
   passed to every command.
4. The 1SD/2SD fork is a flag rather than two files kept in sync by hand.

`classify_gator_windows.py`, `evominingtoGatorGC_perEFDB.py` and the heatmap
scripts are deliberately outside this package — they run after GATOR-GC and
belong with it.
