# EvoMining

Enzyme family expansion analysis for biosynthetic gene cluster discovery — a
Python reimplementation of EvoMining (Selem-Mojica et al. 2019), packaged as a
single conda environment and a single `evomining` command.

Genomes are read directly from GenBank (`.gbff`/`.gbk`); the parser is
vendored/adapted from corason-py, and every artifact is keyed on a composite
`<genome_stem>__<locus_tag>` ID. GATOR-GC, CORASON and antiSMASH stay in their
own environments — they are separate steps in the workflow rather than part of
EvoMining proper.

## Installation

Clone the repository and build the environment from `environment.yml`. This
installs the external tools (`blast`, `muscle` v5, `mafft`, `trimal`,
`fasttree`, `newick_utils`) **and** the `evomining` package in one step:

```bash
git clone https://github.com/bscheep/evomining-py.git
cd evomining-py
micromamba env create -f environment.yml
micromamba activate evomining
evomining check           # verify every external tool resolves
```

If you already have the external tools on your `$PATH` and only want the Python
package in an existing environment:

```bash
pip install "git+https://github.com/bscheep/evomining-py.git@main"
```

Nothing is hardcoded to a conda prefix — tools are taken from `$PATH`. To force
one tool to come from a different environment (e.g. an older BLAST to reproduce
a previous run):

```bash
evomining check --tool-env blastp=/vol/local/conda_envs/blast+
```

## Usage

There is no project file and no `init` step: the directory you run commands from
**is** the workspace. Each step writes to a fixed location beneath it and the
next step reads from there, so after `generate-genome-db` you rarely pass any
paths. Run the six steps in order from one folder:

```bash
evomining generate-genome-db -i <genomes_dir> --list-dir <lists_dir>
evomining generate-enzyme-db --fasta-dir <enzyme_families_dir>
evomining generate-antismash-db --antismash-dir <antismash_dir>   # optional
evomining start
evomining analyze --mibig MiBIG_DB.faa
evomining trees   --mibig MiBIG_DB.faa
```

Every command accepts `--help` for its full flag list. The sections below cover
each step.

### 1. `generate-genome-db` — build the genome protein DB

Reads a folder of GenBank genomes and writes the composite-ID protein database
the rest of the pipeline consumes: `GENOMES.fasta` (one record per CDS, header =
`<genome_stem>__<locus_tag>`), `genome_names.tsv` (id → organism name) and
`genome_functions.tsv` (id → product).

```bash
evomining generate-genome-db -i <genomes_dir> --list-dir <lists_dir>
```

- **`-i, --input-dir`** — directory with one subdirectory per genome. The
  subdirectory name is the genome name and must match the names in your list
  files exactly (case-sensitive). Each subdirectory holds exactly one
  `.gbff`/`.gbk` and one primary `.faa` (see *Genome input layout* below).
- **`-l, --lists`** — one or more `.txt` files, each listing genome stems to
  include (one per line), space-separated on the flag:
  `--lists cladeA.txt cladeB.txt`. Only genomes named in these lists are built.
  A listed stem with no matching GenBank file is warned about and skipped.
- **`--list-dir`** — a directory of such `.txt` files (all `*.txt` in it are
  read); every stem listed across them is included. A single `.txt` or many both
  work — the lists' only required job is to select which genomes enter the DB.
  (The file basenames are additionally reused as clade labels for downstream
  iTOL/heatmap annotation, but that is incidental and optional.)
- **Selection is optional.** If neither `-l` nor `--list-dir` is given, **every**
  GenBank genome under `--input-dir` is built. The lists are only a subset
  filter; use them when you want fewer than all genomes.
- **`--names`** — optional TSV overriding organism display names, as
  `<genome_stem><TAB>organism_name` (the first column is the key, the last column
  is the name, so wider corason-style tables also work). Usually unnecessary:
  in GenBank mode names are taken from each file's SOURCE/ORGANISM record, and
  `--names` only overrides that for genomes present in the table.
- **`--keep-pseudogenes`** — pseudogene CDS are dropped by default (they are not
  functional enzymes and inflate enzyme-family copy counts). Pass this to retain
  them, e.g. for continuity with older RAST-based runs.
- **`-o, --output-dir`** — default `./evomining_db`.

### 2. `generate-enzyme-db` — build the Central (enzyme) DB

Builds the Central DB of seed enzyme families that genomes are searched against.
It accepts four input modes; the usual one is `--fasta-dir` (one file per enzyme
family). Pick the mode matching what you have:

```bash
evomining generate-enzyme-db --fasta-dir <enzyme_families_dir>
```

- **`--fasta-dir DIR`** — a directory of FASTA files (`.faa`, `.fasta` or `.fa`),
  one file per enzyme family (all sequences in a file are that family's seeds).
  The **filename** encodes the family — `SUBSYSTEM__enzyme.faa` groups the family
  under a subsystem (the part before `__`), a bare `enzyme.faa` uses
  `--custom-pathway-name`. The **family grouping comes entirely from the
  filename**, so the sequence headers inside need no special format — they are
  reused only as per-seed labels, and the tool synthesizes the canonical
  `SUBSYSTEM|family|Function|Organism` header itself. This is the mode to use
  when you have plain per-family FASTAs (e.g. seeds pulled from UniProt/NCBI) and
  do not want to hand-format canonical headers; give sequences meaningful header
  IDs if you want readable seed labels later.
- **`--tsv FILE`** — a TSV (`protein_id, family, enzyme_name, [subsystem]`) whose
  `protein_id`s are pulled from the genome FASTAs in `--genomes-dir`. The
  optional 4th column groups families under a subsystem.
- **`--transcripts FILE`** — specialized transcript workflow: a cleaned
  transcript FASTA (`TRINITY_ID|gene|host|description`) where each transcript
  becomes its own single-seed `custom|N` family. Use this to seed from expressed
  transcripts rather than a curated enzyme set.
- **`--genomes-dir`** — genome FASTA directory, required only for `--tsv` mode.
- **`--custom-pathway-name`** — subsystem name for families that have no
  subsystem of their own (bare `--fasta-dir` filenames, `--tsv` rows without a
  subsystem column, and all `--transcripts` families). Replaces the `custom` in
  `custom|N`. Families that already carry a subsystem keep it. Default `custom`.
- **`--fasta FILE`** — a FASTA already in canonical EvoMining header format,
  `SUBSYSTEM|family_number|Function_querynumber|Organism`; headers are passed
  through verbatim. Mostly useful for re-using or verifying the format of an
  existing Central DB (e.g. one produced by a previous run), rather than as the
  normal way to build one.
- **`-o, --output`** — default `./evomining_db/enzymes_db.faa`.

### 3. `generate-antismash-db` — antiSMASH NP mapping (optional)

Builds an optional mapping from antiSMASH results. This is **cosmetic**: it only
sets a leaf label and an `in_antismash` flag downstream and never gates
predictions — the candidate set is identical with and without it.

```bash
evomining generate-antismash-db --antismash-dir <antismash_dir>
```

- **`--antismash-dir`** — antiSMASH results directory (per-genome subdirectories,
  optionally grouped under clade subdirectories).
- **`--genome-names`** — `genome_names.tsv` from step 1, used to recognise valid
  genome stems. Default `./evomining_db/genome_names.tsv`.
- **`-o, --output`** — default `./evomining_db/antismash_db.tsv`.

### 4. `start` — forward + reverse BLAST

Runs the reciprocal BLAST between the Central DB and the genome DB (forward:
Central vs genomes → `thereand.blast`; reverse: genomes vs Central →
`backagain.blast`), the input to expansion and BBH analysis.

```bash
evomining start
```

- **`-g, --genomes`** — genome protein FASTA. Default
  `./evomining_db/GENOMES.fasta`.
- **`-c, --central-db`** — Central DB. Default `./evomining_db/enzymes_db.faa`.
- **`-o, --output-dir`** — default `./evomining_out`.
- **`--threads`** — BLAST threads. Default 32.
- **`--evalue`** — BLAST e-value. Default `1e-4`.
- **`--chunk-size`** — reverse BLAST query block size in proteins. Lower it if
  `blastp` fails to create threads on a very large query. Default 100000.

`start` never touches MIBiG.

### 5. `analyze` — copy counts, expansions, classification

Computes per-family copy counts, expansion thresholds, protein classification
(writing the authoritative `classified_proteins.tsv`), and BLASTs the expanded
proteins against MIBiG.

```bash
evomining analyze --mibig MiBIG_DB.faa
```

- **`--mibig`** — MIBiG protein FASTA for the natural-product BLAST.
- **`--sd`** — expansion threshold, `mean + SD*stdev`. Default `1.0`; pass
  `--sd 2` for the mean + 2SD variant described in the paper.
- **`--ccm {multi,best}`** — copy-count method. `multi` (default, original
  behaviour) counts a protein in every family it matches; `best` assigns each
  protein to its single best-scoring family only.
- **`--antismash`** — the antiSMASH mapping from step 3 (cosmetic; see above).
- **`--blast-dir`** — defaults to the single run under `./evomining_out`.
- **`--central-db`, `--genome-names`, `--genome-functions`, `--genomes-fasta`** —
  default to the step 1/4 artifacts under `./evomining_db` and `./evomining_out`.
- **`-o, --output-dir`** — default `./evomining_analysis`. Re-running with a
  different `--sd` overwrites this directory; use `-o` to keep both variants.

### 6. `trees` — per-family coloured phylogenetic trees

Builds one coloured phylogenetic tree per expanded family, placing MIBiG
reference sequences as expansion seeds and colouring leaves by classification.

```bash
evomining trees --mibig MiBIG_DB.faa
```

- **`--mibig`** — MIBiG FASTA. Must be the **same** file passed to `analyze`;
  `trees` uses it only to place reference sequences keyed off `analyze`'s BLAST
  output.
- **`--mibig-blast`** — `analyze`'s MIBiG BLAST. Default
  `./evomining_analysis/expanded_vs_mibig.blast`.
- **`--classified`** — the classified-proteins table. Default
  `./evomining_analysis/classified_proteins.tsv`.
- **`--families`** — restrict to these families, or a `.txt` with one family per
  line.
- **`--itol`** — also write `itol_colors.txt` (iTOL `TREE_COLORS`) colouring each
  leaf's branch/label by classification.
- **`--microreact`** — also write `microreact.csv` (id, classification, colour)
  for loading into microreact.org alongside the `.nwk` tree.
- **`--max-seqs`** — skip families with more than this many sequences. Default:
  no limit.
- **`--skip-existing`** — resume without recomputing families already done.
- **`--central-db`, `--genomes-fasta`, `--threads`, `-o`** — the usual defaults
  under `./evomining_db` and `./evomining_trees`.

### Genome input layout

`generate-genome-db` expects one subdirectory per genome under `--input-dir`.
The subdirectory name is the genome name and must match the name in your list
files exactly (case-sensitive):

```
input/
  Nostoc_sp_PCC_7107/
    <anything>.gbff   (or .gbk)      required — exactly one
    <anything>.faa                   required — the primary protein FASTA
    <anything>.ffn                   optional — nucleotide source
```

Filenames are not inspected, so this works with Bakta, Prokka, PGAP, or any
annotator producing a GenBank file and a protein FASTA — one genome per
directory. The `.faa` and `.gbff` are joined by `locus_tag` (the first token of
each `.faa` header must equal a `/locus_tag` in the GenBank CDS features), which
Bakta, Prokka and PGAP's `annot.faa` all satisfy.

When an annotator writes more than one protein FASTA, the tool keeps the primary
one and ignores known secondaries by name (markers `hypotheticals`, `_proteins`,
`translated_cds`):

| Annotator | kept | ignored |
|---|---|---|
| Bakta | `<name>.faa` | `<name>.hypotheticals.faa`, `<name>_proteins.faa` |
| PGAP | `annot.faa` | `annot_translated_cds.faa` |
| Prokka | `<name>.faa` | (none) |

If more than one `.faa` (or `.gbff`/`.gbk`) still remains, the tool stops with an
error naming the files rather than guessing.

## Changes from the original EvoMining

The original EvoMining (EvoMining 2.0, Selem-Mojica et al. 2019) is a Perl
pipeline distributed as a Docker image, driven through a web interface and a
single `perl startEvoMining.pl` invocation:

```
perl startEvoMining.pl -g <genome-DB> -r <myRastIds> -c <central-DB> \
                       -n <natural-DB> -a <antismash_db>
```

It required genomes to be annotated through the RAST platform (the `-r
<myRastIds>` argument), bundled its BLAST, MUSCLE, Gblocks, FastTree and Newick
Utilities dependencies inside the container, and included a simplified copy of
CORASON for genomic-neighbourhood visualisation.

This reimplementation keeps the same underlying method — following central
enzyme families into expansions and recruitments toward secondary metabolism,
copy counting, bidirectional best hits, the green-clade expansion, and coloured
trees — but changes how it is run and what it consumes:

1. **GenBank input instead of RAST.** Genomes are parsed straight from GenBank
   (`.gbff`/`.gbk`) with composite `<genome_stem>__<locus_tag>` IDs. There is no
   RAST dependency and no `myRastIds` file; any annotator that emits GenBank and
   a protein FASTA (Bakta, Prokka, PGAP, …) works. Strand is a proper `+1`/`-1`
   field rather than being encoded by coordinate order.
2. **A pip/conda package instead of a Docker image.** The pipeline installs as a
   normal conda environment plus an `evomining` command, rather than running
   inside a container. External tools resolve from `$PATH`; none are hardcoded.
3. **A command-line pipeline instead of a web interface.** Each stage is an
   `evomining` subcommand run from the terminal; there is no web server or
   browser step. Trees can still be exported for iTOL (`--itol`) and Microreact
   (`--microreact`).
4. **Convention over configuration.** The working directory is the workspace and
   stages find each other's outputs automatically, so the long argument string
   of `startEvoMining.pl` is replaced by defaults you only override when you
   deviate from the layout.
5. **Customizable seeds, including transcripts.** Beyond a curated Central DB,
   `generate-enzyme-db` can seed families from a per-family FASTA directory, a
   TSV, or expressed transcripts (`--transcripts`), each transcript becoming its
   own `custom|N` family.

Note on alignment: the original wrapped MUSCLE with Gblocks for alignment
cleanup; this version uses MUSCLE (falling back to MAFFT) with trimAl.

`classify_gator_windows.py`, `evominingtoGatorGC_perEFDB.py` and the heatmap
scripts are deliberately outside this package — they run after GATOR-GC and
belong with it.
