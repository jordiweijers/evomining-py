"""
evomining.workspace
===================
Convention over configuration. There is no project file and no discovery: the
current working directory *is* the workspace, and every step reads and writes
paths at fixed, predictable locations beneath it.

    <cwd>/
      evomining_db/
        GENOMES.fasta              genome protein DB (composite <stem>__<locus_tag> IDs)
        genome_names.tsv           protein_id -> organism display name
        genome_functions.tsv       protein_id -> product
        enzymes_db.faa             Central DB
      antismash_db.tsv             antiSMASH NP mapping (optional)
      evomining_out/<run>/blast/thereand.blast, backagain.blast
      evomining_analysis/          copy counts, expansions, classified proteins
      evomining_trees/             per-family trees

So `analyze` looks for GENOMES.fasta / genome_names.tsv / genome_functions.tsv and
the BLAST output where `generate-genome-db` and `start` put them, and only needs a
flag when you deviate from the layout. Nothing is remembered between commands
except the files the commands themselves produce.

Clade / iTOL metadata is deliberately not part of this layout: it is a
study-specific overlay produced by a separate tool, not by the core pipeline.
"""
from __future__ import annotations
from pathlib import Path

# Conventional layout, relative to the working directory.
DB_DIR = "evomining_db"
GENOMES_FASTA = "evomining_db/GENOMES.fasta"
GENOME_NAMES = "evomining_db/genome_names.tsv"
GENOME_FUNCTIONS = "evomining_db/genome_functions.tsv"
ENZYME_DB = "evomining_db/enzymes_db.faa"
ANTISMASH_DB = "evomining_db/antismash_db.tsv"
BLAST_ROOT = "evomining_out"
ANALYSIS_DIR = "evomining_analysis"
TREES_DIR = "evomining_trees"
MIBIG_BLAST = "evomining_analysis/expanded_vs_mibig.blast"

# Legacy RAST-era location, kept only for `generate-enzyme-db`'s --tsv mode,
# which still consumes a per-genome GENOMES/ directory.
GENOMES_DIR = "evomining_db/GENOMES"


class WorkspaceError(RuntimeError):
    pass


def default(name, root=None):
    """Absolute path for a conventional location under `root` (default: CWD)."""
    root = Path(root or Path.cwd())
    return (root / globals()[name]).resolve()


def resolve(override, conventional_name, root=None):
    """CLI override wins; otherwise the conventional path under the workspace."""
    if override:
        return Path(override).resolve()
    return default(conventional_name, root)


def require_dir(path, *, flag, produced_by, what):
    path = Path(path)
    if not path.is_dir():
        raise WorkspaceError(
            f"{what} not found at:\n  {path}\n"
            f"Run `evomining {produced_by}` first, or pass {flag} explicitly."
        )
    return path


def require_file(path, *, flag, produced_by, what):
    path = Path(path)
    if not path.is_file():
        raise WorkspaceError(
            f"{what} not found at:\n  {path}\n"
            f"Run `evomining {produced_by}` first, or pass {flag} explicitly."
        )
    return path


def find_blast_dir(override, root=None):
    """
    Locate the blast/ directory produced by `start`.
    `start` writes evomining_out/<central>_<genomes>/blast/. We glob for
    exactly one such directory. Zero -> tell the user to run start; more than
    one -> ambiguous, make them choose with --blast-dir.
    """
    if override:
        return require_dir(override, flag="--blast-dir", produced_by="start",
                           what="BLAST output directory")
    root = Path(root or Path.cwd())
    candidates = sorted((root / BLAST_ROOT).glob("*/blast"))
    candidates = [c for c in candidates if list(c.glob("thereand*.blast"))]
    if not candidates:
        raise WorkspaceError(
            f"no BLAST output found under {root / BLAST_ROOT}/.\n"
            f"Run `evomining start` first, or pass --blast-dir explicitly."
        )
    if len(candidates) > 1:
        listing = "\n  ".join(str(c) for c in candidates)
        raise WorkspaceError(
            "multiple BLAST runs found — pass --blast-dir to pick one:\n  "
            + listing
        )
    return candidates[0]
