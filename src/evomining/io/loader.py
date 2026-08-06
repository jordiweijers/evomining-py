"""Genome-set loading and format discovery.

Adapted from corason-py (miguel-mx/corason-py), src/corason/io/loader.py, with the
legacy RAST path removed: EvoMining reads GenBank only, since retiring RAST is the
whole point of this pipeline. ``--genomes`` accepts a directory, a glob, or
explicit paths.
"""

from __future__ import annotations

import logging
from pathlib import Path

from .genbank import GENBANK_SUFFIXES, is_genbank, load_genbank
from .models import Genome
from .names import load_names

logger = logging.getLogger(__name__)


def discover(paths: list[Path]) -> list[Path]:
    """Expand directories into the GenBank files they contain."""
    found: list[Path] = []

    for path in paths:
        if path.is_file():
            found.append(path)
            continue
        if not path.is_dir():
            raise FileNotFoundError(f"{path}: no such file or directory")

        genbanks = sorted(p for p in path.iterdir() if is_genbank(p))
        if genbanks:
            found.extend(genbanks)
            continue

        raise FileNotFoundError(
            f"{path}: no GenBank files found "
            f"({', '.join(sorted(GENBANK_SUFFIXES))})."
        )

    if not found:
        raise FileNotFoundError("no genome files found")
    return found


def load_genomes(paths: list[Path], names_file: Path | None = None,
                 keep_pseudo: bool = False) -> list[Genome]:
    """Load every GenBank genome in ``paths``, applying an optional name override."""
    names = load_names(names_file) if names_file else {}
    genomes: list[Genome] = []

    for path in discover(paths):
        override = names.get(path.stem)
        if not is_genbank(path):
            raise ValueError(
                f"{path}: unrecognized genome format. Expected a GenBank file "
                f"({', '.join(sorted(GENBANK_SUFFIXES))})."
            )
        genomes.append(load_genbank(path, name=override, keep_pseudo=keep_pseudo))

    _disambiguate(genomes)
    return genomes


def _disambiguate(genomes: list[Genome]) -> None:
    """Make organism display names unique.

    Names become tree leaf labels, where duplicates make the tree ambiguous. The
    accession is preferred for disambiguation, the filename stem as fallback.

    NOTE: this touches only ``Genome.name`` (the display name). EvoMining's protein
    IDs are namespaced by ``Genome.id`` (the filename stem), which is unique per
    file regardless of this step, so IDs stay stable across batches.
    """
    groups: dict[str, list[Genome]] = {}
    for genome in genomes:
        groups.setdefault(genome.name, []).append(genome)

    for name, group in groups.items():
        if len(group) == 1:
            continue

        accessions = [g.accessions[0] if g.accessions else None for g in group]
        distinct = all(accessions) and len(set(accessions)) == len(accessions)
        logger.warning(
            "%d genomes share the organism name %r; disambiguating with %s",
            len(group), name, "accession" if distinct else "filename",
        )
        for genome in group:
            suffix = genome.accessions[0] if distinct else genome.id
            genome.name = f"{name} {suffix}"

    _enforce_unique_names(genomes)


def _enforce_unique_names(genomes: list[Genome]) -> None:
    """Last resort if two genomes still collide after accession/filename suffixing."""
    taken: set[str] = set()
    for genome in genomes:
        if genome.name not in taken:
            taken.add(genome.name)
            continue
        ordinal = 2
        while f"{genome.name} {ordinal}" in taken:
            ordinal += 1
        logger.warning("organism name %r is still duplicated; suffixing", genome.name)
        genome.name = f"{genome.name} {ordinal}"
        taken.add(genome.name)
