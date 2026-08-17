"""Genome-set loading and format discovery.

Adapted from corason-py (miguel-mx/corason-py), src/corason/io/loader.py, with the
legacy RAST path removed: EvoMining reads GenBank only, since retiring RAST is the
whole point of this pipeline. ``--genomes`` accepts a directory, a glob, or
explicit paths.
"""

from __future__ import annotations

import logging
from pathlib import Path

from .genbank import (GENBANK_SUFFIXES, genbank_suffix, genome_stem, is_genbank,
                      load_genbank)
from .models import Genome
from .names import load_names

logger = logging.getLogger(__name__)


def _accepted() -> str:
    """The recognised GenBank extensions, for error messages."""
    return ", ".join(sorted(GENBANK_SUFFIXES)) + " (optionally .gz)"


def discover(paths: list[Path]) -> list[Path]:
    """Expand directories into the GenBank files they contain (flat only).

    Kept for callers that already hold a list of files/dirs. New genome-set entry
    points should prefer :func:`resolve_genome_inputs`, which also understands the
    one-folder-per-genome layout and reports the intended stem for each file.
    """
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
            f"{path}: no GenBank files found ({_accepted()})."
        )

    if not found:
        raise FileNotFoundError("no genome files found")
    return found


# Suffix preference when a per-genome folder holds more than one GenBank file.
_GB_PREFERENCE = (".gbff", ".gbk", ".gb", ".genbank")


def _choose_genbank(files: list[Path], label: str) -> Path:
    """Pick the single GenBank file for one genome folder, or fail clearly."""
    for suffix in _GB_PREFERENCE:
        matches = [f for f in files if genbank_suffix(f) == suffix]
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise FileNotFoundError(
                f"{label}: {len(matches)} '{suffix}' files "
                f"({', '.join(sorted(f.name for f in matches))}); "
                f"expected one GenBank file per genome folder."
            )
    if len(files) == 1:
        return files[0]
    raise FileNotFoundError(f"{label}: multiple GenBank files; expected one per genome.")


def resolve_genome_inputs(root: Path) -> list[tuple[Path, str]]:
    """Resolve an input path to ``(genbank_file, genome_stem)`` pairs.

    Accepts, transparently:

    * a single GenBank file                       -> stem = file stem
    * a flat directory of GenBank files           -> stem = file stem
    * a directory of per-genome sub-folders,      -> stem = *folder* name
      each holding that genome's annotation
      (e.g. Bakta output, ``<genome>/<genome>.gbff``)

    and any mix of the last two. In the per-genome-folder layout the stem is the
    folder name -- not the inner filename -- because the folder is the genome's
    identity and the annotation file inside may be named generically or by
    accession. The stem namespaces every protein ID and must line up with the
    antiSMASH directory names, so the resolver enforces stem uniqueness.
    """
    root = Path(root)
    if root.is_file():
        if is_genbank(root):
            return [(root, genome_stem(root))]
        raise FileNotFoundError(f"{root}: not a GenBank file.")
    if not root.is_dir():
        raise FileNotFoundError(f"{root}: no such file or directory")

    pairs: list[tuple[Path, str]] = []

    # flat: GenBank files sitting directly in root -> stem = file stem
    for p in sorted(root.iterdir()):
        if p.is_file() and is_genbank(p):
            pairs.append((p, genome_stem(p)))

    # nested: each immediate sub-folder holding GenBank file(s) -> stem = folder
    for sub in sorted(root.iterdir()):
        if not sub.is_dir():
            continue
        gbs = [q for q in sorted(sub.iterdir()) if q.is_file() and is_genbank(q)]
        if gbs:
            pairs.append((_choose_genbank(gbs, sub.name), sub.name))

    if not pairs:
        raise FileNotFoundError(
            f"{root}: no GenBank files found -- neither directly "
            f"({_accepted()}) nor one level down in per-genome sub-folders."
        )

    seen: dict[str, Path] = {}
    for path, stem in pairs:
        if stem in seen:
            raise ValueError(
                f"duplicate genome stem {stem!r} (from {seen[stem]} and {path}); "
                f"stems namespace protein IDs and must be unique."
            )
        seen[stem] = path
    return pairs


def load_genomes(paths: list[Path], names_file: Path | None = None,
                 keep_pseudo: bool = False,
                 stems: dict[Path, str] | None = None) -> list[Genome]:
    """Load every GenBank genome in ``paths``, applying an optional name override.

    ``stems`` optionally maps a file path to the genome stem it should use (from
    :func:`resolve_genome_inputs`, e.g. the folder name in the per-genome-folder
    layout). When absent, the filename stem is used, preserving flat-layout
    behaviour. The ``--names`` override is keyed on the resolved stem.
    """
    names = load_names(names_file) if names_file else {}
    genomes: list[Genome] = []

    for path in discover(paths):
        stem = stems.get(path) if stems else None
        override = names.get(stem or genome_stem(path))
        if not is_genbank(path):
            raise ValueError(
                f"{path}: unrecognized genome format. Expected a GenBank file "
                f"({_accepted()})."
            )
        genomes.append(load_genbank(path, name=override, keep_pseudo=keep_pseudo,
                                    stem=stem))

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
