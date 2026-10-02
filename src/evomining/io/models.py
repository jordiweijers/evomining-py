"""Genome data model.

Adapted from corason-py (miguel-mx/corason-py), src/corason/models.py. Trimmed to
the parsing-side types EvoMining needs: Gene, Contig, Genome. The neighborhood /
cluster / orthogroup types (CORASON's synteny machinery) are not carried, since
EvoMining does no neighborhood analysis itself.

Coordinates are 0-based half-open ``[start, end)`` throughout, matching Biopython
and Python slicing. Conversion to the 1-based inclusive convention happens only at
the I/O boundary. Strand is ``+1`` / ``-1``.

A gene's identity is a stable qualifier (``locus_tag`` first), never an integer
parsed out of the identifier. ``Gene.index`` carries the gene's ordinal position
along its contig instead, assigned after sorting by coordinate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator


@dataclass(slots=True)
class Gene:
    """A single CDS feature."""

    id: str
    """Stable identifier: locus_tag, else protein_id, else gene, else synthesized."""

    contig: str
    index: int
    """0-based ordinal among the CDS of this contig, sorted by start coordinate."""

    start: int
    end: int
    strand: int
    product: str
    translation: str

    is_compound: bool = False
    """True when the CDS has a spliced/joined location (introns, ribosomal slippage)."""

    parts: list[tuple[int, int]] = field(default_factory=list)
    """Exon extents for compound features; empty when ``is_compound`` is False."""

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def strand_symbol(self) -> str:
        return "+" if self.strand >= 0 else "-"


@dataclass(slots=True)
class Contig:
    """One replicon or assembly contig, holding its CDS in coordinate order."""

    id: str
    genes: list[Gene] = field(default_factory=list)
    length: int | None = None

    def __len__(self) -> int:
        return len(self.genes)


@dataclass(slots=True)
class GenomeMetadata:
    """Metadata for a genome, separate from its contigs and genes.
    
    Split out from Genome so a clear separation exists between the genome's
    structural data (contigs and genes) and its descriptive metadata.
    This is usefull for the downstream loading and processing of genome metadata 
    independently of the structural genome data.
    """
    id: str
    """Filename stem. Used for output paths, as EvoMining's ID namespace, and as
    the fallback display name."""

    name: str
    """Organism name, from the GenBank SOURCE record or the --names table."""

    accessions: list[str] = field(default_factory=list)
    """List of accession numbers associated with this genome."""

    source_path: Path | None = None
    """Path to the source GenBank file, if available."""


@dataclass(slots=True)
class Genome:
    """One annotated genome: its contigs, and the names it answers to (metadata)."""
    metadata: GenomeMetadata
    contigs: dict[str, Contig] = field(default_factory=dict)

    def __init__(
            self, id: str, 
            name: str, 
            contigs: dict[str, Contig] | None = None,
            accessions: list[str] | None = None,
            source_path: Path | None = None,
            ) -> None:
        self.metadata = GenomeMetadata(
            id=id,
            name=name,
            accessions=accessions if accessions is not None else [],
            source_path=source_path,
        )
        self.contigs = contigs if contigs is not None else {}

    def genes(self) -> Iterator[Gene]:
        for contig in self.contigs.values():
            yield from contig.genes

    @property
    def n_genes(self) -> int:
        return sum(len(c.genes) for c in self.contigs.values())


def finalize_contigs(genome: Genome, by_contig: dict[str, list[Gene]]) -> None:
    """Sort each contig's genes by coordinate and assign their positional index."""
    for contig_id, genes in by_contig.items():
        genes.sort(key=lambda g: (g.start, g.end))
        for index, gene in enumerate(genes):
            gene.index = index
        genome.contigs[contig_id] = Contig(
            id=contig_id,
            genes=genes,
            length=max(g.end for g in genes),
        )
