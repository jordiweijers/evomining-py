# analysis.py class refactor — design

Plan for restructuring `src/evomining/analysis.py` from nested dicts threaded between free
functions into a small class hierarchy. **Nothing here is implemented yet.**

## Why

`analyze` currently passes `copy_counts`, `family_members`, `bbh_set`, `expansions` and
`classified` between functions as parallel dicts. Two problems:

- `copy_counts[family][genome]` and `family_members[family][genome]` hold the same
  information twice — the count is just `len(members)`.
- `compute_copy_counts` and `compute_bbh` each iterate the forward BLAST independently and
  re-parse every header.

The goal is a clearer view of the data, not speed.

## The classes

### `BlastHit`
One row of BLAST outfmt 6. **Direction-agnostic**: it does not know whether its `query` is
a Central DB seed or a genome protein, which is why `query`/`subject` stay raw strings
rather than typed IDs. Carries the hit-local filter arithmetic (`bitscore`, and
`send*100/qend` with its `qend == 0` guard) as a predicate that takes thresholds as
arguments — the cutoff values belong to the reader, not to the class.

Only 7 of the 12 outfmt-6 columns are read anywhere in `src/`: `query`, `subject`,
`pident`, `qend`, `send`, `evalue`, `bitscore`. Never read: `length`, `mismatch`,
`gapopen`, `qstart`, `sstart`.

### `Member`
One seed–protein match:

```python
protein_id: str
seed_id:    str
forward:    BlastHit
reverse:    BlastHit | None
```

A protein matching several seeds of the same family produces **several Members in one
Cell**, so `len(cell.members)` is *not* the copy count — `Cell.n_copies()` deduplicates by
`protein_id`. Because no Member elects a representative hit, the question of which forward
hit "counts" never arises; the Cell maxes over every match.

`protein_id` and `seed_id` are **assigned by the reader, never derived from
`query`/`subject`**. The orientation convention is guaranteed only by our own `start.py`
(forward is `-query central_db -db genome_db`, reverse the other way round), while
`analyze --blast-dir` accepts any files. Deriving identity from position would make a
swapped file fail silently — see *Behaviour changes* below.

No family field of any kind. Family belongs to the Cell.

No proxied `bitscore`/`evalue`: a Member holds two hits, so a bare `member.bitscore` would
be ambiguous. Callers say `member.forward.bitscore`.

### `Cell`
```python
family_id: str
genome_id: str
members:   list[Member]
n_copies() -> int
bbh()      -> Member | None
```

Owns the **entire** BBH computation — both the forward-best election and the reciprocity
comparison — because the family lives here. `bbh()` returns the Member rather than a
protein id, so the winning seed and both hits come along.

> **Footgun:** sibling Members share a protein, so callers must compare
> `cell.bbh().protein_id == protein_id`, never object identity.

Deliberately **no `is_expanded`**: expansion needs every genome's count for the family, so
it is a Matrix-level question. Keeping it off the Cell also means a Cell is fully
determined at construction and can be frozen.

### `Matrix`
```python
cells:      dict[(genome_id, family_id), Cell]
genome_ids
family_ids
count(genome_id, family_id) -> int      # 0 when absent
expansion(family_id) -> ExpansionStats
expansions()         -> list[family_id]
```


`genome_ids` **must** be stored — it cannot be derived from the cells. A genome with no
hits has no Cell, but its zeros are what set the expansion threshold. Deriving the axis
from the cells would compute mean and std over only the genomes that hit, raising both and
silently under-calling expansions.

A flat dict beats a list because two writers walk the dense genome × family grid, probing a
possibly-absent cell at every position; with a list you would end up building the dict
internally anyway.

### `ExpansionStats`
```python
family_id, mean, std, threshold, max
expanded: dict[genome_id, int]      # inserted in sorted genome order
n_expanded -> len(expanded)
```

`expanded` is a mapping rather than a sequence of pairs, which is what lets a separate
`is_expanded()` disappear: the summary writer gets ordered pairs from `.items()`, the
classified writer gets O(1) membership, and the count is `len()`.

`min` is dropped — `compute_expansions` writes it today and nothing ever reads it.

## Design principles

1. A calculation belongs to the smallest class that already holds all the data it needs.
2. Context flows **downward as immutable data**, never upward via parent pointers. There
   are no back-references in the hierarchy.
3. Don't denormalise a matrix axis into its payload. A Member is the intersection of a
   protein and a family; the family half of its identity is expressed by which Cell it
   sits in.
4. Identity is assigned by the layer that knows the convention, never parsed downstream.

## The reader

Three jobs: apply the direction-appropriate filter, verify orientation, and label which
side of each hit is which.

The filter is a **parameter, not a constant**. The heatplot.pl cutoffs are only meaningful
on the forward file: in reverse, `query` and `subject` are swapped, so `send*100/qend`
measures something different, and a `bitscore >= 100` floor would discard weak reverse hits
that currently still count for reciprocity.

Filtering at read time also settles where BBH election happens — only counted members are
ever candidates.

Note `send * 100 / qend` mixes a subject coordinate with a query one, so it is not a
coverage measure despite looking like one. Naming it `min_coverage` would be misleading.

## What stays outside the hierarchy

The four metadata loaders, `blast_vs_mibig` (subprocess + FASTA IO), all four writers, and
`run()`.

**antiSMASH membership, `function`, `genome_name` and MIBiG hits are write-time joins, not
fields.** Keeping antiSMASH out of the Matrix turns its "cosmetic only — must never
influence counts, expansion, BBH or MIBiG" rule from a comment into a structural
guarantee.

`trees.py` is not a Matrix consumer; it reads `classified_proteins.tsv` back off disk.

## The expansion algorithm

Two passes over `cells()`, never materialising the dense grid.

Pass one accumulates `Σx`, `Σx²` and `max` per family — O(families) memory regardless of
cell count. Then:

```
mean      = Σx / len(genome_ids)
std       = sqrt(Σx²/len(genome_ids) - mean²)
threshold = mean + sd * std
```

Pass two collects the expanded genomes once thresholds exist.

Two facts make this work: the genome count comes from the axis, so absent genomes count as
zeros without being materialised; and **a zero-count genome can never be expanded**
(`threshold` is always ≥ 0, so `0 > threshold` is never true), so pass two only needs to
touch cells that exist.

`sd` is set at construction. The Matrix already takes `method` there, and nothing in the
codebase varies `sd` within a run.

## Behaviour changes vs today

- **The cross-family BBH leak is fixed structurally.** `compute_bbh` currently returns a
  flat protein-ID set tested by plain membership, so under `multi` counting a protein that
  earns BBH in one family reads `is_bbh=True` in *every* family row it appears in. With BBH
  on the Cell, a protein's reverse hit lands in exactly one family, so at most one of its
  cells can match.
- **`evalue`/`bitscore` in `classified_proteins.tsv`** become the best hit's rather than
  the first-seen hit's. Free: the golden snapshot already drops both columns as
  BLAST-version-dependent.
- **BBH election** moves from all forward hits (unfiltered) to counted members only.

## Open decisions

1. **Family axis dense or sparse?** Today it is sparse (derived from `copy_counts.keys()`),
   so a Central DB family with zero hits vanishes from every output — no summary row, no
   heatplot column, no evidence it was searched for. A dense axis would need `family_ids`
   stored. Adopting it is free against the current golden files, since both fixture
   families have hits.
2. **`BlastHit` columns and representation** — all 12 or the 7 that are used, and plain
   dataclass vs `NamedTuple` vs manual `__slots__`. Note `requires-python = ">=3.9"` rules
   out `@dataclass(slots=True)`.
3. **What to call the filter predicate**, given the coordinate-mixing issue above.

## Verified against the los17 fixture

Measured on this branch (17 genomes, 2 families, 443 forward / 508 reverse hits):

| claim | result |
|---|---|
| member rows (family, genome, protein) | 89 |
| forward hits per member | mean 3.21, max 4, 286 retained |
| members whose first-seen hit is not the best-scoring | 68 of 89 |
| BBH set size / of those, not a counted copy | 34 / 0 |
| cells electing a non-counted protein | 0 |
| reciprocity passed / rejected | 34 / 0 |
| running sums vs numpy mean, std | identical; std within 8.9e-16 |
| swapped BLAST orientation | 508 rows parsed → 0 families, 0 members, **no exception** |

That last row is why `protein_id`/`seed_id` are assigned rather than derived: pointing
`analyze` at files with the wrong orientation produces an empty analysis and a clean exit.

Expansion behaviour on the fixture:

```
3PGA_AMINOACIDS|1  counts 9 8 8 7 6 6 5 4 3 3 3 2 2 2 2 1 1
                   mean 4.2353  std 2.5559
                   sd=1 -> threshold 6.7912 -> 4 expanded
                   sd=2 -> threshold 9.3472 -> 0 expanded

3PGA_AMINOACIDS|2  counts all 1
                   mean 1.0  std 0.0  threshold 1.0 -> 0 expanded
```

`|2` is why the expansion test is `>` and not `>=`: with `>=`, a perfectly flat family
would report every genome as expanded.

## Fixture blind spots

The fixture has only two families (`3PGA_AMINOACIDS|1` PGDH and `|2` PSAT) and **no protein
appears in both**, so it exercises neither the cross-family BBH leak nor the election gap.
Reciprocity is effectively untested too — all 34 elected proteins reciprocate, 0
rejections, because those two enzymes never compete for the same protein. No genome has a
zero count in either family, so the dense-genome-axis requirement is unexercised as well.

Snapshot tests will not catch a regression in any of those. Worth adding targeted cases
alongside the refactor.
