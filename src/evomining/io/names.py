"""Optional organism-name override table.

Adapted from corason-py (miguel-mx/corason-py), src/corason/io/names.py, so the
same --names file works for both tools.

In GenBank mode names come from the SOURCE/ORGANISM records, so a table is usually
unnecessary. When GenBank metadata is poor, a tab-separated file overrides it:

    <genome_stem><TAB><organism name>

Column 0 is the key (the genome file's stem) and the *last* column is the display
name, so the original 3-column ``Example.Ids`` form also works.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def load_names(path: Path) -> dict[str, str]:
    """Parse a names table into ``{genome_stem: organism_name}``."""
    names: dict[str, str] = {}
    with open(path, newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or not row[0].strip():
                continue
            key = row[0].strip().removesuffix(".faa")
            if len(row) < 2:
                logger.warning("%s: no name for %r, ignoring row", path.name, key)
                continue
            names[key] = row[-1].strip()
    if not names:
        raise ValueError(f"{path}: no usable rows found")
    return names
