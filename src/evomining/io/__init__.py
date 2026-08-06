"""Genome input handling for EvoMining.

GenBank flat files in, in-memory Genome objects out. The genome/gene data model
and the GenBank parser are adapted from corason-py (miguel-mx/corason-py), so the
two related tools handle genomes identically. The RAST loader is intentionally
NOT vendored: EvoMining reads .gbff directly and the RAST format is retired.
"""
