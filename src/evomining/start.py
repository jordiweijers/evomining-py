#!/usr/bin/env python3
"""
evomining.start
==============

Forward + reverse BLAST between the Central (enzyme) DB and the genome DB.

  Forward:  Central DB  vs  GENOMES.fasta   -> thereand.blast
  Reverse:  GENOMES.fasta  vs  Central DB   -> backagain.blast

The genome DB is the single GENOMES.fasta written by `evomining generate-genome-db`
(protein records keyed on composite <genome_stem>__<locus_tag> IDs). There is no
longer a directory of per-genome .faa to concatenate, and no RAST ids file: BLAST
operates on the protein FASTA directly and does not care about the header format.

Reverse BLAST is chunked by query. With millions of genome proteins as the query,
blastp keeps the entire query resident and then fails to allocate thread stacks
("CThread::Run() -- error creating thread"), regardless of --num_threads. So the
reverse query is split into blocks of REVERSE_CHUNK_SIZE proteins, each BLASTed
against the (tiny) Central DB and concatenated. outfmt-6 rows are independent per
query, so the concatenation is identical to a single run; `analyze` does not
depend on row order. The forward BLAST queries only the 455-seq Central DB and
needs no chunking.

Each step skips if its output already exists, so re-runs are safe; the reverse
step additionally skips already-finished chunks. `start` does not touch MIBiG or
antiSMASH; the output directory is labelled {central}_{genomes} and discovered by
`analyze` via globbing.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

from .external import Tools


# Reverse-BLAST query block size, in proteins. 100k keeps the resident query
# small enough that thread-stack allocation succeeds while still amortising
# blastp startup over many sequences. Override with --chunk-size.
REVERSE_CHUNK_SIZE = 100_000


def _iter_fasta_chunks(path, chunk_size):
    """Yield successive blocks of `chunk_size` FASTA records as text.

    Streams the file; each yielded string is a valid multi-record FASTA. Chunk
    boundaries are deterministic (fixed size, file order), so a given chunk index
    always covers the same records across runs -- which is what makes per-chunk
    resume correct.
    """
    buf = []
    count = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if count == chunk_size:
                    yield "".join(buf)
                    buf = []
                    count = 0
                count += 1
            buf.append(line)
    if buf:
        yield "".join(buf)


def _reverse_blast_chunked(tools, genome_fasta, central_blastdb, reverse_blast,
                           chunk_dir, threads, evalue, chunk_size):
    """Reverse BLAST (genome proteins vs Central DB), one query block at a time.

    Per-chunk outputs are written under `chunk_dir` and concatenated into
    `reverse_blast`. A finished chunk is a renamed (atomic) .blast file, so an
    interrupted run resumes by re-BLASTing only the chunks that never completed.
    On success the chunk directory is removed.
    """
    chunk_dir.mkdir(parents=True, exist_ok=True)
    out_paths = []

    for idx, records in enumerate(_iter_fasta_chunks(genome_fasta, chunk_size), start=1):
        out_path = chunk_dir / f"backagain_{idx:05d}.blast"
        out_paths.append(out_path)

        if out_path.exists():
            print(f"    chunk {idx:>4}: done, skipping")
            continue

        query_path = chunk_dir / f"query_{idx:05d}.faa"
        query_path.write_text(records)
        tmp_out = chunk_dir / f"backagain_{idx:05d}.blast.part"

        print(f"    chunk {idx:>4}: blastp...")
        tools.run("blastp",
                  ["-query", query_path, "-db", central_blastdb, "-out", tmp_out,
                   "-outfmt", "6", "-evalue", evalue,
                   "-num_threads", threads, "-max_target_seqs", "10000"],
                  check=True)

        os.replace(tmp_out, out_path)   # atomic: only a completed chunk is named .blast
        query_path.unlink(missing_ok=True)

    # Concatenate finished chunks (in index order) into the final reverse file.
    with open(reverse_blast, "w") as out:
        for p in out_paths:
            with open(p) as fh:
                shutil.copyfileobj(fh, out)

    shutil.rmtree(chunk_dir, ignore_errors=True)
    return reverse_blast


def run(args, tools=None):
    """Forward + reverse BLAST between the Central DB and the genome DB."""
    tools = tools or Tools()
    tools.require("makeblastdb", "blastp")

    genome_fasta = Path(args.genomes).resolve()
    central_db = Path(args.central_db).resolve()
    outdir = Path(args.output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if not genome_fasta.is_file():
        raise SystemExit(
            f"ERROR: genome FASTA not found: {genome_fasta}\n"
            f"Run `evomining generate-genome-db` first, or pass --genomes."
        )

    chunk_size = int(getattr(args, "chunk_size", None) or REVERSE_CHUNK_SIZE)

    central_name = central_db.name
    # Token that labels the run. Taken from the genome DB directory name
    # (conventionally 'evomining_db'), now that the DB is a single FASTA.
    genomes_name = genome_fasta.parent.name

    result_dir = outdir / f"{central_name}_{genomes_name}"
    blast_dir = result_dir / "blast"
    blast_dir.mkdir(parents=True, exist_ok=True)

    n_proteins = sum(1 for line in open(genome_fasta) if line.startswith(">"))

    print("=" * 60)
    print("  EvoMining BLAST Pipeline")
    print("=" * 60)
    print(f"  Genome FASTA:      {genome_fasta}  ({n_proteins} proteins)")
    print(f"  Central DB:        {central_db}")
    print(f"  Output:            {result_dir}")
    print(f"  Threads:           {args.threads}")
    print(f"  E-value:           {args.evalue}")
    print(f"  Reverse chunk:     {chunk_size} proteins/block")
    print()

    # ---- Step 1: Forward BLAST (Central DB vs Genome DB) ----
    forward_blast = blast_dir / "thereand.blast"
    genome_db = outdir / "genome_blastdb"

    if forward_blast.exists() and forward_blast.stat().st_size > 0:
        n_hits = sum(1 for _ in open(forward_blast))
        print(f"Step 1: Forward BLAST already exists ({n_hits} hits)")
    else:
        print("Step 1: Forward BLAST (Central DB vs Genome DB)")
        print("  Building genome BLAST database...")
        tools.run("makeblastdb",
                  ["-in", genome_fasta, "-dbtype", "prot", "-out", genome_db],
                  check=True)

        print("  Running blastp...")
        tools.run("blastp",
                  ["-query", central_db, "-db", genome_db, "-out", forward_blast,
                   "-outfmt", "6", "-evalue", args.evalue,
                   "-num_threads", args.threads, "-max_target_seqs", "10000"],
                  check=True)

        n_hits = sum(1 for _ in open(forward_blast))
        print(f"  Forward BLAST: {n_hits} hits")

    # ---- Step 2: Reverse BLAST (Genome DB vs Central DB), chunked by query ----
    reverse_blast = blast_dir / "backagain.blast"
    central_blastdb = outdir / "central_blastdb"

    if reverse_blast.exists() and reverse_blast.stat().st_size > 0:
        n_hits = sum(1 for _ in open(reverse_blast))
        print(f"\nStep 2: Reverse BLAST already exists ({n_hits} hits)")
    else:
        print("\nStep 2: Reverse BLAST (Genome DB vs Central DB)")
        print("  Building Central DB BLAST database...")
        tools.run("makeblastdb",
                  ["-in", central_db, "-dbtype", "prot", "-out", central_blastdb],
                  check=True)

        n_chunks = (n_proteins + chunk_size - 1) // chunk_size
        print(f"  Running blastp in {n_chunks} query chunk(s) of up to {chunk_size}...")
        _reverse_blast_chunked(
            tools, genome_fasta, central_blastdb, reverse_blast,
            chunk_dir=blast_dir / "reverse_chunks",
            threads=args.threads, evalue=args.evalue, chunk_size=chunk_size,
        )

        n_hits = sum(1 for _ in open(reverse_blast))
        print(f"  Reverse BLAST: {n_hits} hits")

    print(f"\n{'='*60}")
    print(f"  BLAST complete. Next: evomining analyze")
    print(f"  Forward BLAST:    {forward_blast}")
    print(f"  Reverse BLAST:    {reverse_blast}")
    print(f"{'='*60}")
    return blast_dir
