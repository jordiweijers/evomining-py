from pathlib import Path

BLAST_DIR = Path("data/blast")
METADATA_DIR = Path("data/metadata")

referenced = set()
for name, col in [("thereand.blast", 1), ("backagain.blast", 0)]:
    with open(BLAST_DIR / name) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 12:
                value = fields[col]
                if "__" in value:
                    referenced.add(value)

referenced_stems = {value.split("__", 1)[0] for value in referenced}

print(f"{len(referenced)} referenced proteins found "
      f"in {len(referenced_stems)} genomes")

# genome_names.tsv is keyed on the bare genome stem (one row per genome),
# genome_functions.tsv on the composite <stem>__<locus_tag> protein id.
for name, wanted in (("genome_names.tsv", referenced_stems),
                     ("genome_functions.tsv", referenced)):
    path = METADATA_DIR / name
    with open(path) as fh:
        lines = fh.readlines()
    header, rows = lines[0], lines[1:]
    kept = [header] + [r for r in rows if r.split("\t", 1)[0] in wanted]
    with open(path, "w") as fh:
        fh.writelines(kept)
    print(f"{name}: {len(rows)} -> {len(kept) - 1} rows")
