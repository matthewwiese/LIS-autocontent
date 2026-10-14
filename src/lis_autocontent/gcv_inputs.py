"""Write the redis_loader inputs the Data Store doesn't publish.

Runs inside the redis_loader container, on the standard library alone:

    gcv_inputs.py chromosomes <genome .fai> <out.gff3> [--rename OLD NEW]
    gcv_inputs.py genes <gene_models_main.gff3.gz> <out.gff3> --rename OLD NEW

Sources are local paths or URLs, gunzipped when their names end in .gz. --rename
replaces the first ``.OLD.`` in every sequence name, so that two annotations of one
assembly load onto distinct chromosomes.
"""

import argparse
import gzip
import io
import sys
import urllib.parse
import urllib.request
from contextlib import contextmanager

USER_AGENT = "lis-autocontent gcv_inputs"
# Characters GFF3 reserves in column 9.
ATTRIBUTE_ESCAPES = str.maketrans({c: f"%{ord(c):02X}" for c in "%;=&,\t"})


@contextmanager
def open_text(source):
    """A local path or URL as lines of text."""
    if urllib.parse.urlparse(source).scheme:
        request = urllib.request.Request(source, headers={"User-Agent": USER_AGENT})
        handle = urllib.request.urlopen(request)  # pylint: disable=consider-using-with
    else:
        handle = open(source, "rb")  # pylint: disable=consider-using-with
    with handle:
        raw = gzip.GzipFile(fileobj=handle) if source.endswith(".gz") else handle
        yield io.TextIOWrapper(raw, encoding="utf-8")


def renamer(rename):
    """A seqid transform for an (old, new) token pair, or the identity for None."""
    if rename is None:
        return lambda seqid: seqid
    old, new = (f".{token}." for token in rename)

    def rename_seqid(seqid):
        if old not in seqid:
            raise ValueError(f"sequence {seqid} has no {old} to rename")
        return seqid.replace(old, new, 1)

    return rename_seqid


def chromosome_rows(lines, rename=None):
    """GFF3 lines with one chromosome per .fai sequence, from 1 to its length.

    redis_loader stores chromosome and supercontig features alike, so all are chromosomes.
    """
    rename_seqid = renamer(rename)
    yield "##gff-version 3\n"
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2:
            continue
        seqid = rename_seqid(fields[0])
        name = seqid.translate(ATTRIBUTE_ESCAPES)
        row = [seqid, ".", "chromosome", "1", fields[1], ".", ".", "."]
        yield "\t".join(row + [f"ID={name};Name={name}"]) + "\n"


def gene_rows(lines, rename=None):
    """The GFF3's gene lines, up to any ##FASTA, with their sequences renamed."""
    rename_seqid = renamer(rename)
    yield "##gff-version 3\n"
    for line in lines:
        if line.startswith("##FASTA"):
            break
        fields = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(fields) != 9 or fields[2] != "gene":
            continue
        fields[0] = rename_seqid(fields[0])
        yield "\t".join(fields) + "\n"


def main(argv=None):
    """Command-line entry; exits non-zero, writing nothing, on a sequence it can't rename."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("kind", choices=("chromosomes", "genes"))
    parser.add_argument("source")
    parser.add_argument("out")
    parser.add_argument("--rename", nargs=2, metavar=("OLD", "NEW"))
    args = parser.parse_args(argv)
    rows = chromosome_rows if args.kind == "chromosomes" else gene_rows
    try:
        with open_text(args.source) as lines:
            text = "".join(rows(lines, args.rename))
    except ValueError as err:
        sys.exit(f"{args.source}: {err}")
    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write(text)


if __name__ == "__main__":
    main()
