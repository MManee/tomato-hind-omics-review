#!/usr/bin/env python3
"""
COI barcode analysis, step 2: extract COI sequences and metadata from the
GenBank flat files downloaded in step 1.

COI is taken from the record's own annotation (CDS, else gene feature) whose
/gene is COI/COX1/CO1/COXI/MT-CO1 or whose /product is "cytochrome (c) oxidase
subunit I/1". Records without such a feature are used whole only if they are
<= 2,500 bp and the definition line names COI; this choice is recorded per
sequence in the 'coi_source' column so it can be checked.

Usage (Python environment with Biopython active, from the repository root):
  python3 scripts/coi_02_parse_genbank.py data/coi_20261008
Output: data/coi_20261008/coi_genbank.tsv and parse_log.txt
Requires Biopython.
"""
import sys, os, re, csv, collections
from Bio import SeqIO
import Bio

D = sys.argv[1] if len(sys.argv) > 1 else "."
FILES = ["genbank_cephalopholis_mito.gb", "genbank_named_refs.gb"]

GENE_RE = re.compile(r"^(COI|COX1|CO1|COXI|MT-CO1|CO-?I)$", re.I)
PROD_RE = re.compile(r"cytochrome\s*c?\s*oxidase\s+(subunit\s+)?(I|1)(?![I\d])", re.I)
DEF_RE  = re.compile(r"(cytochrome\s*c?\s*oxidase\s+subunit\s+(I|1)(?![I\d])|\bCOI\b|\bCOX1\b|\bCO1\b)", re.I)

def q(f, key):
    return f.qualifiers.get(key, [""])[0]

def is_coi(f):
    return GENE_RE.match(q(f, "gene").strip()) or PROD_RE.search(q(f, "product"))

rows, seen = [], set()
counts = collections.Counter()
for fn in FILES:
    path = os.path.join(D, fn)
    if not os.path.exists(path):
        print("missing:", path); continue
    for rec in SeqIO.parse(path, "genbank"):
        acc = rec.id
        if acc in seen:
            counts["duplicate_record_skipped"] += 1; continue
        seen.add(acc); counts["records_read"] += 1
        src = next((f for f in rec.features if f.type == "source"), None)
        sq = lambda k: q(src, k) if src is not None else ""
        cds = [f for f in rec.features if f.type == "CDS" and is_coi(f)]
        gene = [f for f in rec.features if f.type == "gene" and is_coi(f)]
        if cds:
            seq, how = str(cds[0].extract(rec.seq)), "CDS_annotation"
        elif gene:
            seq, how = str(gene[0].extract(rec.seq)), "gene_annotation"
        elif len(rec.seq) <= 2500 and DEF_RE.search(rec.description):
            seq, how = str(rec.seq), "whole_record_definition"
        else:
            counts["no_COI"] += 1; continue
        counts[how] += 1
        loc = sq("geo_loc_name") or sq("country")
        dbx = ";".join(x for x in (src.qualifiers.get("db_xref", []) if src is not None else [])
                       if x.upper().startswith("BOLD"))
        ref_title = rec.annotations.get("references", [None])[0]
        rows.append({
            "accession": acc,
            "organism": rec.annotations.get("organism", ""),
            "country": loc.split(":")[0].strip(),
            "locality": loc,
            "lat_lon": sq("lat_lon"),
            "specimen_voucher": sq("specimen_voucher"),
            "isolate": sq("isolate"),
            "collection_date": sq("collection_date"),
            "bold_xref": dbx,
            "record_length": len(rec.seq),
            "coi_source": how,
            "first_reference": ref_title.title if ref_title is not None else "",
            "definition": rec.description,
            "seq": seq.upper(),
        })

out = os.path.join(D, "coi_genbank.tsv")
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

with open(os.path.join(D, "parse_log.txt"), "w") as f:
    f.write(f"Biopython {Bio.__version__}\n")
    for k, v in sorted(counts.items()):
        f.write(f"{k}\t{v}\n")
    f.write(f"COI sequences written\t{len(rows)}\n")
    sp = collections.Counter(r["organism"] for r in rows)
    f.write("\nCOI sequences per organism name:\n")
    for k, v in sp.most_common():
        f.write(f"{k}\t{v}\n")
print(open(os.path.join(D, "parse_log.txt")).read())
print("Written:", out)
