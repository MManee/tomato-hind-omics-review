#!/usr/bin/env python3
"""Summarise the GenBank records retrieved by kx_records_info.sh (one row per record).
Usage: python3 kx_records_parse.py data/kx2986_check"""
import sys, os, csv
from Bio import SeqIO

D = sys.argv[1]
rows = []
for r in SeqIO.parse(os.path.join(D, "kx_records.gb"), "genbank"):
    src = next((f for f in r.features if f.type == "source"), None)
    q = lambda k: "; ".join(src.qualifiers.get(k, [])) if src is not None else ""
    genes = sorted({(f.qualifiers.get("gene", [""])[0] or f.qualifiers.get("product", [""])[0])
                    for f in r.features if f.type in ("rRNA", "gene", "CDS")} - {""})
    refs = r.annotations.get("references", [])
    xrefs = "; ".join(r.dbxrefs)
    rows.append({
        "accession": r.id,
        "organism": r.annotations.get("organism", ""),
        "definition": r.description,
        "length_bp": len(r.seq),
        "genes": "; ".join(genes),
        "country": q("geo_loc_name") or q("country"),
        "lat_lon": q("lat_lon"),
        "specimen_voucher": q("specimen_voucher"),
        "isolate": q("isolate"),
        "collection_date": q("collection_date"),
        "collected_by": q("collected_by"),
        "identified_by": q("identified_by"),
        "tissue_type": q("tissue_type"),
        "note": q("note"),
        "record_date": r.annotations.get("date", ""),
        "dbxrefs": xrefs,
        "ref_authors": " || ".join(x.authors for x in refs),
        "ref_titles": " || ".join(x.title for x in refs),
        "ref_journals": " || ".join(x.journal for x in refs),  # includes 'Submitted (date) address'
        "ref_pubmed": " || ".join(x.pubmed_id for x in refs if x.pubmed_id),
    })
out = os.path.join(D, "kx_records_summary.tsv")
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(f"{len(rows)} records summarised -> {out}")
for r in rows:
    print(f"{r['accession']}\t{r['organism']}\t{r['genes']}\t{r['length_bp']} bp\t{r['country']}\tvoucher={r['specimen_voucher'] or '-'}")
