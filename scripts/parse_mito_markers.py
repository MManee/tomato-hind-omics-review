#!/usr/bin/env python3
"""
Parse NCBI GenBank XML (efetch -format gbc) that may contain several
concatenated <INSDSet> documents, and summarise mitochondrial records
for Cephalopholis sonnerati by marker and locality.

Usage (from the data folder):
  python3 ../scripts/parse_mito_markers.py ncbi_inventory_20261001/followup/mito_records.xml
Writes, next to the input file:
  mito_markers_by_locality.tsv
  mito_markers_summary.txt
"""
import re, sys, os, csv, collections
import xml.etree.ElementTree as ET

src = sys.argv[1]
outdir = os.path.dirname(os.path.abspath(src))
raw = open(src, encoding="utf-8", errors="replace").read()

# Remove every XML declaration / DOCTYPE and wrap all chunks in one root
raw = re.sub(r"<\?xml[^>]*\?>", "", raw)
raw = re.sub(r"<!DOCTYPE[^>]*>", "", raw)
root = ET.fromstring("<ROOT>" + raw + "</ROOT>")

def marker(defn):
    d = defn.lower()
    if "complete genome" in d:
        return "mitogenome"
    for k, v in [("cytochrome oxidase subunit i", "COI"), ("cytochrome c oxidase subunit 1", "COI"),
                 ("cytochrome c oxidase subunit i", "COI"), ("(coi)", "COI"), ("cox1", "COI"),
                 ("cytochrome b", "CYTB"), ("cytb", "CYTB"), ("16s", "16S"), ("12s", "12S"),
                 ("control region", "D-loop"), ("d-loop", "D-loop"), ("nd2", "ND2"), ("nd4", "ND4")]:
        if k in d:
            return v
    return "other"

rows = []
for s in root.iter("INSDSeq"):
    q = {}
    for f in s.iter("INSDQualifier"):
        n = f.findtext("INSDQualifier_name")
        v = f.findtext("INSDQualifier_value") or ""
        if n and n not in q:
            q[n] = v
    ref = s.find(".//INSDReference/INSDReference_title")
    rows.append({
        "accession": s.findtext("INSDSeq_accession-version", ""),
        "length": s.findtext("INSDSeq_length", ""),
        "marker": marker(s.findtext("INSDSeq_definition", "")),
        "country": q.get("geo_loc_name", q.get("country", "")),
        "lat_lon": q.get("lat_lon", ""),
        "voucher_or_isolate": q.get("specimen_voucher", q.get("isolate", "")),
        "collection_date": q.get("collection_date", ""),
        "first_reference": (ref.text if ref is not None else "").replace("\t", " "),
        "definition": s.findtext("INSDSeq_definition", "").replace("\t", " "),
    })

tsv = os.path.join(outdir, "mito_markers_by_locality.tsv")
with open(tsv, "w", newline="") as h:
    w = csv.DictWriter(h, fieldnames=list(rows[0].keys()) if rows else ["accession"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)

ARAB = ["saudi", "egypt", "sudan", "eritrea", "yemen", "djibouti", "jordan", "israel", "oman",
        "emirates", "uae", "qatar", "bahrain", "kuwait", "iran", "red sea", "persian", "arabian",
        "aden", "somalia"]

lines = []
lines.append(f"Total mitochondrial records: {len(rows)}")
lines.append("\nRecords by marker:")
for k, v in collections.Counter(r["marker"] for r in rows).most_common():
    lines.append(f"  {k}: {v}")
lines.append("\nRecords by country (text before ':'; blank = not given):")
for k, v in collections.Counter((r["country"].split(":")[0].strip() or "(none)") for r in rows).most_common():
    lines.append(f"  {k}: {v}")
lines.append(f"\nRecords with lat_lon but no country: "
             f"{sum(1 for r in rows if r['lat_lon'] and not r['country'])}")
lines.append(f"Records with neither country nor lat_lon: "
             f"{sum(1 for r in rows if not r['lat_lon'] and not r['country'])}")
arab = [r for r in rows if any(w in (r["country"] + " " + r["lat_lon"]).lower() for w in ARAB)]
lines.append(f"\nRecords from Red Sea / Gulf of Aden / Arabian region: {len(arab)}")
for r in arab:
    lines.append(f"  {r['accession']}  {r['marker']}  {r['country']}  {r['lat_lon']}  {r['first_reference'][:80]}")
lines.append("\nCOI records by country:")
for k, v in collections.Counter((r["country"].split(":")[0].strip() or "(none)")
                                for r in rows if r["marker"] == "COI").most_common():
    lines.append(f"  {k}: {v}")

summary = "\n".join(lines)
open(os.path.join(outdir, "mito_markers_summary.txt"), "w").write(summary + "\n")
print(summary)
print(f"\nWrote {tsv}")
