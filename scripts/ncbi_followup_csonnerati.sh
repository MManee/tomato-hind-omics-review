#!/usr/bin/env bash
# =============================================================================
# Follow-up to ncbi_inventory_csonnerati.sh
#  1. Robust COI / mitochondrial-marker table with country + lat_lon (fixed parser)
#  2. Robust BioSample table (locality, tissue, sex, stage, isolate, collector)
#  3. Details of BioProjects/assemblies that need checking
#  4. SRA runs with sample-level metadata (to catch mislabelled library strategies)
#
# Usage (from the data folder):
#   bash ../scripts/ncbi_followup_csonnerati.sh ncbi_inventory_20261001
# =============================================================================
set -uo pipefail
export PATH=${HOME}/bin:${HOME}/edirect:${PATH}
SP="Cephalopholis sonnerati"
OUT="${1:?give the inventory folder, e.g. ncbi_inventory_20261001}"
mkdir -p "$OUT/followup"
F="$OUT/followup"

# ---------------------------------------------------------------------------
# 1. Mitochondrial markers (COI, CYTB, 16S, 12S, control region) -> GenBank XML
#    then parse with Python so missing fields stay empty instead of shifting
# ---------------------------------------------------------------------------
esearch -db nuccore -query "\"$SP\"[Organism] AND mitochondrion[Filter]" \
 | efetch -format gbc > "$F/mito_records.xml"

python3 - "$F/mito_records.xml" > "$F/mito_markers_by_locality.tsv" <<'PY'
import sys, xml.etree.ElementTree as ET
root = ET.parse(sys.argv[1]).getroot()
print("accession\tlength\tmarker\tcountry\tlat_lon\tisolate_or_voucher\tcollection_date\tfirst_reference")
def marker(defn):
    d = defn.lower()
    if "complete genome" in d: return "mitogenome"
    for k, v in [("cytochrome oxidase subunit i", "COI"), ("cytochrome c oxidase subunit 1", "COI"),
                 ("coi", "COI"), ("cox1", "COI"), ("cytochrome b", "CYTB"), ("cytb", "CYTB"),
                 ("16s", "16S"), ("12s", "12S"), ("control region", "D-loop"), ("d-loop", "D-loop")]:
        if k in d: return v
    return "other"
for s in root.iter("INSDSeq"):
    q = {}
    for f in s.iter("INSDQualifier"):
        n = f.findtext("INSDQualifier_name"); v = f.findtext("INSDQualifier_value") or ""
        if n and n not in q: q[n] = v
    ref = s.find(".//INSDReference/INSDReference_title")
    print("\t".join([
        s.findtext("INSDSeq_accession-version", ""),
        s.findtext("INSDSeq_length", ""),
        marker(s.findtext("INSDSeq_definition", "")),
        q.get("country", q.get("geo_loc_name", "")),
        q.get("lat_lon", ""),
        q.get("specimen_voucher", q.get("isolate", "")),
        q.get("collection_date", ""),
        (ref.text if ref is not None else "").replace("\t", " "),
    ]))
PY
echo "[1] mitochondrial marker table done"

# Quick look: records by marker and by country
python3 - "$F/mito_markers_by_locality.tsv" > "$F/mito_markers_summary.txt" <<'PY'
import sys, csv, collections
r = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
print("Records by marker:", dict(collections.Counter(x["marker"] for x in r)))
print("\nRecords by country (blank = not given):")
for k, v in collections.Counter(x["country"].split(":")[0] or "(none)" for x in r).most_common():
    print(f"  {k}: {v}")
print("\nRecords with lat_lon but no country:", sum(1 for x in r if x["lat_lon"] and not x["country"]))
red = [x for x in r if any(w in (x["country"] + x["lat_lon"]).lower()
       for w in ["saudi", "egypt", "sudan", "eritrea", "yemen", "djibouti", "jordan", "israel",
                 "oman", "emirates", "qatar", "bahrain", "kuwait", "iran", "red sea", "persian", "arabian"])]
print("\nRecords from Red Sea / Arabian region:", len(red))
for x in red: print("  ", x["accession"], x["country"], x["lat_lon"])
PY
echo "[1b] summary done"

# ---------------------------------------------------------------------------
# 2. BioSamples, parsed robustly
# ---------------------------------------------------------------------------
esearch -db biosample -query "\"$SP\"[Organism]" | efetch -format xml > "$F/biosamples.xml"
python3 - "$F/biosamples.xml" > "$F/biosamples_clean.tsv" <<'PY'
import sys, xml.etree.ElementTree as ET
root = ET.parse(sys.argv[1]).getroot()
keys = ["geo_loc_name", "lat_lon", "tissue", "sex", "dev_stage", "isolate", "collected_by",
        "collection_date", "specimen_voucher"]
print("accession\tsubmission_date\ttitle\t" + "\t".join(keys))
for b in root.iter("BioSample"):
    a = {x.get("attribute_name"): (x.text or "") for x in b.iter("Attribute")}
    title = b.findtext(".//Description/Title", "")
    print("\t".join([b.get("accession", ""), b.get("submission_date", "")[:10], title] +
                    [a.get(k, "").replace("\t", " ") for k in keys]))
PY
echo "[2] BioSamples done"

# ---------------------------------------------------------------------------
# 3. Projects / assemblies to check
# ---------------------------------------------------------------------------
for p in PRJNA1098793 PRJNA720393 PRJNA1160846 PRJEB123829; do
  echo "=== $p ===" >> "$F/projects_to_check.txt"
  esearch -db bioproject -query "$p" | esummary \
   | xtract -pattern DocumentSummary -element Project_Acc Organization Project_Title Project_Description \
   >> "$F/projects_to_check.txt"
  sleep 0.4
done
# Sample behind the new OceanOmics draft genome
datasets summary genome accession GCA_988286015.1 > "$F/GCA_988286015_summary.json"
esearch -db biosample -query "SAMN49597358" | efetch -format xml > "$F/SAMN49597358.xml"
echo "[3] projects done"

# ---------------------------------------------------------------------------
# 4. SRA runs with experiment titles (to flag strategy/source mismatches)
# ---------------------------------------------------------------------------
esearch -db sra -query "\"$SP\"[Organism]" | efetch -format runinfo > "$F/sra_runinfo.csv"
python3 - "$F/sra_runinfo.csv" > "$F/sra_by_project.tsv" <<'PY'
import sys, csv, collections
r = [x for x in csv.DictReader(open(sys.argv[1])) if x.get("Run") and x["Run"] != "Run"]
print("BioProject\truns\tGb\tstrategies\tsources\tplatforms\tmismatch_flag")
g = collections.defaultdict(list)
for x in r: g[x["BioProject"]].append(x)
for bp, xs in g.items():
    st = collections.Counter(x["LibraryStrategy"] for x in xs)
    so = collections.Counter(x["LibrarySource"] for x in xs)
    pl = collections.Counter(x["Platform"] for x in xs)
    gb = sum(float(x["bases"] or 0) for x in xs) / 1e9
    flag = "WGS_but_TRANSCRIPTOMIC" if any(x["LibraryStrategy"] == "WGS" and x["LibrarySource"] == "TRANSCRIPTOMIC" for x in xs) else ""
    print(f"{bp}\t{len(xs)}\t{gb:.1f}\t{dict(st)}\t{dict(so)}\t{dict(pl)}\t{flag}")
PY
echo "[4] SRA by project done"

echo
cat "$F/mito_markers_summary.txt"
echo
column -t -s $'\t' "$F/sra_by_project.tsv"
echo
cat "$F/projects_to_check.txt"
