#!/usr/bin/env python3
"""
COI barcode analysis, step 1b: download public BOLD records for the genus
Cephalopholis through the BOLD portal API (https://portal.boldsystems.org/api),
following the three documented calls: preprocessor -> query (query_id) ->
documents/<query_id>/download?format=json.

Fail-safe: any error is written to bold_status.txt and the script exits 0, so
the rest of the pipeline continues with GenBank data only.

Output (in OUT):
  bold_raw.json        raw download (BCDM JSON, kept unchanged for the record)
  bold_records.tsv     one row per record that has a COI-5P sequence
  bold_fields.txt      field names present in the download (for checking)
  bold_status.txt      OK / FAILED + reason, date, counts
Standard library only.
"""
import sys, os, json, csv, time, urllib.request, urllib.parse, datetime

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
BASE = "https://portal.boldsystems.org/api"
TERMS = ["tax:genus:Cephalopholis", "tax:Cephalopholis"]   # second is a fallback
status_path = os.path.join(OUT, "bold_status.txt")

def get(url, timeout=300):
    req = urllib.request.Request(url, headers={"User-Agent": "tomato-hind-review/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()

def status(msg):
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    with open(status_path, "w") as f:
        f.write(f"BOLD download ({stamp}): {msg}\n")
    print("BOLD:", msg)

def first(rec, *keys):
    """Return the first non-empty value among candidate field names."""
    for k in keys:
        v = rec.get(k)
        if v not in (None, "", [], {}):
            return v
    return ""

def parse_records(raw):
    """BCDM JSON may arrive as one JSON array/object or as JSON lines."""
    txt = raw.decode("utf-8", errors="replace").strip()
    try:
        obj = json.loads(txt)
        if isinstance(obj, list):
            return obj
        if isinstance(obj, dict):
            for k in ("records", "data", "results"):
                if isinstance(obj.get(k), list):
                    return obj[k]
            return [obj]
    except json.JSONDecodeError:
        pass
    return [json.loads(l) for l in txt.splitlines() if l.strip().startswith("{")]

try:
    raw_path = os.path.join(OUT, "bold_raw.json")
    if not (os.path.exists(raw_path) and os.path.getsize(raw_path) > 0):
        qid, used = None, None
        for term in TERMS:
            try:
                pre = get(f"{BASE}/query/preprocessor?query={urllib.parse.quote(term, safe=':;/')}")
                with open(os.path.join(OUT, "bold_preprocessor.json"), "wb") as f:
                    f.write(pre)
                q = json.loads(get(f"{BASE}/query?query={urllib.parse.quote(term, safe=':;/')}&extent=full"))
                qid = q.get("query_id") or next((v for k, v in q.items() if "query" in k and isinstance(v, str)), None)
                if qid:
                    used = term
                    break
            except Exception as e:
                print(f"BOLD query '{term}' failed: {e}")
        if not qid:
            raise RuntimeError("no query_id returned for any query term")
        time.sleep(1)
        raw = get(f"{BASE}/documents/{urllib.parse.quote(qid)}/download?format=json", timeout=900)
        with open(raw_path, "wb") as f:
            f.write(raw)
        with open(os.path.join(OUT, "bold_query_used.txt"), "w") as f:
            f.write(f"query term: {used}\nquery_id: {qid}\n")
    with open(raw_path, "rb") as f:
        recs = parse_records(f.read())
    if not recs:
        raise RuntimeError("download contained no records")

    fields = sorted({k for r in recs for k in r})
    with open(os.path.join(OUT, "bold_fields.txt"), "w") as f:
        f.write("\n".join(fields) + "\n")

    cols = ["processid", "sampleid", "species", "identification", "identified_by",
            "bin_uri", "marker_code", "country", "province_state", "lat", "lon",
            "insdc_acs", "collection_date", "inst", "nuc_len", "nuc"]
    n_coi = 0
    with open(os.path.join(OUT, "bold_records.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in recs:
            marker = str(first(r, "marker_code", "markercode"))
            seq = str(first(r, "nuc", "nucleotides", "sequence")).replace("-", "").upper()
            if "COI" not in marker.upper() or len(seq) < 100:
                continue
            coord = first(r, "coord", "coordinates")
            lat = lon = ""
            if isinstance(coord, (list, tuple)) and len(coord) == 2:
                lat, lon = coord
            elif isinstance(coord, str) and "," in coord:
                lat, lon = [c.strip(" []") for c in coord.split(",")[:2]]
            w.writerow([
                first(r, "processid"), first(r, "sampleid"),
                first(r, "species", "species_name"),
                first(r, "identification"), first(r, "identified_by"),
                first(r, "bin_uri"), marker,
                first(r, "country/ocean", "country_ocean", "country"),
                first(r, "province/state", "province_state"),
                lat, lon,
                first(r, "insdc_acs", "genbank_accession"),
                first(r, "collection_date_start", "collection_date"),
                first(r, "inst", "institution_storing"),
                len(seq), seq])
            n_coi += 1
    status(f"OK; {len(recs)} records downloaded, {n_coi} with a COI sequence of >= 100 bp")
except Exception as e:
    status(f"FAILED ({type(e).__name__}: {e}); the analysis will use GenBank only")
sys.exit(0)
