#!/usr/bin/env python3
"""
Verify the identity of the Saudi Arabian "Cephalopholis sonnerati" 16S/12S
records (KX298697-KX298699, KX298707-KX298710) by aligning each one against
every Cephalopholis 16S/12S record and every Cephalopholis mitogenome on GenBank.

If they are truly C. sonnerati, their best matches should be C. sonnerati
references at ~99-100% identity, clearly above any congener.

Usage (from the data folder, inside the project's virtual environment):
  python ../scripts/verify_saudi_records.py ncbi_inventory_20261001/verify_saudi
Requires: Entrez Direct on PATH (esearch/efetch) and Biopython.
"""
import os, sys, subprocess, collections, csv
from Bio import SeqIO
from Bio.Align import PairwiseAligner

OUT = sys.argv[1] if len(sys.argv) > 1 else "verify_saudi"
# Optional 2nd argument: comma-separated accession prefixes to EXCLUDE from the
# references (e.g. KX2986 = all records of the unpublished 2016 Saudi study)
EXCLUDE = [p for p in (sys.argv[2].split(",") if len(sys.argv) > 2 else []) if p]
TAG = ("_excl_" + "_".join(EXCLUDE)) if EXCLUDE else ""
os.makedirs(OUT, exist_ok=True)

QUERIES = ["KX298697.1", "KX298698.1", "KX298699.1",
           "KX298707.1", "KX298708.1", "KX298709.1", "KX298710.1"]

def fetch(query, path, by_ids=False):
    """Download sequences as GenBank flat file using Entrez Direct."""
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return
    if by_ids:
        cmd = f'efetch -db nuccore -id {",".join(query)} -format gb > "{path}"'
    else:
        cmd = f'esearch -db nuccore -query \'{query}\' | efetch -format gb > "{path}"'
    print("Downloading:", path)
    subprocess.run(cmd, shell=True, check=True)

q_gb   = os.path.join(OUT, "saudi_queries.gb")
r16_gb = os.path.join(OUT, "ref_cephalopholis_16S_12S.gb")
rmt_gb = os.path.join(OUT, "ref_cephalopholis_mitogenomes.gb")

fetch(QUERIES, q_gb, by_ids=True)
fetch('Cephalopholis[Organism] AND mitochondrion[Filter] AND (16S[Title] OR 12S[Title] OR '
      '"ribosomal RNA"[Title]) NOT "complete genome"[Title]', r16_gb)
fetch('Cephalopholis[Organism] AND mitochondrion[Filter] AND "complete genome"[Title]', rmt_gb)

queries = list(SeqIO.parse(q_gb, "genbank"))
refs = []
for f in (r16_gb, rmt_gb):
    for r in SeqIO.parse(f, "genbank"):
        if r.id in QUERIES:
            continue                      # don't match a query to itself
        if any(r.id.startswith(p) for p in EXCLUDE):
            continue                      # drop references from excluded studies
        country = ""
        for feat in r.features:
            if feat.type == "source":
                country = (feat.qualifiers.get("geo_loc_name") or
                           feat.qualifiers.get("country") or [""])[0]
        refs.append((r.id, r.annotations.get("organism", "?"), country, str(r.seq).upper()))
print(f"{len(queries)} queries, {len(refs)} reference sequences"
      + (f" (excluded prefixes: {', '.join(EXCLUDE)})" if EXCLUDE else ""))
sp_counts = collections.Counter(o for _, o, _, _ in refs)
for sp in ["Cephalopholis sonnerati", "Cephalopholis oligosticta", "Cephalopholis hemistiktos",
           "Cephalopholis miniata", "Cephalopholis urodeta", "Cephalopholis argus"]:
    srcs = sorted({c.split(":")[0] or "no locality" for i, o, c, _ in refs if o == sp})
    print(f"  references for {sp}: {sp_counts.get(sp, 0)}  [{'; '.join(srcs)}]")

aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score, aligner.mismatch_score = 2, -3
aligner.open_gap_score, aligner.extend_gap_score = -5, -2

def identity(q, s):
    """Percent identity and aligned length of the best local alignment."""
    aln = aligner.align(q, s)[0]
    c = aln.counts()
    alen = c.identities + c.mismatches + c.gaps
    return (100.0 * c.identities / alen if alen else 0.0), alen

rows = []
for q in queries:
    qs = str(q.seq).upper()
    gene = "16S" if "16S" in q.description else ("12S" if "12S" in q.description else "rRNA")
    hits = []
    for rid, org, ctry, s in refs:
        pid, alen = identity(qs, s)
        # coverage is judged against the SHORTER of the two sequences, so short
        # references (e.g. 340-380 bp 12S fragments) are still compared over
        # their shared region; require at least 300 aligned bp
        if alen >= max(300, 0.8 * min(len(qs), len(s))):
            hits.append((pid, alen, rid, org, ctry))
    hits.sort(reverse=True)
    best_by_sp = {}
    for pid, alen, rid, org, ctry in hits:
        if org not in best_by_sp:
            best_by_sp[org] = (pid, alen, rid, ctry)
    ranked = sorted(best_by_sp.items(), key=lambda kv: -kv[1][0])
    top_sp, top = ranked[0] if ranked else ("none", (0, 0, "", ""))
    son = best_by_sp.get("Cephalopholis sonnerati", (0, 0, "", ""))
    second = ranked[1] if len(ranked) > 1 else ("none", (0, 0, "", ""))
    rows.append({
        "query": q.id, "gene": gene, "length": len(qs),
        "best_species": top_sp, "best_identity": f"{top[0]:.2f}", "best_ref": top[2],
        "identity_to_C_sonnerati": f"{son[0]:.2f}", "C_sonnerati_ref": son[2],
        "second_species": second[0], "second_identity": f"{second[1][0]:.2f}",
        "n_species_compared": len(ranked),
    })
    with open(os.path.join(OUT, f"{q.id}_best_per_species{TAG}.tsv"), "w") as h:
        h.write("species\tidentity\taligned_len\tref_accession\tref_country\n")
        for sp, (pid, alen, rid, ctry) in ranked:
            h.write(f"{sp}\t{pid:.2f}\t{alen}\t{rid}\t{ctry}\n")

summary = os.path.join(OUT, f"saudi_records_verification{TAG}.tsv")
with open(summary, "w", newline="") as h:
    w = csv.DictWriter(h, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

print("\nquery        gene  len  best_species                     best%   to_sonnerati%  second_species              second%")
for r in rows:
    print(f"{r['query']:<12} {r['gene']:<4} {r['length']:<4} {r['best_species']:<32} "
          f"{r['best_identity']:>6}   {r['identity_to_C_sonnerati']:>6}         "
          f"{r['second_species']:<27} {r['second_identity']:>6}")

# Which Cephalopholis species have NO 16S/12S/mitogenome reference at all?
red_sea_congeners = ["Cephalopholis hemistiktos", "Cephalopholis oligosticta",
                     "Cephalopholis miniata", "Cephalopholis argus", "Cephalopholis sexmaculata"]
have = {o for _, o, _, _ in refs}
missing = [s for s in red_sea_congeners if s not in have]
print("\nRed Sea congeners with NO reference sequence for this comparison:",
      ", ".join(missing) if missing else "none")
print(f"\nPer-query details: {OUT}/<accession>_best_per_species{TAG}.tsv")
print(f"Summary table:     {summary}")
