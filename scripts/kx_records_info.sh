#!/usr/bin/env bash
# =============================================================================
# Retrieve the full GenBank records for the Saudi Arabian Cephalopholis records
# KX298691 and KX298695-KX298699 (16S), together with the 12S records from the
# same submission (KX298707-KX298710) and the neighbouring accessions
# (KX298680-KX298720), so that the whole submission can be identified.
#
# Output (data/kx2986_check/):
#   kx_records.gb          full GenBank flat files (the primary record)
#   kx_records_summary.tsv one row per record: organism, gene, length, country,
#                          voucher/isolate, collector, dates, authors, title,
#                          journal/submitter line, BioProject/BioSample, PopSet
#   kx_links.tsv           PubMed and PopSet links found by elink
#   kx_log.txt             tool versions and run date
#
# Run (conda env with Biopython active):
#   (Python environment with Biopython active)
#   cd tomato-hind-omics-review   # repository root
#   bash scripts/kx_records_info.sh
# =============================================================================
set -uo pipefail
PROJ="${PROJ:-$(pwd)}"   # repository root
OUT="$PROJ/data/kx2986_check"
mkdir -p "$OUT"
LOG="$OUT/kx_log.txt"
{
  echo "Run: $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
  echo "Entrez Direct: $(esearch -version 2>/dev/null)"
  echo "Python: $(python3 --version 2>&1)"
} > "$LOG"

# accessions KX298680-KX298720 (includes the requested ones)
IDS=$(seq 298680 298720 | sed 's/^/KX/' | paste -sd, -)
echo "Requested range: KX298680-KX298720" >> "$LOG"

efetch -db nuccore -id "$IDS" -format gb > "$OUT/kx_records.gb" 2>> "$LOG"
echo "Records retrieved: $(grep -c '^LOCUS' "$OUT/kx_records.gb")" | tee -a "$LOG"

# links to PubMed and PopSet for each record
echo -e "accession\tpubmed_ids\tpopset_ids" > "$OUT/kx_links.tsv"
for acc in $(grep '^VERSION' "$OUT/kx_records.gb" | awk '{print $2}'); do
  pm=$(elink -db nuccore -id "$acc" -target pubmed 2>/dev/null | efetch -format uid 2>/dev/null | paste -sd, -)
  ps=$(elink -db nuccore -id "$acc" -target popset 2>/dev/null | efetch -format uid 2>/dev/null | paste -sd, -)
  echo -e "$acc\t${pm:-none}\t${ps:-none}" >> "$OUT/kx_links.tsv"
  sleep 0.4
done

python3 "$PROJ/scripts/kx_records_parse.py" "$OUT"
echo; echo "Done. Outputs in $OUT"
