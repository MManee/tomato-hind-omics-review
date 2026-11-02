#!/usr/bin/env bash
# =============================================================================
# COI barcode analysis, step 1: download Cephalopholis COI data
#   (a) GenBank (NCBI nuccore) via Entrez Direct, frozen at a publication-date
#       cutoff so the query can be repeated exactly;
#   (b) BOLD Systems public portal API (genus Cephalopholis), fail-safe:
#       if BOLD cannot be reached, the pipeline continues with GenBank only.
#
# Run from Terminal (Python environment with Biopython active, the same one used for
# verify_saudi_records.py, because steps 1b and 2 need Python/Biopython):
#   cd tomato-hind-omics-review   # repository root
#   bash scripts/coi_01_download.sh
# =============================================================================
set -uo pipefail

PROJ="${PROJ:-$(pwd)}"   # repository root
CUTOFF="2026/10/08"                       # GenBank publication-date cutoff (fixed)
OUT="$PROJ/data/coi_20261008"
mkdir -p "$OUT"
LOG="$OUT/download_log.txt"

{
  echo "COI download log"
  echo "Run started      : $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
  echo "Entrez Direct    : $(esearch -version 2>/dev/null || echo 'not found')"
  echo "Python           : $(python3 --version 2>&1)"
  echo "GenBank cutoff   : PDAT <= $CUTOFF"
} > "$LOG"

# ---- (a) GenBank ------------------------------------------------------------
# All Cephalopholis mitochondrial records plus any record that names COI in the
# gene or title fields. RefSeq is excluded (NC_ records duplicate INSDC records,
# e.g. NC_022143 = KC593378). Records of 100-20,000 bp (mitogenomes included).
# COI is extracted from the annotation in step 2, so non-COI records are dropped there.
QUERY="\"Cephalopholis\"[Organism] AND (mitochondrion[filter] OR COI[Gene] OR COX1[Gene] OR CO1[Gene] OR \"cytochrome oxidase subunit\"[Title] OR \"cytochrome c oxidase subunit\"[Title]) AND 1900/01/01:${CUTOFF}[PDAT] AND 100:20000[SLEN] NOT srcdb_refseq[PROP]"
echo "GenBank query    : $QUERY" >> "$LOG"

GB="$OUT/genbank_cephalopholis_mito.gb"
if [ ! -s "$GB" ]; then
  N=$(esearch -db nuccore -query "$QUERY" | xtract -pattern ENTREZ_DIRECT -element Count)
  echo "GenBank hits     : $N" | tee -a "$LOG"
  esearch -db nuccore -query "$QUERY" | efetch -format gb > "$GB"
else
  echo "GenBank file exists, not re-downloaded: $GB" | tee -a "$LOG"
fi
echo "GenBank records  : $(grep -c '^LOCUS' "$GB")" | tee -a "$LOG"

# Named reference records, fetched explicitly so they are present even if the
# search above misses them: KC593378 (C. sonnerati mitogenome; COI reference
# frame for all positions) and the two C. oligosticta COI references.
REF="$OUT/genbank_named_refs.gb"
if [ ! -s "$REF" ]; then
  efetch -db nuccore -id KC593378,KM656840,MH707291 -format gb > "$REF"
fi
echo "Named refs       : $(grep '^VERSION' "$REF" | awk '{print $2}' | tr '\n' ' ')" | tee -a "$LOG"

# ---- (b) BOLD -----------------------------------------------------------------
python3 "$PROJ/scripts/coi_01b_bold.py" "$OUT" || echo "BOLD step returned an error (see bold_status.txt); continuing with GenBank only." | tee -a "$LOG"
cat "$OUT/bold_status.txt" >> "$LOG" 2>/dev/null

echo "Run finished     : $(date -u '+%Y-%m-%d %H:%M:%S UTC')" >> "$LOG"
echo; echo "Done. Log: $LOG"
