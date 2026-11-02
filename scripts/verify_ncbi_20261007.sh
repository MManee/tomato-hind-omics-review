#!/usr/bin/env bash
# =============================================================================
# verify_ncbi_20261007.sh
# Runs the NCBI record checks (V1-V6) used in the C. sonnerati review in one go.
# Requires: Entrez Direct (esearch/efetch/esummary/xtract) and NCBI Datasets
# CLI (datasets).
#
# Usage (from the project folder):
#   cd tomato-hind-omics-review   # repository root
#   bash scripts/verify_ncbi_20261007.sh
# Outputs are written to data/verification_20261007/.
# =============================================================================
set -u
OUT="data/verification_20261007"
mkdir -p "$OUT"
LOG="$OUT/00_run_log.txt"
echo "Run started: $(date)" > "$LOG"
echo "edirect: $(esearch -version 2>&1 | head -1)" >> "$LOG"
echo "datasets: $(datasets --version 2>&1 | head -1)" >> "$LOG"

step () { echo; echo "=== $1 ==="; echo "=== $1 ===" >> "$LOG"; }

# -----------------------------------------------------------------------------
# V1. PRJNA1433391 (gut microbiome, Zhu et al. 2026) -- is it public?
# -----------------------------------------------------------------------------
step "V1 BioProject PRJNA1433391"
esearch -db bioproject -query PRJNA1433391 \
  | esummary \
  | xtract -pattern DocumentSummary -element Project_Acc Project_Title \
      Organism_Name Registration_Date Project_Description \
  > "$OUT/V1_PRJNA1433391_bioproject.tsv" 2>>"$LOG"
esearch -db sra -query PRJNA1433391 | efetch -format runinfo \
  > "$OUT/V1_PRJNA1433391_sra_runinfo.csv" 2>>"$LOG"
echo "BioProject lines: $(wc -l < "$OUT/V1_PRJNA1433391_bioproject.tsv")" | tee -a "$LOG"
echo "SRA run lines (incl. header): $(wc -l < "$OUT/V1_PRJNA1433391_sra_runinfo.csv")" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# V2. BioSample SAMN49597358 (OG354 draft genome) -- locality and voucher
# -----------------------------------------------------------------------------
step "V2 BioSample SAMN49597358"
efetch -db biosample -id SAMN49597358 -format full \
  > "$OUT/V2_SAMN49597358_full.txt" 2>>"$LOG"
grep -iE "geo_loc|lat_lon|collection|voucher|collected|isolate|tissue|sex|organism|identified" \
  "$OUT/V2_SAMN49597358_full.txt" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# V3. PRJNA720393 -- exact title (Table 4 says "Genome skimming ... (NOAA)")
# -----------------------------------------------------------------------------
step "V3 BioProject PRJNA720393"
esearch -db bioproject -query PRJNA720393 \
  | esummary \
  | xtract -pattern DocumentSummary -element Project_Acc Project_Title \
      Organism_Name Registration_Date Project_Description \
  > "$OUT/V3_PRJNA720393_bioproject.tsv" 2>>"$LOG"
cat "$OUT/V3_PRJNA720393_bioproject.tsv" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# V4. NC_022143 vs KC593378 -- is the RefSeq a copy of KC593378?
# -----------------------------------------------------------------------------
step "V4 NC_022143 vs KC593378"
efetch -db nuccore -id NC_022143.1 -format gb > "$OUT/V4_NC_022143.gb" 2>>"$LOG"
sed -n '/^COMMENT/,/^FEATURES/p' "$OUT/V4_NC_022143.gb" | tee -a "$LOG"
esearch -db nuccore -query "KC593378.1 OR NC_022143.1" \
  | esummary | xtract -pattern DocumentSummary -element Caption Slen CreateDate Title \
  | tee "$OUT/V4_lengths_dates.tsv" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# V5. GCA_019142605.1 (JJU_Cson_1.0) -- assembly level, N50s, and whether any
#     chromosome-level (assembled-molecule) sequences were submitted
# -----------------------------------------------------------------------------
step "V5 GCA_019142605.1 assembly stats"
datasets summary genome accession GCA_019142605.1 --as-json-lines \
  > "$OUT/V5_GCA_019142605_summary.jsonl" 2>>"$LOG"
grep -oE '"(assembly_level|contig_n50|scaffold_n50|number_of_contigs|number_of_scaffolds|total_number_of_chromosomes|total_sequence_length)":[^,}]*' \
  "$OUT/V5_GCA_019142605_summary.jsonl" | tee -a "$LOG"
datasets summary genome accession GCA_019142605.1 --report sequence --as-json-lines \
  > "$OUT/V5_GCA_019142605_sequence_report.jsonl" 2>>"$LOG" \
  || echo "sequence report option not supported by this datasets version" | tee -a "$LOG"
echo "Sequence roles (count):" | tee -a "$LOG"
grep -oE '"role":"[^"]*"' "$OUT/V5_GCA_019142605_sequence_report.jsonl" | sort | uniq -c | tee -a "$LOG"

# -----------------------------------------------------------------------------
# V6. DQ107916 / DQ107927 (C. sonnerati COI used by Payet et al. 2016)
# -----------------------------------------------------------------------------
step "V6 DQ107916 DQ107927"
efetch -db nuccore -id DQ107916.1,DQ107927.1 -format gb > "$OUT/V6_DQ107916_DQ107927.gb" 2>>"$LOG"
grep -E "^LOCUS|ORGANISM|/country|/specimen_voucher|/gene=|/product=|^  TITLE" \
  "$OUT/V6_DQ107916_DQ107927.gb" | tee -a "$LOG"

echo; echo "Run finished: $(date)" | tee -a "$LOG"
echo "Outputs in $OUT"
