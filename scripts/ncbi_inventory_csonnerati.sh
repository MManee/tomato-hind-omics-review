#!/usr/bin/env bash
# =============================================================================
# NCBI resource inventory for the tomato hind, Cephalopholis sonnerati
# (plus genus-level counts for Cephalopholis, for comparison)
#
# Requirements (one-time install):
#   conda install -c bioconda -c conda-forge entrez-direct ncbi-datasets-cli
#   (python3 is used for the final summary; no extra packages needed)
#
# Optional: an NCBI API key raises the rate limit from 3 to 10 requests/s
#   export NCBI_API_KEY=your_key_here
#
# Usage:
#   bash ncbi_inventory_csonnerati.sh [output_dir]
# =============================================================================
set -uo pipefail

SP="Cephalopholis sonnerati"
GENUS="Cephalopholis"
OUT="${1:-ncbi_inventory_$(date +%Y%m%d)}"
mkdir -p "$OUT"
LOG="$OUT/run.log"
echo "Inventory started $(date)" | tee "$LOG"

need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing tool: $1" | tee -a "$LOG"; exit 1; }; }
need esearch; need efetch; need esummary; need xtract; need datasets; need dataformat; need python3

count() {  # count <db> <query>
  esearch -db "$1" -query "$2" 2>/dev/null | xtract -pattern ENTREZ_DIRECT -element Count
}

# -----------------------------------------------------------------------------
# 0. Record counts per database (species vs genus)
# -----------------------------------------------------------------------------
echo -e "database\tCephalopholis_sonnerati\tCephalopholis_genus" > "$OUT/00_counts.tsv"
for db in assembly sra bioproject biosample nuccore protein gds pubmed; do
  if [ "$db" = "pubmed" ]; then
    qs="(\"Cephalopholis sonnerati\"[Title/Abstract] OR \"tomato hind\"[Title/Abstract] OR \"tomato grouper\"[Title/Abstract])"
    qg="Cephalopholis[Title/Abstract]"
  else
    qs="\"$SP\"[Organism]"
    qg="$GENUS[Organism]"
  fi
  echo -e "$db\t$(count $db "$qs")\t$(count $db "$qg")" >> "$OUT/00_counts.tsv"
  sleep 0.4
done
echo "[0] counts done" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# 1. Taxonomy
# -----------------------------------------------------------------------------
esearch -db taxonomy -query "\"$SP\"[Scientific Name]" | efetch -format xml \
 | xtract -pattern Taxon -element TaxId ScientificName Rank Lineage \
 > "$OUT/01_taxonomy.tsv"
echo "[1] taxonomy done" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# 2. Genome assemblies (species + genus) with key quality metrics
# -----------------------------------------------------------------------------
FIELDS="accession,organism-name,assminfo-name,assminfo-level,assminfo-release-date,\
assminfo-submitter,assminfo-bioproject,assminfo-biosample-accession,\
assmstats-total-sequence-len,assmstats-number-of-contigs,assmstats-contig-n50,\
assmstats-scaffold-n50,assmstats-gc-percent,annotinfo-name,annotinfo-busco-lineage,\
annotinfo-busco-complete"
# If a field name is rejected by your dataformat version, list valid ones with:
#   dataformat tsv genome --help

datasets summary genome taxon "$SP" --as-json-lines \
 | dataformat tsv genome --fields "$FIELDS" > "$OUT/02_assemblies_species.tsv"
datasets summary genome taxon "$GENUS" --as-json-lines \
 | dataformat tsv genome --fields "$FIELDS" > "$OUT/02_assemblies_genus.tsv"
# Annotation (gene counts) if NCBI has annotated any assembly
datasets summary genome taxon "$SP" --as-json-lines \
 | dataformat tsv genome --fields accession,annotinfo-provider,annotinfo-release-date,\
annotinfo-stats-gene-counts-protein-coding,annotinfo-stats-gene-counts-total \
 > "$OUT/02_annotation_species.tsv" 2>/dev/null
echo "[2] assemblies done" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# 3. SRA runs (raw reads: WGS, RNA-seq, Hi-C, amplicon, ...)
# -----------------------------------------------------------------------------
esearch -db sra -query "\"$SP\"[Organism]" | efetch -format runinfo \
 > "$OUT/03_sra_runinfo.csv"
echo "[3] SRA done" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# 4. BioProjects
# -----------------------------------------------------------------------------
esearch -db bioproject -query "\"$SP\"[Organism]" | esummary \
 | xtract -pattern DocumentSummary -element Project_Acc Registration_Date \
     Project_Data_Type Organization Project_Title \
 > "$OUT/04_bioprojects.tsv"
echo "[4] BioProjects done" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# 5. BioSamples with collection locality (key for the Red Sea / range question)
# -----------------------------------------------------------------------------
esearch -db biosample -query "\"$SP\"[Organism]" | efetch -format docsum \
 | xtract -pattern DocumentSummary -element Accession Date \
     -block Attribute -if @attribute_name -equals geo_loc_name -element Attribute \
     -block Attribute -if @attribute_name -equals tissue -element Attribute \
     -block Attribute -if @attribute_name -equals sex -element Attribute \
     -block Attribute -if @attribute_name -equals dev_stage -element Attribute \
 > "$OUT/05_biosamples.tsv"
echo "[5] BioSamples done" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# 6. Nucleotide records (all), then targeted subsets
# -----------------------------------------------------------------------------
esearch -db nuccore -query "\"$SP\"[Organism]" | esummary \
 | xtract -pattern DocumentSummary -element AccessionVersion Slen Biomol CreateDate Title \
 > "$OUT/06_nuccore_all.tsv"

# Complete mitogenomes
esearch -db nuccore -query "\"$SP\"[Organism] AND mitochondrion[Filter] AND \"complete genome\"[Title]" \
 | esummary | xtract -pattern DocumentSummary -element AccessionVersion Slen CreateDate Title \
 > "$OUT/06a_mitogenomes.tsv"

# COI barcodes, with country of origin (needed for distribution check)
esearch -db nuccore -query "\"$SP\"[Organism] AND (COI[Gene] OR COX1[Gene] OR \"cytochrome oxidase subunit I\"[Title] OR \"cytochrome c oxidase subunit 1\"[Title])" \
 | efetch -format gbc \
 | xtract -pattern INSDSeq -element INSDSeq_accession-version INSDSeq_length \
     -block INSDQualifier -if INSDQualifier_name -equals country -element INSDQualifier_value \
     -block INSDQualifier -if INSDQualifier_name -equals lat_lon -element INSDQualifier_value \
     -block INSDReference -position first -element INSDReference_title \
 > "$OUT/06b_COI_with_country.tsv"

# Other markers commonly used in grouper phylogeography
for g in CYTB 16S 12S "control region"; do
  tag=$(echo "$g" | tr ' ' '_')
  esearch -db nuccore -query "\"$SP\"[Organism] AND \"$g\"[All Fields]" \
   | esummary | xtract -pattern DocumentSummary -element AccessionVersion Slen Title \
   > "$OUT/06c_marker_${tag}.tsv"
  sleep 0.4
done

# Transcriptome Shotgun Assemblies (TSA)
esearch -db nuccore -query "\"$SP\"[Organism] AND tsa-master[Properties]" \
 | esummary | xtract -pattern DocumentSummary -element AccessionVersion Slen Title \
 > "$OUT/06d_TSA.tsv"
echo "[6] nucleotide done" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# 7. Proteins
# -----------------------------------------------------------------------------
esearch -db protein -query "\"$SP\"[Organism]" | esummary \
 | xtract -pattern DocumentSummary -element AccessionVersion Slen Title \
 > "$OUT/07_proteins.tsv"
echo "[7] proteins done" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# 8. GEO expression datasets
# -----------------------------------------------------------------------------
esearch -db gds -query "\"$SP\"[Organism]" | esummary \
 | xtract -pattern DocumentSummary -element Accession PDAT gdsType n_samples title \
 > "$OUT/08_geo.tsv"
echo "[8] GEO done" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# 9. PubMed literature list
# -----------------------------------------------------------------------------
esearch -db pubmed -query "\"Cephalopholis sonnerati\"[Title/Abstract] OR \"tomato hind\"[Title/Abstract] OR \"tomato grouper\"[Title/Abstract]" \
 | esummary | xtract -pattern DocumentSummary -element Id PubDate Source \
     -block ArticleId -if IdType -equals doi -element Value -block DocumentSummary -element Title \
 > "$OUT/09_pubmed.tsv"
echo "[9] PubMed done" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# 10. Genus overview: mitogenomes and assemblies per Cephalopholis species
# -----------------------------------------------------------------------------
esearch -db nuccore -query "$GENUS[Organism] AND mitochondrion[Filter] AND \"complete genome\"[Title]" \
 | esummary | xtract -pattern DocumentSummary -element Organism AccessionVersion Slen \
 | sort > "$OUT/10_genus_mitogenomes.tsv"
echo "[10] genus overview done" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# 11. Human-readable summary
# -----------------------------------------------------------------------------
python3 - "$OUT" <<'PY' > "$OUT/SUMMARY.txt"
import csv, sys, os, collections
out = sys.argv[1]
def rows(f):
    p = os.path.join(out, f)
    if not os.path.exists(p) or os.path.getsize(p) == 0: return []
    with open(p) as h: return [l.rstrip("\n").split("\t") for l in h if l.strip()]

print("=== Record counts (species | genus) ===")
for r in rows("00_counts.tsv"): print("  " + " | ".join(r))

print("\n=== Genome assemblies (species) ===")
for r in rows("02_assemblies_species.tsv"): print("  " + "\t".join(r))

print("\n=== SRA runs by strategy / platform / BioProject ===")
p = os.path.join(out, "03_sra_runinfo.csv")
if os.path.exists(p) and os.path.getsize(p) > 0:
    with open(p) as h:
        recs = [r for r in csv.DictReader(h) if r.get("Run") and r["Run"] != "Run"]
    total_gb = sum(float(r.get("bases") or 0) for r in recs) / 1e9
    print(f"  Runs: {len(recs)}   Total bases: {total_gb:,.1f} Gb")
    for key in ("LibraryStrategy", "LibrarySource", "Platform", "BioProject"):
        c = collections.Counter(r.get(key, "NA") for r in recs)
        print(f"  {key}: " + ", ".join(f"{k}={v}" for k, v in c.most_common()))
    bp = collections.defaultdict(lambda: collections.Counter())
    for r in recs: bp[r.get("BioProject","NA")][r.get("LibraryStrategy","NA")] += 1
    print("  Per BioProject:")
    for k, v in bp.items(): print(f"    {k}: " + ", ".join(f"{s}={n}" for s, n in v.items()))

print("\n=== BioSample localities ===")
c = collections.Counter((r[2] if len(r) > 2 else "not given") for r in rows("05_biosamples.tsv"))
for k, v in c.most_common(): print(f"  {k}: {v}")

print("\n=== COI records by country ===")
c = collections.Counter((r[2] if len(r) > 2 and r[2] else "not given") for r in rows("06b_COI_with_country.tsv"))
for k, v in c.most_common(): print(f"  {k}: {v}")

for f, label in [("06a_mitogenomes.tsv","Mitogenomes"),("06d_TSA.tsv","TSA assemblies"),
                 ("08_geo.tsv","GEO datasets"),("04_bioprojects.tsv","BioProjects")]:
    print(f"\n=== {label} ===")
    for r in rows(f): print("  " + "\t".join(r))

print("\n=== Cephalopholis species with complete mitogenomes ===")
c = collections.Counter(r[0] for r in rows("10_genus_mitogenomes.tsv"))
for k, v in sorted(c.items()): print(f"  {k}: {v}")
PY

echo "Finished $(date). See $OUT/SUMMARY.txt" | tee -a "$LOG"
