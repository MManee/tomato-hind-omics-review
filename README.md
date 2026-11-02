# Tomato hind (*Cephalopholis sonnerati*): public sequence data audit and re-analysis

Scripts, data snapshots and results supporting the review

> Manee, M.M.; Alqahtani, M.S.; Al-Shomrani, B.M.; Alqahtani, F.H. *Omics Research in the Tomato Hind (*Cephalopholis sonnerati*): Current Advances, Knowledge Gaps and Future Perspectives.* Manuscript prepared for *Fishes* (MDPI); citation details will be added on publication.

The repository reproduces the analyses of public data reported in the review:

1. an inventory of NCBI resources for *C. sonnerati* (genome assemblies, SRA runs, BioProjects, BioSamples, mitochondrial records);
2. a re-evaluation of the GenBank 16S and 12S rRNA records deposited as *C. sonnerati* from Saudi Arabia (maximum-likelihood phylogeny and distances);
3. a curated range-wide COI barcode reference for the genus *Cephalopholis* from GenBank and BOLD, with quality control, distance and BIN analyses, a *C. sonnerati* haplotype network, candidate positions distinguishing *C. sonnerati* from *C. oligosticta*, robustness checks and a map of record origins.

All data analysed are public. No new sequences were generated.

## Repository layout

```
scripts/   analysis scripts (Bash, Python, R)
data/      snapshots of the public records used, with retrieval logs
results/   tables and figures produced by the scripts
docs/      software versions and the link between manuscript items and outputs
```

## Software

Analyses were run on macOS 12.7.6 with:

| Tool | Version |
|---|---|
| R | 4.6.0 |
| ape | 5.8.1 |
| Biostrings | 2.80.2 |
| pwalign | 1.8.0 |
| DECIPHER | 3.8.1 |
| phangorn | 2.12.1 |
| ggtree | 4.2.0 |
| pegas | 1.4 |
| ggplot2 | 4.0.3 |
| Python | 3.12.15 |
| Biopython | 1.88 |
| Entrez Direct | 26.3 |
| NCBI Datasets CLI | 18.38.0 |

Full R session information is in `results/coi_analysis/sessionInfo.txt` and `results/box1_tree/sessionInfo.txt`; see also `docs/software_versions.md`. The version of the R package `maps` (country outlines in the map) was not recorded.

Installation (one time):

```r
install.packages(c("ape", "phangorn", "ggplot2", "pegas", "maps", "BiocManager"))
BiocManager::install(c("Biostrings", "pwalign", "DECIPHER", "ggtree"))
```

```bash
pip install biopython==1.88
# Entrez Direct: https://www.ncbi.nlm.nih.gov/books/NBK179288/
# NCBI Datasets CLI: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/
```

## How to run

All commands are run from the repository root unless stated otherwise. R scripts read the environment variable `PROJ` (default: the current directory).

### Re-running the analyses on the archived data (no download)

```bash
Rscript scripts/coi_03_analysis.R            # COI curation, distances, BINs, network, candidate positions
Rscript scripts/coi_04_robustness_checks.R   # record ledger, threshold sensitivity, repeated specimens
Rscript scripts/coi_05_additional_checks.R   # extreme distances, diversity, synonymy, reference provenance
Rscript scripts/coi_06_map.R                 # map of record origins
Rscript scripts/box1_16S_tree.R              # 16S rRNA phylogeny of the Saudi Arabian records
```

`coi_03_analysis.R` and `coi_04_robustness_checks.R` read `data/coi_20261008/bold_records.tsv` if present, and `coi_04`–`coi_06` read the alignment `results/coi_analysis/coi_window_alignment_qc_pass.fasta` written by `coi_03`. Because these files contain BOLD sequences, they are not distributed here (see *Data and licences*). Without `bold_records.tsv` the COI analysis runs on GenBank records only and the results differ from those reported; to reproduce the full analysis, first download the BOLD records with `python3 scripts/coi_01b_bold.py data/coi_20261008`.

**Reproducibility check.** All scripts above were rerun from a fresh copy of this repository on 2 November 2026. All COI outputs and 16S distances matched the original run, apart from run dates, file paths in log messages, and the amino-acid labels at COI positions 441 and 687, which were corrected from Met to Ile after single codons were translated without start-codon rules (`no.init.codon = TRUE`). The 16S tree had the same topology, but bootstrap values differed by up to 2 percentage points (for example, 66% versus 65% for the clade containing all *C. oligosticta* sequences), because `box1_16S_tree.R` runs the bootstrap replicates in parallel (`multicore = TRUE`), and `set.seed()` does not fix the random streams of parallel workers. The tree files in `results/box1_tree/` are those of the original run reported in the manuscript.

### Re-downloading the data

Public databases change over time, so a new download may not return exactly the records used here. The GenBank COI query is fixed by a publication-date cutoff (`PDAT <= 2026/10/08`); BOLD records cannot be frozen in the same way, and the process IDs used are listed in `data/coi_20261008/bold_processids.tsv` for comparison.

```bash
# COI barcodes (GenBank + BOLD), then parse GenBank
bash scripts/coi_01_download.sh                      # runs coi_01b_bold.py for BOLD
python3 scripts/coi_02_parse_genbank.py data/coi_20261008

# NCBI inventory (run from data/)
cd data
bash ../scripts/ncbi_inventory_csonnerati.sh ncbi_inventory_20261001
bash ../scripts/ncbi_followup_csonnerati.sh ncbi_inventory_20261001
python3 ../scripts/verify_saudi_records.py ncbi_inventory_20261001/verify_saudi
python3 ../scripts/verify_saudi_records.py ncbi_inventory_20261001/verify_saudi KX2986
python3 ../scripts/verify_saudi_records.py ncbi_inventory_20261001/verify_saudi KX2986,KX2987   # excludes the whole 2016 submission (used in Table 1)
cd ..

# Full GenBank records of the Saudi Arabian submission, and single-record checks
bash scripts/kx_records_info.sh
bash scripts/verify_ncbi_20261007.sh
```

`scripts/parse_mito_markers.py` is a stand-alone parser for the mitochondrial-record XML written by the follow-up script; the follow-up script contains its own parser and does not call it.

## Link between the manuscript and the outputs

See `docs/manuscript_outputs.md` for each table, figure and supplementary table, the script that produces it and the output file.

## Data and licences

* **Code** (`scripts/`): MIT licence (`LICENSE`).
* **Results produced by this study** (`results/`, `docs/`): CC BY 4.0 (`LICENSE-DATA.md`).
* **Third-party records** (`data/`): GenBank, SRA, BioProject, BioSample and Assembly records retrieved from NCBI, distributed unchanged as retrieved. They remain subject to the terms of their source databases and submitters (NCBI: https://www.ncbi.nlm.nih.gov/home/about/policies/).
* **BOLD records** are not redistributed. Only their process IDs and selected metadata (no sequences) are listed in `data/coi_20261008/bold_processids.tsv`; the records can be downloaded with `scripts/coi_01b_bold.py`.

## Citation

Please cite the review (citation details will be added on publication) and this repository (`CITATION.cff`).

## Contact

Manee M. Manee (malmanee@kacst.gov.sa) and Fahad H. Alqahtani (fqahtani@kacst.gov.sa), King Abdulaziz City for Science and Technology (KACST), Riyadh, Saudi Arabia.
