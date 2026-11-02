# Manuscript items and the files that produce them

Numbering follows the manuscript as prepared for submission. Values were checked against these files.

| Manuscript item | Script | Output file(s) |
|---|---|---|
| Table 1. Distances of the Saudi Arabian 16S/12S records | `box1_16S_tree.R` (16S *p*-distances); `verify_saudi_records.py` with exclusion `KX2986,KX2987` (12S, 100 minus percent identity) | `results/box1_tree/box1_query_distances.tsv`; `data/ncbi_inventory_20261001/verify_saudi/saudi_records_verification_excl_KX2986_KX2987.tsv` and `*_best_per_species_excl_KX2986_KX2987.tsv` |
| Figure 1. 16S rRNA maximum-likelihood phylogeny | `box1_16S_tree.R` | `results/box1_tree/Box1_16S_ML_tree_full.pdf`, `box1_16S_ML_tree.nwk`, `box1_16S_alignment_trimmed.fasta`, `box1_tip_metadata.tsv` |
| Box text: the 2016 Saudi Arabian submission (accessions, vouchers, replicate suffixes) | `kx_records_info.sh`, `kx_records_parse.py` | `data/kx2986_check/kx_records_summary.tsv`, `kx_records.gb` |
| Table 2. Genomic, transcriptomic and mitochondrial resources | `ncbi_inventory_csonnerati.sh`, `ncbi_followup_csonnerati.sh`, `verify_ncbi_20261007.sh`; plus the cited publications | `data/ncbi_inventory_20261001/` (`SUMMARY.txt`, `02_assemblies_species.tsv`, `followup/mito_markers_summary.txt`, `followup/mito_markers_by_locality.tsv`); `data/verification_20261007/` |
| Table 3. Genome assembly statistics | `ncbi_inventory_csonnerati.sh`, `verify_ncbi_20261007.sh` (V2, V5); plus the cited publications | `data/ncbi_inventory_20261001/02_assemblies_species.tsv`; `data/verification_20261007/V5_GCA_019142605_*.jsonl`, `V2_SAMN49597358_full.txt` |
| Table 4. SRA runs by BioProject; 69 runs, 1068.8 Gb; 24 runs with library strategy WGS but transcriptomic source | `ncbi_inventory_csonnerati.sh`, `ncbi_followup_csonnerati.sh` | `data/ncbi_inventory_20261001/03_sra_runinfo.csv`, `followup/sra_by_project.tsv` |
| 16 gut-microbiome runs (PRJNA1433391) | `verify_ncbi_20261007.sh` (V1) | `data/verification_20261007/V1_PRJNA1433391_*` |
| COI dataset: 1522 sequences, 1356 passing QC, 1276 from 23 named species; nearest-neighbour screen; BINs | `coi_01_download.sh`, `coi_01b_bold.py`, `coi_02_parse_genbank.py`, `coi_03_analysis.R` | `data/coi_20261008/download_log.txt`, `parse_log.txt`; `results/coi_analysis/run_log.txt`, `coi_qc_all_records.tsv`, `coi_flagged_sequences.tsv`, `coi_bins_per_species.tsv`, `coi_species_per_bin.tsv`, `coi_genbank_bold_name_conflicts.tsv` |
| Threshold sensitivity (27 combinations), record ledger, 1% rule for all species, repeated specimens | `coi_04_robustness_checks.R` | `results/coi_robustness_checks/` |
| Extreme within-species distances, diversity for 104 and 102 sequences, synonymy check, *C. oligosticta* reference provenance | `coi_05_additional_checks.R` | `results/coi_additional_checks/` |
| Figure 2. Map of record origins | `coi_06_map.R` | `results/coi_map/Fig_Csonnerati_records_map.pdf`, `map_points.tsv`, `map_summary.txt` |
| Figure 3. *C. sonnerati* COI haplotype network | `coi_03_analysis.R` | `results/coi_analysis/Fig_Csonnerati_haplotype_network_BIN.pdf`, `csonnerati_haplotypes.tsv`, `csonnerati_diversity.tsv` |
| Table 5. Candidate positions distinguishing *C. sonnerati* from *C. oligosticta* | `coi_03_analysis.R` (positions, bases, counts); `coi_05_additional_checks.R` (amino acids per sequence) | `results/coi_analysis/diagnostic_sites_sonnerati_vs_oligosticta_curated.tsv` (also `_all_qc_pass`, `_excluding_flagged`); `results/coi_additional_checks/synonymy_check.tsv` |
| Table 6. FAIR assessment | none (synthesis of the results above) | — |
| Table S1. Assembly and quality control of the COI dataset | `coi_03_analysis.R`, `coi_04_robustness_checks.R` | `results/coi_analysis/run_log.txt`, `coi_qc_all_records.tsv`; `results/coi_robustness_checks/coi_record_ledger.tsv` |
| Table S2. Within- and between-species K2P distances | `coi_03_analysis.R` | `results/coi_analysis/coi_species_distances_curated.tsv` (also `coi_species_distances.tsv`, `coi_species_distances_excluding_flagged.tsv`) |
| Table S3. *C. sonnerati* barcodes conflicting with their label | `coi_03_analysis.R`, `coi_04_robustness_checks.R` | `results/coi_analysis/csonnerati_probable_misidentifications.tsv`, `coi_flagged_sequences.tsv`; `results/coi_robustness_checks/french_polynesia_check.tsv` |
| Table S4. *C. sonnerati* COI haplotypes | `coi_03_analysis.R` | `results/coi_analysis/csonnerati_haplotypes.tsv` |

Additional outputs not shown in the manuscript: `Fig_COI_distance_histogram`, `Fig_COI_barcode_gap_*`, `Fig_COI_NJ_tree.pdf` and `coi_NJ_unique_haplotypes.nwk` (`coi_03_analysis.R`); `Box1_16S_ML_tree_focus` (`box1_16S_tree.R`).
