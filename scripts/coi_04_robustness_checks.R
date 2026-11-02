# =============================================================================
# COI barcode analysis, step 4: robustness checks of the curation.
# Reads the outputs of coi_03_analysis.R (no re-alignment, no download) and
# writes results/coi_robustness_checks/:
#
#   1  coi_record_ledger.tsv           one row per record: original name, link to
#                                      GenBank/BOLD, QC outcome, flags, exclusion
#                                      reason and final analytical status
#      bold_linked_records.tsv         BOLD records counted once with their GenBank record
#   2  misid_rule_all_species.tsv      the "within 1% of another named species" rule
#                                      applied to every species, not only C. sonnerati
#      misid_rule_all_species_summary.tsv
#      french_polynesia_check.tsv      why JQ431574/JQ431575/MK566843 were not caught
#                                      by the nearest-neighbour rule
#   3  csonnerati_not_in_network.tsv   the C. sonnerati records omitted from the
#                                      580-bp haplotype network, with missing bases
#   4  qc_threshold_sensitivity.tsv    QC counts under alternative thresholds,
#                                      including the WIO and India BIN groups
#      misid_threshold_sensitivity.tsv C. sonnerati records set aside at 0.5/1/2%
#   5  csonnerati_shared_specimen_check.tsv  C. sonnerati records sharing a voucher,
#                                      isolate, BOLD sample ID or identical sequence
#                                      within the same study
#   6  robustness_checks_summary.txt       all headline numbers
#
# Run from Terminal (after coi_03_analysis.R has been run):
#   cd tomato-hind-omics-review   # repository root
#   Rscript scripts/coi_04_robustness_checks.R
# =============================================================================

suppressPackageStartupMessages({ library(ape); library(Biostrings) })

proj    <- Sys.getenv("PROJ", ".")
in_dir  <- file.path(proj, "data/coi_20261008")
res_dir <- file.path(proj, "results/coi_analysis")
out_dir <- file.path(proj, "results/coi_robustness_checks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sumf <- file.path(out_dir, "robustness_checks_summary.txt"); unlink(sumf)
say <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
wt  <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE, row.names = FALSE)
need <- function(f) { if (!file.exists(f)) stop("Missing input: ", f, "\nRun coi_03_analysis.R first."); f }

focal <- "Cephalopholis sonnerati"
say("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "|", R.version.string,
    "| ape", as.character(packageVersion("ape")), "| Biostrings", as.character(packageVersion("Biostrings")))

## ---- inputs --------------------------------------------------------------------
qc   <- read.delim(need(file.path(res_dir, "coi_qc_all_records.tsv")), colClasses = "character", quote = "")
for (v in c("pid", "cov_window", "ambig", "stops")) qc[[v]] <- as.numeric(qc[[v]])
qc$frameshift <- qc$frameshift == "TRUE"
qc$species[qc$species %in% c("NA", "")] <- NA
scr  <- read.delim(need(file.path(res_dir, "coi_nearest_neighbour_screen_all.tsv")), colClasses = "character", quote = "")
mis  <- read.delim(need(file.path(res_dir, "csonnerati_probable_misidentifications.tsv")), colClasses = "character", quote = "")
drop <- readLines(need(file.path(res_dir, "bold_genbank_mined_unlinked_dropped.txt")))
gb   <- read.delim(need(file.path(in_dir, "coi_genbank.tsv")), colClasses = "character", quote = "")
bold_file <- file.path(in_dir, "bold_records.tsv")
bold <- if (file.exists(bold_file)) read.delim(bold_file, colClasses = "character", quote = "") else NULL

fa  <- readDNAStringSet(need(file.path(res_dir, "coi_window_alignment_qc_pass.fasta")))
hdr <- strsplit(names(fa), "|", fixed = TRUE)
ids <- sapply(hdr, `[`, 1); spp <- gsub("_", " ", sapply(hdr, `[`, 2))
X   <- do.call(rbind, strsplit(as.character(fa), ""))
rownames(X) <- ids
D   <- as.matrix(dist.dna(as.DNAbin(tolower(X)), model = "K80", pairwise.deletion = TRUE))

excl  <- scr$id[grepl("closer_to_another_species|singleton_within", scr$flag)]   # as in step 3
excl2 <- union(excl, mis$id)

## ---- 1. record ledger -------------------------------------------------------
gbx <- gb[, intersect(c("accession", "specimen_voucher", "isolate", "first_reference"), names(gb))]
led <- merge(qc[, c("id", "source", "accession", "processid", "species_raw", "species", "country", "bin",
                    "pid", "cov_window", "ambig", "stops", "frameshift", "qc")],
             gbx, by = "accession", all.x = TRUE, sort = FALSE)
if (!is.null(bold)) {
  bs <- bold[, c("processid", "sampleid")]
  led <- merge(led, bs, by.x = "id", by.y = "processid", all.x = TRUE, sort = FALSE)
} else led$sampleid <- NA
led <- merge(led, scr[, c("id", "min_conspecific_pct", "nearest_other_species", "nearest_other_id",
                          "min_other_species_pct", "flag")], by = "id", all.x = TRUE, sort = FALSE)
led$flag[is.na(led$flag)] <- ""
led$status <- with(led, ifelse(qc != "pass", paste0("qc_fail:", qc),
                   ifelse(is.na(species), "pass_unresolved_name",
                   ifelse(id %in% excl, "excluded_nearest_neighbour",
                   ifelse(id %in% mis$id, "set_aside_probable_misidentification",
                          "included_species_level")))))
led$retained_flag_only <- led$status == "included_species_level" & led$flag != ""
dup <- data.frame(id = drop, source = "BOLD", status = "dropped_duplicate_of_GenBank (GB-mined process ID)")
led <- merge(led, dup, all = TRUE, sort = FALSE)
led <- led[order(led$status, led$species_raw, led$id), ]

say("\n== 1. Record ledger (coi_record_ledger.tsv) ==")
st <- table(led$status); for (k in names(st)) say(sprintf("  %-60s %5d", k, st[[k]]))
npu <- sum(qc$qc == "pass" & is.na(qc$species))
say("  QC-passed records with an unresolved name (for Table S1):", npu)
pu <- qc$qc == "pass" & is.na(qc$species); pr <- qc$qc == "pass" & !is.na(qc$species)
say("  ... by source (GenBank / BOLD): unresolved", sum(pu & qc$source == "GenBank"), "/", sum(pu & qc$source == "BOLD"),
    "| resolved", sum(pr & qc$source == "GenBank"), "/", sum(pr & qc$source == "BOLD"))
say("  Unresolved names among all aligned records (run log value):", sum(is.na(qc$species)))
say("  QC-passed, resolved name:", sum(qc$qc == "pass" & !is.na(qc$species)),
    "in", length(unique(qc$species[qc$qc == "pass" & !is.na(qc$species)])), "species")
inc <- led$status == "included_species_level"
say("  Included in species-level statistics:", sum(inc), "in", length(unique(led$species[inc])), "species")
lost <- setdiff(unique(qc$species[qc$qc == "pass" & !is.na(qc$species)]), unique(led$species[inc]))
say("  Species with no record left after exclusions:", if (length(lost)) paste(lost, collapse = ", ") else "none")
say("  Flagged but retained (far from conspecifics only):",
    if (any(led$retained_flag_only, na.rm = TRUE)) paste(led$id[which(led$retained_flag_only)], collapse = ", ") else "none")
wt(led, "coi_record_ledger.tsv")
if (!is.null(bold) && "processid" %in% names(qc)) {
  lk <- qc[qc$source == "GenBank" & qc$processid != "", c("accession", "processid", "species_raw")]
  wt(lk, "bold_linked_records.tsv")
  say("  BOLD records linked to a GenBank record (bold_linked_records.tsv):", nrow(lk))
}

## ---- 2. 1% rule applied to all species -------------------------------------
keep <- which(!ids %in% excl)
rows <- lapply(keep, function(k) {
  o <- keep[spp[keep] != spp[k]]; c_ <- setdiff(keep[spp[keep] == spp[k]], k)
  dh <- D[k, o]; j <- o[which.min(dh)]
  dc <- if (length(c_)) min(D[k, c_], na.rm = TRUE) else NA
  data.frame(id = ids[k], species = spp[k], nearest_other_species = spp[j], nearest_other_id = ids[j],
             min_other_species_pct = abs(round(100 * min(dh, na.rm = TRUE), 2)),
             min_conspecific_pct = abs(round(100 * dc, 2)), stringsAsFactors = FALSE)
})
ar <- do.call(rbind, rows)
ar1 <- ar[ar$min_other_species_pct < 1, ]
wt(ar1[order(ar1$species, ar1$min_other_species_pct), ], "misid_rule_all_species.tsv")
sm <- if (nrow(ar1)) aggregate(id ~ species + nearest_other_species, ar1, length) else
  data.frame(species = character(), nearest_other_species = character(), id = integer())
names(sm)[3] <- "n_records_within_1pct"
wt(sm[order(-sm$n_records_within_1pct), ], "misid_rule_all_species_summary.tsv")
say("\n== 2. 'Within 1% of another named species' applied to all species (after nearest-neighbour exclusions) ==")
say("  Records affected:", nrow(ar1), "in", length(unique(ar1$species)), "species")
for (i in seq_len(nrow(sm))) say(sprintf("  %-28s vs %-28s %4d", sm$species[i], sm$nearest_other_species[i], sm$n_records_within_1pct[i]))
fp <- ar[ar$id %in% mis$id, ]
fp$nearest_conspecific_id <- sapply(fp$id, function(x) {
  k <- match(x, ids); c_ <- setdiff(keep[spp[keep] == spp[k]], k); ids[c_][which.min(D[k, c_])] })
wt(fp, "french_polynesia_check.tsv")
say("  Records set aside for C. sonnerati: nearest conspecific and distance:")
for (i in seq_len(nrow(fp))) say(sprintf("    %s: nearest conspecific %s at %.2f%%; nearest other species %s (%s) at %.2f%%",
  fp$id[i], fp$nearest_conspecific_id[i], fp$min_conspecific_pct[i], fp$nearest_other_species[i],
  fp$nearest_other_id[i], fp$min_other_species_pct[i]))

## ---- 3. C. sonnerati records not in the network ----------------------------------
fs <- which(spp == focal & !ids %in% excl2)
Xs <- X[fs, , drop = FALSE]
okcol <- colMeans(matrix(Xs %in% c("A","C","G","T"), nrow = nrow(Xs))) >= 0.95   # as in step 3
rr <- rle(okcol); en <- cumsum(rr$lengths); kk <- which(rr$values)[which.max(rr$lengths[rr$values])]
core <- (en[kk] - rr$lengths[kk] + 1):en[kk]
nmiss <- apply(Xs[, core, drop = FALSE], 1, function(x) sum(!x %in% c("A","C","G","T")))
nn <- data.frame(id = ids[fs], n_core_positions_missing_or_ambiguous = nmiss, stringsAsFactors = FALSE)
nn <- merge(nn[nn$n_core_positions_missing_or_ambiguous > 0, ],
            led[, c("id", "country", "bin", "specimen_voucher", "first_reference")], by = "id", all.x = TRUE)
wt(nn, "csonnerati_not_in_network.tsv")
say("\n== 3. C. sonnerati records omitted from the haplotype network ==")
say(sprintf("  Curated C. sonnerati: %d; core = window columns %d-%d (%d bp); omitted: %d",
            length(fs), min(core), max(core), length(core), nrow(nn)))
for (i in seq_len(nrow(nn))) say(sprintf("    %s (%s; %s): %d core positions missing or ambiguous",
  nn$id[i], nn$country[i], ifelse(is.na(nn$bin[i]) | nn$bin[i] == "", "no BIN", nn$bin[i]),
  nn$n_core_positions_missing_or_ambiguous[i]))

## ---- 4. threshold sensitivity ----------------------------------------------------
grid <- expand.grid(min_cov = c(0.80, 0.90, 0.95), max_amb = c(0.005, 0.01, 0.02), min_pid = c(75, 80, 85))
sens <- do.call(rbind, lapply(seq_len(nrow(grid)), function(g) {
  p <- grid[g, ]
  ok <- with(qc, pid >= p$min_pid & cov_window >= p$min_cov & ambig <= p$max_amb & stops == 0 & !frameshift)
  ok[is.na(ok)] <- FALSE
  cs <- ok & qc$species %in% focal
  data.frame(p, pass = sum(ok), pass_named = sum(ok & !is.na(qc$species)),
             n_species = length(unique(qc$species[ok & !is.na(qc$species)])),
             csonnerati_pass = sum(cs),
             csonnerati_countries = length(unique(qc$country[cs & qc$country != "Not recorded"])),
             csonnerati_BIN_AAB5431 = sum(cs & qc$bin == "BOLD:AAB5431"),
             csonnerati_BIN_ACC0647 = sum(cs & qc$bin == "BOLD:ACC0647"),
             csonnerati_BIN_AAD4737 = sum(cs & qc$bin == "BOLD:AAD4737"))
}))
sens$used_in_manuscript <- with(sens, min_cov == 0.90 & max_amb == 0.01 & min_pid == 80)
wt(sens, "qc_threshold_sensitivity.tsv")
say("\n== 4a. QC threshold sensitivity (qc_threshold_sensitivity.tsv) ==")
say("  Ranges across the 27 combinations:")
for (v in c("pass", "pass_named", "n_species", "csonnerati_pass", "csonnerati_countries",
            "csonnerati_BIN_AAB5431", "csonnerati_BIN_ACC0647"))
  say(sprintf("  %-24s %5d - %5d (manuscript setting: %d)", v, min(sens[[v]]), max(sens[[v]]),
              sens[[v]][sens$used_in_manuscript]))
fi <- which(spp == focal & !ids %in% excl); oi <- which(spp != focal & !ids %in% excl)
dmin <- apply(D[fi, oi, drop = FALSE], 1, min, na.rm = TRUE)
mt <- data.frame(threshold_pct = c(0.5, 1, 2),
                 n_set_aside = sapply(c(0.005, 0.01, 0.02), function(t) sum(dmin < t)),
                 records = sapply(c(0.005, 0.01, 0.02), function(t) paste(ids[fi][dmin < t], collapse = ";")))
wt(mt, "misid_threshold_sensitivity.tsv")
say("\n== 4b. C. sonnerati records set aside at alternative thresholds ==")
for (i in seq_len(nrow(mt))) say(sprintf("  < %.1f%%: %d (%s)", mt$threshold_pct[i], mt$n_set_aside[i], mt$records[i]))

## ---- 5. possible repeated specimens among C. sonnerati ---------------------------
cs <- led[led$status == "included_species_level" & led$species %in% focal, ]
cs$seq <- apply(X[match(cs$id, ids), , drop = FALSE], 1, paste, collapse = "")
key <- function(v) ifelse(is.na(v) | v == "", NA, v)
grp <- list(voucher = key(cs$specimen_voucher), isolate = key(cs$isolate), bold_sampleid = key(cs$sampleid),
            same_study_identical_seq = ifelse(is.na(key(cs$first_reference)), NA,
                                              paste(cs$first_reference, cs$seq)))
sh <- do.call(rbind, lapply(names(grp), function(g) {
  v <- grp[[g]]; d <- v[!is.na(v) & duplicated(v)]
  if (!length(d)) return(NULL)
  do.call(rbind, lapply(unique(d), function(x) data.frame(criterion = g,
    ids = paste(cs$id[which(v == x)], collapse = ";"), n = sum(v == x, na.rm = TRUE),
    shared_value = if (g == "same_study_identical_seq") substr(cs$first_reference[which(v == x)[1]], 1, 80) else x)))
}))
say("\n== 5. Possible repeated specimens among curated C. sonnerati ==")
if (is.null(sh)) { say("  none found"); wt(data.frame(criterion = character()), "csonnerati_shared_specimen_check.tsv") } else {
  wt(sh, "csonnerati_shared_specimen_check.tsv")
  for (i in seq_len(nrow(sh))) say(sprintf("  %s (n = %d): %s | %s", sh$criterion[i], sh$n[i], sh$ids[i], sh$shared_value[i]))
}
say("  Records with a voucher:", sum(!is.na(key(cs$specimen_voucher))), "of", nrow(cs),
    "| with a BOLD sample ID:", sum(!is.na(key(cs$sampleid))))
say("\nNote: identical sequences from the same study may be different fish; they are listed for checking, not removed.")
say("\nDone. Summary written to:", sumf)
