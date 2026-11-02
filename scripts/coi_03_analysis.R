# =============================================================================
# COI barcode analysis, step 3: range-wide COI reference set for
# Cephalopholis sonnerati (GenBank + BOLD), with
#   (1) quality control (coverage, ambiguity, stop codons, frameshifts);
#   (2) within- and between-species K2P distances and the barcode gap;
#   (3) a screen of sequences whose nearest neighbour is another species;
#   (4) BOLD BIN concordance;
#   (5) a C. sonnerati haplotype network coloured by country (pegas);
#   (6) fixed nucleotide differences between C. sonnerati and C. oligosticta;
#   (7) a neighbour-joining tree of unique haplotypes (bootstrap).
#
# Positions are given in the coordinates of the COI CDS of KC593378
# (C. sonnerati mitogenome; position 1 = first base of the annotated CDS).
#
# Inputs : data/coi_20261008/coi_genbank.tsv  (step 2)
#          data/coi_20261008/bold_records.tsv (step 1b; optional)
# Output : results/coi_analysis/
#
# Run from Terminal:
#   cd tomato-hind-omics-review   # repository root
#   Rscript scripts/coi_03_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(ape); library(Biostrings); library(ggplot2)
})
if (!requireNamespace("pegas", quietly = TRUE))
  install.packages("pegas", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(library(pegas))
pa <- if (requireNamespace("pwalign", quietly = TRUE)) pwalign::pairwiseAlignment else
  Biostrings::pairwiseAlignment
pid_fun <- if (requireNamespace("pwalign", quietly = TRUE)) pwalign::pid else Biostrings::pid
nsm <- if (requireNamespace("pwalign", quietly = TRUE)) pwalign::nucleotideSubstitutionMatrix else
  Biostrings::nucleotideSubstitutionMatrix
# alignment accessors moved from Biostrings to pwalign in Bioconductor 3.19+
pw <- requireNamespace("pwalign", quietly = TRUE)
aln_pat <- if (pw) pwalign::alignedPattern else Biostrings::alignedPattern
aln_sub <- if (pw) pwalign::alignedSubject else Biostrings::alignedSubject
aln_score <- if (pw) pwalign::score else Biostrings::score
aln_pstart <- function(a) start(if (pw) pwalign::pattern(a) else Biostrings::pattern(a))

## ---- paths & parameters (all reported in parameters.txt) ---------------------
proj    <- Sys.getenv("PROJ", ".")
in_dir  <- file.path(proj, "data/coi_20261008")
out_dir <- file.path(proj, "results/coi_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

P <- list(
  ref_acc      = "KC593378",   # C. sonnerati mitogenome: COI reference frame
  focal        = "Cephalopholis sonnerati",
  sister       = "Cephalopholis oligosticta",
  sister_refs  = c("KM656840", "MH707291"),
  min_pid      = 80,     # % identity to the reference COI over the aligned part
  win_frac     = 0.80,   # window = longest stretch covered by >= 80% of the max coverage
  min_cov      = 0.90,   # sequence must cover >= 90% of the window
  max_amb      = 0.01,   # <= 1% ambiguous bases within the covered window
  thr          = 0.02,   # K2P distance used to flag high divergence (2%)
  core_frac    = 0.95,   # network core: columns complete in >= 95% of C. sonnerati
  n_boot       = 500,    # NJ bootstrap replicates
  seed         = 2026)
set.seed(P$seed)
if (nzchar(Sys.getenv("WIN_START"))) P$win_override <- as.integer(c(Sys.getenv("WIN_START"), Sys.getenv("WIN_END")))
msg <- function(...) { cat(..., "\n"); cat(..., "\n", file = file.path(out_dir, "run_log.txt"), append = TRUE) }
unlink(file.path(out_dir, "run_log.txt"))

## ---- 1. read and combine GenBank + BOLD ---------------------------------------
gb <- read.delim(file.path(in_dir, "coi_genbank.tsv"), colClasses = "character", quote = "")
gb$acc_nov <- sub("\\.\\d+$", "", gb$accession)
gb$bold_pid <- sub("\\..*$", "", sub("^BOLD:", "", sapply(strsplit(gb$bold_xref, ";"), `[`, 1)))
gb$bold_pid[is.na(gb$bold_pid)] <- ""
msg("GenBank COI sequences read:", nrow(gb))

bold_file <- file.path(in_dir, "bold_records.tsv")
bold <- if (file.exists(bold_file)) read.delim(bold_file, colClasses = "character", quote = "") else NULL
bold_status <- if (file.exists(file.path(in_dir, "bold_status.txt")))
  readLines(file.path(in_dir, "bold_status.txt")) else "bold_status.txt not found"
msg("BOLD status:", bold_status)

gb$bin <- ""; gb$bold_species <- ""; gb$bold_processid <- ""
name_conflicts <- NULL
if (!is.null(bold) && nrow(bold) > 0) {
  msg("BOLD COI records read:", nrow(bold))
  bold$insdc1 <- sub("\\.\\d+$", "", sapply(strsplit(bold$insdc_acs, "[,;| ]+"), `[`, 1))
  bold$insdc1[is.na(bold$insdc1)] <- ""
  # link BOLD records to GenBank records (by INSDC accession or by BOLD process ID)
  link <- ifelse(bold$insdc1 != "" & bold$insdc1 %in% gb$acc_nov, match(bold$insdc1, gb$acc_nov),
          ifelse(bold$processid %in% gb$bold_pid[gb$bold_pid != ""], match(bold$processid, gb$bold_pid), NA))
  hit <- which(!is.na(link))
  gb$bin[link[hit]] <- bold$bin_uri[hit]
  gb$bold_species[link[hit]] <- bold$species[hit]
  gb$bold_processid[link[hit]] <- bold$processid[hit]
  msg("BOLD records linked to a GenBank record (counted once, as GenBank):", length(hit))
  nc <- link[hit][gb$bold_species[link[hit]] != "" & gb$bold_species[link[hit]] != gb$organism[link[hit]]]
  name_conflicts <- gb[nc, c("accession", "organism", "bold_processid", "bold_species", "bin", "country")]
  bold_only <- bold[is.na(link), ]
  # BOLD process IDs beginning "GB" are records mined from GenBank; when they could not be
  # linked by accession they are still copies of GenBank records, so they are dropped here.
  gbm <- grepl("^GB", bold_only$processid)
  writeLines(bold_only$processid[gbm], file.path(out_dir, "bold_genbank_mined_unlinked_dropped.txt"))
  msg("BOLD-only records with a GenBank-mined process ID (GB...) dropped as duplicates:", sum(gbm))
  bold_only <- bold_only[!gbm, ]
  msg("BOLD-only records (no GenBank counterpart found):", nrow(bold_only))
} else bold_only <- NULL

rec <- data.frame(
  id = gb$accession, source = "GenBank", accession = gb$accession, processid = gb$bold_processid,
  species_raw = gb$organism, country = gb$country, bin = gb$bin, coi_source = gb$coi_source,
  seq = gb$seq, stringsAsFactors = FALSE)
if (!is.null(bold_only) && nrow(bold_only) > 0)
  rec <- rbind(rec, data.frame(
    id = bold_only$processid, source = "BOLD", accession = "", processid = bold_only$processid,
    species_raw = bold_only$species, country = bold_only$country, bin = bold_only$bin_uri,
    coi_source = "BOLD_COI-5P", seq = bold_only$nuc, stringsAsFactors = FALSE))
rec$id <- make.unique(rec$id)
rec$country[is.na(rec$country) | rec$country == ""] <- "Not recorded"
rec$species <- ifelse(grepl("^Cephalopholis [a-z]+$", rec$species_raw), rec$species_raw, NA)
msg("Records with an unresolved name (sp., cf., hybrids, blank) kept out of species statistics:",
    sum(is.na(rec$species)))

## ---- 2. reference COI and per-sequence alignment to it --------------------------
ri <- grep(paste0("^", P$ref_acc, "\\."), rec$accession)[1]
if (is.na(ri)) stop("Reference ", P$ref_acc, " not found in coi_genbank.tsv")
ref <- DNAString(rec$seq[ri]); L <- length(ref)
msg("Reference COI:", rec$accession[ri], "length", L, "bp; first codon", as.character(subseq(ref, 1, 3)))
smat <- nsm(match = 2, mismatch = -3, baseOnly = FALSE)

clean <- function(s) gsub("[^ACGTRYSWKMBDHVN]", "N", toupper(gsub("[-. ]", "", s)))
place <- function(s) {
  s <- DNAString(clean(s))
  a <- pa(ref, s, type = "local", substitutionMatrix = smat, gapOpening = 5, gapExtension = 2)
  b <- pa(ref, reverseComplement(s), type = "local", substitutionMatrix = smat,
          gapOpening = 5, gapExtension = 2)
  rc <- aln_score(b) > aln_score(a); if (rc) a <- b
  ap <- strsplit(as.character(aln_pat(a)), "")[[1]]
  as <- strsplit(as.character(aln_sub(a)), "")[[1]]
  ps <- aln_pstart(a)
  row <- rep("-", L)
  row[ps:(ps + sum(ap != "-") - 1)] <- as[ap != "-"]
  ins <- rle(ap == "-"); del <- rle(as == "-")
  ins_l <- ins$lengths[ins$values]; del_l <- del$lengths[del$values]
  list(row = row, pid = pid_fun(a, type = "PID1"), rc = rc, start = ps,
       end = ps + sum(ap != "-") - 1,
       frameshift = any(ins_l %% 3 != 0) || any(del_l %% 3 != 0),
       n_indel = length(ins_l) + length(del_l))
}
msg("Aligning", nrow(rec), "sequences to the reference COI ...")
pl <- lapply(rec$seq, place)
M  <- do.call(rbind, lapply(pl, `[[`, "row")); rownames(M) <- rec$id
rec$pid <- sapply(pl, `[[`, "pid"); rec$revcomp <- sapply(pl, `[[`, "rc")
rec$ref_start <- sapply(pl, `[[`, "start"); rec$ref_end <- sapply(pl, `[[`, "end")
rec$frameshift <- sapply(pl, `[[`, "frameshift"); rec$n_indel <- sapply(pl, `[[`, "n_indel")

## ---- 3. common barcode window (data-driven) -------------------------------------
cov <- integer(L)
for (k in which(rec$pid >= P$min_pid)) cov[rec$ref_start[k]:rec$ref_end[k]] <-
  cov[rec$ref_start[k]:rec$ref_end[k]] + 1L
if (!is.null(P$win_override)) { ws <- P$win_override[1]; we <- P$win_override[2] } else {
  r <- rle(cov >= P$win_frac * max(cov)); ends <- cumsum(r$lengths)
  k <- which(r$values)[which.max(r$lengths[r$values])]
  ws <- ends[k] - r$lengths[k] + 1; we <- ends[k]
}
ws <- ws + (3 - (ws - 1) %% 3) %% 3                        # start on codon position 1
we <- we - (we - ws + 1) %% 3                              # end on codon position 3
if (we == L) we <- L - 3                                   # never include the stop codon
msg(sprintf("Barcode window: KC593378 COI positions %d-%d (%d bp; max coverage %d sequences)",
            ws, we, we - ws + 1, max(cov)))
W <- M[, ws:we, drop = FALSE]
colnames(W) <- ws:we

## ---- 4. quality control ------------------------------------------------------------
rec$cov_window <- rowMeans(W != "-")
rec$ambig <- apply(W, 1, function(x) { y <- x[x != "-"]; if (!length(y)) 1 else mean(!y %in% c("A","C","G","T")) })
gc <- getGeneticCode("SGC1")                               # vertebrate mitochondrial code
rec$stops <- apply(W, 1, function(x) {
  x[!x %in% c("A","C","G","T")] <- "N"
  aa <- as.character(translate(DNAString(paste(x, collapse = "")), genetic.code = gc,
                               if.fuzzy.codon = "X"))
  sum(strsplit(aa, "")[[1]] == "*")
})
rec$qc <- with(rec, ifelse(pid < P$min_pid, "low_identity_to_reference",
              ifelse(cov_window < P$min_cov, "short_coverage_of_window",
              ifelse(ambig > P$max_amb, "too_many_ambiguous_bases",
              ifelse(stops > 0, "internal_stop_codon",
              ifelse(frameshift, "frameshift_indel", "pass"))))))
write.table(rec[, setdiff(names(rec), "seq")], file.path(out_dir, "coi_qc_all_records.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
msg("QC outcome:"); print(table(rec$qc, rec$source))
capture.output(table(rec$qc, rec$source), file = file.path(out_dir, "run_log.txt"), append = TRUE)

ok  <- rec$qc == "pass" & !is.na(rec$species)
R   <- rec[ok, ]; X <- W[ok, , drop = FALSE]
X[!X %in% c("A","C","G","T","-")] <- "N"
writeXStringSet(DNAStringSet(setNames(apply(X, 1, paste, collapse = ""),
                paste(R$id, gsub(" ", "_", R$species), gsub(" ", "_", R$country), sep = "|"))),
                file.path(out_dir, "coi_window_alignment_qc_pass.fasta"))
msg("Sequences passing QC with a resolved species name:", nrow(R), "in", length(unique(R$species)), "species")

## ---- 5. distances, barcode gap, nearest-neighbour screen --------------------------
db <- as.DNAbin(tolower(X))
D  <- as.matrix(dist.dna(db, model = "K80", pairwise.deletion = TRUE))
sp <- R$species
species_table <- function(keep) {
  sp <- sp[keep]; D <- D[keep, keep, drop = FALSE]; R <- R[keep, ]
  sp_tab <- do.call(rbind, lapply(sort(unique(sp)), function(s) {
    i <- which(sp == s); j <- which(sp != s)
    intra <- if (length(i) > 1) D[i, i][upper.tri(D[i, i])] else numeric(0)
    inter <- D[i, j, drop = FALSE]
    nn <- which(inter == min(inter, na.rm = TRUE), arr.ind = TRUE)[1, ]
    data.frame(species = s, n = length(i),
               n_countries = length(unique(R$country[i][R$country[i] != "Not recorded"])),
               n_GenBank = sum(R$source[i] == "GenBank"), n_BOLD_only = sum(R$source[i] == "BOLD"),
               intra_mean_pct = if (length(intra)) round(100 * mean(intra, na.rm = TRUE), 2) else NA,
               intra_max_pct  = if (length(intra)) round(100 * max(intra, na.rm = TRUE), 2) else NA,
               nn_species = sp[j][nn[2]], nn_min_pct = round(100 * min(inter, na.rm = TRUE), 2),
               stringsAsFactors = FALSE)
  }))
  sp_tab$barcode_gap_pct <- round(sp_tab$nn_min_pct - sp_tab$intra_max_pct, 2)
  sp_tab <- sp_tab[order(-sp_tab$n), ]
  sp_tab
}
sp_tab <- species_table(rep(TRUE, nrow(R)))
write.table(sp_tab, file.path(out_dir, "coi_species_distances.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

screen <- do.call(rbind, lapply(seq_len(nrow(R)), function(k) {
  d <- D[k, ]; d[k] <- NA
  con <- which(sp == sp[k]); con <- setdiff(con, k); het <- which(sp != sp[k])
  dc <- if (length(con)) min(d[con], na.rm = TRUE) else NA
  dh <- min(d[het], na.rm = TRUE); hn <- het[which.min(d[het])]
  reason <- c(if (!is.na(dc) && dh < dc) "closer_to_another_species",
              if (is.na(dc) && dh < P$thr) "singleton_within_threshold_of_another_species",
              if (!is.na(dc) && dc > P$thr) "far_from_all_conspecifics")
  data.frame(id = R$id[k], source = R$source[k], species = sp[k], country = R$country[k],
             bin = R$bin[k], min_conspecific_pct = round(100 * dc, 2),
             nearest_other_species = sp[hn], nearest_other_id = R$id[hn],
             min_other_species_pct = round(100 * dh, 2),
             flag = if (length(reason)) paste(reason, collapse = ";") else "", stringsAsFactors = FALSE)
}))
write.table(screen, file.path(out_dir, "coi_nearest_neighbour_screen_all.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
flagged <- screen[screen$flag != "", ]
write.table(flagged, file.path(out_dir, "coi_flagged_sequences.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
msg("Sequences flagged by the nearest-neighbour screen:", nrow(flagged),
    "(flags mark records to check; they are not proof of misidentification)")

# pairwise distance distributions (conspecific vs congeneric) for the gap histogram
ut <- which(upper.tri(D), arr.ind = TRUE)
pairs <- data.frame(dist = 100 * D[ut], type = ifelse(sp[ut[, 1]] == sp[ut[, 2]],
                    "Within species", "Between species"))
g1 <- ggplot(pairs, aes(dist, fill = type)) +
  geom_histogram(binwidth = 0.25, position = "identity", alpha = 0.65, colour = NA) +
  scale_fill_manual(values = c("Within species" = "#1F618D", "Between species" = "#C0392B"), name = NULL) +
  scale_y_sqrt() + labs(x = "K2P distance (%)", y = "Number of pairs (square-root scale)") +
  theme_bw(base_size = 9) + theme(legend.position = "top")
ggsave(file.path(out_dir, "Fig_COI_distance_histogram.pdf"), g1, width = 6.5, height = 3.8)
ggsave(file.path(out_dir, "Fig_COI_distance_histogram.png"), g1, width = 6.5, height = 3.8, dpi = 300)

gp <- sp_tab[!is.na(sp_tab$intra_max_pct), ]
gp$focal <- ifelse(gp$species %in% c(P$focal, P$sister), gp$species, "Other")
g2 <- ggplot(gp, aes(intra_max_pct, nn_min_pct)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(colour = focal), size = 2) +
  geom_text(aes(label = sub("^Cephalopholis ", "C. ", species)), fontface = "italic",
            size = 2.4, vjust = -0.8, check_overlap = TRUE) +
  scale_colour_manual(values = setNames(c("#C0392B", "#7D3C98", "grey40"), c(P$focal, P$sister, "Other")),
                      name = NULL) +
  labs(x = "Maximum within-species K2P distance (%)", y = "Distance to nearest other species (%)") +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.08)) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 9) + theme(legend.position = "none", plot.margin = margin(8, 18, 6, 6))
ggsave(file.path(out_dir, "Fig_COI_barcode_gap_all_records.pdf"), g2, width = 5.5, height = 5)
ggsave(file.path(out_dir, "Fig_COI_barcode_gap_all_records.png"), g2, width = 5.5, height = 5, dpi = 300)

## ---- 6. BIN concordance (BOLD) --------------------------------------------------------
bins <- R[R$bin != "" & !is.na(R$bin), c("species", "bin")]
if (nrow(bins)) {
  b1 <- aggregate(bin ~ species, bins, function(x) paste(sort(unique(x)), collapse = ";"))
  b1$n_bins <- sapply(strsplit(b1$bin, ";"), length)
  b2 <- aggregate(species ~ bin, bins, function(x) paste(sort(unique(x)), collapse = ";"))
  b2$n_species <- sapply(strsplit(b2$species, ";"), length)
  write.table(b1, file.path(out_dir, "coi_bins_per_species.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(b2, file.path(out_dir, "coi_species_per_bin.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  msg("BINs:", nrow(b2), "| BINs shared by >1 species name:", sum(b2$n_species > 1),
      "| species split across >1 BIN:", sum(b1$n_bins > 1))
} else msg("No BIN assignments available (BOLD data missing or unlinked).")
if (!is.null(name_conflicts) && nrow(name_conflicts))
  write.table(name_conflicts, file.path(out_dir, "coi_genbank_bold_name_conflicts.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
msg("GenBank records whose linked BOLD record carries a different species name:",
    if (is.null(name_conflicts)) 0 else nrow(name_conflicts))

## sequences whose nearest neighbour is another species are left out of the
## network and of the second diagnostic-site run (listed in run_log.txt)
excl <- screen$id[grepl("closer_to_another_species|singleton_within", screen$flag)]
msg("Excluded from network and from the 'excluding flagged' diagnostic run:",
    if (length(excl)) paste(excl, collapse = ", ") else "none")
write.table(species_table(!R$id %in% excl), file.path(out_dir, "coi_species_distances_excluding_flagged.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## C. sonnerati records within 1% K2P of a record of another named species (that record
## not itself flagged) are set aside as probable misidentifications pending voucher checks.
P$misid_thr <- 0.01
fi <- which(sp == P$focal & !R$id %in% excl)
oi <- which(sp != P$focal & !R$id %in% excl)
dmin <- apply(D[fi, oi, drop = FALSE], 1, min, na.rm = TRUE)
nni  <- oi[apply(D[fi, oi, drop = FALSE], 1, which.min)]
mis  <- data.frame(id = R$id[fi], country = R$country[fi], bin = R$bin[fi],
                   nearest_other_species = sp[nni], nearest_other_id = R$id[nni],
                   min_other_species_pct = round(100 * dmin, 2), stringsAsFactors = FALSE)
mis  <- mis[dmin < P$misid_thr, ]
write.table(mis, file.path(out_dir, "csonnerati_probable_misidentifications.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
msg(sprintf("C. sonnerati records within %.0f%% of another named species (set aside): %d", 100 * P$misid_thr, nrow(mis)),
    if (nrow(mis)) paste0("(", paste(mis$id, collapse = ", "), ")") else "")
excl2 <- union(excl, mis$id)
sp_cur <- species_table(!R$id %in% excl2)
write.table(sp_cur, file.path(out_dir, "coi_species_distances_curated.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
gpc <- sp_cur[!is.na(sp_cur$intra_max_pct), ]
gpc$focal <- ifelse(gpc$species %in% c(P$focal, P$sister), gpc$species, "Other")
g3 <- g2 %+% gpc + geom_text(aes(label = sub("^Cephalopholis ", "C. ", species)), fontface = "italic",
                             size = 2.4, vjust = -0.8)
g3$layers[[3]] <- NULL                                       # drop the overlap-thinned label layer
ggsave(file.path(out_dir, "Fig_COI_barcode_gap_curated.pdf"), g3, width = 5.5, height = 5)
ggsave(file.path(out_dir, "Fig_COI_barcode_gap_curated.png"), g3, width = 5.5, height = 5, dpi = 300)

## ---- 7. C. sonnerati haplotype network and diversity ------------------------------------
fs <- which(sp == P$focal & !R$id %in% excl2)
Xs <- X[fs, , drop = FALSE]
okcol <- colMeans(matrix(Xs %in% c("A","C","G","T"), nrow = nrow(Xs))) >= P$core_frac
rr <- rle(okcol); en <- cumsum(rr$lengths); kk <- which(rr$values)[which.max(rr$lengths[rr$values])]
core <- (en[kk] - rr$lengths[kk] + 1):en[kk]
Xc <- Xs[, core, drop = FALSE]
full <- apply(Xc, 1, function(x) all(x %in% c("A","C","G","T")))
msg(sprintf("C. sonnerati: %d QC-pass sequences; network core = window columns %d-%d (KC593378 COI %s-%s, %d bp); %d complete over the core, %d excluded",
            length(fs), min(core), max(core), colnames(X)[min(core)], colnames(X)[max(core)],
            length(core), sum(full), sum(!full)))
Xc <- Xc[full, , drop = FALSE]; Rs <- R[fs[full], ]
hs <- NULL
if (nrow(Xc) >= 2) {
  dbs <- as.DNAbin(tolower(Xc))
  hs  <- haplotype(dbs)
  idx <- attr(hs, "index")
  hap_tab <- do.call(rbind, lapply(seq_along(idx), function(h) {
    i <- idx[[h]]
    data.frame(haplotype = rownames(hs)[h], n = length(i),
               countries = paste(names(table(Rs$country[i])), table(Rs$country[i]), sep = ":", collapse = "; "),
               bins = paste(sort(unique(Rs$bin[i][Rs$bin[i] != ""])), collapse = ";"),
               records = paste(Rs$id[i], collapse = ";"), stringsAsFactors = FALSE)
  }))
  write.table(hap_tab, file.path(out_dir, "csonnerati_haplotypes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  div <- function(x, lab) data.frame(group = lab, n = nrow(x),
    n_haplotypes = nrow(haplotype(x)), segregating_sites = length(seg.sites(x)),
    hap_diversity = round(hap.div(x), 4), nuc_diversity = round(nuc.div(x), 5))
  dv <- div(dbs, "All C. sonnerati")
  for (cc in names(which(table(Rs$country) >= 3)))
    dv <- rbind(dv, div(dbs[Rs$country == cc, ], cc))
  write.table(dv, file.path(out_dir, "csonnerati_diversity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  print(dv)
  if (nrow(hs) >= 2) {
    net <- haploNet(hs)
    hf  <- haploFreq(dbs, fac = factor(Rs$country), haplo = hs)
    cols <- setNames(hcl.colors(ncol(hf), "Dark 3"), colnames(hf))
    if ("Not recorded" %in% names(cols)) cols["Not recorded"] <- "grey85"
    for (ext in c("pdf", "png")) {
      if (ext == "pdf") pdf(file.path(out_dir, "Fig_Csonnerati_haplotype_network.pdf"), width = 8, height = 7) else
        png(file.path(out_dir, "Fig_Csonnerati_haplotype_network.png"), width = 8, height = 7, units = "in", res = 300)
      par(mar = c(1, 1, 2, 1))
      plot(net, size = sqrt(attr(net, "freq")), pie = hf, bg = cols[colnames(hf)],
           scale.ratio = 2, cex = 0.6, labels = FALSE, show.mutation = 2, legend = FALSE,
           threshold = 0)                                   # hide alternative links
      legend("topright", legend = sprintf("%s (n = %d)", colnames(hf), colSums(hf)), fill = cols[colnames(hf)],
             bty = "n", cex = 0.7, title = "Country (as recorded)")
      title(main = substitute(italic(x) ~ "COI haplotype network (" * n ~ "sequences," ~ L ~ "bp)",
                              list(x = P$focal, n = nrow(Xc), L = length(core))), cex.main = 0.9)
      dev.off()
    }
    bf  <- ifelse(is.na(Rs$bin) | Rs$bin == "", "No BIN", Rs$bin)
    hb  <- haploFreq(dbs, fac = factor(bf), haplo = hs)
    cb  <- setNames(hcl.colors(ncol(hb), "Set 2"), colnames(hb))
    if ("No BIN" %in% names(cb)) cb["No BIN"] <- "grey85"
    for (ext in c("pdf", "png")) {
      if (ext == "pdf") pdf(file.path(out_dir, "Fig_Csonnerati_haplotype_network_BIN.pdf"), width = 8, height = 7) else
        png(file.path(out_dir, "Fig_Csonnerati_haplotype_network_BIN.png"), width = 8, height = 7, units = "in", res = 300)
      par(mar = c(1, 1, 2, 1))
      plot(net, size = sqrt(attr(net, "freq")), pie = hb, bg = cb[colnames(hb)],
           scale.ratio = 2, cex = 0.6, labels = TRUE, show.mutation = 2, legend = FALSE, threshold = 0)
      legend("topright", legend = sprintf("%s (n = %d)", colnames(hb), colSums(hb)), fill = cb[colnames(hb)],
             bty = "n", cex = 0.7, title = "BOLD BIN")
      title(main = substitute(italic(x) ~ "COI haplotypes by BOLD BIN (" * n ~ "sequences," ~ L ~ "bp)",
                              list(x = P$focal, n = nrow(Xc), L = length(core))), cex.main = 0.9)
      dev.off()
    }
  }
}

## ---- 8. fixed differences: C. sonnerati vs C. oligosticta ---------------------------------
## Run twice: with all QC-pass records, and excluding the flagged records above.
codon_aa <- function(cod) if (all(cod %in% c("A","C","G","T")))
  as.character(translate(DNAString(paste(cod, collapse = "")), genetic.code = gc, no.init.codon = TRUE)) else "?"
majority <- function(m, j) { b <- m[, j]; b <- b[b %in% c("A","C","G","T")]; if (length(b)) names(which.max(table(b))) else "N" }
diag_run <- function(drop, tag) {
  fs <- which(sp == P$focal & !R$id %in% drop); fo <- which(sp == P$sister & !R$id %in% drop)
  present <- sapply(P$sister_refs, function(a) any(grepl(paste0("^", a, "\\."), R$accession[fo])))
  msg(sprintf("[%s] C. sonnerati n = %d, C. oligosticta n = %d | named references present & passing QC: %s",
      tag, length(fs), length(fo), paste(names(present), present, sep = "=", collapse = ", ")))
  if (!length(fo) || !length(fs)) { msg("  skipped: one species has no sequence"); return(invisible()) }
  S <- X[fs, , drop = FALSE]; O <- X[fo, , drop = FALSE]
  Oth <- X[sp != P$focal & !R$id %in% drop, , drop = FALSE]
  diag <- do.call(rbind, lapply(seq_len(ncol(X)), function(j) {
    bs <- S[, j][S[, j] %in% c("A","C","G","T")]; bo <- O[, j][O[, j] %in% c("A","C","G","T")]
    if (!length(bs) || !length(bo) || length(intersect(bs, bo))) return(NULL)
    pos <- as.integer(colnames(X)[j]); cp <- (pos - 1) %% 3 + 1
    cj  <- j - cp + 1 + 0:2
    if (min(cj) < 1 || max(cj) > ncol(X)) aaS <- aaO <- "?" else {
      aaS <- codon_aa(sapply(cj, function(c) majority(S, c)))
      aaO <- codon_aa(sapply(cj, function(c) majority(O, c)))
    }
    bt <- Oth[, j][Oth[, j] %in% c("A","C","G","T")]
    data.frame(KC593378_COI_position = pos, codon_position = cp,
               sonnerati_bases = paste(names(table(bs)), table(bs), sep = ":", collapse = ","),
               oligosticta_bases = paste(names(table(bo)), table(bo), sep = ":", collapse = ","),
               n_sonnerati = length(bs), n_oligosticta = length(bo),
               aa_sonnerati = aaS, aa_oligosticta = aaO,
               change = ifelse(aaS == "?" | aaO == "?", "undetermined",
                               ifelse(aaS == aaO, "synonymous", "nonsynonymous")),
               sonnerati_base_absent_from_all_other_species = !any(unique(bs) %in% bt),
               stringsAsFactors = FALSE)
  }))
  nd <- if (is.null(diag)) 0 else nrow(diag)
  if (is.null(diag)) diag <- data.frame(note = "no fixed differences found")
  write.table(diag, file.path(out_dir, sprintf("diagnostic_sites_sonnerati_vs_oligosticta_%s.tsv", tag)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  so <- D[fs, fo, drop = FALSE]
  msg(sprintf("  fixed differences: %d (C. sonnerati base absent from all other species at %d); K2P between species: min %.2f%%, mean %.2f%%, max %.2f%%",
      nd, if (nd) sum(diag$sonnerati_base_absent_from_all_other_species) else 0,
      max(0, 100 * min(so, na.rm = TRUE)), 100 * mean(so, na.rm = TRUE), 100 * max(so, na.rm = TRUE)))
}
diag_run(character(0), "all_qc_pass")
diag_run(excl, "excluding_flagged")
diag_run(excl2, "curated")

## ---- 9. NJ tree of unique haplotypes per species ----------------------------------------------
key <- paste(sp, apply(X, 1, paste, collapse = ""))
u   <- !duplicated(key)
cnt <- table(key)[key[u]]
Xu  <- X[u, , drop = FALSE]
lab <- sprintf("%s | %s | n=%d", sub("^Cephalopholis ", "C. ", sp[u]), R$id[u], as.integer(cnt))
rownames(Xu) <- make.unique(gsub("[ ():,;']", "_", lab))
dbu <- as.DNAbin(tolower(Xu))
treefun <- function(x) { d <- dist.dna(x, model = "K80", pairwise.deletion = TRUE)
  if (any(is.na(d))) njs(d) else bionj(d) }
tr <- treefun(dbu)
bs <- boot.phylo(tr, dbu, treefun, B = P$n_boot, quiet = TRUE, mc.cores = 1)
tr$node.label <- round(100 * bs / P$n_boot)
if (requireNamespace("phangorn", quietly = TRUE)) tr <- phangorn::midpoint(tr, node.labels = "support")
write.tree(tr, file.path(out_dir, "coi_NJ_unique_haplotypes.nwk"))
tipsp <- sub(" \\|.*", "", tr$tip.label)
tcol <- ifelse(grepl("^C._sonnerati", tr$tip.label), "#C0392B",
        ifelse(grepl("^C._oligosticta", tr$tip.label), "#7D3C98", "grey25"))
h <- max(8, 0.11 * Ntip(tr))
pdf(file.path(out_dir, "Fig_COI_NJ_tree.pdf"), width = 8, height = h)
par(mar = c(1, 0.5, 1, 0.5))
plot(tr, cex = 0.45, tip.color = tcol, label.offset = 0.001, no.margin = FALSE)
nl <- suppressWarnings(as.numeric(tr$node.label)); nl[is.na(nl) | nl < 70] <- NA
nodelabels(ifelse(is.na(nl), "", nl), frame = "none", cex = 0.4, adj = c(1.2, -0.4))
add.scale.bar(cex = 0.6); dev.off()
msg("NJ tree:", Ntip(tr), "unique haplotypes;", P$n_boot, "bootstrap replicates")

## ---- 10. parameters, summary and session --------------------------------------------------------
writeLines(c(sprintf("%s = %s", names(P), sapply(P, paste, collapse = ",")),
             sprintf("window = %d-%d", ws, we),
             sprintf("distance = K2P (ape::dist.dna, model K80, pairwise deletion)"),
             sprintf("genetic code = vertebrate mitochondrial (Biostrings SGC1)"),
             sprintf("run date = %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z"))),
           file.path(out_dir, "parameters.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
msg("\nDone. Outputs in", out_dir)
