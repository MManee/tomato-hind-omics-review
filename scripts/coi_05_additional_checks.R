# =============================================================================
# COI barcode analysis, step 5: additional checks of distances, diversity,
# synonymy and reference provenance.
# Reads outputs of coi_03_analysis.R (no download, no re-alignment) and writes
# results/coi_additional_checks/:
#
#   1  extreme_distances.tsv      for every species in the curated set: the pair of
#                                 records giving the maximum within-species distance
#                                 and the pair giving the minimum distance to another
#                                 species, with countries, BINs and overlap
#      outlier_records.tsv        records whose median distance to conspecifics is > 2%
#   2  diversity_102.tsv          haplotypes, segregating sites, Hd and pi for the 104
#                                 network sequences and for 102 after removing one
#                                 record of each same-voucher pair
#   3  synonymy_check.tsv         for each of the 15 candidate positions, the amino acid
#                                 coded by every sequence's own codon in each species
#   4  oligosticta_provenance.tsv GenBank/BOLD metadata of KM656840 and MH707291
#   5  additional_checks_summary.txt         headline results
#
# Run (after coi_03_analysis.R and coi_04_robustness_checks.R):
#   cd tomato-hind-omics-review   # repository root
#   Rscript scripts/coi_05_additional_checks.R
# =============================================================================
suppressPackageStartupMessages({ library(ape); library(Biostrings); library(pegas) })

proj    <- Sys.getenv("PROJ", ".")
in_dir  <- file.path(proj, "data/coi_20261008")
res_dir <- file.path(proj, "results/coi_analysis")
out_dir <- file.path(proj, "results/coi_additional_checks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sumf <- file.path(out_dir, "additional_checks_summary.txt"); unlink(sumf)
say <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
wt  <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE, row.names = FALSE)
need <- function(f) { if (!file.exists(f)) stop("Missing input: ", f); f }
focal <- "Cephalopholis sonnerati"; sister <- "Cephalopholis oligosticta"
WIN_START <- 94                                         # window column 1 = KC593378 COI position 94
say("Run:", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), "|", R.version.string)

## ---- inputs ---------------------------------------------------------------------
qc   <- read.delim(need(file.path(res_dir, "coi_qc_all_records.tsv")), colClasses = "character", quote = "")
scr  <- read.delim(need(file.path(res_dir, "coi_nearest_neighbour_screen_all.tsv")), colClasses = "character", quote = "")
mis  <- read.delim(need(file.path(res_dir, "csonnerati_probable_misidentifications.tsv")), colClasses = "character", quote = "")
gb   <- read.delim(need(file.path(in_dir, "coi_genbank.tsv")), colClasses = "character", quote = "")
dg   <- read.delim(need(file.path(res_dir, "diagnostic_sites_sonnerati_vs_oligosticta_curated.tsv")), colClasses = "character", quote = "")
fa   <- readDNAStringSet(need(file.path(res_dir, "coi_window_alignment_qc_pass.fasta")))
hdr  <- strsplit(names(fa), "|", fixed = TRUE)
ids  <- sapply(hdr, `[`, 1); spp <- gsub("_", " ", sapply(hdr, `[`, 2))
X    <- do.call(rbind, strsplit(as.character(fa), "")); rownames(X) <- ids
excl  <- scr$id[grepl("closer_to_another_species|singleton_within", scr$flag)]
keep  <- !ids %in% union(excl, mis$id)                    # curated set (1247 records)
meta  <- qc[match(ids, qc$id), c("id", "country", "bin")]
acgt  <- matrix(X %in% c("A","C","G","T"), nrow = nrow(X))
D     <- as.matrix(dist.dna(as.DNAbin(tolower(X)), model = "K80", pairwise.deletion = TRUE))
ovl   <- function(i, j) sum(acgt[i, ] & acgt[j, ])
say("Curated records:", sum(keep), "in", length(unique(spp[keep])), "species")

## ---- 1. extreme distances ----------------------------------------------------
K <- which(keep); rows <- list(); outl <- list()
for (s in sort(unique(spp[K]))) {
  i <- K[spp[K] == s]; j <- K[spp[K] != s]
  r <- data.frame(species = s, n = length(i), stringsAsFactors = FALSE)
  if (length(i) > 1) {
    Di <- D[i, i]; Di[lower.tri(Di, diag = TRUE)] <- NA
    w <- which(Di == max(Di, na.rm = TRUE), arr.ind = TRUE)[1, ]
    a <- i[w[1]]; b <- i[w[2]]
    r <- cbind(r, max_within_pct = round(100 * D[a, b], 2), max_pair = paste(ids[a], ids[b], sep = " / "),
               max_pair_countries = paste(meta$country[a], meta$country[b], sep = " / "),
               max_pair_bins = paste(meta$bin[a], meta$bin[b], sep = " / "), max_pair_overlap_bp = ovl(a, b))
    med <- sapply(i, function(k) median(D[k, setdiff(i, k)], na.rm = TRUE))
    o <- i[med > 0.02]
    if (length(o)) outl[[s]] <- data.frame(species = s, id = ids[o], country = meta$country[o], bin = meta$bin[o],
                                           median_to_conspecifics_pct = round(100 * med[med > 0.02], 2),
                                           nearest_other_species = sapply(o, function(k) spp[j][which.min(D[k, j])]),
                                           nearest_other_pct = round(100 * sapply(o, function(k) min(D[k, j], na.rm = TRUE)), 2))
  } else r <- cbind(r, max_within_pct = NA, max_pair = "", max_pair_countries = "", max_pair_bins = "", max_pair_overlap_bp = NA)
  Dij <- D[i, j, drop = FALSE]; w <- which(Dij == min(Dij, na.rm = TRUE), arr.ind = TRUE)[1, ]
  a <- i[w[1]]; b <- j[w[2]]
  rows[[s]] <- cbind(r, min_between_pct = round(100 * D[a, b], 2), min_pair = paste(ids[a], ids[b], sep = " / "),
                     min_pair_species = spp[b], min_pair_overlap_bp = ovl(a, b))
}
ex <- do.call(rbind, rows); ex <- ex[order(-ex$max_within_pct), ]
wt(ex, "extreme_distances.tsv")
ol <- if (length(outl)) do.call(rbind, outl) else data.frame()
wt(ol, "outlier_records.tsv")
say("\n== 1. Records behind the largest within-species distances (curated set) ==")
for (k in seq_len(nrow(ex))) if (!is.na(ex$max_within_pct[k]) && ex$max_within_pct[k] > 2)
  say(sprintf("  %-28s max %6.2f%%: %s | %s | %s | overlap %d bp",
              ex$species[k], ex$max_within_pct[k], ex$max_pair[k], ex$max_pair_countries[k],
              ex$max_pair_bins[k], ex$max_pair_overlap_bp[k]))
say("  Minimum pairwise overlap among all curated pairs:",
    min(sapply(seq_len(length(K) - 1), function(a) min(colSums(t(acgt[K[(a+1):length(K)], , drop = FALSE]) & acgt[K[a], ])))), "bp")
say("  Records whose median distance to conspecifics exceeds 2%:", nrow(ol))
if (nrow(ol)) for (k in seq_len(nrow(ol))) say(sprintf("    %s %s (%s; %s): median %.2f%% to conspecifics; nearest other species %s at %.2f%%",
  ol$species[k], ol$id[k], ol$country[k], ol$bin[k], ol$median_to_conspecifics_pct[k], ol$nearest_other_species[k], ol$nearest_other_pct[k]))

## ---- 2. diversity for 104 and 102 sequences -------------------------------------
fs <- which(spp == focal & keep)
Xs <- X[fs, , drop = FALSE]
okcol <- colMeans(matrix(Xs %in% c("A","C","G","T"), nrow = nrow(Xs))) >= 0.95
rr <- rle(okcol); en <- cumsum(rr$lengths); kk <- which(rr$values)[which.max(rr$lengths[rr$values])]
core <- (en[kk] - rr$lengths[kk] + 1):en[kk]
full <- apply(Xs[, core, drop = FALSE], 1, function(x) all(x %in% c("A","C","G","T")))
Xc <- Xs[full, core, drop = FALSE]
div <- function(m, lab) { d <- as.DNAbin(tolower(m))
  data.frame(set = lab, n = nrow(m), haplotypes = nrow(haplotype(d)), segregating_sites = length(seg.sites(d)),
             hap_diversity = round(hap.div(d), 4), nuc_diversity = round(nuc.div(d), 5)) }
drop2 <- c("PP680539.1", "PP680540.1")                     # one record of each same-voucher pair
dv <- rbind(div(Xc, "network (complete over core)"), div(Xc[!rownames(Xc) %in% drop2, ], "same-voucher duplicates removed"))
wt(dv, "diversity_102.tsv")
say("\n== 2. Diversity over the", length(core), "bp core ==")
for (k in seq_len(nrow(dv))) say(sprintf("  %-34s n = %d; haplotypes %d; S = %d; Hd = %.4f; pi = %.5f",
  dv$set[k], dv$n[k], dv$haplotypes[k], dv$segregating_sites[k], dv$hap_diversity[k], dv$nuc_diversity[k]))
say("  Removed for the second row:", paste(intersect(drop2, rownames(Xc)), collapse = ", "))

## ---- 3. synonymy at the 15 candidate positions, sequence by sequence --------------
gc <- getGeneticCode("SGC1")
fo <- which(spp == sister & keep)
syn <- do.call(rbind, lapply(seq_len(nrow(dg)), function(k) {
  pos <- as.integer(dg$KC593378_COI_position[k]); cp <- as.integer(dg$codon_position[k])
  cols <- (pos - cp + 1):(pos - cp + 3) - WIN_START + 1
  aa <- function(rws) { cod <- apply(X[rws, cols, drop = FALSE], 1, paste, collapse = "")
    cod <- cod[!grepl("[^ACGT]", cod)]
    tab <- table(sapply(cod, function(c) as.character(translate(DNAString(c), genetic.code = gc, no.init.codon = TRUE))))
    list(txt = paste(names(tab), tab, sep = ":", collapse = ","), set = names(tab),
         codons = paste(names(table(cod)), table(cod), sep = ":", collapse = ",")) }
  s <- aa(fs); o <- aa(fo)
  data.frame(position = pos, codon_position = cp, sonnerati_codons = s$codons, sonnerati_aa = s$txt,
             oligosticta_codons = o$codons, oligosticta_aa = o$txt,
             all_synonymous = length(union(s$set, o$set)) == 1, stringsAsFactors = FALSE)
}))
wt(syn, "synonymy_check.tsv")
say("\n== 3. Amino acids coded by each sequence's own codon at the 15 positions ==")
for (k in seq_len(nrow(syn))) say(sprintf("  %3d (codon pos %d): sonnerati %s | oligosticta %s | synonymous for all sequences: %s",
  syn$position[k], syn$codon_position[k], syn$sonnerati_aa[k], syn$oligosticta_aa[k], syn$all_synonymous[k]))
say("  Positions not synonymous for every sequence:", sum(!syn$all_synonymous))

## ---- 4. provenance of the two C. oligosticta COI records ---------------------------
pv <- gb[sub("\\.\\d+$", "", gb$accession) %in% c("KM656840", "MH707291"),
         intersect(c("accession", "organism", "country", "locality", "lat_lon", "specimen_voucher", "isolate",
                     "collection_date", "bold_xref", "first_reference", "definition", "record_length"), names(gb))]
wt(pv, "oligosticta_provenance.tsv")
say("\n== 4. C. oligosticta COI records ==")
for (k in seq_len(nrow(pv))) say(paste0("  ", paste(names(pv), pv[k, ], sep = " = ", collapse = " | ")))
if (nrow(pv) == 2) say("  Same voucher:", identical(pv$specimen_voucher[1], pv$specimen_voucher[2]) && pv$specimen_voucher[1] != "",
  "| K2P between them:", round(100 * D[grep("^KM656840", ids), grep("^MH707291", ids)], 2), "%",
  "| overlap:", ovl(grep("^KM656840", ids), grep("^MH707291", ids)), "bp")
say("\nDone. Summary written to:", sumf)
