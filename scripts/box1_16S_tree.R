# =============================================================================
# Box 1 figure: 16S rRNA maximum-likelihood tree placing the Saudi
# "Cephalopholis sonnerati" records (KX298697-KX298699) among Cephalopholis
# references from GenBank.
#
# Input : GenBank files already downloaded by verify_saudi_records.py
# Output: results/box1_tree/  (PDF + PNG figure, Newick tree, alignment,
#         distance table, metadata)
#
# One-time setup in R:
#   install.packages(c("ape", "phangorn", "ggplot2", "BiocManager"))
#   BiocManager::install(c("Biostrings", "pwalign", "DECIPHER", "ggtree"))
#
# Run from Terminal:
#   cd tomato-hind-omics-review   # repository root
#   Rscript scripts/box1_16S_tree.R
# =============================================================================

suppressPackageStartupMessages({
  library(ape); library(phangorn); library(Biostrings); library(DECIPHER)
  library(ggplot2); library(ggtree)
})

## ---- paths & settings -------------------------------------------------------
proj    <- Sys.getenv("PROJ", ".")   # repository root
in_dir  <- file.path(proj, "data/ncbi_inventory_20261001/verify_saudi")
out_dir <- file.path(proj, "results/box1_tree")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gb_files <- file.path(in_dir, c("saudi_queries.gb", "ref_cephalopholis_16S_12S.gb"))
queries  <- c("KX298697.1", "KX298698.1", "KX298699.1")   # Saudi "C. sonnerati" 16S
anchor   <- "KX298698.1"                                   # defines trimmed region (cleanest 5' end)
same_study_prefix <- "^KX2986|^KX2987"                     # unpublished 2016 Saudi study
min_cov  <- 0.70     # keep references covering >= 70% of the trimmed region
n_boot   <- 1000     # bootstrap replicates (100 is enough for a quick test)
set.seed(2026)

## ---- 1. parse GenBank flat files ---------------------------------------------
parse_gb <- function(file) {
  txt <- readLines(file, warn = FALSE)
  st  <- grep("^LOCUS", txt); en <- grep("^//", txt)
  do.call(rbind, lapply(seq_along(st), function(i) {
    r <- txt[st[i]:en[i]]; j <- paste(r, collapse = "\n")
    acc <- sub("^VERSION\\s+(\\S+).*", "\\1", grep("^VERSION", r, value = TRUE)[1])
    d1  <- grep("^DEFINITION", r)[1]; k <- d1 + 1
    while (k <= length(r) && grepl("^\\s{12}", r[k])) k <- k + 1
    def <- paste(trimws(sub("^DEFINITION", "", r[d1:(k - 1)])), collapse = " ")
    org <- trimws(sub("^\\s+ORGANISM\\s+", "", grep("^\\s+ORGANISM", r, value = TRUE)[1]))
    m   <- regmatches(j, regexpr('/(geo_loc_name|country)="[^"]*"', j))
    ctry <- if (length(m)) sub('.*="([^":]*).*', "\\1", m) else ""
    o   <- grep("^ORIGIN", r)[1]
    sq  <- toupper(gsub("[^A-Za-z]", "", paste(r[(o + 1):(length(r) - 1)], collapse = "")))
    data.frame(acc = acc, organism = org, country = ctry, definition = def,
               seq = sq, stringsAsFactors = FALSE)
  }))
}
gb <- unique(do.call(rbind, lapply(gb_files, parse_gb)))
gb <- gb[!duplicated(gb$acc), ]

## ---- 2. keep 16S records only -------------------------------------------------
is16 <- grepl("16S", gb$definition) | gb$acc %in% queries
d <- gb[is16 & nchar(gb$seq) >= 300, ]
d$group <- ifelse(d$acc %in% queries, "Saudi 'C. sonnerati' (query)",
           ifelse(grepl(same_study_prefix, d$acc), "Same unpublished Saudi study",
           ifelse(d$organism == "Cephalopholis sonnerati", "C. sonnerati (reference)",
           ifelse(d$organism == "Cephalopholis oligosticta", "C. oligosticta (independent)",
                  "Other Cephalopholis"))))
cat(nrow(d), "16S sequences;", sum(d$acc %in% queries), "queries\n")

## ---- 3. orient, then trim EACH sequence to the anchor region by pairwise -----
##        local alignment (robust to fragments from different parts of 16S)
dna <- DNAStringSet(setNames(d$seq, d$acc))
dna <- OrientNucleotides(dna, verbose = FALSE)
anc <- dna[[anchor]]
smat <- pwalign::nucleotideSubstitutionMatrix(match = 2, mismatch = -3, baseOnly = FALSE)
trim <- do.call(rbind, lapply(names(dna), function(id) {
  a  <- pwalign::pairwiseAlignment(anc, dna[[id]], type = "local",
                                   substitutionMatrix = smat, gapOpening = 5, gapExtension = 2)
  ps <- start(pwalign::pattern(a)); pe <- end(pwalign::pattern(a))
  ss <- start(pwalign::subject(a)); se <- end(pwalign::subject(a))
  data.frame(acc = id, cov = (pe - ps + 1) / length(anc),
             pid = pwalign::pid(a, type = "PID1"),
             seg = if (se >= ss) as.character(subseq(dna[[id]], ss, se)) else "",
             stringsAsFactors = FALSE)
}))
keep <- (trim$cov >= min_cov & trim$pid >= 85) | trim$acc %in% queries
cat("Coverage of anchor region (all refs): median",
    round(median(trim$cov[!trim$acc %in% queries]) * 100), "%\n")
cat("Kept", sum(keep), "of", nrow(trim), "sequences (>=", min_cov * 100,
    "% of anchor region, >= 85% identity)\n")
trim <- trim[keep, ]
if (sum(!trim$acc %in% queries) < 3)
  stop("Too few references passed the coverage/identity filter; lower min_cov and re-run.")

## ---- 4. multiple alignment of the trimmed segments (DECIPHER) ----------------
seg <- DNAStringSet(setNames(trim$seg, trim$acc))
for (q in queries) seg[[q]] <- dna[[q]]        # queries kept full-length
aln <- AlignSeqs(seg, verbose = FALSE, processors = NULL)
m   <- as.matrix(aln)
m   <- m[, colSums(m != "-") > 0, drop = FALSE]          # remove all-gap columns
d   <- d[match(rownames(m), d$acc), ]
writeXStringSet(DNAStringSet(apply(m, 1, paste, collapse = "")),
                file.path(out_dir, "box1_16S_alignment_trimmed.fasta"))
cat("Trimmed alignment:", nrow(m), "sequences x", ncol(m), "columns\n")

## ---- 5. distances from each query to each species (p-distance, pairwise del.) -
db <- as.DNAbin(tolower(m))
pd <- as.matrix(dist.dna(db, model = "raw", pairwise.deletion = TRUE)) * 100
dist_tab <- do.call(rbind, lapply(queries[queries %in% rownames(pd)], function(q) {
  others <- setdiff(rownames(pd), c(queries, grep(same_study_prefix, rownames(pd), value = TRUE)))
  sp <- d$organism[match(others, d$acc)]
  agg <- aggregate(pd[q, others], by = list(species = sp), FUN = min)
  names(agg)[2] <- "min_p_distance_pct"
  agg$closest_ref <- sapply(agg$species, function(s)
    others[sp == s][which.min(pd[q, others[sp == s]])])
  cbind(query = q, agg[order(agg$min_p_distance_pct), ])
}))
write.table(dist_tab, file.path(out_dir, "box1_query_distances.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nClosest species to each Saudi query (independent references only):\n")
print(do.call(rbind, lapply(split(dist_tab, dist_tab$query), head, 4)), row.names = FALSE)

## ---- 6. ML tree: model selection + NNI search + bootstrap (phangorn) ---------
##        Set REPLOT=1 to skip this step and only redraw from the saved tree.
tree_file <- file.path(out_dir, "box1_16S_ML_tree.nwk")
if (Sys.getenv("REPLOT") == "1" && file.exists(tree_file)) {
  tr <- read.tree(tree_file); best <- "(read from saved tree)"
  cat("\nREPLOT=1: re-using", tree_file, "\n")
} else {
dat <- phyDat(db, type = "DNA")
mt  <- modelTest(dat, model = c("JC", "K80", "HKY", "GTR"), G = TRUE, I = TRUE,
                 multicore = TRUE, mc.cores = max(1, parallel::detectCores() - 1))
best <- mt$Model[which.min(mt$BIC)]
cat("\nBest model (BIC):", best, "\n")
fit <- pml_bb(mt, rearrangement = "NNI", control = pml.control(trace = 0))  # uses best-BIC model
bs  <- bootstrap.pml(fit, bs = n_boot, optNni = TRUE, multicore = TRUE,
                     mc.cores = max(1, parallel::detectCores() - 1),
                     control = pml.control(trace = 0))
tr  <- plotBS(fit$tree, bs, type = "none", p = 0)
tr  <- midpoint(tr, node.labels = "support")
sup <- suppressWarnings(as.numeric(tr$node.label))
if (all(is.na(sup) | sup <= 1)) sup <- sup * 100          # express support as %
tr$node.label <- ifelse(is.na(sup), "", as.character(round(sup)))
write.tree(tr, tree_file)
}

## ---- 7. figure (ggtree) ------------------------------------------------------
## Tip labels and legend use plotmath so that every scientific name is italic.
abbr <- function(x) sub("^Cephalopholis ", "C. ", x)
esc  <- function(x) gsub("'", "\\\\'", x)                 # escape quotes for plotmath
grp_code <- c("Saudi 'C. sonnerati' (query)" = "query",
              "Same unpublished Saudi study" = "same_study",
              "C. oligosticta (independent)" = "olig_ind",
              "C. sonnerati (reference)"     = "sonn_ref",
              "Other Cephalopholis"          = "other")
acc_short <- sub("\\.\\d+$", "", d$acc)
where     <- ifelse(d$country != "", paste0(" (", d$country, ")"), "")
meta <- data.frame(
  label = d$acc,
  group = unname(grp_code[d$group]),
  tip   = ifelse(d$acc %in% queries,
                 sprintf("bolditalic('%s')~bold('%s%s')", esc(abbr(d$organism)),
                         acc_short, esc(where)),
                 sprintf("italic('%s')~'%s%s'", esc(abbr(d$organism)),
                         acc_short, esc(where))),
  stringsAsFactors = FALSE)

grp_breaks <- c("query", "same_study", "olig_ind", "sonn_ref", "other")
pal <- c(query = "#C0392B", same_study = "#E67E22", olig_ind = "#7D3C98",
         sonn_ref = "#1F618D", other = "grey50")
grp_labels <- c(
  expression("Saudi records deposited as"~italic("C. sonnerati")~"(queries)"),
  expression(italic("C. oligosticta")~"from the same unpublished Saudi study"),
  expression(italic("C. oligosticta")~"(independent voucher)"),
  expression(italic("C. sonnerati")~"(reference sequences)"),
  expression("Other"~italic("Cephalopholis")~"species"))

draw_tree <- function(tr, label_size = 2.3, scale_w = 0.005) {
  p <- ggtree(tr, linewidth = 0.35) %<+% meta
  bsv <- suppressWarnings(as.numeric(p$data$label))         # bootstrap values (internal nodes)
  p$data$bs_lab <- ifelse(!p$data$isTip & !is.na(bsv) & bsv >= 50, as.character(round(bsv)), "")
  p +
    geom_tippoint(aes(colour = group), size = 1.3) +
    geom_tiplab(aes(label = tip, colour = group), parse = TRUE,
                size = label_size, offset = 0.001, show.legend = FALSE) +
    geom_nodelab(aes(label = bs_lab), size = label_size * 0.82, hjust = 1.25, vjust = -0.4,
                 colour = "grey20") +
    scale_colour_manual(values = pal, breaks = grp_breaks, labels = grp_labels,
                        name = "Sequence source") +
    geom_treescale(width = scale_w, fontsize = 2.4, linesize = 0.3) +
    guides(colour = guide_legend(ncol = 1, title.position = "top", override.aes = list(size = 2.6))) +
    theme(legend.position = "bottom",
          legend.justification = "left",
          legend.title = element_text(size = 8, face = "bold"),
          legend.text  = element_text(size = 7.5),
          legend.key.height = unit(0.38, "cm"),
          legend.margin = margin(t = 2, l = 4)) +
    hexpand(0.55)
}

## (a) full tree -> supplementary figure
p <- draw_tree(tr)
h <- max(6, 0.13 * Ntip(tr) + 1.5)
ggsave(file.path(out_dir, "Box1_16S_ML_tree_full.pdf"), p, width = 7.5, height = h, limitsize = FALSE)
ggsave(file.path(out_dir, "Box1_16S_ML_tree_full.png"), p, width = 7.5, height = h, dpi = 300,
       limitsize = FALSE)

## (b) focused clade for the Box: smallest clade containing the Saudi queries,
##     the independent C. oligosticta and the main C. sonnerati references
focus <- intersect(c(queries, "KM656818.1", "KM261673.1", "HQ592260.1", "DQ088037.1"), tr$tip.label)
subtr <- extract.clade(tr, getMRCA(tr, focus))
ps  <- draw_tree(subtr, label_size = 3, scale_w = 0.005)
hs  <- max(4, 0.2 * Ntip(subtr) + 1.5)
ggsave(file.path(out_dir, "Box1_16S_ML_tree_focus.pdf"), ps, width = 7.5, height = hs)
ggsave(file.path(out_dir, "Box1_16S_ML_tree_focus.png"), ps, width = 7.5, height = hs, dpi = 300)
cat("Focused clade:", Ntip(subtr), "tips\n")
write.table(meta, file.path(out_dir, "box1_tip_metadata.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

cat("\nDone. Outputs in", out_dir, "\n")
cat("Model:", best, "| tips:", Ntip(tr), "| bootstrap replicates:", n_boot, "\n")
sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()
