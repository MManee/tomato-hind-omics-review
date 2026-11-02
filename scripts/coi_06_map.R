# =============================================================================
# Map of the geographic origin of public molecular records of Cephalopholis
# sonnerati, for assessment as a possible manuscript figure.
#
# Plotted (country level only): each country's records are placed at one point in
# the adjacent sea, chosen for display; the point does not indicate a sampling site.
# The Saudi Arabian records are deliberately placed on land, because they record the
# country only and do not show whether the fish came from the Red Sea or the Arabian Gulf.
#   - the curated C. sonnerati COI barcodes (n = 113), coloured by BOLD BIN
#     (main group BOLD:AAD4737, western Indian Ocean BOLD:AAB5431, India
#     BOLD:ACC0647, no BIN), point area proportional to the number of barcodes;
#   - the GenBank 16S/12S records deposited as C. sonnerati from Saudi Arabia
#     (country only; Section 2.3);
#   - the C. oligosticta COI reference with a recorded country (MH707291, Egypt);
#   - the source of the two published genome assemblies (hatchery fish, Hainan,
#     China; Xie et al. 2022, Lu et al. 2024).
# Barcodes without a recorded country are counted in the figure note, not plotted.
#
# Country outlines: maps::world (public domain); no download needed.
#
# Inputs : results/coi_analysis/ (outputs of coi_03_analysis.R)
# Outputs: results/coi_map/Fig_Csonnerati_records_map.pdf / .png
#          results/coi_map/map_points.tsv   (every plotted point and its count)
#          results/coi_map/map_summary.txt  (counts behind the figure)
#
# Run:
#   cd tomato-hind-omics-review   # repository root
#   Rscript scripts/coi_06_map.R
# =============================================================================
for (p in c("ggplot2", "maps")) if (!requireNamespace(p, quietly = TRUE))
  install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({ library(ggplot2); library(Biostrings) })

proj    <- Sys.getenv("PROJ", ".")
res_dir <- file.path(proj, "results/coi_analysis")
out_dir <- file.path(proj, "results/coi_map")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sumf <- file.path(out_dir, "map_summary.txt"); unlink(sumf)
say <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
need <- function(f) { if (!file.exists(f)) stop("Missing input: ", f); f }
focal <- "Cephalopholis sonnerati"

## ---- curated C. sonnerati COI set (same rules as coi_03_analysis.R) -------------
qc  <- read.delim(need(file.path(res_dir, "coi_qc_all_records.tsv")), colClasses = "character", quote = "")
scr <- read.delim(need(file.path(res_dir, "coi_nearest_neighbour_screen_all.tsv")), colClasses = "character", quote = "")
mis <- read.delim(need(file.path(res_dir, "csonnerati_probable_misidentifications.tsv")), colClasses = "character", quote = "")
fa  <- readDNAStringSet(need(file.path(res_dir, "coi_window_alignment_qc_pass.fasta")))
hdr <- strsplit(names(fa), "|", fixed = TRUE)
ids <- sapply(hdr, `[`, 1); spp <- gsub("_", " ", sapply(hdr, `[`, 2))
excl <- scr$id[grepl("closer_to_another_species|singleton_within", scr$flag)]
cur  <- ids[spp == focal & !ids %in% union(excl, mis$id)]
cs   <- qc[match(cur, qc$id), c("id", "country", "bin")]
cs$bin[is.na(cs$bin) | cs$bin == ""] <- "No BIN"
cs$group <- c("BOLD:AAD4737" = "BOLD:AAD4737 (main group)",
              "BOLD:AAB5431" = "BOLD:AAB5431 (western Indian Ocean)",
              "BOLD:ACC0647" = "BOLD:ACC0647 (India)")[cs$bin]
cs$group[is.na(cs$group)] <- ifelse(cs$bin[is.na(cs$group)] == "No BIN", "No BIN assigned",
                                    paste("Other BIN:", cs$bin[is.na(cs$group)]))
say("Curated C. sonnerati COI barcodes:", nrow(cs))
say("Without a recorded country (not plotted):", sum(cs$country == "Not recorded"))

## ---- approximate country positions (country level, not sampling sites) ----------
pos <- read.table(header = TRUE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE, text = "
country          | lat   | lon   | water
Australia        | -17.0 | 149.0 | Coral Sea
Bangladesh       |  20.5 |  90.5 | Bay of Bengal
Egypt            |  26.0 |  35.0 | Red Sea
India            |   7.0 |  76.0 | Laccadive Sea, off the southern tip
Indonesia        |  -5.5 | 113.0 | Java Sea
Japan            |  30.0 | 133.0 | Pacific Ocean
Kenya            |  -3.0 |  41.5 | western Indian Ocean
Madagascar       | -20.0 |  50.5 | western Indian Ocean
Malaysia         |   4.5 | 105.0 | South China Sea
Maldives         |   3.2 |  73.2 | Indian Ocean
Mauritius        | -21.0 |  59.0 | western Indian Ocean
Mozambique       | -20.0 |  37.5 | Mozambique Channel
New Caledonia    | -22.5 | 167.5 | Coral Sea
Oman             |  20.5 |  60.0 | Arabian Sea
Papua New Guinea |  -9.0 | 152.0 | Solomon Sea
Philippines      |  11.5 | 121.5 | Sibuyan Sea
Reunion          | -21.1 |  55.5 | western Indian Ocean
Seychelles       |  -4.7 |  56.5 | western Indian Ocean
Solomon Islands  |  -9.6 | 160.2 | Solomon Sea
Somalia          |   4.0 |  50.0 | western Indian Ocean
South Africa     | -31.0 |  34.5 | off KwaZulu-Natal
Sri Lanka        |   7.5 |  82.8 | Bay of Bengal
Taiwan           |  23.0 | 122.5 | Pacific Ocean
Tanzania         |  -6.5 |  41.0 | western Indian Ocean
Thailand         |  11.0 | 101.5 | Gulf of Thailand
Viet Nam         |  12.5 | 110.5 | South China Sea
Vietnam          |  12.5 | 110.5 | South China Sea
Yemen            |  13.5 |  50.0 | Gulf of Aden")
world <- ggplot2::map_data("world")
locate <- function(cn) {
  k <- match(cn, pos$country)
  if (!is.na(k)) return(c(pos$lat[k], pos$lon[k]))
  w <- world[world$region == cn, ]                       # fallback: mean of outline points (may fall on land)
  if (nrow(w)) { say("  NOTE: no offshore position defined for", cn, "- placed at the mean of its outline; add it to the table"); return(c(mean(w$lat), mean(w$long))) }
  say("  WARNING: no position for country", cn, "- not plotted"); c(NA, NA)
}

## ---- points: COI barcodes by country and BIN group ------------------------------------
pl <- cs[cs$country != "Not recorded", ]
agg <- aggregate(id ~ country + group, pl, length); names(agg)[3] <- "n"
agg <- agg[order(agg$country, agg$group), ]
xy  <- t(sapply(agg$country, locate)); agg$lat <- xy[, 1]; agg$lon <- xy[, 2]
for (cn in unique(agg$country)) {                         # side-by-side when a country has >1 group
  k <- which(agg$country == cn); m <- length(k)
  if (m == 1) next
  r <- 3.6 * sqrt(agg$n[k] / max(agg$n))                  # approximate circle radius in degrees (max_size = 10)
  x <- cumsum(c(0, r[-m] + r[-1] + 0.6))                   # circles touching, 0.6 degree gap
  agg$lon[k] <- agg$lon[k] + x - mean(range(x))
}
agg$type <- "COI barcode, C. sonnerati (curated)"

## ---- other records ----------------------------------------------------------------
other <- rbind(
  data.frame(country = "Saudi Arabia", group = NA, n = 7, lat = 24.0, lon = 45.0,
             type = "16S/12S records deposited as C. sonnerati (Section 2.3)"),
  data.frame(country = "Egypt", group = NA, n = 1, lat = 26.0, lon = 35.0,
             type = "COI reference, C. oligosticta (MH707291)"),
  data.frame(country = "China (Hainan)", group = NA, n = 2, lat = 19.2, lon = 109.8,
             type = "Genome assemblies, hatchery fish"))
pts <- rbind(agg[, c("country", "group", "n", "lat", "lon", "type")], other)
write.table(pts, file.path(out_dir, "map_points.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("\nPlotted points (country, group, n):")
for (k in seq_len(nrow(pts))) say(sprintf("  %-16s %-40s %3d", pts$country[k],
  ifelse(is.na(pts$group[k]), pts$type[k], pts$group[k]), pts$n[k]))

## ---- figure ----------------------------------------------------------------------------
cols <- c("BOLD:AAD4737 (main group)" = "#1F618D", "BOLD:AAB5431 (western Indian Ocean)" = "#D35400",
          "BOLD:ACC0647 (India)" = "#7D3C98", "No BIN assigned" = "grey70")
extra <- setdiff(unique(agg$group), names(cols)); if (length(extra)) cols[extra] <- "grey40"
nr <- sum(cs$country == "Not recorded")
g <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group), fill = "grey92", colour = "grey60", linewidth = 0.15) +
  geom_point(data = agg, aes(lon, lat, size = n, fill = group), shape = 21, colour = "black", stroke = 0.3, alpha = 0.9) +
  geom_text(data = agg, aes(lon, lat - 1.2 - 3.6 * sqrt(n / max(n)), label = n), size = 2.4, vjust = 1) +
  geom_point(data = other, aes(lon, lat, shape = type), size = 3, fill = "white", colour = "black", stroke = 0.6) +
  scale_size_area(max_size = 10, breaks = c(1, 5, 10, 30), name = "COI barcodes (n)") +
  scale_fill_manual(values = cols, name = "BOLD BIN") +
  scale_shape_manual(values = c("16S/12S records deposited as C. sonnerati (Section 2.3)" = 24,
                                "COI reference, C. oligosticta (MH707291)" = 23,
                                "Genome assemblies, hatchery fish" = 22), name = "Other records",
                     labels = c("16S/12S records deposited as C. sonnerati (Section 2.3)" =
                                  expression(paste("16S/12S records deposited as ", italic("C. sonnerati"), ", Saudi Arabia (country only)")),
                                "COI reference, C. oligosticta (MH707291)" =
                                  expression(paste("COI reference, ", italic("C. oligosticta"), " (MH707291, Egypt)")),
                                "Genome assemblies, hatchery fish" = "Genome assemblies (hatchery fish, Hainan)")) +
  guides(fill = guide_legend(override.aes = list(size = 4), order = 1, ncol = 2),
         size = guide_legend(order = 2), shape = guide_legend(order = 3, ncol = 1)) +
  annotate("text", x = 17, y = -38, hjust = 0, size = 2.6,
           label = sprintf("COI barcodes without a recorded country: n = %d (not plotted)", nr)) +
  coord_quickmap(xlim = c(15, 180), ylim = c(-40, 40), expand = FALSE) +
  labs(x = "Longitude (degrees E)", y = "Latitude (degrees)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 7, hjust = 0),
        legend.key = element_blank(),
        panel.background = element_rect(fill = "#EAF2F8"), panel.grid = element_blank())
ggsave(file.path(out_dir, "Fig_Csonnerati_records_map.pdf"), g, width = 7.5, height = 6.6)
ggsave(file.path(out_dir, "Fig_Csonnerati_records_map.png"), g, width = 7.5, height = 6.6, dpi = 300)
say("\nDone. Outputs:", file.path(out_dir, "Fig_Csonnerati_records_map.png"), "and", sumf)
