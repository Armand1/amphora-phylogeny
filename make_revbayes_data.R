## make_revbayes_data.R -------------------------------------------------------
## Emit RevBayes input for the 39-tip three-clade test tree, in the exact format
## of Joel's pilot (data/mini_character.nex + data/mini_taxa.tsv).
##
##   characters : coordinate-PCA scores of the 39 Karcher-aligned outlines
##                (self-consistent with three_clade_nj.Rmd; NOT the SRVF-tPCA
##                 basis of the pilot — see README).
##   ages       : min/max in years BP, datum 1950 (radiocarbon convention):
##                min_age = 1950 - date_latest ; max_age = 1950 - date_earliest
##
## Run:  Rscript make_revbayes_data.R           # 10 PCs (default)
##       Rscript make_revbayes_data.R 20        # override PC count
## ---------------------------------------------------------------------------
suppressMessages({library(dplyr)})

args   <- commandArgs(trailingOnly = TRUE)
N_PC   <- if (length(args)) as.integer(args[1]) else 10
HERE   <- "/Users/armandleroi/Documents/Amphora_phylogeny/phylogenies/2026_07_three_clade_tree"
DATA   <- "/Users/armandleroi/Documents/Amphora_phylogeny/Making amphora types/05_2026_amphora_297_types"
OUT    <- file.path(HERE, "data"); dir.create(OUT, showWarnings = FALSE)

sel  <- read.csv(file.path(HERE, "selected_taxa.csv"), stringsAsFactors = FALSE)
prof <- read.csv(file.path(DATA, "2026_amphora_297_types_profiles_with_metadata.csv"),
                 stringsAsFactors = FALSE)

## --- coordinate matrix (39 x 500), points in fixed order -------------------
coord <- prof |>
  filter(type %in% sel$type) |>
  arrange(type, point_order) |>
  group_by(type) |>
  summarise(v = list(c(x, y)), .groups = "drop")
X <- do.call(rbind, coord$v); rownames(X) <- coord$type

## --- PCA (centre, no scale) -> N_PC scores ---------------------------------
pca <- prcomp(X, center = TRUE, scale. = FALSE)
ve  <- pca$sdev^2 / sum(pca$sdev^2)
S   <- pca$x[, 1:N_PC, drop = FALSE]
cat(sprintf("using %d PCs (%.2f%% cumulative variance)\n", N_PC, 100 * cumsum(ve)[N_PC]))

## --- write NEXUS continuous character matrix (pilot format) ----------------
## filename carries the PC count so different sweeps coexist (mini_character_20.nex)
nex <- file.path(OUT, sprintf("mini_character_%d.nex", N_PC))
con <- file(nex, "w")
writeLines(c("#NEXUS",
             sprintf("[Data written by make_revbayes_data.R, %s]", format(Sys.time())),
             "BEGIN DATA;",
             sprintf("  DIMENSIONS NTAX=%d NCHAR=%d;", nrow(S), N_PC),
             "  FORMAT DATATYPE=CONTINUOUS MISSING=? GAP=- INTERLEAVE=NO;",
             "  MATRIX"), con)
pad <- max(nchar(rownames(S))) + 4
for (t in rownames(S))
  writeLines(sprintf("    %-*s%s", pad, t,
                     # fixed-point only: RevBayes' NEXUS reader rejects 1.2e-03 etc.
                     paste(formatC(S[t, ], format = "f", digits = 12), collapse = "\t")), con)
writeLines(c("  ;", "END;"), con); close(con)

## --- write taxa TSV (min/max age, BP, datum 1950) --------------------------
ages <- prof |>
  filter(type %in% sel$type) |>
  distinct(type, date_earliest, date_latest) |>
  transmute(taxon   = type,
            min_age = 1950 - date_latest,     # youngest bound
            max_age = 1950 - date_earliest) |> # oldest bound
  arrange(match(taxon, rownames(S)))
write.table(ages, file.path(OUT, "mini_taxa.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## --- sanity: ages positive & ordered; cross-check vs pilot where types overlap
stopifnot(all(ages$max_age >= ages$min_age), all(ages$min_age > 0))
pilot <- tryCatch(read.delim(
  "/Users/armandleroi/Desktop/simpletree/run_2026-07-22_100k_65PC/data/mini_taxa.tsv"),
  error = function(e) NULL)
if (!is.null(pilot)) {
  ov <- intersect(ages$taxon, pilot$taxon)
  if (length(ov)) {
    m <- merge(ages, pilot, by = "taxon", suffixes = c(".new", ".pilot"))
    cat(sprintf("\ncross-check vs pilot on %d shared taxa: ages identical = %s\n",
                nrow(m), all(m$min_age.new == m$min_age.pilot &
                             m$max_age.new == m$max_age.pilot)))
  } else cat("\n(no taxa shared with pilot to cross-check)\n")
}
cat(sprintf("\nwrote:\n  %s  (NTAX=%d NCHAR=%d)\n  %s  (%d taxa)\n",
            nex, nrow(S), N_PC, file.path(OUT, "mini_taxa.tsv"), nrow(ages)))
cat(sprintf("age range: %d - %d BP  (%.0f - %.0f BCE)\n",
            min(ages$min_age), max(ages$max_age),
            1950 - min(ages$min_age), 1950 - max(ages$max_age)))
