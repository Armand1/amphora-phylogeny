## make_tree_plot.R -----------------------------------------------------------
## Plot the dated MCC tree of the 39-tip three-clade run: pot silhouettes at the
## tips filled by clade, tip labels coloured by origin, posterior support at nodes,
## calendar (BCE) time axis.
## ---------------------------------------------------------------------------
suppressMessages({library(treeio); library(ggtree); library(ape)
                  library(ggplot2); library(dplyr)})

HERE <- "/Users/armandleroi/Documents/Amphora_phylogeny/three_clade_tree"
DATA <- "/Users/armandleroi/Documents/Amphora_phylogeny/2026_amphora_297_types"
FIGW <- 10; FIGH <- 12

trb <- read.beast(file.path(HERE, "hpc_output/simple_run_mcc.tree"))
phy <- as.phylo(trb)

## --- metadata (clade + origin) from the selection ---------------------------
sel  <- read.csv(file.path(HERE, "selected_taxa.csv"), stringsAsFactors = FALSE)
meta <- sel |>
  transmute(label = type,
            clade  = dplyr::recode(sub("_(greek|ital)$", "", group),
                                   loutro = "loutrophoros", pan = "panathenaic",
                                   outgroup = "neck (outgroup)"),
            origin = ifelse(place == "greek_italiote", "Italiote", "Greek"))

## --- absolute chronology (BCE) ---------------------------------------------
## tip ages (BP) = midpoint of the age prior; anchor the tree by them.
ages <- read.delim(file.path(HERE, "data/mini_taxa.tsv"))
tipBP <- setNames((ages$min_age + ages$max_age) / 2, ages$taxon)
depth <- node.depth.edgelength(phy)                 # root->node distance (years)
ntip  <- Ntip(phy)
rootBP <- median(tipBP[phy$tip.label] + depth[1:ntip])   # consistent across tips
nodeBP <- rootBP - depth
nodeYr <- 1950 - nodeBP                              # calendar year (neg = BCE)
root_yr <- nodeYr[ntip + 1]
cat(sprintf("root age: %.0f BCE   (origin/oldest tip span %.0f to %.0f BCE)\n",
            -root_yr, -min(nodeYr), -max(nodeYr)))

## --- base tree (x = years, scaled so tips align to their dates) -------------
p <- ggtree(trb, mrsd = NULL) + theme_tree2()
## ggtree x runs 0..rootheight; convert axis to BCE via nodeYr on the layout
xr <- max(p$data$x)
p$data$year <- root_yr + (p$data$x)                 # x increases toward present
tipxy <- p$data |> filter(isTip) |> transmute(label, tx = x, ty = y)

## --- pot silhouettes, aligned column at the right --------------------------
prof <- read.csv(file.path(DATA, "2026_amphora_297_types_profiles_with_metadata.csv"),
                 stringsAsFactors = FALSE)
xt  <- max(tipxy$tx); gap <- 0.05 * xr
k_asp <- (xr * 1.9 / ntip) * (FIGH / FIGW)
sy <- 0.85 / 3; sx <- sy * k_asp
pots <- prof |>
  filter(type %in% phy$tip.label) |>
  transmute(label = type, point_order, px = x, py = y) |>
  left_join(tipxy, by = "label") |>
  left_join(meta, by = "label") |>
  arrange(label, point_order) |>
  mutate(x = xt + gap + px * sx, y = ty + py * sy)

lab <- meta |> left_join(tipxy, by = "label") |>
  mutate(x = xt + gap + 2.4 * max(abs(prof$x)) * sx,
         txt = sprintf("%s  %s", origin, clade))

pal_clade <- c(pelike = "#E69F00", loutrophoros = "#009E73",
               panathenaic = "#56B4E9", `neck (outgroup)` = "#999999")
pal_origin <- c(Greek = "grey20", Italiote = "#D55E00")

## --- BCE axis breaks --------------------------------------------------------
brks_yr <- seq(-1200, 0, by = 200)
brks_x  <- brks_yr - root_yr
keep <- brks_x >= 0 & brks_x <= xr

g <- p +
  geom_tippoint(aes(colour = origin), size = 1.7,
                data = left_join(p$data, meta, by = "label")) +
  geom_polygon(data = pots, aes(x, y, group = label, fill = clade),
               colour = "grey25", linewidth = 0.2) +
  geom_text(data = lab, aes(x, ty, label = txt, colour = origin),
            hjust = 0, size = 2.6) +
  geom_text2(aes(subset = !isTip & !is.na(as.numeric(posterior)) &
                          as.numeric(posterior) >= 0.5,
                 label = sprintf("%.2f", as.numeric(posterior))),
             hjust = 1.2, vjust = -0.4, size = 2.3, colour = "grey35") +
  scale_fill_manual(values = pal_clade, name = "Clade") +
  scale_colour_manual(values = pal_origin, name = "Origin") +
  scale_x_continuous(breaks = brks_x[keep],
                     labels = paste0(-brks_yr[keep], " BCE")) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(6, 170, 6, 6),
        axis.text.x = element_text(size = 8)) +
  labs(title = "Dated MCC tree — three-clade amphora test (39 tips, 10 PCs)",
       subtitle = sprintf("root %.0f BCE · fill = clade · label colour = origin · posterior >=0.5 shown",
                          -root_yr))

ggsave(file.path(HERE, "tree_dated_pots.png"), g, width = FIGW, height = FIGH, dpi = 130)
cat("wrote tree_dated_pots.png\n")

## --- structure report: where does each clade sit? --------------------------
cat("\n== clade clustering (dated tree) ==\n")
for (cl in c("pelike","loutrophoros","panathenaic","neck (outgroup)")) {
  tips <- meta$label[meta$clade == cl]
  nd <- getMRCA(phy, tips); desc <- extract.clade(phy, nd)$tip.label
  cat(sprintf("  %-16s MRCA subtree %2d tips = %2d of clade + %2d intruders  (mono=%s)\n",
              cl, length(desc), sum(desc %in% tips), sum(!desc %in% tips),
              is.monophyletic(phy, tips)))
}
