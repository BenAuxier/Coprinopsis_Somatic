#!/usr/bin/env Rscript

### Counts of somA-like loci
#############################################################################
# =======================================
# Sandra Lorena Ament Velasquez
# 2026-05-21
# =======================================
library(ggplot2)
library(tidyr)
library(dplyr)

# To handle the tree
library(ggtree)
library(ggtreeExtra)
library(ape)


# ============================
# Data
# ============================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))       # Set working directory to source file location
getwd()

counts <- read.table("../../2_FindSomAs/results/Loci_counts.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  stringsAsFactors = FALSE ) %>% select(-c(NACHTS_locustags, Kinases_locustags))

species <- read.table("../data/Sample_Species_map.txt", 
                      header = TRUE,
                      sep = "\t",
                      quote = "")

counts <- merge(species, counts, by = "Strain") %>% rename(Loci = N_candidate_loci, NACHT = N_NACHTs, Kinases = N_Kinases)

# --- Read tree ---
tree <- read.tree("../data/Mushrooms.tre")

# --- Output ---
path2figures <- "../results/"

# ============================
# Visualizing data
# ============================

# --- Root the tree on Agaricus bisporus ---
outgroup <- "Agaricus_bisporus_var_bisporus__Agabi_varbisH97_2"
tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)

# Optional: ladderize for a tidier display
tree <- ladderize(tree)

# Quick sanity check
is.rooted(tree)   # should be TRUE
# plot(tree); nodelabels(); tiplabels()  # uncomment to eyeball

# --- Split tip labels into Species + Strain ---
tip_df <- data.frame(full_label = tree$tip.label) %>%
  separate(full_label, into = c("Species_tree", "Strain"),
           sep = "__", remove = FALSE) %>%
  mutate(Species_tree = gsub("_", " ", Species_tree))

# Rename tree tips to just the Strain (so they match counts$Strain)
tree$tip.label <- tip_df$Strain

# Sanity check
setdiff(tree$tip.label, counts$Strain)
setdiff(counts$Strain, tree$tip.label)
# Both should be character(0)

# --- Build a tip-data table to attach species names + counts ---

# Map: original Strain -> shortened name to show in parentheses
strain_pretty <- c(
  "Candolleomyces_eurysporus_reassembly_enz"        = "reassembly_enz",
  "Candolleomyces_efflorescens_reassembly_enz"      = "reassembly_enz",
  "GCA_900156845.1_ASM90015684v1"                   = "ASM90015684v1",
  "GCA_951394405.1_gfCopMica1.1"                    = "gfCopMica1.1",
  "Coprinellus_aureogranulatus_MPG-14n_reassembly"  = "MPG-14n_reassembly",
  "Laccaria_amethystina_LaAM-08-1_v1.0"             = "LaAM-08-1_v1.0",
  "C_micaeus_w1"             = "w1",
  "C_micaeus_DM1047"             = "DM1047",
  "Agabi_varbisH97_2"             = "H97_2"
)

tip_info <- tip_df %>%
  select(Strain, Species_tree) %>%
  mutate(Strain_pretty = ifelse(Strain %in% names(strain_pretty),
                           strain_pretty[Strain],
                           Strain),
    display_label = paste0(Species_tree, " (", Strain_pretty, ")")
  ) %>%
  left_join(counts, by = "Strain")

# --- Heatmap data (rownames must match tip labels) ---
heat_data <- tip_info %>%
  select(Loci, NACHT, Kinases) %>%
  as.data.frame()
rownames(heat_data) <- tip_info$Strain

# --- Build the tree plot, attaching tip metadata ---
( treeplot <- ggtree(tree) %<+% tip_info +
  geom_tiplab(aes(label = display_label),
              align = TRUE, size = 4, fontface = "italic") +
  geom_nodelab(aes(label = label,
                   subset = !is.na(suppressWarnings(as.numeric(label))) ),
               size = 2.5, hjust = -0.1, vjust = 0) +
  geom_treescale(x = 0, y = 15, width = 0.1, fontsize = 3) +
  ylim(0, nrow(tip_info) + 1) + # To make space for the heatmap col names
  xlim(0, 4) ) # extend x so labels + heatmap fit
  # xlim(0, 3) ) # extend x so labels + heatmap fit (version 2)

# --- Add the heatmap ---
gheatmap(treeplot, heat_data,
         # offset = 1.6, # v2
         offset = 1.75, 
         width  = 1,
         colnames = TRUE,
         colnames_angle = 0,
         colnames_position = "top",
         hjust = 0.5,
         font.size = 3.5) +
  # scale_fill_viridis_c(name = "Count") +
  theme(legend.position = "right") +
  scale_fill_gradient(low = "white", high = "steelblue",
                      name = "Count", limits = c(0, NA)) +
  theme(legend.position   = "right",
        # legend.key.height = unit(0.4, "cm"),
        # legend.key.width  = unit(0.3, "cm"),
        legend.title      = element_text(size = 9),
        # legend.text       = element_text(size = 7),
        legend.margin     = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, -80),  # pull legend closer to plot
        plot.margin       = margin(0, 0, 0, -22)) # t, r, b, l in pts

ggsave(paste0(path2figures, "tree_heatmap.png"), width = 8, height = 5)

