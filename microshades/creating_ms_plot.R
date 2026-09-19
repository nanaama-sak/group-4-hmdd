## Author: Juliette Kudrick 
## Purpose: Create a Script that will Create Microshades Plot

## load ps object
load("tutorial2_objects.RData")
  
## prep for microshades
library(microshades)
library(tidyverse)

## new tax table with genus + family level
tax_w_fam <- read.csv(file = "taxonomy_temp.csv")
rownames(tax_w_fam) <- tax_w_fam[,1] 
tax_w_fam <- tax_w_fam[, 2:ncol(tax_w_fam)]
tax_w_fam$Family <- trimws(tax_w_fam$Family)

## create a new ps object
ps_upd <- ps
tax_table(ps_upd) <- tax_table(as.matrix(tax_w_fam))

## make ps object smaller (Visit 1 only; people with known IR/IS classification)
ps_upd_v1 <- subset_samples(ps_upd, VisitID == "01") %>%
  subset_samples(., ir_is_classification != "Unknown") %>%
  prune_taxa(taxa_sums(.) > 0, .)

temp_mdf <- microshades::prep_mdf(ps_upd_v1, subgroup_level = "Family", as_relative_abundance = TRUE)


temp_cdf <- create_color_dfs(temp_mdf, selected_groups = c("Clostridiaceae", 
                                                           "Oscillospiraceae",
                                                           "Lachnospiraceae", 
                                                           "Bacteroidaceae"),
                             top_n_subgroups = 4, group_level = "Family", 
                             subgroup_level = "Genus")

ord_df <- reorder_samples_by(temp_cdf$mdf, temp_cdf$cdf, order_tax = "Bacteroidaceae",
                   group_level = "Family", subgroup_level = "Genus")

## create plot
ms_plot <- plot_microshades(ord_df$mdf, ord_df$cdf, group_label = "Family - Genus") +
  facet_wrap(~ir_is_classification, scale = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ms_plot

ggsave(file = "ir_is_ms_plot.png", plot = ms_plot, width = 20, height = 10)
  
