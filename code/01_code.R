# install pkd ####
#install.packages("tidyverse")
#install.packages("bookdown")
#install.packages("xtable")
#install.packages("bench")
#install.packages("readr")
#install.packages("ggcorrplot")

library(tidyverse)
library(bookdown)
library(xtable)
library(bench)
library(readr)
library(ggplot2)
library(ggrepel) # for text repel
library (devtools)
library(ggpubr) # for stat signif
library(ggcorrplot) # for corelation signif
library(dplyr)
library (devtools)
library(readr)
library(ggplot2)
library(dplyr)

# read data frame ----
df_promoter <- read.table("DF_PROMOTER.BED", header = T)
DF_ATAC <- read.table("DF_ATAC.BED", header = T)
DF_ENHANCER <- read.table("DF_ENHANCER.BED", header = T)

# prepare data frame ----
## --- promoter: rebuilt groups
df_promoter <- df_promoter %>%
  mutate(
    class_expr = case_when(
      logFC_expr >  1 & FDR_expr < 0.05 ~ "up",
      logFC_expr < -1 & FDR_expr < 0.05 ~ "down",
      TRUE                              ~ "neutral"
    ),
    # ARID1B-KD class at the looser 0.3 cutoff  <-- VERIFY name/rule vs your file
    class_expr_arid1bKD = case_when(
      logFC_expr_arid1bKD >  0.3 & FDR_expr_arid1bKD < 0.05 ~ "up",
      logFC_expr_arid1bKD < -0.3 & FDR_expr_arid1bKD < 0.05 ~ "down",
      TRUE                                                  ~ "neutral"
    ),
    interest_group = factor(
      case_when(
        interest_group == "MLL_target_interest"                          ~ "MLL_target_interest",
        str_starts(interest_group, "MLL_target") & class_expr == "up"   ~ "MLL_target_up",
        str_starts(interest_group, "MLL_target") & class_expr == "down" ~ "MLL_target_down",
        str_starts(interest_group, "MLL_target")                           ~ "MLL_target_neutral",
        TRUE                                                               ~ "nonMLL_target"
      ),
      levels = c("MLL_target_interest", "MLL_target_up", "MLL_target_neutral",
                 "MLL_target_down", "nonMLL_target")
    )
  )

## --- enhancer: rebuilt groups
DF_ENHANCER <- DF_ENHANCER %>%
  mutate(class_atac = case_when(
    logFC_atac >  1 & FDR_atac < 0.05 ~ "up",
    logFC_atac < -1 & FDR_atac < 0.05 ~ "down",
    TRUE                              ~ "neutral"
  ))

## --- df_selected_4
# NOTE: selects class_atac etc. from RAW DF_ENHANCER (original classes), exactly
#       as your code did -- NOT from DF_ENHANCER_fc1.
df_selected <- dplyr::select(
  DF_ENHANCER,
  name_dhs_dox,
  hic_gene_dox, dist_to_hic_gene_tss_dox, interaction_score_dox,
  interaction_class_dox, closest_gene_dox, dist_to_closest_gene_tss_dox,
  hic_gene_veh, closest_gene_veh, interaction_score_veh,
  dist_to_hic_gene_tss_veh, interaction_class_veh, dist_to_closest_gene_tss_veh,
  logFC_atac, class_atac, novel_atac_peak,
  class_h3k27ac, logFC_h3k27ac, class_h3k4me1, logFC_h3k4me1,
  class_activity
)

expr <- df_promoter %>%
  dplyr::transmute(closest_gene_dox   = name,        # explicit join key (was colnames()[4])
                   class_expr_closest = class_expr,
                   logFC_expr_closest = logFC_expr) %>%
  dplyr::distinct(closest_gene_dox, .keep_all = TRUE)

df_selected_4 <- df_selected %>% dplyr::left_join(expr, by = "closest_gene_dox")
############################ Figure 1B (horizontal) ######################## ---- 
library(tidyverse)
library(clusterProfiler)
library(msigdbr)

# make sure a gene appears only once per condition
df_promoter_clean <- df_promoter %>% 
  group_by(name) %>%
  filter(n_distinct(class_expr_arid1bKD) == 1,
         n_distinct(class_expr) == 1) %>%
  ungroup()

# Get gene lists
up_genes <- df_promoter_clean %>% filter(class_expr == "up") %>% pull(name)
down_genes <- df_promoter_clean %>% filter(class_expr == "down") %>% pull(name)
all_genes <- df_promoter_clean %>% pull(name)

cat("UP genes:", length(up_genes), "\n")
cat("DOWN genes:", length(down_genes), "\n")
cat("All genes:", length(all_genes), "\n")

# 2. Get FULL MSigDB gene sets

# C2: Curated gene sets (CGP, CP, etc.)
msig_c2 <- msigdbr(species = "Mus musculus", category = "C2") %>%
  dplyr::select(gs_name, gene_symbol)

# H: Hallmark
msig_h <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# C5 GO:BP: Biological Process
msig_gobp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)

# Combine ALL
msig_all <- bind_rows(msig_c2, msig_h, msig_gobp)

cat("\nTotal gene sets:", length(unique(msig_all$gs_name)), "\n")

# 3. Prepare gene list for compareCluster

# Create a named list of gene vectors
gene_clusters <- list(
  Upregulated = up_genes,
  Downregulated = down_genes
)

# 4. Run compareCluster
cc_result <- compareCluster(
  geneClusters = gene_clusters,
  fun = "enricher",
  universe = all_genes,
  TERM2GENE = msig_all
)

# Convert to data frame
cc_df <- as.data.frame(cc_result)

# --- Source data for reviewers: full gene-set ORA table (p.adjust < 0.05) ---
ora_out <- as.data.frame(cc_result)
ora_out <- ora_out[!is.na(ora_out$p.adjust) & ora_out$p.adjust < 0.05, ]
write.csv(ora_out, "Fig1B_geneset_ORA.csv", row.names = FALSE)

# 6. Gene Set ORA Visualization 

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# 6.1. Load compareCluster results

df <- cc_df

# 6.2. Define selected gene sets - TOP 5 per category per direction
#    With source annotations for CGP only (author names)
selected_genesets <- tribble(
  ~geneset_id, ~display_name, ~theme,
  
  # 1) LEUKEMIA (Top significant in Upregulated) 
  "HUANG_AML_LSC47", "LSC47 (Huang)", "Leukemia",
  "KEGG_MEDICUS_VARIANT_MLL_ENL_FUSION_TO_TRANSCRIPTIONAL_ACTIVATION", "MLL-ENL fusion targets (KEGG)", "Leukemia",
  "WANG_IMMORTALIZED_BY_HOXA9_AND_MEIS1_DN", "HOXA9/MEIS1 vs HOXA9 UP (Wang)", "Leukemia",
  "HESS_TARGETS_OF_HOXA9_AND_MEIS1_UP", "HOXA9/MEIS1 targets UP (Hess)", "Leukemia",
  # ESC-like stemness programs (co-opted by cancer/leukemia)
  #"WONG_EMBRYONIC_STEM_CELL_CORE", "Embryonic stem cell core (Wong)", "Leukemia",
  #"BHATTACHARYA_EMBRYONIC_STEM_CELL", "Embryonic stem cell (Bhattacharya)", "Leukemia",
  #"RAMALHO_STEMNESS_UP", "Stemness UP (Ramalho)", "Leukemia",
  
  # 2) HEMATOPOIESIS 
  "JAATINEN_HEMATOPOIETIC_STEM_CELL_UP", "HSC UP (Jaatinen)", "Hematopoiesis",
  # Complete Ivanova hierarchy (stem cell to mature)
  #"IVANOVA_HEMATOPOIESIS_STEM_CELL", "Stem cell (Ivanova)", "Hematopoiesis",
  "IVANOVA_HEMATOPOIESIS_STEM_CELL_LONG_TERM", "LT-HSC (Ivanova)", "Hematopoiesis",
  #"IVANOVA_HEMATOPOIESIS_STEM_CELL_SHORT_TERM", "ST-HSC (Ivanova)", "Hematopoiesis",
  "IVANOVA_HEMATOPOIESIS_EARLY_PROGENITOR", "Early progenitor (Ivanova)", "Hematopoiesis",
  #"IVANOVA_HEMATOPOIESIS_INTERMEDIATE_PROGENITOR", "Intermediate progenitor (Ivanova)", "Hematopoiesis",
  "IVANOVA_HEMATOPOIESIS_MATURE_CELL", "Mature cell (Ivanova)", "Hematopoiesis",
  # Other hematopoiesis CGP sets
  "BROWN_MYELOID_CELL_DEVELOPMENT_UP", "Myeloid development UP (Brown)", "Hematopoiesis",
  "KAMIKUBO_MYELOID_CEBPA_NETWORK", "Myeloid CEBPA network (Kamikubo)", "Hematopoiesis",
  # Lymphoid markers (to show myeloid vs lymphoid)
  #"LEE_EARLY_T_LYMPHOCYTE_UP", "Early T lymphocyte UP (Lee)", "Hematopoiesis",
  #"MORI_MATURE_B_LYMPHOCYTE_UP", "Mature B lymphocyte UP (Mori)", "Hematopoiesis",
  
  # 3) BIOLOGICAL PATHWAYS (GO:BP) 
  # Upregulated - Development/Morphogenesis
  "GOBP_EMBRYONIC_ORGAN_DEVELOPMENT", "Embryonic organ development", "Biological pathways",
  "GOBP_RENAL_SYSTEM_DEVELOPMENT", "Skeletal system development", "Biological pathways",
  "GOBP_SENSORY_ORGAN_MORPHOGENESIS","Sensory organ morphogenesis","Biological pathways",
  "GOBP_SKELETAL_SYSTEM_MORPHOGENESIS", "Renal system morphogenesis", "Biological pathways",
  # Upregulated - Metabolism
  "GOBP_AMINO_ACID_METABOLIC_PROCESS", "Amino acid metabolism", "Biological pathways",
  "GOBP_ALPHA_AMINO_ACID_METABOLIC_PROCESS", "Alpha amino acid metabolic process", "Biological pathways",
  "GOBP_ORGANIC_ACID_BIOSYNTHETIC_PROCESS", "Organic acid biosynthesis", "Biological pathways",
  
  # Downregulated (immune & migration)
  "GOBP_CELL_CHEMOTAXIS", "Cell chemotaxis", "Biological pathways",
  "GOBP_LEUKOCYTE_MIGRATION", "Leukocyte migration", "Biological pathways",
  "GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE", "Inflammatory response regulation", "Biological pathways",
  "GOBP_PHAGOCYTOSIS", "Phagocytosis", "Biological pathways",
  "GOBP_REGULATION_OF_LEUKOCYTE_PROLIFERATION", "Regulation of leukocyte proliferation", "Biological pathways",
  
  
  # 4) HALLMARK - Top significant each direction 
  # Upregulated
  "HALLMARK_MYC_TARGETS_V2", "MYC targets V2", "Hallmark",
  "HALLMARK_MTORC1_SIGNALING", "mTORC1 signaling", "Hallmark",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "Unfolded protein response", "Hallmark",
  # Downregulated
  "HALLMARK_INFLAMMATORY_RESPONSE", "Inflammatory response", "Hallmark",
  "HALLMARK_IL2_STAT5_SIGNALING", "IL2/STAT5 signaling", "Hallmark",
  "HALLMARK_KRAS_SIGNALING_UP", "KRAS signaling UP", "Hallmark",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNF/NF-kB signaling", "Hallmark",
  "HALLMARK_COMPLEMENT", "Complement", "Hallmark"
)


# 6.3. Prepare data for visualization

# Filter and join with display names
plot_data <- df %>%
  filter(ID %in% selected_genesets$geneset_id) %>%
  left_join(selected_genesets, by = c("ID" = "geneset_id")) %>%
  mutate(
    # Calculate -log10(p.adjust) for color
    log10_padj = -log10(p.adjust),
    # Cap at 15 for visualization
    log10_padj_capped = pmin(log10_padj, 15),
    # Significance stars
    sig = case_when(
      p.adjust < 0.001 ~ "***",
      p.adjust < 0.01 ~ "**",
      p.adjust < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # Order themes
    theme = factor(theme, levels = c("Leukemia", "Hematopoiesis", "Biological pathways", "Hallmark"))
  )

# Order gene sets by theme
geneset_order <- selected_genesets %>%
  mutate(theme = factor(theme, levels = c("Leukemia", "Hematopoiesis", "Biological pathways", "Hallmark"))) %>%
  arrange(theme) %>%
  pull(display_name)

plot_data <- plot_data %>%
  mutate(display_name = factor(display_name, levels = rev(geneset_order)))

# 6.4. dotploting
# Create a complete scaffold with all gene sets and both clusters
scaffold <- expand.grid(
  display_name = geneset_order,
  Cluster = c("Upregulated", "Downregulated"),
  stringsAsFactors = FALSE
) %>%
  left_join(selected_genesets %>% dplyr::select(display_name, theme), by = "display_name") %>%
  mutate(
    theme = factor(theme, levels = c("Leukemia", "Hematopoiesis", "Biological pathways", "Hallmark")),
    display_name = factor(display_name, levels = rev(geneset_order)),
    Cluster = factor(Cluster, levels = c("Upregulated", "Downregulated"))
  )

# Filter to only significant results for plotting points
plot_data_sig <- plot_data %>%
  filter(p.adjust < 0.05) %>%
  mutate(
    Cluster = factor(Cluster, levels = c("Upregulated", "Downregulated")),
    display_name = factor(display_name, levels = rev(geneset_order))
  )

Figure_2B<-ggplot() +
  # Invisible points to set up the scaffold/grid
  geom_point(data = scaffold, aes(x = Cluster, y = display_name), alpha = 0) +
  # Actual data points (only significant)
  geom_point(data = plot_data_sig, aes(x = Cluster, y = display_name, size = Count, fill = log10_padj_capped), 
             shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient(
    low = "white", high = "darkgreen",
    limits = c(0, 15),
    name = "-log10\n(p.adjust)"
  ) +
  scale_size_continuous(range = c(2, 10), name = "Gene count") +
  scale_x_discrete(labels = c("Upregulated" = "RNA up", "Downregulated" = "RNA down")) +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", face = "bold", angle = 45, hjust = 1),
    axis.title = element_blank(),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 10),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(Figure_2B)
# Figuer_1D ----
df_pie <- DF_ATAC %>%
  filter(peak_hmll1 == "yes") %>%
  count(annot.type) %>%
  mutate(percentage = n / sum(n) * 100,
         label = paste0(annot.type, "\n", n, " (", round(percentage, 1), "%)"))


Figuer_1D<-ggplot(df_pie, aes(x = "", y = n, fill = annot.type)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            size = 3) +
  labs(title = "peak_hmll1 - Pie Chart", x = NULL, y = NULL, fill = "Annotation Type") +
  theme_pubr() +
  theme(text = element_text(size = 10))+
  scale_fill_manual(values = c(  "#74ADD1", "#4575B4", "#D6604D")) 

print(Figuer_1D)
# Figure_1F_bottom----
representives <- filter(df_promoter, name == "Hoxa9" | name == "Meis1" |
                   name == "Arid1b" | name == "Six1" |
                   name == "Eya1" |
                   name == "Ikzf2")

Figure_1F_bottom<-ggplot(df_promoter,
                   aes(x = mean_hmll1_dox,
                       y = logFC_expr,
                       color = interest_group)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_smooth(method = "lm") +
  theme_pubr() +
  theme(text = element_text(size = 10),
        axis.text.y = element_text(angle = 0, vjust = 1, hjust = 1, size = 10)) +
  labs(title = "promoter_hmll1") +
  scale_color_manual(values = c("#984EA3","#E41A1C","#4DAF4A","#377EB8", "#999999")) +
  scale_fill_manual(values  = c("#984EA3","#E41A1C","#4DAF4A","#377EB8", "#999999")) +
  stat_cor(aes(color = interest_group), method = "spearman",
           label.x = 3, size = 4) +
  geom_text_repel(
    data    = representives,
    mapping = aes(label = name),
    size    = 5,
    color   = "#984EA3",
    box.padding       = 0.1,
    point.padding     = 0.2,
    max.overlaps      = 20,
    segment.size      = 0.3,
    segment.color     = "gray40",
    min.segment.length = 0.1,
    force    = 10,
    max.iter = 5000
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  geom_hline(yintercept = -1, linetype = "dashed", color = "black")

print(Figure_1F_bottom)

# Figure_S3F----
representives<- filter(df_promoter, name == "Hoxa9" | name == "Meis1" |
                   name == "Arid1b" | name == "Six1" |
                   name == "Eya1" |
                   name == "Ikzf2")

Figure_S3F<-ggplot(df_promoter,
       aes(x = mean_hmll1_dox,
           y = log2FoldChange_3h,
           color = interest_group)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_smooth(method = "lm") +
  theme_pubr() +
  theme(text = element_text(size = 15),
        axis.text.y = element_text(angle = 0, vjust = 1, hjust = 1, size = 10)) +
  labs(title = "promoter_hmll1") +
  scale_color_manual(values = c("#984EA3","#E41A1C","#4DAF4A","#377EB8", "#999999")) +
  scale_fill_manual(values  = c("#984EA3","#E41A1C","#4DAF4A","#377EB8", "#999999")) +
  stat_cor(aes(color = interest_group), method = "spearman",
           label.x = 6, size = 5) +
  geom_text_repel(
    data    = representives,
    mapping = aes(label = name),
    size    = 5,
    color   = "#984EA3",
    box.padding       = 0.1,
    point.padding     = 0.2,
    max.overlaps      = 20,
    segment.size      = 0.3,
    segment.color     = "gray40",
    min.segment.length = 0.1,
    force    = 10,#3,
    max.iter = 5000
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

print(Figure_S3F)
# Figure_S1E  ----
representives <- filter(df_promoter, name == "Hoxa9" | name == "Meis1" |
                   name == "Arid1b" | name == "Six1" |
                   name == "Eya1" |
                   name == "Ikzf2")

Figure_S1E<-ggplot(df_promoter,
         aes(x = mean_hmll1_dox,
             y = logFC_hmll1,
             color = interest_group)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_smooth(method = "lm") +
  theme_pubr() +
  theme(text = element_text(size = 10),
        axis.text.y = element_text(angle = 0, vjust = 1, hjust = 1, size = 10)) +
  labs(title = "promoter_hmll1") +
  scale_color_manual(values = c("#984EA3","#E41A1C","#4DAF4A","#377EB8", "#999999")) +
  scale_fill_manual(values  = c("#984EA3","#E41A1C","#4DAF4A","#377EB8", "#999999")) +
  stat_cor(aes(color = interest_group), method = "spearman",
           label.x = 3, size = 4) +
  geom_text_repel(
    data    = representives,
    mapping = aes(label = name),
    size    = 1,
    color   = "#984EA3",
    box.padding       = 0.1,
    point.padding     = 0.2,
    max.overlaps      = 20,
    segment.size      = 0.3,
    segment.color     = "gray40",
    min.segment.length = 0.1,
    force    = 10,
    max.iter = 5000
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  geom_hline(yintercept = -1, linetype = "dashed", color = "black")

print(Figure_S1E)
# Figure_1F_top ----
representives <- filter(df_promoter, name == "Hoxa9" | name == "Meis1" |
                   name == "Arid1b" | name == "Six1" |
                   name == "Eya1" |
                   name == "Ikzf2")


Figure_1F_top<- ggplot(df_promoter,
         aes(x = mean_hmll1_dox,
             y = logFC_h3k79me2,
             color = class_h3k79me2)) +#interest_group
  geom_point(size = 1, alpha = 0.8) +
  geom_smooth(method = "lm") +
  theme_pubr() +
  theme(text = element_text(size = 10),
        axis.text.y = element_text(angle = 0, vjust = 1, hjust = 1, size = 10)) +
  labs(title = "promoter_hmll1") +
  scale_color_manual(values = c("#377EB8", "#999999","#E41A1C")) +
  geom_text_repel(
    data    = representives,
    mapping = aes(label = name),
    size    = 2,
    color   = "#984EA3",
    box.padding       = 0.1,
    point.padding     = 0.2,
    max.overlaps      = 20,
    segment.size      = 0.3,
    segment.color     = "gray40",
    min.segment.length = 0.1,
    force    = 10,#3,
    max.iter = 5000
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  geom_hline(yintercept = -1, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 8, linetype = "dashed", color = "black")

print(Figure_1F_top)
######################## Figure_1I, Figure_S1H, Figure_S1J  ############# ----
gene_group <- dplyr::select(df_promoter, name, interest_group, class_expr)
colnames(gene_group)[1] <- "overlapping_gene"

# Factor levels for interest_group
ig_levels <- c("MLL_target_interest", "MLL_target_up", "MLL_target_neutral",
               "MLL_target_down", "nonMLL_target")
ig_colors <- c("#CAB2D6", "#E41A1C", "#4DAF4A", "#377EB8", "#999999")
names(ig_colors) <- ig_levels

# Filter MA9-bound non-promoter regions 
ma9_bound <- filter(DF_ATAC, peak_hmll1 == "yes", annot.type != "promoters")

ma9_bound <- ma9_bound %>%
  mutate(mll_enh = case_when(
    class_h3k79me2 == "up" & class_atac == "up" ~ "yes",
    TRUE ~ "no"
  )) %>%
  left_join(gene_group, by = "overlapping_gene") %>%
  mutate(interest_group = factor(interest_group, levels = ig_levels))

# Figure_1I top ----
Figure_1I<-ggplot(ma9_bound, aes(x = mll_enh, fill = annot.type)) +
  geom_bar(alpha = 1, position = "fill") +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values = c("#74ADD1", "#4575B4")) +
  labs(title = "MLL bound dCREs", x = "Gain both H3K79me2 and ATAC") +
  theme_pubr() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 15))

print(Figure_1I)
# Filter: mll_enh == "yes" 
yes <- filter(ma9_bound, mll_enh == "yes")

#  Figure S1J ----
yes$overlapping_gene <- factor(yes$overlapping_gene,
                               levels = names(sort(table(yes$overlapping_gene),
                                                   decreasing = TRUE)))
# all genes
ggplot(filter(yes, annot.type == "intragenic"),
       aes(x = overlapping_gene, fill = interest_group)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_stack(vjust = 1.1)) +
  scale_fill_manual(values = ig_colors) +
  labs(title = "MLL-bound dCREs", x = "Gene") +
  theme_pubr() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10))

# count > 2 only 
gene_counts  <- table(yes$overlapping_gene)
genes_keep   <- names(gene_counts[gene_counts > 2])
yes_filtered <- yes %>%
  filter(overlapping_gene %in% genes_keep) %>%
  mutate(overlapping_gene = factor(overlapping_gene,
                                   levels = names(sort(table(overlapping_gene),
                                                       decreasing = TRUE))))

Figure_S1J<-ggplot(filter(yes_filtered, annot.type == "intragenic"),
       aes(x = overlapping_gene, fill = interest_group)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_stack(vjust = 1.1)) +
  scale_fill_manual(values = ig_colors) +
  labs(title = "MA9-bound pE count > 2", x = "Gene") +
  theme_pubr() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15))

print(Figure_S1J)
#  Figure_1I_bottom ----
Figure_1I_bottom<-ggplot(filter(yes, annot.type == "intragenic"),
       aes(x = interest_group, fill = interest_group)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_stack(vjust = 1.1)) +
  scale_fill_manual(values = ig_colors) +
  labs(title = "MLL-bound dCREs", x = "Gene group (logFC >= 1)") +
  theme_pubr() +
  theme(legend.position = "right",
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15))

print(Figure_1I_bottom)
############################  Figure_S1I  ############################ ----
library(tidyverse)
library(patchwork)

# 1. Build per-gene summary 

# All genes with their group
all_genes <- df_promoter %>%
  dplyr::select(name, interest_group) %>%
  distinct(name, .keep_all = TRUE) %>%
  mutate(interest_group = factor(interest_group, levels = ig_levels))

# Count of MLL-bound intragenic enhancers per gene
enh_per_gene_counts <- yes %>%
  filter(annot.type == "intragenic") %>%
  count(overlapping_gene, name = "n_enh")

# Join to all genes — genes with no enhancer get 0
all_genes <- all_genes %>%
  left_join(enh_per_gene_counts, by = c("name" = "overlapping_gene")) %>%
  mutate(n_enh   = replace_na(n_enh, 0),
         has_enh = n_enh > 0)

cat("Total enhancers:", sum(all_genes$n_enh), "\n")
cat("Genes with >= 1 enhancer:", sum(all_genes$has_enh), "\n\n")

# Summary per group
group_summary <- all_genes %>%
  group_by(interest_group) %>%
  summarise(
    total_genes   = n(),
    genes_with_enh = sum(has_enh),
    total_enh     = sum(n_enh),
    enh_per_gene  = total_enh / total_genes,   # total enhancers / total genes
    pct_with_enh  = genes_with_enh / total_genes,
    .groups = "drop"
  )

cat("=== Group summary ===\n")
print(group_summary)

# 2. Fisher's exact OR for each group vs all others (background) 
run_or <- function(group_name, df) {
  in_grp  <- df$interest_group == group_name
  has_enh <- df$has_enh
  
  mat <- matrix(
    c(sum( in_grp &  has_enh),  sum(!in_grp &  has_enh),
      sum( in_grp & !has_enh),  sum(!in_grp & !has_enh)),
    nrow = 2,
    dimnames = list(c("in_group", "not_in_group"),
                    c("has_enh", "no_enh"))
  )
  
  ft <- fisher.test(mat)
  
  tibble(
    group  = group_name,
    log2OR = log2(ft$estimate),
    pval   = ft$p.value,
    ci_lo  = log2(ft$conf.int[1]),
    ci_hi  = log2(ft$conf.int[2])
  )
}

or_results <- map_dfr(ig_levels, run_or, df = all_genes) %>%
  mutate(
    padj = p.adjust(pval, method = "BH"),
    sig  = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ "ns"
    )
  ) %>%
  left_join(group_summary, by = c("group" = "interest_group")) %>%
  mutate(group = factor(group, levels = rev(ig_levels)))

cat("\n=== OR results ===\n")
print(or_results %>% dplyr::select(group, total_enh, total_genes,
                                   enh_per_gene, log2OR, padj, sig))

# 3. Colors 
ig_colors <- c(
  "MLL_target_interest" = "#CAB2D6",
  "MLL_target_up"       = "#E41A1C",
  "MLL_target_neutral"  = "#4DAF4A",
  "MLL_target_down"     = "#377EB8",
  "nonMLL_target"       = "#999999"
)

# 4. Forest plot: log2OR with raw counts annotated 
Figure_S1I<-ggplot(or_results,
               aes(x = log2OR, y = group, color = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.25, linewidth = 0.8) +
  geom_point(size = 4, alpha = 0.9) +
  # significance stars to the right of CI
  geom_text(aes(x = ci_hi, label = paste0(" ", sig)),
            hjust = 0, size = 5, fontface = "bold") +
  # genes_with_enh / total_genes — fixed position at left of plot
  geom_text(aes(label = paste0(genes_with_enh, "/", total_genes, " genes")),
            x = -Inf, hjust = -0.05, size = 3.5, color = "grey30") +
  scale_color_manual(values = ig_colors, guide = "none") +
  labs(
    x     = "log2(OR) vs genome background",
    y     = NULL,
    title = "Enrichment of MLL-bound intragenic enhancers\nper gene group (Fisher's exact, vs genome background)",
    caption = "Annotation: genes with ≥1 MLL-bound enhancer / total genes in group"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.caption  = element_text(size = 9, color = "grey40", hjust = 0),
    axis.text.y   = element_text(size = 11),
    axis.text.x   = element_text(size = 11),
    axis.title.x  = element_text(size = 11)
  )

print(Figure_S1I)
########################## Figure_2A_top ######################### ----
# Calculate counts per class
class_counts <- DF_ENHANCER %>%
  group_by(class_atac) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(class_atac, ": n=", n))

label_vec <- setNames(class_counts$label, class_counts$class_atac)

Figure_2A_top<-ggplot(DF_ENHANCER, 
       aes(y=-log10(FDR_atac), x=logFC_atac, color=class_atac))+
  geom_point(size=1,alpha=0.1)+
  # Vertical cutoff lines
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color="black", linewidth=0.5)+
  # Horizontal cutoff line
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="black", linewidth=0.5)+
  # Label for horizontal line only
  annotate("text", x=Inf, y=-log10(0.05), label="FDR=0.05", 
           vjust=-0.5, hjust=1.1, size=4, color="black")+
  #labs(title="log2FC_atac (on vs off)")+
  scale_color_manual(
    values = c("down" = "#377EB8", "neutral" = "#999999", "up" = "#E41A1C"),
    labels = label_vec
  )+
  # Add -1 and 1 as extra breaks on x axis
  scale_x_continuous(breaks = sort(c(-1, 1, scales::extended_breaks()(DF_ENHANCER$logFC_atac))))+
  theme_pubr(base_size = 20)

print(Figure_2A_top)
############################### Figure_2C ############################### ----
library(dplyr)
library(ggplot2)
library(tidyr)

df <- df_selected_4

# Keep unique enhancers with expression data
enh <- df[!duplicated(df$name_dhs_dox), ]
enh <- enh[!is.na(enh$logFC_expr_closest), ]

# STEP 2: Aggregate to gene level

gene_stats <- enh %>%
  group_by(closest_gene_dox) %>%
  summarise(
    n_enhancers        = n(),
    mean_atac_logFC    = mean(logFC_atac,         na.rm = TRUE),
    gene_expr          = logFC_expr_closest[1],
    class_expr         = class_expr_closest[1],
    .groups = "drop"
  )

# STEP 3: Split genes by ATAC sign, bin each half into 2 groups
#
#  negative side  →  bin 1 = ATAC--,  bin 2 = ATAC-
#  positive side  →  bin 1 = ATAC+,   bin 2 = ATAC++
#
# Then re-bin the full range into n_bins for the bar plot x-axis

n_bins <- 50  # number of bars across the full range

gene_stats <- gene_stats %>%
  mutate(
    # Fine bins for x-axis bars
    bin_num = ntile(mean_atac_logFC, n_bins),
    
    # Coarse 4-group classification for dashed-line annotations
    atac_group = case_when(
      mean_atac_logFC <  0 ~ ntile(mean_atac_logFC[mean_atac_logFC < 0], 2)[
        match(mean_atac_logFC, sort(mean_atac_logFC[mean_atac_logFC < 0]))],
      mean_atac_logFC >= 0 ~ 2L + ntile(mean_atac_logFC[mean_atac_logFC >= 0], 2)[
        match(mean_atac_logFC, sort(mean_atac_logFC[mean_atac_logFC >= 0]))]
    )
  )

# Cleaner way to assign the 4 groups (avoids nested ntile complexity):
neg_genes <- gene_stats %>% filter(mean_atac_logFC <  0) %>%
  mutate(atac_group = ntile(mean_atac_logFC, 2))          # 1=ATAC--, 2=ATAC-
pos_genes <- gene_stats %>% filter(mean_atac_logFC >= 0) %>%
  mutate(atac_group = ntile(mean_atac_logFC, 2) + 2L)     # 3=ATAC+,  4=ATAC++

gene_stats <- bind_rows(neg_genes, pos_genes)

# STEP 4: Compute per-bin statistics for the stacked bar

bin_info <- gene_stats %>%
  group_by(bin_num) %>%
  summarise(
    n_genes        = n(),
    mean_atac      = mean(mean_atac_logFC, na.rm = TRUE),
    min_atac       = min(mean_atac_logFC,  na.rm = TRUE),
    max_atac       = max(mean_atac_logFC,  na.rm = TRUE),
    n_up           = sum(class_expr == "up",      na.rm = TRUE),
    n_down         = sum(class_expr == "down",    na.rm = TRUE),
    n_neutral      = sum(class_expr == "neutral", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_up      = 100 * n_up      / n_genes,
    pct_down    = 100 * n_down    / n_genes,
    pct_neutral = 100 * n_neutral / n_genes,
    bin_label   = sprintf("%.2f~%.2f", min_atac, max_atac)
  )

# STEP 5: Find dashed-line positions between the 4 ATAC groups
#
#  The boundary bin numbers are where the atac_group changes:
#    ATAC-- | ATAC-  : between group 1 and 2 (within negatives)
#    ATAC-  | ATAC+  : between negatives and positives  (around 0)
#    ATAC+  | ATAC++ : between group 3 and 4 (within positives)

# Get the max bin_num for each atac_group → place dashed line after that bin
boundary_bins <- gene_stats %>%
  group_by(atac_group) %>%
  summarise(max_bin = max(bin_num), .groups = "drop") %>%
  arrange(atac_group)

# Dashed lines go at x = max_bin + 0.5 (between bars)
dash_positions <- boundary_bins$max_bin[-nrow(boundary_bins)] + 0.5

# Labels for the 4 groups: x position = midpoint of each group's bin range
group_label_x <- gene_stats %>%
  group_by(atac_group) %>%
  summarise(mid_bin = mean(range(bin_num)), .groups = "drop") %>%
  arrange(atac_group) %>%
  pull(mid_bin)

group_labels <- c("ATAC\u2212\u2212", "ATAC\u2212", "ATAC+", "ATAC++")

# STEP 6: Build stacked bar data (long format)

bar_data <- bin_info %>%
  select(bin_num, pct_down, pct_neutral, pct_up) %>%
  pivot_longer(
    cols      = c(pct_down, pct_neutral, pct_up),
    names_to  = "class",
    values_to = "pct"
  ) %>%
  mutate(
    class = factor(class,
                   levels = c("pct_up", "pct_neutral", "pct_down"),
                   labels = c("up", "neutral", "down"))
  )

# STEP 7: Plot

# Show label + tick only every 5th bin
labeled_bins <- c(1, seq(5, n_bins, by = 5), n_bins)
x_labels     <- ifelse(bin_info$bin_num %in% labeled_bins, bin_info$bin_label, "")
tick_colours <- ifelse(bin_info$bin_num %in% labeled_bins, "black", "transparent")

Figure_2C<-ggplot(bar_data, aes(x = factor(bin_num), y = pct, fill = class)) +
  geom_bar(stat = "identity", position = "stack", width = 1,
           colour = "black", linewidth = 0.1) +
  
  # Dashed vertical lines separating ATAC--, ATAC-, ATAC+, ATAC++
  geom_vline(xintercept = dash_positions,
             linetype = "dashed", colour = "black", linewidth = 0.5) +
  
  # Group labels at the top
  annotate("text",
           x     = group_label_x,
           y     = 103,
           label = group_labels,
           size  = 5, fontface = "bold", hjust = 0.5) +
  
  scale_fill_manual(
    values = c("up" = "#E05C44", "neutral" = "#BDBDBD", "down" = "#5B9BD5"),
    name   = "Gene response",
    breaks = c("up", "neutral", "down"),
    labels = c("up", "neutral", "down")
  ) +
  
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(limits = c(0, 107), expand = c(0, 0),
                     breaks = c(0, 25, 50, 75, 100)) +
  
  labs(
    x = "Mean ATAC log2FC of proximal pEs per gene",
    y = "% of Genes"
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, size = 15),
    axis.ticks.x     = element_line(colour = tick_colours),
    axis.line        = element_line(colour = "black"),
    legend.position  = "right",
    legend.title     = element_text(size = 15),
    legend.text      = element_text(size = 15),
    plot.margin      = margin(t = 20, r = 10, b = 40, l = 10)
  )

print(Figure_2C)
############################### Figure_2D, Figure_S2B ########################### ----

# STEP 1: Load and prepare data

df <- df_selected_4

# Get unique enhancers
enh <- df[!duplicated(df$name_dhs_dox), ]
enh <- enh[!is.na(enh$logFC_expr_closest), ]

# Aggregate to gene level
gene_stats <- enh %>%
  group_by(closest_gene_dox) %>%
  summarise(
    n_enhancers = n(),
    mean_atac_logFC = mean(logFC_atac, na.rm = TRUE), #logFC_h3k27ac, #logFC_atac, #logFC_h3k4me1
    gene_expr = logFC_expr_closest[1],
    class_expr = class_expr_closest[1],
    .groups = "drop"
  )

# STEP 2: Filter concordant genes and create bins

n_bins <- 2
n2_bins <- 1

# Concordant UP: H3K27ac > 0 AND class_expr = 'up'
concordant_up <- gene_stats %>%
  filter(mean_atac_logFC > 0 & class_expr == "up",gene_expr>1) %>%
  mutate(bin_num = ntile(mean_atac_logFC, n_bins))

disconcordant_up <- gene_stats %>%
  filter(mean_atac_logFC < 0 & class_expr == "up",gene_expr>1) %>%
  mutate(bin_num = ntile(mean_atac_logFC, n2_bins))

# Concordant DOWN: H3K27ac < 0 AND class_expr = 'down'
concordant_down <- gene_stats %>%
  filter(mean_atac_logFC < 0 & class_expr == "down",gene_expr< -1) %>%
  mutate(bin_num = ntile(mean_atac_logFC, n_bins))

disconcordant_down <- gene_stats %>%
  filter(mean_atac_logFC > 0 & class_expr == "down",gene_expr< -1) %>%
  mutate(bin_num = ntile(mean_atac_logFC, n2_bins))

# Add bin labels
bin_info_up <- concordant_up %>%
  group_by(bin_num) %>%
  summarise(
    min_atac = min(mean_atac_logFC),
    max_atac = max(mean_atac_logFC),
    .groups = "drop"
  ) %>%
  mutate(bin_label = sprintf("UP_bin%02d_atac_%.2f_to_%.2f", bin_num, min_atac, max_atac))

bin_info_up_dis <-disconcordant_up %>%
  group_by(bin_num) %>%
  summarise(
    min_atac = min(mean_atac_logFC),
    max_atac = max(mean_atac_logFC),
    .groups = "drop"
  ) %>%
  mutate(bin_label = sprintf("UP_bin%02d_atac%.2f_to_%.2f", bin_num, min_atac, max_atac))

bin_info_down <- concordant_down %>%
  group_by(bin_num) %>%
  summarise(
    min_atac = min(mean_atac_logFC),
    max_atac = max(mean_atac_logFC),
    .groups = "drop"
  ) %>%
  mutate(bin_label = sprintf("DOWN_bin%02d_atac_%.2f_to_%.2f", bin_num, min_atac, max_atac))

bin_info_down_dis <- disconcordant_down %>%
  group_by(bin_num) %>%
  summarise(
    min_atac = min(mean_atac_logFC),
    max_atac = max(mean_atac_logFC),
    .groups = "drop"
  ) %>%
  mutate(bin_label = sprintf("DOWN_bin%02d_atac_%.2f_to_%.2f", bin_num, min_atac, max_atac))

concordant_up <- concordant_up %>% left_join(bin_info_up, by = "bin_num")
concordant_down <- concordant_down %>% left_join(bin_info_down, by = "bin_num")
concordant_up_dis <- disconcordant_up %>% left_join(bin_info_up_dis, by = "bin_num")
concordant_down_dis <- disconcordant_down %>% left_join(bin_info_down_dis, by = "bin_num")

# STEP 3: Create gene lists for GSEA

# --- CONCORDANT UP GENE LISTS ---
concordant_up_genelist <- split(concordant_up$closest_gene_dox, concordant_up$bin_label)

# --- CONCORDANT DOWN GENE LISTS ---
concordant_down_genelist <- split(concordant_down$closest_gene_dox, concordant_down$bin_label)
# --- dis-CONCORDANT UP GENE LISTS ---
disconcordant_up_genelist <- split(concordant_up_dis$closest_gene_dox, concordant_up_dis$bin_label)

# --- dis-CONCORDANT DOWN GENE LISTS ---
disconcordant_down_genelist <- split(concordant_down_dis$closest_gene_dox, concordant_down_dis$bin_label)

# --- COMBINED LIST ---
df_up<-filter(df_promoter, 
              logFC_expr>1, 
              class_expr=="up")
df_down<-filter(df_promoter, 
                logFC_expr< -1, 
                class_expr=="down")


up_ma9ON<-df_up$name
down_ma9ON<-df_down$name


list<-list(
  "up_ma9ON"=up_ma9ON, 
  "down_ma9ON"=down_ma9ON

)


all_genelists <- c(list,concordant_up_genelist, disconcordant_up_genelist,
                   concordant_down_genelist,disconcordant_down_genelist
)
dplyr::glimpse(all_genelists)

# run geneset enrichment

Mm <- bind_rows(
  msigdbr(species = "Mus musculus") %>%
    filter(gs_collection == "C2") %>%
    dplyr::select(gs_name, ncbi_gene, gene_symbol),
  
  msigdbr(species = "Mus musculus") %>%
    filter(gs_collection == "C5", gs_subcollection == "GO:BP") %>%
    dplyr::select(gs_name, ncbi_gene, gene_symbol),
  
  msigdbr(species = "Mus musculus") %>% 
    filter(gs_collection == "H") %>%
    dplyr::select(gs_name, ncbi_gene, gene_symbol)
)

# RUN
res <- compareCluster(all_genelists, enricher, TERM2GENE=Mm[,c(1,3)])

# ploting
df <- as.data.frame(res@compareClusterResult)

# --- Source data for reviewers: full compareCluster ORA table (p.adjust < 0.05) ---
ora_out <- as.data.frame(res@compareClusterResult)
ora_out <- ora_out[!is.na(ora_out$p.adjust) & ora_out$p.adjust < 0.05, ]
write.csv(ora_out, "Fig2D_S2B_ORA.csv", row.names = FALSE)

N <- 10  # top N (20) per cluster to define the union (set 20 if you can tolerate more rows)

# 1) top N per cluster (ONLY to define union, not to cut the universe globally)
df_topN <- df %>%
  filter(!is.na(Description)) %>%
  group_by(Cluster) %>%
  arrange(p.adjust, desc(Count)) %>%
  slice_head(n = N) %>%
  ungroup()

# 2) union of terms across clusters
union_terms <- unique(df_topN$Description)

# 3) keep ALL results for those union terms (so terms appear across clusters if present)
df_union <- df %>%
  filter(Description %in% union_terms) %>%
  filter(!is.na(Description))

# 4) plot object: copy + replace slot (do not overwrite the original)
res_plot <- res
res_plot@compareClusterResult <- df_union

# 5) dotplot provides the clean layout/order
p <- dotplot(res_plot, showCategory = length(union_terms))

# 6) extract dotplot data and apply your styling
plot_data <- as.data.frame(p$data)
plot_data$log10_padjust <- -log10(plot_data$p.adjust)

ggplot(plot_data, aes(x = Cluster, y = Description, size = Count)) +
  geom_point(shape = 21, aes(fill = log10_padjust), color = "black", stroke = 0.5) +
  scale_fill_gradient(low = "#9ECAE1", high = "#08519C",, name = "-log10(p.adjust)") +
  theme_bw() +
  labs(title = "bin2.ccdt_bin1.dist_6groups0.3_atac_H_C2_BP_top10",
       x = "Cluster", y = "Pathway") +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 8),
        axis.text.y = element_text(size = 6)) +
  scale_size_continuous(name = "Gene Count", range = c(3, 10), guide = "legend")


# Figure_S2B ----

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# 1. Load compareCluster results

df <- as.data.frame(res@compareClusterResult)

# Check clusters
cat("Clusters found in data:\n")
print(unique(df$Cluster))

# 2. Define cluster order and display labels

# Your actual cluster names from compareCluster
cluster_order <- c(
  "up_ma9ON",
  "down_arid1bKD",
  "UP_bin02_atac_0.43_to_3.41",
  "UP_bin01_atac_0.00_to_0.43",
  #"DOWN_bin01_atac_0.01_to_2.63",
  
  "down_ma9ON",
  "up_arid1bKD",
  "DOWN_bin01_atac_-3.64_to_-0.57",
  "DOWN_bin02_atac_-0.57_to_-0.00"
  #"UP_bin01_atac_-2.21_to_-0.00"
)

# Display labels
cluster_labels <- c(
  "up_ma9ON" = "RNA up",
  
  "down_ma9ON" = "RNA down",
  "UP_bin02_atac_0.43_to_3.41" = "ATAC++ RNA up",
  "UP_bin01_atac_0.00_to_0.43" = "ATAC + RNA up",
  #"DOWN_bin01_atac_0.01_to_2.63" = "Enh+\nRNA Down",
  
  "DOWN_bin01_atac_-3.64_to_-0.57" = "ATAC-- RNA down",
  "DOWN_bin02_atac_-0.57_to_-0.00" = "ATAC- RNA down",
  #"UP_bin01_atac_-2.21_to_-0.00" = "Enh-\nRNA Up",
  
  "up_arid1bKD"="RNA up (Arid1b KD)",
  "down_arid1bKD"="RNA down (Arid1b KD)"
)


# Keep only clusters that exist in data
cluster_order <- cluster_order[cluster_order %in% unique(df$Cluster)]
cat("\nClusters to plot:\n")
print(cluster_order)

# 3. Define selected gene sets

selected_genesets <- tribble(
  ~geneset_id, ~display_name, ~theme,
  
  # 1) LEUKEMIA (Top significant in Upregulated) 
  "HUANG_AML_LSC47", "LSC47 (Huang)", "Leukemia",
  "KEGG_MEDICUS_VARIANT_MLL_ENL_FUSION_TO_TRANSCRIPTIONAL_ACTIVATION", "MLL-ENL fusion targets (KEGG)", "Leukemia",
  "WANG_IMMORTALIZED_BY_HOXA9_AND_MEIS1_DN", "HOXA9/MEIS1 vs HOXA9 UP (Wang)", "Leukemia",
  "HESS_TARGETS_OF_HOXA9_AND_MEIS1_UP", "HOXA9/MEIS1 targets UP (Hess)", "Leukemia",
  # ESC-like stemness programs (co-opted by cancer/leukemia)
  #"WONG_EMBRYONIC_STEM_CELL_CORE", "Embryonic stem cell core (Wong)", "Leukemia",
  #"BHATTACHARYA_EMBRYONIC_STEM_CELL", "Embryonic stem cell (Bhattacharya)", "Leukemia",
  #"RAMALHO_STEMNESS_UP", "Stemness UP (Ramalho)", "Leukemia",
  
  # 2) HEMATOPOIESIS 
  "JAATINEN_HEMATOPOIETIC_STEM_CELL_UP", "HSC UP (Jaatinen)", "Hematopoiesis",
  # Complete Ivanova hierarchy (stem cell to mature)
  #"IVANOVA_HEMATOPOIESIS_STEM_CELL", "Stem cell (Ivanova)", "Hematopoiesis",
  "IVANOVA_HEMATOPOIESIS_STEM_CELL_LONG_TERM", "LT-HSC (Ivanova)", "Hematopoiesis",
  #"IVANOVA_HEMATOPOIESIS_STEM_CELL_SHORT_TERM", "ST-HSC (Ivanova)", "Hematopoiesis",
  "IVANOVA_HEMATOPOIESIS_EARLY_PROGENITOR", "Early progenitor (Ivanova)", "Hematopoiesis",
  #"IVANOVA_HEMATOPOIESIS_INTERMEDIATE_PROGENITOR", "Intermediate progenitor (Ivanova)", "Hematopoiesis",
  "IVANOVA_HEMATOPOIESIS_MATURE_CELL", "Mature cell (Ivanova)", "Hematopoiesis",
  # Other hematopoiesis CGP sets
  "BROWN_MYELOID_CELL_DEVELOPMENT_UP", "Myeloid development UP (Brown)", "Hematopoiesis",
  "KAMIKUBO_MYELOID_CEBPA_NETWORK", "Myeloid CEBPA network (Kamikubo)", "Hematopoiesis",
  # Lymphoid markers (to show myeloid vs lymphoid)
  #"LEE_EARLY_T_LYMPHOCYTE_UP", "Early T lymphocyte UP (Lee)", "Hematopoiesis",
  #"MORI_MATURE_B_LYMPHOCYTE_UP", "Mature B lymphocyte UP (Mori)", "Hematopoiesis",
  
  # 3) BIOLOGICAL PATHWAYS (GO:BP) 
  # Upregulated - Development/Morphogenesis
  "GOBP_EMBRYONIC_ORGAN_DEVELOPMENT", "Embryonic organ development", "Biological pathways",
  "GOBP_RENAL_SYSTEM_DEVELOPMENT", "Skeletal system development", "Biological pathways",
  "GOBP_SENSORY_ORGAN_MORPHOGENESIS","Sensory organ morphogenesis","Biological pathways",
  "GOBP_SKELETAL_SYSTEM_MORPHOGENESIS", "Renal system morphogenesis", "Biological pathways",
  # Upregulated - Metabolism
  "GOBP_AMINO_ACID_METABOLIC_PROCESS", "Amino acid metabolism", "Biological pathways",
  "GOBP_ALPHA_AMINO_ACID_METABOLIC_PROCESS", "Alpha amino acid metabolic process", "Biological pathways",
  "GOBP_ORGANIC_ACID_BIOSYNTHETIC_PROCESS", "Organic acid biosynthesis", "Biological pathways",
  
  # Downregulated (immune & migration)
  "GOBP_CELL_CHEMOTAXIS", "Cell chemotaxis", "Biological pathways",
  "GOBP_LEUKOCYTE_MIGRATION", "Leukocyte migration", "Biological pathways",
  "GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE", "Inflammatory response regulation", "Biological pathways",
  "GOBP_PHAGOCYTOSIS", "Phagocytosis", "Biological pathways",
  "GOBP_REGULATION_OF_LEUKOCYTE_PROLIFERATION", "Regulation of leukocyte proliferation", "Biological pathways",
  
  
  # 4) HALLMARK - Top significant each direction 
  # Upregulated
  "HALLMARK_MYC_TARGETS_V2", "MYC targets V2", "Hallmark",
  "HALLMARK_MTORC1_SIGNALING", "mTORC1 signaling", "Hallmark",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "Unfolded protein response", "Hallmark",
  # Downregulated
  "HALLMARK_INFLAMMATORY_RESPONSE", "Inflammatory response", "Hallmark",
  "HALLMARK_IL2_STAT5_SIGNALING", "IL2/STAT5 signaling", "Hallmark",
  "HALLMARK_KRAS_SIGNALING_UP", "KRAS signaling UP", "Hallmark",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNF/NF-kB signaling", "Hallmark",
  "HALLMARK_COMPLEMENT", "Complement", "Hallmark"
)

# 4. Prepare data for visualization
plot_data <- df %>%
  filter(ID %in% selected_genesets$geneset_id) %>%
  filter(Cluster %in% cluster_order) %>%
  left_join(selected_genesets, by = c("ID" = "geneset_id")) %>%
  mutate(
    log10_padj = -log10(p.adjust),
    log10_padj_capped = pmin(log10_padj, 10),
    sig = case_when(
      p.adjust < 0.001 ~ "***",
      p.adjust < 0.01 ~ "**",
      p.adjust < 0.05 ~ "*",
      TRUE ~ ""
    ),
    theme = factor(theme, levels = c("Leukemia", "Hematopoiesis", "Biological pathways", "Hallmark")),
    Cluster = factor(Cluster, levels = cluster_order)
  )

# Order gene sets by theme
geneset_order <- selected_genesets %>%
  mutate(theme = factor(theme, levels = c("Leukemia", "Hematopoiesis", "Biological pathways", "Hallmark"))) %>%
  arrange(theme) %>%
  pull(display_name)

plot_data <- plot_data %>%
  mutate(display_name = factor(display_name, levels = rev(geneset_order)))

# 5. DOTPLOT

# Create scaffold
scaffold <- expand.grid(
  display_name = geneset_order,
  Cluster = cluster_order,
  stringsAsFactors = FALSE
) %>%
  left_join(selected_genesets %>% dplyr::select(display_name, theme), by = "display_name") %>%
  mutate(
    theme = factor(theme, levels = c("Leukemia", "Hematopoiesis", "Biological pathways", "Hallmark")),
    display_name = factor(display_name, levels = rev(geneset_order)),
    Cluster = factor(Cluster, levels = cluster_order)
  )

# Filter significant results
plot_data_sig <- plot_data %>%
  filter(p.adjust < 0.05)
# dotplot
Figure_S2B<-ggplot() +
  geom_point(data = scaffold, aes(x = Cluster, y = display_name), alpha = 0) +
  geom_point(data = plot_data_sig, aes(x = Cluster, y = display_name, size = Count, fill = log10_padj_capped), 
             shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient(
    low = "white", high = "darkgreen",
    limits = c(0, 10),
    name = "-log10\n(p.adjust)"
  ) +
  scale_size_continuous(range = c(2, 10), name = "Gene count") +
  scale_x_discrete(labels = cluster_labels[cluster_order]) +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", face = "bold", angle = 45, hjust = 1),
    axis.title = element_blank(),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 10),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(Figure_S2B)
# Figure_2D  ----

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# 1. Load compareCluster results

df <- as.data.frame(res@compareClusterResult)

# Check clusters
cat("Clusters found in data:\n")
print(unique(df$Cluster))

# 2. Define cluster order and display labels

# Your actual cluster names from compareCluster
cluster_order <- c(
  
  "UP_bin02_atac_0.43_to_3.41",
  "UP_bin01_atac_0.00_to_0.43",
  #"DOWN_bin01_atac_0.01_to_2.63",
 
  "DOWN_bin01_atac_-3.64_to_-0.57",
  "DOWN_bin02_atac_-0.57_to_-0.00"
  #"UP_bin01_atac_-2.21_to_-0.00"
)

# Display labels
cluster_labels <- c(

  "UP_bin02_atac_0.43_to_3.41" = "ATAC++ RNA up",
  "UP_bin01_atac_0.00_to_0.43" = "ATAC + RNA up",

  "DOWN_bin01_atac_-3.64_to_-0.57" = "ATAC-- RNA down",
  "DOWN_bin02_atac_-0.57_to_-0.00" = "ATAC- RNA down"

)

# Keep only clusters that exist in data
cluster_order <- cluster_order[cluster_order %in% unique(df$Cluster)]
cat("\nClusters to plot:\n")
print(cluster_order)

# 3. Define selected gene sets

selected_genesets <- tribble(
  ~geneset_id, ~display_name, ~theme,
  
  # 1) LEUKEMIA (Top significant in Upregulated) 
  "HUANG_AML_LSC47", "LSC47 (Huang)", "Leukemia",
  "KEGG_MEDICUS_VARIANT_MLL_ENL_FUSION_TO_TRANSCRIPTIONAL_ACTIVATION", "MLL-ENL fusion targets (KEGG)", "Leukemia",
  "WANG_IMMORTALIZED_BY_HOXA9_AND_MEIS1_DN", "HOXA9/MEIS1 vs HOXA9 UP (Wang)", "Leukemia",
  "HESS_TARGETS_OF_HOXA9_AND_MEIS1_UP", "HOXA9/MEIS1 targets UP (Hess)", "Leukemia",
  # ESC-like stemness programs (co-opted by cancer/leukemia)
  #"WONG_EMBRYONIC_STEM_CELL_CORE", "Embryonic stem cell core (Wong)", "Leukemia",
  #"BHATTACHARYA_EMBRYONIC_STEM_CELL", "Embryonic stem cell (Bhattacharya)", "Leukemia",
  #"RAMALHO_STEMNESS_UP", "Stemness UP (Ramalho)", "Leukemia",
  
  # 2) HEMATOPOIESIS 
  "JAATINEN_HEMATOPOIETIC_STEM_CELL_UP", "HSC UP (Jaatinen)", "Hematopoiesis",
  # Complete Ivanova hierarchy (stem cell to mature)
  #"IVANOVA_HEMATOPOIESIS_STEM_CELL", "Stem cell (Ivanova)", "Hematopoiesis",
  "IVANOVA_HEMATOPOIESIS_STEM_CELL_LONG_TERM", "LT-HSC (Ivanova)", "Hematopoiesis",
  #"IVANOVA_HEMATOPOIESIS_STEM_CELL_SHORT_TERM", "ST-HSC (Ivanova)", "Hematopoiesis",
  "IVANOVA_HEMATOPOIESIS_EARLY_PROGENITOR", "Early progenitor (Ivanova)", "Hematopoiesis",
  #"IVANOVA_HEMATOPOIESIS_INTERMEDIATE_PROGENITOR", "Intermediate progenitor (Ivanova)", "Hematopoiesis",
  "IVANOVA_HEMATOPOIESIS_MATURE_CELL", "Mature cell (Ivanova)", "Hematopoiesis",
  # Other hematopoiesis CGP sets
  "BROWN_MYELOID_CELL_DEVELOPMENT_UP", "Myeloid development UP (Brown)", "Hematopoiesis",
  "KAMIKUBO_MYELOID_CEBPA_NETWORK", "Myeloid CEBPA network (Kamikubo)", "Hematopoiesis",
  # Lymphoid markers (to show myeloid vs lymphoid)
  #"LEE_EARLY_T_LYMPHOCYTE_UP", "Early T lymphocyte UP (Lee)", "Hematopoiesis",
  #"MORI_MATURE_B_LYMPHOCYTE_UP", "Mature B lymphocyte UP (Mori)", "Hematopoiesis",
  
  # 3) BIOLOGICAL PATHWAYS (GO:BP) 
  # Upregulated - Development/Morphogenesis
  "GOBP_EMBRYONIC_ORGAN_DEVELOPMENT", "Embryonic organ development", "Biological pathways",
  "GOBP_RENAL_SYSTEM_DEVELOPMENT", "Skeletal system development", "Biological pathways",
  "GOBP_SENSORY_ORGAN_MORPHOGENESIS","Sensory organ morphogenesis","Biological pathways",
  "GOBP_SKELETAL_SYSTEM_MORPHOGENESIS", "Renal system morphogenesis", "Biological pathways",
  # Upregulated - Metabolism
  "GOBP_AMINO_ACID_METABOLIC_PROCESS", "Amino acid metabolism", "Biological pathways",
  "GOBP_ALPHA_AMINO_ACID_METABOLIC_PROCESS", "Alpha amino acid metabolic process", "Biological pathways",
  "GOBP_ORGANIC_ACID_BIOSYNTHETIC_PROCESS", "Organic acid biosynthesis", "Biological pathways",
  
  # Downregulated (immune & migration)
  "GOBP_CELL_CHEMOTAXIS", "Cell chemotaxis", "Biological pathways",
  "GOBP_LEUKOCYTE_MIGRATION", "Leukocyte migration", "Biological pathways",
  "GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE", "Inflammatory response regulation", "Biological pathways",
  "GOBP_PHAGOCYTOSIS", "Phagocytosis", "Biological pathways",
  "GOBP_REGULATION_OF_LEUKOCYTE_PROLIFERATION", "Regulation of leukocyte proliferation", "Biological pathways",
  
  
  # 4) HALLMARK - Top significant each direction 
  # Upregulated
  "HALLMARK_MYC_TARGETS_V2", "MYC targets V2", "Hallmark",
  "HALLMARK_MTORC1_SIGNALING", "mTORC1 signaling", "Hallmark",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "Unfolded protein response", "Hallmark",
  # Downregulated
  "HALLMARK_INFLAMMATORY_RESPONSE", "Inflammatory response", "Hallmark",
  "HALLMARK_IL2_STAT5_SIGNALING", "IL2/STAT5 signaling", "Hallmark",
  "HALLMARK_KRAS_SIGNALING_UP", "KRAS signaling UP", "Hallmark",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNF/NF-kB signaling", "Hallmark",
  "HALLMARK_COMPLEMENT", "Complement", "Hallmark"
)

# 4. Prepare data for visualization
plot_data <- df %>%
  filter(ID %in% selected_genesets$geneset_id) %>%
  filter(Cluster %in% cluster_order) %>%
  left_join(selected_genesets, by = c("ID" = "geneset_id")) %>%
  mutate(
    log10_padj = -log10(p.adjust),
    log10_padj_capped = pmin(log10_padj, 10),
    sig = case_when(
      p.adjust < 0.001 ~ "***",
      p.adjust < 0.01 ~ "**",
      p.adjust < 0.05 ~ "*",
      TRUE ~ ""
    ),
    theme = factor(theme, levels = c("Leukemia", "Hematopoiesis", "Biological pathways", "Hallmark")),
    Cluster = factor(Cluster, levels = cluster_order)
  )

# Order gene sets by theme
geneset_order <- selected_genesets %>%
  mutate(theme = factor(theme, levels = c("Leukemia", "Hematopoiesis", "Biological pathways", "Hallmark"))) %>%
  arrange(theme) %>%
  pull(display_name)

plot_data <- plot_data %>%
  mutate(display_name = factor(display_name, levels = rev(geneset_order)))

# 5. DOTPLOT

# Create scaffold
scaffold <- expand.grid(
  display_name = geneset_order,
  Cluster = cluster_order,
  stringsAsFactors = FALSE
) %>%
  left_join(selected_genesets %>% dplyr::select(display_name, theme), by = "display_name") %>%
  mutate(
    theme = factor(theme, levels = c("Leukemia", "Hematopoiesis", "Biological pathways", "Hallmark")),
    display_name = factor(display_name, levels = rev(geneset_order)),
    Cluster = factor(Cluster, levels = cluster_order)
  )

# Filter significant results
plot_data_sig <- plot_data %>%
  filter(p.adjust < 0.05)
# dotplot
Figure_2D<-ggplot() +
  geom_point(data = scaffold, aes(x = Cluster, y = display_name), alpha = 0) +
  geom_point(data = plot_data_sig, aes(x = Cluster, y = display_name, size = Count, fill = log10_padj_capped), 
             shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient(
    low = "white", high = "darkgreen",
    limits = c(0, 10),
    name = "-log10\n(p.adjust)"
  ) +
  scale_size_continuous(range = c(2, 10), name = "Gene count") +
  scale_x_discrete(labels = cluster_labels[cluster_order]) +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", face = "bold", angle = 45, hjust = 1),
    axis.title = element_blank(),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 10),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(Figure_2D)
############################# Figure 2E ############################# ----
# Proportion of E-P loops (PCHi-C) by pE class, MA9-on / MA9-off.
#
# pE classes come straight from DF_ENHANCER$class_activity (values a_2 ... o_-2).
# This replaces the old rank-binned `atac_group` column; the two are the same
# grouping (confirmed block-diagonal), just named differently.
#
# Requires DF_ENHANCER already loaded, with columns:
#   logFC_atac, closest_gene_dox, class_activity, interaction_class

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(patchwork)

df <- DF_ENHANCER

n_bins <- 50

# STEP 1: Reconstruct gene-level ATAC bins (same as figure 2C)
gene_stats <- df %>%
  group_by(closest_gene_dox) %>%
  summarise(
    mean_atac_logFC = mean(logFC_atac, na.rm = TRUE),
    n_enhancers     = n(),
    .groups = "drop"
  )

neg_genes <- gene_stats %>% filter(mean_atac_logFC <  0) %>%
  mutate(atac_bin = ntile(mean_atac_logFC, 2))
pos_genes <- gene_stats %>% filter(mean_atac_logFC >= 0) %>%
  mutate(atac_bin = ntile(mean_atac_logFC, 2) + 2L)

gene_stats <- bind_rows(neg_genes, pos_genes) %>%
  mutate(
    bin_num = ntile(mean_atac_logFC, n_bins),
    atac_bin_label = factor(atac_bin,
                            levels = 1:4,
                            labels = c("ATAC\u2212\u2212", "ATAC\u2212", "ATAC+", "ATAC++"))
  )

# Bin label info for x-axis
bin_info <- gene_stats %>%
  group_by(bin_num) %>%
  summarise(
    mean_atac = mean(mean_atac_logFC),
    min_atac  = min(mean_atac_logFC),
    max_atac  = max(mean_atac_logFC),
    .groups   = "drop"
  ) %>%
  mutate(bin_label = sprintf("[%.2f,%.2f)", min_atac, max_atac))

# Dashed line positions and group label x positions (same as figure 2C)
boundary_bins <- gene_stats %>%
  group_by(atac_bin) %>%
  summarise(max_bin = max(bin_num), .groups = "drop") %>%
  arrange(atac_bin)
dash_positions <- boundary_bins$max_bin[-nrow(boundary_bins)] + 0.5

group_label_x <- gene_stats %>%
  group_by(atac_bin) %>%
  summarise(mid_bin = mean(range(bin_num)), .groups = "drop") %>%
  arrange(atac_bin) %>%
  pull(mid_bin)
group_labels <- c("ATAC\u2212\u2212", "ATAC\u2212", "ATAC+", "ATAC++")

df_bins <- df %>%
  left_join(gene_stats %>% select(closest_gene_dox, atac_bin, atac_bin_label, bin_num),
            by = "closest_gene_dox") %>%
  filter(!is.na(atac_bin))

# pE class levels (a-o, ordered by ATAC logFC).
# Levels are taken straight from the column so they match byte-for-byte.
# (class_activity carries a non-ASCII character, so a hand-typed list of
#  "a_2","b_1",... will silently fail to match and turn every row into NA.)
pe_levels <- sort(unique(as.character(df_bins$class_activity)))
pe_labels <- c("a (2)","b (1)","c (1)","d (1)","e (0)",
               "f (0)","g (0)","h (0)","i (-0)",
               "j (-0)","k (-0)","l (-1)","m (-1)","n (-1)","o (-2)")

stopifnot(length(pe_levels) == length(pe_labels))   # must be 15 classes, a -> o

df_bins <- df_bins %>%
  mutate(
    atac_group_f = factor(class_activity, levels = pe_levels, labels = pe_labels)
  )

stopifnot(sum(is.na(df_bins$atac_group_f)) == 0)     # every row must map to a class

bin_order <- c("ATAC\u2212\u2212", "ATAC\u2212", "ATAC+", "ATAC++")

# STEP 3: Number of E-P loops per gene (boxplot)
# Mirrors top panel of figure 2E
# Count loops per gene per condition
loops_per_gene <- df_bins %>%
  group_by(closest_gene_dox, atac_bin_label) %>%
  summarise(
    n_loops_on  = sum(interaction_class %in% c("HiC_HiC", "HiC_nonHiC")),
    n_loops_off = sum(interaction_class %in% c("HiC_HiC", "nonHiC_HiC")),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(n_loops_on, n_loops_off),
               names_to = "condition", values_to = "n_loops") %>%
  mutate(
    condition = ifelse(condition == "n_loops_on", "MA9-on", "MA9-off"),
    condition = factor(condition, levels = c("MA9-on", "MA9-off")),
    atac_bin_label = factor(atac_bin_label, levels = bin_order)
  )

# Count loops by atac_bin x pE_class x condition
loop_pe_on <- df_bins %>%
  filter(interaction_class %in% c("HiC_HiC", "HiC_nonHiC")) %>%
  group_by(atac_bin_label, atac_group_f) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(atac_bin_label) %>%
  mutate(prop = n / sum(n),
         condition = "MA9-on")

loop_pe_off <- df_bins %>%
  filter(interaction_class %in% c("HiC_HiC", "nonHiC_HiC")) %>%
  group_by(atac_bin_label, atac_group_f) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(atac_bin_label) %>%
  mutate(prop = n / sum(n),
         condition = "MA9-off")

loop_pe <- bind_rows(loop_pe_on, loop_pe_off) %>%
  mutate(
    atac_bin_label = factor(atac_bin_label, levels = bin_order),
    condition      = factor(condition, levels = c("MA9-on", "MA9-off")),
    x_label        = interaction(condition, atac_bin_label, sep = "\n")
  )

# Build x-axis order: MA9-on ATAC-- | MA9-off ATAC-- | ...
x_levels <- as.vector(outer(c("MA9-on", "MA9-off"), bin_order,
                            FUN = function(a, b) paste(a, b, sep = "\n")))
loop_pe <- loop_pe %>%
  mutate(x_label = factor(x_label, levels = x_levels))

# pE class color: warm for gain, grey for neutral, cool for loss
pe_colors <- c(
  colorRampPalette(c("#67000d", "#fc4e2a", "#feb24c", "#ffeda0"))(4),  # a-d gain
  colorRampPalette(c("#f7f7f7", "#cccccc", "#969696", "#636363", "#252525", "#bdbdbd", "#737373"))(7), # e-k neutral
  colorRampPalette(c("#c6dbef", "#6baed6", "#2171b5", "#084594"))(4)   # l-o loss
)
names(pe_colors) <- pe_labels

ggplot(loop_pe,
       aes(x = x_label, y = atac_group_f, fill = prop)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = n), size = 2, color = "black") +
  scale_fill_gradient(low = "white", high = "#08306b",
                      name = "Proportion") +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = seq(2.5, length(x_levels) - 0.5, by = 2),
             linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = c(4.5, 11.5), linetype = "dashed",
             color = "black", linewidth = 0.4) +
  labs(x = NULL, y = "pE class (log2FC)",
       title = "Proportion of E-P loops by pE class") +
  theme_pubr(base_size = 10) +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y    = element_text(size = 8),
    legend.position = "right"
  )

# STEP 4b: 50-bin heatmaps - MA9-on, MA9-off, and difference
loop_pe_50_on <- df_bins %>%
  filter(interaction_class %in% c("HiC_HiC", "HiC_nonHiC")) %>%
  group_by(bin_num, atac_group_f) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(bin_num) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

loop_pe_50_off <- df_bins %>%
  filter(interaction_class %in% c("HiC_HiC", "nonHiC_HiC")) %>%
  group_by(bin_num, atac_group_f) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(bin_num) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Difference: MA9-on minus MA9-off (fill missing combos with 0)
all_combos <- expand.grid(
  bin_num      = 1:n_bins,
  atac_group_f = unique(df_bins$atac_group_f)
)
loop_pe_50_diff <- all_combos %>%
  left_join(loop_pe_50_on  %>% select(bin_num, atac_group_f, prop),
            by = c("bin_num", "atac_group_f")) %>%
  rename(prop_on = prop) %>%
  left_join(loop_pe_50_off %>% select(bin_num, atac_group_f, prop),
            by = c("bin_num", "atac_group_f")) %>%
  rename(prop_off = prop) %>%
  mutate(
    prop_on  = replace_na(prop_on,  0),
    prop_off = replace_na(prop_off, 0),
    diff     = prop_on - prop_off
  )

label_at  <- c(1, seq(5, n_bins, by = 5), n_bins)
label_txt <- bin_info$bin_label[label_at]
vlines    <- boundary_bins$max_bin[-nrow(boundary_bins)] + 0.5
abs_lim   <- max(abs(loop_pe_50_diff$diff), na.rm = TRUE)
prop_lim  <- max(c(loop_pe_50_on$prop, loop_pe_50_off$prop), na.rm = TRUE)

# Shared theme for all three panels
hm_theme <- list(
  scale_y_discrete(limits = rev),
  scale_x_continuous(breaks = label_at, labels = label_txt, expand = c(0, 0)),
  geom_vline(xintercept = vlines, linetype = "dashed",
             color = "black", linewidth = 0.5),
  geom_hline(yintercept = c(4.5, 11.5), linetype = "dashed",
             color = "black", linewidth = 0.4),
  theme_classic(base_size = 10),
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y     = element_text(size = 8),
    legend.position = "right",
    plot.margin     = margin(t = 5, r = 10, b = 5, l = 10)
  )
)

Figure_2E<-ggplot(loop_pe_50_on,
       aes(x = bin_num, y = atac_group_f, fill = prop)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "white", high = "#08306b",
                      limits = c(0, prop_lim), name = "Proportion") +
  annotate("text", x = group_label_x, y = 16.2, label = group_labels,
           size = 3, fontface = "bold", hjust = 0.5) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "pE class", title = "MA9-on") +
  hm_theme

print(Figure_2E)
############################# Figure S2A ############################## ----
## Figure 3: Fisher's exact test heatmap (log2 OR vs neutral)
## Three enhancer marks: ATAC, H3K4me1, H3K27ac
## E-P linking: closest_gene_dox (distance-based, no HiC)
## Groups: core, up, neutral, down (MLL + nonMLL merged, core separate)
## Reference: neutral | Plot columns: core / up / down (neutral excluded)

library(tidyverse)
library(patchwork)

# 2. Define groups: core separate, MLL+nonMLL merged by class_expr 
df_promoter_02 <- df_promoter %>%
  mutate(
    group_fc1 = case_when(
      interest_group == "MLL_target_interest" ~ "core",
      TRUE ~ class_expr
    )
  )

group_levels <- c("core", "up", "neutral", "down")
gene_lookup  <- df_promoter_02 %>% dplyr::select(name, group_fc1)

# 3. Deduplicate enhancers 
df_enh <- DF_ENHANCER %>% distinct(name_dhs_dox, .keep_all = TRUE)

# 4. Create 15-bin groups for each mark 
create_mark_bins <- function(df, logfc_col, group_col) {
  df[[group_col]] <- cut(
    rank(-df[[logfc_col]], ties.method = "first"),
    breaks = quantile(rank(-df[[logfc_col]], ties.method = "first"),
                      probs = seq(0, 1, length.out = 16)),
    labels = paste0(letters[1:15], "_class"),
    include.lowest = TRUE
  )
  return(df)
}

df_enh <- create_mark_bins(df_enh, "logFC_atac",    "atac_bin")
df_enh <- create_mark_bins(df_enh, "logFC_h3k27ac", "h3k27ac_bin")

# 5. Distance-based E-P linking 
df_ep <- df_enh %>%
  filter(!is.na(closest_gene_dox)) %>%
  left_join(gene_lookup, by = c("closest_gene_dox" = "name")) %>%
  filter(!is.na(group_fc1)) %>%
  mutate(group_fc1 = factor(group_fc1, levels = group_levels))

cat("Distance E-P pairs:", nrow(df_ep), "\n")
print(table(df_ep$group_fc1))

# 6. Fisher's exact: log2(OR) vs neutral for each bin x group 
bin_levels   <- paste0(letters[1:15], "_class")
test_groups  <- c("core", "up", "down")   # neutral is reference, not tested

run_fisher <- function(df_ep, mark_col) {
  results <- list()
  
  for (bin in bin_levels) {
    in_bin <- df_ep[[mark_col]] == bin
    
    tab <- df_ep %>%
      mutate(in_bin = in_bin) %>%
      group_by(group_fc1) %>%
      summarise(in_b = sum(in_bin), out_b = sum(!in_bin), .groups = "drop")
    
    neutral_in  <- tab %>% filter(group_fc1 == "neutral") %>% pull(in_b)
    neutral_out <- tab %>% filter(group_fc1 == "neutral") %>% pull(out_b)
    
    for (grp in test_groups) {
      grp_row <- tab %>% filter(group_fc1 == grp)
      if (nrow(grp_row) == 0) next
      
      mat <- matrix(
        c(grp_row$in_b,  neutral_in,
          grp_row$out_b, neutral_out),
        nrow = 2,
        dimnames = list(c("in_bin", "out_bin"), c(grp, "neutral"))
      )
      
      ft <- fisher.test(mat)
      
      results[[length(results) + 1]] <- tibble(
        mark_bin = bin,
        group    = grp,
        log2OR   = log2(ft$estimate),
        pval     = ft$p.value
      )
    }
  }
  
  bind_rows(results) %>%
    group_by(group) %>%
    mutate(padj = p.adjust(pval, method = "BH")) %>%
    ungroup() %>%
    mutate(
      sig = case_when(
        padj < 0.001 ~ "***",
        padj < 0.01  ~ "**",
        padj < 0.05  ~ "*",
        TRUE         ~ ""
      ),
      mark_bin = factor(mark_bin, levels = bin_levels),
      group    = factor(group, levels = test_groups)
    )
}

cat("Running Fisher's tests...\n")
fisher_atac    <- run_fisher(df_ep, "atac_bin")
fisher_h3k27ac <- run_fisher(df_ep, "h3k27ac_bin")

cat("log2OR range ATAC:",    range(fisher_atac$log2OR[is.finite(fisher_atac$log2OR)]), "\n")
cat("log2OR range H3K27ac:", range(fisher_h3k27ac$log2OR[is.finite(fisher_h3k27ac$log2OR)]), "\n")

# 7. Shared symmetric color limits 
all_or   <- c(fisher_atac$log2OR,  fisher_h3k27ac$log2OR)
or_limit <- max(abs(all_or[is.finite(all_or)]), na.rm = TRUE)
or_limit <- round(or_limit + 0.1, 1)
cat("Color scale limit: ±", or_limit, "\n")

#  8. Add OR label 
add_or_label <- function(df) {
  df %>% mutate(or_label = paste0(round(log2OR, 2), "\n", sig))
}

fisher_atac    <- add_or_label(fisher_atac)
fisher_h3k27ac <- add_or_label(fisher_h3k27ac)

bin_labels_gain <- c("a (gain++)", letters[2:14], "o (loss++)")
bin_labels_atac <- c("a (open++)", letters[2:14], "o (close++)")

# 9. Plot functions 
base_heatmap <- function(df_fisher, bin_labels_vec, show_y, cell_label_col) {
  ggplot(df_fisher,
         aes(x = group,
             y = factor(mark_bin, levels = rev(bin_levels)),
             fill = log2OR)) +
    geom_tile(aes(color = log2OR), linewidth = 0.4) +
    geom_tile(color = "grey85", fill = NA, linewidth = 0.3) +
    geom_text(aes(label = .data[[cell_label_col]]), size = 2.8, vjust = 0.8,
              lineheight = 0.85, color = "black") +
    scale_x_discrete(limits = test_groups) +
    scale_y_discrete(labels = if (show_y) rev(bin_labels_vec) else rep("", 15)) +
    scale_fill_gradientn(
      colours  = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                   "#FFFFFF",
                   "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
      limits   = c(-or_limit, or_limit),
      na.value = "grey80",
      name     = "log2(OR)"
    ) +
    scale_color_gradientn(
      colours  = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                   "#FFFFFF",
                   "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
      limits   = c(-or_limit, or_limit),
      guide    = "none"
    ) +
    labs(x = NULL,
         y = if (show_y) "Enhancer bin\n(a=most gain/open → o=most loss/close)" else NULL) +
    theme_minimal() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y     = element_text(size = 11),
      axis.title.y    = element_text(size = 10),
      plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
      legend.position = "right",
      panel.grid      = element_blank()
    )
}

make_heatmap_stars <- function(df_fisher, title_str, bin_labels_vec, show_y = TRUE) {
  base_heatmap(df_fisher, bin_labels_vec, show_y, "sig") + ggtitle(title_str)
}

make_heatmap_or <- function(df_fisher, title_str, bin_labels_vec, show_y = TRUE) {
  base_heatmap(df_fisher, bin_labels_vec, show_y, "or_label") + ggtitle(title_str)
}

# 10. Build and save Version A (stars only) 
p_atac_s    <- make_heatmap_stars(fisher_atac,    "ATAC",    bin_labels_atac, show_y = TRUE)
p_h3k27ac_s <- make_heatmap_stars(fisher_h3k27ac, "H3K27ac", bin_labels_gain, show_y = FALSE)

Figure_S2A <- (p_atac_s |  p_h3k27ac_s) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Fisher's exact vs neutral, log2(OR)",
    subtitle = "Groups: core / up / down (MLL+nonMLL merged) | Reference: neutral | BH-adjusted: * <0.05, ** <0.01, *** <0.001",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0.5)
    )
  )

print (Figure_S2A)

# 11. Build and save Version B (OR value + stars) 
p_atac_or    <- make_heatmap_or(fisher_atac,    "ATAC",    bin_labels_atac, show_y = TRUE)
p_h3k27ac_or <- make_heatmap_or(fisher_h3k27ac, "H3K27ac", bin_labels_gain, show_y = FALSE)

p_or <- (p_atac_or | p_h3k27ac_or) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Fisher's exact vs neutral, log2(OR)",
    subtitle = "Groups: core / up / down (MLL+nonMLL merged) | Reference: neutral | BH-adjusted: * <0.05, ** <0.01, *** <0.001",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0.5)
    )
  )

print (p_or)

############################# Figure_2F ############################## ----
# Figure 2F - E-P loops per gene & pE-class composition, by MLL-target class
#   Gene up/neutral/down groups are (re)defined here from RNA logFC & FDR
#   (|log2FC| > LFC_CUT AND FDR < FDR_CUT), so the cutoff is explicit/editable.
#
# Inputs:
#   df_promoter_fc.csv  (name, interest_group, logFC_expr, FDR_expr)
#   (E-P loop table built in-script below from DF_ENHANCER; no external loop CSV needed)
# Output: figure2F.pdf / .png

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggpubr); library(patchwork)
})

# setwd("/path/to/data")

# ---- RNA classification cutoff (your change) 
LFC_CUT <- 1       # |log2FC| threshold
FDR_CUT <- 0.05    # FDR threshold

# ---- 1. Recompute gene groups from expression 
gene_group <- df_promoter %>%
  mutate(
    class_expr = case_when(
      logFC_expr >  LFC_CUT & FDR_expr < FDR_CUT ~ "up",
      logFC_expr < -LFC_CUT & FDR_expr < FDR_CUT ~ "down",
      TRUE                                       ~ "neutral"
    ),
    group = case_when(
      interest_group == "MLL_target_interest"  ~ "MLL_target_interest",  # core: not reclassified
      str_starts(interest_group, "MLL_target") ~ paste0("MLL_target_", class_expr),
      TRUE                                     ~ paste0("nonMLL_target_", class_expr)
    )
  ) %>%
  select(name, group)

# x-axis order + display labels
group_levels <- c("MLL_target_interest", "MLL_target_up", "MLL_target_neutral",
                  "MLL_target_down", "nonMLL_target_up", "nonMLL_target_neutral",
                  "nonMLL_target_down")
group_labels <- c("core", "up", "neutral", "down", "up", "neutral", "down")
names(group_labels) <- group_levels

col_on  <- "#E05C44"   # MA9-on  (DOX)
col_off <- "#5B9BD5"   # MA9-off (VEH)

# ---- 2. Loop table + attach new groups by gene 
# ---- Build the E-P loop table from DF_ENHANCER (replaces df_enh_comb_select_mod.csv) 
# Correct interaction_gene_class, then split DF_ENHANCER into the six E_P_loop categories.
df_enh <- DF_ENHANCER %>%
  mutate(interaction_gene_class = case_when(
    interaction_class %in% c("HiC_nonHiC", "nonHiC_HiC") ~ "differentGene",
    TRUE ~ interaction_gene_class
  ))

# MA9-on novel loops (DOX-specific)
df_dox_spec <- df_enh %>%
  filter(interaction_class == "HiC_nonHiC", interaction_gene_class == "differentGene") %>%
  mutate(hic_gene_dox = strsplit(as.character(hic_gene_dox), ", *")) %>%
  unnest(cols = hic_gene_dox) %>%
  mutate(hic_gene_veh = NA, E_P_loop = "MA9_on_novel")

# Reorganized in DOX (same enhancer, different gene)
df_dox_diff <- df_enh %>%
  filter(interaction_class == "HiC_HiC", interaction_gene_class == "differentGene") %>%
  mutate(hic_gene_dox = strsplit(as.character(hic_gene_dox), ", *")) %>%
  unnest(cols = hic_gene_dox) %>%
  mutate(hic_gene_veh = NA, E_P_loop = "Reorganized_dox")

# MA9-off novel loops (VEH-specific)  [unnest only hic_gene_veh]
df_veh_spec <- df_enh %>%
  filter(interaction_class == "nonHiC_HiC", interaction_gene_class == "differentGene") %>%
  mutate(hic_gene_veh = strsplit(as.character(hic_gene_veh), ", *")) %>%
  unnest(cols = hic_gene_veh) %>%
  mutate(hic_gene_dox = NA, E_P_loop = "MA9_off_novel")

# Reorganized in VEH  [unnest only hic_gene_veh]
df_veh_diff <- df_enh %>%
  filter(interaction_class == "HiC_HiC", interaction_gene_class == "differentGene") %>%
  mutate(hic_gene_veh = strsplit(as.character(hic_gene_veh), ", *")) %>%
  unnest(cols = hic_gene_veh) %>%
  mutate(hic_gene_dox = NA, E_P_loop = "Reorganized_veh")

# nonHiC enhancers (no loop)
df_enh_nonHiC <- df_enh %>%
  filter(interaction_class == "nonHiC_nonHiC") %>%
  mutate(hic_gene_dox = NA, hic_gene_veh = NA, E_P_loop = "none")

# Unchanged loops (same gene, present in VEH)  [unnest only hic_gene_veh]
df_enh_split_same <- df_enh %>%
  filter(interaction_class_veh == "HiC", interaction_gene_class == "sameGene") %>%
  mutate(hic_gene_veh = strsplit(as.character(hic_gene_veh), ", *")) %>%
  unnest(cols = hic_gene_veh) %>%
  mutate(E_P_loop = "Unchanged")

df_enh_comb <- bind_rows(df_dox_spec, df_dox_diff, df_veh_spec,
                         df_veh_diff, df_enh_nonHiC, df_enh_split_same)

df_enh_comb_select <- df_enh_comb %>%
  dplyr::select(dplyr::any_of(c(
    "name_dhs_dox", "interaction_class", "interaction_gene_class", "E_P_loop",
    "hic_gene_dox", "hic_gene_veh",
    "dist_to_hic_gene_tss_dox", "interaction_class_dox", "closest_gene_dox",
    "dist_to_closest_gene_tss_dox",
    "dist_to_hic_gene_tss_veh", "interaction_class_veh", "closest_gene_veh",
    "dist_to_closest_gene_tss_veh",
    "atac_group", "class_atac", "logFC_atac",
    "interaction_score_diff", "interaction_score_dox", "interaction_score_veh",
    "loop_length", "loop_sample", "annot.type", "CpG_island_length", "CpG_island",
    "novel_atac_peak", "class_activity", "class_atac_arid1bKD")))

df <- df_enh_comb_select   # Fig 2F input, built above from DF_ENHANCER (no external CSV)

df <- df %>%
  left_join(gene_group %>% transmute(hic_gene_dox = name, group_dox = group),
            by = "hic_gene_dox") %>%
  left_join(gene_group %>% transmute(hic_gene_veh = name, group_veh = group),
            by = "hic_gene_veh")

# TOP PANEL - number of E-P loops per gene, MA9-on vs MA9-off
df_dox <- df %>%
  filter(!is.na(hic_gene_dox), !is.na(group_dox)) %>%
  count(gene = hic_gene_dox, group = group_dox, name = "n_loop_dox")

df_veh <- df %>%
  filter(!is.na(hic_gene_veh), !is.na(group_veh)) %>%
  count(gene = hic_gene_veh, group = group_veh, name = "n_loop_veh")

df_long <- full_join(df_dox, df_veh, by = c("gene", "group")) %>%
  mutate(across(c(n_loop_dox, n_loop_veh), ~ replace_na(.x, 0))) %>%
  pivot_longer(c(n_loop_dox, n_loop_veh),
               names_to = "condition", values_to = "n_loop") %>%
  mutate(
    condition = factor(ifelse(condition == "n_loop_dox", "MA9-on", "MA9-off"),
                       levels = c("MA9-on", "MA9-off")),
    group     = factor(group, levels = group_levels)
  ) %>%
  filter(!is.na(group))

p_top <- ggplot(df_long, aes(group, n_loop, fill = condition)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.85,
               linewidth = 0.3) +
  stat_compare_means(aes(group = condition), method = "wilcox.test",
                     label = "p.format", label.y = 46, size = 2.8) +
  geom_vline(xintercept = 4.5, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  scale_fill_manual(values = c("MA9-on" = col_on, "MA9-off" = col_off), name = NULL) +
  scale_x_discrete(labels = group_labels) +
  labs(x = NULL, y = "Number of E-P loops per gene") +
  coord_cartesian(ylim = c(0, 50)) +
  theme_pubr(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

# BOTTOM PANEL - proportion of E-P loops by pE class (a-o), with counts
df_bottom <- df %>%
  filter(!E_P_loop %in% c("none", "MA9_off_novel", "Reorganized_veh"),
         !is.na(group_dox), !is.na(class_activity)) %>%
  mutate(group = factor(group_dox, levels = group_levels),
         pE    = factor(class_activity, levels = sort(unique(class_activity)))) %>%
  filter(!is.na(group))

pe_pal <- c("#B2182B","#D6604D","#F4A582","#FDDBC7",
            "#000000","#404040","#737373","#969696","#BDBDBD","#D9D9D9","#F0F0F0",
            "#D1E5F0","#92C5DE","#4393C3","#2166AC")
pe_txt <- c("white","white","black","black",
            "white","white","white","black","black","black","black",
            "black","black","white","white")
pe_lvls <- levels(df_bottom$pE)
if (length(pe_lvls) != 15)
  warning("class_activity has ", length(pe_lvls), " levels (expected 15): ",
          paste(pe_lvls, collapse = ", "))
names(pe_pal) <- pe_lvls
names(pe_txt) <- pe_lvls

p_bottom <- ggplot(df_bottom, aes(group, fill = pE)) +
  geom_bar(position = "fill", alpha = 0.9, width = 0.85) +
  geom_text(stat = "count", aes(label = after_stat(count), colour = pE),
            position = position_fill(vjust = 0.5), size = 2.2, show.legend = FALSE) +
  geom_vline(xintercept = 4.5, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  scale_fill_manual(values = pe_pal, name = "pE class") +
  scale_colour_manual(values = pe_txt, guide = "none") +
  scale_x_discrete(labels = group_labels) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion of E-P loops by pE class") +
  theme_pubr(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

# ---- group sizes used (handy for sanity-checking the new cutoff) 
cat("\n--- genes per group at LFC>", LFC_CUT, " & FDR<", FDR_CUT, " ---\n", sep = "")
print(gene_group %>% filter(group %in% group_levels) %>% count(group))

fig2F <- p_top / p_bottom + plot_layout(heights = c(1, 1.15))
print(fig2F)
############################# Figure S2E ############################## ----
# Paired lollipop: E-P loops per core gene (only genes WITH HiC loops)
library(tidyverse)

# --> 1. Read data ----
core_genes <- df_promoter %>%
  filter(interest_group == "MLL_target_interest") %>% pull(name)

df_enh <- df_selected_4 %>% distinct(name_dhs_dox, .keep_all = TRUE)

# --> 2. Count HiC loops per gene ----
ct_dox <- df_enh %>%
  filter(interaction_class_dox == "HiC") %>%
  mutate(gene = strsplit(as.character(hic_gene_dox), ",\\s*")) %>%
  unnest(cols = gene) %>%
  filter(gene %in% core_genes) %>%
  count(gene, name = "MA9_on")

ct_veh <- df_enh %>%
  filter(interaction_class_veh == "HiC") %>%
  mutate(gene = strsplit(as.character(hic_gene_veh), ",\\s*")) %>%
  unnest(cols = gene) %>%
  filter(gene %in% core_genes) %>%
  count(gene, name = "MA9_off")

df_loops <- tibble(gene = core_genes) %>%
  left_join(ct_dox, by = "gene") %>%
  left_join(ct_veh, by = "gene") %>%
  replace_na(list(MA9_on = 0, MA9_off = 0)) %>%
  mutate(diff = MA9_on - MA9_off)

# --> 3. Filter to only genes with at least one HiC loop ----
df_loops <- df_loops %>%
  filter(MA9_on > 0 | MA9_off > 0)

cat("Core genes with HiC:", nrow(df_loops), "out of", length(core_genes), "\n")

# Order by MA9_on (descending)
df_loops <- df_loops %>% arrange(MA9_on) %>%
  mutate(gene = factor(gene, levels = gene))

# --> 4. Long format for points ----
df_long <- df_loops %>%
  pivot_longer(cols = c(MA9_on, MA9_off), names_to = "condition", values_to = "loops") %>%
  mutate(condition = factor(condition, levels = c("MA9_on", "MA9_off"),
                            labels = c("MA9-on", "MA9-off")))

# Line color
df_loops <- df_loops %>%
  mutate(line_color = ifelse(MA9_on >= MA9_off, "grey60", "#6BAED6"))

# --> 5. Plot ----
Figure_S2E <- ggplot() +
  geom_segment(data = df_loops,
               aes(x = MA9_off, xend = MA9_on, y = gene, yend = gene, color = line_color),
               linewidth = 0.7, alpha = 0.7) +
  scale_color_identity() +
  geom_point(data = df_long %>% filter(condition == "MA9-on"),
             aes(x = loops, y = gene), color = "#D73027", size = 2.5) +
  geom_point(data = df_long %>% filter(condition == "MA9-off"),
             aes(x = loops, y = gene), color = "#4575B4", size = 2.5) +
  labs(x = "Number of E-P loops", y = NULL,
       title = "3D E-P loops per core MA9 target gene",
       subtitle = paste0("(", nrow(df_loops), " of ", length(core_genes), " core genes with HiC loops)")) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 8, face = "italic"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40")
  )

print(Figure_S2E)
# Manual legend
Figure_S2E <- Figure_S2E +
  annotate("point", x = max(df_loops$MA9_on) - 5, y = 3, color = "#D73027", size = 2.5) +
  annotate("text", x = max(df_loops$MA9_on) - 4, y = 3, label = "MA9-on", size = 3, hjust = 0) +
  annotate("point", x = max(df_loops$MA9_on) - 5, y = 1.5, color = "#4575B4", size = 2.5) +
  annotate("text", x = max(df_loops$MA9_on) - 4, y = 1.5, label = "MA9-off", size = 3, hjust = 0)

print(Figure_S2E)

####################### CRISRP screen ################################### ---- 
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggsignif)
library(patchwork)

# Color palette — distinct, publication-quality, colorblind-safe
# Avoid pure RGB (too saturated) and red-green combos (colorblind issue)
col_hema  <-  "#08306B" # muted red (from RColorBrewer Set1)
col_novel <-   "#4292C6"#6BAED6"#"#9ECAE1"  # muted blue
col_ctrl  <- "grey80"#"#4DAF4A"   # muted green — OK here because not paired with red in same encoding
# NOTE: ctrl is only used for background reference dots; hema vs novel is the key
#       comparison, and red vs blue is fully colorblind-safe
col_essential <-   "#D6604D"  # dark red for essential hits in proportion plot

# Common theme for publication
theme_pub <- function(base_size = 12) {
  theme_pubr(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      axis.title    = element_text(size = base_size),
      axis.text     = element_text(size = base_size - 1),
      legend.position = "right"
    )
}

# Problematic control regions to exclude before plotting
# #113 = low GC (30%) + low complexity
# #78  = high GC (49%) + duplication risk
# #50  = PAM-rich
# #70  = targeting unannotated genes
ctrl_exclude <- c("random_target_113", "random_target_78",
                  "random_target_50",  "random_target_70")

# ── 1. Data Preparation (common to both screens) ─────────────────────────────

# --- 1a. Define sgRNA uniqueness across pools ---
pool2 <- read.csv("pool2_core_sgRNA.csv",   stringsAsFactors = FALSE)
pool3 <- read.csv("pool3_active_sgRNA.csv",  stringsAsFactors = FALSE)

three_pools <- bind_rows(pool2, pool3) %>%
  filter(Target_ID == "enhancer") %>%
  distinct(sgRNA_ID, .keep_all = TRUE)

# Unique enhancer sgRNAs
effect   <- read.csv("enhancers_withEffect_unique.csv",    stringsAsFactors = FALSE)
NOeffect <- read.csv("enhancers_withoutEffect_unique.csv", stringsAsFactors = FALSE)

enhancer_unique <- bind_rows(effect, NOeffect) %>%
  mutate(unique = TRUE)
colnames(enhancer_unique)[colnames(enhancer_unique) == "enh_sgrna"] <- "sgRNA_ID"
enhancer_unique <- enhancer_unique %>%
  dplyr::select(sgRNA_ID, unique)

# Non-unique enhancer sgRNAs
enhancer_nonunique <- three_pools %>%
  dplyr::select(sgRNA_ID) %>%
  anti_join(enhancer_unique, by = "sgRNA_ID") %>%
  mutate(unique = FALSE)

all_enhancer_sgrna <- bind_rows(
  mutate(enhancer_unique, unique = TRUE),
  enhancer_nonunique
)

# Control sgRNAs
ctrl_sgrna <- read.csv("control_sgrna.csv", stringsAsFactors = FALSE) %>%
  group_by(sgRNA_ID) %>%
  mutate(unique = n() == 1) %>%
  ungroup() %>%
  distinct(sgRNA_ID, .keep_all = TRUE) %>%
  dplyr::select(sgRNA_ID, unique)

# Combined: enhancer + control sgRNAs with region label
enh_ctrl <- bind_rows(ctrl_sgrna, all_enhancer_sgrna) %>%
  mutate(region = case_when(
    grepl("REGION", sgRNA_ID) ~ "enhancer",
    grepl("random", sgRNA_ID) ~ "random_closed",
    TRUE ~ NA_character_
  ))


# --- 1b. Process screen results ---
process_screen <- function(rra_file, enhancer_file = "DF_ENHANCER.csv") {
  # Read MAGeCK RRA sgRNA summary
  sgrna_rra <- read.csv(rra_file, stringsAsFactors = FALSE)
  colnames(sgrna_rra)[colnames(sgrna_rra) == "Gene"] <- "enhancer"
  
  # Add uniqueness and region info
  sgrna_rra <- sgrna_rra %>%
    left_join(enh_ctrl, by = c("sgrna" = "sgRNA_ID"))
  
  # Compute per-enhancer stats (using unique sgRNAs only)
  sgrna_rra <- sgrna_rra %>%
    group_by(enhancer) %>%
    mutate(
      nonunique_sgrna = sum(unique == FALSE, na.rm = TRUE),
      unique_sgrna    = sum(unique == TRUE,  na.rm = TRUE),
      mean_lfc        = mean(LFC[unique == TRUE],   na.rm = TRUE),
      median_lfc      = median(LFC[unique == TRUE], na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      class_sgrna      = if_else(LFC < 0 & FDR < 0.05, "depleted", "others"),
      class_mean_lfc   = if_else(mean_lfc < 0,   "depleted", "others"),
      class_median_lfc = if_else(median_lfc < 0,  "depleted", "others")
    )
  
  # Add novel/known annotation from enhancer master table
  df_enh <- read_csv(enhancer_file, show_col_types = FALSE)
  novel_info <- df_enh %>%
    dplyr::select(name_dhs_dox, novel_atac_peak, logCPM_h3k27ac)
  colnames(novel_info)[1] <- "enhancer"
  
  sgrna_rra <- sgrna_rra %>%
    left_join(novel_info, by = "enhancer")
  
  # Distinct enhancer-level summary (only enhancers with all 10 unique sgRNAs)
  res_dist <- sgrna_rra %>%
    distinct(enhancer, .keep_all = TRUE) %>%
    filter(unique_sgrna == 10) %>%
    mutate(
      novel_atac_peak = replace_na(as.character(novel_atac_peak), "ctrl"),
      novel_atac_peak = factor(novel_atac_peak, levels = c("ctrl", "known", "novel"))
    )
  
  # Assign groups
  res_dist <- res_dist %>%
    mutate(group = case_when(
      region == "random_closed"    ~ "ctrl",
      novel_atac_peak == "known"   ~ "hema",
      novel_atac_peak == "novel"   ~ "novel",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(group)) %>%
    mutate(group = factor(group, levels = c("ctrl", "hema", "novel")))
  
  # Rank by median LFC
  res_dist <- res_dist %>%
    # Remove problematic control regions
    filter(!enhancer %in% ctrl_exclude) %>%
    arrange(median_lfc) %>%
    mutate(rank = row_number())
  
  list(sgrna = sgrna_rra, dist = res_dist)
}

# Process screen #1 (pool3) and screen #2 (pool2)
screen1 <- process_screen("mll_af9_crispr_pool3_rra.sgrna_summary.csv")
screen2 <- process_screen("mll_af9_crispr_pool2_rra.sgrna_summary.csv")  # Update filename

# Figure_3B ----

Figure_3B <- ggplot(screen1$dist,
                  aes(x = rank, y = median_lfc, fill = group)) +
  geom_col(width = 1, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  scale_fill_manual(
    values = c("ctrl" = "grey75", "hema" = col_hema, "novel" = col_novel),
    labels = c("ctrl", "hema", "novel"),
    name   = "pEs:"
  ) +
  labs(
    title = "CRISPR screen #1",
    x     = "Ranked by CRISPR median log2FC",
    y     = "CRISPR median log2FC\n(d12 vs d0)"
  ) +
  theme_pub() +
  theme(
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    legend.key.size = unit(0.4, "cm")
  ) +
  coord_cartesian(ylim = c(-1.2, 0.8))

print (Figure_3B)


# Figure_3C ----

Figure_3C <- ggplot(screen1$dist,
                  aes(x = group, y = median_lfc, fill = group)) +
  geom_violin(trim = FALSE, scale = "width", color = "black", linewidth = 0.3, alpha = 0.9) +
  geom_boxplot(width = 0.12, outlier.shape = 18, outlier.size = 1.5, fill = "white") +
  stat_compare_means(
    comparisons = list(c("ctrl", "hema"), c("ctrl", "novel")),
    method      = "wilcox.test",
    label       = "p.format",
    bracket.size = 0.4
  ) +
  scale_fill_manual(values = c("ctrl" = col_ctrl, "hema" = col_hema, "novel" = col_novel)) +
  labs(
    title = "CRISPR screen #1",
    x     = NULL,
    y     = "Median log2FC\n(d12 vs d0)"
  ) +
  theme_pub() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(-1.5, 1.2))

print (Figure_3C)
# Figure_3D ----

# --- Fisher's exact tests (pairwise vs ctrl) ---
fisher_prop <- function(screen_dist) {
  contingency <- screen_dist %>%
    count(group, class_median_lfc) %>%
    pivot_wider(names_from = class_median_lfc, values_from = n, values_fill = 0)
  
  groups <- levels(screen_dist$group)
  pairs  <- combn(groups, 2, simplify = FALSE)
  
  map_dfr(pairs, function(pair) {
    sub <- contingency %>% filter(group %in% pair)
    mat <- as.matrix(sub[, -1])
    rownames(mat) <- sub$group
    ft  <- fisher.test(mat)
    tibble(
      group1     = pair[1],
      group2     = pair[2],
      p_value    = ft$p.value,
      odds_ratio = ft$estimate,
      label      = ifelse(ft$p.value < 0.001, formatC(ft$p.value, format = "e", digits = 1),
                          ifelse(ft$p.value < 0.05, paste0("p = ", round(ft$p.value, 3)),
                                 "n.s."))
    )
  })
}

fisher_s1 <- fisher_prop(screen1$dist)
cat("\n=== Fisher's Exact Test — Screen #1 Proportion ===\n")
print(fisher_s1)

# --- Proportion data ---
prop_data_s1 <- screen1$dist %>%
  group_by(group, class_median_lfc) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(
    total   = sum(count),
    percent = round(count / total * 100),
    label   = paste0(count, "\n(", percent, "%)")
  ) %>%
  ungroup()

# --- Significance bracket positions ---
# Build annotation df for geom_signif on proportion (0-1) scale
sig_annot_s1 <- fisher_s1 %>%
  mutate(
    xmin = as.numeric(factor(group1, levels = c("ctrl","hema","novel"))),
    xmax = as.numeric(factor(group2, levels = c("ctrl","hema","novel"))),
    y_position = seq(1.06, by = 0.09, length.out = n())
  )

Figure_3D <- ggplot(prop_data_s1,
                  aes(x = group, y = count, fill = class_median_lfc)) +
  geom_col(position = "fill", color = "black", width = 0.7) +
  geom_text(
    aes(label = label),
    position = position_fill(vjust = 0.5),
    size = 3.5, color = "black"
  ) +
  scale_fill_manual(
    values = c("depleted" = col_essential, "others" = "white"),
    labels = c("Median log2FC<0\n(essential hit)", "Others"),
    name   = NULL
  ) +
  scale_x_discrete(labels = c("ctrl" = "ctrl", "hema" = "hema", "novel" = "novel")) +
  # Manual significance brackets
  annotate("segment",
           x = sig_annot_s1$xmin, xend = sig_annot_s1$xmax,
           y = sig_annot_s1$y_position, yend = sig_annot_s1$y_position,
           linewidth = 0.4) +
  annotate("segment",
           x = sig_annot_s1$xmin, xend = sig_annot_s1$xmin,
           y = sig_annot_s1$y_position - 0.015, yend = sig_annot_s1$y_position,
           linewidth = 0.4) +
  annotate("segment",
           x = sig_annot_s1$xmax, xend = sig_annot_s1$xmax,
           y = sig_annot_s1$y_position - 0.015, yend = sig_annot_s1$y_position,
           linewidth = 0.4) +
  annotate("text",
           x = (sig_annot_s1$xmin + sig_annot_s1$xmax) / 2,
           y = sig_annot_s1$y_position + 0.02,
           label = sig_annot_s1$label,
           size = 3.2) +
  labs(
    title = "CRISPR screen #1",
    x     = NULL,
    y     = "Proportion"
  ) +
  theme_pub() +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(0, max(sig_annot_s1$y_position) + 0.06))

print (Figure_3D)

# Figure_S3B ----

# Load enhancer-gene interaction data
df_enh <- DF_ENHANCER
df_pro <- df_promoter
colnames(df_pro)[4] <- "hic_gene_dox"

# Prepare enhancer-gene mapping
enh_gene <- df_enh %>%
  dplyr::select(name_dhs_dox, hic_gene_dox, interaction_class_dox, target_gene_dox) %>%
  mutate(hic_gene_dox = if_else(
    interaction_class_dox == "nonHiC", target_gene_dox, hic_gene_dox
  )) %>%
  mutate(hic_gene_dox = strsplit(as.character(hic_gene_dox), ",\\s*")) %>%
  unnest(cols = hic_gene_dox)
colnames(enh_gene)[colnames(enh_gene) == "name_dhs_dox"] <- "enhancer"

# Add gene expression class
expr_info <- df_pro %>%
  dplyr::select(hic_gene_dox, class_expr)
colnames(expr_info)[colnames(expr_info) == "class_expr"] <- "class_expr_hic"

enh_gene <- enh_gene %>%
  left_join(expr_info, by = "hic_gene_dox")

# Get top 30 dropout enhancers from screen #1.
# Strategy for multi-gene pEs:
#   - Dots: one per pE (use distinct on enhancer for geom_point)
#   - Labels: one per enhancer-gene pair (keep all rows for geom_text_repel)
#     so each gene gets its own correctly-colored label next to the dot.

# Join enh_gene first to get interaction_class_dox and gene labels
top30_all <- screen1$dist %>%
  filter(region == "enhancer", median_lfc < 0) %>%
  arrange(median_lfc) %>%
  slice_head(n = 30) %>%
  left_join(enh_gene, by = "enhancer") %>%
  mutate(
    enhancer_label      = gsub("REGION_", "RG_", enhancer),
    interaction_grouped = if_else(
      interaction_class_dox %in% c("HiC", "HiC_neighbor"), "E-P loop", "Linear closest gene"
    )
  )

# top30_all already has one row per enhancer-gene pair from the unnest.
# Use it directly for both geom_point and geom_text — geom_point will naturally
# overdraw duplicate dots for multi-gene pEs at the same position.
top30 <- top30_all %>%
  distinct(enhancer, hic_gene_dox, .keep_all = TRUE) %>%
  mutate(enhancer_label = factor(enhancer_label,
                                 levels = unique(enhancer_label[order(-rank)])))

Figure_S3B <- ggplot(top30,
                  aes(x = enhancer_label,
                      y = median_lfc,
                      fill = group,
                      shape = interaction_grouped)) +
  geom_point(size = 4, color = "black") +
  scale_shape_manual(
    values = c("E-P loop" = 21, "Linear closest gene" = 23),
    name   = "Target gene defined by:"
  ) +
  # geom_text: one label per row, so multi-gene pEs get multiple labels at same position
  geom_text(
    aes(label = hic_gene_dox, color = class_expr_hic),
    size = 3, vjust = 0.5, hjust = -0.3,fontface="italic",
  ) +
  scale_fill_manual(
    values = c("hema" = col_hema, "novel" = col_novel),
    name   = "pEs:"
  ) +
  scale_color_manual(
    values = c("down" = "#377EB8", "neutral" = "#999999", "up" = "#E41A1C"),
    name   = "Target gene expression:"
  ) +
  labs(
    title = "CRISPR screen #1\ntop 30 dropout pEs",
    x     = NULL,
    y     = "CRISPR median log2FC"
  ) +
  theme_pub(base_size = 10) +
  theme(
    axis.text.y  = element_text(size = 8),
    legend.position = "top",
    legend.box = "vertical"
  ) +
  coord_flip(ylim = c(-3, 3))

print (Figure_S3B)


# Figure_3I ----
Figure_3I <- ggplot(screen2$dist,
                  aes(x = rank, y = median_lfc, fill = group)) +
  geom_col(width = 1, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  scale_fill_manual(
    values = c("ctrl" = "grey75", "hema" = col_hema, "novel" = col_novel),
    name   = "pEs:"
  ) +
  labs(
    title = "CRISPR screen #2",
    x     = "Ranked by CRISPR median log2FC",
    y     = "CRISPR median log2FC\n(d12 vs d0)"
  ) +
  theme_pub() +
  theme(
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    legend.key.size = unit(0.4, "cm")
  ) +
  coord_cartesian(ylim = c(-1.5, 1))

print (Figure_3I)

# Figure_3J ----
Figure_3J <- ggplot(screen2$dist,
                  aes(x = group, y = median_lfc, fill = group)) +
  geom_violin(trim = FALSE, scale = "width", color = "black", linewidth = 0.3, alpha = 0.9) +
  geom_boxplot(width = 0.12, outlier.shape = 18, outlier.size = 1.5, fill = "white") +
  stat_compare_means(
    comparisons = list(c("ctrl", "hema"), c("ctrl", "novel")),
    method      = "wilcox.test",
    label       = "p.format",
    bracket.size = 0.4
  ) +
  scale_fill_manual(values = c("ctrl" = col_ctrl, "hema" = col_hema, "novel" = col_novel)) +
  labs(
    title = "CRISPR screen #2",
    x     = NULL,
    y     = "Median log2FC\n(d12 vs d0)"
  ) +
  theme_pub() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(-1.5, 1.2))

print(Figure_3J)
# Figure_3K ----
fisher_s2 <- fisher_prop(screen2$dist)
cat("\n=== Fisher's Exact Test — Screen #2 Proportion ===\n")
print(fisher_s2)

prop_data_s2 <- screen2$dist %>%
  group_by(group, class_median_lfc) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(
    total   = sum(count),
    percent = round(count / total * 100),
    label   = paste0(count, "\n(", percent, "%)")
  ) %>%
  ungroup()

sig_annot_s2 <- fisher_s2 %>%
  mutate(
    xmin = as.numeric(factor(group1, levels = c("ctrl","hema","novel"))),
    xmax = as.numeric(factor(group2, levels = c("ctrl","hema","novel"))),
    y_position = seq(1.06, by = 0.09, length.out = n())
  )

Figure_3K <- ggplot(prop_data_s2,
                  aes(x = group, y = count, fill = class_median_lfc)) +
  geom_col(position = "fill", color = "black", width = 0.7) +
  geom_text(
    aes(label = label),
    position = position_fill(vjust = 0.5),
    size = 3.5, color = "black"
  ) +
  scale_fill_manual(
    values = c("depleted" = col_essential, "others" = "white"),
    labels = c("Median log2FC<0\n(essential hit)", "Others"),
    name   = NULL
  ) +
  scale_x_discrete(labels = c("ctrl" = "ctrl", "hema" = "hema", "novel" = "novel")) +
  # Manual significance brackets
  annotate("segment",
           x = sig_annot_s2$xmin, xend = sig_annot_s2$xmax,
           y = sig_annot_s2$y_position, yend = sig_annot_s2$y_position,
           linewidth = 0.4) +
  annotate("segment",
           x = sig_annot_s2$xmin, xend = sig_annot_s2$xmin,
           y = sig_annot_s2$y_position - 0.015, yend = sig_annot_s2$y_position,
           linewidth = 0.4) +
  annotate("segment",
           x = sig_annot_s2$xmax, xend = sig_annot_s2$xmax,
           y = sig_annot_s2$y_position - 0.015, yend = sig_annot_s2$y_position,
           linewidth = 0.4) +
  annotate("text",
           x = (sig_annot_s2$xmin + sig_annot_s2$xmax) / 2,
           y = sig_annot_s2$y_position + 0.02,
           label = sig_annot_s2$label,
           size = 3.2) +
  labs(
    title = "CRISPR screen #2",
    x     = NULL,
    y     = "Proportion"
  ) +
  theme_pub() +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(0, max(sig_annot_s2$y_position) + 0.06))

print(Figure_3K)

# Figure_S3C-D ----
# Step 1: For nonHiC enhancers, fall back to closest TSS gene
# Step 2: Split comma-separated hic_gene_dox → one row per enhancer-gene pair
df_enh_mod_split <- DF_ENHANCER %>%
  dplyr::select(name_dhs_dox, hic_gene_dox, target_gene_dox, interaction_class_dox) %>%
  mutate(
    hic_gene_dox = if_else(interaction_class_dox == "nonHiC", target_gene_dox, hic_gene_dox),
    hic_gene_dox = strsplit(as.character(hic_gene_dox), ",\\s*")
  ) %>%
  unnest(cols = hic_gene_dox) %>%
  dplyr::rename(enhancer = name_dhs_dox)

# Step 3: Join DF_PROMOTER to get interest_group (e.g. "MLL_target_interest")
df_pro <- df_promoter
colnames(df_pro)[4] <- "hic_gene_dox"   # 4th column is the gene name

df_enh_mod_split <- df_enh_mod_split %>%
  left_join(
    df_pro %>% dplyr::select(hic_gene_dox, interest_group),
    by = "hic_gene_dox"
  )

# Step 4: Join onto screen2 results
res_spec_dist_02_s2 <- screen2$dist %>%
  ungroup() %>%
  arrange(median_lfc) %>%
  mutate(rank_median = row_number()) %>%
  left_join(df_enh_mod_split, by = "enhancer")

# --- 8a. All hema (known) pEs that passed QC in screen1, ranked by H3K27ac ---
# screen1$dist already has novel_atac_peak and logCPM_h3k27ac (joined inside
# process_screen from DF_ENHANCER). We just filter by group and rank.
hema_all <- screen1$dist %>%
  filter(group == "hema") %>%
  arrange(logCPM_h3k27ac) %>%
  mutate(rank_h3k27ac = row_number())

# --- 8b. All novel pEs that passed QC in screen1, ranked by H3K27ac ---
novel_all <- screen1$dist %>%
  filter(group == "novel") %>%
  arrange(logCPM_h3k27ac) %>%
  mutate(rank_h3k27ac = row_number())

# --- 8c. Count chosen = top 100 by H3K27ac within each group ---
# "Chosen" means the top 100 H3K27ac candidates that also passed unique_sgrna==10.
# Since screen1$dist IS the unique_sgrna==10 filtered set, we just count
# how many of those fall in the top 100 ranked by H3K27ac.
n_hema_chosen  <- min(100, nrow(hema_all))   # up to 100
n_novel_chosen <- min(100, nrow(novel_all))

# The dashed threshold line = H3K27ac value at the top-100 cutoff
thresh_hema  <- hema_all  %>% arrange(desc(logCPM_h3k27ac)) %>% slice(n_hema_chosen)  %>% pull(logCPM_h3k27ac)
thresh_novel <- novel_all %>% arrange(desc(logCPM_h3k27ac)) %>% slice(n_novel_chosen) %>% pull(logCPM_h3k27ac)
h3k27ac_threshold <- min(thresh_hema, thresh_novel)

# --- 9a. Prepare data: dropout pEs only (median_lfc < 0), ordered by count per gene ---
# Genes are ordered by count of dropout pE-gene rows (same as your original:
# sort(table(hic_gene_dox), increasing = TRUE)), which places highest-count gene first.
suppD_data <- res_spec_dist_02_s2 %>%
  filter(region == "enhancer",
         interest_group == "MLL_target_interest") %>%
  mutate(
    novel_atac_peak = factor(novel_atac_peak, levels = c("known", "novel")),
    hic_gene_dox    = factor(hic_gene_dox,
                             levels = names(sort(table(hic_gene_dox), decreasing = TRUE)))
  )

# Count unique pEs for legend (multi-gene pEs appear >1 row but are 1 pE)
n_hema_pE_s2  <- suppD_data %>% filter(novel_atac_peak == "known") %>% distinct(enhancer) %>% nrow()
n_novel_pE_s2 <- suppD_data %>% filter(novel_atac_peak == "novel") %>% distinct(enhancer) %>% nrow()
n_genes_s2    <- n_distinct(suppD_data$hic_gene_dox)

Figure_S3C <- ggplot(suppD_data,
                 aes(x = hic_gene_dox, fill = novel_atac_peak)) +
  geom_bar(color = "black", width = 0.75) +
  geom_text(
    aes(label = after_stat(count)),
    stat     = "count",
    position = position_stack(vjust = 0.5),
    size = 8, color = "white", #fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("known" = col_hema, "novel" = col_novel),
    labels = c(paste0("hema (", n_hema_pE_s2,  " pEs)"),
               paste0("novel (", n_novel_pE_s2, " pEs)")),
    name = NULL
  ) +
  labs(
    title    = "CRISPR Screen #2",
    subtitle = paste0(n_hema_pE_s2 + n_novel_pE_s2, " pEs in total (", n_genes_s2, " target genes)"),
    x        = NULL,
    y        = "Count of pEs"
  ) +
  theme_pub() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "top"
  )

print(Figure_S3C)

# --- 10a. Prepare data: dropout pEs only (median_lfc < 0) ---
red_genes     <- c("Meis1")        # genes to label red on y-axis
highlight_ids <- c("RG_16304")     # specific pE labels to show in green

suppE_data <- res_spec_dist_02_s2 %>%
  filter(region == "enhancer",
         median_lfc < 0) %>%
  # Deduplicate: one row per enhancer-gene pair (removes duplicates from multi-gene unnest)
  distinct(enhancer, hic_gene_dox, .keep_all = TRUE) %>%
  mutate(
    enhancer_label  = gsub("REGION_", "RG_", enhancer),
    novel_atac_peak = factor(novel_atac_peak,
                             levels = c("known", "novel"),
                             labels = c("hema", "novel")),
    label_color     = ifelse(enhancer_label %in% highlight_ids, "green4", "black")
  ) %>%
  # Order genes by count of dropout pEs, most pEs on top (same as your original)
  mutate(hic_gene_dox = factor(hic_gene_dox,
                               levels = names(sort(table(hic_gene_dox), decreasing = FALSE))))


# Auto-calculate subtitle numbers
n_dropout_pE   <- n_distinct(suppE_data$enhancer)
n_dropout_gene <- n_distinct(suppE_data$hic_gene_dox)

# --- 10b. Build color vectors for italic gene axis labels ---
gene_order  <- levels(suppE_data$hic_gene_dox)
axis_colors <- ifelse(gene_order %in% red_genes, "#CC0000", "black")

Figure_S3D <- ggplot(suppE_data,
                 aes(x    = rank_median,   # rank_median: arranged by median_lfc ascending,
                     # so most depleted (most negative) = leftmost
                     y    = hic_gene_dox,
                     fill = novel_atac_peak)) +
  geom_point(size = 5, alpha = 1, shape = 21, color = "black") +
  geom_text_repel(
    aes(label = enhancer_label,
        color = label_color),
    size         = 2.5,
    direction    = "x",
    nudge_y      = 0.35,
    segment.size = 0.2,
    max.overlaps = Inf,
    show.legend  = FALSE
  ) +
  scale_fill_manual(
    values = c("hema" = col_hema, "novel" = col_novel),
    name   = NULL
  ) +
  scale_color_identity() +
  labs(
    title    = bquote(italic("CRISPR screen #2")),
    subtitle = paste0("Dropout pEs (", n_dropout_pE, " pEs targeting ", n_dropout_gene, " genes in total)"),
    x        = "pEs ranked by CRISPR median log2FC",
    y        = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.y     = element_text(face  = "italic",
                                   color = axis_colors,
                                   size  = 11),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(override.aes = list(size = 6)))

print(Figure_S3D)

############################# Figure 4B ################################# ----
library(ggplot2)
library(dplyr)

df_pie_hoxa9 <- DF_ATAC %>%
  filter(peak_hoxa9 == "yes") %>%
  count(annot.type) %>%
  mutate(percentage = n / sum(n) * 100,
         label = paste0(annot.type, "\n", n, " (", round(percentage, 1), "%)"))

df_pie_meis1 <- DF_ATAC %>%
  filter(peak_meis1 == "yes") %>%
  count(annot.type) %>%
  mutate(percentage = n / sum(n) * 100,
         label = paste0(annot.type, "\n", n, " (", round(percentage, 1), "%)"))

pie_hoxa9<-ggplot(df_pie_hoxa9, aes(x = "", y = n, fill = annot.type)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            size = 3) +
  labs(title = "peak_hoxa9 - Pie Chart", x = NULL, y = NULL, fill = "Annotation Type") +
  theme_pubr() +
  theme(text = element_text(size = 10))+
  scale_fill_manual(values = c(  "#74ADD1", "#4575B4", "#D6604D")) 

pie_meis1<-ggplot(df_pie_meis1, aes(x = "", y = n, fill = annot.type)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            size = 3) +
  labs(title = "peak_meis1 - Pie Chart", x = NULL, y = NULL, fill = "Annotation Type") +
  theme_pubr() +
  theme(text = element_text(size = 10))+
  scale_fill_manual(values = c(  "#74ADD1", "#4575B4", "#D6604D")) 

print (pie_hoxa9)
print (pie_meis1)
############################# Figure 4G ################################# ----
# --- pie chart for arid1b bound regions ----
library(ggplot2)
library(dplyr)

df_pie_arid1b <- DF_ATAC %>%
  filter(peak_arid1b == "yes") %>%
  count(annot.type) %>%
  mutate(percentage = n / sum(n) * 100,
         label = paste0(annot.type, "\n", n, " (", round(percentage, 1), "%)"))


pie_arid1b<-ggplot(df_pie_arid1b, aes(x = "", y = n, fill = annot.type)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            size = 3) +
  labs(title = "peak_arid1b - Pie Chart", x = NULL, y = NULL, fill = "Annotation Type") +
  theme_pubr() +
  theme(text = element_text(size = 10))+
  scale_fill_manual(values = c(  "#74ADD1", "#4575B4", "#D6604D")) 

print(pie_arid1b)
############################# Figure 4I ################################# ----
library(eulerr)

DF_ENH<- DF_ENHANCER %>%
  mutate(class_atac = case_when(
    logFC_atac < -1 & class_atac == "down" ~ "down",
    logFC_atac > 1 & class_atac == "up" ~ "up",
    TRUE ~ "neutral"
  ))
# filter gene sets

df_arid1b<-filter(DF_ENHANCER, class_atac=="up",peak_arid1b_dox=="yes")
df_meis1<-filter(DF_ENHANCER, class_atac=="up",peak_meis1_dox=="yes")
df_hoxa9<-filter(DF_ENHANCER, class_atac=="up",peak_hoxa9=="yes")
# or
df_arid1b<-filter(DF_ENHANCER, peak_arid1b_dox=="yes")
df_meis1<-filter(DF_ENHANCER, peak_meis1_dox=="yes")
df_hoxa9<-filter(DF_ENHANCER, peak_hoxa9=="yes")
#or DF_ENH
df_arid1b<-filter(DF_ENH, class_atac=="up",peak_arid1b_dox=="yes")
df_meis1<-filter(DF_ENH, class_atac=="up",peak_meis1_dox=="yes")
df_hoxa9<-filter(DF_ENH, class_atac=="up",peak_hoxa9=="yes")
#
ARID1B <- df_arid1b$name_dhs_dox
MEIS1 <- df_meis1$name_dhs_dox
HOXA9 <- df_hoxa9$name_dhs_dox


# Create a named list
set1_3I <- list( HOXA9=HOXA9,MEIS1 = MEIS1,ARID1B = ARID1B)
# Make Venn input
venn_data <- euler(set1_3I)

# Plot
venn_arid1b_hoxa9_meis1<-
  plot(venn_data,
     quantities = TRUE,
     fills = list(fill = c("#377EB8","#E0F3F8", "#A6CEE3")),
     legend = TRUE)

print(venn_arid1b_hoxa9_meis1)

################################## Figure 5D ############################ ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# Reshape data: long format for TF binding status
df_long <- DF_ATAC %>%
  pivot_longer(
    cols = c(peak_arid1b, peak_arid1a, peak_brg1),
    names_to = "TF",
    values_to = "bound"
  ) %>%
  filter(bound == "yes") %>%
  mutate(TF = recode(TF,
                     "peak_arid1b" = "ARID1B",
                     "peak_arid1a" = "ARID1A",
                     "peak_brg1" = "BRG1"))

# Prepare labels
df_labels <- df_long %>%
  count(TF, annot.type) %>%
  group_by(TF) %>%
  mutate(percentage = n / sum(n) * 100,
         label = paste0(n, "\n(", round(percentage, 1), "%)"))

# Set the desired order manually
df_long$TF <- factor(df_long$TF, levels = c("ARID1B", "ARID1A", "BRG1"))
df_labels$TF <- factor(df_labels$TF, levels = c("ARID1B", "ARID1A", "BRG1"))


# Plot
proportion_bound_regions<-
  ggplot(df_long, aes(x = TF, fill = annot.type)) +
  geom_bar(position = "fill", alpha = 0.9) +
  geom_text(data = df_labels,
            aes(x = TF, y = n / sum(n), label = label),
            position = position_fill(vjust = 0.5),
            size = 4, color = "black") +
  labs(x = NULL, y = "Proportion", title = "TF-bound Peaks by Annotation Type") +
  theme_pubr() +
  scale_fill_brewer(palette = "Set2") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_manual(values = c(  "#74ADD1", "#4575B4", "#D6604D")) 

print (proportion_bound_regions)
################################## Figure 5E ############################ ----
library(eulerr)

# filter gene sets
df_arid1b<-filter(DF_ENHANCER, class_atac=="up",logFC_atac>1, peak_arid1b_dox=="yes")
df_arid1a<-filter(DF_ENHANCER, class_atac=="up",logFC_atac>1, peak_arid1a_dox=="yes")
df_brg1<-filter(DF_ENHANCER, class_atac=="up",logFC_atac>1, peak_brg1_dox=="yes")


ARID1B <- df_arid1b$name_dhs_dox
ARID1A <- df_arid1a$name_dhs_dox
BRG1 <- df_brg1$name_dhs_dox


# Create a named list
set_4E <- list(ARID1A = ARID1A,ARID1B = ARID1B, BRG1 = BRG1 )
# Make Venn input
venn_data <- euler(set_4E)

# Plot
venn_arid1a_arid1b_brg1<-
plot(venn_data,
     quantities = TRUE,
     fills = list(fill = c("#E0F3F8", "#A6CEE3","#377EB8")),
     legend = TRUE)
print(venn_arid1a_arid1b_brg1)
################################## Figure 5F ############################ ----
# GSEA: MLL-AF9 up/down expression signature vs ARID1B-KD ranking.
library(fgsea)
library(patchwork)

## 1) MLL-AF9 signature (from MA9 on/off DE)
mllaf9_up   <- unique(na.omit(with(df_promoter, name[logFC_expr >  2 & class_expr == "up"])))
mllaf9_down <- unique(na.omit(with(df_promoter, name[logFC_expr < -2 & class_expr == "down"])))
pathways    <- list(MLLAF9_Up = mllaf9_up, MLLAF9_Down = mllaf9_down)

## 2) Preranked vector: ARID1B-KD log2FC (clean, de-duplicated, sorted high -> low)
ranks <- setNames(df_promoter$logFC_expr_arid1bKD, df_promoter$name)
ranks <- ranks[is.finite(ranks) & !is.na(names(ranks))]
if (any(duplicated(names(ranks)))) ranks <- tapply(ranks, names(ranks), median)
ranks <- setNames(as.numeric(ranks), names(ranks))
ranks <- ranks[ranks != 0]
set.seed(42)
ranks <- ranks + rnorm(length(ranks), sd = 1e-7)
ranks <- sort(ranks, decreasing = TRUE)

## 3) fgsea (multilevel, eps = 0 -> accurate small p-values; matches Fig 6D)
set.seed(42)
res <- fgseaMultilevel(pathways = pathways, stats = ranks,
                       minSize = 15, maxSize = 5000, eps = 0)

## 4) Shared enrichment-panel builder (also used by Fig 6D)
##    3 tracks: running ES + gene-hit ticks + ranked-metric gradient strip.
make_gsea_panel <- function(genes, ranks, res_row, title, subtitle,
                            high_label = "ARID1B-high",   # 6D default
                            low_label  = "ARID1B-low",
                            high_col   = "#d7301f",        # left label colour
                            low_col    = "black",
                            n_bins = 512) {
  pd   <- fgsea::plotEnrichmentData(pathway = genes, stats = ranks)
  NES  <- res_row$NES
  padj <- res_row$padj
  N    <- max(pd$stats$rank)
  curve_col <- if (NES > 0) "#d7301f" else "#3182bd"   # red / blue by NES sign
  if (NES > 0) { ann_x <- 0.50 * N; ann_y <- 0.62 * pd$posES }
  else         { ann_x <- 0.03 * N; ann_y <- 0.55 * pd$negES }
  ann_txt <- sprintf("NES = %.2f\npadj = %s",
                     NES, formatC(padj, format = "e", digits = 1))
  
  p_es <- ggplot(pd$curve, aes(rank, ES)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_line(color = curve_col, linewidth = 1) +
    annotate("text", x = ann_x, y = ann_y, label = ann_txt,
             hjust = 0, vjust = 1, size = 4, lineheight = 0.95) +
    labs(y = "Enrichment score", title = title, subtitle = subtitle) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, N)) +
    theme_classic(base_size = 12) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), axis.line.x = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          plot.margin = margin(4, 8, 0, 4))
  
  p_ticks <- ggplot(pd$ticks) +
    geom_segment(aes(x = rank, xend = rank, y = 0, yend = 1), linewidth = 0.2) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, N)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(0, 8, 0, 4))
  
  strip <- pd$stats %>%
    mutate(bin = cut(rank, breaks = n_bins, labels = FALSE)) %>%
    group_by(bin) %>%
    summarise(rank = mean(rank), stat = mean(stat), .groups = "drop")
  p_bar <- ggplot(strip, aes(x = rank, y = 0, fill = stat)) +
    geom_tile() +
    scale_fill_gradient2(low = "#3182bd", mid = "white", high = "#d7301f",
                         midpoint = 0, guide = "none") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, N)) +
    scale_y_continuous(limits = c(-1.6, 0.5), expand = c(0, 0)) +
    annotate("text", x = 0.01 * N, y = -1.0, label = high_label,
             hjust = 0, fontface = "italic", colour = high_col, size = 3.6) +
    annotate("text", x = 0.99 * N, y = -1.0, label = low_label,
             hjust = 1, fontface = "italic", colour = low_col, size = 3.6) +
    coord_cartesian(clip = "off") +
    theme_void() + theme(plot.margin = margin(0, 8, 4, 4))
  
  p_es / p_ticks / p_bar + plot_layout(heights = c(1, 0.14, 0.16))
}

## 5) Build the two panels and assemble (same look as Fig 6D)
fig_5F_activated <- make_gsea_panel(
  genes      = intersect(pathways$MLLAF9_Up,   names(ranks)),
  ranks      = ranks,
  res_row    = res[res$pathway == "MLLAF9_Up", ],
  title      = "KMT2A-MLLT3-activated genes",
  subtitle   = sprintf("(log2FC > 2, FDR < 0.05, %d genes)", length(mllaf9_up)),
  high_label = "Arid1b-KD", low_label = "Arid1b-WT", high_col = "#3182bd"
)
fig_5F_suppressed <- make_gsea_panel(
  genes      = intersect(pathways$MLLAF9_Down, names(ranks)),
  ranks      = ranks,
  res_row    = res[res$pathway == "MLLAF9_Down", ],
  title      = "KMT2A-MLLT3-suppressed genes",
  subtitle   = sprintf("(log2FC < -2, FDR < 0.05, %d genes)", length(mllaf9_down)),
  high_label = "Arid1b-KD", low_label = "Arid1b-WT", high_col = "#3182bd"
)
# side-by-side, matching the published Figure 5F layout
Figure_5F <- wrap_elements(fig_5F_activated) | wrap_elements(fig_5F_suppressed)
print(Figure_5F)

############################## Figure S4A ############################### ----
# Figure S4A - ARID1B / ARID1A per-class enhancer heatmap (MA9-on / MA9-off)
# Output: ./Figure_S4A/arid_percolumn_heatmap.pdf
#
# NOTE: requires DF_ENHANCER and DF_ATAC to be loaded in the session beforehand.

library(tidyverse)
library(ggpubr)
library(patchwork)

# ---- Data preparation

df_arid <- select(DF_ENHANCER,
                  name_dhs_dox,
                  class_activity,
                  class_atac,
                  novel_atac_peak,
                  logFC_arid1b,
                  logFC_arid1a,
                  logFC_atac_arid1bKD,
                  logFC_arid1a_arid1bKD,
                  logCPM_arid1a,
                  logCPM_arid1b,
                  logFC_atac,
                  peak_arid1b_dox,
                  peak_arid1b_veh,
                  peak_arid1a_dox,
                  peak_arid1a_veh)

DF_ATAC_ENHANCER <- filter(DF_ATAC, annot.type != "promoters")
df_ATAC_arid <- select(DF_ATAC_ENHANCER,
                       name,
                       mean_arid1b_dox,
                       mean_arid1b_veh,
                       mean_arid1a_dox,
                       mean_arid1a_veh)
colnames(df_ATAC_arid)[1] <- "name_dhs_dox"

df_arid <- df_arid %>%
  left_join(df_ATAC_arid, by = "name_dhs_dox")

dir.create("./Figure_S4A", recursive = TRUE, showWarnings = FALSE)


# ---- Labels & helpers
# Levels are taken straight from the column (df_arid, self-contained) so they
# match byte-for-byte; class_activity carries a non-ASCII character that a
# hand-typed list would silently fail to match.
class_levels <- sort(unique(as.character(df_arid$class_activity)))
class_labels <- c("a (2)","b (1)","c (1)","d (1)","e (0)",
                  "f (0)","g (0)","h (0)","i (-0)",
                  "j (-0)","k (-0)","l (-1)","m (-1)","n (-1)","o (-2)")

stopifnot(length(class_levels) == length(class_labels))   # must be 15, ordered a -> o

sig_label <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                                   p < 0.05  ~ "*",   TRUE     ~ "")

# ---- 1. Summarise CPM (fold-over-mean per class, Method A)
summary_df <- df_arid %>%
  group_by(class_activity, novel_atac_peak) %>%
  summarise(
    ARID1B_on   = mean(mean_arid1b_dox, na.rm = TRUE),
    ARID1A_on   = mean(mean_arid1a_dox, na.rm = TRUE),
    ARID1B_off  = mean(mean_arid1b_veh, na.rm = TRUE),
    ARID1A_off  = mean(mean_arid1a_veh, na.rm = TRUE),
    log2fc_1b   = mean(logFC_arid1b,    na.rm = TRUE),
    log2fc_1a   = mean(logFC_arid1a,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ARID1B_on   = ARID1B_on  / mean(ARID1B_on),
    ARID1A_on   = ARID1A_on  / mean(ARID1A_on),
    ARID1B_off  = ARID1B_off / mean(ARID1B_off),
    ARID1A_off  = ARID1A_off / mean(ARID1A_off),
    diff_on     = ARID1B_on  - ARID1A_on,
    diff_off    = ARID1B_off - ARID1A_off,
    log2fc_diff = log2fc_1b  - log2fc_1a,
    pE_class    = factor(class_activity, levels = class_levels, labels = class_labels),
    origin      = factor(
      ifelse(novel_atac_peak == "known", "hema pEs", "novel pEs"),
      levels = c("hema pEs", "novel pEs")
    )
  )

stopifnot(!any(is.na(summary_df$pE_class)))   # every class must map to a label

shared_max <- max(c(summary_df$ARID1B_on, summary_df$ARID1A_on,
                    summary_df$ARID1B_off, summary_df$ARID1A_off))

# ---- 2. Fisher's exact z-score difference (drives the significance stars)
# Per protein: observed vs expected peak count under global rate -> z-score.
# zdiff = z(ARID1B) - z(ARID1A); positive = ARID1B more enriched.
compute_enrichment_zdiff_all <- function(data) {
  N <- nrow(data)
  global_1b_on  <- sum(data$peak_arid1b_dox == "yes") / N
  global_1a_on  <- sum(data$peak_arid1a_dox == "yes") / N
  global_1b_off <- sum(data$peak_arid1b_veh == "yes") / N
  global_1a_off <- sum(data$peak_arid1a_veh == "yes") / N
  
  data %>%
    group_by(class_activity) %>%
    summarise(
      n_class         = n(),
      n_1b_peak_on    = sum(peak_arid1b_dox == "yes"),
      n_1a_peak_on    = sum(peak_arid1a_dox == "yes"),
      n_1b_peak_off   = sum(peak_arid1b_veh == "yes"),
      n_1a_peak_off   = sum(peak_arid1a_veh == "yes"),
      n_1b_nopeak_on  = sum(peak_arid1b_dox == "no"),
      n_1a_nopeak_on  = sum(peak_arid1a_dox == "no"),
      n_1b_nopeak_off = sum(peak_arid1b_veh == "no"),
      n_1a_nopeak_off = sum(peak_arid1a_veh == "no"),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      exp_1b_on  = n_class * global_1b_on,
      exp_1a_on  = n_class * global_1a_on,
      exp_1b_off = n_class * global_1b_off,
      exp_1a_off = n_class * global_1a_off,
      z_1b_on    = (n_1b_peak_on  - exp_1b_on)  / sqrt(exp_1b_on  * (1 - global_1b_on)),
      z_1a_on    = (n_1a_peak_on  - exp_1a_on)  / sqrt(exp_1a_on  * (1 - global_1a_on)),
      z_1b_off   = (n_1b_peak_off - exp_1b_off) / sqrt(exp_1b_off * (1 - global_1b_off)),
      z_1a_off   = (n_1a_peak_off - exp_1a_off) / sqrt(exp_1a_off * (1 - global_1a_off)),
      zdiff_on   = z_1b_on  - z_1a_on,
      zdiff_off  = z_1b_off - z_1a_off,
      fisher_on  = fisher.test(matrix(c(n_1b_peak_on,  n_1b_nopeak_on,
                                        n_1a_peak_on,  n_1a_nopeak_on),  nrow=2))$p.value,
      fisher_off = fisher.test(matrix(c(n_1b_peak_off, n_1b_nopeak_off,
                                        n_1a_peak_off, n_1a_nopeak_off), nrow=2))$p.value
    ) %>%
    ungroup() %>%
    mutate(
      padj_on  = p.adjust(fisher_on,  method = "BH"),
      padj_off = p.adjust(fisher_off, method = "BH"),
      sig_on   = sig_label(padj_on),
      sig_off  = sig_label(padj_off),
      origin   = "all",
      pE_class = factor(class_activity, levels = class_levels, labels = class_labels)
    )
}

zdiff_all <- compute_enrichment_zdiff_all(df_arid)

# ---- 3. Peak count panels
peak_counts <- df_arid %>%
  group_by(class_activity, novel_atac_peak) %>%
  summarise(
    both_on     = sum(peak_arid1b_dox == "yes" & peak_arid1a_dox == "yes"),
    only1b_on   = sum(peak_arid1b_dox == "yes" & peak_arid1a_dox == "no"),
    only1a_on   = sum(peak_arid1b_dox == "no"  & peak_arid1a_dox == "yes"),
    both_off    = sum(peak_arid1b_veh == "yes" & peak_arid1a_veh == "yes"),
    only1b_off  = sum(peak_arid1b_veh == "yes" & peak_arid1a_veh == "no"),
    only1a_off  = sum(peak_arid1b_veh == "no"  & peak_arid1a_veh == "yes"),
    .groups = "drop"
  ) %>%
  mutate(
    pE_class = factor(class_activity, levels = class_levels, labels = class_labels),
    origin   = factor(ifelse(novel_atac_peak == "known", "hema pEs", "novel pEs"),
                      levels = c("hema pEs", "novel pEs"))
  )

peak_counts_all <- peak_counts %>%
  group_by(pE_class) %>%
  summarise(across(c(both_on, only1b_on, only1a_on,
                     both_off, only1b_off, only1a_off), sum),
            .groups = "drop")

excl_max     <- max(c(peak_counts$only1b_on,  peak_counts$only1a_on,
                      peak_counts$only1b_off, peak_counts$only1a_off), na.rm = TRUE)
excl_max_all <- max(c(peak_counts_all$only1b_on, peak_counts_all$only1a_on,
                      peak_counts_all$only1b_off, peak_counts_all$only1a_off), na.rm = TRUE)

# ---- 4. Plot helpers
# Column (top) axis labels are now VERTICAL (angle = 90, hjust = 0) so they sit
# above the panels instead of tilting down into the heatmaps.
common_theme <- theme_pubr(base_size = 11) +
  theme(
    axis.text.x     = element_text(size = 9, angle = 90, hjust = 0, vjust = 0.5),
    axis.ticks.x    = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.title.x    = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.margin   = margin(4, 0, 0, 0),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(2, 2, 15, 2)
  )

# Fold-over-mean heatmap column
make_col <- function(data, value_col, x_label, show_y_text = FALSE) {
  d <- data %>% rename(val = all_of(value_col))
  ggplot(d, aes(x = origin, y = pE_class, fill = val)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(
      colors = c("#f7f7f7", "#d9d9d9", "#969696", "#525252"),
      limits = c(0, shared_max),
      name   = "Fold over mean"
    ) +
    scale_y_discrete(limits = rev) +
    scale_x_discrete(labels = c("hema pEs" = "hema", "novel pEs" = "novel"),
                     position = "top") +
    geom_hline(yintercept = c(4.5, 11.5), linetype = "dashed",
               color = "black", linewidth = 0.4) +
    labs(x = x_label, y = NULL) +
    guides(fill = guide_colorbar(barwidth = unit(3,"cm"), barheight = unit(0.5,"cm"),
                                 direction = "horizontal",
                                 title.position = "top", title.hjust = 0.5)) +
    common_theme +
    theme(axis.text.y = if (show_y_text) element_text(size=9) else element_blank())
}

# ARID1B - ARID1A difference column
make_diff_col <- function(data, diff_col, x_label) {
  d <- data %>%
    rename(val = all_of(diff_col)) %>%
    mutate(origin = factor(origin, levels = c("hema pEs", "novel pEs")))
  dlim <- max(abs(d$val), na.rm = TRUE)
  ggplot(d, aes(x = origin, y = pE_class, fill = val)) +
    geom_tile(color = NA) +
    scale_fill_gradient2(low="#2166ac", mid="white", high="#d6604d",
                         midpoint=0, limits=c(-dlim, dlim),
                         name="ARID1B-ARID1A\n(fold over mean)") +
    scale_y_discrete(limits = rev) +
    scale_x_discrete(labels = c("hema pEs" = "hema", "novel pEs" = "novel"),
                     position = "top") +
    geom_hline(yintercept = c(4.5, 11.5), linetype="dashed",
               color="black", linewidth=0.4) +
    labs(x = x_label, y = NULL) +
    guides(fill = guide_colorbar(barwidth=unit(3,"cm"), barheight=unit(0.5,"cm"),
                                 direction="horizontal",
                                 title.position="top", title.hjust=0.5)) +
    common_theme +
    theme(axis.text.y = element_blank())
}

# log2FC heatmap column (MA9-on vs MA9-off per protein)
make_log2fc_col <- function(data, value_col, x_label) {
  d <- data %>% rename(val = all_of(value_col))
  dlim <- max(abs(d$val), na.rm = TRUE)
  ggplot(d, aes(x = origin, y = pE_class, fill = val)) +
    geom_tile(color = NA) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d",
                         midpoint = 0, limits = c(-dlim, dlim),
                         name = "log2FC\n(MA9-on/off)") +
    scale_y_discrete(limits = rev) +
    scale_x_discrete(labels = c("hema pEs" = "hema", "novel pEs" = "novel"),
                     position = "top") +
    geom_hline(yintercept = c(4.5, 11.5), linetype = "dashed",
               color = "black", linewidth = 0.4) +
    labs(x = x_label, y = NULL) +
    guides(fill = guide_colorbar(barwidth = unit(3,"cm"), barheight = unit(0.5,"cm"),
                                 direction = "horizontal",
                                 title.position = "top", title.hjust = 0.5)) +
    common_theme +
    theme(axis.text.y = element_blank())
}

# log2FC difference column (ARID1B log2FC - ARID1A log2FC)
make_log2fc_diff_col <- function(data, diff_col, x_label) {
  d <- data %>% rename(val = all_of(diff_col))
  dlim <- max(abs(d$val), na.rm = TRUE)
  ggplot(d, aes(x = origin, y = pE_class, fill = val)) +
    geom_tile(color = NA) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d",
                         midpoint = 0, limits = c(-dlim, dlim),
                         name = "log2FC diff\n(1B \u2212 1A)") +
    scale_y_discrete(limits = rev) +
    scale_x_discrete(labels = c("hema pEs" = "hema", "novel pEs" = "novel"),
                     position = "top") +
    geom_hline(yintercept = c(4.5, 11.5), linetype = "dashed",
               color = "black", linewidth = 0.4) +
    labs(x = x_label, y = NULL) +
    guides(fill = guide_colorbar(barwidth = unit(3,"cm"), barheight = unit(0.5,"cm"),
                                 direction = "horizontal",
                                 title.position = "top", title.hjust = 0.5)) +
    common_theme +
    theme(axis.text.y = element_blank())
}

# Per-origin (hema/novel dodged) count bars
make_count_bar <- function(data, value_col, title_label, xlim = NULL) {
  d <- data %>% rename(val = all_of(value_col))
  ggplot(d, aes(x = val, y = pE_class, fill = origin)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, color = NA) +
    scale_fill_manual(values = c("hema pEs" = "#08519c", "novel pEs" = "#9ecae1"),
                      labels = c("hema", "novel"), name = "pE origin") +
    scale_y_discrete(limits = rev) +
    scale_x_continuous(position = "bottom",
                       limits = if (!is.null(xlim)) c(0, xlim) else NULL) +
    geom_hline(yintercept = c(4.5, 11.5), linetype = "dashed",
               color = "black", linewidth = 0.4) +
    labs(x = title_label, y = NULL) +
    guides(fill = guide_legend(position = "bottom")) +
    theme_pubr(base_size = 10) +
    theme(
      axis.title.x.bottom = element_text(size = 9, face = "bold"),
      axis.text.x.bottom  = element_text(size = 7),
      axis.ticks.x.bottom = element_line(),
      axis.text.x.top     = element_blank(),
      axis.ticks.x.top    = element_blank(),
      axis.text.y         = element_blank(),
      axis.ticks.y        = element_blank(),
      legend.position     = "bottom",
      legend.margin       = margin(4, 0, 0, 0),
      legend.text         = element_text(size = 8),
      plot.margin         = margin(2, 2, 15, 2)
    )
}

# Merged (all pEs) count bars
make_count_bar_all <- function(data, value_col, title_label, xlim = NULL, bar_color = "#636363") {
  d <- data %>% rename(val = all_of(value_col))
  ggplot(d, aes(x = val, y = pE_class)) +
    geom_bar(stat = "identity", width = 0.7, color = NA, fill = bar_color) +
    scale_y_discrete(limits = rev) +
    scale_x_continuous(position = "bottom",
                       limits = if (!is.null(xlim)) c(0, xlim) else NULL) +
    geom_hline(yintercept = c(4.5, 11.5), linetype = "dashed",
               color = "black", linewidth = 0.4) +
    labs(x = title_label, y = NULL) +
    theme_pubr(base_size = 10) +
    theme(
      axis.title.x.bottom = element_text(size = 9, face = "bold"),
      axis.text.x.bottom  = element_text(size = 7),
      axis.ticks.x.bottom = element_line(),
      axis.text.x.top     = element_blank(),
      axis.ticks.x.top    = element_blank(),
      axis.text.y         = element_blank(),
      axis.ticks.y        = element_blank(),
      legend.position     = "none",
      plot.margin         = margin(2, 2, 15, 2)
    )
}

# ---- 5. Build panels
# Fold-over-mean heatmaps
p_b_on  <- make_col(summary_df, "ARID1B_on",  "ARID1B", show_y_text = TRUE)
p_a_on  <- make_col(summary_df, "ARID1A_on",  "ARID1A") + theme(legend.position = "none")
p_b_off <- make_col(summary_df, "ARID1B_off", "ARID1B") + theme(legend.position = "none")
p_a_off <- make_col(summary_df, "ARID1A_off", "ARID1A") + theme(legend.position = "none")
p_diff_on  <- make_diff_col(summary_df, "diff_on",  "1B vs 1A")
p_diff_off <- make_diff_col(summary_df, "diff_off", "1B vs 1A") + theme(legend.position = "none")

# log2FC heatmaps
p_lfc_1b   <- make_log2fc_col(summary_df, "log2fc_1b",   "ARID1B")
p_lfc_1a   <- make_log2fc_col(summary_df, "log2fc_1a",   "ARID1A") + theme(legend.position = "none")
p_lfc_diff <- make_log2fc_diff_col(summary_df, "log2fc_diff", "1B vs 1A") + theme(legend.position = "none")

# Per-origin count bars (MA9-on)
p_only1b_on <- make_count_bar(peak_counts, "only1b_on", "1B only", xlim = excl_max) +
  theme(legend.position = "none")
p_only1a_on <- make_count_bar(peak_counts, "only1a_on", "1A only", xlim = excl_max) +
  theme(legend.position = "none")
p_both_on   <- make_count_bar(peak_counts, "both_on",   "Both")

# Per-origin count bars (MA9-off)
p_only1b_off <- make_count_bar(peak_counts, "only1b_off", "1B only", xlim = excl_max) +
  theme(legend.position = "none")
p_only1a_off <- make_count_bar(peak_counts, "only1a_off", "1A only", xlim = excl_max) +
  theme(legend.position = "none")
p_both_off   <- make_count_bar(peak_counts, "both_off",   "Both") +
  theme(legend.position = "none")

# Merged (all pEs) count bars
pa_only1b_on  <- make_count_bar_all(peak_counts_all, "only1b_on",  "1B only", xlim = excl_max_all, bar_color = "#4d4d4d")
pa_only1a_on  <- make_count_bar_all(peak_counts_all, "only1a_on",  "1A only", xlim = excl_max_all, bar_color = "#b0b0b0")
pa_both_on    <- make_count_bar_all(peak_counts_all, "both_on",    "Both",    bar_color = "#4292c6")
pa_only1b_off <- make_count_bar_all(peak_counts_all, "only1b_off", "1B only", xlim = excl_max_all, bar_color = "#4d4d4d")
pa_only1a_off <- make_count_bar_all(peak_counts_all, "only1a_off", "1A only", xlim = excl_max_all, bar_color = "#b0b0b0")
pa_both_off   <- make_count_bar_all(peak_counts_all, "both_off",   "Both",    bar_color = "#4292c6")

# Fisher's exact z-score-diff panel (with significance stars)
or_long <- zdiff_all %>%
  select(pE_class, zdiff_on, zdiff_off, sig_on, sig_off) %>%
  pivot_longer(cols = c(zdiff_on, zdiff_off),
               names_to = "condition", values_to = "zdiff") %>%
  mutate(
    MA9_state = ifelse(condition == "zdiff_on", "MA9-on", "MA9-off"),
    sig       = ifelse(MA9_state == "MA9-on", sig_on, sig_off),
    MA9_state = factor(MA9_state, levels = c("MA9-on", "MA9-off"))
  )

or_lim <- max(abs(or_long$zdiff), na.rm = TRUE)

p_or <- ggplot(or_long, aes(x = MA9_state, y = pE_class, fill = zdiff)) +
  geom_tile(color = NA) +
  geom_text(aes(label = sig), size = 5, vjust = 0.75) +
  scale_fill_gradient2(low="#2166ac", mid="white", high="#d6604d",
                       midpoint=0, limits=c(-or_lim, or_lim),
                       name="z-score diff\nARID1B vs ARID1A") +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  geom_hline(yintercept = c(4.5, 11.5), linetype="dashed",
             color="black", linewidth=0.4) +
  labs(x = NULL, y = NULL, title = "Fisher's exact") +
  guides(fill = guide_colorbar(barwidth=unit(3,"cm"), barheight=unit(0.5,"cm"),
                               direction="horizontal",
                               title.position="top", title.hjust=0.5)) +
  common_theme +
  theme(
    plot.title  = element_text(size=12, hjust=0.5, face="bold"),
    axis.text.x = element_text(size=12, angle=90, hjust=0, vjust=0.5),
    axis.text.y = element_blank()
  )

# ---- 6. Assemble final figure
title_on  <- ggplot() +
  annotate("text", x=0.5, y=0.5, label="MA9-on",
           size=4.5, fontface="bold", color="#cc0000") + theme_void()
title_off <- ggplot() +
  annotate("text", x=0.5, y=0.5, label="MA9-off",
           size=4.5, fontface="bold", color="black") + theme_void()
title_lfc <- ggplot() +
  annotate("text", x=0.5, y=0.5, label="log2FC (MA9-on vs off)",
           size=4.5, fontface="bold", color="black") + theme_void()
title_or  <- ggplot() + theme_void()
spacer    <- plot_spacer()

title_row   <- title_on | spacer | title_off | spacer | title_lfc | spacer | title_or
heatmap_row <- p_b_on | p_a_on | p_diff_on | p_only1b_on | p_only1a_on | p_both_on |
  pa_only1b_on | pa_only1a_on | pa_both_on | spacer |
  p_b_off | p_a_off | p_diff_off | p_only1b_off | p_only1a_off | p_both_off |
  pa_only1b_off | pa_only1a_off | pa_both_off | spacer |
  p_lfc_1b | p_lfc_1a | p_lfc_diff | spacer | p_or

final <- title_row / heatmap_row +
  plot_layout(heights = c(0.08, 1),
              widths  = c(1.4, 1, 1, 1, 1, 1, 0.7, 0.7, 0.7, 0.1,
                          1,   1, 1, 1, 1, 1, 0.7, 0.7, 0.7, 0.1,
                          1,   1, 1, 0.1, 1.2))

print(final)
ggsave("./Figure_S4A/arid_percolumn_heatmap.pdf", final, width = 33, height = 7.5)
cat("Done. Figure saved to ./Figure_S4A/arid_percolumn_heatmap.pdf\n")

############################## Figure S4B ############################### ----
# CONCORDANT GENES FOR GSEA IN R
# Creates gene lists organized by bin group for enrichment analysis
library(dplyr)

# STEP 1: Load and prepare data

df <- df_selected_4

# Get unique enhancers
enh <- df[!duplicated(df$name_dhs_dox), ]
enh <- enh[!is.na(enh$logFC_expr_closest), ]

# Aggregate to gene level
gene_stats <- enh %>%
  group_by(closest_gene_dox) %>%
  summarise(
    n_enhancers = n(),
    mean_atac_logFC = mean(logFC_atac, na.rm = TRUE), 
    gene_expr = logFC_expr_closest[1],
    class_expr = class_expr_closest[1],
    .groups = "drop"
  )

# STEP 2: Filter concordant genes and create bins

n_bins <- 2
n2_bins <- 1

# Concordant UP: atac > 0 AND class_expr = 'up'
concordant_up <- gene_stats %>%
  filter(mean_atac_logFC > 0 & class_expr == "up",gene_expr>1) %>%
  mutate(bin_num = ntile(mean_atac_logFC, n_bins))

disconcordant_up <- gene_stats %>%
  filter(mean_atac_logFC < 0 & class_expr == "up",gene_expr>1) %>%
  mutate(bin_num = ntile(mean_atac_logFC, n2_bins))

# Concordant DOWN: atac < 0 AND class_expr = 'down'
concordant_down <- gene_stats %>%
  filter(mean_atac_logFC < 0 & class_expr == "down",gene_expr< -1) %>%
  mutate(bin_num = ntile(mean_atac_logFC, n_bins))

disconcordant_down <- gene_stats %>%
  filter(mean_atac_logFC > 0 & class_expr == "down",gene_expr< -1) %>%
  mutate(bin_num = ntile(mean_atac_logFC, n2_bins))

# Add bin labels
bin_info_up <- concordant_up %>%
  group_by(bin_num) %>%
  summarise(
    min_atac = min(mean_atac_logFC),
    max_atac = max(mean_atac_logFC),
    .groups = "drop"
  ) %>%
  mutate(bin_label = sprintf("UP_bin%02d_atac_%.2f_to_%.2f", bin_num, min_atac, max_atac))

bin_info_up_dis <-disconcordant_up %>%
  group_by(bin_num) %>%
  summarise(
    min_atac = min(mean_atac_logFC),
    max_atac = max(mean_atac_logFC),
    .groups = "drop"
  ) %>%
  mutate(bin_label = sprintf("UP_bin%02d_atac_%.2f_to_%.2f", bin_num, min_atac, max_atac))

bin_info_down <- concordant_down %>%
  group_by(bin_num) %>%
  summarise(
    min_atac = min(mean_atac_logFC),
    max_atac = max(mean_atac_logFC),
    .groups = "drop"
  ) %>%
  mutate(bin_label = sprintf("DOWN_bin%02d_atac_%.2f_to_%.2f", bin_num, min_atac, max_atac))

bin_info_down_dis <- disconcordant_down %>%
  group_by(bin_num) %>%
  summarise(
    min_atac = min(mean_atac_logFC),
    max_atac = max(mean_atac_logFC),
    .groups = "drop"
  ) %>%
  mutate(bin_label = sprintf("DOWN_bin%02d_atac_%.2f_to_%.2f", bin_num, min_atac, max_atac))

concordant_up <- concordant_up %>% left_join(bin_info_up, by = "bin_num")
concordant_down <- concordant_down %>% left_join(bin_info_down, by = "bin_num")
concordant_up_dis <- disconcordant_up %>% left_join(bin_info_up_dis, by = "bin_num")
concordant_down_dis <- disconcordant_down %>% left_join(bin_info_down_dis, by = "bin_num")

# STEP 3: Create gene lists for GSEA

# --- CONCORDANT UP GENE LISTS ---
concordant_up_genelist <- split(concordant_up$closest_gene_dox, concordant_up$bin_label)

# --- CONCORDANT DOWN GENE LISTS ---
concordant_down_genelist <- split(concordant_down$closest_gene_dox, concordant_down$bin_label)
# --- dis-CONCORDANT UP GENE LISTS ---
disconcordant_up_genelist <- split(concordant_up_dis$closest_gene_dox, concordant_up_dis$bin_label)

# --- dis-CONCORDANT DOWN GENE LISTS ---
disconcordant_down_genelist <- split(concordant_down_dis$closest_gene_dox, concordant_down_dis$bin_label)

# --- COMBINED LIST ---
df_up<-filter(df_promoter, 
              logFC_expr>1, 
              class_expr=="up")
df_down<-filter(df_promoter, 
                logFC_expr< -1, 
                class_expr=="down")

# arid1b kd
df_up_arid1bKD<-filter(df_promoter, 
                       logFC_expr_arid1bKD>0.3,
                       class_expr_arid1bKD=="up"
)

df_down_arid1bKD<-filter(df_promoter, 
                         logFC_expr_arid1bKD < -0.3,
                         class_expr_arid1bKD=="down"
)

up_ma9ON<-df_up$name
down_ma9ON<-df_down$name
up_arid1bKD<-df_up_arid1bKD$name
down_arid1bKD<-df_down_arid1bKD$name

list<-list(
  "up_ma9ON"=up_ma9ON, 
  "down_ma9ON"=down_ma9ON,
  "down_arid1bKD"=down_arid1bKD,
  "up_arid1bKD"=up_arid1bKD
)


all_genelists <- c(list,
                   concordant_up_genelist, disconcordant_up_genelist,
                   concordant_down_genelist,disconcordant_down_genelist
)
dplyr::glimpse(all_genelists)
# --- run geneset enrichment ---
library(msigdbr)
# check what are available in msigdbr
msigdbr_collections()
msigdbr_collections() %>%
  filter(gs_collection == "C2")

# combined 1 (final)
Mm <- bind_rows(
  msigdbr(species = "Mus musculus") %>%
    filter(gs_collection == "C2") %>%
    dplyr::select(gs_name, ncbi_gene, gene_symbol),
  
  msigdbr(species = "Mus musculus") %>%
    filter(gs_collection == "C5", gs_subcollection == "GO:BP") %>%
    dplyr::select(gs_name, ncbi_gene, gene_symbol),
  
  msigdbr(species = "Mus musculus") %>% 
    filter(gs_collection == "H") %>%
    dplyr::select(gs_name, ncbi_gene, gene_symbol)
)

# RUN
res <- compareCluster(all_genelists, enricher, TERM2GENE=Mm[,c(1,3)])

# ploting
df <- as.data.frame(res@compareClusterResult)

# --- Source data for reviewers: full compareCluster ORA table (p.adjust < 0.05) ---
ora_out <- as.data.frame(res@compareClusterResult)
ora_out <- ora_out[!is.na(ora_out$p.adjust) & ora_out$p.adjust < 0.05, ]
write.csv(ora_out, "FigS4B_ORA.csv", row.names = FALSE)
#write.csv(df, "ORA_bin2_atac_compareCluster_full.csv", row.names = FALSE)

N <- 5  # top N (20) per cluster to define the union (set 20 if you can tolerate more rows)

# 1) top N per cluster (ONLY to define union, not to cut the universe globally)
df_topN <- df %>%
  filter(!is.na(Description)) %>%
  group_by(Cluster) %>%
  arrange(p.adjust, desc(Count)) %>%
  slice_head(n = N) %>%
  ungroup()

# 2) union of terms across clusters
union_terms <- unique(df_topN$Description)

# 3) keep ALL results for those union terms (so terms appear across clusters if present)
df_union <- df %>%
  filter(Description %in% union_terms) %>%
  filter(!is.na(Description))

# 4) plot object: copy + replace slot (do not overwrite the original)
res_plot <- res
res_plot@compareClusterResult <- df_union

# 5) dotplot provides the clean layout/order
p <- dotplot(res_plot, showCategory = length(union_terms))

# 6) extract dotplot data and apply your styling
plot_data <- as.data.frame(p$data)
plot_data$log10_padjust <- -log10(plot_data$p.adjust)

ggplot(plot_data, aes(x = Cluster, y = Description, size = Count)) +
  geom_point(shape = 21, aes(fill = log10_padjust), color = "black", stroke = 0.5) +
  scale_fill_gradient(low = "#9ECAE1", high = "#08519C",, name = "-log10(p.adjust)") +
  theme_bw() +
  labs(title = "bin2.ccdt_bin1.dist_6groups0.3_atac_H_C2_BP_top10",
       x = "Cluster", y = "Pathway") +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 8),
        axis.text.y = element_text(size = 6)) +
  scale_size_continuous(name = "Gene Count", range = c(3, 10), guide = "legend")

# --- Gene Set ORA Visualization - Dotplot and Heatmap ---

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# 1. Load compareCluster results
df <- as.data.frame(res@compareClusterResult)

# Check clusters
cat("Clusters found in data:\n")
print(unique(df$Cluster))

# 2. Define cluster order and display labels

# Your actual cluster names from compareCluster
cluster_order <- c(
  "up_ma9ON",
  "down_arid1bKD",
  "UP_bin02_atac_0.43_to_3.41",
  "UP_bin01_atac_0.00_to_0.43",
  #"DOWN_bin01_atac_0.01_to_2.63",
  
  "down_ma9ON",
  "up_arid1bKD",
  "DOWN_bin01_atac_-3.64_to_-0.57",
  "DOWN_bin02_atac_-0.57_to_-0.00"
  #"UP_bin01_atac_-2.21_to_-0.00"
)

# Display labels
cluster_labels <- c(
  "up_ma9ON" = "RNA up",
  
  "down_ma9ON" = "RNA down",
  "UP_bin02_atac_0.43_to_3.41" = "ATAC++ RNA up",
  "UP_bin01_atac_0.00_to_0.43" = "ATAC + RNA up",
  #"DOWN_bin01_atac_0.01_to_2.63" = "Enh+\nRNA Down",
  
  "DOWN_bin01_atac_-3.64_to_-0.57" = "ATAC-- RNA down",
  "DOWN_bin02_atac_-0.57_to_-0.00" = "ATAC- RNA down",
  #"UP_bin01_atac_-2.21_to_-0.00" = "Enh-\nRNA Up",
  
  "up_arid1bKD"="RNA up (Arid1b KD)",
  "down_arid1bKD"="RNA down (Arid1b KD)"
)


# Keep only clusters that exist in data
cluster_order <- cluster_order[cluster_order %in% unique(df$Cluster)]
cat("\nClusters to plot:\n")
print(cluster_order)

# 3. Define selected gene sets

selected_genesets <- tribble(
  ~geneset_id, ~display_name, ~theme,
  
  # 1) LEUKEMIA (Top significant in Upregulated) 
  "HUANG_AML_LSC47", "LSC47 (Huang)", "Leukemia",
  "KEGG_MEDICUS_VARIANT_MLL_ENL_FUSION_TO_TRANSCRIPTIONAL_ACTIVATION", "MLL-ENL fusion targets (KEGG)", "Leukemia",
  "WANG_IMMORTALIZED_BY_HOXA9_AND_MEIS1_DN", "HOXA9/MEIS1 vs HOXA9 UP (Wang)", "Leukemia",
  "HESS_TARGETS_OF_HOXA9_AND_MEIS1_UP", "HOXA9/MEIS1 targets UP (Hess)", "Leukemia",
  # ESC-like stemness programs (co-opted by cancer/leukemia)
  #"WONG_EMBRYONIC_STEM_CELL_CORE", "Embryonic stem cell core (Wong)", "Leukemia",
  #"BHATTACHARYA_EMBRYONIC_STEM_CELL", "Embryonic stem cell (Bhattacharya)", "Leukemia",
  #"RAMALHO_STEMNESS_UP", "Stemness UP (Ramalho)", "Leukemia",
  
  # 2) HEMATOPOIESIS 
  "JAATINEN_HEMATOPOIETIC_STEM_CELL_UP", "HSC UP (Jaatinen)", "Hematopoiesis",
  # Complete Ivanova hierarchy (stem cell to mature)
  #"IVANOVA_HEMATOPOIESIS_STEM_CELL", "Stem cell (Ivanova)", "Hematopoiesis",
  "IVANOVA_HEMATOPOIESIS_STEM_CELL_LONG_TERM", "LT-HSC (Ivanova)", "Hematopoiesis",
  #"IVANOVA_HEMATOPOIESIS_STEM_CELL_SHORT_TERM", "ST-HSC (Ivanova)", "Hematopoiesis",
  "IVANOVA_HEMATOPOIESIS_EARLY_PROGENITOR", "Early progenitor (Ivanova)", "Hematopoiesis",
  #"IVANOVA_HEMATOPOIESIS_INTERMEDIATE_PROGENITOR", "Intermediate progenitor (Ivanova)", "Hematopoiesis",
  "IVANOVA_HEMATOPOIESIS_MATURE_CELL", "Mature cell (Ivanova)", "Hematopoiesis",
  # Other hematopoiesis CGP sets
  "BROWN_MYELOID_CELL_DEVELOPMENT_UP", "Myeloid development UP (Brown)", "Hematopoiesis",
  "KAMIKUBO_MYELOID_CEBPA_NETWORK", "Myeloid CEBPA network (Kamikubo)", "Hematopoiesis",
  # Lymphoid markers (to show myeloid vs lymphoid)
  #"LEE_EARLY_T_LYMPHOCYTE_UP", "Early T lymphocyte UP (Lee)", "Hematopoiesis",
  #"MORI_MATURE_B_LYMPHOCYTE_UP", "Mature B lymphocyte UP (Mori)", "Hematopoiesis",
  
  # 3) BIOLOGICAL PATHWAYS (GO:BP) 
  # Upregulated - Development/Morphogenesis
  "GOBP_EMBRYONIC_ORGAN_DEVELOPMENT", "Embryonic organ development", "Biological pathways",
  "GOBP_RENAL_SYSTEM_DEVELOPMENT", "Skeletal system development", "Biological pathways",
  "GOBP_SENSORY_ORGAN_MORPHOGENESIS","Sensory organ morphogenesis","Biological pathways",
  "GOBP_SKELETAL_SYSTEM_MORPHOGENESIS", "Renal system morphogenesis", "Biological pathways",
  # Upregulated - Metabolism
  "GOBP_AMINO_ACID_METABOLIC_PROCESS", "Amino acid metabolism", "Biological pathways",
  "GOBP_ALPHA_AMINO_ACID_METABOLIC_PROCESS", "Alpha amino acid metabolic process", "Biological pathways",
  "GOBP_ORGANIC_ACID_BIOSYNTHETIC_PROCESS", "Organic acid biosynthesis", "Biological pathways",
  
  # Downregulated (immune & migration)
  "GOBP_CELL_CHEMOTAXIS", "Cell chemotaxis", "Biological pathways",
  "GOBP_LEUKOCYTE_MIGRATION", "Leukocyte migration", "Biological pathways",
  "GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE", "Inflammatory response regulation", "Biological pathways",
  "GOBP_PHAGOCYTOSIS", "Phagocytosis", "Biological pathways",
  "GOBP_REGULATION_OF_LEUKOCYTE_PROLIFERATION", "Regulation of leukocyte proliferation", "Biological pathways",
  
  
  # 4) HALLMARK - Top significant each direction 
  # Upregulated
  "HALLMARK_MYC_TARGETS_V2", "MYC targets V2", "Hallmark",
  "HALLMARK_MTORC1_SIGNALING", "mTORC1 signaling", "Hallmark",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "Unfolded protein response", "Hallmark",
  # Downregulated
  "HALLMARK_INFLAMMATORY_RESPONSE", "Inflammatory response", "Hallmark",
  "HALLMARK_IL2_STAT5_SIGNALING", "IL2/STAT5 signaling", "Hallmark",
  "HALLMARK_KRAS_SIGNALING_UP", "KRAS signaling UP", "Hallmark",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNF/NF-kB signaling", "Hallmark",
  "HALLMARK_COMPLEMENT", "Complement", "Hallmark"
)

# 4. Prepare data for visualization
plot_data <- df %>%
  filter(ID %in% selected_genesets$geneset_id) %>%
  filter(Cluster %in% cluster_order) %>%
  left_join(selected_genesets, by = c("ID" = "geneset_id")) %>%
  mutate(
    log10_padj = -log10(p.adjust),
    log10_padj_capped = pmin(log10_padj, 10),
    sig = case_when(
      p.adjust < 0.001 ~ "***",
      p.adjust < 0.01 ~ "**",
      p.adjust < 0.05 ~ "*",
      TRUE ~ ""
    ),
    theme = factor(theme, levels = c("Leukemia", "Hematopoiesis", "Biological pathways", "Hallmark")),
    Cluster = factor(Cluster, levels = cluster_order)
  )

# Order gene sets by theme
geneset_order <- selected_genesets %>%
  mutate(theme = factor(theme, levels = c("Leukemia", "Hematopoiesis", "Biological pathways", "Hallmark"))) %>%
  arrange(theme) %>%
  pull(display_name)

plot_data <- plot_data %>%
  mutate(display_name = factor(display_name, levels = rev(geneset_order)))

# 5. DOTPLOT

# Create scaffold
scaffold <- expand.grid(
  display_name = geneset_order,
  Cluster = cluster_order,
  stringsAsFactors = FALSE
) %>%
  left_join(selected_genesets %>% dplyr::select(display_name, theme), by = "display_name") %>%
  mutate(
    theme = factor(theme, levels = c("Leukemia", "Hematopoiesis", "Biological pathways", "Hallmark")),
    display_name = factor(display_name, levels = rev(geneset_order)),
    Cluster = factor(Cluster, levels = cluster_order)
  )

# Filter significant results
plot_data_sig <- plot_data %>%
  filter(p.adjust < 0.05)
# dotplot
Fiugre_S4B <- ggplot() +
  geom_point(data = scaffold, aes(x = Cluster, y = display_name), alpha = 0) +
  geom_point(data = plot_data_sig, aes(x = Cluster, y = display_name, size = Count, fill = log10_padj_capped), 
             shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient(
    low = "white", high = "darkgreen",
    limits = c(0, 10),
    name = "-log10\n(p.adjust)"
  ) +
  scale_size_continuous(range = c(2, 10), name = "Gene count") +
  scale_x_discrete(labels = cluster_labels[cluster_order]) +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 7, color = "black", face = "bold", angle = 45, hjust = 1),
    axis.title = element_blank(),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 10),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(Fiugre_S4B)

############################## Figure 6D ############################### ----
# GSEA: KMT2A-MLLT3 dTag signatures vs ARID1B-high/low patient ranking
# Reproduces the two stacked enrichment panels:
#   (1) KMT2A-MLLT3-dependent  transcription  (dTag 24h, log2FC < -1)  [NES > 0]
#   (2) KMT2A-MLLT3-suppressed transcription  (dTag 24h, log2FC >  1)  [NES < 0]
#
# Inputs (place both in the working directory):
#   - preranked_rna_arid1b_HIGH_minus_LOW_mllr.tsv   (columns: gene, score)
#   - human_dtag.csv   (columns: name, class_dTag1D_human, logFC_dTag1D_human, ...)
#
# Output:
#   - gsea_results_dependent_suppressed.csv
#   - gsea_kmt2a_mllt3_arid1b.pdf  (both panels)  + individual PDFs
#
# Requires fgsea >= 1.18 (for plotEnrichmentData()).

library(readr)
library(dplyr)
library(tibble)
library(fgsea)
library(ggplot2)
library(patchwork)

# setwd("/path/to/folder/with/the/two/input/files")   # <- point this at your data

infile_preranked <- "preranked_rna_arid1b_HIGH_minus_LOW_mllr.tsv"
infile_dtag      <- "human_dtag.csv"

# ---- 1. Build the ranked list (named, decreasing) 
df_patient <- read_tsv(infile_preranked, show_col_types = FALSE)
stopifnot(all(c("gene", "score") %in% names(df_patient)))

ranks_df <- df_patient %>%
  transmute(gene = toupper(gene), score = as.numeric(score)) %>%
  filter(!is.na(gene), is.finite(score))

# collapse duplicate symbols by median, sort decreasing
ranks <- tapply(ranks_df$score, ranks_df$gene, median)
ranks <- sort(ranks, decreasing = TRUE)

# tiny deterministic jitter only to break exact ties
set.seed(42)
ranks <- ranks + rnorm(length(ranks), sd = 1e-9)

# ---- 2. Define the two signatures explicitly (no overwrite toggle) 
dtag <- read_csv(infile_dtag, show_col_types = FALSE)
dtag$name <- toupper(dtag$name)

genes_dep <- dtag %>%
  filter(class_dTag1D_human == "down", logFC_dTag1D_human < -1) %>%
  pull(name) %>% unique()

genes_sup <- dtag %>%
  filter(class_dTag1D_human == "up", logFC_dTag1D_human > 1) %>%
  pull(name) %>% unique()

pathways <- list(
  dependent  = intersect(genes_dep, names(ranks)),
  suppressed = intersect(genes_sup, names(ranks))
)

cat(sprintf("dependent set: %d genes (%d in ranking)\n",
            length(genes_dep), length(pathways$dependent)))
cat(sprintf("suppressed set: %d genes (%d in ranking)\n",
            length(genes_sup), length(pathways$suppressed)))

# ---- 3. Run fgsea 
set.seed(42)
res <- fgseaMultilevel(
  pathways = pathways,
  stats    = ranks,
  minSize  = 5,
  maxSize  = 5000,
  eps      = 0          # eps = 0 for accurate small p-values
)

res_out <- res %>%
  dplyr::select(pathway, size, ES, NES, pval, padj, leadingEdge) %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write_csv(res_out, "gsea_results_dependent_suppressed.csv")
print(res %>% dplyr::select(pathway, size, NES, padj))

# ---- 4. Custom enrichment-panel builder 
make_gsea_panel <- function(genes, ranks, res_row, title, subtitle,
                            high_label = "ARID1B-high",
                            low_label  = "ARID1B-low",
                            n_bins = 512) {
  
  pd   <- fgsea::plotEnrichmentData(pathway = genes, stats = ranks)
  NES  <- res_row$NES
  padj <- res_row$padj
  N    <- max(pd$stats$rank)
  
  curve_col <- if (NES > 0) "#d7301f" else "#3182bd"   # red / blue
  
  # NES / padj annotation: upper-right for positive ES, mid-left for negative
  if (NES > 0) {
    ann_x <- 0.50 * N; ann_y <- 0.62 * pd$posES
  } else {
    ann_x <- 0.03 * N; ann_y <- 0.55 * pd$negES
  }
  ann_txt <- sprintf("NES = %.2f\npadj = %s",
                     NES, formatC(padj, format = "e", digits = 1))
  
  # (a) running enrichment score
  p_es <- ggplot(pd$curve, aes(rank, ES)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_line(color = curve_col, linewidth = 1) +
    annotate("text", x = ann_x, y = ann_y, label = ann_txt,
             hjust = 0, vjust = 1, size = 4, lineheight = 0.95) +
    labs(y = "Enrichment score", title = title, subtitle = subtitle) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, N)) +
    theme_classic(base_size = 12) +
    theme(
      axis.title.x  = element_blank(),
      axis.text.x   = element_blank(),
      axis.ticks.x  = element_blank(),
      axis.line.x   = element_blank(),
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      plot.margin   = margin(4, 8, 0, 4)
    )
  
  # (b) gene hit ticks
  p_ticks <- ggplot(pd$ticks) +
    geom_segment(aes(x = rank, xend = rank, y = 0, yend = 1), linewidth = 0.2) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, N)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 8, 0, 4))
  
  # (c) ranked-metric gradient strip (binned to keep the file light) + labels
  strip <- pd$stats %>%
    mutate(bin = cut(rank, breaks = n_bins, labels = FALSE)) %>%
    group_by(bin) %>%
    summarise(rank = mean(rank), stat = mean(stat), .groups = "drop")
  
  p_bar <- ggplot(strip, aes(x = rank, y = 0, fill = stat)) +
    geom_tile() +
    scale_fill_gradient2(low = "#3182bd", mid = "white", high = "#d7301f",
                         midpoint = 0, guide = "none") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, N)) +
    scale_y_continuous(limits = c(-1.6, 0.5), expand = c(0, 0)) +
    annotate("text", x = 0.01 * N, y = -1.0, label = high_label,
             hjust = 0, fontface = "italic", colour = "#d7301f", size = 3.6) +
    annotate("text", x = 0.99 * N, y = -1.0, label = low_label,
             hjust = 1, fontface = "italic", colour = "black", size = 3.6) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 8, 4, 4))
  
  p_es / p_ticks / p_bar + plot_layout(heights = c(1, 0.14, 0.16))
}

# ---- 5. Build the two panels and assemble
panel_dep <- make_gsea_panel(
  genes    = pathways$dependent,
  ranks    = ranks,
  res_row  = res[res$pathway == "dependent", ],
  title    = "KMT2A-MLLT3-dependent transcription",
  subtitle = "(Olsen, human leukemia, dTag 24h, log2FC < -1)"
)

panel_sup <- make_gsea_panel(
  genes    = pathways$suppressed,
  ranks    = ranks,
  res_row  = res[res$pathway == "suppressed", ],
  title    = "KMT2A-MLLT3-suppressed transcription",
  subtitle = "(Olsen, human leukemia, dTag 24h, log2FC > 1)"
)

# wrap_elements() nests each composed panel as one unit so its internal
# plot_layout heights survive; a plain `panel_dep / panel_sup` flattens them
# and triggers the wrap_dims() "need 4 panels" error.
combined <- wrap_elements(panel_dep) / wrap_elements(panel_sup)

print(combined)

############################## Figure 6E ############################### ----
# Figure E - GSEA enrichment panels, ARID1B-high vs low (De novo KMT2A-r AML)
#   Top:    Hematopoietic stem cell up (Jaatinen)   [C2:CGP]
#           AML LSC47 (Huang)                        [C2:CGP]
#   Bottom: G2M Checkpoint (Hallmark)                [H]
#           E2F Targets (Hallmark)                   [H]
#
# Each panel's NES/padj is taken from a FULL-collection fgsea run, so the
# BH-adjusted padj reproduces the reported values (e.g. Jaatinen padj 4.69e-24
# only appears when adjusted across all of C2:CGP). NES is collection-independent.
#
# Input:  rna_arid1b_low_high_mllr.tsv   (columns: Gene, "Log2 Ratio" = LOW-HIGH)
# Output: figureE_gsea_arid1b_4panels.pdf  + figureE_gsea_stats.csv
#
# Requires fgsea >= 1.18 (plotEnrichmentData).

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
})

# setwd("/path/to/folder/with/the/input/file")   # <- point this at your data
infile <- "rna_arid1b_low_high_mllr.tsv"

# ---- 1. Ranked list (flip LOW-HIGH -> HIGH-LOW so NES > 0 = ARID1B-high) 
df <- fread(infile, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)

ranks <- df %>%
  dplyr::select(Gene, `Log2 Ratio`) %>%
  filter(!is.na(Gene), !is.na(`Log2 Ratio`)) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  mutate(`Log2 Ratio (HIGH-LOW)` = -`Log2 Ratio`)

ranks_vec <- ranks$`Log2 Ratio (HIGH-LOW)`
names(ranks_vec) <- toupper(ranks$Gene)

set.seed(42)
ord <- order(ranks_vec, decreasing = TRUE, na.last = NA)
ranks_vec <- ranks_vec[ord]
eps_jit <- rnorm(length(ranks_vec), sd = 1e-9)
ranks_vec <- ranks_vec + eps_jit[rank(ranks_vec, ties.method = "first")] * 1e-12

# ---- 2. Gene-set collections -
msig_h   <- msigdbr(species = "Homo sapiens", category = "H")
msig_cgp <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
pathways_h   <- split(msig_h$gene_symbol,   msig_h$gs_name)
pathways_cgp <- split(msig_cgp$gene_symbol, msig_cgp$gs_name)

# ---- AML LSC47 (Huang): MSigDB C2:CGP set HUANG_AML_LSC47 
lsc_id <- "HUANG_AML_LSC47"
if (!lsc_id %in% names(pathways_cgp)) {
  hits <- grep("HUANG|LSC", names(pathways_cgp), value = TRUE)
  stop("'", lsc_id, "' not found in C2:CGP (msigdbr naming may differ). Candidates: ",
       if (length(hits)) paste(hits, collapse = ", ") else "none")
}

# ---- 3. fgsea per collection (full families -> correct padj) 
set.seed(1234)
fg_h   <- fgsea(pathways_h,   ranks_vec, minSize = 15, maxSize = 500, eps = 0)
set.seed(1234)
fg_cgp <- fgsea(pathways_cgp, ranks_vec, minSize = 15, maxSize = 500, eps = 0)

# --- Source data for reviewers: full GSEA tables (padj < 0.05) ---
gsea_out <- as.data.frame(fg_h)
gsea_out$leadingEdge <- sapply(gsea_out$leadingEdge, paste, collapse = ", ")
gsea_out <- gsea_out[!is.na(gsea_out$padj) & gsea_out$padj < 0.05, ]
write.csv(gsea_out, "Fig6E_Hallmark_GSEA.csv", row.names = FALSE)

gsea_out <- as.data.frame(fg_cgp)
gsea_out$leadingEdge <- sapply(gsea_out$leadingEdge, paste, collapse = ", ")
gsea_out <- gsea_out[!is.na(gsea_out$padj) & gsea_out$padj < 0.05, ]
write.csv(gsea_out, "Fig6E_C2_CGP_GSEA.csv", row.names = FALSE)

get_stat <- function(fg, id) {
  row <- fg[fg$pathway == id, ]
  if (nrow(row) == 0) stop("Pathway not found in results: ", id)
  list(NES = row$NES[1], padj = row$padj[1])
}

# ---- 4. Panel builder (shared styling) 
make_panel <- function(genes, ranks, NES, padj, title,
                       show_xlabels = TRUE,
                       high_label = "ARID1B-high", low_label = "ARID1B-low",
                       curve_col = "#cb181d", n_bins = 512) {
  
  genes <- intersect(toupper(genes), names(ranks))
  pd <- fgsea::plotEnrichmentData(pathway = genes, stats = ranks)
  N  <- max(pd$stats$rank)
  
  padj_txt <- if (padj < 1e-3) formatC(padj, format = "e", digits = 2)
  else             formatC(padj, format = "f", digits = 3)
  ann <- sprintf("NES=%.2f\npadj=%s", NES, padj_txt)
  
  # (a) running enrichment score
  p_es <- ggplot(pd$curve, aes(rank, ES)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_line(colour = curve_col, linewidth = 0.9) +
    annotate("text", x = 0.30 * N, y = 0.55 * pd$posES, label = ann,
             hjust = 0, vjust = 1, size = 3.4, lineheight = 0.95) +
    labs(title = title, y = "Enrichment score") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, N)) +
    theme_classic(base_size = 11) +
    theme(
      axis.title.x = element_blank(), axis.text.x = element_blank(),
      axis.ticks.x = element_blank(), axis.line.x = element_blank(),
      axis.title.y = element_text(size = 10),
      plot.title   = element_text(hjust = 0.5, size = 11),
      plot.margin  = margin(4, 8, 0, 4)
    )
  
  # (b) gene hit ticks
  p_ticks <- ggplot(pd$ticks) +
    geom_segment(aes(x = rank, xend = rank, y = 0, yend = 1), linewidth = 0.2) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, N)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(0, 8, 0, 4))
  
  # (c) ranked-metric gradient strip (+ high/low labels on bottom row only)
  strip <- pd$stats %>%
    mutate(bin = cut(rank, breaks = n_bins, labels = FALSE)) %>%
    group_by(bin) %>%
    summarise(rank = mean(rank), stat = mean(stat), .groups = "drop")
  
  p_bar <- ggplot(strip, aes(rank, 0, fill = stat)) +
    geom_tile() +
    scale_fill_gradient2(low = "#3182bd", mid = "white", high = "#cb181d",
                         midpoint = 0, guide = "none") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, N)) +
    scale_y_continuous(limits = c(-1.8, 0.5), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() + theme(plot.margin = margin(0, 8, 2, 4))
  
  if (show_xlabels) {
    p_bar <- p_bar +
      annotate("text", x = 0.01 * N, y = -1.1, label = high_label,
               hjust = 0, fontface = "italic", colour = "#cb181d", size = 3.2) +
      annotate("text", x = 0.99 * N, y = -1.1, label = low_label,
               hjust = 1, fontface = "italic", colour = "black", size = 3.2)
  }
  
  p_es / p_ticks / p_bar + plot_layout(heights = c(1, 0.13, 0.14))
}

# ---- 5. Build the four panels 
s_jaat <- get_stat(fg_cgp, "JAATINEN_HEMATOPOIETIC_STEM_CELL_UP")
s_g2m  <- get_stat(fg_h,   "HALLMARK_G2M_CHECKPOINT")
s_e2f  <- get_stat(fg_h,   "HALLMARK_E2F_TARGETS")

panel_jaat <- make_panel(pathways_cgp[["JAATINEN_HEMATOPOIETIC_STEM_CELL_UP"]],
                         ranks_vec, s_jaat$NES, s_jaat$padj,
                         "Hematopoietic stem cell up (Jaatinen)", show_xlabels = FALSE)

panel_g2m  <- make_panel(pathways_h[["HALLMARK_G2M_CHECKPOINT"]],
                         ranks_vec, s_g2m$NES, s_g2m$padj,
                         "G2M Checkpoint (Hallmark)", show_xlabels = TRUE)

panel_e2f  <- make_panel(pathways_h[["HALLMARK_E2F_TARGETS"]],
                         ranks_vec, s_e2f$NES, s_e2f$padj,
                         "E2F Targets (Hallmark)", show_xlabels = TRUE)

s_lsc <- get_stat(fg_cgp, lsc_id)
panel_lsc <- make_panel(pathways_cgp[[lsc_id]], ranks_vec, s_lsc$NES, s_lsc$padj,
                        "AML LSC47 (Huang)", show_xlabels = FALSE)

# ---- 6. Assemble 2x2 (wrap_elements avoids patchwork nesting error) 
Figure_6E <- (wrap_elements(panel_jaat) | wrap_elements(panel_lsc)) /
  (wrap_elements(panel_g2m)  | wrap_elements(panel_e2f))

print(Figure_6E)
# ---- 7. Source-data table 
stats_tbl <- bind_rows(
  data.frame(panel = "Hematopoietic stem cell up (Jaatinen)",
             set = "JAATINEN_HEMATOPOIETIC_STEM_CELL_UP", NES = s_jaat$NES, padj = s_jaat$padj),
  data.frame(panel = "AML LSC47 (Huang)", set = lsc_id,
             NES = s_lsc$NES, padj = s_lsc$padj),
  data.frame(panel = "G2M Checkpoint (Hallmark)",
             set = "HALLMARK_G2M_CHECKPOINT", NES = s_g2m$NES, padj = s_g2m$padj),
  data.frame(panel = "E2F Targets (Hallmark)",
             set = "HALLMARK_E2F_TARGETS", NES = s_e2f$NES, padj = s_e2f$padj)
)
fwrite(stats_tbl, "figureE_gsea_stats.csv")


print(stats_tbl)
cat("Done. Saved figureE_gsea_arid1b_4panels.pdf\n")

############################## Figure 6F ############################### ----
# Figure F - Volcano: ARID1B-high vs ARID1B-low DE (De novo KMT2A-r AML patients)
#
# Input:  rna_arid1b_low_high_mllr.tsv  (same file the GSEA scripts rank from)
#         columns: Gene, "Log2 Ratio" (= LOW - HIGH), "p.Value"
# Output: figureF_volcano_arid1b_high_vs_low.pdf / .png
#         figureF_volcano_source_data.csv

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# setwd("/path/to/folder/with/the/input/file")   # <- point this at your data
infile <- "rna_arid1b_low_high_mllr.tsv"

# ---- Thresholds & labels (edit here) 
pval_thr      <- 0.05     # significance cutoff -> dashed horizontal line
lfc_thr       <- 0        # log2FC cutoff       -> dashed vertical line
label_p_cut   <- 1e-5     # auto-label a gene if p < this ...
label_lfc_cut <- 0.5      # ... or |log2FC| > this
pval_floor    <- 1e-300   # floor p-values for plotting stability

# Genes always labelled AND ringed (curated)
forced_genes <- c(
  "CD34","IRAK1BP1","MECOM","PROM1","ERG","HOXA11",
  "CD70","TLR10","CADM1","CCR10","ADAM23","IKBKG","ARID1A","SMARCA4"
)

col_high <- "#D62728"   # up in ARID1B-HIGH (red)
col_low  <- "#1F77B4"   # up in ARID1B-LOW  (blue)
col_ns   <- "grey60"    # not significant
lab_cols <- c("Up in ARID1B-HIGH" = col_high,
              "Up in ARID1B-LOW"  = col_low,
              "NS"                = "grey40")

# ---- 1. DE table (flip LOW-HIGH -> HIGH-LOW) 
df <- fread(infile, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
names(df) <- make.names(names(df), unique = TRUE)

if (!"Log2.Ratio" %in% names(df))
  stop("No 'Log2.Ratio' column after make.names(). Columns: ", paste(names(df), collapse = ", "))
if (!"p.Value" %in% names(df))
  stop("No 'p.Value' column after make.names(). Columns: ", paste(names(df), collapse = ", "))

volc <- df %>%
  transmute(Gene   = Gene,
            log2FC = -`Log2.Ratio`,    # source is LOW-HIGH; negate -> HIGH-LOW
            pval   = `p.Value`) %>%
  filter(!is.na(Gene), !is.na(log2FC), !is.na(pval)) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  mutate(
    pval_plot = pmax(pval, pval_floor),
    negLog10P = -log10(pval_plot),
    sig = case_when(
      pval < pval_thr & log2FC >  lfc_thr ~ "Up in ARID1B-HIGH",
      pval < pval_thr & log2FC < -lfc_thr ~ "Up in ARID1B-LOW",
      TRUE                                ~ "NS"
    ),
    Gene_key = toupper(trimws(Gene))
  )

# ---- 2. Label sets: forced (always) + auto (droppable if crowded) 
wanted        <- toupper(forced_genes)
forced_labels <- volc %>% filter(Gene_key %in% wanted)
auto_labels   <- volc %>%
  filter(pval < label_p_cut | abs(log2FC) > label_lfc_cut) %>%
  filter(!(Gene_key %in% wanted)) %>%
  distinct(Gene_key, .keep_all = TRUE)

missing_forced <- setdiff(wanted, unique(forced_labels$Gene_key))
if (length(missing_forced))
  message("Forced genes not found in data: ", paste(missing_forced, collapse = ", "))

# ---- 3. Base volcano 
p_base <- ggplot(volc, aes(log2FC, negLog10P)) +
  geom_point(data = subset(volc, sig == "NS"),
             size = 1.2, alpha = 0.55, colour = col_ns) +
  geom_point(data = subset(volc, sig == "Up in ARID1B-HIGH"),
             size = 1.6, alpha = 0.90, colour = col_high) +
  geom_point(data = subset(volc, sig == "Up in ARID1B-LOW"),
             size = 1.6, alpha = 0.90, colour = col_low) +
  geom_vline(xintercept = lfc_thr,            linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(pval_thr),   linetype = "dashed", linewidth = 0.4) +
  labs(
    title = "De novo KMT2A-r AML patients",
    x = expression("Gene expression log2FC ("*italic("ARID1B")*"-high vs "*italic("ARID1B")*"-low)"),
    y = expression(-log[10]("p-value"))
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "italic", size = 12))

# ---- 4. Rings on forced genes + two-layer labels (= Figure F) 
p_volcano_final <- p_base +
  geom_point(data = forced_labels, aes(log2FC, negLog10P),
             shape = 21, fill = NA, colour = "black",
             size = 3.8, stroke = 1.1, inherit.aes = FALSE) +
  ggrepel::geom_text_repel(
    data = forced_labels,
    aes(log2FC, negLog10P, label = Gene, colour = sig),
    size = 4.2, fontface = "italic",
    box.padding = 0.25, point.padding = 0.4,
    min.segment.length = 0.6, segment.alpha = 0.8,
    max.overlaps = Inf, force = 4, seed = 123, show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = auto_labels,
    aes(log2FC, negLog10P, label = Gene, colour = sig),
    size = 4.0, fontface = "italic",
    box.padding = 0.25, point.padding = 0.2,
    min.segment.length = 0, segment.alpha = 0.6,
    force = 2, seed = 123, show.legend = FALSE
  ) +
  scale_colour_manual(values = lab_cols) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 20, 5.5, 5.5))   # right margin for edge labels

print(p_volcano_final)
cat("Done. Saved figureF_volcano_arid1b_high_vs_low.pdf\n")

