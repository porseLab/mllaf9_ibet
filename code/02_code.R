##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load custom functions
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
matrix2Heatmap <- function(x, y=NULL, scale=NULL, col=NULL, bias=NULL, clusterRows=FALSE, clusterCols=FALSE, displayN=FALSE, return_rows=FALSE, title=NULL,
                           columnLabels=NULL, rowAnnot=NULL, legendPos="bottom", showRowNames=T, showColNames=T, displayN_customMat=NULL) {
  library(ComplexHeatmap)
  library(ggpubr)
  library(tiff)
  
  # x = matrix or data frame
  # y = cluster information (from cluster_matrix)
  if(is.null(y)) {
    y <- rep(1, nrow(x));
  }
  
  if(is.null(scale)) { scale="none"; }
  else if(scale=="col") {
    scale="column"
  }
  
  if(is.null(bias)) { bias=1; }
  
  if(scale=="row") {
    x <- t(scale(t(x)))
  } else if(scale=="col" | scale=="column") {
    x <- scale(x)
  }
  if(is.null(col)) { col=colorRampPalette(c("dodgerblue4","grey97", "sienna3"), bias=bias)(250); } 
  else if(length(col) > 1) { col=col; }
  else { col=(colorBrewer2palette(name = col, count = 250, bias=bias)); }
  
  if(is.null(title)) { title=NA; }
  
  if(is.null(columnLabels)) { columnLabels = colnames(x); }
  
  mat <- as.matrix(x)
  df <- as.data.frame.matrix(mat)
  if(length(y[grepl("^\\d+$", y)])==nrow(df)) {
    df$clusters <- as.numeric(y)
    df <- df[order(df$clusters),]
  } else {
    df$clusters <- factor(y, levels=unique(y), ordered=T)
  }
  
  mat <- as.matrix(df[,c(1:ncol(df)-1)])

  if(return_rows==TRUE) {
    return(row_order(draw(Heatmap(mat))))
  } else if(displayN==FALSE) {
    grid.grabExpr(draw(Heatmap(mat, cluster_rows=clusterRows, cluster_columns=clusterCols, col=col, use_raster = T, column_title = title,
                               show_row_names=showRowNames, show_column_names=showColNames, raster_device="png", raster_quality = 2, split = df$clusters, gap = unit(2, "mm"),
                               border=F, column_labels=columnLabels, right_annotation = rowAnnot, 
                               heatmap_legend_param = list(direction = "horizontal"), row_title_rot = 0),
                       heatmap_legend_side = legendPos, annotation_legend_side = legendPos))
  } else {
    if(is.null(displayN_customMat)) { displayN_customMat <- x }
    grid.grabExpr(draw(Heatmap(mat, cluster_rows=clusterRows, cluster_columns=clusterCols, col=col, use_raster = T, column_title = title,
                               show_row_names=showRowNames, show_column_names=showColNames, raster_device="png", raster_quality = 2, split = df$clusters, gap = unit(2, "mm"),
                               cell_fun = function(j, i, x, y, width, height, fill) { grid.text(sprintf("%s", displayN_customMat[i, j]), x, y, 
                                                                                                gp = gpar(fontsize = 8)) },
                               rect_gp = gpar(col = "black"), column_labels=columnLabels, right_annotation = rowAnnot, heatmap_legend_param = list(direction = "horizontal")),
                       heatmap_legend_side = legendPos, annotation_legend_side = legendPos))
  }
}

scale_df <- function(x, y) {
  ## de (default center normalized)
  ## zo (Zero-One); 
  ## zs (Z-Score);
  ## qu (quantile normalization)
  ## rank (rank normalization)
  ## cp (contribution percentage)
  if(missing(y)) { y="de"; }
  if(y=="de") {
    res <- apply(x, 2, function(x) scale(x))
  }
  else if(y=="zo") {
    maxs <- apply(x, 2, max)
    mins <- apply(x, 2, min)
    res <- scale(x, center = mins, scale = maxs - mins)
  } else if(y=="zs") {
    maxs <- apply(x, 2, max)
    mins <- apply(x, 2, min)
    res <- scale(x, center=(maxs+mins)/2, scale=(maxs-mins)/2)
  } else if(y=="qu") {
    res <- normalize.quantiles(as.matrix(x))
  } else if(y=="rank") {
    set.seed(1)
    res <- apply(x, 2, function(x) rank(x, ties.method = "last"))
  } else if (y=="cp") {
    res <- sweep(x, 2, colSums(x), FUN = "/") * 100
  } else {
    res <- x
  }
  
  if(!is.null(colnames(x))) {
    colnames(res) <- colnames(x)
  }
  return(res)
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load data
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(preprocessCore)
library(reshape2)
library(plyr)
library(dplyr)

df_atac <- read.table("DF_ATAC.BED", header = T)
df_promoter <- read.table("DF_PROMOTER.BED", header = T)
df_enhancer <- read.table("DF_ENHANCER.BED", header = T)
df_enhancer_human <- read.table("DF_ENHANCER_HUMAN.BED", header = T)
df_ep_atac <- read.table("DF_ENHANCER.BED.hematopoiesis", header=T)
colnames(df_ep_atac) <- gsub("GSE100738_", "", colnames(df_ep_atac))
df_cebpAE_rpkm <- read.table("DF_CEBPAE_RPKM.BED", header = T)

pop_interest <- c("LTHSC.34..BM", "LTHSC.34..BM.1", "STHSC.150..BM", "MPP4.135..BM", "MPP3.48..BM",
                  "proB.CLP.BM", "preT.DN1.Th", "NK.27.11b..BM", "NK.27.11b..Sp",
                  "ILC2.SI", "ILC3.NKp46.CCR6..SI", "ILC3.NKp46..SI", "ILC3.CCR6..SI",
                  "Mo.6C.II..Bl", "Mo.6C.II..Bl.1", "GN.BM",
                  "DC.4..Sp", "DC.8..Sp", "DC.pDC.Sp")

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## Figure 1F
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
ggscatter(df_promoter, x = "mean_hmll1_dox", y = "logFC_h3k79me2", 
          color = "class_h3k79me2",
          palette = get_palette(c("#3182bd", "#d9d9d9", "#de2d26"), 3),
          alpha=0.8) + 
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 8, lty=2) +
  ylab("Changes in H3K79me2 (logFC; dox vs veh)") +
  xlab("Mean MLL-AF9 signal at promoter (dox)") +
  theme(legend.title= element_blank(), legend.position = "none")

ggscatter(df_promoter, x = "mean_hmll1_dox", y = "logFC_expr", color = "interest_group",
          palette = get_palette(c("#54278f", "#de2d26", "#31a354", "#3182bd", "#d9d9d9"), 5),
          alpha=0.8,
          add = "reg.line", conf.int = TRUE) + 
  stat_cor(aes(color = interest_group), label.x = 1.5, method="pearson") +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 8, lty=2) +
  ylab("Changes in gene expr. (logFC; dox vs veh)") +
  xlab("Mean MLL-AF9 signal at promoter (dox)") +
  theme(legend.title= element_blank())

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## Figure 1G
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## select expression of 39 primary target genes
mat <- df_cebpAE_rpkm %>% dplyr::filter(name %in% df_promoter[which(df_promoter$interest_group=="MLL_target_interest"),]$name) %>% 
  dplyr::select("lsk", "pregm", "gmp", "epm", "lpm") %>%
  'row.names<-'(df_promoter[which(df_promoter$interest_group=="MLL_target_interest"),]$name)
## order then by fold change in gene expression between LSK and GMP
mat <- mat[df_promoter[which(df_promoter$interest_group=="MLL_target_interest"),] %>% dplyr::arrange(-logFC_lsk_gmp) %>% dplyr::pull(name),]
## set Six1 expression to NA as it is not expressed
mat[which(rowMax(as.matrix(mat[,c(1:5)]))<=0.1),c(1:5)] <- NA
## plot heatmap
ggarrange(matrix2Heatmap(mat, clusterRows = F, clusterCols = F, scale="row"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## Figure 2A
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
ggscatter(df_enhancer %>% dplyr::mutate(FDR_atac = -log10(FDR_atac)), x="logFC_atac", y="FDR_atac", color="class_atac",
          palette = c("#3182bd", "#d9d9d9", "#de2d26"), xlab="ATAC-seq logFC for pEs (MA9-on/off)", ylab="-log10(FDR)") +
  geom_hline(yintercept = -log10(0.05), lty=2) + 
  geom_vline(xintercept = c(-1, 1), lty=2)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## Figure 2B
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
df <- df_enhancer
df$mean_atac_dox <- rowMeans(df[,grep("^atac_dox", colnames(df), value=T)])
df$mean_atac_veh <- rowMeans(df[,grep("^atac_veh", colnames(df), value=T)])
df$mean_h3k4me1_dox <- (rowMeans(df[,grep("^h3k4me1_dox_", colnames(df), value=T)]))
df$mean_h3k4me1_veh <- (rowMeans(df[,grep("^h3k4me1_veh_", colnames(df), value=T)]))
df$mean_h3k27ac_dox <- (rowMeans(df[,grep("^h3k27ac_dox_", colnames(df), value=T)]))
df$mean_h3k27ac_veh <- (rowMeans(df[,grep("^h3k27ac_veh_", colnames(df), value=T)]))

df <- merge(df, df_ep_atac[,c("name", pop_interest[c(1:3,5,16)])], by.x="name_dhs_dox", by.y="name")
df[,c("mean_atac_lsk", "mean_atac_gmp", "GN.BM")] <- normalize.quantiles(as.matrix(df[,c("mean_atac_lsk", "mean_atac_gmp", "GN.BM")]))
df <- df[,c(2:4,1,5:ncol(df))]

ra1 <- rowAnnotation(na_col = "black",
                     N1 = anno_barplot(as.data.frame.matrix(table(df[,c("class_activity", "novel_atac_peak")])), gp = gpar(fill = c("#6baed6", "#c6dbef"), col = c("#6baed6", "#c6dbef")), add_numbers=T),
                     atac_dox_veh = data.matrix(dcast(data.table::as.data.table(summaryBy(logFC_atac ~ class_activity + novel_atac_peak, df, fun="mean")),
                                                      class_activity ~ novel_atac_peak, value.var="logFC_atac.mean")[,c(2:3)]),
                     atac_dox_known = dcast(as.data.table(summaryBy(mean_atac_dox + mean_atac_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                            class_activity ~ novel_atac_peak, value.var=c("mean_atac_dox.mean", "mean_atac_veh.mean"))$mean_atac_dox.mean_known,
                     atac_veh_known = dcast(as.data.table(summaryBy(mean_atac_dox + mean_atac_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                            class_activity ~ novel_atac_peak, value.var=c("mean_atac_dox.mean", "mean_atac_veh.mean"))$mean_atac_veh.mean_known,
                     atac_dox_novel = dcast(as.data.table(summaryBy(mean_atac_dox + mean_atac_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                            class_activity ~ novel_atac_peak, value.var=c("mean_atac_dox.mean", "mean_atac_veh.mean"))$mean_atac_dox.mean_novel,
                     atac_veh_novel = dcast(as.data.table(summaryBy(mean_atac_dox + mean_atac_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                            class_activity ~ novel_atac_peak, value.var=c("mean_atac_dox.mean", "mean_atac_veh.mean"))$mean_atac_veh.mean_novel,
                     blank1 = anno_empty(border = FALSE, width=unit(0.4, "cm")),
                     h3k4me1_dox_veh = data.matrix(dcast(data.table::as.data.table(summaryBy(logFC_h3k4me1 ~ class_activity + novel_atac_peak, df, fun="mean")),
                                                         class_activity ~ novel_atac_peak, value.var="logFC_h3k4me1.mean")[,c(2:3)]),
                     h3k4me1_dox_known = dcast(as.data.table(summaryBy(mean_h3k4me1_dox + mean_h3k4me1_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                               class_activity ~ novel_atac_peak, value.var=c("mean_h3k4me1_dox.mean", "mean_h3k4me1_veh.mean"))$mean_h3k4me1_dox.mean_known,
                     h3k4me1_veh_known = dcast(as.data.table(summaryBy(mean_h3k4me1_dox + mean_h3k4me1_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                               class_activity ~ novel_atac_peak, value.var=c("mean_h3k4me1_dox.mean", "mean_h3k4me1_veh.mean"))$mean_h3k4me1_veh.mean_known,
                     h3k4me1_dox_novel = dcast(as.data.table(summaryBy(mean_h3k4me1_dox + mean_h3k4me1_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                               class_activity ~ novel_atac_peak, value.var=c("mean_h3k4me1_dox.mean", "mean_h3k4me1_veh.mean"))$mean_h3k4me1_dox.mean_novel,
                     h3k4me1_veh_novel = dcast(as.data.table(summaryBy(mean_h3k4me1_dox + mean_h3k4me1_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                               class_activity ~ novel_atac_peak, value.var=c("mean_h3k4me1_dox.mean", "mean_h3k4me1_veh.mean"))$mean_h3k4me1_veh.mean_novel,
                     blank3 = anno_empty(border = FALSE, width=unit(0.4, "cm")),
                     h3k27ac_dox_veh = data.matrix(dcast(data.table::as.data.table(summaryBy(logFC_h3k27ac ~ class_activity + novel_atac_peak, df, fun="mean")),
                                                         class_activity ~ novel_atac_peak, value.var="logFC_h3k27ac.mean")[,c(2:3)]),
                     h3k27ac_dox_known = dcast(as.data.table(summaryBy(mean_h3k27ac_dox + mean_h3k27ac_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                               class_activity ~ novel_atac_peak, value.var=c("mean_h3k27ac_dox.mean", "mean_h3k27ac_veh.mean"))$mean_h3k27ac_dox.mean_known,
                     h3k27ac_veh_known = dcast(as.data.table(summaryBy(mean_h3k27ac_dox + mean_h3k27ac_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                               class_activity ~ novel_atac_peak, value.var=c("mean_h3k27ac_dox.mean", "mean_h3k27ac_veh.mean"))$mean_h3k27ac_veh.mean_known,
                     h3k27ac_dox_novel = dcast(as.data.table(summaryBy(mean_h3k27ac_dox + mean_h3k27ac_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                               class_activity ~ novel_atac_peak, value.var=c("mean_h3k27ac_dox.mean", "mean_h3k27ac_veh.mean"))$mean_h3k27ac_dox.mean_novel,
                     h3k27ac_veh_novel = dcast(as.data.table(summaryBy(mean_h3k27ac_dox + mean_h3k27ac_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                               class_activity ~ novel_atac_peak, value.var=c("mean_h3k27ac_dox.mean", "mean_h3k27ac_veh.mean"))$mean_h3k27ac_veh.mean_novel,
                     blank2 = anno_empty(border = FALSE, width=unit(0.4, "cm")),
                     atac_lsk_known = (dcast(data.table::as.data.table(summaryBy(mean_atac_lsk ~ class_activity + novel_atac_peak, df, fun="mean")),
                                             class_activity ~ novel_atac_peak, value.var="mean_atac_lsk.mean")[,c(2)])$known,
                     atac_gmp_known = (dcast(data.table::as.data.table(summaryBy(mean_atac_gmp ~ class_activity + novel_atac_peak, df, fun="mean")),
                                             class_activity ~ novel_atac_peak, value.var="mean_atac_gmp.mean")[,c(2)])$known,
                     atac_gran_known = (dcast(data.table::as.data.table(summaryBy(GN.BM ~ class_activity + novel_atac_peak, df, fun="mean")),
                                              class_activity ~ novel_atac_peak, value.var="GN.BM.mean")[,c(2)])$known,
                     atac_lsk_novel = (dcast(data.table::as.data.table(summaryBy(mean_atac_lsk ~ class_activity + novel_atac_peak, df, fun="mean")),
                                             class_activity ~ novel_atac_peak, value.var="mean_atac_lsk.mean")[,c(3)])$novel,
                     atac_gmp_novel = (dcast(data.table::as.data.table(summaryBy(mean_atac_gmp ~ class_activity + novel_atac_peak, df, fun="mean")),
                                             class_activity ~ novel_atac_peak, value.var="mean_atac_gmp.mean")[,c(3)])$novel,
                     atac_gran_novel = (dcast(data.table::as.data.table(summaryBy(GN.BM ~ class_activity + novel_atac_peak, df, fun="mean")),
                                              class_activity ~ novel_atac_peak, value.var="GN.BM.mean")[,c(3)])$novel,
                     col=list(
                       atac_dox_veh = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
                       atac_dox_known = colorRamp2(c(50, 400), c("white", "#737373")),
                       atac_veh_known = colorRamp2(c(50, 400), c("white", "#737373")),
                       atac_dox_novel = colorRamp2(c(50, 400), c("white", "#737373")),
                       atac_veh_novel = colorRamp2(c(50, 400), c("white", "#737373")),
                       h3k4me1_dox_veh = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                       h3k4me1_dox_known = colorRamp2(c(35, 65), c("white", "#737373")),
                       h3k4me1_veh_known = colorRamp2(c(35, 65), c("white", "#737373")),
                       h3k4me1_dox_novel = colorRamp2(c(35, 65), c("white", "#737373")),
                       h3k4me1_veh_novel = colorRamp2(c(35, 65), c("white", "#737373")),
                       h3k27ac_dox_veh = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                       h3k27ac_dox_known = colorRamp2(c(15, 35), c("white", "#737373")),
                       h3k27ac_veh_known = colorRamp2(c(15, 35), c("white", "#737373")),
                       h3k27ac_dox_novel = colorRamp2(c(15, 35), c("white", "#737373")),
                       h3k27ac_veh_novel = colorRamp2(c(15, 35), c("white", "#737373")),
                       atac_lsk_known = colorRamp2(c(0, 14), c("white", "#737373")),
                       atac_gmp_known = colorRamp2(c(0, 14), c("white", "#737373")),
                       atac_gran_known = colorRamp2(c(0, 14), c("white", "#737373")),
                       atac_lsk_novel = colorRamp2(c(0, 14), c("white", "#737373")),
                       atac_gmp_novel = colorRamp2(c(0, 14), c("white", "#737373")),
                       atac_gran_novel = colorRamp2(c(0, 14), c("white", "#737373"))
                     ),
                     annotation_legend_param = list(
                       atac_dox_veh = list(direction = "horizontal"),
                       atac_dox_known = list(direction = "horizontal"),
                       atac_veh_known = list(direction = "horizontal"),
                       atac_dox_novel = list(direction = "horizontal"),
                       atac_veh_novel = list(direction = "horizontal"),
                       h3k4me1_dox_veh = list(direction = "horizontal"),
                       h3k4me1_dox_known = list(direction = "horizontal"),
                       h3k4me1_veh_known = list(direction = "horizontal"),
                       h3k4me1_dox_novel = list(direction = "horizontal"),
                       h3k4me1_veh_novel = list(direction = "horizontal"),
                       h3k27ac_dox_veh = list(direction = "horizontal"),
                       h3k27ac_dox_known = list(direction = "horizontal"),
                       h3k27ac_veh_known = list(direction = "horizontal"),
                       h3k27ac_dox_novel = list(direction = "horizontal"),
                       h3k27ac_veh_novel = list(direction = "horizontal"),
                       atac_lsk_known = list(direction = "horizontal"),
                       atac_gmp_known = list(direction = "horizontal"),
                       atac_gran_known = list(direction = "horizontal"),
                       atac_lsk_novel = list(direction = "horizontal"),
                       atac_gmp_novel = list(direction = "horizontal"),
                       atac_gran_novel = list(direction = "horizontal")
                     )
)

mat <- as.data.frame.matrix(table(df[,c("class_activity", "novel_atac_peak")]))

ggarrange(matrix2Heatmap(mat, scale="none", clusterRows = FALSE, clusterCols = FALSE, title = "Enhancers (mouse)", bias = 3, rowAnnot = c(ra1), displayN = T), ncol=2)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## Figure 4C
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
df <- df_enhancer
df$mean_meis1_dox <- rowMeans(df[,grep("mll_af9_on_meis1", colnames(df), value=T)])
df$mean_meis1_veh <- rowMeans(df[,grep("mll_af9_off_meis1", colnames(df), value=T)])

ra1 <- rowAnnotation(na_col = "black",
                     hoxa9_known = dcast(data.table::as.data.table(summaryBy(hoxa9_signal_rep1 ~ class_activity + novel_atac_peak, df, fun="mean")),
                                         class_activity ~ novel_atac_peak, value.var="hoxa9_signal_rep1.mean")$known,
                     hoxa9_novel = dcast(data.table::as.data.table(summaryBy(hoxa9_signal_rep1 ~ class_activity + novel_atac_peak, df, fun="mean")),
                                         class_activity ~ novel_atac_peak, value.var="hoxa9_signal_rep1.mean")$novel,
                     blank1 = anno_empty(border = FALSE, width=unit(0.4, "cm")),
                     meis1_dox_known = dcast(as.data.table(summaryBy(mean_meis1_dox + mean_meis1_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                             class_activity ~ novel_atac_peak, value.var=c("mean_meis1_dox.mean", "mean_meis1_veh.mean"))$mean_meis1_dox.mean_known,
                     meis1_veh_known = dcast(as.data.table(summaryBy(mean_meis1_dox + mean_meis1_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                             class_activity ~ novel_atac_peak, value.var=c("mean_meis1_dox.mean", "mean_meis1_veh.mean"))$mean_meis1_veh.mean_known,
                     meis1_dox_novel = dcast(as.data.table(summaryBy(mean_meis1_dox + mean_meis1_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                             class_activity ~ novel_atac_peak, value.var=c("mean_meis1_dox.mean", "mean_meis1_veh.mean"))$mean_meis1_dox.mean_novel,
                     meis1_veh_novel = dcast(as.data.table(summaryBy(mean_meis1_dox + mean_meis1_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                             class_activity ~ novel_atac_peak, value.var=c("mean_meis1_dox.mean", "mean_meis1_veh.mean"))$mean_meis1_veh.mean_novel,
                     col=list(
                       hoxa9_known = colorRamp2(c(0, 15), c("white", "#737373")),
                       hoxa9_novel = colorRamp2(c(0, 15), c("white", "#737373")),
                       meis1_dox_known = colorRamp2(c(0, 200), c("white", "#737373")),
                       meis1_veh_known = colorRamp2(c(0, 200), c("white", "#737373")),
                       meis1_dox_novel = colorRamp2(c(0, 200), c("white", "#737373")),
                       meis1_veh_novel = colorRamp2(c(0, 200), c("white", "#737373"))
                     ),
                     annotation_legend_param = list(
                       hoxa9_known = list(direction = "horizontal"),
                       hoxa9_novel = list(direction = "horizontal"),
                       meis1_dox_known = list(direction = "horizontal"),
                       meis1_veh_known = list(direction = "horizontal"),
                       meis1_dox_novel = list(direction = "horizontal"),
                       meis1_veh_novel = list(direction = "horizontal")
                     )
)

mat <- as.data.frame.matrix(table(df[,c("class_activity", "novel_atac_peak")]))

ggarrange(matrix2Heatmap(mat, scale="none", clusterRows = FALSE, clusterCols = FALSE, title = "Enhancers (mouse)", bias = 3, rowAnnot = c(ra1), displayN = T), ncol=2)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## Figure 4H
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
df <- df_enhancer
df$mean_arid1b_dox <- rowMeans(df[,grep("^arid1b_chip_imf9_dox", colnames(df), value=T)])
df$mean_arid1b_veh <- rowMeans(df[,grep("^arid1b_chip_imf9_veh", colnames(df), value=T)])

ra1 <- rowAnnotation(na_col = "black",
                     arid1b_dox_known = dcast(as.data.table(summaryBy(mean_arid1b_dox + mean_arid1b_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                              class_activity ~ novel_atac_peak, value.var=c("mean_arid1b_dox.mean", "mean_arid1b_veh.mean"))$mean_arid1b_dox.mean_known,
                     arid1b_veh_known = dcast(as.data.table(summaryBy(mean_arid1b_dox + mean_arid1b_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                              class_activity ~ novel_atac_peak, value.var=c("mean_arid1b_dox.mean", "mean_arid1b_veh.mean"))$mean_arid1b_veh.mean_known,
                     arid1b_dox_novel = dcast(as.data.table(summaryBy(mean_arid1b_dox + mean_arid1b_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                              class_activity ~ novel_atac_peak, value.var=c("mean_arid1b_dox.mean", "mean_arid1b_veh.mean"))$mean_arid1b_dox.mean_novel,
                     arid1b_veh_novel = dcast(as.data.table(summaryBy(mean_arid1b_dox + mean_arid1b_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                              class_activity ~ novel_atac_peak, value.var=c("mean_arid1b_dox.mean", "mean_arid1b_veh.mean"))$mean_arid1b_veh.mean_novel,
                     col=list(
                       arid1b_dox_known = colorRamp2(c(0, 70), c("white", "#737373")),
                       arid1b_veh_known = colorRamp2(c(0, 70), c("white", "#737373")),
                       arid1b_dox_novel = colorRamp2(c(0, 70), c("white", "#737373")),
                       arid1b_veh_novel = colorRamp2(c(0, 70), c("white", "#737373"))
                     ),
                     annotation_legend_param = list(
                       arid1b_dox_known = list(direction = "horizontal"),
                       arid1b_veh_known = list(direction = "horizontal"),
                       arid1b_dox_novel = list(direction = "horizontal"),
                       arid1b_veh_novel = list(direction = "horizontal")
                     )
)

mat <- as.data.frame.matrix(table(df[,c("class_activity", "novel_atac_peak")]))

ggarrange(matrix2Heatmap(mat, scale="none", clusterRows = FALSE, clusterCols = FALSE, title = "Enhancers (mouse)", bias = 3, rowAnnot = c(ra1), displayN = T), ncol=2)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## Figure 5C
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
df <- df_enhancer
df$mean_arid1a_dox <- rowMeans(df[,grep("^arid1a_chip_imf9_dox", colnames(df), value=T)])
df$mean_arid1a_veh <- rowMeans(df[,grep("^arid1a_chip_imf9_veh", colnames(df), value=T)])
df[,c("mean_atac_arid1bWT", "mean_atac_arid1bKD")] <- normalize.quantiles(as.matrix(df[,c("mean_atac_arid1bWT", "mean_atac_arid1bKD")]))

ra1 <- rowAnnotation(na_col = "black",
                     atac_arid1bWT_known = dcast(as.data.table(summaryBy(mean_atac_arid1bWT + mean_atac_arid1bKD ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                 class_activity ~ novel_atac_peak, value.var=c("mean_atac_arid1bWT.mean", "mean_atac_arid1bKD.mean"))$mean_atac_arid1bWT.mean_known,
                     atac_arid1bKD_known = dcast(as.data.table(summaryBy(mean_atac_arid1bWT + mean_atac_arid1bKD ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                 class_activity ~ novel_atac_peak, value.var=c("mean_atac_arid1bWT.mean", "mean_atac_arid1bKD.mean"))$mean_atac_arid1bKD.mean_known,
                     atac_arid1bWT_novel = dcast(as.data.table(summaryBy(mean_atac_arid1bWT + mean_atac_arid1bKD ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                 class_activity ~ novel_atac_peak, value.var=c("mean_atac_arid1bWT.mean", "mean_atac_arid1bKD.mean"))$mean_atac_arid1bWT.mean_novel,
                     atac_arid1bKD_novel = dcast(as.data.table(summaryBy(mean_atac_arid1bWT + mean_atac_arid1bKD ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                 class_activity ~ novel_atac_peak, value.var=c("mean_atac_arid1bWT.mean", "mean_atac_arid1bKD.mean"))$mean_atac_arid1bKD.mean_novel,
                     logFC_atac_arid1bKD = data.matrix(dcast(data.table::as.data.table(summaryBy(logFC_atac_arid1bKD ~ class_activity + novel_atac_peak,
                                                                                                 df, fun="mean")),
                                                             class_activity ~ novel_atac_peak, value.var="logFC_atac_arid1bKD.mean")[,c(2:3)]),
                     blank1 = anno_empty(border = FALSE, width=unit(0.4, "cm")),
                     arid1a_dox_known = dcast(as.data.table(summaryBy(mean_arid1a_dox + mean_arid1a_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                              class_activity ~ novel_atac_peak, value.var=c("mean_arid1a_dox.mean", "mean_arid1a_veh.mean"))$mean_arid1a_dox.mean_known,
                     arid1a_veh_known = dcast(as.data.table(summaryBy(mean_arid1a_dox + mean_arid1a_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                              class_activity ~ novel_atac_peak, value.var=c("mean_arid1a_dox.mean", "mean_arid1a_veh.mean"))$mean_arid1a_veh.mean_known,
                     arid1a_dox_novel = dcast(as.data.table(summaryBy(mean_arid1a_dox + mean_arid1a_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                              class_activity ~ novel_atac_peak, value.var=c("mean_arid1a_dox.mean", "mean_arid1a_veh.mean"))$mean_arid1a_dox.mean_novel,
                     arid1a_veh_novel = dcast(as.data.table(summaryBy(mean_arid1a_dox + mean_arid1a_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                              class_activity ~ novel_atac_peak, value.var=c("mean_arid1a_dox.mean", "mean_arid1a_veh.mean"))$mean_arid1a_veh.mean_novel,
                     blank2 = anno_empty(border = FALSE, width=unit(0.4, "cm")),
                     arid1a_arid1bWT_known = dcast(as.data.table(summaryBy(mean_arid1a_arid1bWT + mean_arid1a_arid1bKD ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                   class_activity ~ novel_atac_peak, value.var=c("mean_arid1a_arid1bWT.mean", "mean_arid1a_arid1bKD.mean"))$mean_arid1a_arid1bWT.mean_known,
                     arid1a_arid1bKD_known = dcast(as.data.table(summaryBy(mean_arid1a_arid1bWT + mean_arid1a_arid1bKD ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                   class_activity ~ novel_atac_peak, value.var=c("mean_arid1a_arid1bWT.mean", "mean_arid1a_arid1bKD.mean"))$mean_arid1a_arid1bKD.mean_known,
                     arid1a_arid1bWT_novel = dcast(as.data.table(summaryBy(mean_arid1a_arid1bWT + mean_arid1a_arid1bKD ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                   class_activity ~ novel_atac_peak, value.var=c("mean_arid1a_arid1bWT.mean", "mean_arid1a_arid1bKD.mean"))$mean_arid1a_arid1bWT.mean_novel,
                     arid1a_arid1bKD_novel = dcast(as.data.table(summaryBy(mean_arid1a_arid1bWT + mean_arid1a_arid1bKD ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                   class_activity ~ novel_atac_peak, value.var=c("mean_arid1a_arid1bWT.mean", "mean_arid1a_arid1bKD.mean"))$mean_arid1a_arid1bKD.mean_novel,
                     blank3 = anno_empty(border = FALSE, width=unit(0.4, "cm")),
                     brg1_dox_known = dcast(as.data.table(summaryBy(mean_brg1_dox + mean_brg1_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                            class_activity ~ novel_atac_peak, value.var=c("mean_brg1_dox.mean", "mean_brg1_veh.mean"))$mean_brg1_dox.mean_known,
                     brg1_veh_known = dcast(as.data.table(summaryBy(mean_brg1_dox + mean_brg1_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                            class_activity ~ novel_atac_peak, value.var=c("mean_brg1_dox.mean", "mean_brg1_veh.mean"))$mean_brg1_veh.mean_known,
                     brg1_dox_novel = dcast(as.data.table(summaryBy(mean_brg1_dox + mean_brg1_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                            class_activity ~ novel_atac_peak, value.var=c("mean_brg1_dox.mean", "mean_brg1_veh.mean"))$mean_brg1_dox.mean_novel,
                     brg1_veh_novel = dcast(as.data.table(summaryBy(mean_brg1_dox + mean_brg1_veh ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                            class_activity ~ novel_atac_peak, value.var=c("mean_brg1_dox.mean", "mean_brg1_veh.mean"))$mean_brg1_veh.mean_novel,
                     blank4 = anno_empty(border = FALSE, width=unit(0.4, "cm")),
                     logFC_atac_brm014 = data.matrix(dcast(data.table::as.data.table(summaryBy(logFC_atac_brm014 ~ class_activity + novel_atac_peak,
                                                                                               df, fun="mean")),
                                                           class_activity ~ novel_atac_peak, value.var="logFC_atac_brm014.mean")[,c(2:3)]),
                     col=list(
                       atac_arid1bWT_known = colorRamp2(c(0, 50), c("white", "#737373")),
                       atac_arid1bKD_known = colorRamp2(c(0, 50), c("white", "#737373")),
                       atac_arid1bWT_novel = colorRamp2(c(0, 50), c("white", "#737373")),
                       atac_arid1bKD_novel = colorRamp2(c(0, 50), c("white", "#737373")),
                       logFC_atac_arid1bKD = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")),
                       arid1a_dox_known = colorRamp2(c(0, 35), c("white", "#737373")),
                       arid1a_veh_known = colorRamp2(c(0, 35), c("white", "#737373")),
                       arid1a_dox_novel = colorRamp2(c(0, 35), c("white", "#737373")),
                       arid1a_veh_novel = colorRamp2(c(0, 35), c("white", "#737373")),
                       arid1a_arid1bWT_known = colorRamp2(c(0, 45), c("white", "#737373")),
                       arid1a_arid1bKD_known = colorRamp2(c(0, 45), c("white", "#737373")),
                       arid1a_arid1bWT_novel = colorRamp2(c(0, 45), c("white", "#737373")),
                       arid1a_arid1bKD_novel = colorRamp2(c(0, 45), c("white", "#737373")),
                       brg1_dox_known = colorRamp2(c(0, 35), c("white", "#737373")),
                       brg1_veh_known = colorRamp2(c(0, 35), c("white", "#737373")),
                       brg1_dox_novel = colorRamp2(c(0, 35), c("white", "#737373")),
                       brg1_veh_novel = colorRamp2(c(0, 35), c("white", "#737373")),
                       logFC_atac_brm014 = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
                     ),
                     annotation_legend_param = list(
                       atac_arid1bWT_known = list(direction = "horizontal"),
                       atac_arid1bKD_known = list(direction = "horizontal"),
                       atac_arid1bWT_novel = list(direction = "horizontal"),
                       atac_arid1bKD_novel = list(direction = "horizontal"),
                       logFC_atac_arid1bKD = list(direction = "horizontal"),
                       arid1a_dox_known = list(direction = "horizontal"),
                       arid1a_veh_known = list(direction = "horizontal"),
                       arid1a_dox_novel = list(direction = "horizontal"),
                       arid1a_veh_novel = list(direction = "horizontal"),
                       arid1a_arid1bWT_known = list(direction = "horizontal"),
                       arid1a_arid1bKD_known = list(direction = "horizontal"),
                       arid1a_arid1bWT_novel = list(direction = "horizontal"),
                       arid1a_arid1bKD_novel = list(direction = "horizontal"),
                       brg1_dox_known = list(direction = "horizontal"),
                       brg1_veh_known = list(direction = "horizontal"),
                       brg1_dox_novel = list(direction = "horizontal"),
                       brg1_veh_novel = list(direction = "horizontal"),
                       logFC_atac_brm014 = list(direction = "horizontal")
                     )
)

mat <- as.data.frame.matrix(table(df[,c("class_activity", "novel_atac_peak")]))

ggarrange(matrix2Heatmap(mat, scale="none", clusterRows = FALSE, clusterCols = FALSE, title = "Enhancers (mouse)", bias = 3, rowAnnot = c(ra1), displayN = T), ncol=2)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## Figure 6A
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
df <- df_enhancer_human
df$mean_atac_hsc <- rowMeans(df[,grep("donor7256_1", colnames(df), value=T)])
df$mean_atac_gmp <- rowMeans(df[,grep("donor7256_5", colnames(df), value=T)])
df$mean_atac_mono <- rowMeans(df[,grep("donor7256_7", colnames(df), value=T)])
df$mean_h3k27ac_molm13 <- rowMeans(df[,grep("^h3k27ac_molm13_rep1|^h3k27ac_molm13_dmso_rep1", colnames(df), value=T)])
df$mean_h3k27ac_epzvtp <- rowMeans(df[,grep("^h3k27ac_molm13_epzvtp_rep1|^h3k27ac_molm13_vtp_rep1", colnames(df), value=T)])
df$mean_meis1_molm13 <- rowMeans(df[,grep("meis1_molm13_rep", colnames(df), value=T)])
df[,c("mean_atac_hsc", "mean_atac_gmp", "mean_atac_mono")] <- normalize.quantiles(as.matrix(df[,c("mean_atac_hsc", "mean_atac_gmp", "mean_atac_mono")]))
df[,grep("mean_h3k27ac", colnames(df), value=T)] <- matrix2norm(df[,grep("mean_h3k27ac", colnames(df), value=T)], method="qu")

ra1 <- rowAnnotation(na_col = "black",
                     N1 = anno_barplot(as.data.frame.matrix(table(df[,c("class_activity", "novel_atac_peak")])), gp = gpar(fill = c("gray", "black"), col = c("gray", "black")), add_numbers=T),
                     logFC_h3k27ac_epzvtp_known = dcast(as.data.table(summaryBy(logFC_h3k27ac_epzvtp ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                        class_activity ~ novel_atac_peak, value.var=c("logFC_h3k27ac_epzvtp.mean"))$known,
                     logFC_h3k27ac_epzvtp_novel = dcast(as.data.table(summaryBy(logFC_h3k27ac_epzvtp ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                        class_activity ~ novel_atac_peak, value.var=c("logFC_h3k27ac_epzvtp.mean"))$novel,
                     h3k27ac_molm13_known = dcast(as.data.table(summaryBy(mean_h3k27ac_molm13 + mean_h3k27ac_epzvtp ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                  class_activity ~ novel_atac_peak, value.var=c("mean_h3k27ac_molm13.mean", "mean_h3k27ac_epzvtp.mean"))$mean_h3k27ac_molm13.mean_known,
                     h3k27ac_epzvtp_known = dcast(as.data.table(summaryBy(mean_h3k27ac_molm13 + mean_h3k27ac_epzvtp ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                  class_activity ~ novel_atac_peak, value.var=c("mean_h3k27ac_molm13.mean", "mean_h3k27ac_epzvtp.mean"))$mean_h3k27ac_epzvtp.mean_known,
                     h3k27ac_molm13_novel = dcast(as.data.table(summaryBy(mean_h3k27ac_molm13 + mean_h3k27ac_epzvtp ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                  class_activity ~ novel_atac_peak, value.var=c("mean_h3k27ac_molm13.mean", "mean_h3k27ac_epzvtp.mean"))$mean_h3k27ac_molm13.mean_novel,
                     h3k27ac_epzvtp_novel = dcast(as.data.table(summaryBy(mean_h3k27ac_molm13 + mean_h3k27ac_epzvtp ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                  class_activity ~ novel_atac_peak, value.var=c("mean_h3k27ac_molm13.mean", "mean_h3k27ac_epzvtp.mean"))$mean_h3k27ac_epzvtp.mean_novel,
                     blank1 = anno_empty(border = FALSE, width=unit(0.4, "cm")),
                     atac_hsc_known = dcast(as.data.table(summaryBy(mean_atac_hsc + mean_atac_gmp + mean_atac_mono ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                            class_activity ~ novel_atac_peak, value.var=c("mean_atac_hsc.mean", "mean_atac_gmp.mean", "mean_atac_mono.mean"))$mean_atac_hsc.mean_known,
                     atac_gmp_known = dcast(as.data.table(summaryBy(mean_atac_hsc + mean_atac_gmp + mean_atac_mono ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                            class_activity ~ novel_atac_peak, value.var=c("mean_atac_hsc.mean", "mean_atac_gmp.mean", "mean_atac_mono.mean"))$mean_atac_gmp.mean_known,
                     atac_mono_known = dcast(as.data.table(summaryBy(mean_atac_hsc + mean_atac_gmp + mean_atac_mono ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                             class_activity ~ novel_atac_peak, value.var=c("mean_atac_hsc.mean", "mean_atac_gmp.mean", "mean_atac_mono.mean"))$mean_atac_mono.mean_known,
                     atac_hsc_novel = dcast(as.data.table(summaryBy(mean_atac_hsc + mean_atac_gmp + mean_atac_mono ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                            class_activity ~ novel_atac_peak, value.var=c("mean_atac_hsc.mean", "mean_atac_gmp.mean", "mean_atac_mono.mean"))$mean_atac_hsc.mean_novel,
                     atac_gmp_novel = dcast(as.data.table(summaryBy(mean_atac_hsc + mean_atac_gmp + mean_atac_mono ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                            class_activity ~ novel_atac_peak, value.var=c("mean_atac_hsc.mean", "mean_atac_gmp.mean", "mean_atac_mono.mean"))$mean_atac_gmp.mean_novel,
                     atac_mono_novel = dcast(as.data.table(summaryBy(mean_atac_hsc + mean_atac_gmp + mean_atac_mono ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                             class_activity ~ novel_atac_peak, value.var=c("mean_atac_hsc.mean", "mean_atac_gmp.mean", "mean_atac_mono.mean"))$mean_atac_mono.mean_novel,
                     blank4 = anno_empty(border = FALSE, width=unit(0.4, "cm")),
                     meis1_molm13_known = dcast(as.data.table(summaryBy(mean_meis1_molm13 ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                class_activity ~ novel_atac_peak, value.var=c("mean_meis1_molm13.mean"))$known,
                     meis1_molm13_novel = dcast(as.data.table(summaryBy(mean_meis1_molm13 ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                class_activity ~ novel_atac_peak, value.var=c("mean_meis1_molm13.mean"))$novel,
                     blank5 = anno_empty(border = FALSE, width=unit(0.4, "cm")),
                     atac_brm014_72h_known = dcast(as.data.table(summaryBy(logFC_brm014_72h ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                   class_activity ~ novel_atac_peak, value.var=c("logFC_brm014_72h.mean"))$known,
                     atac_brm014_72h_novel = dcast(as.data.table(summaryBy(logFC_brm014_72h ~ class_activity + novel_atac_peak, df, FUN = c(mean))),
                                                   class_activity ~ novel_atac_peak, value.var=c("logFC_brm014_72h.mean"))$novel,
                     col=list(
                       logFC_h3k27ac_epzvtp_known = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                       logFC_h3k27ac_epzvtp_novel = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                       h3k27ac_molm13_known = colorRamp2(c(0, 100), c("white", "#737373")),
                       h3k27ac_epzvtp_known = colorRamp2(c(0, 100), c("white", "#737373")),
                       h3k27ac_molm13_novel = colorRamp2(c(0, 100), c("white", "#737373")),
                       h3k27ac_epzvtp_novel = colorRamp2(c(0, 100), c("white", "#737373")),
                       atac_hsc_known = colorRamp2(c(15, 50), c("white", "#737373")),
                       atac_gmp_known = colorRamp2(c(15, 50), c("white", "#737373")),
                       atac_mono_known = colorRamp2(c(15, 50), c("white", "#737373")),
                       atac_hsc_novel = colorRamp2(c(0, 50), c("white", "#737373")),
                       atac_gmp_novel = colorRamp2(c(0, 50), c("white", "#737373")),
                       atac_mono_novel = colorRamp2(c(0, 50), c("white", "#737373")),
                       meis1_molm13_known = colorRamp2(c(80, 200), c("white", "#737373")),
                       meis1_molm13_novel = colorRamp2(c(80, 200), c("white", "#737373")),
                       atac_brm014_72h_known = colorRamp2(c(-0.8, -0.4, 0, 0.4, 0.8), c("blue", "#d1e5f0", "#f7f7f7", "#fddbc7", "red")),
                       atac_brm014_72h_novel = colorRamp2(c(-0.8, -0.4, 0, 0.4, 0.8), c("blue", "#d1e5f0", "#f7f7f7", "#fddbc7", "red"))
                     ), 
                     annotation_legend_param = list(
                       logFC_h3k27ac_epzvtp_known = list(direction = "horizontal"),
                       logFC_h3k27ac_epzvtp_novel = list(direction = "horizontal"),
                       h3k27ac_molm13_known = list(direction = "horizontal"),
                       h3k27ac_epzvtp_known = list(direction = "horizontal"),
                       h3k27ac_molm13_novel = list(direction = "horizontal"),
                       h3k27ac_epzvtp_novel = list(direction = "horizontal"),
                       atac_hsc_known = list(direction = "horizontal"),
                       atac_gmp_known = list(direction = "horizontal"),
                       atac_mono_known = list(direction = "horizontal"),
                       atac_hsc_novel = list(direction = "horizontal"),
                       atac_gmp_novel = list(direction = "horizontal"),
                       atac_mono_novel = list(direction = "horizontal"),
                       meis1_molm13_known = list(direction = "horizontal"),
                       meis1_molm13_novel = list(direction = "horizontal"),
                       atac_brm014_72h_known = list(direction = "horizontal"),
                       atac_brm014_72h_novel = list(direction = "horizontal")
                     )
)

mat <- as.data.frame.matrix(table(df[,c("class_activity", "novel_atac_peak")]))

ggarrange(matrix2Heatmap(mat, scale="none", clusterRows = FALSE, clusterCols = FALSE, title = "Enhancers (human)", bias = 3, rowAnnot = c(ra1), displayN = T), ncol=2)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## Figure 7A
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
df <- df_enhancer
df[,c("mean_h3k27ac_sen", "mean_h3k27ac_res")] <- normalize.quantiles(as.matrix(df[,c("mean_h3k27ac_sen", "mean_h3k27ac_res")]))

ra1 <- rowAnnotation(na_col = "black",
                     h3k27ac_sensitive_known = data.matrix(dcast(data.table::as.data.table(summaryBy(mean_h3k27ac_sen ~ class_activity + novel_atac_peak,
                                                                                                     df, fun="mean")),
                                                                 class_activity ~ novel_atac_peak, value.var="mean_h3k27ac_sen.mean")[,c(2)]),
                     h3k27ac_resistant_known = data.matrix(dcast(data.table::as.data.table(summaryBy(mean_h3k27ac_res ~ class_activity + novel_atac_peak,
                                                                                                     df, fun="mean")),
                                                                 class_activity ~ novel_atac_peak, value.var="mean_h3k27ac_res.mean")[,c(2)]),
                     h3k27ac_sensitive_novel = data.matrix(dcast(data.table::as.data.table(summaryBy(mean_h3k27ac_sen ~ class_activity + novel_atac_peak,
                                                                                                     df, fun="mean")),
                                                                 class_activity ~ novel_atac_peak, value.var="mean_h3k27ac_sen.mean")[,c(3)]),
                     h3k27ac_resistant_novel = data.matrix(dcast(data.table::as.data.table(summaryBy(mean_h3k27ac_res ~ class_activity + novel_atac_peak,
                                                                                                     df, fun="mean")),
                                                                 class_activity ~ novel_atac_peak, value.var="mean_h3k27ac_res.mean")[,c(3)]),
                     logFC_h3k27ac_ibet = data.matrix(dcast(data.table::as.data.table(summaryBy(logFC_h3k27ac_ibet ~ class_activity + novel_atac_peak,
                                                                                                df, fun="mean")),
                                                            class_activity ~ novel_atac_peak, value.var="logFC_h3k27ac_ibet.mean")[,c(2:3)]),
                     col=list(
                       h3k27ac_sensitive_known = colorRamp2(c(0, 40), c("white", "#737373")),
                       h3k27ac_resistant_known = colorRamp2(c(0, 40), c("white", "#737373")),
                       h3k27ac_sensitive_novel = colorRamp2(c(0, 40), c("white", "#737373")),
                       h3k27ac_resistant_novel = colorRamp2(c(0, 40), c("white", "#737373")),
                       logFC_h3k27ac_ibet = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
                     ),
                     annotation_legend_param = list(
                       h3k27ac_sensitive_known = list(direction = "horizontal"),
                       h3k27ac_resistant_known = list(direction = "horizontal"),
                       h3k27ac_sensitive_novel = list(direction = "horizontal"),
                       h3k27ac_resistant_novel = list(direction = "horizontal"),
                       logFC_h3k27ac_ibet = list(direction = "horizontal")
                     )
)

mat <- as.data.frame.matrix(table(df[,c("class_activity", "novel_atac_peak")]))

ggarrange(matrix2Heatmap(mat, scale="none", clusterRows = FALSE, clusterCols = FALSE, title = "Enhancers (mouse)", bias = 3, rowAnnot = c(ra1), displayN = T), ncol=2)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## Supplementary Figure S1E, F and H
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

ggscatter(df_promoter, x = "mean_hmll1_dox", y = "logFC_hmll1", color = "interest_group",
          palette = get_palette(c("#54278f", "#de2d26", "#31a354", "#3182bd", "#d9d9d9"), 5),
          alpha=0.8,
          label="name", label.select=c("Six1", "Ikzf2", "Eya1", "Meis1", "Hoxa9", "Arid1b"),
          repel=TRUE, font.label = c(14, "plain", "blue"),
          add = "reg.line", conf.int = TRUE) + 
  stat_cor(aes(color = interest_group), label.x = 1.5, method="pearson") +
  ylab("Changes in MLL-AF9 binding (logFC; dox vs veh)") +
  xlab("Mean MLL-AF9 signal at promoter (dox)") +
  geom_hline(yintercept = 0, lty=2) +
  theme(legend.title= element_blank())

ggscatter(df_promoter, x = "mean_hmll1_dox", y = "log2FoldChange_3h", color = "interest_group",
          palette = get_palette(c("#54278f", "#de2d26", "#31a354", "#3182bd", "#d9d9d9"), 5),
          alpha=0.8,
          label="name", label.select=c("Six1", "Ikzf2", "Eya1", "Meis1", "Hoxa9", "Arid1b"),
          repel=TRUE, font.label = c(14, "plain", "blue"),
          add = "reg.line", conf.int = TRUE) + 
  stat_cor(aes(color = interest_group), label.x = 1.5, method="pearson") +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 8, lty=2) +
  ylab("Changes in gene expr. (log2FC; 3h dTAG)") +
  xlab("Mean MLL-AF9 signal at promoter (dox)") +
  geom_hline(yintercept = 0, lty=2) +
  theme(legend.title= element_blank(), legend.position = "none")

ggboxplot(df_enhancer, x="peak_mll1_dox", y="dist_to_closest_gene_tss_dox", yscale="log2")
