# Cluster analysis done elsewhere, this script  only meant for producing final figures

# load packages
library(tidyverse)
library(edgeR) # loads limma as a dependency
library(DESeq2)
library(pheatmap) # alternative to ComplexHeatmap
library(cowplot)
theme_set(theme_cowplot()) #white background instead of grey -> don't load if want grey grid
behandeling <- c("dodgerblue2","firebrick3","black") 

threshold <- 0.01 #FDR threshold I want to use ####

# load data and results
load("analysis/_RNAseq/_results/limma_mod_allweeks_filtn2.RData") # all weeks analyzed together

head(dfilt$counts) # raw counts filtered by cutoff ~20/median(library size) in at least 2 replicates
head(dfilt$samples) # sample information incl. library sizes, normalisation factors

rm(res_list.3h, res_list.24h, results.3h, results.24h, v.24h, v.3h, cutoff, filtwhich.n2) # remove object not going to use

DEGs <- read.table(file=paste("analysis/_RNAseq/_results/DEGs_annot_n2_p", threshold,".tsv", sep=""), header=T)
head(DEGs) # DEGs_annot

hm_clust <- read.csv("analysis/_RNAseq/_results/DEGS_hm_clusters.csv")
head(hm_clust)


#-------------------------------------------------------
# HEAT MAPS ####
#-------------------------------------------------------
vst <- vst(dfilt$counts) # all genes

# get means
vst.means <- matrix(c(rowMeans(vst[, c(colnames(vst)[grepl("BEF_",colnames(vst))==T &
                                                       (gsub("BEF_(\\d+)","\\1",colnames(vst)) %in% c("136", "326", "473"))])]), # w2
                      rowMeans(vst[, c(colnames(vst)[grepl("C_3",colnames(vst))==T &
                                                       (gsub("C_3_(\\d+)","\\1",colnames(vst)) %in% c("136", "326", "473"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("C_24",colnames(vst))==T &
                                                       (gsub("C_24_(\\d+)","\\1",colnames(vst)) %in% c("136", "326", "473"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("N_3",colnames(vst))==T &
                                                       (gsub("N_3_(\\d+)","\\1",colnames(vst)) %in% c("136", "326", "473"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("N_24",colnames(vst))==T &
                                                       (gsub("N_24_(\\d+)","\\1",colnames(vst)) %in% c("136", "326", "473"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("W_3",colnames(vst))==T &
                                                       (gsub("W_3_(\\d+)","\\1",colnames(vst)) %in% c("136", "326", "473"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("W_24",colnames(vst))==T &
                                                       (gsub("W_24_(\\d+)","\\1",colnames(vst)) %in% c("136", "326", "473"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("BEF_",colnames(vst))==T &
                                                       (gsub("BEF_(\\d+)","\\1",colnames(vst)) %in% c("102", "219", "367"))])]), # w4
                      rowMeans(vst[, c(colnames(vst)[grepl("C_3",colnames(vst))==T &
                                                       (gsub("C_3_(\\d+)","\\1",colnames(vst)) %in% c("102", "219", "367"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("C_24",colnames(vst))==T &
                                                       (gsub("C_24_(\\d+)","\\1",colnames(vst)) %in% c("102", "219", "367"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("N_3",colnames(vst))==T &
                                                       (gsub("N_3_(\\d+)","\\1",colnames(vst)) %in% c("102", "219", "367"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("N_24",colnames(vst))==T &
                                                       (gsub("N_24_(\\d+)","\\1",colnames(vst)) %in% c("102", "219", "367"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("W_3",colnames(vst))==T &
                                                       (gsub("W_3_(\\d+)","\\1",colnames(vst)) %in% c("102", "219", "367"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("W_24",colnames(vst))==T &
                                                       (gsub("W_24_(\\d+)","\\1",colnames(vst)) %in% c("102", "219", "367"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("BEF_",colnames(vst))==T &
                                                       (gsub("BEF_(\\d+)","\\1",colnames(vst)) %in% c("128", "390", "471"))])]), # w6
                      rowMeans(vst[, c(colnames(vst)[grepl("C_3",colnames(vst))==T &
                                                       (gsub("C_3_(\\d+)","\\1",colnames(vst)) %in% c("128", "390", "471"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("C_24",colnames(vst))==T &
                                                       (gsub("C_24_(\\d+)","\\1",colnames(vst)) %in% c("128", "390", "471"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("N_3",colnames(vst))==T &
                                                       (gsub("N_3_(\\d+)","\\1",colnames(vst)) %in% c("128", "390", "471"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("N_24",colnames(vst))==T &
                                                       (gsub("N_24_(\\d+)","\\1",colnames(vst)) %in% c("128", "390", "471"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("W_3",colnames(vst))==T &
                                                       (gsub("W_3_(\\d+)","\\1",colnames(vst)) %in% c("128", "390", "471"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("W_24",colnames(vst))==T &
                                                       (gsub("W_24_(\\d+)","\\1",colnames(vst)) %in% c("128", "390", "471"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("BEF_",colnames(vst))==T &
                                                       (gsub("BEF_(\\d+)","\\1",colnames(vst)) %in% c("94", "407", "411"))])]), # w8
                      rowMeans(vst[, c(colnames(vst)[grepl("C_3",colnames(vst))==T &
                                                       (gsub("C_3_(\\d+)","\\1",colnames(vst)) %in% c("94", "407", "411"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("C_24",colnames(vst))==T &
                                                       (gsub("C_24_(\\d+)","\\1",colnames(vst)) %in% c("94", "407", "411"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("N_3",colnames(vst))==T &
                                                       (gsub("N_3_(\\d+)","\\1",colnames(vst)) %in% c("94", "407", "411"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("N_24",colnames(vst))==T &
                                                       (gsub("N_24_(\\d+)","\\1",colnames(vst)) %in% c("94", "407", "411"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("W_3",colnames(vst))==T &
                                                       (gsub("W_3_(\\d+)","\\1",colnames(vst)) %in% c("94", "407", "411"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("W_24",colnames(vst))==T &
                                                       (gsub("W_24_(\\d+)","\\1",colnames(vst)) %in% c("94", "407", "411"))])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("BEF_",colnames(vst))==T])]), # means as ref
                      rowMeans(vst[, c(colnames(vst)[grepl("C_3",colnames(vst))==T])]), 
                      rowMeans(vst[, c(colnames(vst)[grepl("C_24",colnames(vst))==T])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("N_3",colnames(vst))==T])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("N_24",colnames(vst))==T])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("W_3",colnames(vst))==T])]),
                      rowMeans(vst[, c(colnames(vst)[grepl("W_24",colnames(vst))==T])]),
                      rowMeans(vst[, c(colnames(vst)[gsub(".+_(\\d+)","\\1",colnames(vst)) %in% c("136", "326", "473")])]), # w2
                      rowMeans(vst[, c(colnames(vst)[gsub(".+_(\\d+)","\\1",colnames(vst)) %in% c("102", "219", "367")])]), # w4
                      rowMeans(vst[, c(colnames(vst)[gsub(".+_(\\d+)","\\1",colnames(vst)) %in% c("128", "390", "471")])]), # w6
                      rowMeans(vst[, c(colnames(vst)[gsub(".+_(\\d+)","\\1",colnames(vst)) %in% c("94", "407", "411")])]),  # w8
                      rowMeans(vst)), ncol=40) # average all samples
colnames(vst.means) <- c("w2_BEF", "w2_C3", "w2_C24", "w2_N3", "w2_N24", "w2_W3", "w2_W24",
                         "w4_BEF", "w4_C3", "w4_C24", "w4_N3", "w4_N24", "w4_W3", "w4_W24",
                         "w6_BEF", "w6_C3", "w6_C24", "w6_N3", "w6_N24", "w6_W3", "w6_W24",
                         "w8_BEF", "w8_C3", "w8_C24", "w8_N3", "w8_N24", "w8_W3", "w8_W24",
                         "av_BEF", "av_C3", "av_C24", "av_N3", "av_N24", "av_W3", "av_W24",
                         "w2", "w4", "w6", "w8", "AV")
rownames(vst.means) <- rownames(vst)
head(vst.means)

vst.10 <- vst.means[,which(colnames(vst.means) %in% c("w2_BEF", "w2_N3", "w2_N24",
                                                      "w4_BEF", "w4_N3", "w4_N24", 
                                                      "w6_BEF", "w6_N3", "w6_N24",
                                                      "w8_BEF", "w8_N3", "w8_N24",
                                                      "av_BEF", "av_N3", "av_N24"))]
head(vst.10)

vst.5 <- vst.means[,which(colnames(vst.means) %in% c("w2_BEF", "w2_C3", "w2_C24",
                                                     "w4_BEF", "w4_C3", "w4_C24",
                                                     "w6_BEF", "w6_C3", "w6_C24",
                                                     "w8_BEF", "w8_C3", "w8_C24",
                                                     "av_BEF", "av_C3", "av_C24"))]

vst.5.rel <- vst.5 / vst.10
head(vst.5.rel)

vst.15 <- vst.means[,which(colnames(vst.means) %in% c("w2_BEF", "w2_W3", "w2_W24",
                                                      "w4_BEF", "w4_W3", "w4_W24",
                                                      "w6_BEF", "w6_W3", "w6_W24",
                                                      "w8_BEF", "w8_W3", "w8_W24",
                                                      "av_BEF", "av_W3", "av_W24"))]

vst.15.rel <- vst.15 / vst.10
head(vst.15.rel)

vst.rel <- cbind(vst.15.rel[,which(colnames(vst.15.rel) %in% c("w2_W3", "w2_W24", "w4_W3", "w4_W24",
                                                               "w6_W3", "w6_W24", "w8_W3", "w8_W24",
                                                               "av_W3", "av_W24"))], 
                 vst.5.rel[,which(colnames(vst.5.rel) %in% c("w2_C3", "w2_C24", "w4_C3", "w4_C24",
                                                             "w6_C3", "w6_C24", "w8_C3", "w8_C24",
                                                             "av_C3", "av_C24"))])
rm(vst.5.rel, vst.15.rel)
head(vst.rel)


# Heat map per contrast set
sets <- unique(DEGs$Comparison)

hms_3h <- list()
hms_24h <- list()

# Pheatmap aesthetics ####
ann_colors <- list(Week=c(Week2="lavender", Week4="lightsteelblue2", Week6="lavenderblush3", Week8="lightsteelblue4"),
                   Treatment=c(Cold="dodgerblue2", Warm="firebrick3"))
# max 5 clusters
cl.col <- list(Cluster=c(cl1="grey10", cl2="grey50", cl3="grey30", cl4="grey90", cl5="grey70"))
names(cl.col) <- "Cluster"
cl.col[[1]][1]

for(set in 1:length(sets)){
  setname <- sets[[set]]
  
  # for after 3h
  #-------------------------------------------------------
  clustinf.3h <- filter(hm_clust, set==setname & Comparison=="3h")
  nrclust.3h <- clustinf.3h$Cluster %>% unique() %>% length()
  head(clustinf.3h)
  
  vst.3h.relmean <- vst.rel[rownames(vst.rel) %in% unique(filter(DEGs, Comparison==sets[set] & (C3==1 | W3==1))$GeneID),]
  head(vst.3h.relmean)
  
  # Get z-scores
  z.3h.mean <- t(scale(t(vst.3h.relmean), center=TRUE, scale=TRUE)) # pick counts or counts relative to BEF
  range(z.3h.mean)
  head(z.3h.mean)
  
  
  # get column order
  colDat <- subset(dfilt$samples, Timepoint=="3" | Timepoint=="BEF")
  colorder <- as.data.frame(colDat) %>% 
    mutate(WeekTim = ifelse(Timepoint=="BEF", paste(Treat_week, Timepoint, sep="_"), # In case want to plot means instead of all samples
                            paste(Treat_week, "_", Treatment,Timepoint, sep="")),
           Treatment=factor(Treatment, levels=c("C", "N", "W"))) %>% 
    arrange(Trw_num, Treatment, Timepoint, Tube)
  
  colorder.3h.mean <- subset(colorder, Timepoint=="3" & Treatment!="N") %>% select(WeekTim, Treat_week, Timepoint, Treatment) %>% unique # to plot means
  annotation.3h.mean <- colorder.3h.mean %>% mutate(Week=gsub("w(\\d)", "Week\\1", Treat_week), Treatment=ifelse(Treatment=="C", "Cold", "Warm")) %>% 
    select(Treatment, Week)# without BEF
  rownames(annotation.3h.mean) <- colorder.3h.mean$WeekTim
  annotation.3h.mean
  
  row.annot.3h <- data.frame(Cluster=paste("cl", clustinf.3h$Cluster, sep=""))
  rownames(row.annot.3h) <- clustinf.3h$GeneID
  head(row.annot.3h)
  
  ann_colors.3h <- c(ann_colors, list(c(cl.col[[1]][c(1:nrclust.3h)])))
  names(ann_colors.3h) <- c(names(ann_colors), names(cl.col))
  
  # Make heatmap
  hm.3h.mean <- pheatmap::pheatmap(z.3h.mean[,c(colorder.3h.mean$WeekTim)], 
                                   #color=my.colors,
                                   annotation_colors = ann_colors.3h,
                                   cluster_rows=TRUE, show_rownames=FALSE, annotation_row = row.annot.3h,
                                   cluster_cols=FALSE, annotation_col=annotation.3h.mean, show_colnames=FALSE, gaps_col=cumsum(c(2,2,2,2)),
                                   clustering_distance_rows="correlation", #Pearson correlation
                                   cutree_rows=nrclust.3h,
                                   legend_breaks=c(min(z.3h.mean)+0.5,0,max(z.3h.mean)-0.5), legend_labels=c("", "", ""),
                                   fontsize=12, border_color=NA) 
  hms_3h[[set]] <-  hm.3h.mean
  names(hms_3h)[set] <- sets[set]
  plot_name <- paste(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/pretty/Hmap_DEGs_3h_means_", sets[set], sep=""), ".jpeg", sep="")# name for file
  jpeg(plot_name, width=200, height=200, units="mm", res=300, quality=100)
  print(hm.3h.mean)
  dev.off()  
  
  
  # for after 24h 
  #-------------------------------------------------------
  clustinf.24h <- filter(hm_clust, set==setname & Comparison=="24h")
  nrclust.24h <- clustinf.24h$Cluster %>% unique() %>% length()
  head(clustinf.24h)
  
  vst.24h.relmean <- vst.rel[rownames(vst.rel) %in% unique(filter(DEGs, Comparison==sets[set] & (C24==1 | W24==1))$GeneID),]
  head(vst.24h.relmean)
  
  
  # Get z-scores
  z.24h.mean <- t(scale(t(vst.24h.relmean), center=TRUE, scale=TRUE)) # pick counts or counts relative to BEF
  range(z.24h.mean)
  head(z.24h.mean)
  
  
  # get column order
  colDat <- subset(dfilt$samples, Timepoint=="24" | Timepoint=="BEF")
  colorder <- as.data.frame(colDat) %>% 
    mutate(WeekTim = ifelse(Timepoint=="BEF", paste(Treat_week, Timepoint, sep="_"), # In case want to plot means instead of all samples
                            paste(Treat_week, "_", Treatment,Timepoint, sep="")),
           Treatment=factor(Treatment, levels=c("C", "N", "W"))) %>% 
    arrange(Trw_num, Treatment, Timepoint, Tube)

  colorder.24h.mean <- subset(colorder, Timepoint=="24" & Treatment!="N") %>% select(WeekTim, Treat_week, Timepoint, Treatment) %>% unique # to plot means
  annotation.24h.mean <- colorder.24h.mean %>% mutate(Week=gsub("w(\\d)", "Week\\1", Treat_week), Treatment=ifelse(Treatment=="C", "Cold", "Warm")) %>% 
    select(Treatment, Week)# without BEF
  rownames(annotation.24h.mean) <- colorder.24h.mean$WeekTim
  head(annotation.24h.mean)

  row.annot.24h <- data.frame(Cluster=paste("cl", clustinf.24h$Cluster, sep=""))
  rownames(row.annot.24h) <- clustinf.24h$GeneID
  head(row.annot.24h)
  
  ann_colors.24h <- c(ann_colors, list(c(cl.col[[1]][c(1:nrclust.24h)])))
  names(ann_colors.24h) <- c(names(ann_colors), names(cl.col))
  
  # With pheatmap
  hm.24h.mean <- pheatmap::pheatmap(z.24h.mean[,c(colorder.24h.mean$WeekTim)], 
                                    annotation_colors = ann_colors.24h, #color=my.colors,
                                    cluster_rows=TRUE, show_rownames=FALSE, annotation_row = row.annot.24h,
                                    cluster_cols=FALSE, annotation_col=annotation.24h.mean, show_colnames=FALSE, gaps_col=cumsum(c(2,2,2,2)),
                                    clustering_distance_rows="correlation", #Pearson correlation
                                    cutree_rows=nrclust.24h,
                                    legend_breaks=c(min(z.24h.mean)+0.5,0,max(z.24h.mean)-0.5), legend_labels=c("", "", ""),
                                    fontsize=12, border_color=NA)
  
  hms_24h[[set]] <-  hm.24h.mean
  names(hms_24h)[set] <- sets[set]
  plot_name <- paste(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/pretty/Hmap_DEGs_24h_means_", sets[set], sep=""), ".jpeg", sep="")# name for file
  jpeg(plot_name, width=200, height=200, units="mm", res=300, quality=100)
  print(hm.24h.mean)
  dev.off()  
  
}

names(hms_3h) # stored heatmap trees
names(hms_24h)
