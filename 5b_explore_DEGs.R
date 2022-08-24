#-------------------------------------------------------
#-------------------------------------------------------
# Explore RNAseq analysis results of egg samples from Wintermoth Transfer experiment
# Analysis Method: DE with limma
#-------------------------------------------------------
#-------------------------------------------------------

# Open R project in top folder

#-------------------------------------------------------
# Set up environment ####
#-------------------------------------------------------

# load packages
library(tidyverse)
library(edgeR) # loads limma as a dependency
library(DESeq2)
library(pheatmap) # alternative to ComplexHeatmap
library(cowplot)
theme_set(theme_cowplot()) #white background instead of grey -> don't load if want grey grid
behandeling <- c("dodgerblue2","firebrick3","black") 


# load data and results
#-------------------------------------------------------
load("analysis/_RNAseq/_results/limma_mod_allweeks_filtn2.RData") # all weeks analyzed together

head(dfilt$counts) # raw counts filtered by cutoff ~20/median(library size) in at least 2 replicates
head(dfilt$samples) # sample information incl. library sizes, normalisation factors
lapply(res_list.3h, head) # results for 22 contrasts 3h model
lapply(res_list.24h, head) # results for 22 contrasts 24h model

res_lists <- list(res_list.3h, res_list.24h)
names(res_lists) <- c("3h", "24h")
rm(res_list.3h, res_list.24h)

# annotation
annot <- read.table("analysis/_RNAseq/_annotation/genesAnnot_wblastp_rows_exp.tsv", header=T)

# DEGs annotated
DEGs_annot <- read.table(file="analysis/_RNAseq/_results/DEGs_annot_n2_p0.01.tsv", header=T, sep="\t")

# WGCNA analysis: description 10C samples
load("analysis/_RNAseq/_results/WGCNA_Results_CoExpr10C.Rdata")
head(modNames) # ANOVA results for each module
head(SiteInfo.sign) # genes assigned to significant modules with significant membership 
head(SiteInfo) # results for all genes

# Save WGCNA results for significant modules and whether overlap with DE or not
# write.csv(modNames, file="analysis/_RNAseq/_results/WGCNA_ANOVA-results.csv", row.names = F )

SiteInfo.sign <- SiteInfo.sign %>% 
  mutate(DE_Overlap=ifelse(gene %in% unique(DEGs_annot$GeneID), "Yes", "No"), SignMM="Yes")

WGCNA <- filter(SiteInfo[,c(colnames(SiteInfo) %in% colnames(SiteInfo.sign))], Module %in% filter(modNames, signMod==1)$Module) 
WGCNA <- merge(WGCNA, SiteInfo.sign[,c("gene", "Module", "DE_Overlap", "SignMM")], by=c("gene", "Module"), all.x=T)
WGCNA$SignMM <- ifelse(is.na(WGCNA$SignMM)==T, "No", WGCNA$SignMM)
head(WGCNA)

table(WGCNA$SignMM, useNA="always")
table(WGCNA$DE_Overlap, useNA="always")

# write.csv(WGCNA[,c(1:3,13,4:12)], file="analysis/_RNAseq/_results/WGCNA_MMgenes-results.csv", row.names = F)


# Overlap analyses ####
#-------------------------------------------------------
overview <- list(unique(SiteInfo.sign$gene), unique(DEGs_annot$GeneID))
names(overview) <- c("WGCNA analysis", "DE analysis")

venn <- ggVennDiagram::ggVennDiagram(overview, # needs to be a list, with each item = vector of genes names
                                     set_size=7, 
                                     label="count", label_alpha = 0, label_size=6)+
  scale_fill_gradient(low="white",high ="white") + # "mediumorchid4")
  theme(legend.position="none", text = element_text(vjust=6)) + scale_y_continuous(expand = expansion(mult = .1))
png(paste("analysis/_RNAseq/_plots/Venn_overlap_analyses.png", sep=""), width=100, height=60, units="mm", res=300)
print(venn)
dev.off()

WGCNAresults <- list(filter(SiteInfo.sign, Module=="Mod10")$gene, filter(SiteInfo.sign, Module=="Mod8")$gene,
                 filter(SiteInfo.sign, Module=="Mod13")$gene, filter(SiteInfo.sign, Module=="Mod3")$gene)
names(WGCNAresults) <- c("Mod10.high.w2", "Mod8.low.w2", "Mod13.down", "Mod3.up")

splitview <- c(WGCNAresults, list(unique(DEGs_annot$GeneID)))
names(splitview) <- c("Mod10.high.w2", "Mod8.low.w2", "Mod13.down", "Mod3.up", "DEGs")

vennsplit <- ggVennDiagram::ggVennDiagram(splitview, # needs to be a list, with each item = vector of genes names
                                     label_alpha = 0)+
  scale_fill_gradient(low="white",high = "mediumorchid4")
png(paste("analysis/_RNAseq/_plots/Venn_overlap_analyses_split.png", sep=""), width=250, height=250, units="mm", res=300)
print(vennsplit)
dev.off()

# Overall_Temp DEGS
results1 <-c(WGCNAresults, list(unique(filter(DEGs_annot, Comparison=="Overall_Temp")$GeneID)))
names(results1) <- c("Mod10.high.w2", "Mod8.low.w2", "Mod13.down", "Mod3.up", "DEGs_Overall.Temp")

venn1 <- ggVennDiagram::ggVennDiagram(results1, # needs to be a list, with each item = vector of genes names
                                         label_alpha = 0)+
  scale_fill_gradient(low="white",high = "mediumorchid4")
png(paste("analysis/_RNAseq/_plots/Venn_overlap_analyses1.png", sep=""), width=250, height=150, units="mm", res=300)
print(venn1)
dev.off()

# Within-week DEGS
results2 <-c(WGCNAresults, list(unique(filter(DEGs_annot, Comparison=="Within_week")$GeneID)))
names(results2) <- c("Mod10.high.w2", "Mod8.low.w2", "Mod13.down", "Mod3.up", "DEGs_WithinWeek")

venn2 <- ggVennDiagram::ggVennDiagram(results2, # needs to be a list, with each item = vector of genes names
                                      label_alpha = 0)+
  scale_fill_gradient(low="white",high = "mediumorchid4")
png(paste("analysis/_RNAseq/_plots/Venn_overlap_analyses2.png", sep=""), width=250, height=150, units="mm", res=300)
print(venn2)
dev.off()

results2a <-c(WGCNAresults, list(unique(filter(DEGs_annot, origin=="w2")$GeneID)), list(unique(filter(DEGs_annot, origin=="w4")$GeneID)))
names(results2a) <- c("Mod10.high.w2", "Mod8.low.w2", "Mod13.down", "Mod3.up", "DEGs_Week2", "DEGs_Week4")

venn2a <- ggVennDiagram::ggVennDiagram(results2a, # needs to be a list, with each item = vector of genes names
                                      label_alpha = 0)+
  scale_fill_gradient(low="white",high = "mediumorchid4")
png(paste("analysis/_RNAseq/_plots/Venn_overlap_analyses2early.png", sep=""), width=300, height=250, units="mm", res=300)
print(venn2a)
dev.off()

results2b <-c(WGCNAresults, list(unique(filter(DEGs_annot, origin=="w6")$GeneID)), list(unique(filter(DEGs_annot, origin=="w8")$GeneID)))
names(results2b) <- c("Mod10.high.w2", "Mod8.low.w2", "Mod13.down", "Mod3.up", "DEGs_Week6", "DEGs_Week8")

venn2b <- ggVennDiagram::ggVennDiagram(results2b, # needs to be a list, with each item = vector of genes names
                                       label_alpha = 0)+
  scale_fill_gradient(low="white",high = "mediumorchid4")
png(paste("analysis/_RNAseq/_plots/Venn_overlap_analyses2late.png", sep=""), width=300, height=250, units="mm", res=300)
print(venn2b)
dev.off()

# Unexpected overlap
w2 <- filter(DEGs_annot, origin=="w2") %>% filter(GeneID %in% filter(SiteInfo.sign, Module=="Mod3")$gene) %>% arrange(GeneID, Evalue)
w4 <- filter(DEGs_annot, origin=="w4") %>% filter(GeneID %in% filter(SiteInfo.sign, Module=="Mod3")$gene) %>% arrange(GeneID, Evalue)
early <- rbind(w2, w4)
early$ModuleOverlap <- "Module3"

w6 <- filter(DEGs_annot, origin=="w6") %>% filter(GeneID %in% filter(SiteInfo.sign, Module=="Mod13")$gene)%>% arrange(GeneID, Evalue)
w8 <- filter(DEGs_annot, origin=="w8") %>% filter(GeneID %in% filter(SiteInfo.sign, Module=="Mod13")$gene)%>% arrange(GeneID, Evalue)
late1 <- rbind(w6, w8)
late1$ModuleOverlap <- "Module13"

w6 <- filter(DEGs_annot, origin=="w6") %>% filter(GeneID %in% filter(SiteInfo.sign, Module=="Mod10")$gene)%>% arrange(GeneID, Evalue)
w8 <- filter(DEGs_annot, origin=="w8") %>% filter(GeneID %in% filter(SiteInfo.sign, Module=="Mod10")$gene)%>% arrange(GeneID, Evalue)
late2 <- rbind(w6, w8)
late2$ModuleOverlap <- "Module10"

unxoverlap <- rbind(early, late1, late2)
rm(early, late1, late2)
head(unxoverlap)

hmDEGs <- read.csv("analysis/_RNAseq/_results/DEGS_hm_clusters_annot.csv")
hmDEGs.sub <- filter(hmDEGs, Set=="Within_week" & GeneID %in% unxoverlap$GeneID)
hmDEGs.sub <- merge(hmDEGs.sub, unique(unxoverlap[,c("GeneID", "ModuleOverlap")]), by="GeneID")
head(hmDEGs.sub)
# write.table(hmDEGs.sub, file="analysis/_RNAseq/_results/unexpected_overlap_genes.tsv", sep="\t", row.names=F)
filter(hmDEGs.sub, grepl("w8_", contr)==T & (ModuleOverlap=="Module13"|ModuleOverlap=="Module10" ))$GeneID %>% unique %>% length()

# Between-week DEGS
results3 <- list(filter(SiteInfo.sign, Module=="Mod10")$gene, filter(SiteInfo.sign, Module=="Mod8")$gene,
                 filter(SiteInfo.sign, Module=="Mod13")$gene, filter(SiteInfo.sign, Module=="Mod3")$gene,
                 unique(filter(DEGs_annot, Comparison=="Between_week")$GeneID))
names(results3) <- c("Mod10.high.w2", "Mod8.low.w2", "Mod13.down", "Mod3.up", "DEGs_BetweenWeek")

venn3 <- ggVennDiagram::ggVennDiagram(results3, # needs to be a list, with each item = vector of genes names
                                      label_alpha = 0)+
  scale_fill_gradient(low="white",high = "mediumorchid4")
png(paste("analysis/_RNAseq/_plots/Venn_overlap_analyses3.png", sep=""), width=250, height=150, units="mm", res=300)
print(venn3)
dev.off()


# Extract all info I need for plotting ####
#-------------------------------------------------------
threshold <- 0.01 #FDR threshold I want to use ####

# Want plots per contrast set (3x) per comparison (CvsN and WvsN) per timepoint (3h and 24h)

# Denote for each result mod + contrast origin
for(list in 1:2){
  res_list <- res_lists[[list]]
  
  for(item in 1:length(res_list)){
    res <- res_list[[item]]
    res$mod <- names(res_lists)[list]
    res$contr <- names(res_list)[item]
    res$GeneID <- rownames(res)
    res_list[[item]] <- res
  }
  
  res_lists[[list]] <- data.table::rbindlist(res_list)
}
lapply(res_lists, head)
res <- data.table::rbindlist(res_lists)
#res <- filter(res, padj<threshold) # I will also need the non sign ones for many plots

# Get baseMean = log2 normalized mean for all samples
lcpm <- cpm(dfilt, log=TRUE) # log2 library size normalized counts
log2means <- as.data.frame(rowMeans(lcpm))
colnames(log2means) <- "log2Mean"
log2means$GeneID <- rownames(lcpm)
head(log2means)

res <- merge(res, log2means, by="GeneID")
head(res)

# Denote set and comparison
contr <- DEGs_annot[,c("origin", "Comparison")] %>% unique
contr

res <- res %>% mutate(Set=ifelse(grepl("\\w\\.", contr)==T,
                                    "Between_week", ifelse(grepl("w\\d_\\wvsN",contr)==T, "Within_week", "Overall_Temp")),
                      Comparison=ifelse(grepl("\\w\\.", contr)==T, gsub("(\\w)\\.w\\dvs\\d", "\\1vsN",contr), 
                                        ifelse(grepl("w\\d_\\wvsN",contr)==T, gsub("w\\d_(\\wvsN)", "\\1",contr), contr))) %>%
  arrange(padj, desc(abs(log2FoldChange)))
table(res$Comparison)
table(res$Set)
head(res)

res_annot <- merge(unique(annot[,c("GeneID", "Read_frame", "GeneName", "Evalue", "Prot_descr", "Species")]), res, by=c("GeneID"))
res_annot <- res_annot %>% select(GeneID, Set, mod, Comparison, contr, log2Mean, log2FoldChange, padj, GeneName, Prot_descr, Evalue, Read_frame, Species) %>%
  arrange(padj, desc(abs(log2FoldChange)), GeneID, Evalue)
res_annot <- res_annot[!duplicated(res_annot[,c(1:8,10)]),]
head(res_annot)

res_annot_DEGs <- filter(res_annot, padj<threshold)
res_ot <-filter(res, Set=="Overall_Temp")
res_ww <-filter(res, Set=="Within_week")
res_bw <-filter(res, Set=="Between_week")

# write.table(res_annot_DEGs, file=paste("analysis/_RNAseq/_results/all_results_n2_p", threshold,".tsv", sep=""), row.names=F, sep="\t")
# write.table(res_ot, file=paste("analysis/_RNAseq/_results/all_results_n2_OverallTemp.tsv", sep=""), row.names=F, sep="\t")
# write.table(res_ww, file=paste("analysis/_RNAseq/_results/all_results_n2_WithinWeek.tsv", sep=""), row.names=F, sep="\t")
# write.table(res_bw, file=paste("analysis/_RNAseq/_results/all_results_n2_BetweenWeek.tsv", sep=""), row.names=F, sep="\t")


#-------------------------------------------------------
# MA PLOTS ####
#-------------------------------------------------------

# log2 normalized mean for all samples vs. Log2FoldChange
# NB: takes into account how highly the gene is expressed for which you observe 
#  a particular logFoldchange.

# Check if results behave as you would expect, so make per contrast!
for(c in 1:length(unique(res$contr))){
  contrast <- unique(res$contr)[c]
  
  for(m in 1:length(unique(res$mod))){
    model <- unique(res$mod)[m]
    
    res1 <- filter(res, mod==model & contr==contrast) %>%
      mutate(colour=ifelse(padj > threshold, "NS", ifelse(log2FoldChange>0, "UP", "DOWN"))) %>%
      mutate(colour=factor(colour, levels=c("UP", "DOWN", "NS")))
    head(res1)
    table(res1$colour)
    table(res1$contr)
    
    # with ggplot
    # leave out NS for now, super messy + bit double
    
    plot_title <- paste(model, contrast, sep="_")
    
    plot <- ggplot(res1, aes(x = log2Mean, y=log2FoldChange, col=colour)) + # looks better with shrunken results
      #scale_colour_manual(values=c("coral", "royalblue3", "grey"))+
      geom_point(data=filter(res1, colour=="NS"),size=2, col="darkgrey")+
      geom_point(data=filter(res1, colour!="NS"),size=3)+
      labs(x="Log2 Mean of normalised counts", y="Log fold change") +
      ggtitle(plot_title)+
      geom_hline(yintercept=0, col="grey", size=2)+
      #coord_cartesian(xlim = c(-5, 20), ylim=c(-30,30)) +
      scale_y_continuous(breaks=seq(-30,30, by=2))+ scale_x_continuous(breaks=seq(-20,20, by=2))+
      guides(shape = "none")+
      theme(legend.title= element_blank(), legend.text=element_text(size=18), axis.title = element_text(size=20), 
            axis.text=element_text(size=16), plot.title = element_text(size=18))
    
    plot_name <- paste(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold, "/MA/MA_", plot_title, sep=""), ".png", sep="") # name for file
    png(plot_name, width=250, height=150, units="mm", res=300)
    print(plot)
    dev.off() #close png writer
  }
}
rm(s, set, m, model, c, comp, res1, plot_title, plot, plot_name)


#-------------------------------------------------------
# VOLCANO PLOTS ####
#-------------------------------------------------------
# Log2FoldChange vs. -log10(padj)

for(s in 1:length(unique(res$Set))){
  set <- unique(res$Set)[s]
  
  for(m in 1:length(unique(res$mod))){
    model <- unique(res$mod)[m]
    
    for(c in 1:length(unique(res$Comparison))){
      comp <- unique(res$Comparison)[c]
      
      res1 <- filter(res, Set==set & mod==model & Comparison==comp) %>%
        mutate(log10=-log10(padj), colour=ifelse(padj > threshold, "NS", ifelse(log2FoldChange>0, "UP", "DOWN"))) %>%
        mutate(colour=factor(colour, levels=c("UP", "DOWN", "NS")))
      head(res1)
      table(res1$colour)
      
      # with ggplot
      # leave out NS for now, super messy + bit double
      
      plot_title <- paste(set, model, comp, sep="_")
      
      plot <- ggplot(res1, aes(x = log2FoldChange, y=log10, col=colour)) + 
        scale_colour_manual(values=c("coral", "royalblue3", "grey"))+
        geom_point(size=2) +
        geom_hline(yintercept=-log10(threshold), linetype="dashed", col="grey")+
        labs(x="Log fold change", y="-log10(FDR)")+
        ggtitle(plot_title)+
        coord_cartesian(ylim = c(0, 12.5), xlim=c(-24,24)) +
        scale_x_continuous(breaks=seq(-30,30, by=5))+ scale_y_continuous(breaks=seq(0,12.5, by=1))+
        guides(shape = "none")+
        theme(legend.title= element_blank(), legend.text=element_text(size=18), axis.title = element_text(size=20), 
              axis.text=element_text(size=16), plot.title = element_text(size=22))
      
      plot_name <- paste(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold, "/Volcano_", plot_title, sep=""), ".png", sep="") # name for file
      png(plot_name, width=250, height=150, units="mm", res=300)
      print(plot)
      dev.off() #close png writer
    }
  }
}
rm(s, set, m, model, c, comp, res1, plot_title, plot, plot_name)


#-------------------------------------------------------
# Per gene expression patterns ####
#-------------------------------------------------------
# Loop to plot all 837 DEGs, z-scores of vst transformed expression
# Plot 4 panel figure; expression for each gene in each week
head(res)
head(dfilt$counts)
vst <- vst(dfilt$counts)
head(vst)

# use z-scores for expression value ####
z <- t(scale(t(vst), center=TRUE, scale=TRUE))
range(z)
head(z)

z.ind <- t(z[rownames(z) %in% unique(DEGs_annot$GeneID),]) %>% as.data.frame
z.ind$Sample <- rownames(z.ind)
z.ind <- merge(z.ind, dfilt$samples[,c("Sample", "Trw_num", "Treatment","Timepoint")], by="Sample")
z.ind <- z.ind %>% mutate(Treat_week=paste("Week", Trw_num, sep=""),
                          Hour=ifelse(Timepoint=="BEF", "0", Timepoint),
                          Hour_num=ifelse(Timepoint=="BEF", 1, ifelse(Timepoint=="3", 2, 6))) %>%
  mutate(Treatment=factor(Treatment, levels=c("C", "W", "N")), Hour=factor(Hour, levels=c("0", "3", "24")))
levels(z.ind$Treatment)
levels(z.ind$Hour)
unique(z.ind$Treat_week)
unique(z.ind$Hour_num)

# get mean values per treatment group per week
z.means <- matrix(c(rowMeans(z[, c(colnames(z)[grepl("BEF_",colnames(z))==T &
                                                       (gsub("BEF_(\\d+)","\\1",colnames(z)) %in% c("136", "326", "473"))])]), # w2
                      rowMeans(z[, c(colnames(z)[grepl("C_3",colnames(z))==T &
                                                       (gsub("C_3_(\\d+)","\\1",colnames(z)) %in% c("136", "326", "473"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("C_24",colnames(z))==T &
                                                       (gsub("C_24_(\\d+)","\\1",colnames(z)) %in% c("136", "326", "473"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("N_3",colnames(z))==T &
                                                       (gsub("N_3_(\\d+)","\\1",colnames(z)) %in% c("136", "326", "473"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("N_24",colnames(z))==T &
                                                       (gsub("N_24_(\\d+)","\\1",colnames(z)) %in% c("136", "326", "473"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("W_3",colnames(z))==T &
                                                       (gsub("W_3_(\\d+)","\\1",colnames(z)) %in% c("136", "326", "473"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("W_24",colnames(z))==T &
                                                       (gsub("W_24_(\\d+)","\\1",colnames(z)) %in% c("136", "326", "473"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("BEF_",colnames(z))==T &
                                                       (gsub("BEF_(\\d+)","\\1",colnames(z)) %in% c("102", "219", "367"))])]), # w4
                      rowMeans(z[, c(colnames(z)[grepl("C_3",colnames(z))==T &
                                                       (gsub("C_3_(\\d+)","\\1",colnames(z)) %in% c("102", "219", "367"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("C_24",colnames(z))==T &
                                                       (gsub("C_24_(\\d+)","\\1",colnames(z)) %in% c("102", "219", "367"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("N_3",colnames(z))==T &
                                                       (gsub("N_3_(\\d+)","\\1",colnames(z)) %in% c("102", "219", "367"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("N_24",colnames(z))==T &
                                                       (gsub("N_24_(\\d+)","\\1",colnames(z)) %in% c("102", "219", "367"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("W_3",colnames(z))==T &
                                                       (gsub("W_3_(\\d+)","\\1",colnames(z)) %in% c("102", "219", "367"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("W_24",colnames(z))==T &
                                                       (gsub("W_24_(\\d+)","\\1",colnames(z)) %in% c("102", "219", "367"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("BEF_",colnames(z))==T &
                                                       (gsub("BEF_(\\d+)","\\1",colnames(z)) %in% c("128", "390", "471"))])]), # w6
                      rowMeans(z[, c(colnames(z)[grepl("C_3",colnames(z))==T &
                                                       (gsub("C_3_(\\d+)","\\1",colnames(z)) %in% c("128", "390", "471"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("C_24",colnames(z))==T &
                                                       (gsub("C_24_(\\d+)","\\1",colnames(z)) %in% c("128", "390", "471"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("N_3",colnames(z))==T &
                                                       (gsub("N_3_(\\d+)","\\1",colnames(z)) %in% c("128", "390", "471"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("N_24",colnames(z))==T &
                                                       (gsub("N_24_(\\d+)","\\1",colnames(z)) %in% c("128", "390", "471"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("W_3",colnames(z))==T &
                                                       (gsub("W_3_(\\d+)","\\1",colnames(z)) %in% c("128", "390", "471"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("W_24",colnames(z))==T &
                                                       (gsub("W_24_(\\d+)","\\1",colnames(z)) %in% c("128", "390", "471"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("BEF_",colnames(z))==T &
                                                       (gsub("BEF_(\\d+)","\\1",colnames(z)) %in% c("94", "407", "411"))])]), # w8
                      rowMeans(z[, c(colnames(z)[grepl("C_3",colnames(z))==T &
                                                       (gsub("C_3_(\\d+)","\\1",colnames(z)) %in% c("94", "407", "411"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("C_24",colnames(z))==T &
                                                       (gsub("C_24_(\\d+)","\\1",colnames(z)) %in% c("94", "407", "411"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("N_3",colnames(z))==T &
                                                       (gsub("N_3_(\\d+)","\\1",colnames(z)) %in% c("94", "407", "411"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("N_24",colnames(z))==T &
                                                       (gsub("N_24_(\\d+)","\\1",colnames(z)) %in% c("94", "407", "411"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("W_3",colnames(z))==T &
                                                       (gsub("W_3_(\\d+)","\\1",colnames(z)) %in% c("94", "407", "411"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("W_24",colnames(z))==T &
                                                       (gsub("W_24_(\\d+)","\\1",colnames(z)) %in% c("94", "407", "411"))])]),
                      rowMeans(z[, c(colnames(z)[grepl("BEF_",colnames(z))==T])]), # means as ref
                      rowMeans(z[, c(colnames(z)[grepl("C_3",colnames(z))==T])]), 
                      rowMeans(z[, c(colnames(z)[grepl("C_24",colnames(z))==T])]),
                      rowMeans(z[, c(colnames(z)[grepl("N_3",colnames(z))==T])]),
                      rowMeans(z[, c(colnames(z)[grepl("N_24",colnames(z))==T])]),
                      rowMeans(z[, c(colnames(z)[grepl("W_3",colnames(z))==T])]),
                      rowMeans(z[, c(colnames(z)[grepl("W_24",colnames(z))==T])]),
                      rowMeans(z[, c(colnames(z)[gsub(".+_(\\d+)","\\1",colnames(z)) %in% c("136", "326", "473")])]), # w2
                      rowMeans(z[, c(colnames(z)[gsub(".+_(\\d+)","\\1",colnames(z)) %in% c("102", "219", "367")])]), # w4
                      rowMeans(z[, c(colnames(z)[gsub(".+_(\\d+)","\\1",colnames(z)) %in% c("128", "390", "471")])]), # w6
                      rowMeans(z[, c(colnames(z)[gsub(".+_(\\d+)","\\1",colnames(z)) %in% c("94", "407", "411")])]),  # w8
                      rowMeans(z)), ncol=40) # average all samples
colnames(z.means) <- c("w2_BEF", "w2_C3", "w2_C24", "w2_N3", "w2_N24", "w2_W3", "w2_W24",
                         "w4_BEF", "w4_C3", "w4_C24", "w4_N3", "w4_N24", "w4_W3", "w4_W24",
                         "w6_BEF", "w6_C3", "w6_C24", "w6_N3", "w6_N24", "w6_W3", "w6_W24",
                         "w8_BEF", "w8_C3", "w8_C24", "w8_N3", "w8_N24", "w8_W3", "w8_W24",
                         "av_BEF", "av_C3", "av_C24", "av_N3", "av_N24", "av_W3", "av_W24",
                         "w2", "w4", "w6", "w8", "AV")
rownames(z.means) <- rownames(z)
head(z.means)

# use z-scores for expression value ####
# z.mat <- t(scale(t(vst.means), center=TRUE, scale=TRUE)) # base R functions expect rows to be samples, while all RNAseq tools expect them as columns
# range(z.mat)
# head(z.mat)


plotDat <- z.means %>%
  as.data.frame() %>%
  mutate(w2_CBEF=w2_BEF, w2_WBEF=w2_BEF, 
         w4_CBEF=w4_BEF, w4_WBEF=w4_BEF, 
         w6_CBEF=w6_BEF, w6_WBEF=w6_BEF, 
         w8_CBEF=w8_BEF, w8_WBEF=w8_BEF) %>%
  select(w2_BEF, w2_CBEF, w2_WBEF, w2_N3, w2_N24, w2_C3, w2_C24, w2_W3, w2_W24,
         w4_BEF, w4_CBEF, w4_WBEF, w4_N3, w4_N24, w4_C3, w4_C24, w4_W3, w4_W24,
         w6_BEF, w6_CBEF, w6_WBEF, w6_N3, w6_N24, w6_C3, w6_C24, w6_W3, w6_W24,
         w8_BEF, w8_CBEF, w8_WBEF, w8_N3, w8_N24, w8_C3, w8_C24, w8_W3, w8_W24) %>%
  t() %>% as.data.frame()

plotDat <- plotDat[, c(colnames(plotDat) %in% unique(DEGs_annot$GeneID))]
plotDat$Timepoint <- rownames(plotDat)

# pull Treatment and Time info from row.names
plotDat <- plotDat %>%  mutate(Treatment=gsub("w\\d_(.+)","\\1",rownames(plotDat)), Hour=gsub("w\\d_(.+)","\\1",rownames(plotDat)),
                               Treat_week=gsub("w(\\d)_.+","Week\\1",rownames(plotDat)))
plotDat <- plotDat %>%  mutate(Treatment=gsub("(\\w)\\d+","\\1",Treatment), Hour=gsub("\\w(\\d+)","\\1",Hour))
plotDat <- plotDat %>%  mutate(Treatment=as.factor(ifelse(grepl("C", Treatment)==T, "C", ifelse(grepl("W", Treatment)==T, "W","N"))),
                               Hour_num=as.numeric(ifelse(grepl("BEF", Hour)==T, 1, ifelse(grepl("3", Hour)==T, 2, 6))),
                               Hour=as.factor(ifelse(grepl("BEF", Hour)==T, "0", Hour)))
plotDat <- plotDat %>%  mutate(Treatment=factor(plotDat$Treatment, levels=c("C", "W", "N")), Hour=factor(plotDat$Hour, levels=c("0", "3", "24"))) %>%
  arrange(Treatment) # use this order,otherwise 0h point colored red in graph, don't know why...
levels(plotDat$Treatment)
levels(plotDat$Hour)
unique(plotDat$Treat_week)
unique(plotDat$Hour_num)

for(gene in 1:837){
  g <- colnames(plotDat)[gene]
  Dat <- plotDat[,c("Treat_week", "Treatment", "Hour", "Hour_num", g)]
  colnames(Dat)[5] <- "y"
  head(Dat)
  
  ind.Dat <- z.ind[,c("Treat_week", "Treatment", "Hour", "Hour_num", g)]
  colnames(ind.Dat)[5] <- "y"
  
  # Add sign info
  g_info <- filter(DEGs_annot, GeneID==g) %>% select(GeneID, origin, C3, W3, C24, W24, Comparison) %>% unique()
  info <- paste(g_info[1,"Comparison"], " ", g_info[1, "origin",], ", Sign: ", sep="")
  
  if(g_info[1,"C3"]==1){info <- paste(info, "C3")}
  if(g_info[1,"W3"]==1){info <- paste(info, "W3")}
  if(g_info[1,"C24"]==1){info <- paste(info, "C24")}
  if(g_info[1,"W24"]==1){info <- paste(info, "W24")}
  
  if(nrow(g_info)>1){
    for(row in 2:nrow(g_info)){
      grab <- paste(g_info[row,"Comparison"], " ", g_info[row, "origin",], ", Sign: ", sep="")
      
      if(g_info[row,"C3"]==1){grab <- paste(grab, "C3")}
      if(g_info[row,"W3"]==1){grab <- paste(grab, "W3")}
      if(g_info[row,"C24"]==1){grab <- paste(grab, "C24")}
      if(g_info[row,"W24"]==1){grab <- paste(grab, "W24")}
      
      info <- paste(info, "\n", grab)
      
    }
  }
  # would like to have LogFoldChange and pvalue too, get from res df?
  
  # plot
  plot <- ggplot(Dat, aes(x=Hour_num, y=y, col=Treatment, shape=Treatment))+
    scale_colour_manual(values=c("dodgerblue2","firebrick3", "black"))+
    scale_shape_manual(values=c(19,17,18)) + # use different shapes to prevent overlap
    facet_wrap(~Treat_week)+
    # geom_point(data=ind.Dat, size=2.5, alpha=0.5)+
    geom_jitter(data=ind.Dat, size=2.5, alpha=0.5, width=0.3, height=0)+
    geom_point(size=4) +
    geom_line()+
    scale_x_continuous(breaks=c(1, 2, 6), labels=paste(levels(Dat$Hour), "h", sep=""))+
    scale_y_continuous(breaks=seq(-10,10, by=1))+
    labs(title= paste(g, filter(DEGs_annot, GeneID==g)$Prot_descr[1], sep=": "),
         subtitle = info, 
         x="Timepoint",
         y="GeneExpression z-scores")+
    theme(strip.text = element_text(size=8), axis.text = element_text(size=8))
  
  plot_name <- paste(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/",threshold, "/DEGs/",g,  sep=""), ".png", sep="")# name for file
  png(plot_name, width=250, height=150, units="mm", res=300)
  print(plot)
  dev.off()
  
}

sessionInfo()