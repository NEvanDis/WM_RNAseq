#-------------------------------------------------------
#-------------------------------------------------------
# Visualisation of RNAseq analysis results of egg samples from Wintermoth Transfer experiment
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
load("analysis/_RNAseq/_results/limma_mod_allweeks_filtn2.RData") # all weeks analyzed together

head(dfilt$counts) # raw counts filtered by cutoff ~20/median(library size) in at least 2 replicates
head(dfilt$samples) # sample information incl. library sizes, normalisation factors
lapply(res_list.3h, head) # results for 22 contrasts 3h model
lapply(res_list.24h, head) # results for 22 contrasts 24h model

res_lists <- list(res_list.3h, res_list.24h)
names(res_lists) <- c("3h", "24h")
rm(res_list.3h, res_list.24h)


# ANNOTATION ####
#-------------------------------------------------------
annot <- read.table("analysis/_RNAseq/_annotation/genesAnnot_wblastp_rows_exp.tsv", header=T)
str(annot)
head(annot[,c("GeneID", "TranscriptID", "Evalue")])# ordered by GeneID and Evalue, most significant hit for each GeneID first
length(unique(annot$GeneID)) # this should be 29113
table(is.na(annot$Evalue)) #4142 genes without hit


#-------------------------------------------------------
#-------------------------------------------------------
# Make list of significant genes ####
#-------------------------------------------------------
#-------------------------------------------------------
threshold <- 0.01 #FDR threshold I want to use ####

DEGs_lists <- list()

for(list in 1:length(res_lists)){
  res_list <- res_lists[[list]]
  names(res_list)

  # contrasts of interest
  names(res_list)[c(1:2)] <- c("temp_CvsN", "temp_WvsN")
  contr_within <- c(names(res_list)[c(3:10)])
  contr_between <- c(names(res_list)[c(11:22)])

  # Create lists to store results in
  sigGenes_all_down <- list()
  sigGenes_all_up <- list()
  sigGenes_all <- list()

  weeks <- unique(gsub("(w\\d_)\\wvs\\w", "\\1", names(res_list[contr_within])))
  week_comp <- unique(gsub("\\w.(w\\dvs\\d)", "\\1", names(res_list[contr_between])))
  contrasts <- c(weeks, week_comp, "temp_")

  mod <- names(res_lists)[list]

  for(contrast in 1:length(contrasts)){
    origin <- ifelse(grepl("w\\d_", contrasts[contrast])==T, gsub("(w\\d)_", "\\1", contrasts[contrast]),
                     ifelse(grepl("temp_", contrasts[contrast])==T, "Overall_Temp",contrasts[contrast]))

    contr <- names(res_list)[grepl(contrasts[contrast], names(res_list))==T] # should be 2 contrasts

    res <- res_list[contr]
    lapply(res, head)

    # pull contrasts of interest
    ResC <- res[[contr[1]]]
    ResW <- res[[contr[2]]]

    head(ResC)
    summary(ResC)

    # pull significant down-regulated genes from contrast
    sigGenesC <- subset(ResC, padj<threshold & log2FoldChange<0)
    sigGenesW <- subset(ResW, padj<threshold & log2FoldChange<0)

    # combine sigGenes from all contrast, remove duplicates and denote from which contrast each gene comes
    sigGenes_down <- c(rownames(sigGenesC), rownames(sigGenesW))
    sigGenes_down <- as.data.frame(unique(sigGenes_down))
    colnames(sigGenes_down) <- "GeneID"
    sigGenes_down <- as.data.frame(unique(sigGenes_down)) %>% mutate(C=ifelse(sigGenes_down$GeneID %in% rownames(sigGenesC), 1, 0),
                                                                     W=ifelse(sigGenes_down$GeneID %in% rownames(sigGenesW), 1, 0),
                                                                     origin=origin, EffDir="down")
    head(sigGenes_down)

    sigGenes_all_down[[contrast]] <- sigGenes_down
    names(sigGenes_all_down)[contrast] <- contrasts[contrast]


    # pull significant up-regulated genes from contrast
    sigGenesC <- subset(ResC, padj<threshold & log2FoldChange>0)
    sigGenesW <- subset(ResW, padj<threshold & log2FoldChange>0)

    # combine sigGenes from all contrast, remove duplicates and denote from which contrast each gene comes
    sigGenes_up <- c(rownames(sigGenesC), rownames(sigGenesW))
    sigGenes_up <- as.data.frame(unique(sigGenes_up))
    colnames(sigGenes_up) <- "GeneID"
    sigGenes_up <- as.data.frame(unique(sigGenes_up)) %>% mutate(C=ifelse(sigGenes_up$GeneID %in% rownames(sigGenesC), 1, 0),
                                                                 W=ifelse(sigGenes_up$GeneID %in% rownames(sigGenesW), 1, 0),
                                                                 origin=origin, EffDir="up")
    sigGenes_all_up[[contrast]] <- sigGenes_up
    names(sigGenes_all_up)[contrast] <- contrasts[contrast]

    # list all sigGenes together
    #sigGenes <- rbind(sigGenes_down, sigGenes_up) # causes duplicates if in one contrast Gene is upregulated and in another down

    # pull all significant genes from contrast
    sigGenesC <- subset(ResC, padj<threshold)
    sigGenesW <- subset(ResW, padj<threshold)

    # combine sigGenes from all contrast, remove duplicates and denote from which contrast each gene comes
    sigGenes <- c(rownames(sigGenesC), rownames(sigGenesW))
    sigGenes <- as.data.frame(unique(sigGenes))
    colnames(sigGenes) <- "GeneID"
    sigGenes <- as.data.frame(unique(sigGenes)) %>% mutate(C=ifelse(sigGenes$GeneID %in% rownames(sigGenesC), 1, 0),
                                                           W=ifelse(sigGenes$GeneID %in% rownames(sigGenesW), 1, 0),
                                                           origin=origin)
    head(sigGenes)

    sigGenes_all[[contrast]] <- sigGenes
    names(sigGenes_all)[contrast] <- contrasts[contrast]

  }
  rm(res, sigGenesC, sigGenesW, ResC, ResW)

  lapply(sigGenes_all, head)
  lapply(sigGenes_all, nrow)

  sigGenes_down <- data.table::rbindlist(sigGenes_all_down)
  sigGenes_up <- data.table::rbindlist(sigGenes_all_up)
  sigGenes <- data.table::rbindlist(sigGenes_all)
  #sigGenes <- unique(subset(sigGenes, select=-c(origin)))

  DEGs <- list(sigGenes, sigGenes_down, sigGenes_up)
  names(DEGs) <- c("sigGenes", "sigGenes_down", "sigGenes_up")
  DEGs_lists[[list]] <- DEGs
  names(DEGs_lists)[list] <- mod
}
rm(sigGenes, sigGenes_down, sigGenes_up, sigGenes_all, sigGenes_all_down, sigGenes_all_up,
   contrast, weeks, res_list, contr, list, mod, DEGs, contr_between, contr_within,
   origin, week_comp)

lapply(DEGs_lists, head)
####

# Rename columns so can combine 3h and 24h results
lapply(DEGs_lists, head)
for(list in 1:length(DEGs_lists)){
  DEGs_list <- DEGs_lists[[list]]
  names(DEGs_list)

  sigGenes_down <- DEGs_list$sigGenes_down
  sigGenes_up <- DEGs_list$sigGenes_up
  sigGenes <- DEGs_list$sigGenes

  mod <- names(DEGs_lists)[list]

  if(mod=="3h"){
    colnames(sigGenes_down)[c(2,3)] <- c("C3", "W3")
    colnames(sigGenes_up)[c(2,3)] <- c("C3", "W3")
    colnames(sigGenes)[c(2,3)] <- c("C3", "W3")
  }else{
    colnames(sigGenes_down)[c(2,3)] <- c("C24", "W24")
    colnames(sigGenes_up)[c(2,3)] <- c("C24", "W24")
    colnames(sigGenes)[c(2,3)] <- c("C24", "W24")
  }

  DEGs <- list(sigGenes, sigGenes_down, sigGenes_up)
  names(DEGs) <- c("sigGenes", "sigGenes_down", "sigGenes_up")
  DEGs_lists[[list]] <- DEGs
  names(DEGs_lists)[list] <- mod
}
rm(list, DEGs_list, sigGenes_down, sigGenes_up, sigGenes, mod, DEGs)
lapply(DEGs_lists, head)

#-------------------------------------------------------
# SAVE sigGenes set ####
#-------------------------------------------------------

# contrasts
contr_all <- ifelse(grepl("w\\d_", contrasts)==T, gsub("(w\\d)_", "\\1", contrasts),
                    ifelse(grepl("temp_", contrasts)==T, "Overall_Temp",contrasts))
contr_within <- contr_all[1:4]
contr_between <- contr_all[5:10]
contr_temp <- contr_all[11]

# Combine results 3h and 24h model ####
DEGs <- merge(DEGs_lists$`3h`$sigGenes, DEGs_lists$`24h`$sigGenes, by=c("GeneID", "origin"), all=T)
DEGs <- DEGs %>% mutate(C3=ifelse(is.na(C3)==T, 0, C3), W3=ifelse(is.na(W3)==T, 0, W3),
                        C24=ifelse(is.na(C24)==T, 0, C24), W24=ifelse(is.na(W24)==T, 0, W24),
                        Comparison=ifelse(origin %in% contr_within, "Within_week",
                                          ifelse(origin %in% contr_between, "Between_week", origin))) %>%
  arrange(GeneID, origin, C3, W3, C24, W24)
head(DEGs)
length(unique(DEGs$GeneID))

DEGs_down <- merge(DEGs_lists$`3h`$sigGenes_down, DEGs_lists$`24h`$sigGenes_down, by=c("GeneID", "origin", "EffDir"), all=T)
DEGs_down <- DEGs_down %>% mutate(C3=ifelse(is.na(C3)==T, 0, C3), W3=ifelse(is.na(W3)==T, 0, W3),
                        C24=ifelse(is.na(C24)==T, 0, C24), W24=ifelse(is.na(W24)==T, 0, W24),
                        Comparison=ifelse(origin %in% contr_within, "Within_week",
                                          ifelse(origin %in% contr_between, "Between_week", origin))) %>%
  arrange(GeneID, origin, C3, W3, C24, W24)
head(DEGs_down)
length(unique(DEGs_down$GeneID))

DEGs_up <- merge(DEGs_lists$`3h`$sigGenes_up, DEGs_lists$`24h`$sigGenes_up, by=c("GeneID", "origin", "EffDir"), all=T)
DEGs_up <- DEGs_up %>% mutate(C3=ifelse(is.na(C3)==T, 0, C3), W3=ifelse(is.na(W3)==T, 0, W3),
                        C24=ifelse(is.na(C24)==T, 0, C24), W24=ifelse(is.na(W24)==T, 0, W24),
                        Comparison=ifelse(origin %in% contr_within, "Within_week",
                                          ifelse(origin %in% contr_between, "Between_week", origin))) %>%
  arrange(GeneID, origin, C3, W3, C24, W24)
head(DEGs_up)
length(unique(DEGs_up$GeneID))


# Add annotation and save sign genes ####
head(annot[,c("GeneID", "Read_frame", "GeneName", "Evalue", "Prot_descr", "Species")])

DEGs_all <- rbind(DEGs_down, DEGs_up)
length(unique(DEGs_all$GeneID)) # some genes double, because differs per model in which contrast they were significant

DEGs_annot <- merge(DEGs_all, annot[,c("GeneID", "Read_frame", "GeneName", "Evalue", "Prot_descr", "Species")],
                    by="GeneID", all.x=T)
DEGs_annot <- DEGs_annot[!duplicated(DEGs_annot[,c(1:8,12)]),] %>% arrange(Comparison, GeneID, origin, C3, W3, C24, W24, Evalue)
head(DEGs_annot)
# write.table(DEGs_annot, file=paste("analysis/_RNAseq/_results/DEGs_annot_n2_p", threshold,".tsv", sep=""), row.names=F, sep="\t")
# write.table(paste(unique(DEGs_annot$GeneID), ".", sep=""), file=paste("analysis/_RNAseq/_results/DEGsIDs_n2_p", threshold,".txt", sep=""), row.names=F, col.names=F)

table(is.na(annot$Evalue)) #4142 genes without hit in total
table(unique(DEGs_annot$GeneID) %in% annot[is.na(annot$Evalue),"GeneID"]) #84 DEGs without hit


#-------------------------------------------------------
#-------------------------------------------------------
# VISUALIZATION ####
#-------------------------------------------------------
#-------------------------------------------------------
DEGs <- read.table(file=paste("analysis/_RNAseq/_results/DEGs_annot_n2_p", threshold,".tsv", sep=""), header=T)
head(DEGs) # DEGs_annot

# Overlap 3h and 24h model
table(unique(filter(DEGs, C24==1 | W24==1)$GeneID) %in% unique(filter(DEGs, C3==1 | W3==1)$GeneID))


#-------------------------------------------------------
# VENN DIAGRAMS ####
#-------------------------------------------------------
length(unique(DEGs$GeneID))

# Overlap between contrasts ####
#-------------------------------------------------------
genes <- list(unique(filter(DEGs, Comparison=="Overall_Temp")$GeneID), 
              unique(filter(DEGs, Comparison=="Within_week")$GeneID),
              unique(filter(DEGs, Comparison=="Between_week")$GeneID))
names(genes) <- c("Overall_Temp", "Within_week", "Between_week")
lapply(genes, length)

venn_all <- ggVennDiagram::ggVennDiagram(genes, # needs to be a list, with each item = vector of genes names
                                         label_alpha = 0)+
  scale_fill_gradient(low="white",high = "mediumorchid4")
png(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold, "/Venn_DEG_contrasts.png", sep=""), width=250, height=150, units="mm", res=300)
print(venn_all)
dev.off()

shared_bw_contrasts <- filter(DEGs_annot, GeneID %in% filter(DEGs, Comparison=="Overall_Temp")$GeneID &
  GeneID %in% filter(DEGs, Comparison=="Within_week")$GeneID &
  GeneID %in% filter(DEGs, Comparison=="Between_week")$GeneID)
shared_bw_contrasts[,c("GeneID", "Prot_descr")] %>% unique


# Overall Temp effects ####
#-------------------------------------------------------
length(unique(filter(DEGs, Comparison=="Overall_Temp")$GeneID))

# all genes
sigGenes_all2 <- list(unique(filter(DEGs, origin=="Overall_Temp" & C3==1)$GeneID), unique(filter(DEGs, origin=="Overall_Temp" & C24==1)$GeneID), 
                      unique(filter(DEGs, origin=="Overall_Temp" & W24==1)$GeneID), unique(filter(DEGs, origin=="Overall_Temp" & W3==1)$GeneID))
names(sigGenes_all2) <- c("Cold 3h", "Cold 24h", "Warm 24h", "Warm 3h")
lapply(sigGenes_all2, length)

venn_all <- ggVennDiagram::ggVennDiagram(sigGenes_all2, # needs to be a list, with each item = vector of genes names
                                         label_alpha = 0)+
  scale_fill_gradient(low="white",high = "grey") + 
  scale_colour_manual(values=c("dodgerblue2", "dodgerblue4", "firebrick4", "firebrick3"))
png(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold, "/Venn_DEG_OverallTemp_all.png", sep=""), width=250, height=150, units="mm", res=300)
print(venn_all)
dev.off()


# Within week comparisons ####
#-------------------------------------------------------
length(unique(filter(DEGs, Comparison=="Within_week")$GeneID))

# all genes
sigGenes_all2 <- list(unique(filter(DEGs, origin=="w2")$GeneID), unique(filter(DEGs, origin=="w4")$GeneID), 
                      unique(filter(DEGs, origin=="w6")$GeneID), unique(filter(DEGs, origin=="w8")$GeneID))
names(sigGenes_all2) <- c("Week2", "Week4", "Week6", "Week8")
lapply(sigGenes_all2, length)

venn_all <- ggVennDiagram::ggVennDiagram(sigGenes_all2, # needs to be a list, with each item = vector of genes names
                                         label_alpha = 0)+
  scale_fill_gradient(low="white",high = "mediumorchid4")
png(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold, "/Venn_DEG_weeks_all.png", sep=""), width=250, height=150, units="mm", res=300)
print(venn_all)
dev.off()
# majority of genes development week specific, only 695 from 4875 differentially expressed in all developmental stages

# down regulated genes
sigGenes_all_down2 <- list(unique(filter(DEGs_down, origin=="w2")$GeneID), unique(filter(DEGs_down, origin=="w4")$GeneID), 
                           unique(filter(DEGs_down, origin=="w6")$GeneID), unique(filter(DEGs_down, origin=="w8")$GeneID))
names(sigGenes_all_down2) <- c("w2", "w4", "w6", "w8")
lapply(sigGenes_all_down2, length)

png(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold, "/Venn_DEG_weeks_down.png", sep=""), width=250, height=150, units="mm", res=300)
venn_down <- ggVennDiagram::ggVennDiagram(sigGenes_all_down2, # needs to be a list, with each item = vector of genes names
                                          label_alpha = 0) + 
  scale_fill_gradient(low="white",high = "royalblue3")
print(venn_down)
dev.off()

# up regulated genes
sigGenes_all_up2 <- list(unique(filter(DEGs_up, origin=="w2")$GeneID), unique(filter(DEGs_up, origin=="w4")$GeneID), 
                         unique(filter(DEGs_up, origin=="w6")$GeneID), unique(filter(DEGs_up, origin=="w8")$GeneID))
names(sigGenes_all_up2) <- c("w2", "w4", "w6", "w8")
lapply(sigGenes_all_up2, length)

png(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold, "/Venn_DEG_weeks_up.png", sep=""), width=250, height=150, units="mm", res=300)
venn_up <- ggVennDiagram::ggVennDiagram(sigGenes_all_up2, # needs to be a list, with each item = vector of genes names
                                        label_alpha = 0) + 
  scale_fill_gradient(low="white",high = "coral")
print(venn_up)
dev.off()

# if numbers don't add up, 
# means that for some genes that are sign DE in >1 week are upregulated in one week and down regulated in the other
rm(sigGenes_all2, sigGenes_all_up2, sigGenes_all_down2)


# Between week comparisons ####
#-------------------------------------------------------
length(unique(filter(DEGs, Comparison=="Between_week")$GeneID))

# Overlap between and within week comparisons
table(unique(filter(DEGs, Comparison=="Between_week")$GeneID) %in% unique(filter(DEGs, Comparison=="Within_week")$GeneID))

# all genes
sigGenes_all2 <- list(unique(filter(DEGs_all, origin=="w2vs4")$GeneID), unique(filter(DEGs_all, origin=="w2vs6")$GeneID), unique(filter(DEGs_all, origin=="w2vs8")$GeneID), 
                     unique(filter(DEGs_all, origin=="w4vs6")$GeneID), unique(filter(DEGs_all, origin=="w4vs8")$GeneID), unique(filter(DEGs_all, origin=="w6vs8")$GeneID))
names(sigGenes_all2) <- c(contr_between)
lapply(sigGenes_all2, length)

venn_betw_all <- ggVennDiagram::ggVennDiagram(sigGenes_all2, # needs to be a list, with each item = vector of genes names
                                              label_alpha = 0)+
  scale_fill_gradient(low="white",high = "mediumorchid4")
png(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold, "/Venn_DEG_between_weeks.png", sep=""), width=250, height=150, units="mm", res=300)
print(venn_betw_all)
dev.off()


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

for(set in 1:length(sets)){
  
  # for after 3h
  #-------------------------------------------------------
  vst.3h <- vst[rownames(vst) %in% unique(filter(DEGs, Comparison==sets[set] & (C3==1 | W3==1))$GeneID),]
  vst.3h <- vst.3h[,c("BEF_136", "BEF_326", "BEF_473", # put in order so can substract BEF
                      "BEF_102", "BEF_219", "BEF_367", 
                      "BEF_128", "BEF_390", "BEF_471", 
                      "BEF_407", "BEF_411","BEF_94",
                      "N_3_136", "N_3_326", "N_3_473",
                      "N_3_102", "N_3_219", "N_3_367", 
                      "N_3_128", "N_3_390", "N_3_471", 
                      "N_3_407", "N_3_411", "N_3_94",
                      "C_3_136", "C_3_326", "C_3_473",
                      "C_3_102", "C_3_219", "C_3_367",
                      "C_3_128", "C_3_390", "C_3_471",
                      "C_3_407", "C_3_411", "C_3_94",
                      "W_3_136", "W_3_326", "W_3_473", 
                      "W_3_102", "W_3_219", "W_3_367",
                      "W_3_128", "W_3_390", "W_3_471",
                      "W_3_407", "W_3_411", "W_3_94")]
  vst.3h.relN <- vst.3h[,c(13:24)] - vst.3h[,c(1:12)]
  vst.3h.relC <- vst.3h[,c(25:36)] - vst.3h[,c(1:12)]
  vst.3h.relW <- vst.3h[,c(37:48)] - vst.3h[,c(1:12)]
  
  vst.3h.rel <- cbind(vst.3h.relN, vst.3h.relC, vst.3h.relW)
  rm(vst.3h.relN, vst.3h.relC, vst.3h.relW)
  head(vst.3h.rel) # counts relative to BEF
  
  vst.3h.relmean <- vst.rel[rownames(vst.rel) %in% unique(filter(DEGs, Comparison==sets[set] & (C3==1 | W3==1))$GeneID),]
  head(vst.3h.relmean)
  
  
  # Get z-scores
  z.3h <- t(scale(t(vst.3h.rel), center=TRUE, scale=TRUE)) # pick counts or counts relative to BEF
  range(z.3h)
  head(z.3h)
  
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
  colorder.3h <- subset(colorder, Timepoint=="3") %>% select(Sample, Treat_week, Treatment, Timepoint) # do not show BEF timepoint
  colorder.3h.mean <- subset(colorder, Timepoint=="3" & Treatment!="N") %>% select(WeekTim, Treat_week, Timepoint) %>% unique # to plot means
  
  annotation.3h <- colorder.3h %>% select(Treat_week) # without BEF
  rownames(annotation.3h) <- colorder.3h$Sample
  annotation.3h.mean <- colorder.3h.mean %>% select(Treat_week) # with BEF
  rownames(annotation.3h.mean) <- colorder.3h.mean$WeekTim
  
  # With pheatmap
  my.breaks <- seq(-3, 3, by=0.1) 
  my.colors <- c(colorRampPalette(colors = c("blue", "black"))(length(my.breaks)/2), colorRampPalette(colors = c("black", "yellow"))(length(my.breaks)/2))
  
  hm.3h <- pheatmap::pheatmap(z.3h[,c(colorder.3h$Sample)], 
                              color=my.colors,
                              cluster_rows=TRUE, show_rownames=FALSE,
                              cluster_cols=FALSE, annotation_col=annotation.3h,
                              clustering_distance_rows="correlation") #Pearson correlation
  plot_name <- paste(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Hmap_DEGs_3h_all_", sets[set], sep=""), ".png", sep="")# name for file
  png(plot_name, width=350, height=250, units="mm", res=300)
  print(hm.3h)
  dev.off()

  hm.3h.mean <- pheatmap::pheatmap(z.3h.mean[,c(colorder.3h.mean$WeekTim)], 
                              color=my.colors,
                              cluster_rows=TRUE, show_rownames=FALSE,
                              cluster_cols=FALSE, annotation_col=annotation.3h.mean,
                              clustering_distance_rows="correlation") #Pearson correlation
  hms_3h[[set]] <-  hm.3h.mean
  names(hms_3h)[set] <- sets[set]
  plot_name <- paste(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Hmap_DEGs_3h_means_", sets[set], sep=""), ".png", sep="")# name for file
  png(plot_name, width=350, height=250, units="mm", res=300)
  print(hm.3h.mean)
  dev.off()  
  
  
  # for after 24h 
  #-------------------------------------------------------
  vst.24h <- vst[rownames(vst) %in% unique(filter(DEGs, Comparison==sets[set] & (C24==1 | W24==1))$GeneID),]
  vst.24h <- vst.24h[,c("BEF_136", "BEF_326", "BEF_473", # put in order so can substract BEF
                        "BEF_102", "BEF_219", "BEF_367", 
                        "BEF_128", "BEF_390", "BEF_471", 
                        "BEF_407", "BEF_411","BEF_94",
                        "N_24_136", "N_24_326", "N_24_473",
                        "N_24_102", "N_24_219", "N_24_367", 
                        "N_24_128", "N_24_390", "N_24_471", 
                        "N_24_407", "N_24_411", "N_24_94",
                        "C_24_136", "C_24_326", "C_24_473",
                        "C_24_102", "C_24_219", "C_24_367",
                        "C_24_128", "C_24_390", "C_24_471",
                        "C_24_407", "C_24_411", "C_24_94",
                        "W_24_136", "W_24_326", "W_24_473", 
                        "W_24_102", "W_24_219", "W_24_367",
                        "W_24_128", "W_24_390", "W_24_471",
                        "W_24_407", "W_24_411", "W_24_94")]
  vst.24h.relN <- vst.24h[,c(13:24)] - vst.24h[,c(1:12)]
  vst.24h.relC <- vst.24h[,c(25:36)] - vst.24h[,c(1:12)]
  vst.24h.relW <- vst.24h[,c(37:48)] - vst.24h[,c(1:12)]
  
  vst.24h.rel <- cbind(vst.24h.relN, vst.24h.relC, vst.24h.relW)
  rm(vst.24h.relN, vst.24h.relC, vst.24h.relW)
  head(vst.24h.rel) # counts relative to BEF
  
  vst.24h.relmean <- vst.rel[rownames(vst.rel) %in% unique(filter(DEGs, Comparison==sets[set] & (C24==1 | W24==1))$GeneID),]
  head(vst.24h.relmean)
  
  
  # Get z-scores
  z.24h <- t(scale(t(vst.24h.rel), center=TRUE, scale=TRUE)) # pick counts or counts relative to BEF
  range(z.24h)
  head(z.24h)
  
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
  colorder.24h <- subset(colorder, Timepoint=="24") %>% select(Sample, Treat_week, Treatment, Timepoint) # do not show BEF timepoint
  colorder.24h.mean <- subset(colorder, Timepoint=="24" & Treatment!="N") %>% select(WeekTim, Treat_week, Timepoint) %>% unique # to plot means
  
  annotation.24h <- colorder.24h %>% select(Treat_week) # without BEF
  rownames(annotation.24h) <- colorder.24h$Sample
  annotation.24h.mean <- colorder.24h.mean %>% select(Treat_week) # with BEF
  rownames(annotation.24h.mean) <- colorder.24h.mean$WeekTim
  
  # With pheatmap
  my.breaks <- seq(-3, 3, by=0.1) 
  my.colors <- c(colorRampPalette(colors = c("blue", "black"))(length(my.breaks)/2), colorRampPalette(colors = c("black", "yellow"))(length(my.breaks)/2))
  
  hm.24h <- pheatmap::pheatmap(z.24h[,c(colorder.24h$Sample)], 
                               color=my.colors,
                               cluster_rows=TRUE, show_rownames=FALSE,
                               cluster_cols=FALSE, annotation_col=annotation.24h,
                               clustering_distance_rows="correlation") #Pearson correlation
  plot_name <- paste(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Hmap_DEGs_24h_all_", sets[set], sep=""), ".png", sep="")# name for file
  png(plot_name, width=350, height=250, units="mm", res=300)
  print(hm.24h)
  dev.off()
 
  hm.24h.mean <- pheatmap::pheatmap(z.24h.mean[,c(colorder.24h.mean$WeekTim)], 
                                   color=my.colors,
                                   cluster_rows=TRUE, show_rownames=FALSE,
                                   cluster_cols=FALSE, annotation_col=annotation.24h.mean,
                                   clustering_distance_rows="correlation") #Pearson correlation
  hms_24h[[set]] <-  hm.24h.mean
  names(hms_24h)[set] <- sets[set]
  plot_name <- paste(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Hmap_DEGs_24h_means_", sets[set], sep=""), ".png", sep="")# name for file
  png(plot_name, width=350, height=250, units="mm", res=300)
  print(hm.24h.mean)
  dev.off()  
   
}


# Get heat map clusters to do GO analysis on 6x ####
#-------------------------------------------------------
names(hms_3h) # stored heatmap trees
names(hms_24h)

clusters <- list()


# Find cut heights for 3h ####
## 1
dev.off() # clean plot window
hms_3h$Overall_Temp # hm

# save Tree dendrogram with where I cut
pdf(file = paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Overall_temp_3h_cut.pdf", sep=""), width = 12, height = 9)

plot(hms_3h$Overall_Temp$tree_row)
abline(h=1.75, col="red", lty=2, lwd=2)

dev.off()

cl.3h <- as.data.frame(sort(cutree(hms_3h$Overall_Temp$tree_row, h=1.75)))
names(cl.3h) <- "Cluster"
cl.3h$GeneID <- row.names(cl.3h)
cl.3h$Set <- "Overall_Temp"
cl.3h$Comparison <- "3h"
cl.3h$plotColors <- WGCNA::labels2colors(cl.3h$Cluster)
head(cl.3h) 
table(cl.3h$Cluster)
table(cl.3h$plotColors)

clusters[[1]] <- cl.3h # save clusters with hm they belong to

# Color + legend for cluster determination
plotColors <- as.data.frame(cl.3h$plotColors) # get colors
rownames(plotColors) <- cl.3h$GeneID

cols <- unique(cl.3h[,c("Cluster","plotColors")]) # get legend
name <- "Clusters: "
for(clust in 1:length(unique(cl.3h$Cluster))){
  name1 <- paste(clust, filter(cols, Cluster==clust)$plotColors, sep="=")
  name <- paste(name, name1, sep=" ")
}
name

# Plot colors in order of heat map labels
pdf(file = paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Overall_temp_3h_ClustCol.pdf", sep=""), width = 12, height = 9)
WGCNA::plotDendroAndColors(hms_3h$Overall_Temp$tree_row, plotColors[c(hms_3h$Overall_Temp$tree_row$labels),],
                    groupLabels = "Cluster",
                    main = name)
dev.off()


## 2
dev.off() # clean plot window
hms_3h$Within_week # hm

# save Tree dendrogram with where I cut
pdf(file = paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Within_week_3h_cut.pdf", sep=""), width = 12, height = 9)

plot(hms_3h$Within_week$tree_row)
abline(h=1.75, col="red", lty=2, lwd=2)

dev.off()

cl.3h <- as.data.frame(sort(cutree(hms_3h$Within_week$tree_row, h=1.75)))
names(cl.3h) <- "Cluster"
cl.3h$GeneID <- row.names(cl.3h)
cl.3h$Set <- "Within_week"
cl.3h$Comparison <- "3h"
cl.3h$plotColors <- WGCNA::labels2colors(cl.3h$Cluster)
head(cl.3h) 
table(cl.3h$Cluster)
table(cl.3h$plotColors)

clusters[[2]] <- cl.3h # save clusters with hm they belong to

# Color + legend for cluster determination
plotColors <- as.data.frame(cl.3h$plotColors) # get colors
rownames(plotColors) <- cl.3h$GeneID

cols <- unique(cl.3h[,c("Cluster","plotColors")]) # get legend
name <- "Clusters: "
for(clust in 1:length(unique(cl.3h$Cluster))){
  name1 <- paste(clust, filter(cols, Cluster==clust)$plotColors, sep="=")
  name <- paste(name, name1, sep=" ")
}
name

# Plot colors in order of heat map labels
pdf(file = paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Within_week_3h_ClustCol.pdf", sep=""), width = 12, height = 9)
WGCNA::plotDendroAndColors(hms_3h$Within_week$tree_row, plotColors[c(hms_3h$Within_week$tree_row$labels),],
                           groupLabels = "Cluster",
                           main = name)
dev.off()


## 3
dev.off() # clean plot window
hms_3h$Between_week # hm

# save Tree dendrogram with where I cut
pdf(file = paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Between_week_3h_cut.pdf", sep=""), width = 12, height = 9)

plot(hms_3h$Between_week$tree_row)
abline(h=1.75, col="red", lty=2, lwd=2)

dev.off()

cl.3h <- as.data.frame(sort(cutree(hms_3h$Between_week$tree_row, h=1.75)))
names(cl.3h) <- "Cluster"
cl.3h$GeneID <- row.names(cl.3h)
cl.3h$Set <- "Between_week"
cl.3h$Comparison <- "3h"
cl.3h$plotColors <- WGCNA::labels2colors(cl.3h$Cluster)
head(cl.3h) 
table(cl.3h$Cluster)
table(cl.3h$plotColors)

clusters[[3]] <- cl.3h # save clusters with hm they belong to

# Color + legend for cluster determination
plotColors <- as.data.frame(cl.3h$plotColors) # get colors
rownames(plotColors) <- cl.3h$GeneID

cols <- unique(cl.3h[,c("Cluster","plotColors")]) # get legend
name <- "Clusters: "
for(clust in 1:length(unique(cl.3h$Cluster))){
  name1 <- paste(clust, filter(cols, Cluster==clust)$plotColors, sep="=")
  name <- paste(name, name1, sep=" ")
}
name

# Plot colors in order of heat map labels
pdf(file = paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Between_week_3h_ClustCol.pdf", sep=""), width = 12, height = 9)
WGCNA::plotDendroAndColors(hms_3h$Between_week$tree_row, plotColors[c(hms_3h$Between_week$tree_row$labels),],
                           groupLabels = "Cluster",
                           main = name)
dev.off()


# Find cut heights for 24h ####
## 1
dev.off() # clean plot window
hms_24h$Overall_Temp # hm

# save Tree dendrogram with where I cut
pdf(file = paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Overall_temp_24h_cut.pdf", sep=""), width = 12, height = 9)

plot(hms_24h$Overall_Temp$tree_row)
abline(h=1.75, col="red", lty=2, lwd=2)

dev.off()

cl.24h <- as.data.frame(sort(cutree(hms_24h$Overall_Temp$tree_row, h=1.75)))
names(cl.24h) <- "Cluster"
cl.24h$GeneID <- row.names(cl.24h)
cl.24h$Set <- "Overall_Temp"
cl.24h$mod <- "24h"
cl.24h$plotColors <- WGCNA::labels2colors(cl.24h$Cluster)
head(cl.24h) 
table(cl.24h$Cluster)
table(cl.24h$plotColors)

clusters[[4]] <- cl.24h # save clusters with hm they belong to

# Color + legend for cluster determination
plotColors <- as.data.frame(cl.24h$plotColors) # get colors
rownames(plotColors) <- cl.24h$GeneID

cols <- unique(cl.24h[,c("Cluster","plotColors")]) # get legend
name <- "Clusters: "
for(clust in 1:length(unique(cl.24h$Cluster))){
  name1 <- paste(clust, filter(cols, Cluster==clust)$plotColors, sep="=")
  name <- paste(name, name1, sep=" ")
}
name

# Plot colors in order of heat map labels
pdf(file = paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Overall_temp_24h_ClustCol.pdf", sep=""), width = 12, height = 9)
WGCNA::plotDendroAndColors(hms_24h$Overall_Temp$tree_row, plotColors[c(hms_24h$Overall_Temp$tree_row$labels),],
                           groupLabels = "Cluster",
                           main = name)
dev.off()


## 2
dev.off() # clean plot window
hms_24h$Within_week # hm

# save Tree dendrogram with where I cut
pdf(file = paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Within_week_24h_cut.pdf", sep=""), width = 12, height = 9)

plot(hms_24h$Within_week$tree_row)
abline(h=1.8, col="red", lty=2, lwd=2)

dev.off()

cl.24h <- as.data.frame(sort(cutree(hms_24h$Within_week$tree_row, h=1.8)))
names(cl.24h) <- "Cluster"
cl.24h$GeneID <- row.names(cl.24h)
cl.24h$Set <- "Within_week"
cl.24h$mod <- "24h"
cl.24h$plotColors <- WGCNA::labels2colors(cl.24h$Cluster)
head(cl.24h) 
table(cl.24h$Cluster)
table(cl.24h$plotColors)

clusters[[5]] <- cl.24h # save clusters with hm they belong to

# Color + legend for cluster determination
plotColors <- as.data.frame(cl.24h$plotColors) # get colors
rownames(plotColors) <- cl.24h$GeneID

cols <- unique(cl.24h[,c("Cluster","plotColors")]) # get legend
name <- "Clusters: "
for(clust in 1:length(unique(cl.24h$Cluster))){
  name1 <- paste(clust, filter(cols, Cluster==clust)$plotColors, sep="=")
  name <- paste(name, name1, sep=" ")
}
name

# Plot colors in order of heat map labels
pdf(file = paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Within_week_24h_ClustCol.pdf", sep=""), width = 12, height = 9)
WGCNA::plotDendroAndColors(hms_24h$Within_week$tree_row, plotColors[c(hms_24h$Within_week$tree_row$labels),],
                           groupLabels = "Cluster",
                           main = name)
dev.off()


## 3
dev.off() # clean plot window
hms_24h$Between_week # hm

# save Tree dendrogram with where I cut
pdf(file = paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Between_week_24h_cut.pdf", sep=""), width = 12, height = 9)

plot(hms_24h$Between_week$tree_row)
abline(h=1.75, col="red", lty=2, lwd=2)

dev.off()

cl.24h <- as.data.frame(sort(cutree(hms_24h$Between_week$tree_row, h=1.75)))
names(cl.24h) <- "Cluster"
cl.24h$GeneID <- row.names(cl.24h)
cl.24h$Set <- "Between_week"
cl.24h$mod <- "24h"
cl.24h$plotColors <- WGCNA::labels2colors(cl.24h$Cluster)
head(cl.24h) 
table(cl.24h$Cluster)
table(cl.24h$plotColors)

clusters[[6]] <- cl.24h # save clusters with hm they belong to

# Color + legend for cluster determination
plotColors <- as.data.frame(cl.24h$plotColors) # get colors
rownames(plotColors) <- cl.24h$GeneID

cols <- unique(cl.24h[,c("Cluster","plotColors")]) # get legend
name <- "Clusters: "
for(clust in 1:length(unique(cl.24h$Cluster))){
  name1 <- paste(clust, filter(cols, Cluster==clust)$plotColors, sep="=")
  name <- paste(name, name1, sep=" ")
}
name

# Plot colors in order of heat map labels
pdf(file = paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/hms/Between_week_24h_ClustCol.pdf", sep=""), width = 12, height = 9)
WGCNA::plotDendroAndColors(hms_24h$Between_week$tree_row, plotColors[c(hms_24h$Between_week$tree_row$labels),],
                           groupLabels = "Cluster",
                           main = name)
dev.off()


# Combine clusters into one dataframe ####
lapply(clusters, head)
clusters <- data.table::rbindlist(clusters)
table(clusters$Cluster, clusters$set, clusters$mod)
head(clusters)
# write.csv(clusters, file="analysis/_RNAseq/_results/DEGS_hm_clusters.csv", row.names=F)


# Add annotation ####
annotfull <- read.table("analysis/_RNAseq/_results/all_results_n2_p0.01.tsv", header=T) # table made in next script 5b ####
head(annotfull) #incl padj & log2FoldChange

cl.annot <- merge(clusters, annotfull, by=c("GeneID", "Set", "mod"))
cl.annot <- cl.annot[,c("Set", "mod", "Cluster", "plotColors", "GeneID", "Comparison", "contr",
                        "log2Mean", "log2FoldChange", "padj", "GeneName", "Prot_descr",    
                        "Evalue", "Read_frame", "Species")] # reorder
cl.annot <- cl.annot %>% arrange(Set, mod, Cluster, Comparison, contr, padj, desc(abs(log2FoldChange)),GeneID)
head(cl.annot)
# write.csv(cl.annot, file="analysis/_RNAseq/_results/DEGS_hm_clusters_annot.csv", row.names=F)


sessionInfo() %>% capture.output(file="analysis/_RNAseq/_src/env_visualization.txt")
