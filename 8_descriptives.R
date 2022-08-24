# Open R project in top folder

#-------------------------------------------------------
# Set up environment ####
#-------------------------------------------------------
# load packages
library(tidyverse)
library(PCAtools)
library(DESeq2)
library(edgeR) # loads limma as a dependency
library(circlize) #to generate a colour scale
library(pheatmap) # alternative to ComplexHeatmap
# library(ComplexHeatmap) # to make a heatmap

behandeling <- c("dodgerblue2","firebrick3","black") 
behandeling2 <- c("grey", "black", "dodgerblue2","firebrick3") 


# Load and prep data ####
#-------------------------------------------------------
annot <- read.table("analysis/_RNAseq/_annotation/genesAnnot_wblastp_rows_exp.tsv", header=T, sep="\t")
head(annot)

load("analysis/_RNAseq/_data/preprocessed2_v3cor_filt5.RData") # filt n2 count>5
counts <- countdata_unfilt # use unfiltered set, will filter below
dim(counts)
head(sampleinfo)
rm(countdata_unfilt, filtwhich.n2, filtwhich.n3)

load("analysis/_RNAseq/_results/limma_mod_allweeks_filtn2.RData") # all weeks analyzed together 
# contains used cut-off and filtered data set used for analysis
dim(dfilt$counts) # raw counts filtered by cutoff ~20/median(library size) in at least 2 replicates
str(dfilt$samples) # sample information incl. library sizes, normalisation factors

rm(res_list.24h, res_list.3h, results.24h, results.3h, v.24h, v.3h) # not needed here

# Make sure factor levels are coded correctly ####
levels(dfilt$samples$Treatment)
levels(dfilt$samples$Treatment2)
levels(dfilt$samples$Timepoint)
levels(dfilt$samples$Treat_week)
levels(dfilt$samples$Tube.n) # Hack individual effects nested in group


#-------------------------------------------------------
# Visualize descriptives ####
#-------------------------------------------------------

# Number of genes shared between weeks and/or treatments ####
#-------------------------------------------------------
filter <- cpm(counts)> cutoff # note for each clutch for each timepoint/treatment combo which genes expressed with at least >expression value (filter out lowly expressed genes)
head(filter)
filt <- filt <- matrix(c(rowCounts(filter[, c(colnames(filter)[grepl("BEF_",colnames(filter))==T &
                                                                 (gsub("BEF_(\\d+)","\\1",colnames(filter)) %in% c("136", "326", "473"))])]), # w2
                         rowCounts(filter[, c(colnames(filter)[grepl("C_3",colnames(filter))==T &
                                                                 (gsub("C_3_(\\d+)","\\1",colnames(filter)) %in% c("136", "326", "473"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("C_24",colnames(filter))==T &
                                                                 (gsub("C_24_(\\d+)","\\1",colnames(filter)) %in% c("136", "326", "473"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("N_3",colnames(filter))==T &
                                                                 (gsub("N_3_(\\d+)","\\1",colnames(filter)) %in% c("136", "326", "473"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("N_24",colnames(filter))==T &
                                                                 (gsub("N_24_(\\d+)","\\1",colnames(filter)) %in% c("136", "326", "473"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("W_3",colnames(filter))==T &
                                                                 (gsub("W_3_(\\d+)","\\1",colnames(filter)) %in% c("136", "326", "473"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("W_24",colnames(filter))==T &
                                                                 (gsub("W_24_(\\d+)","\\1",colnames(filter)) %in% c("136", "326", "473"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("BEF_",colnames(filter))==T &
                                                                 (gsub("BEF_(\\d+)","\\1",colnames(filter)) %in% c("102", "219", "367"))])]), # w4
                         rowCounts(filter[, c(colnames(filter)[grepl("C_3",colnames(filter))==T &
                                                                 (gsub("C_3_(\\d+)","\\1",colnames(filter)) %in% c("102", "219", "367"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("C_24",colnames(filter))==T &
                                                                 (gsub("C_24_(\\d+)","\\1",colnames(filter)) %in% c("102", "219", "367"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("N_3",colnames(filter))==T &
                                                                 (gsub("N_3_(\\d+)","\\1",colnames(filter)) %in% c("102", "219", "367"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("N_24",colnames(filter))==T &
                                                                 (gsub("N_24_(\\d+)","\\1",colnames(filter)) %in% c("102", "219", "367"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("W_3",colnames(filter))==T &
                                                                 (gsub("W_3_(\\d+)","\\1",colnames(filter)) %in% c("102", "219", "367"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("W_24",colnames(filter))==T &
                                                                 (gsub("W_24_(\\d+)","\\1",colnames(filter)) %in% c("102", "219", "367"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("BEF_",colnames(filter))==T &
                                                                 (gsub("BEF_(\\d+)","\\1",colnames(filter)) %in% c("128", "390", "471"))])]), # w6
                         rowCounts(filter[, c(colnames(filter)[grepl("C_3",colnames(filter))==T &
                                                                 (gsub("C_3_(\\d+)","\\1",colnames(filter)) %in% c("128", "390", "471"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("C_24",colnames(filter))==T &
                                                                 (gsub("C_24_(\\d+)","\\1",colnames(filter)) %in% c("128", "390", "471"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("N_3",colnames(filter))==T &
                                                                 (gsub("N_3_(\\d+)","\\1",colnames(filter)) %in% c("128", "390", "471"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("N_24",colnames(filter))==T &
                                                                 (gsub("N_24_(\\d+)","\\1",colnames(filter)) %in% c("128", "390", "471"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("W_3",colnames(filter))==T &
                                                                 (gsub("W_3_(\\d+)","\\1",colnames(filter)) %in% c("128", "390", "471"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("W_24",colnames(filter))==T &
                                                                 (gsub("W_24_(\\d+)","\\1",colnames(filter)) %in% c("128", "390", "471"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("BEF_",colnames(filter))==T &
                                                                 (gsub("BEF_(\\d+)","\\1",colnames(filter)) %in% c("94", "407", "411"))])]), # w8
                         rowCounts(filter[, c(colnames(filter)[grepl("C_3",colnames(filter))==T &
                                                                 (gsub("C_3_(\\d+)","\\1",colnames(filter)) %in% c("94", "407", "411"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("C_24",colnames(filter))==T &
                                                                 (gsub("C_24_(\\d+)","\\1",colnames(filter)) %in% c("94", "407", "411"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("N_3",colnames(filter))==T &
                                                                 (gsub("N_3_(\\d+)","\\1",colnames(filter)) %in% c("94", "407", "411"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("N_24",colnames(filter))==T &
                                                                 (gsub("N_24_(\\d+)","\\1",colnames(filter)) %in% c("94", "407", "411"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("W_3",colnames(filter))==T &
                                                                 (gsub("W_3_(\\d+)","\\1",colnames(filter)) %in% c("94", "407", "411"))])]),
                         rowCounts(filter[, c(colnames(filter)[grepl("W_24",colnames(filter))==T &
                                                                 (gsub("W_24_(\\d+)","\\1",colnames(filter)) %in% c("94", "407", "411"))])])), ncol=28)
colnames(filt) <- c("w2_BEF", "w2_C3", "w2_C24", "w2_N3", "w2_N24", "w2_W3", "w2_W24",
                    "w4_BEF", "w4_C3", "w4_C24", "w4_N3", "w4_N24", "w4_W3", "w4_W24",
                    "w6_BEF", "w6_C3", "w6_C24", "w6_N3", "w6_N24", "w6_W3", "w6_W24",
                    "w8_BEF", "w8_C3", "w8_C24", "w8_N3", "w8_N24", "w8_W3", "w8_W24")
rownames(filt) <- rownames(counts)
head(filt) # per timepoint/treatment combo, count number of times genes expressed >20 out of 3 replicates

filt.w2 <- filt[,colnames(filt)[grepl("w2", colnames(filt))==T]]
filt.w4 <- filt[,colnames(filt)[grepl("w4", colnames(filt))==T]]
filt.w6 <- filt[,colnames(filt)[grepl("w6", colnames(filt))==T]]
filt.w8 <- filt[,colnames(filt)[grepl("w8", colnames(filt))==T]]

table(rowCounts(filt.w2>=2)>=1) # only keep genes expressed in >=2 replicates at least in one timepoint/treatment combi

filtwhich.w2 <- row.names(filt.w2[rowCounts(filt.w2>=2)>=1,]) # get gene names to keep
filtwhich.w4 <- row.names(filt.w4[rowCounts(filt.w4>=2)>=1,]) # get gene names to keep
filtwhich.w6 <- row.names(filt.w6[rowCounts(filt.w6>=2)>=1,]) # get gene names to keep
filtwhich.w8 <- row.names(filt.w8[rowCounts(filt.w8>=2)>=1,]) # get gene names to keep

countsfiltw2 <- counts[rownames(counts) %in% filtwhich.w2,] # select these genes
countsfiltw4 <- counts[rownames(counts) %in% filtwhich.w4,] # select these genes
countsfiltw6 <- counts[rownames(counts) %in% filtwhich.w6,] # select these genes
countsfiltw8 <- counts[rownames(counts) %in% filtwhich.w8,] # select these genes


# Which genes week specific? ####
#-------------------------------------------------------
genesw468 <- unique(c(rownames(countsfiltw4), rownames(countsfiltw6), rownames(countsfiltw8)))
genesw268 <- unique(c(rownames(countsfiltw2), rownames(countsfiltw6), rownames(countsfiltw8)))
genesw248 <- unique(c(rownames(countsfiltw2), rownames(countsfiltw4), rownames(countsfiltw8)))
genesw246 <- unique(c(rownames(countsfiltw2), rownames(countsfiltw4), rownames(countsfiltw6)))

wkspecif.w2 <- countsfiltw2[!(row.names(countsfiltw2) %in% genesw468),] %>% rownames() %>% unique() %>% 
  as.data.frame() %>% mutate(week="w2")
wkspecif.w4 <- countsfiltw4[!(row.names(countsfiltw4) %in% genesw268),] %>% rownames() %>% unique()%>% 
  as.data.frame() %>% mutate(week="w4")
wkspecif.w6 <- countsfiltw6[!(row.names(countsfiltw6) %in% genesw248),] %>% rownames() %>% unique()%>% 
  as.data.frame() %>% mutate(week="w6")
wkspecif.w8 <- countsfiltw8[!(row.names(countsfiltw8) %in% genesw246),] %>% rownames() %>% unique()%>% 
  as.data.frame() %>% mutate(week="w8")

wkspecif <- rbind(wkspecif.w2, wkspecif.w4, wkspecif.w6, wkspecif.w8) # genes that meet the filter threshold on week level
colnames(wkspecif) <- c("GeneID", "wkspecif")
wkspecif <- merge(wkspecif, annot[,c("GeneID", "GeneName", "Prot_descr", "Evalue", "Read_frame", "Species")], by="GeneID")
wkspecif <- wkspecif %>% arrange(wkspecif, GeneID, Evalue, Prot_descr)
wkspecif <- wkspecif[!duplicated(wkspecif[,c(1,2,4)]),]
head(wkspecif)
# write.csv(wkspecif, file="analysis/_RNAseq/_results/weekspecific_genes_annot.csv", row.names=F)

# Venn diagram 
week_genes <- list(row.names(countsfiltw2), row.names(countsfiltw4), row.names(countsfiltw6), row.names(countsfiltw8))
names(week_genes) <- c("w2", "w4", "w6", "w8")
lapply(week_genes, length)

venn_week <- ggVennDiagram::ggVennDiagram(week_genes, # needs to be a list, with each item = vector of genes names
                                          label_alpha = 0)+
  scale_fill_gradient(low="white",high = "mediumorchid4")
png("analysis/_RNAseq/_plots/descriptive/Venn_Week_genesets.png", width=250, height=150, units="mm", res=300)
print(venn_week)
dev.off()


# Number of genes shared between temperatures ####
#-------------------------------------------------------
filt.10 <- filt[,colnames(filt)[grepl("_N", colnames(filt))==T | grepl("_BEF", colnames(filt))==T]]
filt.5 <- filt[,colnames(filt)[grepl("_C", colnames(filt))==T]]
filt.15 <- filt[,colnames(filt)[grepl("_W", colnames(filt))==T]]

table(rowCounts(filt.10>=2)>=1) # only keep genes expressed in >=2 replicates at least in one timepoint/treatment combi
table(rowCounts(filt.5>=2)>=1)
table(rowCounts(filt.15>=2)>=1)

filtwhich.10 <- row.names(filt.10[rowCounts(filt.10>=2)>=1,]) # get gene names to keep
filtwhich.5 <- row.names(filt.5[rowCounts(filt.5>=2)>=1,]) # get gene names to keep
filtwhich.15 <- row.names(filt.15[rowCounts(filt.15>=2)>=1,]) # get gene names to keep

ddsfilt10 <- counts[rownames(counts) %in% filtwhich.10,] # select these genes
ddsfilt5 <- counts[rownames(counts) %in% filtwhich.5,] # select these genes
ddsfilt15 <- counts[rownames(counts) %in% filtwhich.15,] # select these genes

# Venn diagram 
all_genes <- list(row.names(ddsfilt10), row.names(ddsfilt5), row.names(ddsfilt15))
names(all_genes) <- c("Ref.10C", "Cold.5C", "Warm.15C")
lapply(all_genes, length)

venn_temp <- ggVennDiagram::ggVennDiagram(all_genes, # needs to be a list, with each item = vector of genes names
                                          label_alpha = 0)+
  scale_fill_gradient(low="white",high = "grey")+
  scale_colour_manual(values=c("black", "dodgerblue2", "firebrick3"))
png("analysis/_RNAseq/_plots/descriptive/Venn_temperature.png", width=250, height=150, units="mm", res=300)
print(venn_temp)
dev.off()

# temp specific genes
tempspecif.10 <- filtwhich.10[!(filtwhich.10 %in% filtwhich.15)]
tempspecif.10 <- tempspecif.10[!(tempspecif.10 %in% filtwhich.5)]
tempspecif.10 <- data.frame(geneID=tempspecif.10, tempspecif="Control.10C")

tempspecif.5 <- filtwhich.5[!(filtwhich.5 %in% filtwhich.10)]
tempspecif.5 <- tempspecif.5[!(tempspecif.5 %in% filtwhich.15)]
tempspecif.5 <- data.frame(geneID=tempspecif.5, tempspecif="Cold.5C")

tempspecif.15 <- filtwhich.15[!(filtwhich.15 %in% filtwhich.10)]
tempspecif.15 <- tempspecif.15[!(tempspecif.15 %in% filtwhich.5)]
tempspecif.15 <- data.frame(geneID=tempspecif.15, tempspecif="Warm.15C")

tempspecif.shift <- filtwhich.5[filtwhich.5 %in% filtwhich.15 & !(filtwhich.5 %in% filtwhich.10)]
tempspecif.shift <- data.frame(geneID=tempspecif.shift, tempspecif="TempShift")

tempspecif <- rbind(tempspecif.10, tempspecif.5, tempspecif.15, tempspecif.shift)
rm(tempspecif.10, tempspecif.15, tempspecif.5, tempspecif.shift)

tempspecif <- merge(tempspecif, annot[,c("GeneID", "GeneName", "Prot_descr", "Evalue", "Read_frame", "Species")], by.x="geneID", by.y="GeneID")
tempspecif <- tempspecif %>% arrange(tempspecif, geneID, Evalue, Prot_descr)
tempspecif <- tempspecif[!duplicated(tempspecif[,c(1,2,4)]),]
head(tempspecif)
table(tempspecif$tempspecif)
# write.csv(tempspecif, file="analysis/_RNAseq/_results/tempspecific_genes_annot.csv", row.names=F)


# How many DEGs only expressed after temp change or week specific? How many of those not WGCNA overlap? ####
#-------------------------------------------------------------------------------------------------
DEGs_annot <- read.csv(file="analysis/_RNAseq/_results/DEGS_hm_clusters_annot.csv")
head(DEGs_annot)

WGCNAres <- read.csv("analysis/_RNAseq/_results/WGCNA_MMgenes-results.csv")
head(WGCNAres)

notOverlap <- filter(DEGs_annot, !(GeneID %in% filter(WGCNAres, SignMM=="Yes")$gene)) %>% select(GeneID) %>%  unique
head(notOverlap)

getGenes1 <- filter(DEGs_annot, GeneID %in% unique(tempspecif$geneID)) %>% # temp specific
  arrange(Set, mod, GeneID, Comparison, contr, log2FoldChange, padj, Evalue, Prot_descr)
unique(getGenes1$GeneID) %>% length
head(getGenes1)

getGenes2 <- filter(DEGs_annot, GeneID %in% unique(wkspecif$GeneID)) %>%
  arrange(Set, mod, GeneID, Comparison, contr, log2FoldChange, padj, Evalue, Prot_descr)
unique(getGenes2$GeneID) %>% length
head(getGenes2)

getGenes <- rbind(getGenes1, getGenes2) %>% unique
rm(getGenes1, getGenes2)
table(unique(getGenes2[,c("Set", "mod", "GeneID")])$Set)
getGenes <- merge(getGenes, unique(wkspecif[,c("GeneID", "wkspecif")]), all.x=T, by="GeneID")
getGenes <- merge(getGenes, unique(tempspecif[,c("geneID", "tempspecif")]), all.x=T, by.x="GeneID", by.y="geneID")
head(getGenes)

table(notOverlap$GeneID %in% all_genes$Ref.10C)
table(notOverlap$GeneID %in% unique(wkspecif$GeneID)) # ok also some week specific genes!!

getGenes$Overlap <- ifelse(getGenes$GeneID %in% notOverlap$GeneID, "No", "Yes")
getGenes <- getGenes %>% arrange(desc(Overlap), Set, desc(mod), log2FoldChange, padj, Evalue, Prot_descr) %>% unique
head(getGenes)
unique(getGenes$GeneID) %>% length
# write.table(getGenes, file="analysis/_RNAseq/_results/tempandwkspecific_genes_annot.tsv", row.names=F, sep="\t")

rm(wkspecif.w2, wkspecif.w4, wkspecif.w6, wkspecif.w8, genesw468, genesw268, genesw248, genesw246, 
   filtwhich.w2, filtwhich.w4, filtwhich.w6, filtwhich.w8, filter, filt, filt.w2, filt.w4, filt.w6, filt.w8,
   countsfiltw2, countsfiltw4, countsfiltw6, countsfiltw8)
####-------------------------------------------------------


# Can normalize together or not? ####
#-------------------------------------------------------
# Checked with PCA for full gene set of 26453 genes
#  and with Heatmaps for subset of DEGs from DE analysis

# Observations:
# - No outliers if control for sequencing depth
# - Big proportion of genes expressed in all weeks (see Venn diagram above)
# - Separate vst() per week doesn't look good
# - Don't see binary patterns when normalized together (or for raw data)

# Conclusion: Yes can normalize together ####

sessionInfo()