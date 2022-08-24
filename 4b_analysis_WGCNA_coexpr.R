# Script by Lindner et al. 2021 Mol Ecol used as template ####

# Open R project in top folder

# load packages
library(tidyverse)
library(WGCNA);
options(stringsAsFactors = FALSE);  ## this setting is important, do not omit.
allowWGCNAThreads()
library(dendextend)
library(edgeR)
library(lme4)
library(lmerTest)
library(car)
library(cowplot)
theme_set(theme_cowplot()) #white background instead of grey -> don't load if want grey grid

# WGCNA with RNAseq data: work on library size normalized counts and 
#                         do variance-stabilizing transformation (e.g. log-transform them using log2(x+1))

# -----------------------------------------------------------------------------
### 1. Prep dataset -----------------------------------------------------------
#------------------------------------------------------------------------------

# Load data ####
load("analysis/_RNAseq/_results/limma_mod_allweeks_filtn2.RData") # all weeks analyzed together

head(dfilt$counts) # raw counts filtered by cutoff ~20/median(library size) in at least 2 replicates
head(dfilt$samples) # sample information incl. library sizes, normalisation factors

# use expression values calculated by voom in DE analysis
head(v.24h$E) # numeric matrix of normalized expression values on the log2 scale for 24h samples
head(v.3h$E) # and for 3h samples

v <- cbind(v.3h$E, v.24h$E[,grepl("BEF", colnames(v.24h$E))==F]) # recombine, BEF samples only once!

rm(res_list.24h, res_list.3h, results.24h, results.3h, v.24h, v.3h) # clean up the rest


# Transpose counts #### 
# so that genes are columns and rows are samples
datExpr <- t(v)
dim(datExpr)
rownames(datExpr)
colnames(datExpr)
datExpr[1,c(1:20)]

datExpr <- datExpr[grepl("N_", rownames(datExpr))==T | grepl("BEF", rownames(datExpr))==T, ] # only look at constant 10C samples
dim(datExpr) # 36 samples


# Check if contains genes or samples with too many missing values ####
gsg <- goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK # if returns FALSE, remove the offending genes/samples from the data


# Cluster samples to look for outliers ####
sampleTree <- hclust(dist(datExpr), method = "average");
labels2colors(sampleTree$labels)

# Plot the sample tree
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)#, col.axis=labels2colors(sampleTree$labels))


# Convert traits to a color representation: white means low, red means high, grey means missing entry
sampleinfo <- filter(dfilt$samples, Sample %in% rownames(datExpr))
head(sampleinfo)

sampleinfo$Timenum <- as.numeric(ifelse(sampleinfo$Timepoint=="BEF", 0, ifelse(sampleinfo$Timepoint=="3", 3, 24)))
numColors <- numbers2colors(sampleinfo[,c("dev_median", "Timenum")])
colnames(numColors) <-c("dev_median", "Timenum")

infoColors <- as.data.frame(labels2colors(sampleinfo[,c("Treat_week")]))#, "Site", "Tree", "NovemberDate")])
colnames(infoColors) <-c("Treat_week") #, "Site", "Tree", "NovemberDate")
plotColors1 <- cbind(numColors, infoColors)
plotColors1 <- plotColors1[,c("Treat_week", "dev_median", "Timenum")]#, "Site", "Tree", "NovemberDate")] #change order
rm(infoColors, numColors)

# Plot the sample dendrogram with sample info colors underneath
# pdf(file = "analysis/_RNAseq/_plots/descriptive/WGCNA/sampleClustering.pdf", width = 12, height = 9)
plotDendroAndColors(sampleTree, plotColors1,
                    groupLabels = names(plotColors1),
                    main = "Sample dendrogram and trait heatmap")
# dev.off()

# No obvious outliers; samples not perfectly seggregated by Treat_week; some mixing between w4 and w6; and w6 and w8
# Data now ready for analysis


# -----------------------------------------------------------------------------
### 2. Call network -----------------------------------------------------------
#------------------------------------------------------------------------------

# Choose a set of soft-thresholding powers ####
powers <- c(c(1:10), seq(from = 12, to=30, by=2))

sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed") # other RNAseq analyses also used signed network
save(sft, file="analysis/_RNAseq/_results/WGCNA_sft_CoExpr10C.RData") ## time consuming --> save result to avoid running it again!
# load(file="analysis/_RNAseq/_results/WGCNA_sft_CoExpr10C.RData")

# access result
sft$fitIndices

# plot result: (1) scale free topology and (2) mean connectivity as a function of the soft-thresholding power
#sizeGrWindow(9, 5)
par(mar = c(5,6,3,2))
#par(mfrow = c(1,2));
#cex1 <- 1

# pdf(file = "analysis/_RNAseq/_plots/descriptive/WGCNA/power_threshold.pdf", width = 12, height = 9)
# (1)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab=expression(paste("Scale Free Topology, R "^"2")),type="n",
     main = "", cex.lab=2, cex.axis=1.8);
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1.5,col="red");
abline(h=0.91,col="red") ## line corresponds to using an R^2 cut-off of h, change if needed
# (2)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = "", cex.lab=2, cex.axis=1.8)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=1.5,col="red") 
# dev.off()


# Module making ####
# Here use soft-thresholding power=12 (lowest power for which the scale-free topology fit index reaches 0.90, see plots above)
# Run on server so can do all genes at once
net <- blockwiseModules(datExpr, power = 12, networkType = "signed",
                        TOMType = "signed", minModuleSize = 30, maxBlockSize = 23961, ## maxBlockSize is the number of genes in your dataset
                        reassignThreshold = 0, mergeCutHeight = 0.15, ## mergecutheight defines threshold for mergeing modules, 0.15=default
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "analysis/_RNAseq/_results/WGCNA_datExprAfterFilterTOM0.15",
                        verbose = 3)
save(net, file = "analysis/_RNAseq/_results/WGCNA_net_CoExpr10C.RData") ## time consuming --> save result to avoid running it again!
# load(file="analysis/_RNAseq/_results/WGCNA_net_CoExpr10C.RData")

# Decide on cutHeight ####
# Convert labels to colors for plotting
moduleColors <- as.data.frame(net$colors)
colnames(moduleColors) <- "Module"
moduleColors$colors <- labels2colors(moduleColors)

length(net$colors[moduleColors==0]) ## here, 3517 genes without module
unique(net$colors) %>% length # 29 modules
unique(moduleColors) # turquoise=0

# calculate module eigensites (or use net$MEs (net_update$newMEs if 'mergeCloseModules' used) for module eigensites)
MEList1 <- moduleEigengenes(datExpr, colors = moduleColors$colors, nPC=10) # calculate % of variance explained up to PC10

head(MEList1$varExplained) # var explained by PC1-10
dim(MEList1$varExplained) # cols=modules, rows = PCs

# isolate module eigensites and calculate dissimilarity
MEs1 <- MEList1$eigengenes # 1 per sample per module
MEDiss <- 1-cor(MEs1);
METree <- hclust(as.dist(MEDiss), method = "average"); ## cluster module eigengenes

sizeGrWindow(7, 6)
# pdf(file = "analysis/_RNAseq/_plots/descriptive/WGCNA/MEDendrogram_CH0.15.pdf", width = 12, height = 9)
plot(METree, main = "Clustering of module eigengenes", ## plot clustering of module eigengenes
     xlab = "", sub = "")
abline(h=0.4, col="red") # cutHeight we want?
# dev.off()

# if you want to increase the mergeCutHeight, update network using: ####
net_update <- mergeCloseModules(datExpr, net$colors,
                                net$MEs, cutHeight = 0.4, ## insert new cutHeight here!!
                                verbose = 3)
save(net_update, file = "analysis/_RNAseq/_results/WGCNA_net_CoExpr10C_0.4.RData") ## time consuming --> save result to avoid running it again!
# load(file="analysis/_RNAseq/_results/WGCNA_net_CoExpr10C_0.4.RData")


# Check module clustering ####
# Convert labels to colors for plotting
mergedColors <- as.data.frame(net_update$colors)
colnames(mergedColors) <- "Module"
mergedColors$colors <- labels2colors(mergedColors)
length(net_update$colors[mergedColors==0])
unique(mergedColors$colors) %>% length # 14 modules left

# calculate module eigensites (or use net$MEs (net_update$newMEs if 'mergeCloseModules' used) for module eigensites)
MEList2 <- moduleEigengenes(datExpr, colors = mergedColors$colors, nPC=10) # calculate % of variance explained up to PC10

head(MEList2$varExplained) # var explained by PC1-10
dim(MEList2$varExplained) # cols=modules, rows = PCs

# isolate module eigensites and calculate dissimilarity
MEs2 <- MEList2$eigengenes # 1 per sample per module
MEDiss <- 1-cor(MEs2);
METree <- hclust(as.dist(MEDiss), method = "average"); ## cluster module eigengenes

par(mar = c(5,6,3,2))
sizeGrWindow(7, 6)
# pdf(file = "analysis/_RNAseq/_plots/descriptive/WGCNA/MEDendrogram_CH0.4.pdf", width = 12, height = 9)
plot(METree, main = "Clustering of module eigengenes", ## plot clustering of module eigengenes
     xlab = "", sub = "", cex.main=2.5, cex.lab=2, cex.axis=1.8, cex=1.5)
#abline(h=0.55, col="red") # cutheight we want to try
# dev.off()


# FINAL CUTHEIGHT ####
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(orderMEs(MEs2), "", plotAdjacency=T, colorLabels = F, ## Plot the heatmap matrix (NOTE: this plot will overwrite the dendrogram plot!)
                      marDendro = c(0,4,2,4), marHeatmap = c(3,4,1,1.5),
                      plotDendrograms = T, xLabelsAngle = 90)

# get correlation between modules
cor(orderMEs(MEs2))


#------------------------------------------------------------------------------
### 3. Visualize network ------------------------------------------------------
#------------------------------------------------------------------------------
head(net$colors) #contains the module assignment per gene
head(net$MEs) #contains the module eigengenes of the modules = PC1; representative of gene expression profile
head(net_update$newMEs) 
head(net_update2$newMEs) 

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], moduleColors$colors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Plot the dendrogram and the module colors underneath (incl. for reduced model set)
# pdf(file = "analysis/_RNAseq/_plots/descriptive/WGCNA/GeneClusters_Dendrogram.pdf", width = 12, height = 9)
plotDendroAndColors(net$dendrograms[[1]], 
                    cbind(moduleColors$colors[net$blockGenes[[1]]], mergedColors$colors[net$blockGenes[[1]]], mergedColors2$colors[net$blockGenes[[1]]]),
                    c("Module colors", "ModCol merged 0.4", "ModCol merged 0.55"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# dev.off()


# get sites per module
SitesPerModule <- as.data.frame(table(net_update$colors))
colnames(SitesPerModule) <- c("Module", "NoSites")
SitesPerModule <- merge(SitesPerModule, unique(mergedColors), by="Module") # get same colors as used in dendrogram
SitesPerModule


#------------------------------------------------------------------------------
### 4. Relation between modules and Treatment Week ----------------------
#------------------------------------------------------------------------------

# Define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# Sampleinfo available
rownames(sampleinfo) <- sampleinfo$Sample
sampleinfo <- sampleinfo[rownames(MEs2),] # make sure Samples are in the same order in both dataframes
colnames(sampleinfo) # Trw_num = Treat_week numeric
str(sampleinfo)

datTrait <- MEs2
datTrait$Sample <- rownames(MEs2)
datTrait <- merge(datTrait, sampleinfo[,c("Sample", "Treat_week", "Trw_num", "dev_median", "Tube", "Tube.n", "Timepoint", "Timenum")], by=c("Sample"))
head(datTrait)

modNames <- data.frame(ModuleName=names(MEs2))
modNames <- modNames %>% mutate(Module=ifelse(ModuleName=="MEturquoise", "0", "Mod1")) %>% # turquoise is group of genes without Module
        arrange(Module)
modNames$Module[2:nrow(modNames)] <- paste("Mod", rep(1:(nrow(modNames)-1), each=1), sep="")
head(modNames)

# Plot module eigen genes
#------------------------------------------------------------------------------
for(mod in 1:nrow(modNames)){
        module <- modNames$ModuleName[mod]
        dat <- datTrait[,c(module, "Sample", "Treat_week", "Trw_num", "dev_median", "Tube", "Tube.n", "Timepoint", "Timenum")]
        colnames(dat)[1] <- "Module"
        color <- gsub("ME(\\w+)","\\1",module)
        plot <- ggplot(data=dat, aes(x=Trw_num, y=Module))+
                scale_size_manual(values=c(2,3,4))+
                geom_jitter(aes(size=Tube.n), shape=21, fill=color, alpha=0.7, height=0, width=0.3)+ #size=3 
                #geom_smooth()+
                labs(x="Treatment week", y="Module PC1", title=modNames$Module[mod]) +
                scale_y_continuous(breaks=seq(-1,1, by=0.1))+
                guides(shape = "none")+
                theme(legend.title= element_blank(), legend.text=element_text(size=18), axis.title = element_text(size=20), 
                      axis.text=element_text(size=16), plot.title = element_text(size=18))
        
        plot_name <- paste("_plots/WGCNA/modules/", modNames$Module[mod], "_", color, ".png", sep="")# name for file
        png(plot_name, width=250, height=150, units="mm", res=300)
        print(plot)
        dev.off()
        
}
# Save objects needed for plotting
modNames1 <- modNames
rownames(modNames1) <- modNames1$ModuleName

datTrait_plot <- datTrait
colnames(datTrait_plot)[c(2:15)] <- paste(colnames(datTrait)[c(2:15)], modNames1[colnames(datTrait)[c(2:15)], "Module"], sep="_") 
head(datTrait_plot) # module eigen genes per sample with color and Module number as colnames

# write.csv(datTrait_plot, file="analysis/_RNAseq/_results/WGCNA_0.4_MEs_to_plot.csv", row.names=F)
rm(modNames1, datTrait_plot)


# LMM to Test for module relation to TreatmentWeek ####
#------------------------------------------------------------------------------
threshold <- 0.05 # threshold to use for testing module significance and module membership significance ####

head(MEs2) # Module eigen genes = data to do tests on

head(modNames)
for (modul in 1:nrow(modNames)){
        df <- as.data.frame(MEs2[,modNames$ModuleName[modul]])
        colnames(df) <- "MEs" #colnames(MEs2)[modul]
        df$Sample <- rownames(MEs2)
        df <- merge(df, sampleinfo[,c("Sample", "Treat_week", "Trw_num", "dev_median", "Timepoint", "Tube")], by="Sample")
        head(df)
        
        m <- lmer(MEs ~ Treat_week + (1|Tube), data=df)
        test <- anova(m)
        modNames$FTreatWeek[modul] <- test$`F value`
        modNames$SignTreatWeek[modul] <- test$`Pr(>F)`
}

modNames$FDRTreatWeek <- p.adjust(modNames$SignTreatWeek, method = "BH") #use BH to correct for multiple testing
modNames$signMod <- as.numeric(modNames$FDRTreatWeek<threshold)
table(modNames$signMod) # four significant modules
head(modNames)


# Visualize Treat_week significance per module ####
#------------------------------------------------------------------------------
# prepare data for plotting
PlotSign <- as.matrix(modNames$signMod)
row.names(PlotSign) <- modNames$Module

# plot correlation between modules and reproductive state
CorrPlot <- ggcorrplot::ggcorrplot(PlotSign, p.mat=as.matrix(modNames[,"FDRTreatWeek"]), 
                                   method = "circle", colors=c("royalblue", "white", "gold1"), 
                                   legend.title="F-value", show.legend = FALSE) +
        labs(title="Modules relation to Treatment Week")+
        theme(legend.text=element_text(size=13),
              legend.title=element_text(size=16),
              axis.text.x=element_text(size=14, color="black"),
              axis.text.y=element_text(size=14, color="black"),)+
        scale_y_discrete(labels=c("Treatment Week"))
png("_plots/WGCNA/moduleSign.png", width=250, height=150, units="mm", res=300)
print(CorrPlot)
dev.off()


#--------------------------------------------------------------------
### 4. Per gene correlations with module ----------------------------
#--------------------------------------------------------------------

# Create the starting data frame
SiteInfo0 <- data.frame(gene = colnames(datExpr),
                        moduleColor = mergedColors$colors) #,
# SSTrait,
# SSPvalue)
# SiteInfo0$padj.SS_Treatweek <- p.adjust(SiteInfo0$p.SS_Treatweek, method = "BH") # BH correction for multiple testing
head(SiteInfo0)

# Get correlation of all genes with each module eigensite ---> gene module membership (MM)
SiteMM <- as.data.frame(cor(datExpr, MEs2, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(SiteMM), nSamples))
MMPadj <- sapply(MMPvalue, p.adjust, method = "BH") # correct for multiple testing with BH

# Add module membership information into one dataframe
head(modNames)
SiteInfo <- SiteInfo0
for (mod in 1:ncol(SiteMM)) {
        module <- modNames$ModuleName[mod]
        oldNames <- names(SiteInfo)
        SiteInfo <- data.frame(SiteInfo, SiteMM[, module], 
                               MMPadj[, module]);
        names(SiteInfo) = c(oldNames, paste("MM_", substring(module, 3), sep=""),
                            paste("padj.MM_", substring(module, 3), sep=""))
}
head(SiteInfo)


modNames$moduleColor <- gsub("ME(\\w+)", "\\1", modNames$ModuleName)
SiteInfo <- merge(SiteInfo, modNames[,c("moduleColor", "Module")], by="moduleColor") # add Module number

# Get list of genes belonging to significant modules
signMod <- gsub("ME(\\w+)","\\1", filter(modNames, signMod==1)$ModuleName)

SiteInfo.sign <-  cbind(SiteInfo[,c("gene", "Module","moduleColor")],
                        SiteInfo[,grepl("_brown",colnames(SiteInfo))==T],
                        SiteInfo[,grepl("_pink",colnames(SiteInfo))==T], 
                        SiteInfo[,grepl("_red",colnames(SiteInfo))==T], 
                        SiteInfo[,grepl("_yellow",colnames(SiteInfo))==T])

SiteInfo.sign <- SiteInfo.sign %>% filter(moduleColor %in% c(signMod)) %>%
        filter((moduleColor=="brown" & padj.MM_brown<threshold) | 
                       (moduleColor=="pink" & padj.MM_pink<threshold) | 
                       (moduleColor=="red" & padj.MM_red<threshold) |
                       (moduleColor=="yellow" & padj.MM_yellow<threshold)) #%>% # only keep genes with significant module membership to these modules
# filter(padj.SS_Treatweek<0.05) # only keep genes with significant trait-based site significance (think this a bit too conservative for my interests)

head(SiteInfo.sign) 
table(SiteInfo.sign$Module)
# ~7445 genes that significantly belong to the modules that relate to Treat_week; 
# Could still run LMM for each gene, and only select those with sign Treat_week fixed effect
# But seems bit overkill? When module membership sign. already enough ####


# Save Results ####
#------------------------------------------------------------------------------
# save(modNames, SiteInfo, SiteInfo.sign, file="analysis/_RNAseq/_results/WGCNA_Results_CoExpr10C.RData")
# write.csv(SiteInfo, file="analysis/_RNAseq/_results/WGCNA_Results_CoExpr10C.csv", row.names=F)

load(file="analysis/_RNAseq/_results/WGCNA_Results_CoExpr10C.RData")
nrgenes <- table(SiteInfo$Module) %>% as.data.frame
colnames(nrgenes) <- c("Module", "NrGenes")
modNames <- merge(modNames, nrgenes, by="Module")
# write.csv(modNames, file="analysis/_RNAseq/_results/WGCNA_ANOVA-results.csv", row.names=F)

sessionInfo()

