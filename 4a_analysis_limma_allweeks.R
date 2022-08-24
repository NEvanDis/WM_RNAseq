#-------------------------------------------------------
#-------------------------------------------------------
# RNAseq analysis of egg samples from Wintermoth Transfer experiment
# Following: Law CW, Alhamdoosh M, Su S et al. RNA-seq analysis is easy as 1-2-3 with 
#               limma, Glimma and edgeR [version 3; peer review: 3 approved] F1000Research 2018, 
#               5:1408 https://doi.org/10.12688/f1000research.9005.3
#-------------------------------------------------------
#-------------------------------------------------------


# Open R project in top folder

#-------------------------------------------------------
# Set up environment ####
#-------------------------------------------------------
# BiocManager::install("edgeR")

# load packages
library(tidyverse)
library(edgeR) # loads limma as a dependency
library(DESeq2) # needed for rowCounts


# Load and prep data ####
load("analysis/_RNAseq/_data/preprocessed2_v3cor_filt5.RData") # filt n2 count>5
counts <- countdata_unfilt # use unfiltered set, because will filter below
dim(counts)
head(sampleinfo)
rm(countdata_unfilt, filtwhich.n2, filtwhich.n3)

# Make sure factor levels are coded correctly
sampleinfo <- sampleinfo %>% mutate(Tube=as.factor(Tube), Treatment=as.factor(Treatment), Timepoint=as.factor(Timepoint), NovemberDate=as.factor(NovemberDate),
                                    Treat_week=as.factor(Treat_week), Trw_num= as.numeric(gsub("w(\\d)","\\1",Treat_week)), Tr_num=ifelse(Treatment=="N", 10, ifelse(Treatment=="W", 15, 5)),
                                    Treatment2=as.factor(ifelse(Timepoint=="BEF", "BEF", as.character(Treatment))))
sampleinfo <- sampleinfo %>% mutate(Treatment=factor(Treatment, levels=c("N", "C", "W")), Timepoint=factor(Timepoint, levels=c("BEF", "3", "24")),
                                    Treat_week=factor(Treat_week, levels=c("w2", "w4", "w6", "w8")),
                                    Treatment2=factor(Treatment2, levels=c("BEF","N","C", "W")))
levels(sampleinfo$Treatment)
levels(sampleinfo$Treatment2)
levels(sampleinfo$Timepoint)
levels(sampleinfo$Treat_week)
head(sampleinfo)
str(sampleinfo)

# Hack for individual effects nested in group
sampleinfo <- sampleinfo %>% arrange(Treat_week, Tube, Treatment, Timepoint)
sampleinfo$Tube.n <- factor(rep(rep(1:3, each=7),4)) # for 3 replicates in each group, 7 times same number, 4 times repeat
head(sampleinfo, n=28)
table(sampleinfo$Tube.n)
levels(sampleinfo$Tube.n)


# Create DGEList object
d <- DGEList(counts) # library sizes automatically calculated for each sample and normalisation factors are set to 1
d$samples$Sample <- row.names(d$samples)

sampleinfo <- merge(d$samples, sampleinfo, by="Sample")
head(sampleinfo)

d$samples <- sampleinfo # add sampleinfo 
head(d$samples)


# Prefilter ####
#-------------------------------------------------------
# keep genes with about ~20 read counts in at least 2 of 3 replicates on one timepoint
# Do actual filtering on CPM values to avoid giving preference to samples with large library sizes

cutoff <- 20/(median(d$samples$lib.size) * 1e-6) # 20 counts/median library size in million

filter <- cpm(d)> cutoff # note for each clutch for each timepoint/treatment combo which genes expressed with at least >5 (filter out lowly expressed genes)
head(filter)
filt <- matrix(c(rowCounts(filter[, c(colnames(filter)[grepl("BEF_",colnames(filter))==T &
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
rownames(filt) <- rownames(d)
head(filt) # per timepoint/treatment combo, count number of times genes expressed >20 out of 3 replicates

table(rowCounts(filt>=2)>=1) # only keep genes expressed in >=2 replicates at least in one timepoint/treatment combi
filtwhich.n2 <- rownames(filt[rowCounts(filt>=2)>=1,]) # get gene names to keep

dfilt <- d[rownames(d) %in% filtwhich.n2,] # select these genes
dim(d)
dim(dfilt) 
# cutoff ~10: 25569 genes expressed at least once in >=2 replicates at least in one timepoint/treatment combi
# cutoff ~20: 23961 genes expressed at least once in >=2 replicates at least in one timepoint/treatment combi


# Normalisation ####
dfilt <- calcNormFactors(dfilt, method="TMM")
head(dfilt$samples) # norm.factors updated

plotMDS(cpm(dfilt, log=TRUE)) # groups nicely according to mother

lcpm <- cpm(dfilt, log=TRUE)
boxplot(lcpm, las=2, main="")

#-------------------------------------------------------
# Model ####
#-------------------------------------------------------

# Design for Treatment effects within and between weeks, all weeks together
wm_fixed <- ~ -1 + Treat_week + Treatment2 + Treatment2:Treat_week
# wm_fixed <- ~ -1 + Treat_week + Treat_week:Tube.n + Treatment2 + Treatment2:Treat_week # with DESeq2 hack
m1 <- model.matrix(wm_fixed, data = sampleinfo)

# Use voom to remove variance dependency on mean
v <- voom(dfilt, design=m1, plot=TRUE)
dim(v)


# Fit model 3h ####
#-------------------------------------------------------
v.3h <- v[, -which(grepl("_24_",colnames(v)))]
m.3h <- model.matrix(wm_fixed, data = filter(sampleinfo, Timepoint!="24"))
dim(v.3h)

# Fit a random effect
corfit_main.3h <- duplicateCorrelation(v.3h, design=m.3h, block = filter(sampleinfo, Timepoint!="24")$Tube) # fitting a random effect?!
# Estimate the intra-block correlation given a block structure for the arrays or samples

# fit.3h <- lmFit(v.3h, design=m.3h) # with DESeq2 individual within group hack
fit.3h <- lmFit(v.3h, design=m.3h, block = filter(sampleinfo, Timepoint!="24")$Tube, correlation = corfit_main.3h$consensus)


# Fit model 24h ####
#-------------------------------------------------------
v.24h <- v[, -which(grepl("_3_",colnames(v)))]
m.24h <- model.matrix(wm_fixed, data = filter(sampleinfo, Timepoint!="3"))


# Fit a random effect
corfit_main.24h <- duplicateCorrelation(v.24h, design=m.24h, block = filter(sampleinfo, Timepoint!="3")$Tube) # fitting a random effect?!
# Estimate the intra-block correlation given a block structure for the arrays or samples

# fit.24h <- lmFit(v.24h, design=m.24h) # with DESeq2 individual within group hack
fit.24h <- lmFit(v.24h, design=m.24h, block = filter(sampleinfo, Timepoint!="3")$Tube, correlation = corfit_main.24h$consensus)


#-------------------------------------------------------
# Contrasts 
#-------------------------------------------------------

# for 3h ####
#-------------------------------------------------------
mod3h <- filter(sampleinfo, Timepoint!="24")

mod_mat <- fit.3h$design
mod_mat

# Model parameters to make contrasts
TreatN <- colMeans(mod_mat[mod3h$Treatment2 == "N",])
TreatC <- colMeans(mod_mat[mod3h$Treatment2 == "C",])
TreatW <- colMeans(mod_mat[mod3h$Treatment2 == "W",])

TreatN.w2 <- colMeans(mod_mat[mod3h$Treatment2 == "N" & mod3h$Treat_week=="w2", ])
TreatC.w2 <- colMeans(mod_mat[mod3h$Treatment2 == "C" & mod3h$Treat_week=="w2", ])
TreatW.w2 <- colMeans(mod_mat[mod3h$Treatment2 == "W" & mod3h$Treat_week=="w2", ])

TreatN.w4 <- colMeans(mod_mat[mod3h$Treatment2 == "N" & mod3h$Treat_week=="w4", ])
TreatC.w4 <- colMeans(mod_mat[mod3h$Treatment2 == "C" & mod3h$Treat_week=="w4", ])
TreatW.w4 <- colMeans(mod_mat[mod3h$Treatment2 == "W" & mod3h$Treat_week=="w4", ])

TreatN.w6 <- colMeans(mod_mat[mod3h$Treatment2 == "N" & mod3h$Treat_week=="w6", ])
TreatC.w6 <- colMeans(mod_mat[mod3h$Treatment2 == "C" & mod3h$Treat_week=="w6", ])
TreatW.w6 <- colMeans(mod_mat[mod3h$Treatment2 == "W" & mod3h$Treat_week=="w6", ])

TreatN.w8 <- colMeans(mod_mat[mod3h$Treatment2 == "N" & mod3h$Treat_week=="w8", ])
TreatC.w8 <- colMeans(mod_mat[mod3h$Treatment2 == "C" & mod3h$Treat_week=="w8", ])
TreatW.w8 <- colMeans(mod_mat[mod3h$Treatment2 == "W" & mod3h$Treat_week=="w8", ])


# Overal Treatment effect ####
CvsN <- as.matrix(TreatC - TreatN)
WvsN <- as.matrix(TreatW - TreatN)


# Within week ####
# N vs. C
w2_CvsN <- as.matrix(TreatC.w2 - TreatN.w2)
w4_CvsN <- as.matrix(TreatC.w4 - TreatN.w4)
w6_CvsN <- as.matrix(TreatC.w6 - TreatN.w6)
w8_CvsN <- as.matrix(TreatC.w8 - TreatN.w8)

# N vs. W
w2_WvsN <- as.matrix(TreatW.w2 - TreatN.w2)
w4_WvsN <- as.matrix(TreatW.w4 - TreatN.w4)
w6_WvsN <- as.matrix(TreatW.w6 - TreatN.w6)
w8_WvsN <- as.matrix(TreatW.w8 - TreatN.w8)


# Between weeks ####
# individual week comparisons
C.w2vs4 <- as.matrix((TreatC.w2 - TreatN.w2) - (TreatC.w4 - TreatN.w4)) # w2 (N vs. C)  vs. w4 (N vs. C)
C.w2vs6 <- as.matrix((TreatC.w2 - TreatN.w2) - (TreatC.w6 - TreatN.w6)) # w2 (N vs. C)  vs. w6 (N vs. C)
C.w2vs8 <- as.matrix((TreatC.w2 - TreatN.w2) - (TreatC.w8 - TreatN.w8)) # w2 (N vs. C)  vs. w8 (N vs. C)
C.w4vs6 <- as.matrix((TreatC.w4 - TreatN.w4) - (TreatC.w6 - TreatN.w6)) # w4 (N vs. C)  vs. w6 (N vs. C)
C.w4vs8 <- as.matrix((TreatC.w4 - TreatN.w4) - (TreatC.w8 - TreatN.w8)) # w4 (N vs. C)  vs. w8 (N vs. C)
C.w6vs8 <- as.matrix((TreatC.w6 - TreatN.w6) - (TreatC.w8 - TreatN.w8)) # w6 (N vs. C)  vs. w8 (N vs. C)

C.earlyvslate <- as.matrix(((TreatC.w2 - TreatN.w2)+(TreatC.w4 - TreatN.w4)) - 
                             ((TreatC.w6 - TreatN.w6) +(TreatC.w8 - TreatN.w8))) # w2+w4(N vs. C)  vs. w6+w8 (N vs. C), not used

W.w2vs4 <- as.matrix((TreatW.w2 - TreatN.w2) - (TreatW.w4 - TreatN.w4)) # w2 (N vs. W)  vs. w4 (N vs. W)
W.w2vs6 <- as.matrix((TreatW.w2 - TreatN.w2) - (TreatW.w6 - TreatN.w6)) # w2 (N vs. W)  vs. w6 (N vs. W)
W.w2vs8 <- as.matrix((TreatW.w2 - TreatN.w2) - (TreatW.w8 - TreatN.w8)) # w2 (N vs. W)  vs. w8 (N vs. W)
W.w4vs6 <- as.matrix((TreatW.w4 - TreatN.w4) - (TreatW.w6 - TreatN.w6)) # w4 (N vs. W)  vs. w6 (N vs. W)
W.w4vs8 <- as.matrix((TreatW.w4 - TreatN.w4) - (TreatW.w8 - TreatN.w8)) # w4 (N vs. W)  vs. w8 (N vs. W)
W.w6vs8 <- as.matrix((TreatW.w6 - TreatN.w6) - (TreatW.w8 - TreatN.w8)) # w6 (N vs. W)  vs. w8 (N vs. W)

W.earlyvslate <- as.matrix(((TreatW.w2 - TreatN.w2)+(TreatW.w4 - TreatN.w4)) - 
                             ((TreatW.w6 - TreatN.w6) +(TreatW.w8 - TreatN.w8))) # w2+w4(N vs. W)  vs. w6+w8 (N vs. W), not used


# Gather results ####
contr.matrix <- cbind(CvsN, WvsN, # Overall Treatment effects
                      w2_CvsN, w4_CvsN, w6_CvsN, w8_CvsN, # Cold within week
                      w2_WvsN, w4_WvsN, w6_WvsN, w8_WvsN, # Warm  within week
                      C.w2vs4, C.w2vs6, C.w2vs8, C.w4vs6, C.w4vs8, C.w6vs8, # Cold between weeks
                      W.w2vs4, W.w2vs6, W.w2vs8, W.w4vs6, W.w4vs8, W.w6vs8, # Warm between weeks
                      C.earlyvslate, W.earlyvslate # Between weeks: not sens vs. sens, not used
)

colnames(contr.matrix) <- c("CvsN", "WvsN",
                            "w2_CvsN", "w4_CvsN", "w6_CvsN", "w8_CvsN",
                            "w2_WvsN", "w4_WvsN", "w6_WvsN", "w8_WvsN", 
                            "C.w2vs4", "C.w2vs6", "C.w2vs8", "C.w4vs6", "C.w4vs8", "C.w6vs8", # Cold between weeks
                            "W.w2vs4", "W.w2vs6", "W.w2vs8", "W.w4vs6", "W.w4vs8", "W.w6vs8", # Warm between weeks
                            "C.earlyvslate", "W.earlyvslate" # Between weeks: not sens vs. sens, not used
)
contr.matrix


results.3h <- contrasts.fit(fit.3h, contrasts=contr.matrix)
results.3h <- eBayes(results.3h, robust=T)
plotSA(results.3h)

# Check p-value histograms ####
#-------------------------------------------------------
for(reslist in 1:length(colnames(results.3h$p.value))){
  plotname <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/pvalues/mod.3h/hist_", colnames(results.3h$p.value)[reslist], ".png", sep="")
  png(plotname)
  hist(results.3h$p.value[,reslist], col="dodgerblue", main=colnames(results.3h$p.value)[reslist]) # p-value histogram
  dev.off()
}
# how to check pattern: http://varianceexplained.org/statistics/interpreting-pvalue-histogram/ 
# many of histograms are either bimodal or hill shaped

# pvalue corrections for histograms that are off: https://www.huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html
# "Very often, if the assumed variance of the null distribution is too high, 
# we see hill–shaped p–value histogram. If the variance is too low, we get a U–shaped histogram, 
# with peaks at both ends."

# correct pvalues ####
#-------------------------------------------------------
results.3h$padj <- results.3h$p.value # create new list to store correct pvalues in

for(reslist in 1:length(colnames(results.3h$p.value))){

  # use t statistics returned by limma as input to fdrtool to re–estimate the p–values
  FDR.Res <- fdrtool::fdrtool(results.3h$t[,reslist], statistic= "normal", plot = T)
  
  # estimated null model variance
  print(FDR.Res$param[1, "sd"]) # theoretical one used =1
  
  # add the new BH–adjusted p–values add values to the results data frame
  results.3h$padj[,reslist]  <- p.adjust(FDR.Res$pval, method = "BH") #with Benjamini-Hochberg method
  
  # plot corrected p-values
  plotname <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/pvalues/mod.3h/hist_", colnames(results.3h$p.value)[reslist],"_corr", ".png", sep="")
  png(plotname)
  
  hist(FDR.Res$pval, col = "royalblue4", 
       main = paste(colnames(results.3h$p.value)[reslist], "correct null model", sep=", "), xlab = "CORRECTED p-values")
  
  dev.off()
}


# for 24h ####
#-------------------------------------------------------
mod24h <- filter(sampleinfo, Timepoint!="3")

mod_mat <- fit.24h$design
mod_mat

# Model parameters to make contrasts
TreatN <- colMeans(mod_mat[mod24h$Treatment2 == "N",])
TreatC <- colMeans(mod_mat[mod24h$Treatment2 == "C",])
TreatW <- colMeans(mod_mat[mod24h$Treatment2 == "W",])

TreatN.w2 <- colMeans(mod_mat[mod24h$Treatment2 == "N" & mod24h$Treat_week=="w2", ])
TreatC.w2 <- colMeans(mod_mat[mod24h$Treatment2 == "C" & mod24h$Treat_week=="w2", ])
TreatW.w2 <- colMeans(mod_mat[mod24h$Treatment2 == "W" & mod24h$Treat_week=="w2", ])

TreatN.w4 <- colMeans(mod_mat[mod24h$Treatment2 == "N" & mod24h$Treat_week=="w4", ])
TreatC.w4 <- colMeans(mod_mat[mod24h$Treatment2 == "C" & mod24h$Treat_week=="w4", ])
TreatW.w4 <- colMeans(mod_mat[mod24h$Treatment2 == "W" & mod24h$Treat_week=="w4", ])

TreatN.w6 <- colMeans(mod_mat[mod24h$Treatment2 == "N" & mod24h$Treat_week=="w6", ])
TreatC.w6 <- colMeans(mod_mat[mod24h$Treatment2 == "C" & mod24h$Treat_week=="w6", ])
TreatW.w6 <- colMeans(mod_mat[mod24h$Treatment2 == "W" & mod24h$Treat_week=="w6", ])

TreatN.w8 <- colMeans(mod_mat[mod24h$Treatment2 == "N" & mod24h$Treat_week=="w8", ])
TreatC.w8 <- colMeans(mod_mat[mod24h$Treatment2 == "C" & mod24h$Treat_week=="w8", ])
TreatW.w8 <- colMeans(mod_mat[mod24h$Treatment2 == "W" & mod24h$Treat_week=="w8", ])


# Overal Treatment effect ####
CvsN <- as.matrix(TreatC - TreatN)
WvsN <- as.matrix(TreatW - TreatN)


# Within week ####
# N vs. C
w2_CvsN <- as.matrix(TreatC.w2 - TreatN.w2)
w4_CvsN <- as.matrix(TreatC.w4 - TreatN.w4)
w6_CvsN <- as.matrix(TreatC.w6 - TreatN.w6)
w8_CvsN <- as.matrix(TreatC.w8 - TreatN.w8)

# N vs. W
w2_WvsN <- as.matrix(TreatW.w2 - TreatN.w2)
w4_WvsN <- as.matrix(TreatW.w4 - TreatN.w4)
w6_WvsN <- as.matrix(TreatW.w6 - TreatN.w6)
w8_WvsN <- as.matrix(TreatW.w8 - TreatN.w8)


# Between weeks ####
# individual week comparisons
C.w2vs4 <- as.matrix((TreatC.w2 - TreatN.w2) - (TreatC.w4 - TreatN.w4)) # w2 (N vs. C)  vs. w4 (N vs. C)
C.w2vs6 <- as.matrix((TreatC.w2 - TreatN.w2) - (TreatC.w6 - TreatN.w6)) # w2 (N vs. C)  vs. w6 (N vs. C)
C.w2vs8 <- as.matrix((TreatC.w2 - TreatN.w2) - (TreatC.w8 - TreatN.w8)) # w2 (N vs. C)  vs. w8 (N vs. C)
C.w4vs6 <- as.matrix((TreatC.w4 - TreatN.w4) - (TreatC.w6 - TreatN.w6)) # w4 (N vs. C)  vs. w6 (N vs. C)
C.w4vs8 <- as.matrix((TreatC.w4 - TreatN.w4) - (TreatC.w8 - TreatN.w8)) # w4 (N vs. C)  vs. w8 (N vs. C)
C.w6vs8 <- as.matrix((TreatC.w6 - TreatN.w6) - (TreatC.w8 - TreatN.w8)) # w6 (N vs. C)  vs. w8 (N vs. C)

C.earlyvslate <- as.matrix(((TreatC.w2 - TreatN.w2)+(TreatC.w4 - TreatN.w4)) - 
                             ((TreatC.w6 - TreatN.w6) +(TreatC.w8 - TreatN.w8))) # w2+w4(N vs. C)  vs. w6+w8 (N vs. C), not used

W.w2vs4 <- as.matrix((TreatW.w2 - TreatN.w2) - (TreatW.w4 - TreatN.w4)) # w2 (N vs. W)  vs. w4 (N vs. W)
W.w2vs6 <- as.matrix((TreatW.w2 - TreatN.w2) - (TreatW.w6 - TreatN.w6)) # w2 (N vs. W)  vs. w6 (N vs. W)
W.w2vs8 <- as.matrix((TreatW.w2 - TreatN.w2) - (TreatW.w8 - TreatN.w8)) # w2 (N vs. W)  vs. w8 (N vs. W)
W.w4vs6 <- as.matrix((TreatW.w4 - TreatN.w4) - (TreatW.w6 - TreatN.w6)) # w4 (N vs. W)  vs. w6 (N vs. W)
W.w4vs8 <- as.matrix((TreatW.w4 - TreatN.w4) - (TreatW.w8 - TreatN.w8)) # w4 (N vs. W)  vs. w8 (N vs. W)
W.w6vs8 <- as.matrix((TreatW.w6 - TreatN.w6) - (TreatW.w8 - TreatN.w8)) # w6 (N vs. W)  vs. w8 (N vs. W)

W.earlyvslate <- as.matrix(((TreatW.w2 - TreatN.w2)+(TreatW.w4 - TreatN.w4)) - 
                             ((TreatW.w6 - TreatN.w6) +(TreatW.w8 - TreatN.w8))) # w2+w4(N vs. W)  vs. w6+w8 (N vs. W), not used

# Gather results ####
contr.matrix <- cbind(CvsN, WvsN, # Overall Treatment effects
                      w2_CvsN, w4_CvsN, w6_CvsN, w8_CvsN, # Cold within week
                      w2_WvsN, w4_WvsN, w6_WvsN, w8_WvsN, # Warm  within week
                      C.w2vs4, C.w2vs6, C.w2vs8, C.w4vs6, C.w4vs8, C.w6vs8, # Cold between weeks
                      W.w2vs4, W.w2vs6, W.w2vs8, W.w4vs6, W.w4vs8, W.w6vs8, # Warm between weeks
                      C.earlyvslate, W.earlyvslate # Between weeks: not sens vs. sens, not used
)

colnames(contr.matrix) <- c("CvsN", "WvsN",
                            "w2_CvsN", "w4_CvsN", "w6_CvsN", "w8_CvsN",
                            "w2_WvsN", "w4_WvsN", "w6_WvsN", "w8_WvsN", 
                            "C.w2vs4", "C.w2vs6", "C.w2vs8", "C.w4vs6", "C.w4vs8", "C.w6vs8", # Cold between weeks
                            "W.w2vs4", "W.w2vs6", "W.w2vs8", "W.w4vs6", "W.w4vs8", "W.w6vs8", # Warm between weeks
                            "C.earlyvslate", "W.earlyvslate" # Between weeks: not sens vs. sens, not used
)
contr.matrix


results.24h <- contrasts.fit(fit.24h, contrasts=contr.matrix)
results.24h <- eBayes(results.24h, robust=TRUE)
plotSA(results.24h)


# Check p-value histograms ####
#-------------------------------------------------------
for(reslist in 1:length(colnames(results.24h$p.value))){
  plotname <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/pvalues/mod.24h/hist_", colnames(results.24h$p.value)[reslist], ".png", sep="")
  png(plotname)
  hist(results.24h$p.value[,reslist], col="dodgerblue", main=colnames(results.24h$p.value)[reslist]) # p-value histogram
  dev.off()
}
# how to check pattern: http://varianceexplained.org/statistics/interpreting-pvalue-histogram/ 
# many of histograms are either bimodal or hill shaped

# pvalue corrections for histograms that are off: https://www.huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html
# "Very often, if the assumed variance of the null distribution is too high, 
# we see hill–shaped p–value histogram. If the variance is too low, we get a U–shaped histogram, 
# with peaks at both ends."


# Correct pvalues ####
#-------------------------------------------------------
results.24h$padj <- results.24h$p.value

for(reslist in 1:length(colnames(results.24h$p.value))){

  # use t–statistics returned by limma as input to fdrtool to re–estimate the p–values
  FDR.Res <- fdrtool::fdrtool(results.24h$t[,reslist], statistic= "normal", plot = T)
  
  # estimated null model variance
  print(FDR.Res$param[1, "sd"]) # theoretical one used by DESeq2=1
  
  # add the new BH–adjusted p–values add values to the results data frame
  results.24h$padj[,reslist]  <- p.adjust(FDR.Res$pval, method = "BH") #Benjamini-Hochberg
  
  # plot corrected p-values
  plotname <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/pvalues/mod.24h/hist_", colnames(results.24h$p.value)[reslist],"_corr", ".png", sep="")
  png(plotname)
  
  hist(FDR.Res$pval, col = "royalblue4", 
       main = paste(colnames(results.24h$p.value)[reslist], "correct null model", sep=", "), xlab = "CORRECTED p-values")
  
  dev.off()
}


#-------------------------------------------------------
# Results
#-------------------------------------------------------
threshold <- 0.01 #FDR threshold I want to use ####


# for 3h model
#-------------------------------------------------------
topTable(results.3h, n=Inf, sort.by="none")[1:22] %>% head 
# gives log2FoldChanges for each contrast, same as $coefficients in results.3h
names(results.3h)
results.3h %>% head #padj = corrected FDR, coefficients=log2FoldChange

# Quick overview numbers of DEGs per contrast
sign.results.3h <- results.3h$padj<threshold
dim(sign.results.3h) # 24 contrasts

nr.sign.3h <- colCounts(sign.results.3h)
names(nr.sign.3h) <- colnames(sign.results.3h)
nr.sign.3h
sum(nr.sign.3h)

# Make list of 24 contrasts, with each list item containing logFC and FDR for each contrast
res_list.3h <- list() 

for(contrast in 1:length(colnames(results.3h$coefficients))){
  res <- as.data.frame(results.3h$coefficients[,contrast])
  colnames(res) <- c("log2FoldChange")
  res$padj <- results.3h$padj[,contrast]
  head(res)
  
  res_list.3h[[contrast]] <- res
  names(res_list.3h)[contrast] <- colnames(results.3h$coefficients)[contrast]
  
}
lapply(res_list.3h, head)
names(res_list.3h)

try <- res_list.3h[["w4_CvsN"]]
filter(try, padj<threshold) %>% row.names


# for 24h model
#-------------------------------------------------------

# Quick overview numbers of DEGs per contrast
sign.results.24h <- results.24h$padj<threshold
dim(sign.results.24h) # 24 contrasts

nr.sign.24h <- colCounts(sign.results.24h)
names(nr.sign.24h) <- colnames(sign.results.24h)
nr.sign.24h
sum(nr.sign.24h)

# Make list of 22 contrasts, with each list item containing logFC and FDR for each contrast
res_list.24h <- list() 

for(contrast in 1:length(colnames(results.24h$coefficients))){
  res <- as.data.frame(results.24h$coefficients[,contrast])
  colnames(res) <- c("log2FoldChange")
  res$padj <- results.24h$padj[,contrast]
  head(res)
  
  res_list.24h[[contrast]] <- res
  names(res_list.24h)[contrast] <- colnames(results.24h$coefficients)[contrast]
  
}
lapply(res_list.24h, head)
names(res_list.24h)


# Save models and results ####
# save(dfilt, v.3h, v.24h, fit.3h, fit.24h, results.3h, results.24h,
#      cutoff, filtwhich.n2, res_list.3h, res_list.24h, file="analysis/_RNAseq/_results/limma_mod_allweeks_filtn2.RData")
# save a few data objects to use later so we don’t have to rerun everything

load("analysis/_RNAseq/_results/limma_mod_allweeks_filtn2.RData") # all weeks analyzed together
# load("analysis/_RNAseq/_results/limma_mod_allweeks_filtn2_newcontr.RData") # with extra 2 between week contrasts, not used

sessionInfo() %>% capture.output(file="analysis/_RNAseq/_src/env_limma.txt")
