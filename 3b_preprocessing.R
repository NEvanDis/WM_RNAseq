# Have a look at TransfExp2019 RNAseq data and see if patterns follow my expectations ####

# Open R project in top folder

# Set up environment ####
#-------------------------------------------------------
# Load packages
library(DESeq2)
library(tidyverse)
library(ggfortify) # for autoplot
library(readxl)
library(PCAtools)


# Load data ####
d_count <- read.csv("analysis/_RNAseq/_data/TransfExp2019_gene_count_matrix_v3cor.csv")
class(d_count)  
dim(d_count) # v1: 46373 gene_ids, v2: 29095 gene_ids (17278 gene_ids novel transcripts without homologs excluded by BLAST), 
# v3: 30201 gene_ids (filtered out optical duplicates)
# v3cor: 29113 gene_ids (made mistake with grepping novel transcripts, this is correct version)
head(d_count) #first column = gene_id, other 84 cols = each sample, ordered by sample_name


# Format sample meta data ####
info <- read.csv("./analysis/RNA_available-clutches.csv") #quantiles per available RNA sample
info <- filter(info, Year==2019) %>% mutate(Tube=as.character(Tube), dev_median=dst50_before.50.) # dev_stage at BEF time point
head(info)

# treatment, time point and female_ID are in sample_name
sampleinfo <- as.data.frame(names(d_count[,-1]))
names(sampleinfo) <- "Sample"
head(sampleinfo)

sampleinfo <- sampleinfo %>%
  mutate(Tube=gsub(".*_(\\d+)$", "\\1", Sample), Treatment=gsub("(.*)_\\d+$", "\\1", Sample))

sampleinfo <- sampleinfo %>%
  mutate(Timepoint=ifelse(Treatment=="BEF", Treatment, gsub("\\w_(\\d+)$", "\\1", Treatment)), 
         Treatment=ifelse(Treatment=="BEF", "N", gsub("(\\w)_\\d+$", "\\1", Treatment))) # for now BEF time point = 10C treatment only

# get treatment week, and female info (Area, Site, Tree, Nov_date) from freezer list
sampleinfo <- merge(sampleinfo, unique(info[,c("Tube", "Treat_week", "dev_median", "NovemberDate", "Area", "Site", "Tree")]), by="Tube")

sampleinfo <- sampleinfo %>%
  mutate(Treat_week=as.factor(paste("w", Treat_week, sep="")), Area=as.factor(Area), Site=as.factor(Site), Tree=as.factor(Tree), NovemberDate=as.factor(NovemberDate))

sampleinfo <- sampleinfo[,c("Sample", "Treat_week", "dev_median","Tube", "Treatment", "Timepoint", "Area", "Site", "Tree", "NovemberDate")] %>%
  arrange(Sample)

# incubator info
# only separation of treatment weeks, not treatments or time points

# 5C: Inc1 (w2, w4, w6, w8)
# 15C: Inc4 (w2, w4, w6, w8)
# 10C: Inc3 (w2), Inc9 (w6, w8), Inc11 (w4) -> only 10C samples from different treatment weeks were in different incubators

head(sampleinfo)
tail(sampleinfo)


# Format data ####
countdata <- d_count %>% 
  column_to_rownames('gene_id') %>% #give rows gene names
  select(sampleinfo$Sample) %>% # drop gene_id column
  as.matrix() # matrix=all cells have same variable class

class(countdata)  
dim(countdata)  
head(countdata)  # each row = one gene, each column = one sample


# Specify model and create DESeq2 object ####
levels(sampleinfo$Treat_week)
wm_fixed <- ~ 1 + Treat_week + NovemberDate + Site
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData=sampleinfo,
                              design = wm_fixed)


# Prefilter ####
#-------------------------------------------------------
# Filter out genes not expressed in >=2 replicates or more
filter <- counts(dds)>5 # note for each clutch for each timepoint/treatment combo which genes expressed with at least >5 (filter out lowly expressed genes)
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
rownames(filt) <- rownames(dds)
head(filt) # per timepoint/treatment combo, count number of times genes expressed >5 out of 3 replicates

table(rowCounts(filt>=2)>=1) # only keep genes expressed in >=2 replicates at least in one timepoint/treatment combi
table(rowCounts(filt>=3)>=1) # only keep genes expressed in all replicates at least in one timepoint/treatment combi

filtwhich.n2 <- rownames(filt[rowCounts(filt>=2)>=1,]) # get gene names to keep
filtwhich.n3 <- rownames(filt[rowCounts(filt>=3)>=1,]) # get gene names to keep
head(filtwhich.n2)

ddsfilt <- dds[rownames(dds) %in% filtwhich.n2,] # select these genes
dim(dds)
dim(ddsfilt) # 26453 genes expressed at least once in >=2 replicates at least in one timepoint/treatment combi, 2660 genes dropped

# store formatted raw data
countdata_unfilt <- countdata

rm(filt, filter, dds, info, countdata) # clean up


#-------------------------------------------------------
# QUALITY CHECKS ####
#-------------------------------------------------------
# count distributions and visualisation

# Check Raw reads ####
#-------------------------------------------------------
raw <- assay(ddsfilt)

hist(raw[,1]) # distribution gene 1 raw counts
summary(raw)# most samples similar pattern, very high max compared to median/mean/Q3

boxplot(raw, main="Raw counts", las=2) #only see outliers, las to make xlabels vertical

plot(rowMeans(raw), rowSds(raw),
     main="Raw counts: sd vs. mean", 
     xlim=c(0, 10000), # if leave out limits, see a few extreme outliers
     ylim=c(0, 5000))
abline(lm(rowSds(raw)~rowMeans(raw)), col="green") # most points around green line
# RNAseq data often the higher the mean, the higher the variance
# patterns will be driven by weight of genes with highest means
# *heteroscedasticity= dependency of variance on the mean


# Check transformed reads #### 
#-------------------------------------------------------
# make sure bias can be dealt with as expected for RNAseq data
vst <- assay(vst(ddsfilt, blind = FALSE))
summary(vst)
# Sequencing depth correction is done automatically for the vst and rlog via DESeq2

# blind = FALSE means that differences because of variables in the design
# will not contribute to the expected variance-mean trend of the experiment. 
# The experimental design is not used directly in the transformation, only in 
#  estimating the global amount of variability in the counts.

# Alternative: rlog= regularized log transformation (takes long time)
# rlog_counts <- rlog(countdata)
# summary(rlog_counts)

statusCol <- match(sampleinfo$Treatment, c("W", "N", "C")) + 1 #assign a number to each treatment
boxplot(vst, main="vst(counts)",
        xlab="", 
        ylab="vst counts",
        las=2,
        col=statusCol)
abline(h=median(vst), col="red") # all sample medians close to overall median
# No extreme outliers

plot(rowMeans(vst), rowSds(vst), 
     main='vst counts: sd vs mean')
abline(lm(rowSds(vst)~rowMeans(vst)), col="green")
# points close together, variance no longer as dependent on mean, but still a bit of a bump for low mean expression genes


#-------------------------------------------------------
# PCA tools ####
#-------------------------------------------------------
# user manual: https://www.rdocumentation.org/packages/PCAtools/versions/2.5.13

# check if samples cluster according to expectations
# e.g. replicate samples should cluster together, while greatest differences should be between different treatments
# Also check for outliers (and batch effects if any)

# Want:
# - SCREE plot
# - Table with all PCAs 
# - Test correlation between PCAs and variables of interest

# Run PCA
row.names(sampleinfo) <- sampleinfo$Sample # needed to run PCA
sampleinfo <- sampleinfo %>% mutate(Trw_num= as.numeric(gsub("w(\\d)","\\1",Treat_week)), 
                                    Tr_num=ifelse(Treatment=="N", 10, ifelse(Treatment=="W", 15, 5)))

p <- pca(vst, metadata = sampleinfo, removeVar = 0.1) ## -- removing the lower 10% of variables based on variance


# How many PCAs explain X % of the variance in expression counts
screeplot(p, axisLabSize = 18, titleLabSize = 22)+
  geom_hline(yintercept=80)

# Correlate PCAs to sample variables (only for continuous!)
eigencorplot(p, components= getComponents(p, 1:19), # 80% of variance
             metavars = c('dev_median', 'NovemberDate', 'Trw_num', 'Tr_num'))

# Plot PCAs in 2dimensional space
biplot(p, showLoadings = FALSE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5,
       colby = 'dev_median',
       shape='Treat_week',
       legendPosition = 'right') +
  scale_colour_gradient(low = 'blue', high = 'red2')+
  theme(axis.title = element_text(size=25), axis.text=element_text(size=20), 
        legend.title=element_text(size=20), legend.text=element_text(size=20))

biplot(p, showLoadings = FALSE,
       x = 'PC4', y = 'PC3',
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5,
       colby = 'NovemberDate',
       shape='Treat_week',
       legendPosition = 'right')

biplot(p, showLoadings = FALSE,
       x = 'PC15', y = 'PC14',
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5,
       colby = 'Treatment', colkey = c('C'='dodgerblue2', 'N'='black', 'W'='firebrick3'),
       shape='Treat_week',
       legendPosition = 'right')


# Combine 2dim plots into one plane
pairsplot(p,
          colby='Treatment', colkey = c('C'='dodgerblue2', 'N'='black', 'W'='firebrick3'),
          shape='Treat_week')

  
# PCA plot per Treat_week ####
#-------------------------------------------------------

# w2 ####
dds.w2 <- ddsfilt[, which(grepl("_136",colnames(ddsfilt)) | grepl("_326",colnames(ddsfilt)) | grepl("_473",colnames(ddsfilt)))]
vst.w2 <- assay(vst(dds.w2))
p.w2 <- pca(vst.w2, metadata = dds.w2@colData, removeVar = 0.1) ## -- removing the lower 10% of variables based on variance

biplot(p.w2, showLoadings = FALSE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5,
       colby = 'dev_median',
       shape='Treat_week',
       legendPosition = 'right') +
  scale_colour_gradient(low = 'blue', high = 'red2')


# w4 ####
dds.w4 <- ddsfilt[, which(grepl("_102",colnames(ddsfilt)) | grepl("_219",colnames(ddsfilt)) | grepl("_367",colnames(ddsfilt)))]
vst.w4 <- assay(vst(dds.w4))
p.w4 <- pca(vst.w4, metadata = dds.w4@colData, removeVar = 0.1) ## -- removing the lower 10% of variables based on variance

biplot(p.w4, showLoadings = FALSE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5,
       colby = 'dev_median',
       shape='Treat_week',
       legendPosition = 'right') +
  scale_colour_gradient(low = 'blue', high = 'red2')


# w6 ####
dds.w6 <- ddsfilt[, which(grepl("_128",colnames(ddsfilt)) | grepl("_390",colnames(ddsfilt)) | grepl("_471",colnames(ddsfilt)))]
vst.w6 <- assay(vst(dds.w6))
p.w6 <- pca(vst.w6, metadata = dds.w6@colData, removeVar = 0.1) ## -- removing the lower 10% of variables based on variance

biplot(p.w6, showLoadings = FALSE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5,
       colby = 'dev_median',
       shape='Treat_week',
       legendPosition = 'right') +
  scale_colour_gradient(low = 'blue', high = 'red2')


# w8 ####
dds.w8 <- ddsfilt[, which(grepl("_94",colnames(ddsfilt)) | grepl("_407",colnames(ddsfilt)) | grepl("_411",colnames(ddsfilt)))]
vst.w8 <- assay(vst(dds.w8))
p.w8 <- pca(vst.w8, metadata = dds.w8@colData, removeVar = 0.1) ## -- removing the lower 10% of variables based on variance

biplot(p.w8, showLoadings = FALSE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5,
       colby = 'dev_median',
       shape='Treat_week',
       legendPosition = 'right') +
  scale_colour_gradient(low = 'blue', high = 'red2')


#-------------------------------------------------------
# save formatted data ####
#-------------------------------------------------------
# save(countdata_unfilt, sampleinfo, # unfiltered data and sample meta info
#      filtwhich.n2, filtwhich.n3, # genes expressed in =>2 replicates or all 3 replicates at least at one timepoint/treatment combi
#      geneset1, geneset2, # genes that pop up in PCA analysis 1. treat_week related, 2. treatment related
#      file="analysis/_RNAseq/_data/preprocessed2_v3cor_filt5.RData")
# save a few data objects to use later so we don’t have to rerun everything


sessionInfo()