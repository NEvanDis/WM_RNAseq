# BiocManager::install("topGO")

# load required libraries ####
library(tidyverse)
library(DESeq2)
library(edgeR)
library(topGO) # loads three environments: GOBPTerm, GOMFTerm and GOCCTerm
library(Rgraphviz) # needed for plotting GO graphs
library(ViSEAGO) # cluster GO terms according to semantic similarity

library(cowplot)
theme_set(theme_cowplot()) #white background instead of grey -> don't load if want grey grid
library(htmlwidgets) # to save viewer objects
library(webshot) # to save viewer objects

#-------------------------------------------------------------------------------------------------
#### Set up environment
#-------------------------------------------------------------------------------------------------

threshold <- 0.01 #FDR threshold I want to use ####

# Define colMap, ref. https://github.com/Bioconductor-mirror/topGO/blob/master/vignettes/topGO.Rnw
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}


# Load results ####
#-------------------------------------------------------------------------------------------------
load("analysis/_RNAseq/_results/limma_mod_allweeks_filtn2.RData") # all weeks analyzed together with limma

head(dfilt$counts) # raw counts filtered by cutoff ~20/median(library size) in at least 2 replicates
head(dfilt$samples) # sample information incl. library sizes, normalisation factors

lapply(res_list.24h, head) # Log2FoldChange and padj for each gene in each contrast tested

rm(results.24h, results.3h, v.24h, v.3h) # remove data not needed here


# WGCNA analysis: description 10C samples
load("analysis/_RNAseq/_results/WGCNA_Results_CoExpr10C.Rdata")
head(modNames) # ANOVA results for each module
head(SiteInfo.sign) # genes assigned to significant modules with significant membership 
head(SiteInfo) # results for all genes


# Get annotation ####
#-------------------------------------------------------------------------------------------------
# Load custom annotation from file
# format of annotation file: "gene_ID<TAB>GO_ID1, GO_ID2, GO_ID3, ...." 1 gene per line
# also see: https://gist.github.com/slavailn/dcff753cc32fd9e053590308c831b057 
geneID2GO <- readMappings(file="analysis/_RNAseq/_annotation/TransfExp_WM_gene2GO_v2.txt")
str(head(geneID2GO))

# gene universe = all genes in dataset
geneUniverse <- row.names(dfilt$counts)
head(geneUniverse)

# Load DEGS with gene level annotation
DEGs_annot <- read_tsv("analysis/_RNAseq/_results/DEGs_annot_n2_p0.01.tsv", col_names=T) # table of all DEGs
length(unique(DEGs_annot$GeneID))
DEGs_IDs <- unique(DEGs_annot$GeneID) # store GeneIDs of DEGs


#-------------------------------------------------------------------------------------------------
#### ENRICHMENT: Overrepresentation analysis
#-------------------------------------------------------------------------------------------------

## ENRICHMENT PER CONTRAST SET ####
#-------------------------------------------------------------------------------------------------
table(DEGs_annot$Comparison)

# Explore topGO environments
BPterms <- ls(GOBPTerm) # list of GO terms from BP ontology
head(BPterms)

# Choose ontology
ont <- "MF" #ontologies Biological Processes (BP), Molecular Function (MF)

# create lists to store results in
GOdfs_list <- as.list(unique(DEGs_annot$Comparison))
names(GOdfs_list) <- unique(DEGs_annot$Comparison)

FisherRes_list <- as.list(unique(DEGs_annot$Comparison))
names(FisherRes_list) <- unique(DEGs_annot$Comparison)


# Loop through contrast sets 3x
for (set in 1:length(GOdfs_list)){
  comp <- names(GOdfs_list)[set]
  res <- filter(DEGs_annot, Comparison==comp)$GeneID %>% unique
  
  # Prefiltering ####
  # "The number of probes have a direct effect on the multiple testing adjustment of p-values. Too many probes
  # will result in too conservative adjusted p-values which can bias the result of tests like Fisher’s exact test."
  
  # Recommendation: filter out genes with low expression value and/or with very small variability across samples
  # = use filtered gene set as geneUniverse 23 961 genes
  
  # build topGOdata container object ####
  #-------------------------------------------------------------------------------------------------
  # contains: gene identifiers (and optionally their scores e.g. Pvalues); GO annotations; GO hierarchical structure; 
  # additional information needed for analysis. Additional info here =  list of differentially expressed genes
  
  # Note for each gene in the universe whether it's a DEG or not
  table(geneUniverse %in% res)
  geneList <- factor(as.integer(geneUniverse %in% res))
  names(geneList) <- c(geneUniverse)
  str(geneList) # named gene list DEG yes or no
  table(geneList)
  
  sampleGOdata <- new("topGOdata",
                      description = comp, # optional 
                      ontology = ont, # GO graph
                      allGenes = geneList, # gene universe as logical DEG yes or no
                      nodeSize = 5, #prune GO hierarchy from terms with less than X annotated genes (5-10 recommended)
                      annot = annFUN.gene2GO, gene2GO=geneID2GO) #function used to extract gene-to-GO mappings from affyLib object
  GOdfs_list[[set]] <- sampleGOdata
  GOdfs_list[[set]] # see summary of object
  
  
  # Overrepresentation analysis ####
  #-------------------------------------------------------------------------------------------------
  # testing the over-representation of GO terms within the group of differentially expressed genes
  
  resultFisher <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher") 
  # resultFisher.cl <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher") 
  # method=classic, each GO category is tested independently
  # not adjusted for multiple testing
  # method=elim, more conservative
  FisherRes_list[[set]] <- resultFisher # see summary of object
  FisherRes_list[[set]]

  
  # Down stream analysis tools ####
  #-------------------------------------------------------------------------------------------------
  signNodes <- table(resultFisher@score<0.01)["TRUE"] # number of significant GO terms P<0.01
  if(is.na(signNodes)==TRUE) {signNod <- 1
  } else signNod <- signNodes
  
  # gather results in a table ####
  allRes <- GenTable(sampleGOdata, elimFisher = resultFisher,
                     #classicFisher = resultFisher.cl,
                     orderBy = "elimFisher", ranksOf = "elimFisher", topNodes = signNod) # only keep top 10 results of elimKS results
  allRes$set <- comp
  colnames(allRes)[6] <- "elimFisher"
  allRes
  
  
  # Expand table with gene level info ####
  # pull GO2Gene annotation for significant GOterms
  
  if(length(subset(allRes, elimFisher<=0.01)$GO.ID)>0){
    # retrieve genes2GO list from the "expanded" annotation in GOdata
    allGO <- genesInTerm(sampleGOdata) # pull the GO2Gene annotation
    head(allGO)
    
    # pull GO2Gene annotation for DEGs of set
    DEGs_GOannot <- lapply(allGO,function(x) x[x %in% names(geneList[geneList==1])])
    DEGs_GOannot <- DEGs_GOannot[lapply(DEGs_GOannot,length)>0] # remove GO terms without hits
    head(DEGs_GOannot)
    
    DEGs_GOannot[[allRes$GO.ID[1]]]
    DEGs_signGO <- DEGs_GOannot[subset(allRes, elimFisher<=0.01)$GO.ID]
    DEGs_signGO <- DEGs_signGO[lapply(DEGs_signGO,length)>0] # remove GO terms without hits
    head(DEGs_signGO)
    
    gen <- lapply(DEGs_signGO, as.data.frame) # format as data frame with one row per gene
    gen <- mapply(cbind, gen, "GO.ID"=names(gen), SIMPLIFY=F)
    gen <- lapply(names(gen), function(x) setNames(gen[[x]], c("GeneID", "GO.ID")) )
    gen <- data.table::rbindlist(gen)
    head(gen)
    
    head(allRes)
    allRes_ext <- merge(allRes, gen, by="GO.ID") # complete table with one row per gene
    allRes_ext <- merge(allRes_ext, unique(DEGs_annot[,c("GeneID", "GeneName", "Prot_descr")]), by="GeneID") # add annotation to the genes
    head(allRes_ext)
    
  } else{
    allRes_ext <- character(0)
  }
  
  # Visualize distribution significant GO terms over GO graph ####
  plotname <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/GOanalysis/", comp, "_GOgraph_", ont, ".pdf", sep="")
  pdf(plotname, width = 16, height = 9, family="sans", colormodel="srgb")
  showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = signNod, useInfo = 'all') # based on 5 most sign GO terms of results
  dev.off()
  
  # Store results
  head(allRes)
  head(allRes_ext)
  
  allRes <- subset(allRes, elimFisher<=0.01) %>% arrange(elimFisher)
  tablename <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/GOanalysis/", comp, "_GOtable_", ont, ".tsv", sep="")
  write.table(allRes, file=tablename, row.names=FALSE, sep="\t")
  
  if(is.null(allRes_ext$elimFisher)==F){
    allRes_ext <- allRes_ext %>% dplyr::select(set, GO.ID, Term, Annotated, Significant, Expected, elimFisher, GeneID, GeneName, Prot_descr) %>%
      arrange(elimFisher, GO.ID, GeneID, Prot_descr)
    tablename <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/GOanalysis/", comp, "_GOtable_", ont, "_ext", ".tsv", sep="")
    write.table(allRes_ext, file=tablename, row.names=FALSE, sep="\t")
    }
  
  # end loop
}

head(GOdfs_list)
names(GOdfs_list)

# save data ####
save(GOdfs_list, FisherRes_list, file=paste("analysis/_RNAseq/_results/limma_GOanalysis_sets_", ont, ".RData", sep=""))
# save a few data objects to use later so we don’t have to rerun everything
#-------------------------------------------------------------------------------------------------


# ENRICHMENT OVERLAP Y/N WITH WGCNA ANALYSIS ####
#-------------------------------------------------------------------------------------------------
table(unique(DEGs_annot$GeneID) %in% SiteInfo.sign$gene) 
overlap <- list()
overlap[[1]] <- unique(filter(DEGs_annot, GeneID %in% SiteInfo.sign$gene)$GeneID) # overlapping
overlap[[2]] <- unique(filter(DEGs_annot, !GeneID %in% SiteInfo.sign$gene)$GeneID) # not overlapping
names(overlap) <- c("overlapping", "not_overlapping")

# Choose ontology
ont <- "MF" #ontologies Biological Processes (BP), Molecular Function (MF)

# create lists to store results in
GOdfs_list <- overlap
FisherRes_list <- overlap


# Overrepresentation for 2 gene sets
for (set in 1:length(GOdfs_list)){
  comp <- names(GOdfs_list)[set]
  res <- overlap[[set]]
  
  # Prefiltering ####
  # "The number of probes have a direct effect on the multiple testing adjustment of p-values. Too many probes
  # will result in too conservative adjusted p-values which can bias the result of tests like Fisher’s exact test."
  
  # Recommendation: filter out genes with low expression value and/or with very small variability across samples
  # = use filtered gene set as geneUniverse 23 961 genes
  
  # build topGOdata container object ####
  #-------------------------------------------------------------------------------------------------
  # contains: gene identifiers (and optionally their scores e.g. Pvalues); GO annotations; GO hierarchical structure; 
  # additional information needed for analysis. Additional info here =  list of differentially expressed genes
  
  # Note for each gene in the universe whether it's a DEG or not
  table(geneUniverse %in% res)
  geneList <- factor(as.integer(geneUniverse %in% res))
  names(geneList) <- c(geneUniverse)
  str(geneList) # named gene list DEG yes or no
  table(geneList)
  
  sampleGOdata <- new("topGOdata",
                      description = comp, # optional 
                      ontology = ont, # GO graph
                      allGenes = geneList, # gene universe as logical DEG yes or no
                      nodeSize = 5, #prune GO hierarchy from terms with less than X annotated genes (5-10 recommended)
                      annot = annFUN.gene2GO, gene2GO=geneID2GO) #function used to extract gene-to-GO mappings from affyLib object
  GOdfs_list[[set]] <- sampleGOdata
  GOdfs_list[[set]] # see summary of object
  
  
  # Overrepresentation analysis ####
  #-------------------------------------------------------------------------------------------------
  # testing the over-representation of GO terms within the group of differentially expressed genes
  
  resultFisher <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher") 
  # resultFisher.cl <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher") 
  # method=classic, each GO category is tested independently
  # not adjusted for multiple testing
  # method=elim, more conservative
  FisherRes_list[[set]] <- resultFisher # see summary of object
  FisherRes_list[[set]]
  
  
  # Down stream analysis tools ####
  #-------------------------------------------------------------------------------------------------
  signNodes <- table(resultFisher@score<0.01)["TRUE"] # number of significant GO terms P<0.01
  if(is.na(signNodes)==TRUE) {signNod <- 1
  } else signNod <- signNodes
  
  # gather results in a table ####
  allRes <- GenTable(sampleGOdata, elimFisher = resultFisher,
                     #classicFisher = resultFisher.cl,
                     orderBy = "elimFisher", ranksOf = "elimFisher", topNodes = signNod) # only keep top 10 results of elimKS results
  allRes$set <- comp
  colnames(allRes)[6] <- "elimFisher"
  allRes
  
  
  # Expand table with gene level info ####
  # pull GO2Gene annotation for significant GOterms
  
  if(length(subset(allRes, elimFisher<=0.01)$GO.ID)>0){
    # retrieve genes2GO list from the "expanded" annotation in GOdata
    allGO <- genesInTerm(sampleGOdata) # pull the GO2Gene annotation
    head(allGO)
    
    # pull GO2Gene annotation for DEGs of set
    DEGs_GOannot <- lapply(allGO,function(x) x[x %in% names(geneList[geneList==1])])
    DEGs_GOannot <- DEGs_GOannot[lapply(DEGs_GOannot,length)>0] # remove GO terms without hits
    head(DEGs_GOannot)
    
    DEGs_GOannot[[allRes$GO.ID[1]]]
    DEGs_signGO <- DEGs_GOannot[subset(allRes, elimFisher<=0.01)$GO.ID]
    DEGs_signGO <- DEGs_signGO[lapply(DEGs_signGO,length)>0] # remove GO terms without hits
    head(DEGs_signGO)
    
    gen <- lapply(DEGs_signGO, as.data.frame) # format as data frame with one row per gene
    gen <- mapply(cbind, gen, "GO.ID"=names(gen), SIMPLIFY=F)
    gen <- lapply(names(gen), function(x) setNames(gen[[x]], c("GeneID", "GO.ID")) )
    gen <- data.table::rbindlist(gen)
    head(gen)
    
    head(allRes)
    allRes_ext <- merge(allRes, gen, by="GO.ID") # complete table with one row per gene
    allRes_ext <- merge(allRes_ext, unique(DEGs_annot[,c("GeneID", "GeneName", "Prot_descr")]), by="GeneID") # add annotation to the genes
    head(allRes_ext)
    
  } else{
    allRes_ext <- character(0)
  }
  
  # Visualize distribution significant GO terms over GO graph ####
  plotname <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/GOanalysis/", comp, "_GOgraph_", ont, ".pdf", sep="")
  pdf(plotname, width = 16, height = 9, family="sans", colormodel="srgb")
  showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = signNod, useInfo = 'all') # based on 5 most sign GO terms of results
  dev.off()
  
  # Store results
  head(allRes)
  head(allRes_ext)
  
  allRes <- subset(allRes, elimFisher<=0.01) %>% arrange(elimFisher)
  tablename <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/GOanalysis/", comp, "_GOtable_", ont, ".tsv", sep="")
  write.table(allRes, file=tablename, row.names=FALSE, sep="\t")
  
  if(is.null(allRes_ext$elimFisher)==F){
    allRes_ext <- allRes_ext %>% dplyr::select(set, GO.ID, Term, Annotated, Significant, Expected, elimFisher, GeneID, GeneName, Prot_descr) %>%
      arrange(elimFisher, GO.ID, GeneID, Prot_descr)
    tablename <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/GOanalysis/", comp, "_GOtable_", ont, "_ext", ".tsv", sep="")
    write.table(allRes_ext, file=tablename, row.names=FALSE, sep="\t")
  }
  # end loop
}

names(GOdfs_list)

# save data ####
save(GOdfs_list, FisherRes_list, file=paste("analysis/_RNAseq/_results/limma_GOanalysis_overlap_", ont, ".RData", sep=""))
# save a few data objects to use later so we don’t have to rerun everything
#-------------------------------------------------------------------------------------------------


# ENRICHMENT WEEK SPECIFIC GENE SETS ####
#-------------------------------------------------------------------------------------------------
withinweek <- list(unique(filter(DEGs_annot, origin=="w2")$GeneID), unique(filter(DEGs_annot, origin=="w4")$GeneID), 
                      unique(filter(DEGs_annot, origin=="w6")$GeneID), unique(filter(DEGs_annot, origin=="w8")$GeneID))
names(withinweek) <- c("Week2", "Week4", "Week6", "Week8")
lapply(withinweek, length)

# Choose ontology
ont <- "MF" #ontologies Biological Processes (BP), Molecular Function (MF)

# create lists to store results in
GOdfs_list <- withinweek
FisherRes_list <- withinweek


# Overrepresentation for 4 week specific gene sets from Within-week contrasts
for (set in 1:length(GOdfs_list)){
  comp <- names(GOdfs_list)[set]
  res <- withinweek[[set]]
  
  # Prefiltering ####
  # "The number of probes have a direct effect on the multiple testing adjustment of p-values. Too many probes
  # will result in too conservative adjusted p-values which can bias the result of tests like Fisher’s exact test."
  
  # Recommendation: filter out genes with low expression value and/or with very small variability across samples
  # = use filtered gene set as geneUniverse 23 961 genes
  
  # build topGOdata container object ####
  #-------------------------------------------------------------------------------------------------
  # contains: gene identifiers (and optionally their scores e.g. Pvalues); GO annotations; GO hierarchical structure; 
  # additional information needed for analysis. Additional info here =  list of differentially expressed genes
  
  # Note for each gene in the universe whether it's a DEG or not
  table(geneUniverse %in% res)
  geneList <- factor(as.integer(geneUniverse %in% res))
  names(geneList) <- c(geneUniverse)
  str(geneList) # named gene list DEG yes or no
  table(geneList)
  
  sampleGOdata <- new("topGOdata",
                      description = comp, # optional 
                      ontology = ont, # GO graph
                      allGenes = geneList, # gene universe as logical DEG yes or no
                      nodeSize = 5, #prune GO hierarchy from terms with less than X annotated genes (5-10 recommended)
                      annot = annFUN.gene2GO, gene2GO=geneID2GO) #function used to extract gene-to-GO mappings from affyLib object
  GOdfs_list[[set]] <- sampleGOdata
  GOdfs_list[[set]] # see summary of object
  
  
  # Overrepresentation analysis ####
  #-------------------------------------------------------------------------------------------------
  # testing the over-representation of GO terms within the group of differentially expressed genes
  
  resultFisher <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher") 
  # resultFisher.cl <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher") 
  # method=classic, each GO category is tested independently
  # not adjusted for multiple testing
  # method=elim, more conservative
  FisherRes_list[[set]] <- resultFisher # see summary of object
  FisherRes_list[[set]]
  
  
  # Down stream analysis tools ####
  #-------------------------------------------------------------------------------------------------
  signNodes <- table(resultFisher@score<0.01)["TRUE"] # number of significant GO terms P<0.01
  if(is.na(signNodes)==TRUE) {signNod <- 1
  } else signNod <- signNodes
  
  # gather results in a table ####
  allRes <- GenTable(sampleGOdata, elimFisher = resultFisher,
                     #classicFisher = resultFisher.cl,
                     orderBy = "elimFisher", ranksOf = "elimFisher", topNodes = signNod) # only keep top 10 results of elimKS results
  allRes$set <- comp
  colnames(allRes)[6] <- "elimFisher"
  allRes
  
  
  # Expand table with gene level info ####
  # pull GO2Gene annotation for significant GOterms
  
  if(length(subset(allRes, elimFisher<=0.01)$GO.ID)>0){
    # retrieve genes2GO list from the "expanded" annotation in GOdata
    allGO <- genesInTerm(sampleGOdata) # pull the GO2Gene annotation
    head(allGO)
    
    # pull GO2Gene annotation for DEGs of set
    DEGs_GOannot <- lapply(allGO,function(x) x[x %in% names(geneList[geneList==1])])
    DEGs_GOannot <- DEGs_GOannot[lapply(DEGs_GOannot,length)>0] # remove GO terms without hits
    head(DEGs_GOannot)
    
    DEGs_GOannot[[allRes$GO.ID[1]]]
    DEGs_signGO <- DEGs_GOannot[subset(allRes, elimFisher<=0.01)$GO.ID]
    DEGs_signGO <- DEGs_signGO[lapply(DEGs_signGO,length)>0] # remove GO terms without hits
    head(DEGs_signGO)
    
    gen <- lapply(DEGs_signGO, as.data.frame) # format as data frame with one row per gene
    gen <- mapply(cbind, gen, "GO.ID"=names(gen), SIMPLIFY=F)
    gen <- lapply(names(gen), function(x) setNames(gen[[x]], c("GeneID", "GO.ID")) )
    gen <- data.table::rbindlist(gen)
    head(gen)
    
    head(allRes)
    allRes_ext <- merge(allRes, gen, by="GO.ID") # complete table with one row per gene
    allRes_ext <- merge(allRes_ext, unique(DEGs_annot[,c("GeneID", "GeneName", "Prot_descr")]), by="GeneID") # add annotation to the genes
    head(allRes_ext)

  } else{
    allRes_ext <- character(0)
  }
  
  # Visualize distribution significant GO terms over GO graph ####
  plotname <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/GOanalysis/", comp, "_GOgraph_", ont, ".pdf", sep="")
  pdf(plotname, width = 16, height = 9, family="sans", colormodel="srgb")
  showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = signNod, useInfo = 'all') # based on 5 most sign GO terms of results
  dev.off()
  
  # Store results
  head(allRes)
  head(allRes_ext)
  
  allRes <- subset(allRes, elimFisher<=0.01) %>% arrange(elimFisher)
  tablename <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/GOanalysis/", comp, "_GOtable_", ont, ".tsv", sep="")
  write.table(allRes, file=tablename, row.names=FALSE, sep="\t")
  
  if(is.null(allRes_ext$elimFisher)==F){
    allRes_ext <- allRes_ext %>% dplyr::select(set, GO.ID, Term, Annotated, Significant, Expected, elimFisher, GeneID, GeneName, Prot_descr) %>%
      arrange(elimFisher, GO.ID, GeneID, Prot_descr)
    tablename <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/GOanalysis/", comp, "_GOtable_", ont, "_ext", ".tsv", sep="")
    write.table(allRes_ext, file=tablename, row.names=FALSE, sep="\t")
  }
  # end loop
}

names(GOdfs_list)

# save data ####
save(GOdfs_list, FisherRes_list, file=paste("analysis/_RNAseq/_results/limma_GOanalysis_weekspecific_", ont, ".RData", sep=""))
# save a few data objects to use later so we don’t have to rerun everything
#-------------------------------------------------------------------------------------------------


# ENRICHMENT BETWEEN-WEEK SPECIFIC GENE SETS ####
#-------------------------------------------------------------------------------------------------
betweenweek <- list(unique(filter(DEGs_annot, origin=="w4vs6")$GeneID), unique(filter(DEGs_annot, origin=="w6vs8")$GeneID))
names(betweenweek) <- c("w4vs6", "w6vs8")
lapply(betweenweek, length)

# Choose ontology
ont <- "MF" #ontologies Biological Processes (BP), Molecular Function (MF)

# create lists to store results in
GOdfs_list <- betweenweek
FisherRes_list <- betweenweek


# Overrepresentation for 2 BetweenWeek gene sets
for (set in 1:length(GOdfs_list)){
  comp <- names(GOdfs_list)[set]
  res <- betweenweek[[set]]
  
  # Prefiltering ####
  # "The number of probes have a direct effect on the multiple testing adjustment of p-values. Too many probes
  # will result in too conservative adjusted p-values which can bias the result of tests like Fisher’s exact test."
  
  # Recommendation: filter out genes with low expression value and/or with very small variability across samples
  # = use filtered gene set as geneUniverse 23 961 genes
  
  # build topGOdata container object ####
  #-------------------------------------------------------------------------------------------------
  # contains: gene identifiers (and optionally their scores e.g. Pvalues); GO annotations; GO hierarchical structure; 
  # additional information needed for analysis. Additional info here =  list of differentially expressed genes
  
  # Note for each gene in the universe whether it's a DEG or not
  table(geneUniverse %in% res)
  geneList <- factor(as.integer(geneUniverse %in% res))
  names(geneList) <- c(geneUniverse)
  str(geneList) # named gene list DEG yes or no
  table(geneList)
  
  sampleGOdata <- new("topGOdata",
                      description = comp, # optional 
                      ontology = ont, # GO graph
                      allGenes = geneList, # gene universe as logical DEG yes or no
                      nodeSize = 5, #prune GO hierarchy from terms with less than X annotated genes (5-10 recommended)
                      annot = annFUN.gene2GO, gene2GO=geneID2GO) #function used to extract gene-to-GO mappings from affyLib object
  GOdfs_list[[set]] <- sampleGOdata
  GOdfs_list[[set]] # see summary of object
  
  
  # Overrepresentation analysis ####
  #-------------------------------------------------------------------------------------------------
  # testing the over-representation of GO terms within the group of differentially expressed genes
  
  resultFisher <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher") 
  # resultFisher.cl <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher") 
  # method=classic, each GO category is tested independently
  # not adjusted for multiple testing
  # method=elim, more conservative
  FisherRes_list[[set]] <- resultFisher # see summary of object
  FisherRes_list[[set]]
  
  
  # Down stream analysis tools ####
  #-------------------------------------------------------------------------------------------------
  signNodes <- table(resultFisher@score<0.01)["TRUE"] # number of significant GO terms P<0.01
  if(is.na(signNodes)==TRUE) {signNod <- 1
  } else signNod <- signNodes
  
  # gather results in a table ####
  allRes <- GenTable(sampleGOdata, elimFisher = resultFisher,
                     #classicFisher = resultFisher.cl,
                     orderBy = "elimFisher", ranksOf = "elimFisher", topNodes = signNod) # only keep top 10 results of elimKS results
  allRes$set <- comp
  colnames(allRes)[6] <- "elimFisher"
  allRes
  
  
  # Expand table with gene level info ####
  # pull GO2Gene annotation for significant GOterms
  
  if(length(subset(allRes, elimFisher<=0.01)$GO.ID)>0){
    # retrieve genes2GO list from the "expanded" annotation in GOdata
    allGO <- genesInTerm(sampleGOdata) # pull the GO2Gene annotation
    head(allGO)
    
    # pull GO2Gene annotation for DEGs of set
    DEGs_GOannot <- lapply(allGO,function(x) x[x %in% names(geneList[geneList==1])])
    DEGs_GOannot <- DEGs_GOannot[lapply(DEGs_GOannot,length)>0] # remove GO terms without hits
    head(DEGs_GOannot)
    
    DEGs_GOannot[[allRes$GO.ID[1]]]
    DEGs_signGO <- DEGs_GOannot[subset(allRes, elimFisher<=0.01)$GO.ID]
    DEGs_signGO <- DEGs_signGO[lapply(DEGs_signGO,length)>0] # remove GO terms without hits
    head(DEGs_signGO)
    
    gen <- lapply(DEGs_signGO, as.data.frame) # format as data frame with one row per gene
    gen <- mapply(cbind, gen, "GO.ID"=names(gen), SIMPLIFY=F)
    gen <- lapply(names(gen), function(x) setNames(gen[[x]], c("GeneID", "GO.ID")) )
    gen <- data.table::rbindlist(gen)
    head(gen)
    
    head(allRes)
    allRes_ext <- merge(allRes, gen, by="GO.ID") # complete table with one row per gene
    allRes_ext <- merge(allRes_ext, unique(DEGs_annot[,c("GeneID", "GeneName", "Prot_descr")]), by="GeneID") # add annotation to the genes
    head(allRes_ext)
    
  } else{
    allRes_ext <- character(0)
  }
  
  # Visualize distribution significant GO terms over GO graph ####
  plotname <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/GOanalysis/", comp, "_GOgraph_", ont, ".pdf", sep="")
  pdf(plotname, width = 16, height = 9, family="sans", colormodel="srgb")
  showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = signNod, useInfo = 'all') # based on 5 most sign GO terms of results
  dev.off()
  
  # Store results
  head(allRes)
  head(allRes_ext)
  
  allRes <- subset(allRes, elimFisher<=0.01) %>% arrange(elimFisher)
  tablename <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/GOanalysis/", comp, "_GOtable_", ont, ".tsv", sep="")
  write.table(allRes, file=tablename, row.names=FALSE, sep="\t")
  
  if(is.null(allRes_ext$elimFisher)==F){
    allRes_ext <- allRes_ext %>% dplyr::select(set, GO.ID, Term, Annotated, Significant, Expected, elimFisher, GeneID, GeneName, Prot_descr) %>%
      arrange(elimFisher, GO.ID, GeneID, Prot_descr)
    tablename <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/GOanalysis/", comp, "_GOtable_", ont, "_ext", ".tsv", sep="")
    write.table(allRes_ext, file=tablename, row.names=FALSE, sep="\t")
  }
  # end loop
}

names(GOdfs_list)

# save data ####
save(GOdfs_list, FisherRes_list, file=paste("analysis/_RNAseq/_results/limma_GOanalysis_BetweenWeek_", ont, ".RData", sep=""))
# save a few data objects to use later so we don’t have to rerun everything
#-------------------------------------------------------------------------------------------------



# ENRICHMENT HM CLUSTERS ####
#-------------------------------------------------------------------------------------------------
clusters <- read.csv("analysis/_RNAseq/_results/DEGS_hm_clusters.csv")
clusters$subset <- paste(clusters$set, clusters$Comparison, sep="_")
head(clusters)
table(clusters$Cluster, clusters$subset)

# Choose ontology
ont <- "MF" #ontologies Biological Processes (BP), Molecular Function (MF)

# create lists to store results in
GOdfs_list <- list()
FisherRes_list <- list()


# Overrepresentation for identified hm clusters
for (set in 1:6){
  comp <- unique(clusters$subset)[set]
  fullres <- filter(clusters, subset==comp)
  
  GOdfs <- list()
  FishRes <- list()
  
  for (clust in 1:length(unique(fullres$Cluster))){
    res <- filter(fullres, Cluster==clust)$GeneID
    head(res)
    
    # Prefiltering ####
    # "The number of probes have a direct effect on the multiple testing adjustment of p-values. Too many probes
    # will result in too conservative adjusted p-values which can bias the result of tests like Fisher’s exact test."
    
    # Recommendation: filter out genes with low expression value and/or with very small variability across samples
    # = use filtered gene set as geneUniverse 23 961 genes
    
    # build topGOdata container object ####
    #-------------------------------------------------------------------------------------------------
    # contains: gene identifiers (and optionally their scores e.g. Pvalues); GO annotations; GO hierarchical structure; 
    # additional information needed for analysis. Additional info here =  list of differentially expressed genes
    
    # Note for each gene in the universe whether it's a DEG or not
    table(geneUniverse %in% res)
    geneList <- factor(as.integer(geneUniverse %in% res))
    names(geneList) <- c(geneUniverse)
    str(geneList) # named gene list DEG yes or no
    table(geneList)
    
    sampleGOdata <- new("topGOdata",
                        description = paste(comp, clust, sep=" "), # optional 
                        ontology = ont, # GO graph
                        allGenes = geneList, # gene universe as logical DEG yes or no
                        nodeSize = 5, #prune GO hierarchy from terms with less than X annotated genes (5-10 recommended)
                        annot = annFUN.gene2GO, gene2GO=geneID2GO) #function used to extract gene-to-GO mappings from affyLib object
    GOdfs[[clust]] <- sampleGOdata
    GOdfs[[clust]] # see summary of object
    
    
    # Overrepresentation analysis ####
    #-------------------------------------------------------------------------------------------------
    # testing the over-representation of GO terms within the group of differentially expressed genes
    
    resultFisher <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher") 
    # resultFisher.cl <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher") 
    # method=classic, each GO category is tested independently
    # not adjusted for multiple testing
    # method=elim, more conservative
    FishRes[[clust]] <- resultFisher # see summary of object
    FishRes[[clust]]
    
    
    # Down stream analysis tools ####
    #-------------------------------------------------------------------------------------------------
    signNodes <- table(resultFisher@score<0.01)["TRUE"] # number of significant GO terms P<0.01
    if(is.na(signNodes)==TRUE) {signNod <- 1
    } else signNod <- signNodes
    
    # gather results in a table ####
    allRes <- GenTable(sampleGOdata, elimFisher = resultFisher,
                       #classicFisher = resultFisher.cl,
                       orderBy = "elimFisher", ranksOf = "elimFisher", topNodes = signNod) # only keep top 10 results of elimKS results
    allRes$cluster <- paste(comp, clust, sep=" ")
    colnames(allRes)[6] <- "elimFisher"
    allRes
    
    
    # Expand table with gene level info ####
    # pull GO2Gene annotation for significant GOterms
    
    if(length(subset(allRes, elimFisher<=0.01)$GO.ID)>0){
      # retrieve genes2GO list from the "expanded" annotation in GOdata
      allGO <- genesInTerm(sampleGOdata) # pull the GO2Gene annotation
      head(allGO)
      
      # pull GO2Gene annotation for DEGs of set
      DEGs_GOannot <- lapply(allGO,function(x) x[x %in% names(geneList[geneList==1])])
      DEGs_GOannot <- DEGs_GOannot[lapply(DEGs_GOannot,length)>0] # remove GO terms without hits
      head(DEGs_GOannot)
      
      DEGs_GOannot[[allRes$GO.ID[1]]]
      DEGs_signGO <- DEGs_GOannot[subset(allRes, elimFisher<=0.01)$GO.ID]
      DEGs_signGO <- DEGs_signGO[lapply(DEGs_signGO,length)>0] # remove GO terms without hits
      head(DEGs_signGO)
      
      gen <- lapply(DEGs_signGO, as.data.frame) # format as data frame with one row per gene
      gen <- mapply(cbind, gen, "GO.ID"=names(gen), SIMPLIFY=F)
      gen <- lapply(names(gen), function(x) setNames(gen[[x]], c("GeneID", "GO.ID")) )
      gen <- data.table::rbindlist(gen)
      head(gen)
      
      head(allRes)
      allRes_ext <- merge(allRes, gen, by="GO.ID") # complete table with one row per gene
      allRes_ext <- merge(allRes_ext, unique(DEGs_annot[,c("GeneID", "GeneName", "Prot_descr")]), by="GeneID") # add annotation to the genes
      head(allRes_ext)
      
    } else{
      allRes_ext <- NULL
    }
    
    # Visualize distribution significant GO terms over GO graph ####
    plotname <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/GOanalysis/hm_clusters/hm_", comp, "_", clust, "_GOgraph_", ont, ".pdf", sep="")
    pdf(plotname, width = 16, height = 9, family="sans", colormodel="srgb")
    showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = signNod, useInfo = 'all') # based on 5 most sign GO terms of results
    dev.off()
    
    # Store results
    head(allRes)
    head(allRes_ext)
    
    allRes <- subset(allRes, elimFisher<=0.01) %>% arrange(elimFisher)
    tablename <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/GOanalysis/hm_clusters/hm_", comp, "_", clust, "_GOtable_", ont, ".tsv", sep="")
    write.table(allRes, file=tablename, row.names=FALSE, sep="\t")
    
    if(is.null(allRes_ext)==F){
      allRes_ext <- allRes_ext %>% dplyr::select(set, GO.ID, Term, Annotated, Significant, Expected, elimFisher, GeneID, GeneName, Prot_descr) %>%
        arrange(elimFisher, GO.ID, GeneID, Prot_descr)
      tablename <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/", threshold,"/GOanalysis/hm_clusters/hm_", comp, "_", clust, "_GOtable_", ont, "_ext", ".tsv", sep="")
      write.table(allRes_ext, file=tablename, row.names=FALSE, sep="\t")
    }  
  }
  GOdfs_list[[set]] <- GOdfs
  names(GOdfs_list)[set] <- comp
  
  FisherRes_list[[set]] <- FishRes
  names(FisherRes_list)[set] <- comp
  
  # end loop
}

head(GOdfs_list)
names(GOdfs_list)
names(FisherRes_list)
FisherRes_list[[1]]

# save data ####
save(GOdfs_list, FisherRes_list, file=paste("analysis/_RNAseq/_results/limma_GOanalysis_HMclusts_", ont, ".RData", sep=""))
# save a few data objects to use later so we don’t have to rerun everything
#-------------------------------------------------------------------------------------------------


sessionInfo() %>% capture.output(file="analysis/_RNAseq/_src/env_topGO.txt")
