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

# Define colMap, ref. https://github.com/Bioconductor-mirror/topGO/blob/master/vignettes/topGO.Rnw
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}


# Load results ####
#-------------------------------------------------------------------------------------------------
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
geneUniverse <- SiteInfo$gene
head(geneUniverse)

# gene level annotation
annot <- read.table("analysis/_RNAseq/_annotation/genesAnnot_wblastp_rows_exp.tsv", header=T)
str(annot)
head(annot[,c("GeneID", "TranscriptID", "Evalue")])# ordered by GeneID and Evalue, most significant hit for each GeneID first
length(unique(annot$GeneID)) # this should be 29113
table(is.na(annot$Evalue)) #4142 genes without hit

# Load DEGS with gene level annotation
# DEGs_annot <- read_tsv("analysis/_RNAseq/_results/DEGs_annot_n2_p0.01.tsv", col_names=T) # table of all DEGs
# length(unique(DEGs_annot$GeneID))
# DEGs_IDs <- unique(DEGs_annot$GeneID) # store GeneIDs of DEGs


#-------------------------------------------------------------------------------------------------
#### ENRICHMENT: Overrepresentation analysis
#-------------------------------------------------------------------------------------------------
table(SiteInfo.sign$Module) # per Module

# Explore topGO environments
BPterms <- ls(GOBPTerm) # list of GO terms from BP ontology
head(BPterms)

# Choose ontology
ont <- "MF" #ontologies Biological Processes (BP), Molecular Function (MF)

# create lists to store results in ####
GOdfs_list <- as.list(unique(SiteInfo.sign$Module))
names(GOdfs_list) <- unique(SiteInfo.sign$Module)

FisherRes_list <- as.list(unique(SiteInfo.sign$Module))
names(FisherRes_list) <- unique(SiteInfo.sign$Module)


# LOOP THROUGH MODULES = 4x ####
for (mod in 1:length(GOdfs_list)){
  module <- names(GOdfs_list)[mod]
  res <- filter(SiteInfo.sign, Module==module)$gene %>% unique
  
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
                      description = module, # optional 
                      ontology = ont, # GO graph
                      allGenes = geneList, # gene universe as logical DEG yes or no
                      nodeSize = 5, #prune GO hierarchy from terms with less than X annotated genes (5-10 recommended)
                      annot = annFUN.gene2GO, gene2GO=geneID2GO) #function used to extract gene-to-GO mappings from affyLib object
  GOdfs_list[[mod]] <- sampleGOdata
  GOdfs_list[[mod]] # see summary of object
  
  
  # Overrepresentation analysis ####
  #-------------------------------------------------------------------------------------------------
  # testing the over-representation of GO terms within the group of differentially expressed genes
  
  resultFisher <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher") 
  # resultFisher.cl <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher") 
  # method=classic, each GO category is tested independently
  # not adjusted for multiple testing
  # method=elim, more conservative
  FisherRes_list[[mod]] <- resultFisher # see summary of object
  FisherRes_list[[mod]]

  
  # DOWN STREAM ANALYSIS TOOLS ####
  #-------------------------------------------------------------------------------------------------
  signNodes <- table(resultFisher@score<0.01)["TRUE"] # number of significant GO terms P<0.01
  if(is.na(signNodes)==TRUE) {signNod <- 1
  } else signNod <- signNodes
  
  # gather results in a table ####
  allRes <- GenTable(sampleGOdata, elimFisher = resultFisher,
                     #classicFisher = resultFisher.cl,
                     orderBy = "elimFisher", ranksOf = "elimFisher", topNodes = signNod) # only keep top 10 results of elimKS results
  allRes$Module <- module
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
    allRes_ext <- merge(allRes_ext, unique(annot[,c("GeneID", "GeneName", "Prot_descr")]), by="GeneID") # add annotation to the genes
    head(allRes_ext)
    
  } else{
    allRes_ext <- character(0)
  }
  
  # Visualize distribution significant GO terms over GO graph ####
  plotname <- paste("analysis/_RNAseq/_plots/descriptive/WGCNA/CutHeight0.40/GOanalysis/", module, "_GOgraph_", ont, ".pdf", sep="")
  pdf(plotname, width = 16, height = 9, family="sans", colormodel="srgb")
  showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = signNod, useInfo = 'all') # based on 5 most sign GO terms of results
  dev.off()
  
  # Store results
  head(allRes)
  head(allRes_ext)
  
  allRes <- subset(allRes, elimFisher<=0.01) %>% arrange(elimFisher)
  tablename <- paste("analysis/_RNAseq/_plots/descriptive/WGCNA/CutHeight0.40/GOanalysis/", module, "_GOtable_", ont, ".tsv", sep="")
  write.table(allRes, file=tablename, row.names=FALSE, sep="\t")
  
  if(is.null(allRes_ext$elimFisher)==F){
    allRes_ext <- allRes_ext %>% dplyr::select(Module, GO.ID, Term, Annotated, Significant, Expected, elimFisher, GeneID, GeneName, Prot_descr) %>%
      arrange(elimFisher, GO.ID, GeneID, Prot_descr)
    tablename <- paste("analysis/_RNAseq/_plots/descriptive/WGCNA/CutHeight0.40/GOanalysis/", module, "_GOtable_", ont, "_ext", ".tsv", sep="")
    write.table(allRes_ext, file=tablename, row.names=FALSE, sep="\t")
    }
  
  # end loop ####
}

names(GOdfs_list)

# save data ####
save(GOdfs_list, FisherRes_list, file=paste("analysis/_RNAseq/_results/WGCNA_GOanalysis_", ont, ".RData", sep=""))
# save a few data objects to use later so don’t have to rerun everything
#-------------------------------------------------------------------------------------------------

sessionInfo()