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


# Load DEGS with gene level annotation ####
DEGs_annot <- read_tsv("analysis/_RNAseq/_results/DEGs_annot_n2_p0.01.tsv", col_names=T) # table of all DEGs


#-------------------------------------------------------------------------------------------------
# Cluster GO terms with ViSEAGO ####
#-------------------------------------------------------------------------------------------------

# Reformat geneID2GO file ####
#---------------------------------------------------------------------------------------
geneID2GO.txt <- read.table(file="analysis/_RNAseq/_annotation/TransfExp_WM_gene2GO_v2.txt", sep="\t")
names(geneID2GO.txt) <- c("GeneID", "GO")
head(geneID2GO.txt)

# format of text file needed:
# taxid       gene_id	gene_symbol	GOID	      evidence
# myspecies1	idA	    geneA	      GO:0000001	IEA=Inferred from Electronic Annotation

GOtbl <- geneID2GO.txt %>% separate(GO, sep = ",", into=c(sprintf("GO%03d", seq(1:135)))) # 1 row with 135 GO terms
GOtbl <- GOtbl %>% gather(key="key", value="GOterm", colnames(GOtbl)[c(2:136)]) %>%
  arrange(GeneID, key)
GOtbl <- GOtbl %>% mutate(taxid="O.Brumata", gene_id=GeneID, gene_symbol=GeneID, GOID=GOterm, evidence="IEA") %>%
  dplyr::select(taxid, gene_id, gene_symbol, GOID, evidence)

table(is.na(GOtbl$GOID)) # make sure only NAs are kept for genes with 0 GO_IDs
GOtbl <- GOtbl[!duplicated(GOtbl[,c("gene_id", "GOID")]),]
length(unique(GOtbl$gene_id))

table(grepl(" GO:",GOtbl$GOID)) # remove spaces before GO term in Notepad++
head(GOtbl)
# write.table(GOtbl, file="analysis/_RNAseq/_annotation/TransfExp_WM_gene2GO_rows_v2.txt", row.names=F, sep="\t")

GOtbl <- read.table(file="analysis/_RNAseq/_annotation/TransfExp_WM_gene2GO_rows_v2.txt", sep="\t", header=T)
table(grepl(" GO:",GOtbl$GOID)) # spaces are gone
head(GOtbl)

viseago_gos <- select(GO.db,columns=columns(GO.db),keys=keys(GO.db))
head(viseago_gos)

table(filter(GOtbl, !is.na(GOID))$GOID %in% viseago_gos$GOID)
check <- subset(unique(filter(GOtbl, !is.na(GOID))), !(GOID %in% viseago_gos$GOID))
head(check)
length(unique(check$GOID)) #10 GO terms not included in ViSEAGO GO database ####
unique(check$GOID) # don't know why they are not included, but the GO terms do exist when I google them
length(unique(filter(GOtbl, !is.na(GOID))$GOID)) # from 7736, so that's <0.01% excluded then

GOtbl_sub <- subset(GOtbl, GOID %in% viseago_gos$GOID) # can't have NAs in here
table(GOtbl_sub$GOID %in% viseago_gos$GOID)
# write.table(GOtbl_sub, file="analysis/_RNAseq/_annotation/TransfExp_WM_gene2GO_rows_sub.txt", row.names=F, sep="\t")


# Load annotation format needed for ViSEAGO ####
#---------------------------------------------------------------------------------------
geneID2GOb <-ViSEAGO::Custom2GO(file="analysis/_RNAseq/_annotation/TransfExp_WM_gene2GO_rows_sub.txt") 
str(geneID2GOb)

myGENE2GO<-ViSEAGO::annotate( # load GO annotations from Custom annotation table
  id="O.Brumata",
  object=geneID2GOb
)
class(myGENE2GO)
myGENE2GO

rm(geneID2GOb) # keep environment clean from unnecessary files


# Choose ontology ####
#---------------------------------------------------------------------------------------
ont <- "MF" # MF or BP


# Load TopGO objects and results for further analysis ####
#---------------------------------------------------------------------------------------
# load(paste("analysis/_RNAseq/_results/limma_GOanalysis_sets_", ont, ".RData", sep="")) # enrichment on contrast sets Overall_Temp, Within_Week, Between_Week
# load(paste("analysis/_RNAseq/_results/limma_GOanalysis_overlap_", ont, ".RData", sep="")) # enrichment on overlap Y/N with WGCNA analysis
# load(paste("analysis/_RNAseq/_results/limma_GOanalysis_weekspecific_", ont, ".RData", sep="")) # enrichment on Within_week contrasts per week
# load(paste("analysis/_RNAseq/_results/limma_GOanalysis_BetweenWeek_", ont, ".RData", sep="")) # enrichment on 2 Between_week contrasts
# load(paste("analysis/_RNAseq/_results/limma_GOanalysis_HMclusts_", ont, ".RData", sep="")) # enrichment on identified HM clusters
load(paste("analysis/_RNAseq/_results/WGCNA_GOanalysis_", ont, ".RData", sep="")) # enrichment per module

names(GOdfs_list)
head(GOdfs_list[[1]]) # HM clusters are lists within list
GOset <- "WGCNA.Modules" # Make sure to set this variable!! choose from: "ContrastSets", "Overlap", "WeekSpecific", "BetweenWeek", "HMclusts", "WGCNA.Modules"

for(item in 1:length(GOdfs_list)){ # ONLY RUN FULL LOOP FOR HM CLUSTERS ####
  list <- GOdfs_list[[item]]
  ResList <- FisherRes_list[[item]]
  
  # Need to recreate GOdfs via ViSEAGO... ####
  #---------------------------------------------------------------------------------------
  GOdfs_list2 <- GOdfs_list #IF NOT HM CLUSTERS, START HERE ####
  if(is.null(names(GOdfs_list2))==TRUE){names(GOdfs_list2) <- as.character(seq(1, length(GOdfs_list2), by=1)) }
  
  for (set in 1:length(names(GOdfs_list2))){
    #name <- names(GOdfs_list)[set]
    GOdf <- GOdfs_list2[[set]]
    
    # pull info from old GOdf and recreate GOdf
    geneList <- GOdf@allGenes
    
    # get DEG_IDs
    # selection <- geneList[geneList %in% unique(filter(DEGs_annot, Comparison==name)$GeneID)]
    selection <- GOdf@allScores
    names(selection) <- geneList
    selection <- names(selection[selection==1])
    
    GOdf_2 <-ViSEAGO::create_topGOdata(
      geneSel=selection,
      allGenes=geneList,
      geneList=NULL, 
      gene2GO=myGENE2GO, # requires specific class object gene2GO-class object created by annotate
      ont=ont,
      nodeSize=5
    )
    
    GOdfs_list2[[set]] <- GOdf_2
    
  }
  names(GOdfs_list2)
  
  rm(GOdf, GOdf_2, geneList, selection, set) # clean up files no longer needed
  
  
  # Merge all GO terms from different contrasts into one dataframe and cluster them ####
  #---------------------------------------------------------------------------------------
  if(GOset=="ContrastSets"){try({ # if contrast sets
    GO_Overall_Temp <-GOdfs_list2$Overall_Temp
    Res_Overall_Temp <-FisherRes_list$Overall_Temp
    
    GO_Within_week <-GOdfs_list2$Within_week
    Res_Within_week <-FisherRes_list$Within_week
    
    GO_Between_week <-GOdfs_list2$Between_week
    Res_Between_week <-FisherRes_list$Between_week
    
    sResults<-ViSEAGO::merge_enrich_terms(
      Input=list(
        Overall_Temp=c("GO_Overall_Temp","Res_Overall_Temp"),
        Within_week=c("GO_Within_week","Res_Within_week"),
        Between_week=c("GO_Between_week","Res_Between_week")
      )
    )
    rm(GO_Overall_Temp, Res_Overall_Temp, GO_Within_week, Res_Within_week, GO_Between_week, Res_Between_week)
  })
  }
  
  if(GOset=="Overlap"){try({ # if overlap Y/N WGCNA
    GO_Y <-GOdfs_list2$overlapping
    Res_Y <-FisherRes_list$overlapping
    
    GO_N <-GOdfs_list2$not_overlapping
    Res_N <-FisherRes_list$not_overlapping
    
    sResults<-ViSEAGO::merge_enrich_terms(
      Input=list(
        Overlapping=c("GO_Y","Res_Y"),
        Not_Overlapping=c("GO_N","Res_N")
      )
    )
    rm(GO_Y, Res_Y, GO_N, Res_N)
  })
  }
  
  if(GOset=="WeekSpecific"){try({ # if week-specific sets
    GO_w2 <-GOdfs_list2$Week2
    Res_w2 <-FisherRes_list$Week2
    
    GO_w4 <-GOdfs_list2$Week4
    Res_w4 <-FisherRes_list$Week4
    
    GO_w6 <-GOdfs_list2$Week6
    Res_w6 <-FisherRes_list$Week6
    
    GO_w8 <-GOdfs_list2$Week8
    Res_w8 <-FisherRes_list$Week8
    
    sResults<-ViSEAGO::merge_enrich_terms(
      Input=list(
        Week2=c("GO_w2","Res_w2"),
        Week4=c("GO_w4","Res_w4"),
        Week6=c("GO_w6","Res_w6"),
        Week8=c("GO_w8","Res_w8")
      )
    )
    rm(GO_w2, Res_w2, GO_w4, Res_w4, GO_w6, Res_w6, GO_w8, Res_w8)
  })
  }
  
  if(GOset=="BetweenWeek"){try({ # if between-week specific sets
    GO_w4vs6 <-GOdfs_list2$w4vs6
    Res_w4vs6 <-FisherRes_list$w4vs6
    
    GO_w6vs8 <-GOdfs_list2$w6vs8
    Res_w6vs8 <-FisherRes_list$w6vs8
    
    sResults<-ViSEAGO::merge_enrich_terms(
      Input=list(
        w4vs6=c("GO_w4vs6","Res_w4vs6"),
        w6vs8=c("GO_w6vs8","Res_w6vs8")
      )
    )
    rm(GO_w4vs6, Res_w4vs6,GO_w6vs8, Res_w6vs8)
  })
  }
  
  if(GOset=="WGCNA.Modules"){try({ # if WGCNA modules
    GO_mod3 <-GOdfs_list2$Mod3
    Res_mod3 <-FisherRes_list$Mod3
    
    GO_mod8 <-GOdfs_list2$Mod8
    Res_mod8 <-FisherRes_list$Mod8
    
    GO_mod10 <-GOdfs_list2$Mod10
    Res_mod10 <-FisherRes_list$Mod10
    
    GO_mod13 <-GOdfs_list2$Mod13
    Res_mod13 <-FisherRes_list$Mod13
    
    sResults<-ViSEAGO::merge_enrich_terms(
      Input=list(
        Module3=c("GO_mod3","Res_mod3"),
        Module8=c("GO_mod8","Res_mod8"),
        Module10=c("GO_mod10","Res_mod10"),
        Module13=c("GO_mod13","Res_mod13")
      )
    )
    rm(GO_mod3, Res_mod3, GO_mod8, Res_mod8, GO_mod10, Res_mod10, GO_mod13, Res_mod13)
  })
  }
  
  if(GOset=="HMclusts"){try({ # if HM clusters, can be 2 to 5 clusters
    GO_1 <-GOdfs_list2$`1`
    Res_1 <-ResList[[1]]
    
    GO_2 <-GOdfs_list2$`2`
    Res_2 <-ResList[[2]]
    
    if(length(GOdfs_list2)>2){
      GO_3 <-GOdfs_list2$`3`
      Res_3 <-ResList[[3]]}
    
    if(length(GOdfs_list2)>3){
      GO_4 <-GOdfs_list2$`4`
      Res_4 <-ResList[[4]]}
    
    if(length(GOdfs_list2)>4){
      GO_5 <-GOdfs_list2$`5`
      Res_5 <-ResList[[5]]}
    
    if(length(GOdfs_list2)==2){
      sResults<-ViSEAGO::merge_enrich_terms(
      Input=list(
        Cluster1=c("GO_1","Res_1"),
        Cluster2=c("GO_2","Res_2")
      )
    )
    }
    
    if(length(GOdfs_list2)==3){sResults<-ViSEAGO::merge_enrich_terms(
      Input=list(
        Cluster1=c("GO_1","Res_1"),
        Cluster2=c("GO_2","Res_2"),
        Cluster3=c("GO_3","Res_3")
      )
    )
    }
    
    if(length(GOdfs_list2)==4){sResults<-ViSEAGO::merge_enrich_terms(
      Input=list(
        Cluster1=c("GO_1","Res_1"),
        Cluster2=c("GO_2","Res_2"),
        Cluster3=c("GO_3","Res_3"),
        Cluster4=c("GO_4","Res_4")
      )
    )
    }
    
    if(length(GOdfs_list2)==5){sResults<-ViSEAGO::merge_enrich_terms(
      Input=list(
        Cluster1=c("GO_1","Res_1"),
        Cluster2=c("GO_2","Res_2"),
        Cluster3=c("GO_3","Res_3"),
        Cluster4=c("GO_4","Res_4"),
        Cluster5=c("GO_5","Res_5")
      )
    )
    }
    
    #rm(GO_Overall_Temp, Res_Overall_Temp, GO_Within_week, Res_Within_week, GO_Between_week, Res_Between_week)
  })
  }
  
  
  # Cluster significant GO terms based on their semantic similarity ####
  myGOs<-ViSEAGO::build_GO_SS( # initialize
    gene2GO=myGENE2GO,
    enrich_GO_terms=sResults
  )
  
  myGOs<-ViSEAGO::compute_SS_distances(# compute all available Semantic Similarity (SS) measures
    myGOs,
    distance="Wang"
  )
  
  rm(sResults)
  
  
  # Visualize ####
  #---------------------------------------------------------------------------------------
  # ViSEAGO::GOcount(sResults) # Visualize how many GO terms from total number enriched for each contrast
  # ViSEAGO::Upset(sResults) # visualize if enriched GO terms are contrast specific or shared between contrasts
  
  try(
    {Wang_clusters<-ViSEAGO::GOterms_heatmap(# GOterms heatmap with the default parameters
      myGOs,
      showIC=FALSE,
      showGOlabels=TRUE,
      GO.tree=list(
        tree=list(
          distance="Wang",
          aggreg.method="ward.D2"
        ),
        cut=list(
          dynamic=list(
            pamStage=TRUE,
            pamRespectsDendro=TRUE,
            deepSplit=2,
            minClusterSize =2
          )
        )
      ),
      samples.tree=NULL
    )
    })
  rm(myGOs)
  
  # Access and change visual attributes of hm ####
  if(GOset=="HMclusts"){
    title <- paste(ont," GO Overrepresentation Results ", GOset, " ", names(GOdfs_list)[item],"\n Clustered with Wang distance", sep="")
  }else{ title <- paste(ont," GO Overrepresentation Results ", GOset,"\n Clustered with Wang distance", sep="")}
  
  Wang_clusters@heatmap$GOterms$x$layoutAttrs[[2]]$title <- title
  Wang_clusters@heatmap$GOterms$x$layoutAttrs[[2]]$font$size <- 20
  Wang_clusters@heatmap$GOterms$x$layoutAttrs[[3]]$yaxis$tickfont$size <- 12 
  Wang_clusters@heatmap$GOterms$x$layoutAttrs[[4]]$xaxis$tickfont$size <- 18
  
  
  # Display heatmap ####
  try({
    
  if(GOset=="HMclusts"){
      plot_name <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/0.01/GOanalysis/ViSEAGO/_hm_clusters/GO_hm_", names(GOdfs_list)[item],"_", ont,"_heatmap", sep="")
    }else{ plot_name <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/0.01/GOanalysis/ViSEAGO/GO_", GOset,"_", ont,"_heatmap", sep="")}
    
    hm <- ViSEAGO::show_heatmap( # Display the clusters-heatmap
      Wang_clusters,
      "GOterms"
    )
    saveWidget(hm, paste(plot_name, ".html", sep=""))
    webshot(paste(plot_name, ".html", sep=""), paste(plot_name, ".png", sep=""), vwidth = 1500, vheight = 1000)
    
  })
  rm(hm, plot_name, title)
  
  # plot(Wang_clusters@dendrograms$GO)
  # abline(h=1.7, col="red", lty=2, lwd=2)
  # #new.clusters <- as.data.frame(sort(cutree(Wang_clusters@terms_dist$Wang, h=1.7))) # can't reassign clusters this way
  
  
  # Display the clusters-heatmap table ####
  try({
    
    if(GOset=="HMclusts"){
      tbl_name <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/0.01/GOanalysis/ViSEAGO/_hm_clusters/GO_hm_", names(GOdfs_list)[item],"_", ont,"_hmtable", sep="")
    }else{ tbl_name <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/0.01/GOanalysis/ViSEAGO/GO_", GOset,"_", ont,"_hmtable", sep="")}
    
    
    tbl <- ViSEAGO::show_table(Wang_clusters)
    saveWidget(tbl, paste(tbl_name, ".html", sep=""))
    
    tbl <- tbl$x$data #extract the table with values, then can just save it as .csv or .tsv
    tbl$GOID <- gsub(".*term/(GO:\\d+).*","\\1",tbl$`GO ID`)
    tbl <- tbl[,c(2,3, length(colnames(tbl)),5:(length(colnames(tbl))-1))]
    str(tbl)
    write.csv(tbl, file=paste(tbl_name, ".csv", sep=""), row.names=F)
    #filter(tbl, `GO cluster`==1)$GOID %>% head
  })
  rm(tbl, tbl_name)
  
  
  # # display MDSplot colored by assigned cluster ####
  # try(
  #   {plot_name <- paste("analysis/_RNAseq/_plots/", week, "_mod/GOanalysis/", week, "_GO_", ont,"_MDS", sep="")
  #   mds <- ViSEAGO::MDSplot(
  #     Wang_clusters_wardD2,
  #     "GOterms")
  #   saveWidget(mds, paste(plot_name, ".html", sep=""))
  #   webshot(paste(plot_name, ".html", sep=""), paste(plot_name, ".png", sep=""), vwidth = 1200, vheight = 700)}
  # )
  
  
  # Save Wang_clusters objects ####
  #---------------------------------------------------------------------------------------
  # Can use them to still try and cut the trees into different clusters if I want
  if(GOset=="HMclusts"){
    save(Wang_clusters, file=paste("analysis/_RNAseq/_results/limma_ViSEAGO_", ont, "_hm_", names(GOdfs_list)[item],".RData", sep=""))
  }
  if(GOset=="WGCNA.Modules"){
    save(Wang_clusters, file=paste("analysis/_RNAseq/_results/WGCNA_ViSEAGO_", ont,".RData", sep=""))
  }else{ save(Wang_clusters, file=paste("analysis/_RNAseq/_results/limma_ViSEAGO_", ont, "_", GOset,".RData", sep=""))}
  
     
}


#---------------------------------------------------------------------------------------
# Get parent terms with ViSEAGO (only possible if >2 clusters) ####
#---------------------------------------------------------------------------------------
ont <- "MF" # MF or BP
GOset <- "BetweenWeek" # Make sure to set this variable!! choose from: "ContrastSets", "Overlap", "WeekSpecific", "BetweenWeek", "HMclusts", "WGCNA.Modules"

# for HM clusters  
set <- "Overall_Temp" #or Within_week or Between_week
mod <- "24h" #or 3h

# load(paste("analysis/_RNAseq/_results/limma_ViSEAGO_", ont, "_hm_", set, "_", mod, ".RData", sep="")) # hm clusters
# load(paste("analysis/_RNAseq/_results/WGCNA_ViSEAGO_", ont,".RData", sep="")) # WGCNA results
load(paste("analysis/_RNAseq/_results/limma_ViSEAGO_", ont, "_", GOset,".RData", sep="")) # rest

hm <- ViSEAGO::show_heatmap( # Display the clusters-heatmap
  Wang_clusters,
  "GOterms"
)
hm

# calculate semantic similarites between clusters of GO terms
Wang_clusters <-ViSEAGO::compute_SS_distances(
  Wang_clusters,
  distance="BMA" # or use (combination of) "max", "avg","rcmax"
)

# GOclusters heatmap
Wang_clusters<-ViSEAGO::GOclusters_heatmap(
  Wang_clusters,
  tree=list(
    distance="BMA",
    aggreg.method="ward.D2"
  )
)

clusthm <- ViSEAGO::show_heatmap(
  Wang_clusters,
  "GOclusters"
)
clusthm
# "definition of the cluster of GO terms corresponds to the first common GO term ancestor 
# followed by the cluster label in brackets"

# Get GO clusters defining GO term:
BMAres <- Wang_clusters@clusters_dist$BMA %>% as.matrix
rownames(BMAres)

parents <- data.frame(GO.cluster=gsub("(\\d+)_GO:\\d+_.*","\\1",rownames(BMAres)),
                      parent_GOID=gsub("\\d+_(GO:\\d+)_.*","\\1",rownames(BMAres)), 
                      parent_term=gsub("\\d+_GO:\\d+_(.*)","\\1",rownames(BMAres)))
parents

if(GOset=="HMclusts"){
  tblname <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/0.01/GOanalysis/ViSEAGO/_hm_clusters/GO_hm_", set,"_", mod, "_", ont,"_hmtable.csv", sep="")
}else{ tblname <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/0.01/GOanalysis/ViSEAGO/GO_", GOset,"_", ont,"_hmtable.csv", sep="")}
tbl <- read.csv(tblname)
head(tbl)

tbl1 <- merge(tbl, parents, by="GO.cluster")
str(tbl1)

# write.csv(tbl1, file=tblname, row.names=F)


sessionInfo() %>% capture.output(file="analysis/_RNAseq/_src/env_ViSEAGO.txt")
