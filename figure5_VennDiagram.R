# load packages
library(tidyverse)
library(ggVennDiagram)
library(cowplot)
theme_set(theme_cowplot()) #white background instead of grey -> don't load if want grey grid
behandeling <- c("dodgerblue2","firebrick3","black") 

## For tweeking aesthetics ggVennDiagram, see: https://venn.bio-spring.top/using-ggvenndiagram


# Load data ####

# DEGs annotated
DEGs_annot <- read.table(file="analysis/_RNAseq/_results/DEGs_annot_n2_p0.01.tsv", header=T, sep="\t")

# WGCNA analysis: description 10C samples
load("analysis/_RNAseq/_results/WGCNA_Results_CoExpr10C.Rdata")
head(modNames) # ANOVA results for each module
head(SiteInfo.sign) # genes assigned to significant modules with significant membership 
head(SiteInfo) # results for all genes


# Venn diagram: Overlap analyses ####
#-------------------------------------------------------
overview <- list(unique(SiteInfo.sign$gene), unique(DEGs_annot$GeneID))
names(overview) <- c("WGCNA", "DE")
# names(overview) <- c("WGCNA analysis", "DE analysis")

# Use ggVennDiagram wrapper:
venn1 <- ggVennDiagram::ggVennDiagram(overview, # needs to be a list, with each item = vector of genes names
                                     set_size=7, 
                                     label="count", label_alpha = 0, label_size=6)+
  scale_fill_gradient(low="white",high ="white") + # "mediumorchid4")
  theme(legend.position="none", text = element_text(vjust=6)) + scale_y_continuous(expand = expansion(mult = .1))
#png(paste("analysis/_RNAseq/_plots/Venn_overlap_analyses.png", sep=""), width=100, height=60, units="mm", res=300)
print(venn1)
#dev.off()

# Tweek all aesthetics yourself: 
venn <- Venn(overview)
data <- process_data(venn)

venn2 <- ggplot() +
  geom_sf(data= venn_region(data), fill="white") +
  geom_sf(data=venn_setedge(data), aes(color=id), size=1, show.legend=FALSE) +
  geom_sf_text(data=venn_setlabel(data), aes(label=name), size=6, vjust=-0.2) +
  geom_sf_text(data=venn_region(data), aes(label=count), size=5, alpha=1) +
  theme_void() + # get rid of x and y axes
  scale_y_continuous(expand = expansion(mult = .15))
venn2

png(paste("analysis/_RNAseq/_plots/Venn_overlap_analyses.png", sep=""), width=100, height=60, units="mm", res=300)
print(venn2)
dev.off()
