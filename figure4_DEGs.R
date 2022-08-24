# Plot individual expression pattern for DEGs, z-scores of vst transformed expression
# Plot 4 panel figure; expression for each gene in each week

# load packages
library(tidyverse)
library(edgeR) # loads limma as a dependency
library(DESeq2)
library(cowplot)
theme_set(theme_cowplot()) #white background instead of grey -> don't load if want grey grid
behandeling <- c("dodgerblue2","firebrick3","black") 

threshold <- 0.01 #FDR threshold I want to use ####


# load data and results
load("analysis/_RNAseq/_results/limma_mod_allweeks_filtn2.RData") # all weeks analyzed together

head(dfilt$counts) # raw counts filtered by cutoff ~20/median(library size) in at least 2 replicates
head(dfilt$samples) # sample information incl. library sizes, normalisation factors

rm(res_list.3h, res_list.24h, results.3h, results.24h, v.24h, v.3h, cutoff, filtwhich.n2) # remove object not going to use

DEGs_annot <- read.table(file=paste("analysis/_RNAseq/_results/DEGs_annot_n2_p", threshold,".tsv", sep=""), header=T)
head(DEGs_annot)

hm_clust <- read.csv("analysis/_RNAseq/_results/DEGS_hm_clusters.csv")
head(hm_clust) # maybe useful here?


#-------------------------------------------------------
# Plot per gene expression patterns ####
#-------------------------------------------------------
head(dfilt$counts)
vst <- vst(dfilt$counts)
head(vst)

# use z-scores for expression value ####
z <- t(scale(t(vst), center=TRUE, scale=TRUE)) # take means of z-scores, otherwise individual points don't match mean line
range(z)
head(z)

z.ind <- t(z[rownames(z) %in% unique(DEGs_annot$GeneID),]) %>% as.data.frame
z.ind$Sample <- rownames(z.ind)
z.ind <- merge(z.ind, dfilt$samples[,c("Sample", "Trw_num", "Treatment","Timepoint")], by="Sample")
z.ind <- z.ind %>% mutate(Treat_week=paste("Week", Trw_num, sep=""),
                          Hour=ifelse(Timepoint=="BEF", "0", Timepoint),
                          Hour_num=ifelse(Timepoint=="BEF", 1, ifelse(Timepoint=="3", 2, 6)),
                          Treatment=ifelse(Treatment=="C", "Cold", ifelse(Treatment=="W", "Warm", "Ref 10C"))) %>%
  mutate(Treatment=factor(Treatment, levels=c("Cold", "Warm", "Ref 10C")), Hour=factor(Hour, levels=c("0", "3", "24")))
levels(z.ind$Treatment)
levels(z.ind$Hour)
unique(z.ind$Treat_week)
unique(z.ind$Hour_num) # all info per sample per gene

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
plotDat <- z.means %>%
  as.data.frame() %>%
  mutate(w2_CBEF=w2_BEF, w2_WBEF=w2_BEF, #add 0h point for C and W
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
plotDat <- plotDat %>%  mutate(Treatment=as.factor(ifelse(grepl("C", Treatment)==T, "Cold", ifelse(grepl("W", Treatment)==T, "Warm","Ref 10C"))),
                               Hour_num=as.numeric(ifelse(grepl("BEF", Hour)==T, 1, ifelse(grepl("3", Hour)==T, 2, 6))),
                               Hour=as.factor(ifelse(grepl("BEF", Hour)==T, "0", Hour)))
plotDat <- plotDat %>%  mutate(Treatment=factor(plotDat$Treatment, levels=c("Cold", "Warm", "Ref 10C")), Hour=factor(plotDat$Hour, levels=c("0", "3", "24"))) %>%
  arrange(Treatment) # use this order,otherwise 0h point colored red in graph, don't know why...
levels(plotDat$Treatment)
levels(plotDat$Hour)
unique(plotDat$Treat_week)
unique(plotDat$Hour_num) # mean z-scores per treatment-timepoint


# Plot Period ####
#-------------------------------------------------------
gene <- "MSTRG.9278"

Dat <- plotDat[,c("Treat_week", "Treatment", "Hour", "Hour_num", gene)]
colnames(Dat)[5] <- "y"
head(Dat)

ind.Dat <- z.ind[,c("Treat_week", "Treatment", "Hour", "Hour_num", gene)]
colnames(ind.Dat)[5] <- "y"
head(ind.Dat)

plot <- ggplot(Dat, aes(x=Hour_num, y=y, col=Treatment, shape=Treatment))+
  scale_colour_manual(values=c("dodgerblue2","firebrick3", "black"))+
  scale_shape_manual(values=c(19,17,18)) + # use different shapes to prevent overlap
  facet_wrap(~Treat_week)+
  #geom_point(data=ind.Dat, size=2.5, alpha=0.5)+ # no differences in sample timepoint
  geom_jitter(data=ind.Dat, size=2.5, alpha=0.5, width=0.3, height=0)+ # so jitter, but only in width, not height
  geom_point(size=4) +
  geom_line()+
  scale_x_continuous(breaks=c(1, 2, 6), labels=paste(levels(Dat$Hour), "h", sep=""))+
  scale_y_continuous(breaks=seq(-10,10, by=1))+
  labs(title= paste(gene, filter(DEGs_annot, GeneID==gene)$Prot_descr[1], sep=": "),
       x="Timepoint",
       y="Gene expression z-scores")+
  theme(strip.text = element_text(size=18), axis.text = element_text(size=15), plot.title = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_text(size=20, vjust=2.5), 
        legend.title = element_text(size=18, vjust=-1), legend.text = element_text(size=16))
plot

plot_name <- paste(paste("analysis/_RNAseq/_plots/mod_allweeks/limma/",threshold, "/DEGs/pretty/",gene,  sep=""), ".png", sep="")# name for file
#jpeg(plot_name, width=200, height=115, units="mm", res=300, quality=100) # quality of 100 means uncompressed
png(plot_name, width=200, height=115, units="mm", res=300)
print(plot)
dev.off()


# Plot candidate genes Table1 ####
#-------------------------------------------------------
genes <- c("MSTRG.6256", "MSTRG.19976", "MSTRG.30412", "MSTRG.39480", "OBRU01_205030", "MSTRG.37241", "OBRU01_214986")
annot <- c("pdfr", "orcokinin", "asator", "elys", "JNK interacting protein", "cenG1A", "alp4")

Dat <- plotDat[,c("Treat_week", "Treatment", "Hour", "Hour_num", genes)]
Dat <- Dat %>% pivot_longer(c(5:11), names_to="GeneID", values_to="y")
head(Dat)

ind.Dat <- z.ind[,c("Treat_week", "Treatment", "Hour", "Hour_num", genes)]
ind.Dat <- ind.Dat %>% pivot_longer(c(5:11), names_to="GeneID", values_to="y")
head(ind.Dat)

plots <- list()
for(g in 1:length(genes)){
  d <- filter(Dat, GeneID==genes[g])
  ind.d <- filter(ind.Dat, GeneID==genes[g])
  
  plot <- ggplot(d, aes(x=Hour_num, y=y, col=Treatment, shape=Treatment))+
    scale_colour_manual(values=c("dodgerblue2","firebrick3", "black"))+
    scale_shape_manual(values=c(19,17,18)) + # use different shapes to prevent overlap
    facet_wrap(~Treat_week)+
    #geom_point(data=ind.Dat, size=2.5, alpha=0.5)+ # no differences in sample timepoint
    geom_jitter(data=ind.d, size=2.5, alpha=0.5, width=0.3, height=0)+ # so jitter, but only in width, not height
    geom_point(size=4) +
    geom_line()+
    scale_x_continuous(breaks=c(1, 2, 6), labels=paste(levels(Dat$Hour), "h", sep=""))+
    scale_y_continuous(breaks=seq(-10,10, by=1))+
    labs(title= paste("\n", genes[g],": ",annot[g], sep=""),
         x="Timepoint",
         y="  ")+
    theme(strip.text = element_text(size=18), axis.text.x = element_blank(), axis.text.y = element_text(size=15),
          plot.title = element_text(size=15, face="bold", hjust = 0.5),
          axis.title.x = element_blank(), axis.title.y = element_text(size=30, vjust=2),
          legend.title = element_text(size=18, vjust=-1), legend.text = element_text(size=16),
          legend.position="none")
  plot
  
  if(g==3){
    plot <- plot + labs(y="Gene expression z-scores")
  }
  
  if(g==6 | g==7){
    plot <- plot + theme(axis.text.x = element_text(size=15))
    #if(g==7){plot <- plot + theme(legend.position="right")}
  }
  
  plots[[g]] <- plot
}

plot_t1 <- gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]],
                                   nrow = 4)

plot_name <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/",threshold, "/DEGs/pretty/Table1_genes.png",  sep="")# name for file
ggsave(plot=plot_t1, filename=plot_name, device="png", width=350, height=500, units="mm", dpi="print")


# Plot candidate genes Table2 ####
#-------------------------------------------------------
genes <- c("MSTRG.20613", "MSTRG.12645", "MSTRG.3087", "MSTRG.22817", "MSTRG.8284", "MSTRG.32147", 
           "MSTRG.13770", "MSTRG.11979", "MSTRG.25766", "MSTRG.11006", "MSTRG.10930")
annot <- c("tao", "fdl", "tey", "gld", "cyp6a13", "sp7", "atlastin", "nanos", "cycC", "psc", "polo")

Dat <- plotDat[,c("Treat_week", "Treatment", "Hour", "Hour_num", genes)]
Dat <- Dat %>% pivot_longer(c(5:15), names_to="GeneID", values_to="y")
head(Dat)

ind.Dat <- z.ind[,c("Treat_week", "Treatment", "Hour", "Hour_num", genes)]
ind.Dat <- ind.Dat %>% pivot_longer(c(5:15), names_to="GeneID", values_to="y")
head(ind.Dat)

plots <- list()
for(g in 1:length(genes)){
  d <- filter(Dat, GeneID==genes[g])
  ind.d <- filter(ind.Dat, GeneID==genes[g])
  
  plot <- ggplot(d, aes(x=Hour_num, y=y, col=Treatment, shape=Treatment))+
    scale_colour_manual(values=c("dodgerblue2","firebrick3", "black"))+
    scale_shape_manual(values=c(19,17,18)) + # use different shapes to prevent overlap
    facet_wrap(~Treat_week)+
    #geom_point(data=ind.Dat, size=2.5, alpha=0.5)+ # no differences in sample timepoint
    geom_jitter(data=ind.d, size=2.5, alpha=0.5, width=0.3, height=0)+ # so jitter, but only in width, not height
    geom_point(size=4) +
    geom_line()+
    scale_x_continuous(breaks=c(1, 2, 6), labels=paste(levels(Dat$Hour), "h", sep=""))+
    scale_y_continuous(breaks=seq(-10,10, by=1))+
    labs(title= paste("\n", genes[g],": ",annot[g], sep=""),
         x="Timepoint",
         y="  ")+
    theme(strip.text = element_text(size=18), axis.text.x = element_blank(), axis.text.y = element_text(size=15),
          plot.title = element_text(size=15, face="bold", hjust = 0.5),
          axis.title.x = element_blank(), axis.title.y = element_text(size=30, vjust=2),
          legend.title = element_text(size=18, vjust=-1), legend.text = element_text(size=16),
          legend.position="none")
  plot
  
  if(g==3 | g==9){
    plot <- plot + labs(y="Gene expression z-scores")
  }
  
  if(g==5| g==6 | g==10 | g==11){
    plot <- plot + theme(axis.text.x = element_text(size=15))
    #if(g==7){plot <- plot + theme(legend.position="right")}
  }
  
  plots[[g]] <- plot
}

plot_t2 <- gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]],
                                   plots[[7]], plots[[8]], plots[[9]], plots[[10]], plots[[11]],
                                   nrow = 6)

plot_name <- paste("analysis/_RNAseq/_plots/mod_allweeks/limma/",threshold, "/DEGs/pretty/Table2_genes.png",  sep="")# name for file
ggsave(plot=plot_t2, filename=plot_name, device="png", width=350, height=680, units="mm", dpi="print")
