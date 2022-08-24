# Open R project in top folder

# load packages
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot()) #white background instead of grey -> don't load if want grey grid

ME <- read.csv(file="analysis/_RNAseq/_results/WGCNA_0.4_MEs_to_plot.csv") # module eigen genes for CutHeight0.40
ME <- ME %>% mutate(Tube.n=as.character(Tube.n), Tube=as.character(Tube))
head(ME)
str(ME)

sign <- c("MEred_Mod10", "MEpink_Mod8", "MEyellow_Mod13", "MEbrown_Mod3")

# Plot all modules' eigen genes ####
#------------------------------------------------------------------------------
ME_tidy <- ME %>% pivot_longer(c(2:15), names_to="Module", values_to="PC1") %>%
  mutate(color=gsub("ME(\\w+)_Mod\\d+","\\1",Module), modul=gsub("ME\\w+_Mod(\\d+)","Module \\1",Module)) %>%
  mutate(color=gsub("ME(\\w+)_\\d+","\\1",color), modul=gsub("ME\\w+_(\\d+)","Module \\1",modul))
ME_tidy$sign <- ifelse(ME_tidy$Module %in% sign, "yes", "no")
ME_tidy$modul <- ifelse(ME_tidy$sign=="yes", paste(ME_tidy$modul, "*", sep=""), ME_tidy$modul)
ME_tidy$mod_num <- ifelse(ME_tidy$sign=="no", gsub("Module (\\d+)","\\1",ME_tidy$modul),
                          gsub("Module (\\d+).","\\1",ME_tidy$modul)) %>% as.numeric()
head(ME_tidy)

cols <- unique(ME_tidy$color)
plot_all <- ggplot(data=ME_tidy, aes(x=Trw_num, y=PC1, fill=color))+
  scale_fill_manual(values=cols) +
  facet_wrap(~modul, ncol=4)+
  scale_size_manual(values=c(4,5,6))+
  geom_jitter(aes(size=Tube.n), shape=21, alpha=0.7, height=0, width=0.3, stroke=1.5)+
  #geom_smooth()+
  labs(y="Module PC1", size="Week-specific\nClutch") + #, x=paste("Module", gsub("ME\\w+_Mod(\\d+)","\\1",module), sep=" ")) +
  coord_cartesian(ylim=c(-0.6, 0.6), xlim=c(1,9)) + scale_y_continuous(breaks=seq(-1,1, by=0.2))+
  scale_x_continuous(breaks=seq(2,8, by=2), labels=c("Week2", "Week4", "Week6", "Week8"))+
  guides(fill = "none")+
  theme(strip.text= element_text(size=20, face="bold"),
        #legend.position="none", 
        legend.title=element_text(size=20, hjust=0.5), legend.text=element_text(size=20),
        plot.tag=element_text(size=25, face="plain", hjust=0.25), axis.title.y = element_text(size=25, vjust=2), 
        axis.title.x = element_blank(), axis.text.x=element_text(size=15), axis.text.y=element_text(size=20))
plot_all

# ggsave(plot=plot_all, filename="analysis/_RNAseq/_plots/descriptive/WGCNA/CutHeight0.40/all_modules.jpg", device="jpeg", width=350, height=300, units="mm", dpi="print")
ggsave(plot=plot_all, filename="analysis/_RNAseq/_plots/descriptive/WGCNA/CutHeight0.40/all_modules_legend.jpg", device="jpeg", width=350, height=300, units="mm", dpi="print")


# Plot module eigen genes for significant modules
#------------------------------------------------------------------------------
# Significant modules: Mod3, Mod8, Mod10, Mod13
# Pairs: Mod10+Mod8; Mod13+Mod3

ME_sign <- ME[,c("Sample", "Trw_num", "Tube.n", "MEred_Mod10", "MEpink_Mod8", "MEyellow_Mod13", "MEbrown_Mod3")]
head(ME_sign)

plots <- list()

for(mod in c(1:4)){ # 4 significant modules in total
  module <- colnames(ME_sign)[mod+3]
  dat <- ME_sign[,c(module, "Sample", "Trw_num", "Tube.n")]
  colnames(dat)[1] <- "Module"
  color <- gsub("ME(\\w+)_Mod\\d+","\\1",module)
  
  tag <- ifelse(mod==1 | mod==3, "(a)", "(b)") # for two figures, each with 2 panels
  nrgenes <- ifelse(module=="MEred_Mod10", "n=460",ifelse(module=="MEpink_Mod8","n=208", ifelse(module=="MEyellow_Mod13", "n=3280", "n=3497")))
  
  plot <- ggplot(data=dat, aes(x=Trw_num, y=Module))+
    labs(tag= tag)+
    scale_size_manual(values=c(4,5,6))+
    geom_jitter(aes(size=Tube.n), shape=21, fill=color, alpha=0.7, height=0, width=0.3, stroke=1.5)+
    #geom_smooth()+
    labs(y="Module PC1") + #, x=paste("Module", gsub("ME\\w+_Mod(\\d+)","\\1",module), sep=" ")) +
    coord_cartesian(ylim=c(-0.4, 0.5)) + scale_y_continuous(breaks=seq(-1,1, by=0.1))+
    scale_x_continuous(breaks=seq(2,8, by=2), labels=c("Week2", "Week4", "Week6", "Week8"))+
    guides(shape = "none")+
    annotate(geom="text", label=nrgenes, size=7.5, x=7.6, y=0.5)+
    theme(legend.position="none", plot.tag=element_text(size=25, face="plain", hjust=0.25), axis.title.y = element_text(size=25, vjust=2), 
          axis.title.x = element_blank(), axis.text.x=element_text(size=20), axis.text.y=element_text(size=20)) #axis.title.x = element_text(size=30, face="bold", vjust=-0.5)
  
  if(tag=="(b)"){
    plot <- plot + theme(axis.title.y=element_text(colour="white"))
  }
  
  print(plot)
  
  plots[[mod]] <- plot
}

# Pairs: Mod10+Mod8; Mod13+Mod3
pair1 <- gridExtra::grid.arrange(plots[[1]], plots[[2]],
                    nrow = 1)
pair2 <- gridExtra::grid.arrange(plots[[3]], plots[[4]],
                                 nrow = 1)

ggsave(plot=pair1, filename="analysis/_RNAseq/_plots/descriptive/WGCNA/CutHeight0.40/SignMod_pair1_Mod10_8.jpg", device="jpeg", width=300, height=150, units="mm", dpi="print")
ggsave(plot=pair2, filename="analysis/_RNAseq/_plots/descriptive/WGCNA/CutHeight0.40/SignMod_pair2_Mod13_3.jpg", device="jpeg", width=300, height=150, units="mm", dpi="print")

