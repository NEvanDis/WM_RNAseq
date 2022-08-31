#Open R project file in top folder
#load required packages
library(tidyverse)
library(readxl)
library(cowplot)
theme_set(theme_cowplot())#white background instead of grey -> don't load if want grey grid
pd <- position_dodge(width = 1) #prevent overlap

# Data
d_2019 <- read.csv("./analysis/_TransfExp2019/fitted_values/d_im2019_tidy.csv") #from 2019 analysis script, data not included ####
head(d_2019)
table(d_2019$treatment) # fix are RNA samples

RNA_sampl <- read.csv(file="./analysis/RNA_available-clutches2019.csv")
head(RNA_sampl)

final_sel <- c("136", "326", "473", # w2: 136, 326, 473
               "102", "219", "367", # w4: 102, 219, 367
               "128", "390", "471", # w6: 128, 390, 471
               "94", "407", "411")  # w8: 94, 407, 411

RNA_im <- filter(d_2019, Tube %in% final_sel)
unique(RNA_im$Tube) %>% length

RNA_sel <- filter(RNA_sampl, Tube %in% final_sel) %>% mutate(Tube=as.factor(Tube))
RNA_sel$Tube <- factor(RNA_sel$Tube, levels(devcatnum$Tube)[c(4,6,12, #w2
                                                                  2,5,7, #4
                                                                  3,8,11, #6
                                                                  1,9,10)]) #8 # order to pick order of plot colors

head(RNA_sel)
str(RNA_sel)


# Get clutch size range used for experiment
clsize2019 <- read_xlsx("data/_TransfExp2019/20191115_expsamples.xlsx", sheet=7) # data not included ####
head(clsize2019)
table(clsize2019$treat_groupshort)
table(clsize2019$sample)

str(clsize2019)
clsize <- aggregate(as.numeric(`#eggs_sub`)~as.factor(`Tube nbr`), filter(clsize2019 , sample!="shed"&sample!="shed2"), sum)
names(clsize) <- c("Tube", "Nr_eggs")
range(clsize$Nr_eggs) # without eggs that went into the shed

clsize2019b <- clsize2019[!duplicated(clsize2019[,"Tube nbr"]),]
range(clsize2019b$Eggs) # with eggs that went into the shed


#plot dev stage available samples per week ####
devcatnum <- aggregate(imaged_on~sample_week + dst_num + Tube, data=RNA_im, length) %>% # number of embryos per stage per Treat_week
  arrange(sample_week, Tube) %>% mutate(Tube=as.factor(Tube))
devcatnum$Tube <- factor(devcatnum$Tube, levels(devcatnum$Tube)[c(4,6,12, #w2
                                                                  2,5,7, #4
                                                                  3,8,11, #6
                                                                  1,9,10)]) #8 # order to pick order of plot colors
head(devcatnum)
range(devcatnum$dst_num)
table(devcatnum$imaged_on)

col <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdGy"))(55)[c(1,5,10,15,20,25,30,37,45,40,47,55)] #choose colors that are far apart from palet
col <- rev(col)

col2 <- c("skyblue4","steelblue1","gray10","royalblue1", "darkslateblue", "grey","sienna3","firebrick4","orangered2","tomato2","chocolate4","red3")

# Graph with clear discrete values
p_dst50 <- ggplot(devcatnum, aes(y=dst_num, x=sample_week, fill=Tube, group=Tube))+ #
  scale_size_continuous(range=c(3,18))+
  scale_fill_manual(values = c(col2))+
  geom_point(aes(size=imaged_on), position=pd, shape=21, col="black", stroke=1.5, alpha=0.8)+
  geom_point(data=RNA_sel, aes(y=dst50_before.50., x=Treat_week), size=8, position=pd, shape=4, col="grey8", stroke=3)+ # medians
  labs(y="Development stage", size='N observed')+
  scale_y_continuous(breaks=seq(6,12,by=1)) + scale_x_continuous(breaks=seq(2,8,by=2), labels=c("Week2", "Week4", "Week6", "Week8")) +
  theme(axis.title.x = element_blank(), axis.text.x  = element_blank(), axis.ticks.length=unit(0.2,"inch"))+ #changes y axis label size
  theme(axis.title.y = element_blank(), axis.text.y=element_text(angle=0, size=35), legend.title=element_text(size=20), legend.text=element_text(size=20))+ #changes x axis label size
  #theme(legend.position="none")+
  coord_flip(ylim=c(5.5,12.5))+
  guides(fill = "none")+
  background_grid(major = "xy") #get grid lines
p_dst50
ggsave(filename="./analysis/RNAsamp_dst50_sel2a.png", plot=p_dst50, device="png", width=550, height=150, units="mm", dpi="print") # save with legend
ggsave(filename="./analysis/RNAsamp_dst50_sel2b.png", plot=p_dst50, device="png", width=550, height=150, units="mm", dpi="print") # save without legend

