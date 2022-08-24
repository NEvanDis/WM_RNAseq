# Open R project in top folder

# load packages
library(tidyverse)

# load data ####
d_count <- read.csv("analysis/_RNAseq/_data/TransfExp2019_gene_count_matrix_v3cor.csv")
class(d_count)  
dim(d_count) # v1: 46373 gene_ids, v2: 29095 gene_ids (17278 gene_ids novel transcripts without homologs excluded by BLAST), 
# v3: 30201 gene_ids (filtered out optical duplicates)
# v3cor: 29113 gene_ids (made mistake with grepping novel transcripts, this is correct version)
head(d_count) #first column = gene_id, other 84 cols = each sample, ordered by sample_name

# Gene list ####
genes <- as.data.frame(d_count$gene_id)
names(genes) <- "GeneID"
head(genes)
tail(genes)

# blast results at level of transcript! ####
# NB: BLAST output not included ####

# Swiss prot Drosophila melanogaster ####
sprot_Dros <- read.table("analysis/_RNAseq/_data/uniprot_sprot.Dmelanogaster.blastp.tophits.tsv", sep="\t", quote="")
names(sprot_Dros) <- c("TranscriptID", "Strand", "Read_frame", "Length",
                       "Sequence", "Hit", "Evalue") #qseqid qstrand qframe qlen full_qseq stitle evalue
sprot_Dros <- sprot_Dros %>% mutate(GeneID=ifelse(grepl("OBRU", TranscriptID)==T, gsub("(OBRU01_\\d+)\\-PA_\\d", "\\1", TranscriptID),
                                                  gsub("(MSTRG\\.\\d+)\\.\\d+_\\d","\\1",TranscriptID)), 
                                    TranscriptIDcor=ifelse(grepl("OBRU", TranscriptID)==T, gsub("(OBRU01_\\d+\\-PA)_\\d", "\\1", TranscriptID),
                                                           gsub("(MSTRG\\.\\d+\\.\\d+)_\\d","\\1",TranscriptID)),
                                    Read_frame=ifelse(grepl("OBRU", TranscriptID)==T, gsub("OBRU01_\\d+\\-PA_(\\d)", "\\1", TranscriptID),
                                                      gsub("MSTRG\\.\\d+\\.\\d+_(\\d)","\\1",TranscriptID)),
                                    ProtName=gsub("sp\\|\\w+\\|(\\w+) .+", "\\1", Hit), 
                                    Prot_descr=gsub("sp\\|\\w+\\|\\w+ (.+) OS\\=.+", "\\1", Hit), 
                                    GeneName=ifelse(grepl("GN\\=", Hit)==T,gsub("sp\\|\\w+\\|\\w+ .+GN\\=(.+) PE.+", "\\1", Hit), NA),
                                    Species="Drosophila melanogaster",
                                    Origin="sprot_Dros")
sprot_Dros <- sprot_Dros[!duplicated(sprot_Dros[,c('GeneID', 'GeneName')]),] # keep only unique gene hits
sprot_Dros <- sprot_Dros[!duplicated(sprot_Dros[,c('GeneID', 'Prot_descr')]),] # keep only unique gene hits
head(sprot_Dros[,c(1:4,6,7:13)])
tail(sprot_Dros[,c(1:4,6,7:13)])
length(unique(sprot_Dros$GeneID)) # 10 645 genes of 29 113 map to Drosophila proteins

# Swiss prot all Insects ####
sprot_Insect <- read.table("analysis/_RNAseq/_data/uniprot_sprot.Insects.blastp.tophits.tsv", sep="\t", quote="")
names(sprot_Insect) <- c("TranscriptID", "Strand", "Read_frame", "Length",
                         "Sequence", "Hit", "Evalue") #qseqid qstrand qframe qlen full_qseq stitle evalue
sprot_Insect <- sprot_Insect %>% mutate(GeneID=ifelse(grepl("OBRU", TranscriptID)==T, gsub("(OBRU01_\\d+)\\-PA_\\d", "\\1", TranscriptID),
                                                      gsub("(MSTRG\\.\\d+)\\.\\d+_\\d","\\1",TranscriptID)),
                                        TranscriptIDcor=ifelse(grepl("OBRU", TranscriptID)==T, gsub("(OBRU01_\\d+\\-PA)_\\d", "\\1", TranscriptID),
                                                               gsub("(MSTRG\\.\\d+\\.\\d+)_\\d","\\1",TranscriptID)),
                                        Read_frame=ifelse(grepl("OBRU", TranscriptID)==T, gsub("OBRU01_\\d+\\-PA_(\\d)", "\\1", TranscriptID),
                                                          gsub("MSTRG\\.\\d+\\.\\d+_(\\d)","\\1",TranscriptID)),
                                    ProtName=gsub("sp\\|\\w+\\|(\\w+) .+", "\\1", Hit), 
                                    Prot_descr=gsub("sp\\|\\w+\\|\\w+ (.+) OS\\=.+", "\\1", Hit), 
                                    GeneName=ifelse(grepl("GN\\=", Hit)==T,gsub("sp\\|\\w+\\|\\w+ .+GN\\=(.+) PE.+", "\\1", Hit), NA),
                                    Species=gsub("sp\\|\\w+\\|\\w+ .+OS\\=(.+) OX.+", "\\1", Hit),
                                    Origin="sprot_Insects")
sprot_Insect <- sprot_Insect[!duplicated(sprot_Insect[,c('GeneID', 'GeneName')]),] # keep only unique gene hits
sprot_Insect <- sprot_Insect[!duplicated(sprot_Insect[,c('GeneID', 'Prot_descr')]),] # keep only unique gene hits
head(sprot_Insect[,c(1:4,6,7:13)])
tail(sprot_Insect[,c(1:4,6,7:13)])
length(unique(sprot_Insect$GeneID)) # 11 599 genes of 29 113 map to Insect proteins in Swissprot
table(sprot_Insect$Species)

# TrEMBL prot all Insects ####
tprot_Insect <- read.table("analysis/_RNAseq/_data/uniprot_trembl.Insects.blastp.tophits_wnames.tsv", sep="\t", quote="")
names(tprot_Insect) <- c("TranscriptID", "Strand", "Read_frame", "Length",
                         "Sequence", "Hit", "Evalue") #qseqid qstrand qframe qlen full_qseq stitle evalue
tprot_Insect <- tprot_Insect %>% mutate(GeneID=ifelse(grepl("OBRU", TranscriptID)==T, gsub("(OBRU01_\\d+)\\-PA_\\d", "\\1", TranscriptID),
                                                      gsub("(MSTRG\\.\\d+)\\.\\d+_\\d","\\1",TranscriptID)), 
                                        TranscriptIDcor=ifelse(grepl("OBRU", TranscriptID)==T, gsub("(OBRU01_\\d+\\-PA)_\\d", "\\1", TranscriptID),
                                                               gsub("(MSTRG\\.\\d+\\.\\d+)_\\d","\\1",TranscriptID)),
                                        Read_frame=ifelse(grepl("OBRU", TranscriptID)==T, gsub("OBRU01_\\d+\\-PA_(\\d)", "\\1", TranscriptID),
                                                          gsub("MSTRG\\.\\d+\\.\\d+_(\\d)","\\1",TranscriptID)),
                                        ProtName=gsub("tr\\|\\w+\\|(\\w+) .+", "\\1", Hit), 
                                        Prot_descr=gsub("tr\\|\\w+\\|\\w+ (.+) OS\\=.+", "\\1", Hit), 
                                        GeneName=ifelse(grepl("GN\\=", Hit)==T,gsub("tr\\|\\w+\\|\\w+ .+GN\\=(.+) PE.+", "\\1", Hit), NA),
                                        Species=gsub("tr\\|\\w+\\|\\w+ .+OS\\=(.+) OX.+", "\\1", Hit),
                                        Origin="tprot_Insects")
tprot_Insect <- tprot_Insect[!duplicated(tprot_Insect[,c('GeneID', 'GeneName')]),] # keep only unique gene hits
tprot_Insect <- tprot_Insect[!duplicated(tprot_Insect[,c('GeneID', 'Prot_descr')]),] # keep only unique gene hits
head(tprot_Insect[,c(1:4,6,7:13)])
tail(tprot_Insect[,c(1:4,6,7:13)], n=20)

# Remove identical hits ####
prot <- rbind(sprot_Dros, sprot_Insect, tprot_Insect)
prot <- prot[!duplicated(prot[,c('GeneID', 'GeneName', 'Species')]),] # keep only unique gene hits
prot <- arrange(prot, GeneID, Evalue)
prot$TranscriptID <- prot$TranscriptIDcor
head(prot[,c(1:4,6,7:13)])
tail(prot[,c(1:4,6,7:13)], n=20)
# all individual blast results in separate row

# Alternative format: combine all BLAST results into one table ####
transcripts <- unique(prot[,c("GeneID", "TranscriptID")])

sprot_Dros <- subset(prot, Origin=="sprot_Dros")
names(sprot_Dros)[!colnames(sprot_Dros) %in% c("GeneID", "TranscriptID", "Origin")] <- paste("sprotDros_", names(sprot_Dros)[!colnames(sprot_Dros) %in% c("GeneID", "TranscriptID", "Origin")], sep="")
head(sprot_Dros[,c(1:4,6,7:13)])

sprot_Insect <- subset(prot, Origin=="sprot_Insects")
names(sprot_Insect)[!colnames(sprot_Insect) %in% c("GeneID", "TranscriptID", "Origin")] <- paste("sprotInsect_", names(sprot_Insect)[!colnames(sprot_Insect) %in% c("GeneID", "TranscriptID", "Origin")], sep="")
head(sprot_Insect[,c(1:4,6,7:13)])

tprot_Insect <- subset(prot, Origin=="tprot_Insects")
names(tprot_Insect)[!colnames(tprot_Insect) %in% c("GeneID", "TranscriptID", "Origin")] <- paste("tprotInsect_", names(tprot_Insect)[!colnames(sprot_Insect) %in% c("GeneID", "TranscriptID", "Origin")], sep="")
head(tprot_Insect[,c(1:4,6,7:13)])

prot2 <- merge(transcripts, sprot_Dros[,!colnames(sprot_Dros) %in% c("Origin")], by=c("GeneID", "TranscriptID"), all.x=T)
prot2 <- merge(prot2, sprot_Insect[,!colnames(sprot_Insect) %in% c("Origin")], by=c("GeneID", "TranscriptID"), all.x=T)
prot2 <- merge(prot2, tprot_Insect[,!colnames(tprot_Insect) %in% c("Origin")], by=c("GeneID", "TranscriptID"), all.x=T)
prot2 <- prot2[!duplicated(prot2),] %>% arrange(GeneID, TranscriptID)
head(prot2[,c(1:5,7:15, 17:25, 27:32)])
head(prot2[,c("GeneID", "TranscriptID", "sprotDros_GeneName", "sprotDros_Prot_descr",
             "sprotInsect_GeneName", "sprotInsect_Prot_descr", "tprotInsect_GeneName", "tprotInsect_Prot_descr")], n=50)
###############################################################################################################################

# Get StringTie assigned names ####
# StringTie assigned names might include OBRU transcripts!!!
Str_geneIDs <- read.table("analysis/_RNAseq/_data/stringtie_v2cor_transcripts.tsv", sep="\t", quote="")
Str_geneIDs$GeneID_cor <- gsub("gene_id \"(.+)\"\\; transcript_id.+","\\1",Str_geneIDs$V9)
Str_geneIDs$TranscriptID <- gsub("gene_id .+\\; transcript_id \"(.+)\"\\;","\\1",Str_geneIDs$V9)
Str_geneIDs$TranscriptID <- gsub("(.+)\"\\; ref_gene_id .+","\\1", Str_geneIDs$TranscriptID)
Str_geneIDs <- Str_geneIDs[,c("GeneID_cor", "TranscriptID")]
head(Str_geneIDs)
tail(Str_geneIDs)

str(Str_geneIDs) # there's a space behind each MSTRG transcriptID, no idea why but remove
Str_geneIDs$TranscriptID <- ifelse(grepl("MSTRG",Str_geneIDs$TranscriptID)==T, 
                                   gsub("(MSTRG\\.\\d+\\.\\d+) ", "\\1",Str_geneIDs$TranscriptID)
                                   ,Str_geneIDs$TranscriptID)
str(Str_geneIDs) # gone now
Str_geneIDs <- arrange(Str_geneIDs, TranscriptID)
length(unique(Str_geneIDs$GeneID_cor))

# correct GeneID
prot <- merge(prot, Str_geneIDs, by="TranscriptID") # row lay out, novel transcripts excl in v3cor dropped
prot2 <- merge(prot2, Str_geneIDs, by="TranscriptID") # column lay out, novel transcripts excl in v3cor dropped
###############################################################################################################################

############################
# FUNCTIONAL ANNOTATION ####
############################

# load made annotation tables ####
genesAnnot <- read.table(file="analysis/_RNAseq/_annotation/genesAnnot_wblastp.tsv", header=T)
genesAnnot2 <- read.table(file="analysis/_RNAseq/_annotation/genesAnnot_wblastp_rows.tsv", header=T)
genesAnnot2_ext <- read.table(file="analysis/_RNAseq/_annotation/genesAnnot_wblastp_rows_exp.tsv", header=T)
notAnnot<- read.table(file="analysis/_RNAseq/_annotation/genesNotAnnot_wblastp.txt")

intrp <- read.table(file="analysis/_RNAseq/_annotation/intrp_annot_v2.tsv", header=T)
pannz_sub <- read.table(file="analysis/_RNAseq/_annotation/pannz0.7_annot_v2.tsv", header=T)

#########################################
# Table formatted like Martijn Derks ####
#########################################
genesAnnot <- merge(genes, prot2, by.x="GeneID", by.y="GeneID_cor", all.x=T)
length(unique(genesAnnot$GeneID))
# length(subset(prot2, !(TranscriptID %in% genesAnnot$TranscriptID))$GeneID) # all blast hits in there

notAnnot <- as.data.frame(unique(subset(genesAnnot, is.na(sprotDros_Hit)&is.na(sprotInsect_Hit)&is.na(tprotInsect_Hit))$GeneID))
names(notAnnot) <- "GeneID"
head(notAnnot)
# 4898 genes without annotation -> check what hits found in nr DB for these genes and maybe annotate them with that
# write.table(notAnnot, file="analysis/_RNAseq/_annotation/genesNotAnnot_wblastp.txt",row.names=F, col.names=F)

4898/29113*100 # for 17% didn't find a BLAST hit 

# OBRU genes accidently dropped before in there?
leftout <- subset(Str_geneIDs, TranscriptID %in% c("OBRU01_204069-PA", "OBRU01_205414-PA", "OBRU01_201409-PA"))

subset(genesAnnot, GeneID %in% leftout$GeneID_cor)[,c(1,2,8,18,28)]
subset(notAnnot, GeneID %in% leftout$GeneID_cor) # yes one OBRU gene not BLASTED now! redo nr BLAST whole transcriptome and include there

str(genesAnnot)
genesAnnot <-genesAnnot[,!colnames(genesAnnot) %in% c("GeneID.y", "sprotDros_Strand", "sprotDros_TranscriptIDcor",
                                             "sprotInsect_Strand", "sprotInsect_TranscriptIDcor",
                                             "tprotInsect_Strand", "tprotInsect_TranscriptIDcor")]
# write.table(genesAnnot, file="analysis/_RNAseq/_annotation/genesAnnot_wblastp.tsv", row.names=F)
###############################################################################################################################

###################################################
# table with each BLAST result on separate row ####
###################################################
# genesAnnot2 <- merge(genes, prot, by.x="GeneID", by.y="GeneID_cor", all.x=T)
genesAnnot2 <- genesAnnot2 %>% arrange(GeneID, Evalue)
genesAnnot2 <- genesAnnot2 %>% mutate(UniprotID=gsub("\\w+\\|(.+)\\|.+","\\1",Hit))
str(genesAnnot2)
head(genesAnnot2[,c(1:4, 7:13)])
length(unique(genesAnnot2$GeneID))
max(genesAnnot2$Evalue, na.rm=T) # don't need Evalue cut-off! Highest Evalue below 0.001
genesAnnot2 <- genesAnnot2[,!colnames(genesAnnot2) %in% c("GeneID.y", "Strand", "TranscriptIDcor")]
# write.table(genesAnnot2, file="analysis/_RNAseq/_annotation/genesAnnot_wblastp_rows.tsv", row.names=F)
###############################################################################################################################

#########################################################################
# Try to annotate 4898 without annotation with NCBI nr BLAST results ####
#########################################################################
nr_notAnn <- read.table("analysis/_RNAseq/_data/nr_all.blastp.tophits_sub.tsv", sep="\t", quote="")
# nr BLAST results for not annotated genes, if want to excl. hits "uncharacterized protein" use nr_all.blast.tophits_wnames_sub.tsv
head(nr_notAnn[,c(1:4,6,7)])

names(nr_notAnn) <- c("TranscriptID", "Strand", "Read_frame", "Length",
                       "Sequence", "Hit", "Evalue") #qseqid qstrand qframe qlen full_qseq stitle evalue
nr_notAnn <- nr_notAnn %>% mutate(GeneID=ifelse(grepl("OBRU", TranscriptID)==T, gsub("(OBRU01_\\d+)\\-PA_\\d", "\\1", TranscriptID),
                                                  gsub("(MSTRG\\.\\d+)\\.\\d+_\\d","\\1",TranscriptID)),
                                  TranscriptIDcor=ifelse(grepl("OBRU", TranscriptID)==T, gsub("(OBRU01_\\d+\\-PA)_\\d", "\\1", TranscriptID),
                                                         gsub("(MSTRG\\.\\d+\\.\\d+)_\\d","\\1",TranscriptID)),
                                  Read_frame=ifelse(grepl("OBRU", TranscriptID)==T, gsub("OBRU01_\\d+\\-PA_(\\d)", "\\1", TranscriptID),
                                                    gsub("MSTRG\\.\\d+\\.\\d+_(\\d)","\\1",TranscriptID)),
                                    ProtName=gsub("(.+\\.\\d+) .+", "\\1", Hit), 
                                    Prot_descr=gsub(".+\\.\\d+ (.+) \\[.+\\]", "\\1", Hit), 
                                    GeneName=gsub(".+ (.+) \\[.+\\]","\\1",Hit),
                                    Species=gsub(".+\\[(.+)\\]","\\1",Hit),
                                    Origin="nr",
                                  UniprotID=NA)
nr_notAnn <- nr_notAnn[!duplicated(nr_notAnn[,c('GeneID', 'GeneName')]),] # keep only unique gene hits
nr_notAnn <- nr_notAnn[!duplicated(nr_notAnn[,c('GeneID', 'Prot_descr')]),] # keep only unique gene hits
nr_notAnn$TranscriptID <- nr_notAnn$TranscriptIDcor
head(nr_notAnn[,c(1:4,6,7:13)]) # most of it is hypothetical/uncharacterized protein, but not all
tail(nr_notAnn[,c(1:4,6,7:13)])
str(nr_notAnn)

str(notAnnot)
still_notAnnot <- unique(subset(notAnnot, !(GeneID %in% nr_notAnn$GeneID)))
head(still_notAnnot, n=10) # still 4142 genes not annotated... how can that be? --> I checked some, many Genes that match ref genes/transcripts
# also used Blastp now, for excluding novel sites used Blastx
# write.table(still_notAnnot, file="analysis/_RNAseq/_annotation/genes_stillNotAnnot.txt",row.names=F, col.names=F)

subset(still_notAnnot, GeneID %in% leftout$GeneID_cor) # hmm still not annotated with nr, but I guess wouldn't have found a match with other dbs either then

# correct GeneID
nr_Annot <- merge(nr_notAnn, Str_geneIDs, by="TranscriptID") # row lay out

# Add extra annotation
genesAnnot2_ext <- as.data.frame(subset(genesAnnot2, GeneID %in% nr_Annot$GeneID)$GeneID)
colnames(genesAnnot2_ext) <- "GeneID"
genesAnnot2_ext <- merge(genesAnnot2_ext, nr_Annot[,!colnames(nr_Annot) %in% c("GeneID")], by.x="GeneID", by.y="GeneID_cor", all.x=T)

genesAnnot2b <- subset(genesAnnot2, !(GeneID %in% nr_Annot$GeneID))
genesAnnot2b <- rbind(genesAnnot2b, subset(genesAnnot2_ext, select=-c(Strand,TranscriptIDcor)))
genesAnnot2b <- genesAnnot2b %>% arrange(GeneID, Evalue)
head(genesAnnot2b[,c(1:4, 7:13)], n=10)
length(unique(genesAnnot2b$GeneID))
table(is.na(genesAnnot2b$Evalue))
# write.table(genesAnnot2b, file="analysis/_RNAseq/_annotation/genesAnnot_wblastp_rows_exp.tsv", row.names=F, sep="\t") # 4142/29113 = 14%
###############################################################################################################################

uniprIDs <- unique(genesAnnot2_ext[,c("GeneID", "UniprotID")])
head(uniprIDs)
# write.table(subset(uniprIDs, !is.na(UniprotID)), file="analysis/_RNAseq/_annotation/uniprotIDs.txt", row.names=F, col.names=F)

#########################################################
# Add Interpro results = protein domains and GO terms####
#########################################################
# intrp <- read.table("analysis/_RNAseq/_data/interpro_stringtiev2_sub_cor_nopa.tsv", sep="\t", quote="", fill=T) # pathway info left out
# I don't know why but at least 1 line that's missing a field...

names(intrp) <- c("ProteinID", "SeqMD5", "Length", "Analysis", "Hit", "Prot_descr", "Start", "End",
                       "Evalue", "Match", "RunDate", "IntrpID", "Intrp_descr", "GO")
intrp <- intrp %>% mutate(GeneID=ifelse(grepl("OBRU", ProteinID)==T, gsub("(OBRU01_\\d+)\\-PA_\\d", "\\1", ProteinID),
                                                  gsub("(MSTRG\\.\\d+)\\.\\d+_\\d","\\1",ProteinID)), 
                          TranscriptID= ifelse(grepl("OBRU", ProteinID)==T, gsub("(OBRU01_\\d+\\-PA)_\\d", "\\1", ProteinID),
                                              gsub("(MSTRG\\.\\d+\\.\\d+)_\\d","\\1",ProteinID)), 
                          ORF=ifelse(grepl("OBRU", ProteinID)==T, gsub("OBRU01_\\d+\\-PA_(\\d)", "\\1", ProteinID),
                                     gsub("MSTRG\\.\\d+\\.\\d+_(\\d)","\\1",ProteinID)),
                                    Origin="intrp")
intrp <- intrp[!duplicated(intrp[,c('GeneID', 'IntrpID')]),] # keep only unique gene hits
intrp <- intrp[,c("GeneID", "TranscriptID", "ORF", "Length", "Analysis", "Hit", "IntrpID", "Evalue", "Prot_descr", "Intrp_descr", "GO")]
intrp <- merge(Str_geneIDs, intrp, by="TranscriptID", all.x=T) #correct GeneNames, keep all genes
intrp <- arrange(intrp, GeneID_cor, Evalue) # sort by GeneID and Evalue
head(intrp)
tail(intrp)
str(intrp)
# write.table(intrp, file="analysis/_RNAseq/_annotation/intrp_annot_v2.tsv", row.names=F)

intrp_1hit <- intrp[!duplicated(intrp[,c('GeneID_cor')]),] %>% arrange(GeneID_cor, TranscriptID) # keep only most significant hit
intrp_1hit$GO_cor <- ifelse(is.na(intrp_1hit$GO)==F & intrp_1hit$GO!="-", intrp_1hit$GO, NA)
intrp_1hit$IntrpID_cor <- ifelse(is.na(intrp_1hit$IntrpID)==F & intrp_1hit$IntrpID!="-", intrp_1hit$IntrpID, NA)
table(is.na(intrp_1hit$GO_cor)) #only 6210 Genes out of 29113 have GO terms atm....21%

table(is.na(intrp_1hit$GeneID)) # see for how many genes I have Interpro protein domain results = 16194 genes ~56%
head(intrp_1hit)

###############################################################################################################################


########################################
# GO terms annotation from Pannzer2 ####
########################################
# pannz <- read.table("analysis/_RNAseq/_data/pannzer2_stringtiev2_sub_cor_GO.out", sep="\t", header=T)
head(pannz)

table(pannz$ARGOT_score>=0.7) # example paper Kahilainen et al. only used GOs with scores above 0.7

pannz_sub <- subset(pannz, ARGOT_score>=0.7) %>% mutate(TranscriptID=ifelse(grepl("OBRU", qpid)==T, gsub("(OBRU01_\\d+\\-PA)_\\d", "\\1", qpid),
                                                                            gsub("(MSTRG\\.\\d+\\.\\d+)_\\d","\\1",qpid)), 
                                                        GeneID=ifelse(grepl("OBRU", qpid)==T, gsub("(OBRU01_\\d+)\\-PA_\\d", "\\1", qpid),
                                                                      gsub("(MSTRG\\.\\d+)\\.\\d+_\\d","\\1",qpid)))
pannz_sub <- pannz_sub[!duplicated(pannz_sub[,c('GeneID', 'ontology', 'goid')]),] # keep only unique gene hits
table(pannz_sub$ontology) # 5875 BP GO terms

pannz_sub <- merge(Str_geneIDs, pannz_sub, by="TranscriptID", all.x=T) #correct GeneNames, keep all genes
pannz_sub <- arrange(pannz_sub, GeneID_cor, desc(ARGOT_score)) # sort by GeneID and ARGOT_score
head(pannz_sub)
tail(pannz_sub)
str(pannz_sub)
# write.table(pannz_sub, file="analysis/_RNAseq/_annotation/pannz0.7_annot_v2.tsv", row.names=F)

pannz_1hit <- pannz_sub[!duplicated(pannz_sub[,c('GeneID_cor')]),] %>% arrange(GeneID_cor, TranscriptID) # keep only most significant hit
table(is.na(pannz_1hit$goid)) #only 4603 Genes out of 29113 have GO terms atm.... 16%

# pannzBP_1hit <- subset(pannz_sub, ontology=="BP")
# pannzBP_1hit <- pannzBP_1hit[!duplicated(pannzBP_1hit[,c('GeneID_cor')]),] %>% arrange(GeneID_cor, TranscriptID) # keep only most significant hit
# head(pannzBP_1hit) #only for 4473 genes BP GO term
###############################################################################################################################


##############################
# GO - gene mapping table ####
##############################
# Format needed for topGO: "gene_ID<TAB>GO_ID1, GO_ID2, GO_ID3, ...." 1 gene per line in txt file

# format interpro output ####
head(intrp)

intrp_GO <- intrp[,c("GeneID_cor", "GO")] %>% mutate(GOcor=ifelse(grepl("|", GO)==T, gsub("\\|", ", ", GO), GO)) # GO terms go up to 5, only want to put a comma instead of | 
intrp_GO$GOcor <- ifelse(is.na(intrp_GO$GOcor)==T | intrp_GO$GOcor=="", "-", intrp_GO$GOcor)
intrp_GO <- subset(intrp_GO, GOcor!="-")
intrp_GO <- intrp_GO[!duplicated(intrp_GO[,c("GeneID_cor", "GOcor")]),] # remove duplicates

intrp_GO$key <- NA
ungenes <- unique(intrp_GO$GeneID_cor)

for (gene in 1:length(ungenes)){
  g <- ungenes[gene]
  
  for (hit in 1:length(intrp_GO$key[intrp_GO$GeneID_cor==g])){
  intrp_GO$key[intrp_GO$GeneID_cor==g][hit] <- paste("GO", hit, sep="")
  }
}
head(intrp_GO)
table(is.na(intrp_GO$key))

intrp_GO <- spread(intrp_GO[,c("GeneID_cor", "GOcor", "key")], key="key", value="GOcor")
intrp_GO <- intrp_GO %>% mutate(GO=paste(GO1, GO2, GO3, GO4, GO5, GO6, sep=", "))
intrp_GO$GO <- ifelse(grepl("NA", intrp_GO$GO)==T, gsub(", NA", "", intrp_GO$GO), intrp_GO$GO)
intrp_GO <- intrp_GO[, c("GeneID_cor", "GO")]
head(intrp_GO)

# format pannzer output ####
head(pannz_sub)
pannz_GO <- subset(pannz_sub[,c("GeneID_cor", "goid")], !is.na(goid)) %>% mutate(GO=sprintf("%07d", goid))
pannz_GO$GOcor <- paste("GO:", pannz_GO$GO, sep="") # is this transformation correct?!?! ####
pannz_GO <- pannz_GO[!duplicated(pannz_GO[,c("GeneID_cor", "GOcor")]),] # remove duplicates

pannz_GO$key <- NA
ungenes <- unique(pannz_GO$GeneID_cor)

for (gene in 1:length(ungenes)){
  g <- ungenes[gene]
  
  for (hit in 1:length(pannz_GO$key[pannz_GO$GeneID_cor==g])){
    pannz_GO$key[pannz_GO$GeneID_cor==g][hit] <- paste("GO", hit, sep="")
  }
}
head(pannz_GO)
table(is.na(pannz_GO$key))

pannz_GO <- spread(pannz_GO[,c("GeneID_cor", "GOcor", "key")], key="key", value="GOcor")
subset(pannz_GO,!is.na(GO50)) # for one gene 50 GO terms

cols <- paste(colnames(pannz_GO)[2:51],collapse=', ')

pannz_GO <- pannz_GO %>% mutate(GO=paste(GO1, GO10, GO11, GO12, GO13, GO14, GO15, GO16, GO17, GO18, GO19, GO2, GO20, GO21, GO22, GO23, GO24, GO25, GO26, GO27, GO28, GO29, GO3, GO30, GO31, GO32, GO33, GO34, GO35, GO36, GO37, GO38, GO39, GO4, GO40, GO41, GO42, GO43, GO44, GO45, GO46, GO47, GO48, GO49, GO5, GO50, GO6, GO7, GO8, GO9, sep=", "))
pannz_GO$GO <- ifelse(grepl("NA", pannz_GO$GO)==T, gsub(", NA", "", pannz_GO$GO), pannz_GO$GO)
pannz_GO <- pannz_GO[, c("GeneID_cor", "GO")]
head(pannz_GO)

# UniprotID to GO mapping ####
head(uniprIDs)

uniprGOmap <- read.table("analysis/_RNAseq/_data/uniprot.blastp.gos.tsv", fill=T)
colnames(uniprGOmap)[c(2,5)] <- c("UniprotID", "GO")
head(uniprGOmap)

unipr_GO <- merge(uniprIDs, uniprGOmap[,c("UniprotID", "GO")], by="UniprotID", all.x=T)
unipr_GO <- unique(subset(unipr_GO, !is.na(UniprotID)))
table(is.na(unipr_GO$GO)) # for 12438/54019 UniprotIDs no GO term ~23%
unipr_GO <- subset(unipr_GO, !is.na(GO))
head(unipr_GO)

unipr_GO$key <- NA
ungenes <- unique(unipr_GO$GeneID)

for (gene in 1:length(ungenes)){ # takes really long time to run
  g <- ungenes[gene]
  
  for (hit in 1:length(unipr_GO$key[unipr_GO$GeneID==g])){
    unipr_GO$key[unipr_GO$GeneID==g][hit] <- paste("GO", hit, sep="")
  }
}
head(unipr_GO)
table(is.na(unipr_GO$key))

unipr_GO <- spread(unipr_GO[,c("GeneID", "GO", "key")], key="key", value="GO")
subset(unipr_GO,!is.na(GO195)) # for one gene 195 GO terms

cols <- paste(colnames(unipr_GO)[2:196],collapse=', ')

unipr_GO <- unipr_GO %>% mutate(GO=paste(GO1, GO10, GO100, GO101, GO102, GO103, GO104, GO105, GO106, GO107, GO108, GO109, GO11, GO110, GO111, GO112, GO113, GO114, GO115, GO116, GO117, GO118, GO119, GO12, GO120, GO121, GO122, GO123, GO124, GO125, GO126, GO127, GO128, GO129, GO13, GO130, GO131, GO132, GO133, GO134, GO135, GO136, GO137, GO138, GO139, GO14, GO140, GO141, GO142, GO143, GO144, GO145, GO146, GO147, GO148, GO149, GO15, GO150, GO151, GO152, GO153, GO154, GO155, GO156, GO157, GO158, GO159, GO16, GO160, GO161, GO162, GO163, GO164, GO165, GO166, GO167, GO168, GO169, GO17, GO170, GO171, GO172, GO173, GO174, GO175, GO176, GO177, GO178, GO179, GO18, GO180, GO181, GO182, GO183, GO184, GO185, GO186, GO187, GO188, GO189, GO19, GO190, GO191, GO192, GO193, GO194, GO195, GO2, GO20, GO21, GO22, GO23, GO24, GO25, GO26, GO27, GO28, GO29, GO3, GO30, GO31, GO32, GO33, GO34, GO35, GO36, GO37, GO38, GO39, GO4, GO40, GO41, GO42, GO43, GO44, GO45, GO46, GO47, GO48, GO49, GO5, GO50, GO51, GO52, GO53, GO54, GO55, GO56, GO57, GO58, GO59, GO6, GO60, GO61, GO62, GO63, GO64, GO65, GO66, GO67, GO68, GO69, GO7, GO70, GO71, GO72, GO73, GO74, GO75, GO76, GO77, GO78, GO79, GO8, GO80, GO81, GO82, GO83, GO84, GO85, GO86, GO87, GO88, GO89, GO9, GO90, GO91, GO92, GO93, GO94, GO95, GO96, GO97, GO98, GO99, sep=", "))
unipr_GO$GO <- ifelse(grepl("NA", unipr_GO$GO)==T, gsub(", NA", "", unipr_GO$GO), unipr_GO$GO)
unipr_GO <- unipr_GO[, c("GeneID", "GO")]
head(unipr_GO)
length(unique(unipr_GO$GeneID)) # 18562/29113 annotated with at least one GO term = 64%
# write.table(unipr_GO, file="analysis/_RNAseq/_annotation/uniprot_GOannot.tsv", row.names=F)

# combine the three ####
GOtbl <- merge(intrp_GO, pannz_GO, by="GeneID_cor", all=T)
GOtbl <- merge(GOtbl, unipr_GO, by.x="GeneID_cor", by.y="GeneID", all=T)
GOtbl$GO <- paste(GOtbl$GO.x, GOtbl$GO.y, GOtbl$GO, sep=", ")
GOtbl$GO <- ifelse(grepl("NA", GOtbl$GO)==T, gsub(", NA", "", GOtbl$GO), GOtbl$GO)
GOtbl$GO <- ifelse(grepl("NA", GOtbl$GO)==T, gsub("NA, ", "", GOtbl$GO), GOtbl$GO)
table(grepl("NA", GOtbl$GO))
GOtbl <- GOtbl[, c("GeneID_cor", "GO")]

# check if strings do not contain doubles ####
head(GOtbl)

GOtbl1 <- GOtbl %>% separate(GO, sep = ",", into=c(sprintf("GO%03d", seq(1:198))))
GOtbl1 <- GOtbl1 %>% gather(key="key", value="GOterm", colnames(GOtbl1)[c(2:199)]) %>% 
  arrange(GeneID_cor, key)
GOtbl1 <- GOtbl1[!duplicated(GOtbl1[,c("GeneID_cor", "GOterm")]),] # lots of duplicates, but prob mostly NAs
head(GOtbl1)
table(is.na(GOtbl1$GOterm))

GOtbl1 <- subset(GOtbl1, !is.na(GOterm)) %>% spread(key="key", value="GOterm")
subset(GOtbl1,!is.na(GO189)) # for one gene 189 GO terms

cols <- paste(colnames(GOtbl1)[2:190],collapse=', ')

GOtbl1 <- GOtbl1 %>% mutate(GO=paste(GO001, GO002, GO003, GO004, GO005, GO006, GO007, GO008, GO009, GO010, GO011, GO012, GO013, GO014, GO015, GO016, GO017, GO018, GO019, GO020, GO021, GO022, GO023, GO024, GO025, GO026, GO027, GO028, GO029, GO030, GO031, GO032, GO033, GO034, GO035, GO036, GO037, GO038, GO039, GO040, GO041, GO042, GO043, GO044, GO045, GO046, GO047, GO048, GO049, GO050, GO051, GO052, GO053, GO054, GO055, GO056, GO057, GO058, GO059, GO060, GO061, GO062, GO063, GO064, GO065, GO066, GO067, GO068, GO069, GO070, GO071, GO072, GO073, GO074, GO075, GO076, GO077, GO078, GO079, GO080, GO081, GO082, GO083, GO084, GO085, GO086, GO087, GO088, GO089, GO090, GO091, GO092, GO093, GO094, GO095, GO096, GO097, GO098, GO099, GO100, GO101, GO102, GO103, GO104, GO105, GO106, GO107, GO108, GO109, GO110, GO111, GO112, GO113, GO114, GO115, GO116, GO117, GO118, GO119, GO120, GO121, GO122, GO123, GO124, GO125, GO126, GO127, GO128, GO129, GO130, GO131, GO132, GO133, GO134, GO135, GO136, GO137, GO138, GO139, GO140, GO141, GO142, GO143, GO144, GO145, GO146, GO147, GO148, GO149, GO150, GO151, GO152, GO153, GO154, GO155, GO156, GO157, GO158, GO159, GO160, GO161, GO162, GO163, GO164, GO165, GO166, GO167, GO168, GO169, GO170, GO171, GO172, GO173, GO174, GO175, GO176, GO178, GO180, GO181, GO183, GO184, GO187, GO189, GO192, GO193, GO194, GO195, GO197, GO198, sep=", "))
GOtbl1$GO <- ifelse(grepl("NA", GOtbl1$GO)==T, gsub(", NA", "", GOtbl1$GO), GOtbl1$GO)
table(grepl("NA", GOtbl1$GO))
GOtbl <- GOtbl1[, c("GeneID_cor", "GO")]
colnames(GOtbl) <- c("GeneID", "GO")
rm(GOtbl1)
head(GOtbl)

# make sure all genes are in final GO table ####
head(genes)

GOtbl <- merge(genes, GOtbl, by="GeneID", all.x=T)
table(is.na(GOtbl$GO)) # before: 12895/29113 = 44%, v2 pannzer + interpo 10828/29113 = 37%; v2 + uniprot mapping 19153/29113= 66%
# write.table(GOtbl, file="analysis/_RNAseq/_annotation/TransfExp_WM_gene2GO_v2.txt", row.names=F, sep="\t", col.names=F)
# need to remove \" and double space between GO terms in Notepad++ !!

###############################################################################################################################


####################################
# Stats novel transcripts/genes ####
####################################
gff <- read.table(file="analysis/_RNAseq/_data/gff.annotations.txt", sep=",", header=F, fill=T)
head(gff)
str(gff)

# Get relevant info
gff <- gff %>% mutate(TranscriptID=gsub("transcript_id (.+); gene_id .+","\\1",V1), GeneID=gsub("transcript_id .+; gene_id (.+); gene_name.+","\\1",V1), 
                      GeneName=ifelse(grepl("gene_name", V1)==T, gsub(".+; gene_id .+; gene_name (.+); xloc .+","\\1",V1), NA), 
                      Class=gsub(".+\\; class_code (.+)\\; .+","\\1",V1)) %>%
  mutate(GeneID=ifelse(grepl("xloc", GeneID)==T, gsub("transcript_id .+; gene_id (.+); xloc .+","\\1",V1), GeneID)) %>%
  select(TranscriptID, GeneID, GeneName, Class)
head(gff)
table(gff$Class, useNA="always") 
unique(gff$GeneID) %>% length # 29,113 genes
unique(gff$GeneName) %>% length # part of novel transcripts part of ref genes, others completely new
unique(gff$TranscriptID) %>% length # all unique transcripts, or some double? All unique it seems

# Number of transcripts per class per gene
ClassCounts <- as.data.frame(table(gff$GeneID, gff$Class))
colnames(ClassCounts) <- c("GeneID", "Class", "Freq")
head(ClassCounts)

ClassCounts2 <- ClassCounts %>% spread(key="Class", value="Freq")
head(ClassCounts2) # '=' complete match with ref; 'j' potential novel isoform; 'u' novel intergenic transcript
ClassCounts2$NrTranscripts <- rowSums(ClassCounts2[,c(2:13)]) # Number of transcripts per gene
summary(ClassCounts2$NrTranscripts )
head(ClassCounts2)
# write.table(ClassCounts2, file="analysis/_RNAseq/_results/gff_stats.tsv", sep="\t", row.names=F)

ClassCounts3 <- ClassCounts2
ClassCounts3[,c(2:13)] <- ClassCounts3[,c(2:13)]!=0
head(ClassCounts3)
colSums(ClassCounts3[,c(2:13)]) # number of TRUES per column = number of genes for each category

table(ClassCounts2$`=`!=0) #check

# How many novel genes?
novelgenes <- filter(ClassCounts2, u!=0)
colSums(novelgenes[,c(2:10,12,13)]) # col 11 = class u, not all genes with u transcripts are novel genes
novelgenes <- filter(novelgenes, `=`==0, i==0, j==0, k==0, m==0, n==0, o==0, p==0, s==0, x==0, y==0)
head(novelgenes)
colSums(novelgenes[,c(2:13)])

## Stats for DEGs
DEGs <- read.table(file="analysis/_RNAseq/_results/DEGs_annot_n2_p0.01.tsv", header=T) # produced by later script 5a ####
head(DEGs)
unique(DEGs$GeneID) %>% length

DEG_class <- filter(ClassCounts3, GeneID %in% unique(DEGs$GeneID))
colSums(DEG_class[,c(2:13)]) # number of TRUES per column = number of genes for each category

DEGs_novel <- filter(novelgenes, GeneID %in% unique(DEGs$GeneID))
head(DEGs_novel)

DEGs_novel_annot <- filter(DEGs, GeneID %in% DEGs_novel$GeneID)
filter(DEGs_novel_annot, Species=="Drosophila melanogaster")
filter(DEGs, GeneID=="MSTRG.20613") # interesting novel gene as example, Tao
filter(DEGs, GeneID=="MSTRG.1055") # interesting novel gene as example, takeout
filter(DEGs, GeneID=="MSTRG.9103") # interesting novel gene as example, wingless
