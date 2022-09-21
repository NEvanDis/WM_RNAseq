### Working directory: "/home/NIOO.INT/natalied/WM_RNAseq"

### Conda environment location: "/home/NIOO.INT/natalied/.conda/envs/RNAseq"
```
conda activate RNAseq 		# to activate the environment containing software I need
conda list 			# check what software is installed in environment
conda install <package name> 	# install extra software while in environment
conda deactivate 		# to quite the environment again
conda env export > src/env_WM_RNAseq_pipeline.yml	# save software versions for reproducibility
``` 

### more info on conda: 
 - https://github.com/bioconda/bioconda-recipes
 - https://bioconda.github.io/index.html #

## Get data
``` 
ln -s <path name> 		#make symbolic link to data location, so that can access it easily from my analysis folder
```

# FUNCTIONAL ANNOTATION: Step 1
Do diamond blast search on scratch/TransfExp2019/stringtie_merged_v2_sub_cor.gtf

```
## Get .fasta file from Stringtie
## conda_env="gffconversion"
gffread -w scratch/TransfExp2019/stringtie_merged_v2_sub_cor.fa -g scratch/ref/Obru1.fsa scratch/TransfExp2019/stringtie_merged_v2_sub_cor.gtf
## translate .gtf file Stringtie used for count table to .fasta file

## Could use blastp instead of blastx, which is faster, but then need to translate to proteins myself, tool: http://emboss.sourceforge.net/apps/cvs/emboss/apps/transeq.html
## Right frame should be = 1, because data come from transcriptome
## BUT only small part of the assembled transcripts actually starts with start codon: 
cat scratch/TransfExp2019/stringtie_merged_v2_sub_cor.fa | fasta_formatter  | paste - - | cut -f 2 | grep -P '^ATG' | wc -l
## And they are not divisible by 3:
cat scratch/TransfExp2019/stringtie_merged_v2_sub_cor.fa | fasta_formatter  | paste - - | cut -f 2 | awk '{print length($0)%3}' |  grep -v 0 | less
## So use reading frames 1 to 3 (assuming strand info from Stringtie is ok) 

##---------------------------------------
## UPDATE 23/9/2021 
##---------------------------------------
## Change to blastp instead of blastx

## Translate nucleotide to protein seq
transeq -sequence scratch/TransfExp2019/stringtie_merged_v2_sub_cor.fa -outseq scratch/TransfExp2019/stringtie_merged_v2_sub_cor.pep -frame=F -clean
## Translate reading frame 1-3, because at transcript level (should already be in correct order i.e. forward/reverse) 
## -clean to have an 'X' as STOP codon instead of '*' to serve as Interproscan input

## Do diamond blast search, outfmt (-f) as tabular output
## conda_env="RNAseq"
diamond makedb --in uniprot_trembl.fasta.gz -d uniprot_trembl-diamond-2.0.6.dmnd
## make diamond blast database

## 1. Swissprot ##SKIP
#diamond blastx -q scratch/TransfExp2019/stringtie_merged_v2_sub.fa -d /scratch/blast/uniprot/20210210/uniprot_sprot-diamond-2.0.6.dmnd -p 8 -f 6 qseqid qstrand qframe qlen full_qseq stitle evalue -o results/blast/uniprot_sprot.blast.tsv

## 2. TrEMBL ##SKIP
#diamond blastx -q scratch/TransfExp2019/stringtie_merged_v2_sub.fa -d /scratch/blast/uniprot/20210210/uniprot_trembl-diamond-2.0.6.dmnd -b12 -c1 -p 8 -f 6 qseqid qstrand qframe qlen full_qseq stitle evalue -o results/blast/uniprot_trembl.blast.tsv

## 3. Swissprot only Drosophila
diamond blastp -q scratch/TransfExp2019/stringtie_merged_v2_sub_cor.pep -d /scratch/blast/uniprot/20210210/uniprot_sprot.Dmelanogaster-diamond-2.0.6.dmnd -p 8 -f 6 qseqid qstrand qframe qlen full_qseq stitle evalue -o results/blast/uniprot_sprot.Dmelanogaster.blastp.tsv

## 4. Swissprot only Insects
diamond blastp -q scratch/TransfExp2019/stringtie_merged_v2_sub_cor.pep -d /scratch/blast/uniprot/20210210/uniprot_sprot.Insects-diamond-2.0.6.dmnd -p 8 -f 6 qseqid qstrand qframe qlen full_qseq stitle evalue -o results/blast/uniprot_sprot.Insects.blastp.tsv

## 5. TrEMBL only Insects
diamond blastp -q scratch/TransfExp2019/stringtie_merged_v2_sub_cor.pep -d /scratch/blast/uniprot/20210210/uniprot_trembl.Insects-diamond-2.0.6.dmnd -b12 -c1 -p 8 -f 6 qseqid qstrand qframe qlen full_qseq stitle evalue -o results/blast/uniprot_trembl.Insects.blastp.tsv

```

Next, make table with best hit, blast hits are already sorted by transcript and e-value (ascending)
```
#sort -u -k1,1 results/blast/uniprot_sprot.blast.tsv > results/blast/uniprot_sprot.blast.tophits.tsv
#sort -u -k1,1 results/blast/uniprot_trembl.blast.tsv > results/blast/uniprot_trembl.blast.tophits.tsv
sort -u -k1,1 results/blast/uniprot_sprot.Dmelanogaster.blastp.tsv > results/blast/uniprot_sprot.Dmelanogaster.blastp.tophits.tsv
sort -u -k1,1 results/blast/uniprot_sprot.Insects.blastp.tsv > results/blast/uniprot_sprot.Insects.blastp.tophits.tsv
sort -u -k1,1 results/blast/uniprot_trembl.Insects.blastp.tsv > results/blast/uniprot_trembl.Insects.blastp.tophits.tsv
## keep only first line for each transcript

head -n 50 results/blast/uniprot_sprot.Dmelanogaster.blastp.tophits.tsv | cut -f 1-4,6,7
## have look at results, cut out sequence column

#grep -v 'ncharacterized protein' results/blast/uniprot_trembl.blast.tsv | sort -u -k1,1 > results/blast/uniprot_trembl.blast.tophits_wnames.tsv
grep -v 'ncharacterized protein' results/blast/uniprot_trembl.Insects.blastp.tsv | sort -u -k1,1 > results/blast/uniprot_trembl.Insects.blastp.tophits_wnames.tsv
## exclude uncharacterized proteins
```

Need to match StringTie GeneIDs to REF to use for making annotation table in R
```
# grep -P 'transcript\t' scratch/TransfExp2019/stringtie_merged_v2_sub.gtf > scratch/TransfExp2019/stringtie_v2_transcripts.tsv ## OLD
grep -P 'transcript\t' scratch/TransfExp2019/stringtie_merged_v2_sub_cor.gtf > scratch/TransfExp2019/stringtie_v2cor_transcripts.tsv
```

See if I can annotate genes without hits in Swissprot/TrEMBL, with results from nr blast
I just rerun nr blast and have it produce tabular output, because grep with flag -f doesn't work well in XML files

```
## 6. NCBI nr ALL
diamond blastp -q scratch/TransfExp2019/stringtie_merged_v2_sub_cor.pep -d /data/shared/db/blast/nr/20200128/nr-diamond-0.9.29.dmnd -b12 -c1 -p 8 -f 6 qseqid qstrand qframe qlen full_qseq stitle evalue -o results/blast/nr_all.blastp.tsv

# grep -v 'ncharacterized protein' results/blast/nr_all.blast.tsv | sort -u -k1,1 > results/blast/nr_all.blast.tophits_wnames.tsv
sort -u -k1,1 results/blast/nr_all.blastp.tsv > results/blast/nr_all.blastp.tophits.tsv
grep -f results/blast/genesNotAnnot_wblastp_regexp.txt results/blast/nr_all.blastp.tophits.tsv > results/blast/nr_all.blastp.tophits_sub.tsv
```

---- OLD ------------------------------------------------------------
Check what's with the still not annotated genes
```
grep -f results/blast/genes_stillNotAnnot_regexp.txt results/blast/nr_all.blast.tophits.tsv # not in here
grep -f results/blast/genes_stillNotAnnot_regexp.txt results/blast/nr_all.blast.tsv #not in here

grep -f results/blast/genes_stillNotAnnot_regexp.txt results/blast/novel_sites.blast_v2.xml  # why are they in here if I used the exact same db for blast like before?! Something different about translation??

## for the novel_sites.blast I used:
# grep -P 'StringTie\ttranscript' scratch/TransfExp2019/stringtie_merged_v2.gtf | grep -v '_id "OBRU01_' > results/novel_sites.stringtie_v2.gtf
# gffread -w results/novel_sites.stringtie_v2.fa -g scratch/ref/Obru1.fsa results/novel_sites.stringtie_v2.gtf
# diamond blastx -q results/novel_sites.stringtie_v2.fa -d /data/db/blast/nr/20200128/nr-diamond-0.9.29.dmnd -p 8 -f 5 -o results/blast/novel_sites.blast_v2.xml

## for entire genome I used:
# gffread -w scratch/TransfExp2019/stringtie_merged_v2_sub.fa -g scratch/ref/Obru1.fsa scratch/TransfExp2019/stringtie_merged_v2_sub.gtf
# diamond blastx -q scratch/TransfExp2019/stringtie_merged_v2_sub.fa -d /data/shared/db/blast/nr/20200128/nr-diamond-0.9.29.dmnd -b12 -c1 -p 8 -f 6 qseqid qstrand qframe qlen full_qseq stitle evalue -o results/blast/nr_all.blast.tsv

## difference = I also kept exon info in when translating to .fa file + other output format
## Does something weird!!!

sed -n -e '/MSTRG\.10\.2$/,/>MSTRG/ p' scratch/TransfExp2019/stringtie_merged_v2_sub.fa | wc
sed -n -e '/MSTRG\.10\.2$/,/>MSTRG/ p' results/novel_sites.stringtie_v2.fa | wc
## same transcript translated has different length in the two files

grep -E '"MSTRG\.10\.2"' results/novel_sites.stringtie_v2.gtf
grep -E '"MSTRG\.10\.2"' scratch/TransfExp2019/stringtie_merged_v2_sub.gtf
## while the transcript info line is exactly the same, only stringtie_merged_v2_sub.gtf retained exon lines info

## redone conversion novel transcripts to .fa
sed -n -e '/MSTRG\.10\.2$/,/>MSTRG/ p' results/novel_sites.stringtie_v2_cor.fa | wc
# FIXED, same result as scratch/TransfExp2019/stringtie_merged_v2_sub.fa now
```
------------------------------------------------------------


# FUNCTIONAL ANNOTATION: Step 2
Get GO terms and Protein domains with Interproscan
See: https://interproscan-docs.readthedocs.io/en/latest/
Run with conda=RNAseq activated
```
## Use small part to test Interproscan
#head scratch/TransfExp2019/stringtie_merged_v2_sub.pep -n 80000 > scratch/TransfExp2019/TEST_stringtie_merged_v2_sub.pep
## 3532 peptides

## TEST RUN INTERPRO SCAN
#./src/my_interproscan/interproscan-5.51-85.0/interproscan.sh -cpu 8 -appl Pfam,SUPERFAMILY -goterms -i scratch/TransfExp2019/TEST_stringtie_merged_v2_sub.pep -b results/interpro/interproTEST_stringtiev2_sub
## WORKS FINE

## Have memory problems if I run the entire file at once = 253,524 proteins (507,048 before)
## I changed number of embedded workers in interproscan.properties to 4 (instead of 6, max 8)
## Suggested to chunk up file in 80,000 proteins per file to improve performance =~ 3 files
head scratch/TransfExp2019/stringtie_merged_v2_sub_cor.pep -n 1800014 > scratch/TransfExp2019/stringtie_merged_v2_sub_cor_1.pep # 85465 proteins
head scratch/TransfExp2019/stringtie_merged_v2_sub_cor_1.pep -n 900015 > scratch/TransfExp2019/stringtie_merged_v2_sub_cor_1a.pep
tail scratch/TransfExp2019/stringtie_merged_v2_sub_cor_1.pep -n +1800014 > scratch/TransfExp2019/stringtie_merged_v2_sub_cor_1b.pep

tail scratch/TransfExp2019/stringtie_merged_v2_sub_cor.pep -n +1800015 | head -n 1700002 > scratch/TransfExp2019/stringtie_merged_v2_sub_cor_2.pep
# 84297 proteins

tail scratch/TransfExp2019/stringtie_merged_v2_sub_cor.pep -n +3500017 > scratch/TransfExp2019/stringtie_merged_v2_sub_cor_3.pep
# 83762 proteins

grep ">" scratch/TransfExp2019/stringtie_merged_v2_sub_cor.pep | wc
grep ">" scratch/TransfExp2019/stringtie_merged_v2_sub_cor_*.pep | wc
253,524 proteins

## RUN INTERPRO SCAN	
## Analyses: Pfam-33.1, SUPERFAMILY-1.75
(./src/my_interproscan/interproscan-5.52-86.0/interproscan.sh -cpu 8 -appl Pfam,SUPERFAMILY -goterms -pa -i scratch/TransfExp2019/stringtie_merged_v2_sub_cor_1.pep -b results/interpro/interpro_stringtiev2_sub_cor_1) 2> scratch/logs/interproscan1.log
(./src/my_interproscan/interproscan-5.52-86.0/interproscan.sh -cpu 8 -appl Pfam,SUPERFAMILY -goterms -pa -i scratch/TransfExp2019/stringtie_merged_v2_sub_cor_2.pep -b results/interpro/interpro_stringtiev2_sub_cor_2) 2> scratch/logs/interproscan2.log
(./src/my_interproscan/interproscan-5.52-86.0/interproscan.sh -cpu 8 -appl Pfam,SUPERFAMILY -goterms -pa -i scratch/TransfExp2019/stringtie_merged_v2_sub_cor_3.pep -b results/interpro/interpro_stringtiev2_sub_cor_3) 2> scratch/logs/interproscan3.log
## each run ~3,5h 

## where:
##     -cpu = number of cores to use, default=8
##     -t = seqtype if other than protein (n=nucleic acid sequences)
##     -appl = list of analyses I want to run
##     -goterms = also provide Gene Ontology (GO) mappings 
##     -i = input fasta file
##     -b = outpuf-file-base file name
##     -pa = provides mappings from matches to pathway information KEGG/MetaCyc/Reactome

## Append results into one file
cat results/interpro/*.tsv >> results/interpro/interpro_stringtiev2_sub_cor.tsv

## Drop pathway column for now, KEGG never in there
cut -f 1-14 results/interpro/interpro_stringtiev2_sub_cor.tsv > results/interpro/interpro_stringtiev2_sub_cor_nopa.tsv

```

Get GO terms with Pannzer2 web service, download results into results/pannzer2/
```
## v2_sub_cor_1a.pep
wget http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp//Nft1HEWoGle/GO.out -O v2_sub_cor_pep1a_GO.out
wget http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp//Nft1HEWoGle/anno.out -O v2_sub_cor_pep1a_anno.out
wget http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp//Nft1HEWoGle/DE.out -O v2_sub_cor_pep1a_DE.out

## v2_sub_cor_1b.pep
wget http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp//iSCvxTBznwr/GO.out -O v2_sub_cor_pep1b_GO.out

## v2_sub_cor_2.pep
wget http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp//XvevLaxXqC7/GO.out -O v2_sub_cor_pep2_GO.out

## v2_sub_cor_3.pep
wget http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp//iJ7kjCWZlGz/GO.out -O v2_sub_cor_pep3_GO.out


## Paste all GO.out files together into one
head v2_sub_cor_pep1a_GO.out -n 1 > pannzer2_stringtiev2_sub_cor_GO.out
for RESULT in *GO.out ; do tail $RESULT -n +2 >> pannzer2_stringtiev2_sub_cor_GO.out; done
```

Get GO terms and KEGG mapping with UniprotIDs from Blast results
```
zcat src/goa_uniprot_all.gaf.gz | grep -f results/blast/uniprotIDs_clean.txt > results/blast/uniprot.blastp.gos.tsv

wget https://www.uniprot.org/mapping/M20210927A084FC58F6BBA219896F365D15F2EB4420796AM.tab -O uniprot.blastp.KEGG.tsv
```

# FUNCTIONAL ANNOTATION: Step 3
Add functional annotation to Stringtie produced transcriptome according to convention (author: Judith E. Risse).
```
# convert gtf to gff and create mRNA and gene features
singularity exec agat_0.8.0--pl5262hdfd78af_0.sif agat_convert_sp_gxf2gxf.pl -g stringtie_merged_v2_sub_cor.gtf > stringtie_merged_v2_sub_cor2.gff

# remove agat log lines from gff output
# annotate with interproscan results
singularity exec agat_0.8.0--pl5262hdfd78af_0.sif agat_sp_manage_functional_annotation.pl -f stringtie_merged_v2_sub_cor2.gff -i interpro_stringtiev2_sub_cor_noframe_noReactomeMetaCyc.tsv > stringtie_merged_v2_sub_cor_ipr.gff

# remove agat log lines from gff output
# annotate with blast results from swissprot Drosophila
singularity exec agat_0.8.0--pl5262hdfd78af_0.sif agat_sp_manage_functional_annotation.pl -f stringtie_merged_v2_sub_cor_ipr.gff -b uniprot_sprot.Dmelanogaster.blastp.tophits_noframe_outfmt6.tsv -d swissprot_drosophila.fa > stringtie_merged_v2_sub_cor_ipr_sd.gff

# remove agat log lines from gff output
# annotate with blast results from swissprot Insects
singularity exec agat_0.8.0--pl5262hdfd78af_0.sif agat_sp_manage_functional_annotation.pl -f stringtie_merged_v2_sub_cor_ipr_sd.gff -b uniprot_sprot.Insects.blastp.tophits_noframe_outfmt6.tsv -d swissprot_insects.fa > stringtie_merged_v2_sub_cor_ipr_sd_si.gff

# remove agat log lines from gff output
# annotate with blast results from trembl Insect
singularity exec agat_0.8.0--pl5262hdfd78af_0.sif agat_sp_manage_functional_annotation.pl -f stringtie_merged_v2_sub_cor_ipr_sd_si.gff -b uniprot_trembl.Insects.blastp.tophits_wnames_noframe_outfmt6.tsv -d trembl_insects.fa > stringtie_merged_v2_sub_cor_ipr_sd_si_ti.gff

mv stringtie_merged_v2_sub_cor_ipr_sd_si_ti.gff stringtie_merged_v2_sub_cor_annotated.gff
```
