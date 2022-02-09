### Whole genome sequencing data analysis ###
# code updated at: https://github.com/sungminhwang-duke/breseq_Yarrowia

### Quality check for the sequencing data
# gunzip *.gz
# trim_galore -o ./trimming *.fastq --paired --fastqc --illumina

### Revision of the reference genome 
#CLIB89 (W29): https://www.ncbi.nlm.nih.gov/genome/194?genome_assembly_id=285090
#Representative genome of Yarrowia: GCF_000002525.2_ASM252v1_genomic.gbff
# *Removed "Repeat feature" (YalifMp02, YalifMp26) due to the error below:
#      Repeat feature must have only one location: YalifMp02
#      Please rename the duplicates so that they have unique locus_tags & names, e.g., gene_1 and gene_2

### Run breseq by manual
# (done) breseq -p -r GCA_001761485.1_ASM176148v1_genomic.gbff Ctrl_1_val_1.fq Ctrl_2_val_2.fq -o breseq_output      #took 3h 10m
# (done) breseq -p -r GCA_001761485.1_ASM176148v1_genomic.gbff 7-1_1_val_1.fq 7-1_2_val_2.fq -o breseq_output_7-1    #took 6h 41m
# (done) breseq -p -r GCA_001761485.1_ASM176148v1_genomic.gbff 51_1_val_1.fq 51_2_val_2.fq -o breseq_output_51       #took 5h 38m

#for f1 in *1_val_1.fq ; do for f2 in ${f1%%_1_val_1.fq}"_2_val_2.fq" ; do breseq -p -r GCA_001761485.1_ASM176148v1_genomic.gbff $f1 $f2 -o $f1.output ; done; done


### Identify the genomic differences by comparing each genome
library(dplyr)
library(tidyverse)

res.ctrl.0 <- read.csv("z_output_1_notcutwithpairedend/Raw_breseq_result_Ctrl.csv", header = T, fileEncoding = "latin1")
res.ctrl.1 <- res.ctrl.0 %>% filter(freq > 0.8)  #frequency higher than 80%
res.ctrl.1 <- res.ctrl.1 %>% unite(Ctrl, seq_id, position, remove=F) 
write.csv(res.ctrl.1, "z_output_1_notcutwithpairedend/Raw_breseq_result_Ctrl_add.csv")
res.ctrl.2 <- res.ctrl.1 %>% select(Ctrl)
res.ctrl.2$number <- row.names(res.ctrl.2)
  
  
res.7_1.0 <- read.csv("z_output_1_notcutwithpairedend/Raw_breseq_result_7-1.csv", header = T, fileEncoding = "latin1")
res.7_1.1 <- res.7_1.0 %>% filter(freq > 0.8)
res.7_1.1 <- res.7_1.1 %>% unite(S7_1, seq_id, position, remove=F) 
write.csv(res.7_1.1, "z_output_1_notcutwithpairedend/Raw_breseq_result_7-1_add.csv")
res.7_1.2 <- res.7_1.1 %>% select(S7_1)
res.7_1.2$number <- row.names(res.7_1.2)


res.51.0 <- read.csv("z_output_1_notcutwithpairedend/Raw_breseq_result_51.csv", header = T, fileEncoding = "latin1")
res.51.1 <- res.51.0 %>% filter(freq > 0.8)
res.51.1 <- res.51.1 %>% unite(S51, seq_id, position, remove=F) 
write.csv(res.51.1, "z_output_1_notcutwithpairedend/Raw_breseq_result_51_add.csv")
res.51.2 <- res.51.1 %>% select(S51)
res.51.2$number <- row.names(res.51.2)


two.data <- merge(res.ctrl.2, res.7_1.2, by="number", all=T)
three.data <- merge(two.data, res.51.2, by="number", all=T)

write.csv(three.data, "combined.csv") #Draw venn diagram and find the overlaps in here: http://bioinformatics.psb.ugent.be/webtools/Venn/


### Identify the different loci/genes discovered by the genome comparison
diff.data <- read.csv("z_output_1_notcutwithpairedend/combined_2.csv")
unique(diff.data$Names) #"Ctrl S51 S7-1"   "Ctrl S7-1"   "Ctrl S51"   "S51 S7-1"   "Ctrl"   "S7-1"   "S51"
 
Ctrl.S51.S7 <- diff.data %>% filter( Names == "Ctrl S51 S7-1") 
Ctrl.S51.S7$gene <- res.ctrl.1$gene[match(Ctrl.S51.S7$elements, res.ctrl.1$Ctrl)]
Ctrl.S51.S7$description <- res.ctrl.1$description[match(Ctrl.S51.S7$elements, res.ctrl.1$Ctrl)]

Ctrl.S7 <- diff.data %>% filter( Names == "Ctrl S7-1") 
Ctrl.S7$gene <- res.ctrl.1$gene[match(Ctrl.S7$elements, res.ctrl.1$Ctrl)]
Ctrl.S7$description <- res.ctrl.1$description[match(Ctrl.S7$elements, res.ctrl.1$Ctrl)]

Ctrl.S51 <- diff.data %>% filter( Names == "Ctrl S51") 
Ctrl.S51$gene <- res.ctrl.1$gene[match(Ctrl.S51$elements, res.ctrl.1$Ctrl)]
Ctrl.S51$description <- res.ctrl.1$description[match(Ctrl.S51$elements, res.ctrl.1$Ctrl)]

Ctrl. <- diff.data %>% filter( Names == "Ctrl") 
Ctrl.$gene <- res.ctrl.1$gene[match(Ctrl.$elements, res.ctrl.1$Ctrl)]
Ctrl.$description <- res.ctrl.1$description[match(Ctrl.$elements, res.ctrl.1$Ctrl)]

S51.S7 <- diff.data %>% filter( Names == "S51 S7-1") 
S51.S7$gene <- res.51.1$gene[match(S51.S7$elements, res.51.1$S51)]
S51.S7$description <- res.51.1$description[match(S51.S7$elements, res.51.1$S51)]

S51. <- diff.data %>% filter( Names == "S51") 
S51.$gene <- res.51.1$gene[match(S51.$elements, res.51.1$S51)]
S51.$description <- res.51.1$description[match(S51.$elements, res.51.1$S51)]

S7. <- diff.data %>% filter( Names == "S7-1") 
S7.$gene <- res.7_1.1$gene[match(S7.$elements, res.7_1.1$S7_1)]
S7.$description <- res.7_1.1$description[match(S7.$elements, res.7_1.1$S7_1)]

gene.description.1 <- rbind(Ctrl.S51.S7, Ctrl.S7, Ctrl.S51, Ctrl., S51.S7, S51., S7.)
write.csv(gene.description.1, "z_output_1_notcutwithpairedend/combined_3_all.csv")

gene.description.2 <- rbind(Ctrl.S7, Ctrl.S51, Ctrl., S51.S7, S51., S7.)
write.csv(gene.description.2, "z_output_1_notcutwithpairedend/combined_3_removed_alloverlap.csv")

gene.description.3 <- rbind(S51., S7.)
write.csv(gene.description.3, "z_output_1_notcutwithpairedend/combined_3_specific_S51_S7.csv")






### Pair wp number with locus and old locus IDs, along with start and stop positions
library(tidyverse)
library(plyr)
library(magrittr)

rm(list = ls())

gff.parse <- function(x) {
  ids <- list()
  x %<>% filter(type == "gene" |type == "pseudogene" | type == "CDS" | type == "tRNA" | type == "rRNA" | type == "RNase_P_RNA" | type == "SRP_RNA")
  x %<>% arrange(desc(start))
  if ("old_locus_tag" %in% colnames(x)) {
    for (i in 1:(nrow(x)-1)) {
      if (all(x[i,2:4] == x[i+1,2:4])) {
        vec <- data.frame("chr"=x[["seqnames"]][i],
                          "acc"= x[["protein_id"]][(i+1)], 
                          "locus_tag" = x[["locus_tag"]][i], 
                          "old_locus_tag" = x[["old_locus_tag"]][i], 
                          "gene"= x[["gene"]][i],
                          "length_nt"= x[["width"]][i], 
                          "start"= x[["start"]][i], 
                          "end"= x[["end"]][i],
                          "strand" = x[["strand"]][i],
                          "type"= x[["gene_biotype"]][i],
                          "annotation"= x[["product"]][(i+1)],
                          stringsAsFactors = F)
        ids[[i]] <- vec
      }}} 
  else if ("gene" %in% colnames(x)) {
    for (i in 1:(nrow(x)-1)) {
      if (all(x[i,1:2] == x[i+1,1:2])) {
        vec <- data.frame("chr"=x[["seqnames"]][i],
                          "acc"= x[["protein_id"]][(i+1)], 
                          "locus_tag" = x[["locus_tag"]][i], 
                          "gene_name" = x[["gene"]][i], 
                          "length_nt"= x[["width"]][i], 
                          "start"= x[["start"]][i], 
                          "end"= x[["end"]][i],
                          "strand" = x[["strand"]][i],
                          "type"= x[["gene_biotype"]][i],
                          "annotation"= x[["product"]][(i+1)],
                          stringsAsFactors = F)
        ids[[i]] <- vec
      }}}
  
  df <- bind_rows(ids)
  if (length(na.omit(unique(x[["protein_id"]]))) == length(na.omit(unique(df[["acc"]])))) {
    print("all protein ids accounted for")
  }
  return(df)
}

getwd() #Path to find the .gff files
(filenames <- paste("0_materials/", list.files(path = "0_materials/", ".gff$"), sep = ""))

#Execute function for all gff files in working directory and write results out to new file
for (filename in filenames) {
  
  # read data:
  gtf <- rtracklayer::import.gff(filename)
  gtf_df <- as.data.frame(gtf)
  gtf_df %>% mutate_if(is.factor, as.character)-> gtf_df 
  
  # use fxn to make key
  key <- gff.parse(gtf_df)
  
  #write out to new .csv
  write_csv(key, paste(filename, "_key.csv", sep = ""))
}

anno.data <- read.csv("0_materials/GCF_000007805.1_ASM780v1_genomic.gff_key_3geneadded.csv")
#anno.data <- read.csv("0_materials/GCF_000007805.1_ASM780v1_genomic.gff_key.csv")
substi.data <- read.csv("6_DESeq2/0_Data.csv")

substi.data$locus <- anno.data$locus_tag[match(substi.data$gene_id, anno.data$gene)]
write.csv(substi.data, "anno_2.csv")
