#gtf2 <- read.table("/home/celia/Downloads/Inversions_2022/gencode.v40.annotation.gtf", header = TRUE, sep = '\t')
getwd()
setwd("/home/celia/Downloads/Inversions_2022")
library("rtracklayer")
library("dplyr")
library("bedr")

gtf <- rtracklayer::import("/home/celia/Downloads/Inversions_2022/gencode.v40.annotation.gtf")
gtf_df=as.data.frame(gtf)

##keep only gene_annotation
gtf_df1<-gtf_df[gtf_df$type == "gene", ]

#keep only specific columns
gene_annotation<-gtf_df1 %>%
  select(seqnames, start, end, gene_id)

#gene_annotation<-gtf_df1 %>%
  #select(seqnames, start, end, gene_name)



##read the inversions bed file
bed <- as.data.frame(read.table("/home/celia/Downloads/Inversions_2022/inversions_recurrent_nonrecurrent/inversions_recurrent_nonrecurrent/processed/recurrent_invs.bed",
                                header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

bed_non_rec <- as.data.frame(read.table("/home/celia/Downloads/Inversions_2022/inversions_recurrent_nonrecurrent/inversions_recurrent_nonrecurrent/processed/nonrecurrent_invs.bed",
                                header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))


bed$inversion <- NA
bed$inversion <- paste0('inversion', c(1:length(bed$V2)))

bed_non_rec$inversion <- NA
bed_non_rec$inversion <- paste0('inversion', c(1:length(bed_non_rec$V2)))

names(bed_non_rec)[names(bed_non_rec) == 'V1'] <- 'chr'
names(bed_non_rec)[names(bed_non_rec) == 'V2'] <- 'start'
names(bed_non_rec)[names(bed_non_rec) == 'V3'] <- 'end'

names(gene_annotation)[names(gene_annotation) == 'seqnames'] <- 'chr'

#intersect function, seperately for recurrent and non-recurrent inversions
#bedtools_intersect <- function(inversions, gene_annotation) {
  #bedtools intersect -a inversions -b gene_annotation
#}


bedtools_intersect <- function(gene_annotation, inversions) {
 if (check.binary("bedtools")) {
  gene_annotation[,1]<-as.character(gene_annotation[,1])
  is.valid.region(gene_annotation)
  is.valid.region(bed_non_rec)
  gene_annotation.sort <- bedr.sort.region(gene_annotation)
  bed.sort <- bedr.sort.region(bed_non_rec)
  gene_annotation.merge <- bedr.merge.region(gene_annotation.sort)
  #a.sub1 <- bedr.subtract.region(gene_annotation.merge, bed.sort)
  #inversions$inversion <- NA
  #inversions$inversion <- paste0('inversion', c(1:length(inversions)))
  gene_overlap <- bedr.join.region(bed.sort, gene_annotation.sort)
  return(gene_overlap);
 }
 }


#keep only gene ids to pass them to GO
gene_annotation_overlap_only<- gene_overlap[, c("gene_id")]
gene_annotation_overlap_only<-as.data.frame(gene_annotation_overlap_only)
write.csv(gene_annotation_overlap_only, file="/home/celia/Downloads/Inversions_2022/non_recurrent.csv", row.names = FALSE)

#READ FILE WITH ALL INVERSIONS
inversions<-read.table("/home/celia/Downloads/IGV/igv_4celia/igv_4celia/invs/invs_all.bed", header = FALSE, sep = '\t')
names(inversions)[names(inversions) == 'V1'] <- 'chr'
names(inversions)[names(inversions) == 'V2'] <- 'start'
names(inversions)[names(inversions) == 'V3'] <- 'end'
inversions$inversion <- NA
inversions$inversion <- paste0('inversion', c(1:length(inversions$start)))
#remove column names
names(inversions)<-NULL
inversions <- inversions[ -c(4) ]
#write.csv(inversions, file="/home/celia/Downloads/Inversions_2022/inv_for_system.csv", row.names = FALSE)
write.table(inversions, "/home/celia/Downloads/Inversions_2022/inv_for_system.bed",sep = '\t', row.names=FALSE)


#turn dataframe into bed file
#names(gene_annotation)<-NULL
#write.table(gene_annotation, "/home/celia/Downloads/Inversions_2022/gen_annot.bed")

#read genome csv with chromosome sizes
hg38 <- read.csv(file = "/home//celia/Downloads/Inversions_2022/hg38.csv", header=FALSE)
write.table(hg38, "/home/celia/Downloads/Inversions_2022/hg38.bed", sep = '\t', row.names=FALSE)


#run bedtools
system("bedtools shuffle -i /home/celia/Downloads/Inversions_2022/inv_for_system.bed -g /home/celia/Downloads/Inversions_2022/hg38.bed -chrom > output")
  

  
  