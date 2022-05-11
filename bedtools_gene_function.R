#gtf2 <- read.table("/home/celia/Downloads/Inversions_2022/gencode.v40.annotation.gtf", header = TRUE, sep = '\t')
library("rtracklayer")
library("dplyr")
library("bedr")

gtf <- rtracklayer::import("/home/celia/Downloads/Inversions_2022/gencode.v40.annotation.gtf")
gtf_df=as.data.frame(gtf)

##keep only genes
gtf_df1<-gtf_df[gtf_df$type == "gene", ]

#keep only specific columns
gene_annotation<-gtf_df1 %>%
  select(seqnames, start, end, gene_id)

#turn dataframe into bed file
write.table(gene_annotation, "/home/celia/Downloads/Inversions_2022/gen_annot.bed")




##read the inversions bed file
bed <- as.data.frame(read.table("/home/celia/Downloads/Inversions_2022/inversions_recurrent_nonrecurrent/inversions_recurrent_nonrecurrent/processed/recurrent_invs.bed",
                                header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

bed$inversion <- NA
bed$inversion <- paste0('inversion', c(1:length(bed)))

names(bed)[names(bed) == 'V1'] <- 'seqnames'
names(bed)[names(bed) == 'V2'] <- 'start'
names(bed)[names(bed) == 'V3'] <- 'end'


#intersect function, seperately for recurrent and non-recurrent inversions
#bedtools_intersect <- function(inversions, genes) {
  #bedtools intersect -a inversions -b genes
#}


bedtools_intersect <- function(inversions, genes) {
 if (check.binary("bedtools")) {
  #inversions$inversion <- NA
  #inversions$inversion <- paste0('inversion', c(1:length(inversions)))
  a.int1 <- bedr.join.region(inversions, genes)
  return(a.int1);
 }
}
