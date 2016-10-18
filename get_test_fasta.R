source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings", suppressUpdates=TRUE)
biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates=TRUE)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
biocLite("biomaRt", suppressUpdates=TRUE)
library(biomaRt)
install.packages("seqinr", repos = "http://cran.r-project.org")
library(seqinr)
hg38 = BSgenome.Hsapiens.UCSC.hg38
seqnames(hg38)[1:25]
hg38.chrM = DNAStringSet(hg38$chrM)
names(hg38.chrM) = c("chrM")
writeXStringSet(hg38.chrM, file = "hg38_chrM.fa", format = "fasta")

mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
tss_map = getBM(attributes = c("hgnc_symbol", "chromosome_name", "strand", "transcript_start"), mart = mart)
tss_map = tss_map[which(tss_map$chromosome_name == "MT"), ]
start_relative2tss = -300
end_relative2tss = 0
start_pos = tss_map$transcript_start + tss_map$strand*start_relative2tss + 1;
end_pos = tss_map$transcript_start + tss_map$strand*end_relative2tss;
chroms = rep("chrM", times = length(start_pos))
wanted_ranges = GRanges(chroms, IRanges(apply(cbind(start_pos, end_pos), 1, min) , apply(cbind(start_pos, end_pos), 1, max)))
seqs = c()
for(i in 1:length(start_pos)){
	seqs = c(seqs, getSeq(hg38, wanted_ranges[i], as.character=TRUE))
}
tss_seqs = data.frame(genes = tss_map$hgnc_symbol, seqs = seqs, tss = tss_map$transcript_start)

write_seqs <- function(seqs, gene_names, tss, filename){
	stopifnot(dim(seqs)[1] == length(gene_names))
	write.fasta(file.out = filename, sequences = seqs[1], names = gene_names[1], open = "w", nbchar = 80, as.string = TRUE)
	if(length(gene_names) > 1){
	    for(i in 2:length(gene_names)){
	        write.fasta(file.out = filename, sequences = seqs[i], names = gene_names[i], open = "a", nbchar = 80, as.string = TRUE)
	   	}
	}
}

write_seqs(seqs = tss_seqs$seqs, gene_names = tss_seqs$genes, tss = tss_seqs$tss, filename = "hg38_chrM_minus300toTSS.fa")

	