args = commandArgs(trailingOnly=TRUE)
sraid = args[1]
species = args[2]

library("plyr")
library("stringr")

load("supportfiles/gene_mapping.rda")

getKallisto <- function(sra, mapping){
    file = paste0("quant/kallisto/",sra,"/abundance.tsv")
    input = read.table(file, sep="\t", stringsAsFactors=F, header=T)[,c(1,4)]
    input[,1] = gsub("\\.[0-9]$","", input[,1])
    
    gene_match = match(input[,1], mapping[,1])
    
    transcript_counts = data.frame(mapping[gene_match,2], round(as.numeric(input[,2])))
    colnames(transcript_counts) = c("gene", "value")
    
    dd = ddply(transcript_counts,.(gene),summarize,sum=sum(value),number=length(gene))
    ge = dd[,2]
    names(ge) = dd[,1]
    return(ge)
}

getSalmon <- function(sra, mapping){
    file = paste0("quant/salmon/",sra,"/quant.sf")
    input = read.table(file, sep="\t", stringsAsFactors=F, header=T)[,c(1,5)]
    input[,1] = gsub("\\.[0-9]$","", input[,1])
    
    gene_match = match(input[,1], mapping[,1])
    
    transcript_counts = data.frame(mapping[gene_match,2], round(as.numeric(input[,2])))
    colnames(transcript_counts) = c("gene", "value")
    
    dd = ddply(transcript_counts,.(gene),summarize,sum=sum(value),number=length(gene))
    ge = dd[,2]
    names(ge) = dd[,1]
    return(ge)
}

getHISAT2 <- function(sra, mapping){
    file = paste0("quant/hisat2/",sra,"/",sra,".tsv")
    input = read.table(file, sep="\t", stringsAsFactors=F, header=T)[,c(1,7)]
    
    gene_map = match(input[,1], mapping[,3])
    ww = which(!is.na(gene_map))
    ge = input[ww,2]
    names(ge) = mapping[gene_map[ww],2]
    return(ge)
}

getSTAR <- function(sra, mapping){
    file = paste0("quant/star/",sra,"/",sra,"ReadsPerGene.out.tab")
    input = read.table(file, sep="\t", stringsAsFactors=F, skip=4)[,c(1,2)]
    
    gene_map = match(input[,1], mapping[,3])
    ww = which(!is.na(gene_map))
    ge = input[ww,2]
    names(ge) = mapping[gene_map[ww],2]
    return(ge)
}

getBWA <- function(sra, mapping){
    file = paste0("quant/bwa/",sra,"/",sra,".tsv")
    input = read.table(file, sep="\t", stringsAsFactors=F, header=T)[,c(1,7)]
    
    gene_map = match(input[,1], mapping[,3])
    ww = which(!is.na(gene_map))
    ge = input[ww,2]
    names(ge) = mapping[gene_map[ww],2]
    return(ge)
}

gene_map = human_map
if(species == "mouse"){
    gene_map = mouse_map
}

gene_count_kallisto = getKallisto(sraid, gene_map)
gene_count_salmon = getSalmon(sraid, gene_map)
gene_count_hisat2 = getHISAT2(sraid, gene_map)
gene_count_star = getSTAR(sraid, gene_map)
gene_count_bwa = getBWA(sraid, gene_map)

inter1 = intersect(names(gene_count_kallisto), names(gene_count_salmon))
inter2 = intersect(names(gene_count_hisat2), names(gene_count_star))
inter3 = intersect(inter2, names(gene_count_bwa))
inter = sort(intersect(inter1, inter3))

kallisto = gene_count_kallisto[inter]
salmon = gene_count_salmon[inter]
hisat2 = gene_count_hisat2[inter]
star = gene_count_star[inter]
bwa = gene_count_bwa[inter]

expression = do.call(cbind, list(kallisto, salmon, hisat2, star, bwa))
colnames(expression) = c("kallisto","Salmon","HISAT2","STAR", "BWA")

write.table(expression, file=paste0("quant/combined/",sraid,".tsv"), sep="\t", quote=F)

