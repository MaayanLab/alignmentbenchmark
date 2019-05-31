library("plyr")
library("stringr")
setwd("~/OneDrive/alignmentbench/")

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


gene_count_kallisto = getKallisto("SRR827478", human_map)
gene_count_salmon = getSalmon("SRR886587", human_map)
gene_count_hisat2 = getHISAT2("SRR901183", human_map)
gene_count_star = getSTAR("SRR827478", human_map)

inter1 = intersect(names(gene_count_kallisto), names(gene_count_salmon))
inter2 = intersect(names(gene_count_hisat2), names(gene_count_star))
inter = sort(intersect(inter1, inter2))

kallisto = gene_count_kallisto[inter]
salmon = gene_count_salmon[inter]
hisat2 = gene_count_hisat2[inter]
star = gene_count_star[inter]

expression = do.call(cbind, list(kallisto, salmon, hisat2, star))
colnames(expression) = c("kallisto","Salmon","HISAT2","STAR")


gg = grep("A.*\\.[0-9]$", rownames(expression))
expression = expression[-gg,]


scatterplotMatrix(log2(expression[sample(1:nrow(expression), 2000),]+1), pch=".",
    regLine = list(method=lm, lty=1, lwd=2, col="red"),
    col="black")


genequant<-function(abundance){
    
    #args <- commandArgs(TRUE)
    
    #mapping = args[1]
    #res = load(mapping)
    
    kf = "quant/kallisto/SRR1002568/abundance.tsv"
    sf = "quant/salmon/SRR1002568/quant.sf"
    
    kallisto = read.table(kf, sep="\t", stringsAsFactors=F, header=T)
    salmon = read.table(sf, sep="\t", stringsAsFactors=F, header=T)
    salmon = salmon[,c(1,2,3,5,4)]
    
    inter = intersect(kallisto[,1], salmon[,1])
    sdd = setdiff(kallisto[,1], salmon[,1])
    
    ugene = cb[,2]
    
    m3 = match(abu[,1], cb[,1])
    cco = cbind(abu,ugene[m3])[-1,]
    co = cco[,c(6,4)]
    co[,1] = as.character(co[,1])
    df = data.frame(co[,1], as.numeric(co[,2]))
    colnames(df) = c("gene", "value")
    
    dd = ddply(df,.(gene),summarize,sum=sum(value),number=length(gene))
    ge = dd[,2]
    names(ge) = dd[,1]
    
}

library(RColorBrewer)
darkcols <- brewer.pal(8, "Set1")


times = read.table("supportfiles/SRR7460337_time.tsv", sep="\t", stringsAsFactors=F)
rownames(times) = times[,1]
times = times[,-1]
colnames(times) = c(1,2,3,4,6,8,12,16)

oo = order(rowSums(times))
oo = c(2,1,4,3,5)

pdf("figures/speed_align_SRR7460337.pdf")
par(mar=c(6,6,6,6))
barplot(as.matrix(times[oo,]), beside=TRUE, col=darkcols[1:5], xlab="threads", ylab="seconds", cex.lab=1.8, cex.axis=2, cex.names=1.8)
legend("topright",legend=rownames(times)[oo], fill=darkcols[1:5], bty = "n", cex=2)
dev.off()


times = read.table("supportfiles/SRR7460337_time.tsv", sep="\t", stringsAsFactors=F)
rownames(times) = times[,1]
times = times[,-1]
colnames(times) = c(1,2,3,4,6,8,12,16)

times = times[-5,]

oo = c(2,1,4,3)


pdf("figures/speed_align_SRR7460337_nobwa.pdf")
par(mar=c(6,6,6,6))
barplot(as.matrix(times[oo,]), beside=TRUE, col=darkcols[1:4], xlab="threads", ylab="seconds", cex.lab=1.8, cex.axis=2, cex.names=1.8)
legend("topright",legend=rownames(times)[oo], fill=darkcols[1:4], bty = "n", cex=2)
dev.off()








times = read.table("supportfiles/SRR2972202_time.tsv", sep="\t", stringsAsFactors=F)
rownames(times) = times[,1]
times = times[,-1]
colnames(times) = c(1,2,3,4,6,8,12,16)

oo = order(rowSums(times))
oo = c(2,1,4,3,5)

pdf("figures/speed_align_SRR2972202.pdf")
par(mar=c(6,6,6,6))
barplot(as.matrix(times[oo,]), beside=TRUE, col=darkcols[1:5], xlab="threads", ylab="seconds", cex.lab=1.8, cex.axis=2, cex.names=1.8)
legend("topright",legend=rownames(times)[oo], fill=darkcols[1:5], bty = "n", cex=2)
dev.off()


times = read.table("supportfiles/SRR2972202_time.tsv", sep="\t", stringsAsFactors=F)
rownames(times) = times[,1]
times = times[,-1]
colnames(times) = c(1,2,3,4,6,8,12,16)

times = times[-5,]

oo = c(2,1,4,3)


pdf("figures/speed_align_SRR2972202_nobwa.pdf")
par(mar=c(6,6,6,6))
barplot(as.matrix(times[oo,]), beside=TRUE, col=darkcols[1:4], xlab="threads", ylab="seconds", cex.lab=1.8, cex.axis=2, cex.names=1.8)
legend("topright",legend=rownames(times)[oo], fill=darkcols[1:4], bty = "n", cex=2)
dev.off()




pdf("figures/speed_align_SRR7460337_SRR2972202.pdf", 14, 10)
par(mfrow=c(1,2))


par(mar=c(5,5,4,1))

times = read.table("supportfiles/SRR7460337_time.tsv", sep="\t", stringsAsFactors=F)
rownames(times) = times[,1]
times = times[,-1]
colnames(times) = c(1,2,3,4,6,8,12,16)

oo = order(rowSums(times))
oo = c(2,1,4,3,5)

barplot(as.matrix(times[oo,]), beside=TRUE, col=darkcols[1:5], xlab="threads", ylab="seconds", cex.lab=2, cex.axis=2, cex.names=1.8)
legend("topright",legend=rownames(times)[oo], fill=darkcols[1:5], bty = "n", cex=2)


times = read.table("supportfiles/SRR2972202_time.tsv", sep="\t", stringsAsFactors=F)
rownames(times) = times[,1]
times = times[,-1]
colnames(times) = c(1,2,3,4,6,8,12,16)

oo = order(rowSums(times))
oo = c(2,1,4,3,5)

barplot(as.matrix(times[oo,]), beside=TRUE, col=darkcols[1:5], xlab="threads", ylab="seconds", cex.lab=2, cex.axis=2, cex.names=1.8)
legend("topright",legend=rownames(times)[oo], fill=darkcols[1:5], bty = "n", cex=2)

dev.off()





