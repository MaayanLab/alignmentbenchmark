setwd("~/OneDrive/alignmentbench")

library("R.utils")
library("RCurl")
library("DESeq2")
library("limma")
library("statmod")
library("edgeR")
library("enrichR")

diffExpression <- function(controls, samples){
    expression = do.call(cbind, append(controls, samples))
    
    design = c(rep(1,length(controls)), rep(-1, length(samples)))
    dm = model.matrix(~as.factor(design))
    
    dge <- DGEList(counts=expression)
    dge <- calcNormFactors(dge)
    
    v <- voom(dge, dm, plot=TRUE)
    
    fit <- lmFit(v, dm)
    fit <- eBayes(fit)
    #topTable(fit, coef=ncol(dm))
    
    r1 = fit$p.value[,2]
    r2 = fit$t[,2]
    
    return(r2)
}

signatures = read.table("supportfiles/signatures.tsv", stringsAsFactors=F, sep="\t", fill=T, quote="", header=T)
sampleMap = read.table("supportfiles/samplemapping.tsv", stringsAsFactors=F, sep="\t", fill=T, quote="", header=T)

whuman = which(signatures[,14] == "homo sapiens")

controlList = list()
perturbationList = list()

for(co in whuman){
    controlList[[length(controlList)+1]] = sort(unlist(c(strsplit(signatures[co, 6], ","), strsplit(signatures[co,10], ","))))
    perturbationList[[length(perturbationList)+1]] = sort(unlist(c(strsplit(signatures[co,7], ","), strsplit(signatures[co,9], ","))))
}

ww = which((lapply(perturbationList,length) == 3) & (lapply(controlList,length)  == 3))

controlList = controlList[ww]
perturbationList = perturbationList[ww]
signaturesHuman = signatures[whuman,][ww,]

perturbation_type_count = table(signaturesHuman[,11])



samples = unqiue(unlist(c(controlList, perturbationList)))
sampleCounts = list()

for(sample in samples){
    ws = which(sampleMap[,3] == sample)
    sras = sampleMap[ws,4]
    temp = matrix(0, 21933, 5)
    for(sra in sras){
        if(url.exists(paste0("https://s3.amazonaws.com/mssm-genecount-combined/",sra,".tsv"))){
            download.file(paste0("https://s3.amazonaws.com/mssm-genecount-combined/",sra,".tsv"), paste0("downloads/", sra, ".tsv"))
            rr = read.table(paste0("downloads/", sra, ".tsv"), stringsAsFactor=F, sep="\t")
            temp = temp + rr
        }
        else if(url.exists(paste0("https://s3.amazonaws.com/mssm-genecount-combined/",sra,".tsv.gz"))){
            download.file(paste0("https://s3.amazonaws.com/mssm-genecount-combined/",sra,".tsv.gz"), paste0("downloads/", sra, ".tsv.gz"))
            gunzip(paste0("downloads/", sra, ".tsv.gz"), overwrite=TRUE)
            rr = read.table(paste0("downloads/", sra, ".tsv"), stringsAsFactor=F, sep="\t")
            temp = temp + rr
        }
        else{
            print("file missing")
        }
    }
    
    if(sum(rr) > 10000){
        sampleCounts[[length(sampleCounts)+1]] = temp
        names(sampleCounts)[length(sampleCounts)] = sample
    }
}

countsKallisto = list()
countsSalmon = list()
countsHisat2 = list()
countsBwa = list()
countsStar = list()

for(i in 1:length(sampleCounts)){
    countsKallisto[[i]] = sampleCounts[[i]][,1]
    countsSalmon[[i]] = sampleCounts[[i]][,2]
    countsHisat2[[i]] = sampleCounts[[i]][,3]
    countsStar[[i]] = sampleCounts[[i]][,4]
    countsBwa[[i]] = sampleCounts[[i]][,5]
}

kallisto = do.call(cbind, countsKallisto)
rownames(kallisto) = rownames(rr)
colnames(kallisto) = names(sampleCounts)

salmon = do.call(cbind, countsSalmon)
rownames(salmon) = rownames(rr)
colnames(salmon) = names(sampleCounts)

hisat2 = do.call(cbind, countsHisat2)
rownames(hisat2) = rownames(rr)
colnames(hisat2) = names(sampleCounts)

star = do.call(cbind, countsStar)
rownames(star) = rownames(rr)
colnames(star) = names(sampleCounts)

bwa = do.call(cbind, countsBwa)
rownames(bwa) = rownames(rr)
colnames(bwa) = names(sampleCounts)

ll = list(kallisto, salmon, hisat2, star, bwa)
names(ll) = c("kallisto", "Salmon", "HISAT2", "STAR", "BWA")


pdf("correlation_genecount.pdf")
par(mfrow=c(2,5))
par(mar=c(1,1,1,1))
par(oma=c(3,3,3,3))

cx = cor(ll[[1]], method="spearman")
plot(density(cx[upper.tri(cx)], na.rm=T), lwd=3, xlab="Spearman correlation", cex.axis=2, cex.lab=2, main="")

cx = cor(ll[[2]], method="spearman")
plot(density(cx[upper.tri(cx)], na.rm=T), lwd=3, xlab="Spearman correlation", cex.axis=2, cex.lab=2, main="")

cx = cor(ll[[3]], method="spearman")
plot(density(cx[upper.tri(cx)], na.rm=T), lwd=3, xlab="Spearman correlation", cex.axis=2, cex.lab=2, main="")

cx = cor(ll[[4]], method="spearman")
plot(density(cx[upper.tri(cx)], na.rm=T), lwd=3, xlab="Spearman correlation", cex.axis=2, cex.lab=2, main="")

cx = cor(ll[[5]], method="spearman")
plot(density(cx[upper.tri(cx)], na.rm=T), lwd=3, xlab="Spearman correlation", cex.axis=2, cex.lab=2, main="")

cx = cor(diffList[[1]], diffList[[1]], method="pearson")
cx[cx > 0.999] = 0
plot(density(cx[upper.tri(cx)], na.rm=T), lwd=3, xlab="Pearson correlation", cex.axis=2, cex.lab=2, main="")

cx = cor(diffList[[2]], diffList[[2]], method="pearson")
cx[cx > 0.999] = 0
plot(density(cx[upper.tri(cx)], na.rm=T), lwd=3, xlab="Pearson correlation", cex.axis=2, cex.lab=2, main="")

cx = cor(diffList[[3]], diffList[[3]], method="pearson")
cx[cx > 0.999] = 0
plot(density(cx[upper.tri(cx)], na.rm=T), lwd=3, xlab="Pearson correlation", cex.axis=2, cex.lab=2, main="")

cx = cor(diffList[[4]], diffList[[4]], method="pearson")
cx[cx > 0.999] = 0
plot(density(cx[upper.tri(cx)], na.rm=T), lwd=3, xlab="Pearson correlation", cex.axis=2, cex.lab=2, main="")

cx = cor(diffList[[5]], diffList[[5]], method="pearson")
cx[cx > 0.999] = 0
plot(density(cx[upper.tri(cx)], na.rm=T), lwd=3, xlab="Pearson correlation", cex.axis=2, cex.lab=2, main="")


dev.off()



cc = cor(kallisto, salmon)
ks = colSums(kallisto)

dd = diag(cc)
ww = which(!is.na(dd))

plot(log(ks[ww]), dd[ww])

plot(density(dd[ww]), lwd=3)
plot((dd[ww]), lwd=2)


dex = do.call(cbind, ll)
wwc = which(rowSums(dex) > 50)

mc = c(0.770, 0.770, 0.814, 0.808, 0.717)

pdf("figures/absolute_correlation.pdf")
par(oma=c(3,2,3,2))
par(mfrow=c(5,5))
par(mar=c(0.2,0.2,0.2,0.2))

for(i in c(1,2,4,3,5)){
    
    ddx = c()
    
    for(j in c(1,2,4,3,5)){
        if(i == j){
            plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
            text(x = 0.5, y = 0.4, names(ll)[i], cex=2)
            text(x = 0.5, y = 0.22, mc[i], cex=1.5)
        }
        else{
            cc = cor(ll[[i]][wwc,], ll[[j]][wwc,], method="spearman")
            dd = diag(cc)
            ww = which(!is.na(dd))
            ddx = c(ddx, dd[ww])
            de = density(dd[ww])
            de$y = de$y/max(de$y)
            de2 = density(cc[upper.tri(cc)], na.rm=T)
            de2$y = de2$y/max(de2$y)
            plot(de, lwd=3, main="", xlab="", ylab="", ylim=c(0,1.3), xlim=c(0,1), yaxt = 'n', frame.plot=F, cex.axis=1.4)
            polygon(de, col="grey", border="black")
            lines(de2, col="aquamarine", lwd=5)
        }
    }
    
    print(paste(names(ll)[i], mean(ddx)))
}
dev.off()

diffList = list()

for(i in c(1,2,4,3,5)){
    
    tmat = ll[[i]]
    wk = which(colSums(tmat) > 0)
    tmat = tmat[,wk]
    
    diffy = list()
    
    for(i in 1:nrow(signaturesHuman)){
        controlCounts = list()
        for(sample in controlList[[i]]){
            if(sample %in% colnames(tmat)){
                controlCounts[[length(controlCounts)+1]] = tmat[,sample]
            }
        }
        
        perturbationCounts = list()
        for(sample in perturbationList[[i]]){
            if(sample %in% colnames(tmat)){
                perturbationCounts[[length(perturbationCounts)+1]] = tmat[,sample]
            }
        }
        if(length(controlCounts) > 1 & length(perturbationCounts) > 1){
            diffy[[length(diffy)+1]] = diffExpression(controlCounts, perturbationCounts)
            names(diffy)[length(diffy)] = i
        }
    }
    
    diffList[[length(diffList)+1]] = do.call(cbind, diffy)
}

names(diffList) = names(ll)[c(1,2,4,3,5)]



mc = c(0.748, 0.744, 0.761, 0.771, 0.697)

pdf("figures/diff_expr.pdf")
par(oma=c(3,2,3,2))
par(mfrow=c(5,5))
par(mar=c(0.2,0.2,0.2,0.2))

for(i in c(1,2,3,4,5)){
    
    ddx = c()
    
    for(j in c(1,2,3,4,5)){
        if(i == j){
            plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
            text(x = 0.5, y = 0.4, names(diffList)[i], cex=2)
            text(x = 0.5, y = 0.22, mc[i], cex=1.5)
        }
        else{
            cc = cor(diffList[[i]], diffList[[j]], method="pearson")
            dd = diag(cc)
            ww = which(!is.na(dd))
            ddx = c(ddx, dd[ww])
            de = density(dd[ww])
            de$y = de$y/max(de$y)
            de2 = density(cc[upper.tri(cc)], na.rm=T)
            de2$y = de2$y/max(de2$y)
            plot(de, lwd=3, main="", xlab="", ylab="", ylim=c(0,1.3), xlim=c(-0.2,1), yaxt = 'n', frame.plot=F, cex.axis=1.4)
            polygon(de, col="grey", border="black")
            lines(de2, col="aquamarine", lwd=5)
        }
    }
    
    print(paste(names(ll)[i], mean(ddx)))
}
dev.off()

updown = list()

for(i in 1:5){
    upgc = c()
    downgc = c()
    for(j in 1:ncol(diffList[[i]])){
        wwu = which(diffList[[i]][,j] > 5)
        wwd = which(diffList[[i]][,j] < -5)
        print(paste(length(wwu), length(wwd)))
        upgc = c(upgc, length(wwu))
        downgc = c(downgc, length(wwd))
    }
    
    updown[[length(updown)+1]] = downgc
    updown[[length(updown)+1]] = upgc
}

bs = do.call(cbind, updown)

pdf("figures/diff_gene_count.pdf", 9, 7)
par(mar=c(6,6,6,6))
boxplot(updown, lwd=2, outline=F, col=c("aquamarine", "lightpink"), ylab="differntial expressed genes", labels="", cex.lab=1.6, cex.axis=1.6, xaxt="n", ylim=c(-100, 1600))

text(x=1.5, y=-80, labels="kallisto", cex=1.6)
text(x=3.5, y=-80, labels="Salmon", cex=1.6)
text(x=5.5, y=-80, labels="STAR", cex=1.6)
text(x=7.5, y=-80, labels="HISAT2", cex=1.6)
text(x=9.5, y=-80, labels="BWA", cex=1.6)

abline(v=(0.5+2*(1:4)), lwd=1)
legend("topleft", legend=c("up-regulated", "down-regulated"), fill=c("aquamarine", "lightpink"), cex=1.4, bg="white")
dev.off()

sample_count = list()
for(i in 1:5){
    cs = colSums(ll[[i]])
    cs = cs[cs !=0]
    sample_count[[length(sample_count)+1]] = cs
    
}
names(sample_count) = names(ll)

sample_count = sample_count[c(1,2,4,3,5)]

scm = do.call(cbind, sample_count)

log2mean<- function(x) {
    x = x[x != 0]
    x <- log2(1+x)
    return = mean(x)
}

pdf("figures/sample_count_trend.pdf", 9, 7)
par(mar=c(6,6,6,6))
plot(1:5, log2(scm[1,]+1), type="b", pch=20, ylim=c(20, 26), cex.lab=1.6, cex.axis=1.6, xlab="", ylab="log2 sum of gene counts", xaxt="n")
for(i in 2:nrow(scm)){
    lines(1:5, log2(scm[i,]+1), type="b", pch=20, col=rgb(0,0,0,0.2))
}
axis(1, 1:5, names(sample_count), cex.axis=1.6)
lines(1:5, unlist(lapply(sample_count, log2mean)), type="b", col="red", lwd=4, pch=20, cex=2)
legend("topleft", legend=c("sample", "mean trend"), pch=20, col=c("grey", "red"), lwd=c(1,3), border=F, cex=2)
dev.off()



wwk = which(scm[,1] > scm[,5])
plot(1:5, log2(scm[wwk,][1,]+1), type="b", pch=20, ylim=c(20, 26))

for(i in 2:length(wwk)){
    lines(1:5, log2(scm[wwk,][i,]+1), type="b", pch=20)
}



nka = normalize.quantiles(kallisto)
dimnames(nka) =dimnames(kallisto)
nka = log2(1+nka)

nhi = normalize.quantiles(hisat2)
dimnames(nhi) = dimnames(kallisto)
nhi = log2(1+nhi)

nbw = normalize.quantiles(bwa)
dimnames(nbw) = dimnames(kallisto)
nbw = log2(1+nbw)

nst = normalize.quantiles(star)
dimnames(nst) = dimnames(kallisto)
nst = log2(1+nst)

nsa = normalize.quantiles(salmon)
dimnames(nsa) = dimnames(kallisto)
nsa = log2(1+nsa)

nexp = list(nka, nsa, nst, nhi, nbw)
nnexp = do.call(cbind, nexp)

nnnexp = normalize.quantiles(nnexp)

wval = which(rowMeans(nnnexp) > 0.1)
nnnexp = nnnexp[wval, ]

qq = quantile(probs = seq(0, 1, 0.1), rowMeans(nnnexp))
wlow = which(rowMeans(nnnexp) > qq[2] & rowMeans(nnnexp) < qq[5])
whigh = which(rowMeans(nnnexp) > qq[8] & rowMeans(nnnexp) < qq[11])


h = hist(rowMeans(nnnexp), breaks=qq)

pdf("figures/percentilerange.pdf", 8,6)
par(mar=c(6,6,6,6))
plot(density(rowMeans(nnnexp)), lwd=1, cex.lab=1.6, cex.axis=1.6, xlab="mean log2 gene count", ylab="Density", main="")
rect(qq[2],0,qq[5],100, col=rgb(1,0,0,0.3))
rect(qq[8],0,qq[11],100, col=rgb(0,1,0,0.3))
lines(density(rowMeans(nnnexp)), lwd=6)
legend("topright",cex=1, bg="white", legend=c("low expression genes  (10% to 40% percentile range)", "high expression genes (70% to 100% percentile range)"), fill=c(rgb(1,0,0,0.3), rgb(0,1,0,0.3)))

dev.off()


pdf("figures/gene_correlation.pdf")
par(mfrow=c(5,5))
par(mar=c(0.2, 0.2, 0.2, 0.2))
par(oma=c(4,4,4,4))
for(i in 1:length(nexp)){
    for(j in 1:length(nexp)){
        print(paste(i,j))
        if(i == j){
            plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
            text(x = 0.5, y = 0.4, names(diffList)[i], cex=2)
            text(x = 0.5, y = 0.22, "", cex=1.5)
        }
        else{
            clow = cor(t(nexp[[i]][wval,][wlow,]), t(nexp[[j]][wval,][wlow,]))
            chigh = cor(t(nexp[[i]][wval,][whigh,]), t(nexp[[j]][wval,][whigh,]))
            plot(density(diag(clow), na.rm=T), ann = F, col="red", bty = 'n', lwd=4, ylim=c(0,4.5), cex.lab=1.6, cex.axis=1.6, ylab="", yaxt="n")
            lines(density(diag(chigh), na.rm=T), lwd=4, col="green")
        }
    }
}
dev.off()



wwl = which(rowSums(nhi) > 5 & rowSums(nhi) < 100)
clow = cor(t(nhi[wwl,]), t(nka[wwl,]))

wwh = which(rowSums(nhi) > 1000)
chigh = cor(t(nhi[wwh,]), t(nka[wwh,]))

pdf("correlation_low_high.pdf")
plot(density(diag(clow), na.rm=T), lwd=4, ylim=c(0,4.5), cex.lab=1.6, cex.axis=1.6, ylab="", yaxt="n")
lines(density(diag(chigh), na.rm=T), lwd=4, col="red")


dev.off()


pdf("figures/sample_gene_count.pdf", 9, 7)

par(mar=c(6,6,6,6))
boxplot(sample_count, lwd=2, outline=F, col=c("aquamarine", "lightpink"), ylab="gene count sum", labels="", cex.lab=1.6, cex.axis=1.6, xaxt="n")
axis(1, 1:5, names(sample_count), cex.axis=1.6)
legend("topleft", legend=c("up-regulated", "down-regulated"), fill=c("aquamarine", "lightpink"), cex=1.4, bg="white")
dev.off()

boxplot(updown)
abline(v=(0.5+2*(1:4)))


dbs <- c("GO_Biological_Process_2018", "KEGG_2019_Human", "ChEA_2016", "Human_Phenotype_Ontology")

enrichment = list()

for(j in 1:length(diffList)){
    enrichList = list()
    
    for(i in 1:ncol(diffList[[j]])){
        print(i)
        oo = order(diffList[[j]][,i])
        
        up_genes = rownames(diffList[[j]])[oo[1:250]]
        
        enriched = tryCatch({
            enrichr(up_genes, dbs)
        }, warning = function(w) {
        }, error = function(e) {
            print("try again")
            enrichr(up_genes, dbs)
        }, finally = {
        })
        
        processes = c()
        for(en in enriched){
            ww = which(en[,4] < 0.05)
            processes = c(processes, en[ww,1])
        }
        
        enrichList[[length(enrichList)+1]] = processes
    }
    enrichment[[length(enrichment)+1]] = enrichList
}

overl = matrix(0, 5, 5)
colnames(overl) = names(diffList)
rownames(overl) = names(diffList)

for(i in 1:length(enrichment)){
    for(j in 1:length(enrichment)){
        
        s1 = 0
        
        for(k in 1:length(enrichment[[i]])){
            inter = intersect(enrichment[[i]][[k]], enrichment[[j]][[k]])
            s1 = s1 + length(inter)
        }
        
        overl[i,j] = s1
    }
}

xtable(overl, digits=0)





