setwd("~/OneDrive/alignmentbench")

signatures = read.table("supportfiles/signatures.tsv", stringsAsFactors=F, sep="\t", fill=T, quote="", header=T)
sampleMap = read.table("supportfiles/samplemapping.tsv", stringsAsFactors=F, sep="\t", fill=T, quote="", header=T)

whuman = which(signatures[,14] == "homo sapiens")
wmouse = which(signatures[,14] == "Mus musculus")

controlList = list()
perturbationList = list()

for(co in whuman){
    controlList[[length(controlList)+1]] = sort(unlist(c(strsplit(signatures[co, 6], ","), strsplit(signatures[co,10], ","))))
    perturbationList[[length(perturbationList)+1]] = sort(unlist(c(strsplit(signatures[co,7], ","), strsplit(signatures[co,9], ","))))
}

ww = which((lapply(perturbationList,length) == 3) & (lapply(controlList,length)  == 3))

samples = sort(unique(unlist(c(controlList[ww], perturbationList[ww]))))
sraids = sort(unique(sampleMap[sampleMap[,3] %in% samples, 4]))
writeLines(sraids, con="supportfiles/sra_files_human.txt")


srasplit = split(sraids, cut(1:length(sraids), breaks=10, include.lowest=TRUE))

counter = 0
for(s in srasplit){
    counter = counter+1
    writeLines(s, con=paste("supportfiles/split/sra_files_human_",counter,".txt"))
}


controlList = list()
perturbationList = list()

for(co in wmouse){
    controlList[[length(controlList)+1]] = sort(unlist(c(strsplit(signatures[co, 6], ","), strsplit(signatures[co,10], ","))))
    perturbationList[[length(perturbationList)+1]] = sort(unlist(c(strsplit(signatures[co,7], ","), strsplit(signatures[co,9], ","))))
}

ww = which((lapply(perturbationList,length) == 3) & (lapply(controlList,length) == 3))

samples = sort(unique(unlist(c(controlList[ww], perturbationList[ww]))))
sraids = sort(unique(sampleMap[sampleMap[,3] %in% samples, 4]))
writeLines(sraids, con="supportfiles/sra_files_mouse.txt")


