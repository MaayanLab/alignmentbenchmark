library("plyr")
library("stringr")


getMapping <- function(species, mapping){
    file = paste0("reference/",species,"_cdna_96.fa")
    ll = readLines(file)
    ll = ll[grep("^>", ll)]
    
    gsym = str_extract(ll, "gene_symbol:.+?\\s|gene_symbol:.+?$")
    gsym = gsub("gene_symbol:|\\s", "", gsym)
    
    gene = str_extract(ll, "gene:.+?\\s|gene:.+?$")
    gene = gsub("gene:|\\s|\\.[0-9]\\s$", "", gene)
    
    transcript = str_extract(ll, "^>.+?\\s")
    transcript = gsub("^>|\\s|\\.[0-9]\\s$", "", transcript)
    
    return(do.call(cbind, list(transcript, gsym, gene)))
}

human_map = getMapping("human")
mouse_map = getMapping("mouse")

save(human_map, mouse_map, file="supportfiles/gene_mapping.rda")

