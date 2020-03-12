# Install packages
# if(!"RCy3" %in% installed.packages()){
#   install.packages("BiocManager")
#   BiocManager::install("RCy3")
# }
pacman::p_load("data.table", "scales", "tidyverse", "stringr", 
               "Biostrings","DECIPHER","igraph","RColorBrewer", "RCy3")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/guanylurea_analysis/")

# Read in the datasets
gu1 <- readAAStringSet("data/RutB_homologs.fasta")
names(gu1) <- paste0(gsub(" ", "_", names(gu1)), "_guanylurease")
gu2 <- readAAStringSet("data/RutB_homologs_genbank.fasta")
names(gu2) <- paste0(gsub(" ", "_", names(gu2)), "_guanylurease")
bh <- readAAStringSet("data/biuret_hydrolases_of_interest_for_clustering.faa")
bh2 <- bh[-grep("atzE", names(bh))]
names(bh2) <- paste0(names(bh2), "_biuretase")

fam <- list.files("data/CDD_by_subfamily/")
sqs<-AAStringSet()
for(i in 1:length(fam)){
  rd <- readAAStringSet(paste0("data/CDD_by_subfamily/", fam[i]))
  names(rd) <- paste0(names(rd), "_", gsub("\\.fasta", "", fam[i]))
  sqs <- c(rd,sqs)
}
names(sqs)

# Read in the PDB ids
pdb <- readAAStringSet("data/35_unique_CH_PDB.fa")

# Combine everything
comb <- c(pdb, gu1, gu2, bh2, sqs)
dedup <- comb[!duplicated(comb)]
names(dedup) <- gsub(" |\\[|\\]|\\(|\\)|\\/|\\=|:|\\.", "_", names(dedup))
names(dedup)
writeXStringSet(dedup, "data/151_isochorismatases_unaligned.fasta")
al <- AlignSeqs(dedup)
BrowseSeqs(al)

# Run an all-v-all BLAST
# makeblastdb -in 151_isochorismatases_unaligned.fasta -dbtype prot -out db/151_prot_blast_db
# blastp -db db/151_prot_blast_db -query 151_isochorismatases_unaligned.fasta -outfmt 6 -out 151_isochorismatase_superfamily_all_v_all.tsv

# Read in the all-vs-all BLAST table
allvall <- fread("data/151_isochorismatase_superfamily_all_v_all.tsv", data.table = F)
shrt <- allvall[, c(1:2,11)]
colnames(shrt) <- c("prot1","prot2","eval")
noprs <- shrt[!shrt$prot1==shrt$prot2,] # Remove all pairs
noprs <- noprs[order(noprs$eval),]
head(noprs)

###  Make cluster diagram

# Color by PFAM class
namvec <- unique(stringr::word(names(dedup), -1, sep = "_"))
namvec <- gsub("like", "YcaC_like", namvec)
namvec <- rev(c("pdb",  "guanylurease",  "biuretase", "CSHase", "nicotinamidase", "YcaC_like", "isochorismatase", "other"))

# Color
# pal <- c("gray", "#E41A1C", "#92D050", "#377EB8", "#984EA3", "#FF7F00",   
#          "goldenrod", "#A65628", "#F781BF", "blue1", 
#          "gray68", "darkorchid1", "navy", "plum1",
#          "deepskyblue", "gold", "deeppink2", "lightslateblue",
#          "lightblue2", "darkseagreen1")
pal <- rev(c("gray", "red", "turquoise", "mediumorchid1","#F2EA91", "#3288BD", "#FA9D59", "purple4"))




for (i in 101:120) {
  source("lib/make_cluster_diagram.r")
  shrt <- allvall[,(c(1:2,11))]
  colnames(shrt) <- c("prot1", "prot2", "eval")
  
  l <- layout_components(g)
  thresh <- i
  gr <- make_cluster(shrt, thresh = thresh, bynum = 1, namvec, pal)
  
  net <- gr[[1]]
  g <- gr[[2]]
  
  pdf(file = paste0("output/clustering_diagram_eval_10e-", thresh, ".pdf"), width = 8, height = 8)
  par(mar=c(1,1,1,1))
  
  pl<-plot(g, vertex.label = NA,
           vertex.label.cex= 0.25,
           vertex.size = 3,
           layout=l, edge.color = "gray60", edge.width=0.3)
  title(main = paste0("BLAST e-value cut-off: 1e-", thresh))
  dev.off()
}


thresh <- 103
gr <- make_cluster(shrt, thresh = thresh, bynum = 1, namvec, pal)

net <- gr[[1]]
g <- gr[[2]]
l <- layout_components(g)
pdf(file = paste0("output/clustering_diagram_eval_10e-", thresh, "_labeled.pdf"), width = 8, height = 8)
par(mar=c(1,1,1,1))

pl<-plot(g, vertex.label = ifelse(grepl("guanylurease", V(g)$name), V(g)$name, NA),
         vertex.label.cex= 0.25,
         vertex.size = 3,
         layout=l, edge.color = "gray60", edge.width=0.3)
title(main = paste0("BLAST e-value cut-off: 1e-", thresh))
dev.off()

 # Create legend
cleg <- gr[[3]]
namvec
# pal2 <- pal[1:length(namvec)]

pdf(file=paste0("output/clustering_legend.pdf"), width = 5, height = 5)
plot.new()
legend("center", legend = namvec, fill = pal, 
       bty="n", ncol = 1)
dev.off() 

net <- createNetworkFromIgraph(g)# "output/isochorismatase_like_hydrolases_for_cytoscape")
