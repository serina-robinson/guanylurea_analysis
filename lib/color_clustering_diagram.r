color_nodes <- function(noprs, namvec, pal) {
  
  meta<-data.frame(noprs$prot1,stringsAsFactors=F) 
  meta$fam<-rep("unknown",nrow(noprs))
  meta$color<-rep("gray",nrow(noprs))
  meta$size<-rep(2,nrow(noprs))
  meta$shape <-rep("circle", nrow(noprs))
  
  for(i in 1:length(namvec)){
      ind <- grep(paste0("_", namvec[i]), meta$noprs.prot1)
      meta$fam[ind] <- namvec[i]
      meta$size[ind] <- 2
      meta$size[ind] <- ifelse(grepl("_guanylurease", meta$noprs.prot1[ind]), 2, 3)
      meta$color[ind] <- ifelse(grepl("_pdb", meta$noprs.prot1[ind]), "gray", pal[i])
      meta$shape[ind]<- ifelse(grepl("_pdb", meta$noprs.prot1[ind]), "square", "circle")
      
    }
  
  return(meta)
}