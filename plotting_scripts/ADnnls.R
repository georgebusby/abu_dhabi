################################################################
## SCRIPT TO RUN THE NNLS PROCEDURE ON THE ABU DHABI CLUSTERS ##
################################################################

## list of IDs
ids <- read.table("final_analysis/ind_id_pop_id_include_0or1_GB.txt",
                  header = F, as.is = T, stringsAsFactors = F, sep ="\t")
colnames(ids) <- c("ind","pop","include")
copyvectors <- read.table("final_analysis/GB_combined.chunklengths.out", 
                          header = T, row.names = 1)

## list of IDs
ids <- read.table("final_analysis/ind_id_pop_id_include_0or1_KE.txt",
                  header = F, as.is = T, stringsAsFactors = F, sep ="\t")
colnames(ids) <- c("ind","pop","include")
copyvectors <- read.table("final_analysis/KE_combined.chunklengths.out", 
                          header = T, row.names = 1)


## get pop designation for each inds in copyvectors
indpops <- ids$pop[match(rownames(copyvectors),ids$ind)]
indpops[indpops == "Bosnia-Herzegovina"] <- "Bosnia.Herzegovina"
indpops[grep("Abu_Dhabi",indpops)] <- "Abu_Dhabi_Emirati"
indpops[grep("Ras_al_Khaimah",indpops)] <- "Ras_al_Khaimah_Emirati"

## now make a matrix of copying vectors for all non Abu Dhabi pops
## let's call these donors
adpops <- c(grep("Abu_Dhabi",unique(indpops), value = T),
            grep("Ras_al",unique(indpops),value = T),
            grep("Emirati",unique(indpops),value = T))
adpops <- sort(unique(adpops))
donpops <- unique(indpops)[!unique(indpops)%in%adpops]
donmat <- matrix(0,ncol=length(donpops),nrow=length(donpops))
colnames(donmat) <- rownames(donmat) <- sort(donpops)
for(i in rownames(donmat)){
  donmat[i,colnames(copyvectors)] <- colSums(copyvectors[indpops==i,])
}

## now normalise this matrix
donmat <- donmat/rowSums(donmat)


## NOW FOR THE MAGIC
## we want to fit each AD pop as a mixture of all these donors
## using Simon Myer's approach
source("~/repos/plasmodium/nnlsfunctions.R")

## now compute nnls for each AD pop
nnlsmat <- matrix(0,nc = length(donpops),nr = length(adpops))
colnames(nnlsmat) <- sort(donpops)
rownames(nnlsmat) <- sort(adpops)

for(i in adpops){
  popvector <- colSums(copyvectors[indpops==i,colnames(donmat)])
  popvector <- popvector/sum(popvector)
  predmat <- donmat
  ## this is the function
  ourmix <- getoverallfit(predmat,popvector)$x
  nnlsmat[i,names(ourmix)] <- ourmix
}

## remove donors that are not involved
nnlsmat <- nnlsmat[,!colSums(nnlsmat)==0]


popord <- scan("final_analysis/poporder.txt", what = "char")
popord <- scan("final_analysis/poporder2.txt", what = "char")

pdf("final_analysis/ADnnlsKE.pdf",height = 20)
layout(c(1:8))
par(mar = c(8,4,4,1))
for(i in rownames(nnlsmat)){
  barplot((nnlsmat[i,popord]), 
          horiz = F, las = 2,
          ylim=c(0,0.25),
          main = paste0(i," n=(",sum(indpops==i),")"))
}
dev.off()






