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

weights <- c()
for(i in colnames(copyvectors)){
  weights <- c(weights,sum(ids$pop==i))
}

for(i in 1:nrow(copyvectors)){
  tmp <- copyvectors[i,]/weights
  tmp <- tmp/sum(tmp)
  copyvectors[i,] <- tmp
}

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

## now compute nnls for each AD pop - byt only if it doesn't exist...
rrun <- 1
if(rrun == 1){

  nnlsmat <- matrix(0,nc = length(donpops),nr = nrow(copyvectors))
  colnames(nnlsmat) <- sort(donpops)
  rownames(nnlsmat) <- rownames(copyvectors)
  
  for(i in rownames(copyvectors)){
    popvector <- copyvectors[i,colnames(donmat)]
    popvector <- popvector/sum(popvector)
    predmat <- donmat
    ## this is the function
    ourmix <- getoverallfit(predmat,popvector)$x
    nnlsmat[i,names(ourmix)] <- ourmix
  }
  
  ## remove donors that are not involved
  nnlsmat <- nnlsmat[,!colSums(nnlsmat)==0]
}

popord <- read.table(file="final_analysis/indpops_order.txt",
                     header=F)
colnames(popord) <- c("pop","reg")
#popord <- popord[popord$pop!="Egypt",]
colrs <- brewer.pal(length(unique(popord$reg)), "Accent")
popscol <- colrs[popord$reg]



plotmat <- matrix(0,nc=ncol(nnlsmat),nr=length(labels(tdend)))
colnames(plotmat) <- colnames(nnlsmat)
## this puts individuals in same order as tree
# rownames(plotmat) <- labels(tdend)
## we want individuals from emirate together
reg_order <- c("Abu_Dhabi","Dubai","Sharjah","Ajman",
               "Umm_al_Quwain", "Ras_al_Khaimah","Fujairah","Unknown")
ind_order <- ids[ids$ind%in%labels(tdend),]
ind_order[,4] <- factor(gsub("_Emirati","",gsub("_[0-9]","",gsub("_[0-9][0-9]","",ind_order$pop))),levels = reg_order)

ind_order <- ind_order[order(ind_order$V4),]


pdf("figures/NNLSbyEmirate.pdf")

for(i in reg_order){
  inds_in_ord <- ind_order$ind[ind_order$V4==i]
  plotmat <- nnlsmat[rownames(nnlsmat)%in%inds_in_ord,]
  plotmat <- plotmat[,popord$pop]/rowSums(plotmat[,popord$pop])
  
  ## reverse order to match the tree??
  #plotmat <- plotmat[rev(rownames(plotmat)),]
  
  ### REORDER BY LARGEST FRACTION WITHIN EACH EMIRATE
  orderby <- which.max(colMeans(plotmat))
  plotmat <- plotmat[order(plotmat[,orderby]),]
  
  par(mar=c(4,4,2,2))
  barplot(t(plotmat),col = popscol,
          horiz = T, las = 2, border = NA, space = 0,
          axes = F, names.arg = rep("",nrow(plotmat)),
          main=i)

}

# ## 08 PLOT NNLS LEGEND
par(mar=c(0,0,0,0))
plot(0,0,type = "n", xlab = "", ylab = "", axes = F)
legend("top",legend = levels(popord$reg),
       fill =  colrs, ncol = 2, bty = "n")

dev.off()

# 
