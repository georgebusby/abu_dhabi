### SCRIPT TO CALL THE VARIOUS PLOTTING SCRIPTS TO PRODUCE A SINGLE PLOT ###
############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
main_dir <- "~/repos/" ## where the code is kept
source("popgen/packages_dev/functionWriter.R")
setwd(paste0(main_dir,"abu_dhabi/"))
library(ape)
source("popgen/packages_ext/FinestructureLibrary.R")
source("popgen/packages_ext/FinestructureLibrary_GB.R")


###########################################################
## SWITCH IDS TO FAMILY NAMES AND CLEAN UP ##
# ids <- read.table("final_analysis/ind_id_pop_id_include_0or1_KE_GB_comb.txt",
#                   header = T, as.is = T, stringsAsFactors = F, sep ="\t")


###########################################################

fsanalyname <- paste0('710K_Emirati_painted_combined_fs_Emirati')
mcmc_file <- paste0("final_analysis/",fsanalyname,".mcmc.xml")
tree_file <- paste0("final_analysis/",fsanalyname,".tree.xml")
mat_file <- paste0("final_analysis/710K_Emirati_painted_combined_all.chunkcounts.out")
## read in the mcmc files
mcmcxml <- xmlTreeParse(mcmc_file)
## convert this into a data frame
mcmcdata <- as.data.frame.myres(mcmcxml) 
## read the tree as xml format
treexml <- xmlTreeParse(tree_file)
## extract the tree into ape's phylo format
ttree <- extractTree(treexml)
## now is a good time to remove them via:
ttree$node.label<-NULL
## Will will instead remove "perfect" node labels
#ttree$node.label[ttree$node.label=="1"] <-""
## And reduce the amount of significant digits printed:
#ttree$node.label[ttree$node.label!=""] <-format(as.numeric(ttree$node.label[ttree$node.label!=""]),digits=2)
    
## convert to dendrogram format
tdend <- myapetodend(ttree,factor=1, tol=1) 
#tdend2 <- sort(tdend)
##
popdend<-makemydend(tdend,mapstatelist) # use NameSummary to make popdend
popdend<-fixMidpointsComplete(popdend) # needed for obscure dendrogram reasons
# 
popdendclear<-makemydend(tdend,mapstatelist,TRUE)# use NameMoreSummary to make popdend
popdendclear<-fixMidpointsComplete(popdendclear) # needed for obscure dendrogram reasons


ids <- read.table("final_analysis/ind_id_pop_id_include_0or1_KE_GB_comb_corrected3.txt",
                  header = T, as.is = T, stringsAsFactors = F, sep ="\t", comment.char = "$")

slabels <- ids$KE_correct
snames <- slabels[match(labels(tdend),ids$IIDx)]
snames[is.na(snames)] <- "Unknown"

slabels <- paste(ids$KE_correct,ids$name_code,sep="_")
snames <- slabels[match(labels(tdend),ids$IIDx)]
snames[is.na(snames)] <- "Unknown"



## sort out ordering within clusters
mapstate <- extractValue(treexml,"Pop")
mapstatelist <- popAsList(mapstate)
mapstatelistnames <- mapstatelist
for(i in 1:length(mapstatelist)){
  print(i)
  tmpinds <- c()
  for(j in mapstatelist[[i]]){
    tmp <- snames[ids$IIDx==j]
    tmpinds <- c(tmpinds,tmp)
  }  
  mapstatelistnames[[i]] <- tmpinds[order(tmpinds)]
  mapstatelist[[i]] <- mapstatelist[[i]][order(tmpinds)]
}



########
## PLAY WITH TREE
library(dendextend)
library(RColorBrewer)

K <- 14
tdend %>%
  raise.dendrogram (1) %>%
  set("labels_cex", 0.01)   %>%
  set("branches_k_color", value = brewer.pal(K,"Set3"), k = K) %>% 
  set("branches_lwd",2) %>% 
  plot(horiz = T)
    
 ## USE THIS TO FIND WHERE THE BREAKS BETWEEN CLUSTERS ARE
ktree <- cutree(tdend,k = K, order_clusters_as_data = FALSE) 
kclusts <- kclustnames <- vector("list",K)
for(i in 1:K) kclusts[[i]] <- names(which(ktree == i))
for(i in 1:length(kclusts)){
  print(i)
  tmpinds <- c()
  for(j in kclusts[[i]]){
    tmp <- slabels[ids$IIDx==j]
    if(length(tmp) == 0) tmp <- "Unknown"
    tmpinds <- c(tmpinds,tmp)
  }  
  kclustnames[[i]] <- tmpinds
}
    
names(mapstatelist) <- names(mapstatelistnames) <- popnames
    

tree_cols <- "Dark2"
make_palette <- colorRampPalette(brewer.pal(8,tree_cols))
snames_cols <- make_palette(K)

###########################################################
###########################################################
###########################################################
## PLOT HEATMAP / TREE / SURNAMES
pdf(paste0("figures/",fsanalyname,"TreeMatrixEmirates.pdf"), height = 9.5, width = 14)
  par(xaxs="i", yaxs = "i")
  layout(matrix(c(5,2,10,11,
                  1,3,6,7,
                  5,4,9,8),3,4, byrow = T),
         widths = c(1,8,5,0), heights = c(1,8,0.5))
      
###########################################################
## 01 FS TREE
  plot_vert <- TRUE
  plot_rev <- FALSE
  plot_tree_labels <- FALSE
  #par(mar=c(0.5,0.5,0,0))
  par(mar=c(5.5,0.5,0,0))
  tdend %>%
    raise.dendrogram (1) %>%
    set("labels_cex", 0.01)   %>%
    set("branches_k_color",value = snames_cols, k = K) %>% 
    set("branches_lwd",2) %>% 
    plot(horiz = plot_vert, axes = F)
      
###########################################################
## 02 FS TREE_2
  plot_vert <- FALSE
  plot_rev <- TRUE
  par(mar=c(0,0,1,0))
  tdend %>%
    raise.dendrogram (1) %>%
    set("labels_cex", 0.01)   %>%
    set("branches_k_color",value = snames_cols, k = K) %>% 
    set("branches_lwd",2) %>% 
    dendextend::rotate(length(labels(tdend)):1) %>%
    plot(horiz = plot_vert, axes = F)
      
##########################################################
## 03/04 FS COPYING MATRIX
  tmatmax <- 100 # cap the heatmap
  tmatmin <- 0
  plot_scale <- TRUE
  par(mar=c(5.5,0,0,0))
  plot_axis_labels <- FALSE
  ## CODE TO SWITCH LABELS ##
  new_labels <- c()
  for(i in tree_labels){
    tmp <- snames[ids$IIDx==j]
    new_labels <- c(new_labels,tmp)
  }
      
      
  popnames <- 1:length(mapstatelist)
  axis_labels <- names(mapstatelist)
  axis_tmp <- as.vector(unlist(lapply(mapstatelist,length)))
  axis_pos <- nrow(tmpmat) - (cumsum(axis_tmp) - (axis_tmp/2))
  axis_tmp <- (nrow(tmpmat) - cumsum(axis_tmp))
  axis_tmp <- 0
  for(i in 1:K) axis_tmp <- c(axis_tmp,max(which(ktree == i)))
  cluster_labs <- nrow(tmpmat) - (axis_tmp[1:K] + diff(axis_tmp)/2)
  source("plotting_scripts/copyingmatrix.R")
###########################################################
## 05 PLOT EMPTY PANEL TOP LEFT ABOVE TREE
  par(mar=c(0,0,1,0))    
  plot(0,0,type="n",axes=F,xlab="",ylab="")

###########################################################
## 06 SURNAME/EMIRATI GROUPS
  par(mar=c(5.5,0,0,0))  
  numsnames <- length(unique(snames))
  newmat <- matrix(0,nr = length(labels(tdend)),nc = numsnames)
  colnames(newmat) <- sort(unique(snames))
  newmat2 <- newmat
  for(i in 1:K){
    tmp <- table(kclustnames[[i]])
    for(j in which(ktree == i))  newmat[j, names(tmp)] <- i
    tmp <- unlist(kclustnames)[which(ktree == i)]
    k <- 1
    for(j in which(ktree == i)){
      newmat2[j, tmp[k]] <- i 
      k <- k + 1
    }
  }
      
  hc <- hclust(dist(t(newmat[,colnames(newmat2)!="OTHER"])))
## PLOT HCLUST OBJECT AT TOP
  ktree2 <- cutree(as.dendrogram(hc),k = numsnames, order_clusters_as_data = FALSE) 
  col_order <- c(hc$order)
  col_order <- order(colSums(newmat2), decreasing = F)
  
  ## geographically order
  col_order <- c("Abu Dhabi","Dubai","Sharjah","Ajman",
                 "Umm Al Quwain", "Ras Al Khaimah","Fujairah","Unknown")
  
  ## for surnames by emirate
  # nameemir <- unique(ids[,c("KE_correct","name_code")])
  # nameemir[,3] <- paste(nameemir$KE_correct,nameemir$name_code,sep="_")
  # nameemir$name_code[nameemir$name_code == "Name_"] <- "Unknown"
  # nameemir$KE_correct <- factor(nameemir$KE_correct,levels = col_order)
  # nameemir <- nameemir[order(nameemir$KE_correct),]
  
  # col_order <- nameemir$name_code[which(colnames(newmat2)%in%nameemir$name_code)]
  # newmat2 <- newmat2[,col_order]
  
  image(1:ncol(newmat2),
        1:nrow(newmat2),
        t(newmat2),
        col = c("white",snames_cols), axes = F,xlab = "")
  ## add horizontal lines
  abline(h=axis_tmp, col = "grey",lwd = 0.5)
  
  ## 
  # for(i in 1:max(ktree2)){
  #   tmp <- max(which(rev(ktree2) == i))
  #   abline(v = tmp+0.5)
  # }
  
  ## for emirates only
  # axis(1,at = 1:numsnames,
  #      labels = gsub("_Emirati","",col_order), las = 2,
  #      lwd = 0)
  
  sname_labs <- as.character(sapply(colnames(newmat2),function(x){strsplit(x,split="_")[[1]][2]}))
  emir_labs <- as.character(sapply(colnames(newmat2),function(x){strsplit(x,split="_")[[1]][1]}))
  
  
  axis(1,at = 1:ncol(newmat2),
        labels = sname_labs, las = 2,
        lwd = 0)
  
  emir_order <- emir_labs
  for(i in 2:length(emir_order)){
    if(emir_order[(i-1)]!=emir_order[i]){
      abline(v=i-0.5, xpd=T)
    }
  }
  for(i in unique(emir_order)){
    text(x=median(which(emir_order==i)),y=nrow(newmat2),
         label = i, xpd = T, srt = 20, adj = c(1,1))
  }
  
  
  # ###########################################################
  # ## 07/08 PLOT NNLS BAROPLOT FROM ADnnlsByInd.R and legend
  # par(mar=c(5.5,0,0,0))
  # 
  # source("plotting_scripts/ADnnlsByInd.R")
  # 
  # ###########################################################
 
  
dev.off()
    
