### SCRIPT TO CALL THE VARIOUS PLOTTING SCRIPTS TO PRODUCE A SINGLE PLOT ###
############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
main_dir <- "~/repos/" ## where the code is kept
source(paste0(main_dir,"popgen/packages_dev/functionWriter.R"))
setwd(paste0(main_dir,"abu_dhabi/"))
###########################################################
## SWITCH IDS TO FAMILY NAMES AND CLEAN UP ##
snames <- read.table("data/AD_FID_Surnames.txt", header = T, as.is = T,
                     stringsAsFactors = F, comment.char = "",sep = "\t")
ids <- read.table("data/AbuDhabi.inds", header = F, as.is = T, stringsAsFactors = F)

## LET'S TIDY SOME OF THESE UP
snames$Surnames <- gsub("Mohammed","Mohamed",snames$Surnames)
snames$Surnames <- gsub("al ","Al ",snames$Surnames)
snames$Surnames <- gsub("AL ","Al ",snames$Surnames)
snames$Surnames <- gsub("BAl ","Bal ",snames$Surnames)
snames$Surnames <- gsub(" ameri"," Ameri",snames$Surnames)
snames$Surnames <- gsub("Al Balooshi","Al Baloushi",snames$Surnames)
snames$Surnames <- gsub("Al Belooshi","Al Baloushi",snames$Surnames)
snames$Surnames <- gsub("Al Blooshi","Al Baloushi",snames$Surnames)
snames$Surnames <- gsub("Al bloushi","Al Baloushi",snames$Surnames)
snames$Surnames <- gsub("Al Bloushi","Al Baloushi",snames$Surnames)
snames$Surnames <- gsub("Al dhaheri","Al Dhaheri",snames$Surnames)
snames$Surnames <- gsub("Al Ghefli","Al Ghafli",snames$Surnames) 
snames$Surnames <- gsub("Al Ghufli","Al Ghafli",snames$Surnames)
snames$Surnames <- gsub("Al Hajri","Al Hajeri",snames$Surnames)
snames$Surnames <- gsub("Al kindi","Al Kindi",snames$Surnames)
snames$Surnames <- gsub("Al Mansouri","Al Mansoori",snames$Surnames)
snames$Surnames <-  gsub(" ","_",snames$Surnames)
#snames$Surnames <-  gsub("__","_",snames$Surnames)
## FOR ANY SURNAME WITH FEWER THAN TWO INDS REPRESENTING, 
## BLANK OUT TO OTHER
names2keep <- names(table(snames$Surnames))[table(snames$Surnames) > 5]
snames$Surnames[!snames$Surnames%in%names2keep] <- "OTHER"
###########################################################
run <- "A"

#for(run in c("B","C"))
#{
  fsanalyname <- paste0('AbuDhabiPhasedFull',run)
  # leginfo_file <- "data/MalariaGenAdmixturePopulationKey.txt"
  # mapinfo_file <- "data/MalariaGenAdmixturePopulationKey2.txt"
  # latlong_file <- "data/MalariaGenAdmixturePopulationKeyLatLongs.txt"
  # poppos_file <- "data/MalariaGenAdmixturePopulationKeyMapPositions.txt"
  # popkey_file <- "data/MalariaGenAdmixturePopulationOverview.txt"
  # popnums_file <- "data/MalariaGenPopulationKeyCPanalysisPopNumbers.txt"
  # pca_file <- "data/Africa300KPCS.txt"
  mcmc_file <- paste0("data/",fsanalyname,".mcmc.xml")
  tree_file <- paste0("data/",fsanalyname,".tree.xml")
  mat_file <- paste0("data/AbuDhabiPhased.chunkcounts.out")
  mat_file2 <- paste0("data/AbuDhabiPhased.chunklengths.out")
  ## LOAD POPKEY FILE ##
  # popkey <- read.table(popkey_file,header=T,as.is=T)
  # popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
  
  
  
  pdf(paste0("figures/",fsanalyname,"FS_2018.pdf"), height = 9.5, width = 9)
  
      layout(matrix(c(5,2,
                      1,3,
                      5,4),3,2, byrow = T),
             widths = c(1,8), heights = c(1,8,0.5))
  
      ###########################################################
      ## 03 FS TREE
      plot_vert <- TRUE
      plot_rev <- FALSE
      plot_tree_labels <- FALSE
      #par(mar=c(0.5,0.5,0,0))
      par(mar=c(5.5,0.5,0,0))
      source("plotting_scripts/fstree.R")
      ###########################################################
      ## 05 FS TREE_2
      plot_vert <- FALSE
      plot_rev <- TRUE
      par(mar=c(0,0,1,5))
      source("plotting_scripts/fstree.R")
      ##########################################################
      ## 04 FS COPYING MATRIX
      tmatmax <- 150 # cap the heatmap
      tmatmin <- 0
      plot_scale <- TRUE
      # par(mar=c(0.5,0,0,0.5))
      # plot_axis_labels <- FALSE
      par(mar=c(5.5,0,0,5))
      plot_axis_labels <- TRUE
      ## CODE TO SWITCH LABELS ##
      new_labels <- c()
      for(i in tree_labels)
      {
        tmp <- snames$Surnames[snames$FID==i]
        if(length(tmp) == 0) tmp <- snames$Surnames[snames$FID==gsub("a","",i)]
        if(length(tmp) > 1 ) tmp <- tmp[length(tmp)]
        if(is.na(tmp)) tmp <- "NA"
        new_labels <- c(new_labels,tmp)
      }
      
      mapstate <- extractValue(treexml,"Pop")
      mapstatelist <- popAsList(mapstate)
      mapstatelistnames <- mapstatelist
      for(i in 1:length(mapstatelist))
      {
        print(i)
        tmpinds <- c()
        for(j in mapstatelist[[i]])
        {
          tmp <- snames$Surnames[names$ID==j]
          if(length(tmp) == 0) tmp <- snames$Surnames[snames$FID==gsub("a","",j)]
          if(length(tmp) > 1 ) tmp <- tmp[length(tmp)]
          if(is.na(tmp)) tmp <- "NA"
          tmpinds <- c(tmpinds,tmp)
        }  
        mapstatelistnames[[i]] <- tmpinds
      }
      
      popnames <- 1:length(mapstatelist)
      
      clusterNames <- function(cl) { 
        tb <- sort(table(gsub("[0-9]","",unlist(cl))),decreasing=T)
        if("OTHER"%in%names(tb)) tb <- c(tb[!names(tb)%in%"OTHER"],tb[names(tb)%in%"OTHER"])
        nn <- paste0("{",sum(tb),"}", gsub(" ","",paste(names(tb),"[",tb,"]",collapse="/")),sep = " ")
        nn
      }
      popnames <- lapply(mapstatelistnames,clusterNames)
      names(mapstatelist) <- names(mapstatelistnames) <- popnames
      
      axis_labels <- names(mapstatelist)
      axis_tmp <- as.vector(unlist(lapply(mapstatelist,length)))
      #axis_pos <- cumsum(axis_tmp) - apply(cbind(cummin(axis_tmp),cummax(axis_tmp)),1,mean)
      #axis_pos <- nrow(tmpmat) - c(axis_pos[2:length(axis_pos)], nrow(tmpmat))
      axis_pos <- nrow(tmpmat) - (cumsum(axis_tmp) - (axis_tmp/2))
      axis_tmp <- (nrow(tmpmat) - cumsum(axis_tmp))
      
      source("plotting_scripts/copyingmatrix.R")
      ###########################################################
      ## 06 PLOT PANEL LETTER TOP LEFT ABOVE TREE
      par(mar=c(0,0,1,0))    
      plot(0,0,type="n",axes=F,xlab="",ylab="")
      #legend("topleft",legend="D",cex=2,bty="n")
  
  
  dev.off()
    
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
    ##
    popdend<-makemydend(tdend,mapstatelist) # use NameSummary to make popdend
    popdend<-fixMidpointsComplete(popdend) # needed for obscure dendrogram reasons
    # 
    popdendclear<-makemydend(tdend,mapstatelist,TRUE)# use NameMoreSummary to make popdend
    popdendclear<-fixMidpointsComplete(popdendclear) # needed for obscure dendrogram reasons
    
    
    ########
    ## PLAY WITH TREE
    
    library(dendextend)
    library(RColorBrewer)
    
    K <- 8
    tdend %>%
      raise.dendrogram (1) %>%
      set("labels_cex", 0.01)   %>%
      set("branches_k_color",value = brewer.pal(K,"Dark2"), k = K) %>% 
      set("branches_lwd",2) %>% 
      plot(horiz = T)
    
    ## USE THIS TO FIND WHERE THE BREAKS BETWEEN CLUSTERS ARE
    ktree <- cutree(tdend,k = K, order_clusters_as_data = FALSE) 
    kclusts <- kclustnames <- vector("list",K)
    for(i in 1:K) kclusts[[i]] <- names(which(ktree == i))
    
    
    for(i in 1:length(kclusts))
    {
      print(i)
      tmpinds <- c()
      for(j in kclusts[[i]])
      {
        tmp <- snames$Surnames[names$ID==j]
        if(length(tmp) == 0) tmp <- snames$Surnames[snames$FID==gsub("a","",j)]
        if(length(tmp) > 1 ) tmp <- tmp[length(tmp)]
        if(is.na(tmp)) tmp <- "NA"
        tmpinds <- c(tmpinds,tmp)
      }  
      kclustnames[[i]] <- tmpinds
    }
    
    clusterNames <- function(cl) { 
      tb <- sort(table(gsub("[0-9]","",unlist(cl))),decreasing=T)
      if("OTHER"%in%names(tb)) tb <- c(tb[!names(tb)%in%"OTHER"],tb[names(tb)%in%"OTHER"])
      nn <- paste0("{",sum(tb),"}", gsub(" ","",paste(names(tb),"[",tb,"]",collapse="/")),sep = " ")
      nn
    }
    regnames <- lapply(kclustnames,clusterNames)
    names(mapstatelist) <- names(mapstatelistnames) <- popnames
    
    
    
    pdf(paste0("figures/",fsanalyname,"treematrix_2018.pdf"), height = 9.5, width = 11)
      par(xaxs="i", yaxs = "i")
      layout(matrix(c(5,2,8,
                      1,3,6,
                      5,4,7),3,3, byrow = T),
             widths = c(1,8,4), heights = c(1,8,0.5))
      
      ###########################################################
      ## 03 FS TREE
      plot_vert <- TRUE
      plot_rev <- FALSE
      plot_tree_labels <- FALSE
      #par(mar=c(0.5,0.5,0,0))
      par(mar=c(5.5,0.5,0,0))
      K <- 8
      tdend %>%
        raise.dendrogram (1) %>%
        set("labels_cex", 0.01)   %>%
        set("branches_k_color",value = brewer.pal(K,"Dark2"), k = K) %>% 
        set("branches_lwd",2) %>% 
        plot(horiz = plot_vert, axes = F)
      
      ###########################################################
      ## 05 FS TREE_2
      plot_vert <- FALSE
      plot_rev <- TRUE
      par(mar=c(0,0,1,0))
      K <- 8
      tdend %>%
        raise.dendrogram (1) %>%
        set("labels_cex", 0.01)   %>%
        set("branches_k_color",value = brewer.pal(K,"Dark2"), k = K) %>% 
        set("branches_lwd",2) %>% 
        dendextend::rotate(length(labels(tdend)):1) %>%
        plot(horiz = plot_vert, axes = F)
      
      ##########################################################
      ## 04 FS COPYING MATRIX
      tmatmax <- 150 # cap the heatmap
      tmatmin <- 0
      plot_scale <- TRUE
      # par(mar=c(0.5,0,0,0.5))
      # plot_axis_labels <- FALSE
      par(mar=c(5.5,0,0,0))
      plot_axis_labels <- FALSE
      ## CODE TO SWITCH LABELS ##
      new_labels <- c()
      for(i in tree_labels)
      {
        tmp <- snames$Surnames[snames$FID==i]
        if(length(tmp) == 0) tmp <- snames$Surnames[snames$FID==gsub("a","",i)]
        if(length(tmp) > 1 ) tmp <- tmp[length(tmp)]
        if(is.na(tmp)) tmp <- "NA"
        new_labels <- c(new_labels,tmp)
      }
      
      mapstate <- extractValue(treexml,"Pop")
      mapstatelist <- popAsList(mapstate)
      mapstatelistnames <- mapstatelist
      for(i in 1:length(mapstatelist))
      {
        print(i)
        tmpinds <- c()
        for(j in mapstatelist[[i]])
        {
          tmp <- snames$Surnames[names$ID==j]
          if(length(tmp) == 0) tmp <- snames$Surnames[snames$FID==gsub("a","",j)]
          if(length(tmp) > 1 ) tmp <- tmp[length(tmp)]
          if(is.na(tmp)) tmp <- "NA"
          tmpinds <- c(tmpinds,tmp)
        }  
        mapstatelistnames[[i]] <- tmpinds
      }
      
      popnames <- 1:length(mapstatelist)
      
      clusterNames <- function(cl) { 
        tb <- sort(table(gsub("[0-9]","",unlist(cl))),decreasing=T)
        if("OTHER"%in%names(tb)) tb <- c(tb[!names(tb)%in%"OTHER"],tb[names(tb)%in%"OTHER"])
        nn <- paste0("{",sum(tb),"}", gsub(" ","",paste(names(tb),"[",tb,"]",collapse="/")),sep = " ")
        nn
      }
      popnames <- lapply(mapstatelistnames,clusterNames)
      names(mapstatelist) <- names(mapstatelistnames) <- popnames
      
      axis_labels <- names(mapstatelist)
      axis_tmp <- as.vector(unlist(lapply(mapstatelist,length)))
      #axis_pos <- cumsum(axis_tmp) - apply(cbind(cummin(axis_tmp),cummax(axis_tmp)),1,mean)
      #axis_pos <- nrow(tmpmat) - c(axis_pos[2:length(axis_pos)], nrow(tmpmat))
      axis_pos <- nrow(tmpmat) - (cumsum(axis_tmp) - (axis_tmp/2))
      axis_tmp <- (nrow(tmpmat) - cumsum(axis_tmp))
      
      axis_tmp <- 0
      for(i in 1:K) axis_tmp <- c(axis_tmp,max(which(ktree == i)))
      cluster_labs <- nrow(tmpmat) - (axis_tmp[1:K] + diff(axis_tmp)/2)
      source("plotting_scripts/copyingmatrix.R")
      ###########################################################
      ## 06 PLOT PANEL LETTER TOP LEFT ABOVE TREE
      par(mar=c(0,0,1,0))    
      plot(0,0,type="n",axes=F,xlab="",ylab="")
      #legend("topleft",legend="D",cex=2,bty="n")
      
      ## SURNAME GROUPS
      par(mar=c(5.5,0,0,0))  
      numsnames <- length(unique(snames$Surnames))
      newmat <- matrix(0,nr = length(labels(tdend)),nc = numsnames)
      colnames(newmat) <- sort(unique(snames$Surnames))
      newmat2 <- newmat
      for(i in 1:K)
      {
        tmp <- table(kclustnames[[i]])
        for(j in which(ktree == i))  newmat[j, names(tmp)] <- i#tmp
        tmp <- unlist(kclustnames)[which(ktree == i)]
        k <- 1
        for(j in which(ktree == i))
        {
          newmat2[j, tmp[k]] <- i#tmp
          k <- k + 1
        }
      }
      
      hc <- hclust(dist(t(newmat[,colnames(newmat2)!="OTHER"])))
      
      ## PLOT HCLUST OBJECT AT TOP
      #as.dendrogram(hc) %>% dendextend::rotate(30:1) %>% plot(axes = F)
      ktree2 <- cutree(as.dendrogram(hc),k = 7, order_clusters_as_data = FALSE) 
      
      col_order <- c(hc$order,30)
      newmat2 <- newmat2[,col_order]
      
      image(1:ncol(newmat2),
            1:nrow(newmat2),
            t(newmat2),
            col = c("white",brewer.pal(K,"Dark2")), axes = F,xlab = "")
      for(i in 1:max(ktree2))
      {
        tmp <- max(which(rev(ktree2) == i))
        abline(v = tmp+0.5)
      }
      axis(1,at = 1:numsnames,labels = colnames(newmat)[col_order], las = 2,
           lwd = 0)
    dev.off()
    
#}


counts <- read.table(mat_file,header=T,row.names=1)
counts <- counts[tree_labels,rev(tree_labels)]
counts <- as.matrix(counts)

lengths <- read.table(mat_file2,header=T,row.names=1)
lengths <- lengths[tree_labels,rev(tree_labels)]
lengths <- as.matrix(lengths)
