##########################################################
## PLOT FINESTRUCTURE TREE
##########################################################
library(ape)
main_dir <- "~/repos/"
source(paste0(main_dir,"popgen/packages_ext/FinestructureLibrary.R"))
source(paste0(main_dir,"popgen/packages_ext/FinestructureLibrary_GB.R"))

##########################################################
## LOAD DATA AND TREE
#leginfo <-read.table(leginfo_file,header=T,comment.char="")
treexml <- xmlTreeParse(tree_file)
ttree <- extractTree(treexml)
## MANIPULATE TREE ##
ttree$tip.label <- tidyNames(ttree$tip.label)
tree_labels <- ttree$tip.label #getPopSymbols(tree_labels,leginfo)

tree2plot <- ttree
tree2plot <- as.phylo(tree2plot)
###########################################################
## ROTATE SOME OF THE NODES
#tree2plot <- reorder.phylo(rotate(tree2plot,c("JOLA_S1","YORUBA2")))
#tree2plot <- rotate(tree2plot,2211)

if(plot_tree_labels == TRUE)
{
    tree2plot$tip.label <- tree_labels$pop_vec
} else 
{
    tree2plot$tip.label <- rep("",length(tree_labels)) 
}


###########################################################
edge_cols <- c()
for(i in 1:nrow(tree2plot$edge))
{
    # test <- tree_labels$pop_vec[tree2plot$edge[i,]]
    cl <- "#000000"
    # if(sum(!is.na(test))>0) cl <- as.character(leginfo$Colour[leginfo$EthnicGroup==test[!is.na(test)]])
    edge_cols <- c(edge_cols,cl)
}
###########################################################
## PLOT THE TREE
if(plot_vert == TRUE)
{
  ylim <- c(1,length(tree2plot$tip.label))
  if(plot_rev == TRUE) ylim <- rev(ylim)
  plot(tree2plot,lab4ut="axial",type="phylogram",xaxs="i",yaxs="i",
       edge.color=edge_cols,use.edge.length=F, y.lim = ylim) #tip.color=tree_labels$col2plot,
    
    # ###########################################################
    # ## ADD TO THE TREE
    # ## split the tree into major branches
    # for(i in 2:length(rect_splits))
    # {
    #   rect(rect_split_x,rect_splits[(i-1)]-1,2400,rect_splits[i],
    #        col=makeTransparent(rect_cols[(i-1)],alpha=50),border=NA)
    #   rect(-100,rect_splits[(i-1)]-1,rect_split_x,rect_splits[i],
    #        col=makeTransparent(rect_cols[(i-1)],alpha=100),border=NA)
    #   text(550,y=rect_splits[(i)]-((rect_splits[(i)]-rect_splits[(i-1)])/2),
    #        labels=gsub("\\_","\n",rect_labels[(i-1)]),cex=0.5,xpd=T)
    # }
}

if(plot_vert == FALSE)
{
  xlim <- c(1,length(tree2plot$tip.label))
  if(plot_rev == TRUE) xlim <- rev(xlim)
  plot(tree2plot,lab4ut="axial",type="phylogram",xaxs="i",yaxs="i",
       edge.color=edge_cols,use.edge.length=F,direction="downwards",
       x.lim = xlim) #tip.color=tree_labels$col2plot,
}