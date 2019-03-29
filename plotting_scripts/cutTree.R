####function that cuts a tree of XX clusters, where the number of clusters equal to
# exclude.n (e.g.= 1) is < then fract(e.g 0.05)
# It extracts the heights from the popdendclear tree
cutTree=function(popdendclear,tdend,fract=0.05,exclude.n=1){
  
  perc.sing=1  #number of singleton
  hh=1
  
  while (perc.sing > fract){
    tmp.tree=cut.dendrogram(tdend,h=h.list[hh])
    print(paste("Cutting number",hh))
    
    nK=length(unlist(tmp.tree$upper))
    
    all.members=vector(length=nK,mode="numeric")
    for (x in 1:nK){
      n.members=length(unlist(dendrapply(tmp.tree$lower[[x]],function(f) attr(f,"members"))))
      all.members[x]=n.members
      sing.number=length(all.members[all.members<=exclude.n])
    }
    perc.sing=sing.number/nK
    print (paste("Number of unwanted clusters=",perc.sing))
    hh=hh+1
  }
  ####From George####
  clustergroups <- list()
  for(i in 1:nK){
    clustergroups[[i]] <- labels(tmp.tree$lower[[i]])
  }
  clusternames <- unlist(lapply(clustergroups,clusterNames))
  ######
  my.cutted.dend=makemydend(tdend,lablist=clustergroups,TRUE)
  my.cutted.dend=fixMidpointsComplete(my.cutted.dend)
  return(my.cutted.dend)
}

################

