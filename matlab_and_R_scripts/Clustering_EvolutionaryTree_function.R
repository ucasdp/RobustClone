rm(list=ls())

require(igraph)
require(gmodels)
require(sva)
require(RANN)
require(reshape)
library(shape)
library(vegan)


############################## Louvain-Jaccard clustering ##############################
# Input:
#  AA: the genotype matrix recovered by RPCA.
# Output: 
#  robust_clone: a list variable, where each component represents a cluster with labels of the cells contained.

LJClustering <- function(AA){
  m <- dim(AA)[1]
  if(m<100){
    nn <- 10
  }
  if(m<=100){
    nn <- 10
  }
  if(m>100&m<=1000){
    nn <- 30
  }
  if(m>1000){
    nn <- 150
  }
  
  nearest <- nn2(AA,AA,k=nn+1, searchtype="priority")
  nearest$nn.idx <- nearest$nn.idx[,-1]
  nearest$nn.dists <- nearest$nn.dists[,-1] #Convert to a similarity score
  nearest$nn.sim <-  1*(nearest$nn.dists >= 0)
  
  edges <- melt(t(nearest$nn.idx))
  colnames(edges) <-  c("B", "A", "C")
  edges <-  edges[,c("A","B","C")]
  edges$B <- edges$C; 
  edges$C <- 1
  
  #Remove repetitions
  edges <-  unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))
  
  Adj <-  matrix(0, nrow=nrow(AA), ncol=nrow(AA))
  rownames(Adj) <-  rownames(AA)
  colnames(Adj) <-  rownames(AA)
  Adj[cbind(edges$A,edges$B)] <-  edges$C
  Adj[cbind(edges$B,edges$A)] <-  edges$C
  for(i in 1:nrow(Adj)){
    Adj[i,i] <- 0
  }
  
  g <- graph.adjacency(Adj, mode = "undirected")
  graph_out <- cluster_louvain(g)
  
  clust_assign <- factor(graph_out$membership, levels=sort(unique(graph_out$membership)))
  names(clust_assign) <- graph_out$names
  k <- order(table(clust_assign), decreasing = TRUE)
  new_levels <- rep(1,length(unique(graph_out$membership)))
  new_levels[k] <- 1:length(unique(graph_out$membership))
  levels(clust_assign) <- new_levels
  clust_assign <- factor(clust_assign, levels=1:length(unique(graph_out$membership)))
  robust_clone <- list()
  for(i in 1:length(unique(clust_assign))){
    robust_clone <- c(robust_clone,list(which(clust_assign==i)))
  }
  return(robust_clone)
}


############################## subclonal GTM ###########################################
# Input:
#  robust_clone: the clustering result output by LJClustering function;
#  type: the data type ('SNV' or 'CNV') of input.
# Output: 
#  clone_gety: the inferred clonal genotype based on the Louvain-Jaccard clustering.

subclone_GTM <- function(robust_clone, type){
  clone_gety <- matrix(0,length(robust_clone),ncol(AA))

  #SNV data
  if(type=='SNV'){
    unexit <- list()
    for(i in 1:nrow(clone_gety)){
      clone_cells_gety <- AA[robust_clone[[i]],]
      for(j in 1:ncol(clone_gety)){
        aa <- mode_num(clone_cells_gety[,j])
        if(length(aa)==1){
          clone_gety[i,j] <- mode_num(clone_cells_gety[,j])
        }
        if(length(aa)>1){
          unexit <- c(unexit,list(j))
        }
      }
    }
    
    if(length(unexit)>0){
      clone_gety <- clone_gety[,-unlist(unexit)]
    }
  }
  
  #CNV data
  if(type=='CNV'){
    for(i in 1:nrow(clone_gety)){
      clone_cells_gety <- AA[robust_clone[[i]],]
      for(j in 1:ncol(clone_gety)){
        clone_gety[i,j] <- round(median(clone_cells_gety[,j]))
      }
    }
  }
  return(clone_gety)
}


mode_num <- function(x) #the mutation state with highest frequency
{
  return(as.numeric(names(table(x))[table(x) == max(table(x))]))
}


############################## clonal minimum spanning tree ############################
# Input:
#  clone_gety: the inferred clonal genotype based on the Louvain-Jaccard clustering output by subclone_GTM function;
#  robust_clone: the clustering result output by LJClustering function;
#  pdf_name: the name of pdf with the clonal MST graph.
# Output: 
#  the pdf with the clonal MST graph;
#  el: the connected edges in the MST.

plot_MST <- function(clone_gety, robust_clone, pdf_name){
  
  clone_gety_root <- matrix(2,(length(robust_clone)+1),ncol(clone_gety))
  clone_gety_root[2:(length(robust_clone)+1),] <- clone_gety
  
  rownames(clone_gety_root) <- c('Root',paste('subclone',c(1:length(robust_clone)),sep=''))
  rownames(clone_gety) <- paste('subclone',c(1:length(robust_clone)),sep='')
  
  dis1<- dist(clone_gety_root,p=2)
  dis1 <- as.matrix(dis1)
  dis<- dist(clone_gety,p=2)
  
  spanningtree <- spantree(dis)
  rootid <- which.min(dis1[2:nrow(dis1),1])
  
  idspantree<-spanningtree$kid
  numedge<- nrow(clone_gety)-1
  spantree<-matrix(1,numedge,2)
  spantree[,1]<-as.matrix(2:nrow(clone_gety),nrow(clone_gety),1)
  spantree[,2]<-t(idspantree)
  
  adj<- matrix(0,nrow(clone_gety),nrow(clone_gety))
  adj[spantree]<-1
  adj<-(adj+t(adj))
  
  allpath<-findpath(rootid,rootid,info = adj)
  k1 <- 0
  for(i in 1:length(allpath)){
    k1 <- k1+(length(allpath[[i]])-1)
  }
  el <- matrix(0,k1,2)
  k2 <- 1
  for(i in 1:length(allpath)){
    subpath <- allpath[[i]]
    for(j in 1:(length(subpath)-1)){
      el[k2,1] <- subpath[j]
      el[k2,2] <- subpath[j+1]
      k2 <- k2+1
    }
  }
  el <- unique(el)
  colnames(el) <- c('from','to')
  minspantree <- igraph::graph.edgelist(el, directed = TRUE)
  labelname <- rownames(clone_gety)
  col <- intpalette(c('magenta','red','limegreen','skyblue','pink','brown','gold','blue','cyan'),length(robust_clone))
  size <- c(1:length(robust_clone))
  for(i in 1:length(size)){
    size[i] <- length(robust_clone[[i]])/length(clust_assign)*250
  }
  
  pdf(paste('MST_', pdf_name, '.pdf', sep=''))
  plot(minspantree,layout=layout_as_tree,vertex.color=col,edge.color='black',vertex.label.color='black',alpha=0.5,
       edge.width=4,vertex.size=size,vertex.shape= "sphere",vertex.label=labelname)
  dev.off()
  
  return(el)
}


findpath<-function(node,prev_node,info){
  adj<-which(info[node,]!=0)
  kid<-setdiff(adj,prev_node)
  if(length(kid)==0)
    return(list(node))
  t<-list()
  l<-1
  for (i in kid) {
    path <- findpath(i, node, info)
    for (j in 1:length(path)) {
      t[[l]]<-c(node,path[[j]])
      l<-l+1
    }
  }
  return(t)
}


############################## newly mutated genotypes of each subclone ################
# Input:
#  clone_gety: the inferred clonal genotype based on the Louvain-Jaccard clustering output by subclone_GTM function;
#  robust_clone: the clustering result output by LJClustering function;
#  el: the connected edges in the MST output by plot_MST function.
# Output: 
#  clones_mt_change: a list variable where each component contains the newly mutated genotypes of each subclone.

new_mutation <- function(clone_gety, robust_clone, el){
  clones_mutation <- list()
  for(i in 1:length(robust_clone)){
    clone_mt <- which(clone_gety[i,]!=0)
    clones_mutation <- c(clones_mutation,list(clone_mt))
  }
  
  clones_mt_change <- list() #newly mutated genotypes of each subclone
  root_mt <- clones_mutation[[el[1,1]]]
  clones_mt_change <- c(clones_mt_change,list(root_mt))
  for(i in 1:nrow(el)){
    clone1 <- el[i,1]
    clone2 <- el[i,2]
    clone1_mt <- clones_mutation[[clone1]]
    clone2_mt <- clones_mutation[[clone2]]
    new_mt <- setdiff(clone2_mt,clone1_mt)
    clones_mt_change <- c(clones_mt_change,list(new_mt))
  }
  return(clones_mt_change)
}






