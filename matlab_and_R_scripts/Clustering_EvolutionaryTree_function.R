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
  if(m <= 100){
    nn <- 10
  }
  if(m > 100 & m <= 1000){
    nn <- 30
  }
  if(m > 1000){
    nn <- 150
  }
  
  nearest <- nn2(AA, AA, k=nn+1, searchtype="priority")
  nearest$nn.idx <- nearest$nn.idx[,-1]
  nearest$nn.dists <- nearest$nn.dists[,-1] #Convert to a similarity score
  nearest$nn.sim <-  1*(nearest$nn.dists >= 0)
  
  edges <- melt(t(nearest$nn.idx))
  colnames(edges) <-  c("B", "A", "C")
  edges <-  edges[,c("A","B","C")]
  edges$B <- edges$C; 
  edges$C <- 1
  
  #Remove repetitions
  edges <-  unique(transform(edges, A = pmin(A,B), B = pmax(A,B)))
  
  Adj <-  matrix(0, nrow=nrow(AA), ncol=nrow(AA))
  rownames(Adj) <-  rownames(AA)
  colnames(Adj) <-  rownames(AA)
  Adj[cbind(edges$A, edges$B)] <-  edges$C
  Adj[cbind(edges$B, edges$A)] <-  edges$C
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
#  AA: the genotype matrix recovered by RPCA.
#  robust_clone: the clustering result output by LJClustering function;
#  type: the data type ('SNV' or 'CNV') of input. 
#        Here, 'CNV' data element is the number of copies of each chromosome segment in each cell, such as: 0, 1, 2, 3, 4, 5 , ……, 
#            where copy number 2 is normal,
#       'SNV' data element is binary or ternary, that is 0, 1 or 0, 1, 2, 
#            where 0 represents normal, 1 in binary data represents mutation under the hypothesis of infinite site, and 1,2 in ternary data represent mutation under the hypothesis of finite site;
# Output: 
#  clone_gety: the inferred clonal genotype based on the Louvain-Jaccard clustering.

subclone_GTM <- function(AA, robust_clone, type){
  clone_gety <- matrix(0,length(robust_clone),ncol(AA))

  #SNV data
  if(type == 'SNV'){
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
  if(type == 'CNV'){
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
#  type: the data type ('SNV' or 'CNV') of input.
#        Here, 'CNV' data element is the number of copies of each chromosome segment in each cell, such as: 0, 1, 2, 3, 4, 5 , ……, 
#            where copy number 2 is normal,
#       'SNV' data element is binary or ternary, that is 0, 1 or 0, 1, 2, 
#            where 0 represents normal, 1 in binary data represents mutation under the hypothesis of infinite site, and 1,2 in ternary data represent mutation under the hypothesis of finite site;
#  pdf_name: the name of pdf with the clonal MST graph.
# Output: 
#  the pdf with the clonal MST graph;
#  el: the connected edges in the MST.

plot_MST <- function(clone_gety, robust_clone, type, pdf_name){
  if(type == 'SNV'){
    clone_gety_root <- matrix(0,(length(robust_clone)+1),ncol(clone_gety))
  }
  if(type == 'CNV'){
    clone_gety_root <- matrix(2,(length(robust_clone)+1),ncol(clone_gety))
  }
  
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
    size[i] <- length(robust_clone[[i]])/sum(lengths(robust_clone))*250
  }
  
  pdf(paste('RobustClone_MST_', pdf_name, '.pdf', sep=''))
  plot(minspantree, layout=layout_as_tree, vertex.color=col, edge.color='black', vertex.label.color='black', alpha=0.5,
       edge.width=4, vertex.size=size, vertex.shape= "sphere", vertex.label=labelname)
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


########## obtain the new variant SNV loci or CNV genome fragments of each subclone compared with its parent subclone or normal cell or the root/first subclone of evolutionary tree ################
# Input:
#  clone_gety: the inferred clonal genotype based on the Louvain-Jaccard clustering output by subclone_GTM function;
#  el: the connected edges in the MST output by plot_MST function.
#  type: the data type ('SNV' or 'CNV') of input.
#        Here, 'CNV' data element is the number of copies of each chromosome segment in each cell, such as: 0, 1, 2, 3, 4, 5 , ……, 
#            where copy number 2 is normal,
#       'SNV' data element is binary or ternary, that is 0, 1 or 0, 1, 2, 
#            where 0 represents normal, 1 in binary data represents mutation under the hypothesis of infinite site, and 1,2 in ternary data represent mutation under the hypothesis of finite site;
#  compare.to: the comparison objects ('parent' or 'normal') of input.
#        Here, 'parent'represents that the new variant SNV loci or CNV genome fragments of each subclone are obtained compared with its parent subclone;
#       'normal' represents that the new variant SNV loci or CNV genome fragments of each subclone are obtained compared with normal cell;
#       'root' represents that the new variant SNV loci or CNV genome fragments of each subclone are obtained compared with the root/first subclone of evolutionary tree;
# Output: 
#  clones_vt: a list variable where each component contains the gene loci with SNV or genome fragments with CNV of each subclone compared with its parent subclone or normal cell or the root/first subclone of evolutionary tree.

new_variant_site <- function(clone_gety, el, type, compare.to){
  if(type == 'SNV'){
    normal_state <- 0
  }
  if(type == 'CNV'){
    normal_state <- 2
  }
  root_vt <- which(clone_gety[el[1,1],]!= normal_state)
  clones_vt <- list()
  clones_vt <- c(clones_vt, list(root_vt))
  names(clones_vt)[1] <- paste('subclone', el[1,1], sep='')
  for(i in 1:nrow(el)){
    clone1 <- el[i,1]
    clone2 <- el[i,2]
    clone1_clone_gety <- clone_gety[clone1,]
    clone2_clone_gety <- clone_gety[clone2,]
    if(compare.to == 'parent'){
      new_vt <- which((clone2_clone_gety - clone1_clone_gety)!=0)
    }
    if(compare.to == 'normal'){
      new_vt <- which(clone2_clone_gety!= normal_state)
    }
    if(compare.to == 'root'){
      new_vt <- which((clone2_clone_gety - clone_gety[el[1,1],])!=0)
    }
    clones_vt <- c(clones_vt, list(new_vt))
    names(clones_vt)[i+1] <- paste('subclone', as.numeric(clone2), sep='')
  }
  return(clones_vt)
}


##################### obtain the variant chromosomes of each subclone compared with its parent subclone or normal cell or the root/first subclone of evolutionary tree ##################################
# Input:
#  clones_vt: a list variable where each component contains the genome fragments with CNV of each subclone compared with its parent subclone or normal cell or the root/first subclone of evolutionary tree, which is output by new_variant_site function;
#  chr: the chromosomes in which each genome fragment is located;
# Output: 
#  clones_CNV_chr: a list variable where each component contains the variant chromosomes in which all genome fragments with CNV in clones_vt variable are located.

clonal_CNV_chr <- function(clones_vt, chr){
  clones_CNV_chrs <- list()
  for(i in 1:length(clones_vt)){
    clones_CNV_chrs <- c(clones_CNV_chrs, list(chr[clones_vt[[i]]]))
  }
  names(clones_CNV_chrs) <- names(clones_vt)
  
  clones_CNV_chr <- list()
  for(i in 1:length(clones_vt)){
    clones_CNV_chr <- c(clones_CNV_chr, list(unique(clones_CNV_chrs[[i]])))
  }
  names(clones_CNV_chr) <- names(clones_vt)
  return(clones_CNV_chr)
}

############### obtain the variant chromosome sates (loss/gain) of each subclone compared with its parent subclone or normal cell or the root/first subclone of evolutionary tree ################
# Input:
#  clone_gety: the inferred clonal genotype based on the Louvain-Jaccard clustering output by subclone_GTM function;
#  el: the connected edges in the MST output by plot_MST function;
#  chr: a factor variable represents the chromosomes in which each genome fragment is located;
#  clones_CNV_chr: a list variable where each component contains the variant chromosomes in which all genome fragments with CNV in clones_vt variable are located, which is output by clonal_CNV_chr function;
#  compare.to: the comparison objects ('parent' or 'normal') of input.
#        Here, 'parent'represents that the variant chromosome sates (loss/gain) of each subclone are obtained compared with its parent subclone;
#       'normal' represents that the variant chromosome sates (loss/gain) of each subclone are obtained compared with normal cell;
#       'root' represents that the variant chromosome sates (loss/gain) of each subclone are obtained compared with the root/first subclone of evolutionary tree;
# Output: 
#  clones_vt_state: a list variable where each component is still a list variable, which contains the statistical frequency of how many copy number of each genome segment have changed for each chromosome, when compared with its parent subclone or normal cell or the root/first subclone of evolutionary treeor normal cell or the root/first subclone of evolutionary tree.
#                    If the number of genome fragments with increased copy number is more than the number of genome fragments with reduced copy number, the  the state of the chromosome is defined as gain, labeled as '+'. Conversely, it is defined as loss, labeled as '-'. 
#                    the label ('+'/'-') of the chromosome state are reflected on the name of each component.

new_CNV_chr_state <- function(clone_gety, el, chr, clones_CNV_chr, compare.to){
  normal_st <- matrix(2, 1, dim(clone_gety)[2])
  root_clone_st <- clone_gety[el[1,1],] - normal_st
  root_clone_chr_st <- list()
  root_vt_chr <- clones_CNV_chr[[1]]
  for(i in 1:length(root_vt_chr)){
    root_clone_chr_st <- c(root_clone_chr_st, list(table(root_clone_st[which(chr==root_vt_chr[i])])))
  }
  root_state_label <- c(1:length(root_clone_chr_st))
  for(i in 1:length(root_clone_chr_st)){
    chr_table <- root_clone_chr_st[[i]]
    s <- sum(chr_table[which(as.numeric(names(chr_table))>0)])-sum(chr_table[which(as.numeric(names(chr_table))<0)])
    if(s > 0){
      root_state_label[i] = '+'
    }
    if(s < 0){
      root_state_label[i] = '-'
    }
    if(s == 0){
      root_state_label[i] = '='
    }
  }
  names(root_clone_chr_st) <- paste(paste('chr', root_vt_chr, sep=''), root_state_label, sep=',')
  
  clones_vt_state <- list()
  clones_vt_state <- c(clones_vt_state, list(root_clone_chr_st))
  names(clones_vt_state)[1] <- paste('subclone', el[1,1], sep='')
  
  for(i in 1:nrow(el)){
    clone1 <- el[i,1]
    clone2 <- el[i,2]
    clone1_clone_gety <- clone_gety[clone1,]
    clone2_clone_gety <- clone_gety[clone2,]
    if(compare.to == 'parent'){
      clone_st <- clone2_clone_gety - clone1_clone_gety
    }
    if(compare.to == 'normal'){
      clone_st <- clone2_clone_gety - normal_st
    }
    if(compare.to == 'root'){
      clone_st <- clone2_clone_gety - clone_gety[el[1,1],]
    }

    clone_chr_st <- list()
    vt_chr <- clones_CNV_chr[[i+1]]
    for(j in 1:length(vt_chr)){
      clone_chr_st <- c(clone_chr_st, list(table(clone_st[which(chr == vt_chr[j])])))
    }
    state_label <- c(1:length(clone_chr_st))
    for(j in 1:length(clone_chr_st)){
      chr_table <- clone_chr_st[[j]]
      s <- sum(chr_table[which(as.numeric(names(chr_table))>0)])-sum(chr_table[which(as.numeric(names(chr_table))<0)])
      # s represents the difference between the number of genome fragments with increased copy number and the number of genome fragments with reduced copy number on this chromosome. 
      #    If s > 0, the state of this chromosome is defined as gain, that is, '+', if s < 0, it is defined as loss, that is, '-'.
      if(s > 0){
        state_label[j] = '+'
      }
      if(s < 0){
        state_label[j] = '-'
      }
      if(s == 0){
        state_label[j] = '='
      }
    }
    names(clone_chr_st) <- paste(paste('chr', vt_chr, sep=''), state_label, sep=',')
    
    clones_vt_state <- c(clones_vt_state, list(clone_chr_st))
    names(clones_vt_state)[i+1] <- paste('subclone', as.numeric(clone2), sep='')
  }
  return(clones_vt_state)
}



