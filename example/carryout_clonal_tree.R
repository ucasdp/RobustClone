
setwd('./RobustClone/example')
source('./RobustClone/matlab_and_R_scripts/Clustering_EvolutionaryTree_function.R')

## example for SNV data
AA1_example <- R.matlab:::readMat('example_RobustClone.mat')
AA <- AA1_example[[1]] # Read GTM recovered by RPCA model or extended RPCA model 

robust_clone <- LJClustering(AA) # Louvain-Jaccard clustering
clone_gety <- subclone_GTM(AA, robust_clone, 'SNV') # obtain subclonal GTM
MST <- plot_MST(clone_gety, robust_clone, 'SNV', 'exampledata') # calculate and plot clonal MST
clones_mt <- new_variant_site(clone_gety, MST, 'SNV') # obtain the variant SNV loci each subclone compared with its parent subclone


## example for CNV data
AA1_example <- R.matlab:::readMat('SA501X3F_RobustClone.mat')
AA <- AA1_example[[1]] # Read copy number profiles recovered by RPCA model or extended RPCA model 

robust_clone <- LJClustering(AA) # Louvain-Jaccard clustering
clone_gety <- subclone_GTM(AA, robust_clone, 'CNV') # obtain subclonal copy number profiles
MST <- plot_MST(clone_gety, robust_clone, 'CNV', 'exampledata') # calculate and plot clonal MST
clones_vt <- new_variant_site(clone_gety, MST, 'CNV') # obtain the CNV genome fragments of each subclone compared with its parent subclone

data_info <- read.csv('SA501X3F.integer_copy_number.csv',header = TRUE) # the chromosome and copy number profiles information of SA501X3F data
chr <- data_info$chr # the chromosomes in which each genome fragment is located

clones_change_chr <- clonal_CNV_chr(clones_vt, chr) # obtain the variant chromosomes in which all genome fragments with CNV in clones_vt variable are located
clones_vt_state <- new_CNV_chr_state(clone_gety, MST, chr, clones_change_chr) # get the statistical frequency of how many copy number of each genome segment have changed for each chromosome and chromosome states 



