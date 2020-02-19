
source('./matlab_and_R_scripts/Clustering_EvolutionaryTree_function.R')

AA1_example <- R.matlab:::readMat('example_RobustClone.mat')
AA <- AA1_example[[1]] # Read GTM recovered by RPCA model or extended RPCA model 

robust_clone <- LJClustering(AA) # Louvain-Jaccard clustering
clone_gety <- subclone_GTM(robust_clone, 'SNV') # obtain subclonal GTM
MST <- plot_MST(clone_gety, robust_clone) # calculate and plot clonal MST
clones_mt_change <- new_mutation(robust_clone, clone_gety, MST) # get newly mutated genotypes of each subclone
