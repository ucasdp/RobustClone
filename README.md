# RobustClone
A robust PCA method of tumor clone and evolution inference from single-cell sequencing data.

RobustClone runs as follows: 
1. The input data is SNV data, either binary or ternary data. If it is binary data, 0 represents non-mutation site, 1 represents mutation site, 2 represents missing; if ternary data, 0 represents non-mutation site, 1 represents mutation heterozygous site, 2 represents mutation homozygous site and 3 represents missing，for example, ‘example.csv’ in the above list. 
2. Run the matlab script file, named "carryout.m" to recover the genotype matrix.
3. Run the R language script file, named "Cluster_EvolutionaryTree.m" to cluster cells and reconstruct the subclonal evolutionary tree.
