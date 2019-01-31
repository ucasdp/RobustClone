# RobustClone
A robust PCA method of tumor clone and evolution inference from single-cell sequencing data.
RobustClone runs as follows: 
1. The input data is SNV data, either binary or ternary data. If it is binary data, 0 represents non-mutation site, 1 represents mutation site, 2 represents missing; if ternary data, 0 represents non-mutation site, 1 represents mutation heterozygous site, 2 represents mutation homozygous site and 3 represents missing，for example, ‘example.csv’ in the above list. 
2. The variable AA1 output by ‘robubstclone.m’ is the low rank genotype of cells. Variable AA_2 is clone's genotype, variable AA1_3 is a simplified and ordered matrix, clone_genes stores the mutation indexes contained in each mutation class, clone_cells stores the cell indexes contained in each clone, clone_mutation stores the mutation sites for each clone.
3. Under the ininite site model, if the clonal matrix is fully recoverded, we provide a huffman coding schim, such that, the mutational order canbe traced. 
