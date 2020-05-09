# RobustClone
A robust PCA method of tumor clone and evolution inference from single-cell sequencing data. We many thank to the [article](http://arxiv.org/abs/1109.0367) for their available source code.

## Usage
RobustClone runs as follows: 
1. The input data is SNV data or CNV data. For SNV data, it can be either binary or ternary data. If it is binary data, 0 represents non-mutation site, 1 represents mutation site, 3 represents missing; if ternary data, 0 represents non-mutation site, 1 represents mutation heterozygous site, 2 represents mutation homozygous site and 3 represents missing，for example, ‘example.csv’ in the above list. 
2. Run the matlab script file, named "carryout_RPCA.m" to recover the genotype matrix.
3. Run the R language script file, named "carryout_clonal_tree.R" to cluster cells and reconstruct the subclonal evolutionary tree.

## Citation
Please cite the RobustClone in your publications if it helps your research.
```
@article{chen2019robustclone,
  title={RobustClone: A robust PCA method for tumor clone and evolution inference from single-cell sequencing data},
  author={Chen, Ziwei and Gong, Fuzhou and Wan, Lin and Ma, Liang},
  journal={Bioinformatics},
  year={2019}
}

```
Please cite this project in your publications if it helps your research.
```
@misc{robustclone,
    author = {Chen, Ziwei and Gong, Fuzhou and Ma, Liang and Wan, Lin},
    title = {RobustClone},
    howpublished = {\url{https://github.com/ucasdp/RobustClone}},
    year ={2019}
}
```

