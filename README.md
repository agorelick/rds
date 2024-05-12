# Root diversity score (RDS)

R implementation of root diversity score (RDS) calculation (Reiter et al, Nature Genetics, 2020). Translated from the original python by Martin Blohmer. Packaged into an R library and maintained by Alexander Gorelick.

**Please cite**: Reiter et al., Nature Genetics, 2020. https://doi.org/10.1038/s41588-020-0633-2

## Quick installation from GitHub

You can quickly install this package directly from the github repository using the 'devtools' R package:
```r
## Installing directly from github with the 'devtools' package
install.packages('devtools')
devtools::install_github('agorelick/rds')
```

## Alternative installation from source code

You can also install 'rds' the package from the source code. First, in your terminal, clone the github repository:
```
## clone the repo in bash
git clone https://github.com/agorelick/rds
```

Then start R and install the package:
```r
## open R and install the package from the cloned repo
install.packages('rds/rds_0.1.0.tar.gz',type='src',repos=NULL)
```

# Using the RDS package

## Calculating RDS given k,l,m values
RDS values for example trees in Figure 2a-d of the publication can be calculated as follows:

![alt text](https://github.com/agorelick/rds/blob/main/etc/fig2a-d.png "Examples of RDS calculations from Figure 2a-d.")

```r
## Fig. 2a
calculate_rds(k=9,l=2,m=2) # 0.05263158

## Fig. 2b
calculate_rds(k=6,l=5,m=5) # 0.001667064

## Fig. 2c
calculate_rds(k=6,l=2,m=2) # 0.07692308

## Fig. 2d
calculate_rds(k=4,l=2,m=5) # 0.5897436
```

## Calculate RDS values given a cancer phylogeny (phylo object)

```r
## create a random tree with 20 tumor samples (+ 1 normal)
set.seed(123)
size <- 21
tree <- rtree(size)
tips <- tree$tip.label

## randomly assign tips to be the normal (N1), primary (PT) or metastasis (M) samples
primary_samples <- tips[sample(2:(size-3))]
remaining_samples <- tips[!tips %in% primary_samples]
normal_sample <- sample(remaining_samples,1)
metastasis_samples <- remaining_samples[!remaining_samples %in% normal_sample]
tree$tip.label[tree$tip.label==normal_sample] <- 'N1'
tree$tip.label[tree$tip.label %in% primary_samples] <- paste0('PT',1:length(primary_samples))
tree$tip.label[tree$tip.label %in% metastasis_samples] <- paste0('M',1:length(metastasis_samples))

## calculate RDS values for this tree
rds(tree) # RDS for M samples: 0.002316602
```



