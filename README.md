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

# Calculating RDS values

RDS values for example trees in Figure 2a-d of the publication can be calculated as follows:

![alt text](https://github.com/agorelick/rds/blob/main/etc/fig2a-d.png "Examples of RDS calculations from Figure 2a-d.")

```r
## Fig. 2a
rds(k=9,l=2,m=2) # 0.05263158

## Fig. 2b
rds(k=6,l=5,m=5) # 0.001667064

## Fig. 2c
rds(k=6,l=2,m=2) # 0.07692308

## Fig. 2d
rds(k=4,l=2,m=5) # 0.5897436
```

