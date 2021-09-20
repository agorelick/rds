# rds
Root diversity score (RDS)

R implementation of root diversity score (RDS) calculation (Reiter et al, Nature Genetics, 2020). Translated from the original python by Martin Blohmer. Packaged into an R library and maintained by Alexander Gorelick.

**Please cite**: Reiter et al., Nature Genetics, 2020. https://doi.org/10.1038/s41588-020-0633-2

# Installation

You can quickly install this package directly from the github repository using the 'devtools' R package:
```r
## Installing directly from github with the 'devtools' package
install.packages('devtools')
devtools::install_github('agorelick/rds')
```


Alternatively, you install 'rds' the package by first cloning the github repository:
```
## clone the repo in bash
git clone https://github.com/agorelick/rds
```

Then start R and install the package:
```r
## open R and install the package from the cloned repo
install.packages('rds/rds_0.1.0.tar.gz',type='src',repos=NULL)
```



