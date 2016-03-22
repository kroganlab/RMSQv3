# RMSQv3
RMSQ for MSstats v3.3.8



### Installing MSStats
```
# step 1: install dependency packages
install.packages(c("gplots","lme4","ggplot2","ggrepel","reshape","reshape2","data.table","Rcpp","survival"))

source("http://bioconductor.org/biocLite.R")

biocLite(c("limma","marray","preprocessCore","MSnbase"))

# step 2: install MSstats - please chage the 'LocalPath' with your location.
install.packages(pkgs="LocalPath/MSstats_3.3.8.tar.gz",repos=NULL,type="source")

# step 3: load the library of MSstats
library(MSstats)

# step 4: getting started
?MSstats
```

