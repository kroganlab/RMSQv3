RMSQv3
===

RMSQ for MSstats v3.3.10

#### Installing MSStats
```
# step 1: install dependency packages
install.packages(c("gplots","lme4","ggplot2","ggrepel","reshape","reshape2","data.table","Rcpp","survival"))

source("http://bioconductor.org/biocLite.R")

biocLite(c("limma","marray","preprocessCore","MSnbase"))

# step 2: install MSstats - Use the version available in this repository `MSstats_3.3.10.tar.gz`. Install it like this:
install.packages(pkgs="MSstats_3.3.10.tar.gz",repos=NULL,type="source")

# step 3: load the library of MSstats
library(MSstats)

# step 4: getting started. Check that it works
?MSstats
```

### Running MSstats v3, `MSstats_main.R`

```
MSstats_main.R -c configuration_file.yaml
```

Check the folder `test` for a sample configuration file

### `MaxQ_utilities.R`
These functions are designed to work with MaxQuant evidence files. The functions include a wide variety of options listed below.

The `Arguments` section lists the arguments needed for each function. The arcuments (proceeded by the short flag alias -x) should be entered in the same line of the terminal, separated by spaces.


### Typical workflows
The RMSQ pipeline was designed to run in a certain order. The following is the propper order to perform the analysis with RMSQ for SILAC, PTM, and AMPS datasets.

##### PTM analysis
0. MaxQ_utilities -> convert-silac (for SILAC data only)
1. MaxQ_utilities -> convert-sites
2. MSstats
3. MaxQ_utilities -> results-wide
4. MaxQ_utilities -> mapback-sites
5. MaxQ_utilities -> annotate

##### APMS Global Analysis
1. MSstats
2. MaxQ_utilities -> results-wide
4. MaxQ_utilities -> annotate




### Functions
#####concat  
#####convert-silac  
Converts SILAC data to a format compatible with MSstats. 

```
Arguments:
-c convert-silac
-f evidence_file_path
-o output_file_path
```


#####keys  
#####convert-sites
This function is a preprocessing function that is used with PTM data. If you want to run a site specific analysis, before running  MSStats on a UB/PH/AC set, you need to convert the evidence file into a format that MSStats will be able to diffentiate the sites with. This outputs a new file that should be used as the input file for your MSStats group conparison analysis.

```
Arguments:
-c convert-sites
-f evidence_file_path
-o output_file_path
-p reference_proteom_file_path
-t mod_type (ph|ub|ac)
```

##### annotate
Annotates the proteins in the results or results-wide files after they've been through the MSStats pipeline. Multiple species can be searched at once, simply separate them by a "-". (eg. human-mouse)

```
Arguments:
-c annotate
-f results_file_path
-o output_file_path
-s species (human|mouse)
-d uniprot_dir (directory containing uniprot file with protein masses)
```

##### results-wide
Converts the normal MSStats output file into "wide" format where each row represents a protein's results, and each column represents the comparison made by MSStats. The fold change and p-value of each comparison will be it's own column.

```
Arguments:
-c results-wide
-f results_file_path
-o output_file_path
```

##### mapback-sites
Used with PTM datasets. Map back the sites to correct proteins after MSStats analysis. This file is created previously when running the `conver-sites` function.

```
Arguments:
-c mapback-sites
-f results_file_path
-o output_file_path
-m mapping_file_path
```

##### heatmap
Outputs a heatmap of the MSStats results created using the log2 fold changes.

```
Arguments:
-c heatmap
-f results_file_path
-o output_file_path
-l log2FC lower bound
-u log2FC upper bound
-q FDR significance cutoff
```

##### simplify
##### saint-format
Converts the MaxQuant evidence file to the 3 required files for SAINTexpress. One can choose to either use the `spectral counts` or the `intensities` for the analysis. 

```
Arguments:
-c saint-format
-f evidence_file_path
-k keys_file_path
-o output_file_directory
-p reference_proteome
-i identifier_column (sepctral_count|ms1)
```

##### data-plots
##### spectral-counts
Outputs the spectral counts from the MaxQuant evidence file. 

```
Arguments:
-c spectral-counts
-f evidence_file_path
-k keys_file_path
-o output_file_path
```

##### mist
Converts MaxQuant evidence file into a file format compatible with the MiST pipeline. Note that this is the MiST *data* file, and that an additional *keys* file will have to be constructed before running MiST. Multiple species can be searched at once, simply separate them by a "-". (eg. `HUMAN-MOUSE`)

```
Arguments:
-c mist
-f evidence_file_path
-k keys_file_path
-o output_file_path
-s species (e.g.: `-s HUMAN-DENGUE`)
-d uniprot_dir (directory containing uniprot file with protein length)
```












