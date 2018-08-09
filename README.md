RMSQv3
===

RMSQ for MSstats v3.3.10


# Installing MSStats

```
# step 1: install dependency packages
install.packages(c("gplots","lme4","ggplot2","ggrepel","reshape2","data.table","Rcpp","survival", "getopt","yaml", "pheatmap"))

# Install package from bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite(c("limma","marray","preprocessCore","MSnbase","biomaRt"))

# step 2: Install MSstats
# - Recommended option: from bioconductor like this:
source("https://bioconductor.org/biocLite.R")
biocLite("MSstats")

# - Alternatively, use the version available in this repository `MSstats_3.3.10.tar.gz`. Install it like this:
install.packages(pkgs="MSstats_3.3.10.tar.gz",repos=NULL,type="source")

# step 3: load the library of MSstats
library(MSstats)

# step 4: getting started. Check that it works
?MSstats
```

## Input files

3 tab delimited files are required to run RMSQv3:

### `evidence.txt`

The output of the quantitative proteomics software package `MaxQuant`. 
It combines all the information about the identified peptides and normally 
is the only file required for processing the results.

### `keys.txt`

It contains the experimental design. 
This file will be merged with the `evidence.txt`
file (i.e., the output of MaxQuant with the peptides identified) through the 
"Raw.file" column. 
Each raw file corresponds to a unique individual technical replicate / 
biological replicate / Condition / Run.

Example (see `test/example-keys.txt`)

**RawFile**|**IsotopeLabelType**|**Condition**|**BioReplicate**|**Run**
:-----:|:-----:|:-----:|:-----:|:-----:
FU20170922-17|L|H1N1\_03H|H1N1\_03H-1|9
FU20170922-19|L|H1N1\_03H|H1N1\_03H-2|10
FU20170922-21|L|H1N1\_06H|H1N1\_06H-1|11
FU20170922-23|L|H1N1\_06H|H1N1\_06H-2|12
FU20170922-35|L|H1N1\_12H|H1N1\_12H-1|13
FU20170922-37|L|H1N1\_12H|H1N1\_12H-2|14
FU20170922-39|L|H1N1\_18H|H1N1\_18H-1|15
FU20170922-41|L|H1N1\_18H|H1N1\_18H-2|16
FU20170922-01|L|MOCK\_03H|MOCK\_03H-1|1
FU20170922-03|L|MOCK\_03H|MOCK\_03H-2|2
FU20170922-05|L|MOCK\_06H|MOCK\_06H-1|3
FU20170922-07|L|MOCK\_06H|MOCK\_06H-2|4
FU20170922-09|L|MOCK\_12H|MOCK\_12H-1|5
FU20170922-11|L|MOCK\_12H|MOCK\_12H-2|6
FU20170922-13|L|MOCK\_18H|MOCK\_18H-1|7
FU20170922-15|L|MOCK\_18H|MOCK\_18H-2|8

### `contrast.txt`

The comparisons between conditions that we want to quantified. The written
comparisons must follow the following consensus:

```
Condition_A-Condition_B_mutant
```

- The two conditions to be compared are separated by a dash symbol (`-`)
- The condition on the left will take the positive log2FC sign (if it is more abundant) and the one on the right the negative log2FC (if it is more abundant)
- The only special character allowed for the condition names is the underscore (`_`)

Example (see `test/example-contrast.txt`)

```
H1N1_03H-MOCK_03H
H1N1_06H-MOCK_06H
H1N1_12H-MOCK_12H
H1N1_18H-MOCK_18H
```

### The configuration file (`.yaml`)

The configuration file in `yaml` format contains the details of the analyses 
performed by `RMSQv3`. Check the folder `test` for a sample configuration file 
depending on the experiment. It currently covers the quantification of 
protein abundance, phosphorylation, ubiquitination and acetylation.

The configuration file contains the following sections:

```
files :
  keys : /path/to/project/data/keys.txt
  data : /path/to/project/data/data/evidence.txt
  contrasts : /path/to/project/data/contrast.txt
  output : /path/to/project/results/20160621-results.txt
  sample_plots : 1
```

The file path / name of the required files. **sample_plots** creates quality
control plots, including heatmaps of the peptide features based on intensity
values, and peptide counts.


```
data:
  enabled : 1
```

- 1 to pre-process the data provided in the *files* section.
- 0 won't process the data (and a pre-generated MSstats file will be expected)

```
fractions: 
  enabled : 0 # 1 for protein fractions
  aggregate_fun : sum
```
Multiple fractionation or separation methods are often combined in proteomics 
to improve signal-to-noise and proteome coverage and to reduce interference
between peptides in quantitative proteomics.
Use 1 to enable processing protein fractionation datasets. See the
*Special case: Protein fractionation* section below for details.

```
silac: 
  enabled : 0 # 1 for SILAC experiments
```

Mark 1 if the files belong to a SILAC experiment. See *Special case: SILAC*
below for details

```
filters: 
  enabled : 1 # Enables filtering
  contaminants : 1 # Removes contaminants (CON__ and REV__)
  protein_groups : remove # remove or keep protein groups
  modifications :  # empty for all, PH, UB, or AC
```

Filtering the datasets:

- `contaminants` : 1 to remove contaminants (CON__ and REV__)
- `protein_groups` : `remove` or `keep` protein groups
- `modifications` :  `empty` for all, `PH` to select phospho-peptides, `UB`
ubiquitinated peptides, or `AC` to select acetylated peptides.


```
msstats :
  enabled : 1
  msstats_input : 
  version :  # blank = R library version, MSstats.daily = a location where to find it
  profilePlots : before-after # before, after, before-after, none
  normalization_method : equalizeMedians # globalStandards (include a reference protein(s) ), equalizeMedians, quantile, 0
  normalization_reference :  #should be a value in the Protein column
  summaryMethod : TMP # "TMP"(default) means Tukey's median polish, which is robust estimation method. "linear" uses linear mixed model. "logOfSum" conducts log2 (sum of intensities) per run.
  censoredInt : NA  # Missing values are censored or at random. 'NA' (default) assumes that all 'NA's in 'Intensity' column are censored. '0' uses zero intensities as censored intensity. In this case, NA intensities are missing at random. The output from Skyline should use '0'. Null assumes that all NA intensites are randomly missing.
  cutoffCensored : minFeature  # Cutoff value for censoring. only with censoredInt='NA' or '0'. Default is 'minFeature', which uses minimum value for each feature.'minFeatureNRun' uses the smallest between minimum value of corresponding feature and minimum value of corresponding run. 'minRun' uses minumum value for each run.
  MBimpute : 1 # only for summaryMethod="TMP" and censoredInt='NA' or '0'. TRUE (default) imputes 'NA' or '0' (depending on censoredInt option) by Accelated failure model. FALSE uses the values assigned by cutoffCensored.
  feature_subset: all # all|highQuality  : highQuality seems to be buggy right now
```

If enables (1), it will run MSstats with all the specified options 
(read MSstats manual to find out more)

```
output_extras :
  enabled : 1
  msstats_output : 
  annotate : 0 # 1|0 whether to annotate the proteins in the results or not
  species : HUMAN # can use multiple species, but separate with a "-" eg. HUMAN-MOUSE-HIV-...
  annotation_dir : /path/to/the/files
  comparisons : all # or any grep expression that returns a subset of the contrasts file
  LFC : -1 1
  FDR : 0.05
  heatmap : 1 
  heatmap_cluster_cols : 0
  heatmap_display : log2FC #or pvalue
  volcano : 1
```

Extra actions to perform based on the MSstats results.


#### Special case: Protein fractionation

To handle protein fractionation experiments, two options need to be activated

1. the keys file must contain and additional column named "FractionKey" with 
the information about fractionation. For example:

**Raw.file**|**IsotopeLabelType**|**Condition**|**BioReplicate**|**Run**|**FractionKey**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
S9524\_Fx1|L|AB|AB-1|1|1
S9524\_Fx2|L|AB|AB-1|1|2
S9524\_Fx3|L|AB|AB-1|1|3
S9524\_Fx4|L|AB|AB-1|1|4
S9524\_Fx5|L|AB|AB-1|1|5
S9524\_Fx6|L|AB|AB-1|1|6
S9524\_Fx7|L|AB|AB-1|1|7
S9524\_Fx8|L|AB|AB-1|1|8
S9524\_Fx9|L|AB|AB-1|1|9
S9524\_Fx10|L|AB|AB-1|1|10
S9525\_Fx1|L|AB|AB-2|2|1
S9525\_Fx2|L|AB|AB-2|2|2
S9525\_Fx3|L|AB|AB-2|2|3
S9525\_Fx4|L|AB|AB-2|2|4
S9525\_Fx5|L|AB|AB-2|2|5
S9525\_Fx6|L|AB|AB-2|2|6
S9525\_Fx7|L|AB|AB-2|2|7
S9525\_Fx8|L|AB|AB-2|2|8
S9525\_Fx9|L|AB|AB-2|2|9
S9525\_Fx10|L|AB|AB-2|2|10
S9526\_Fx1|L|AB|AB-3|3|1
S9526\_Fx2|L|AB|AB-3|3|2
S9526\_Fx3|L|AB|AB-3|3|3
S9526\_Fx4|L|AB|AB-3|3|4
S9526\_Fx5|L|AB|AB-3|3|5
S9526\_Fx6|L|AB|AB-3|3|6
S9526\_Fx7|L|AB|AB-3|3|7
S9526\_Fx8|L|AB|AB-3|3|8
S9526\_Fx9|L|AB|AB-3|3|9
S9526\_Fx10|L|AB|AB-3|3|10

Internally, the function `getMSstatsFormat` handles the key step 
(just a simple `sum` aggregation)

2. Enable *fractions* in the configuration file as follow:

```
fractions: 
  enabled : 1 # 1 for protein fractions
  aggregate_fun : sum
```

#### Special case: SILAC

One of the most widely used techniques that enable relative protein 
quantitation is stable isotope labeling by amino acids in cell culture (SILAC). 
The keys file will capture the typical SILAC experiment. 
For example, let's show a SILAC experiment with two conditions, 
two biological replicates and two technical replicates:

**RawFile**|**IsotopeLabelType**|**Condition**|**BioReplicate**|**Run**
:-----:|:-----:|:-----:|:-----:|:-----:
QE20140321-01|H|iso|iso-1|1
QE20140321-02|H|iso|iso-1|2
QE20140321-04|L|iso|iso-2|3
QE20140321-05|L|iso|iso-2|4
QE20140321-01|L|iso\_M|iso\_M-1|1
QE20140321-02|L|iso\_M|iso\_M-1|2
QE20140321-04|H|iso\_M|iso\_M-2|3
QE20140321-05|H|iso\_M|iso\_M-2|4

It is also required to activate the *silac* option in the yaml file to be 
activated as follow:

```
silac: 
  enabled : 1 # 1 for SILAC experiments
```

## Running MSstats v3

### `MSstats_main.R`

```
MSstats_main.R -c configuration_file.yaml
```

### `MaxQ_utilities.R`

These functions are designed to work with MaxQuant evidence files. The functions include a wide variety of options listed below.

The `Arguments` section lists the arguments needed for each function. The arcuments (proceeded by the short flag alias -x) should be entered in the same line of the terminal, separated by spaces.


## Typical workflows

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

##### concat  

##### convert-silac  
Converts SILAC data to a format compatible with MSstats. 

```
Arguments:
-c convert-silac
-f evidence_file_path
-o output_file_path
```


##### getRawFiles

Get the list of unique Raw.files from an evidence file.

```
Arguments:
-c getRawFiles
-f evidence_file_path
-o output_file_path
```

##### convert-sites

This function is a preprocessing function that is used with PTM data. 
If you want to run a site specific analysis, before running 
MSStats on a PH/UB/AC set, you need to convert the evidence file into a format 
that MSStats will be able to diffentiate the sites with. 
This function outputs a new file that should be used as the input file for your 
MSStats group comparison analysis.

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

##### replicateplots

Outputs a replicate plots based on a user provied file containing the replicates to be compared. Values are based on the log2 value of the maximum intensities per modified sequence. The "replicate file" should describe which replicates of which conditions should be compared against eachother. Each row represents a replicate plot to be created. The file should be structured using the following format and column names:

| condition1 | rep1_1 | rep1_2 | condition2 | rep2_1 | rep2_2 |
|---|---|---|---|---|---|
|  Infected | Infected_Rep1_name | Infected_Rep2_name | Negative | Negative_Rep1_name | Negative_Rep2_name |
| etc... |   |   |   |   |   |
| etc... |   |   |   |   |   |

The arguments for this optionare as follows:

```
Arguments:
-c replicateplots
-f evidence_file_path
-k keys_file_path
-r replicate_plot_info_file_path
-o output_file_path
```

##### simplify

##### saint-format

Converts the MaxQuant evidence file to the 3 required files for SAINTexpress. 
One can choose to either use the `spectral counts` or the `intensities` for the analysis. 

```
Arguments:
-c saint-format
-f evidence_file_path
-k keys_file_path
-o output_file_path
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

Converts MaxQuant evidence file into a file format compatible with the MiST pipeline using `MS/MS.Count`. Note that this is the MiST *data* file, and that an additional *keys* file will have to be constructed before running MiST. Multiple species can be searched at once, simply separate them by a "-". (eg. `HUMAN-MOUSE`)

```
Arguments:
-c mist
-f evidence_file_path
-k keys_file_path
-o output_file_path
-s species (e.g.: `-s HUMAN-DENGUE`)
-d uniprot_dir (directory containing uniprot file with protein length)
```

##### mistint

Very similar to the `mist`, but instead of using `MS/MS.Count`, uses `Intensity` values. Converts MaxQuant evidence file into a file format compatible with the MiST pipeline. Note that this is the MiST *data* file, and that an additional *keys* file will have to be constructed before running MiST. Multiple species can be searched at once, simply separate them by a "-". (eg. `HUMAN-MOUSE`)

```
Arguments:
-c mistint
-f evidence_file_path
-k keys_file_path
-o output_file_path
-s species (e.g.: `-s HUMAN-DENGUE`)
-d uniprot_dir (directory containing uniprot file with protein length. E.g: ~/Box Sync/db/mist/)
```


##### samplequant

Aggregates the normalized abundance and replicate data from the samples. Uses the MSstat output file  `...mss-sampleQuant.txt` for the aggregations, and is applied directly to the MSstats results in ***wide*** format. The resulting file will have "abundance" appended to the end of the file name.

```
Arguments:
-f sampleQuant_file_path
-n contrast_file_path
-r results_wide_file_path
```







