# TE_feature_density
A script for computing density or number of transposable elements (TEs) in protein coding features.

---

## Dependencies
  * dplyr
  * stringr
  * reshape2
  * bedtoolsr
  * AnnotationHub
  * rtracklayer
  * biomartr
  * UCSCRepeatMasker
  * data.table
  * tools

`bedtoolsr` is an R package that uses internally `bedtools` so [this package](https://bedtools.readthedocs.io/en/latest/) has to be installed previously.

**TE_feature_density** has only been tested in Unix/Linux.

---

## How does it work?


**TE_feature_density** computes the number of TEs or overlapping TE density in different gene regions (gene, exons, introns, 5'UTR, 3'UTR, downstream). You can use your own gene and RepeatMasker annotation or let **TE_feature_density** download them from Ensembl. Furthermore, a subset from gene annotation (defined by *tag*; "basic", "Ensembl_canonical", "MANE_select"...) can be chosen. 


### Automatic annotation download
Gene annotation is downloaded from Ensemble database using function `getGTF` from `biomartr`package. You need to define species (*Canis lupus familiaris*) and Ensembl release (*113*). You can check available species and releases in [this website](https://www.ensembl.org/info/website/archives/assembly.html). You can further filter annotation by tag.

TE annotation is downloaded from UCSC through `AnnotationHub` package using metadata from `UCSCRepeatMasker`package. For this purpose you must know the corresponding code to a UCSC genome version. 

To know what tags are available for a specific species and release in a gene annotation, and code for UCSC RepeatMasker annotation, you can use `checkAnnotations.R` (see **Usage**).


### Using your own annotations
Gene annotation can be downloaded from Ensembl.org. 

Remember to change in `TE_feature_density.R` these variables:

* `fileGene <- TRUE`  
* `gene_annot_file  <- "data/your_genome_annotation.gtf"`
* `gene_annot_format <- "gff"`


TE annotation can be downloaded from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables):
choose Clade, Genome and Assembly of interest, then 

* Group: Variations and Repeats,
* Track: RepeatMasker
* Output format: All fields from selected table.
* Output field separator: tsv

and Get output, saving file in Data folder. 
Remember to change in `TE_feature_density.R` these variables:

* `fileTE <- TRUE`
* `TE_annot_file <- "data/your_rmsk.txt"`
* `TE_annot_format <- "rmsk"`

---

## Usage

### Check genome and TE annotation
Use script `checkAnnotations.R` to get possible filter tags for gene annotation
(or NULL for no tag filtering) and UCSC RepeatMasker code for online TE annotation
downloading. Change variables in `Parameters`section:

* `organism <- "Canis lupus familiaris"`
* `release <- "113"`
* `fileTE <- TRUE`: TRUE for using your own TE annotation. FALSE for automatic downloading.
* `fileGene <- TRUE`: TRUE for using your own gene annotation. FALSE for automatic downloading.
* `gene_annot_file  <- "data/Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf"`


### Change parameters
Before you run `TE_feature_density.R` script you have to change R variables in section 
**Parameters**:

* `organism`: Species name. E.g "Homo sapiens", "Mus musculus", "Danio rerio"
* `UCSC_TE_annot`: UCSC code for TE annotation. You can get this code running first `checkAnnotations.R` script.
* `release`: Ensembl gene annotation version.
* `interest_TEs`: a list of TE classes to analyze. E.g. `c("LINE", "SINE")`. If set to `NULL` all TE classes are analyzed.
* `interest_subF`: a list of TE families to analyze. E.g. `c("Alu", "L1")`. If ser to `NULL` all TE families are considered. Please, notice that if this variable is not NULL overrides `interest_TEs' variable (you can only filter classes or families).
* `tag`: Filter gene annotation according to gene selection ("Ensembl_canonical", "basic", "MANE_select"...). Check tags available running first `checkAnnotations.R` script.
* `minOverlap`: only consider TEs that overlap at least `minOverlap` bp.
   In density analysis this applies to overlapping TE clusters, not individual TEs.
* `downstream`: Number of bp defining downstream region.
* `OUTPUT_DIR`: Results folder.
* `analysis`: You can choose to analyze number of TEs (`number`) or TE density (`density`). In density analysis, overlapping TEs are merged so that common nucleotides are not counted several times. 

You may also consider to change some variables in **Advanced parameters** section, particularly if you want to use your own downloaded annotations.

* `fileGene`: `TRUE` for reading file annotation from `gene_annot_file`. `FALSE` for automatic downloading.
* `gene_annot_file`: Path to gene annotation file.
* `gene_annot_format`: Parameter to import file function. Don't change it if you are not really sure!
* `fileTE`: `TRUE` for reading file annotation from `TE_annot_file`. `FALSE` for automatic downloading.
* `TE_annot_file`: Path to TE annotation file.
* `TE_annot_format`: Parameter to import file function. Don't change it if you are not really sure!
* `feature_types`: List of gene features to analyze. Choose from  c("gene", "five_prime_utr", "three_prime_utr", "exon", "intron", "downstream"). 


Please, notice that using your own annotation files can take longer time than expected!

---

## Session info
```
R version 4.3.0 (2023-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

locale:
 [1] LC_CTYPE=es_ES.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=es_ES.UTF-8        LC_COLLATE=es_ES.UTF-8    
 [5] LC_MONETARY=es_ES.UTF-8    LC_MESSAGES=es_ES.UTF-8   
 [7] LC_PAPER=es_ES.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=es_ES.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Madrid
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets
[7] methods   base     

other attached packages:
 [1] data.table_1.16.4       BRGenomics_1.14.1      
 [3] rtracklayer_1.62.0      GenomicRanges_1.54.1   
 [5] biomartr_1.0.7          UCSCRepeatMasker_3.15.2
 [7] GenomeInfoDb_1.38.8     IRanges_2.36.0         
 [9] S4Vectors_0.40.2        AnnotationHub_3.10.1   
[11] BiocFileCache_2.10.2    dbplyr_2.5.0           
[13] BiocGenerics_0.48.1     dplyr_1.1.4            
[15] stringr_1.5.1           reshape2_1.4.4         
[17] bedr_1.0.7             

loaded via a namespace (and not attached):
  [1] DBI_1.2.3                     bitops_1.0-9             
  [3] formatR_1.14                  testthat_3.2.1.1         
  [5] biomaRt_2.58.2                rlang_1.1.4               
  [7] magrittr_2.0.3                matrixStats_1.4.1         
  [9] compiler_4.3.0                RSQLite_2.3.8                
 [11] png_0.1-8                     vctrs_0.6.5                  
 [13] pkgconfig_2.0.3               crayon_1.5.3                 
 [15] fastmap_1.2.0                 XVector_0.42.0               
 [17] Rsamtools_2.18.0              promises_1.3.0               
 [19] rmarkdown_2.29                tzdb_0.4.0                   
 [21] purrr_1.0.2                   bit_4.5.0.1                  
 [23] xfun_0.49                     zlibbioc_1.48.2              
 [25] cachem_1.1.0                  jsonlite_1.8.9               
 [27] progress_1.2.3                blob_1.2.4                   
 [29] later_1.3.2                   DelayedArray_0.28.0          
 [31] BiocParallel_1.36.0           interactiveDisplayBase_1.40.0
 [33] parallel_4.3.0                prettyunits_1.2.0            
 [35] R6_2.5.1                      stringi_1.8.4                
 [37] brio_1.1.5                    knitr_1.49                   
 [39] Rcpp_1.0.13-1                 SummarizedExperiment_1.32.0  
 [41] downloader_0.4                R.utils_2.12.3               
 [43] readr_2.1.5                   VennDiagram_1.7.3            
 [45] httpuv_1.6.15                 Matrix_1.6-5                 
 [47] tidyselect_1.2.1              rstudioapi_0.17.1            
 [49] abind_1.4-8                   yaml_2.3.10                  
 [51] codetools_0.2-20              curl_6.0.1                   
 [53] lattice_0.22-6                tibble_3.2.1                 
 [55] plyr_1.8.9                    withr_3.0.2                  
 [57] Biobase_2.62.0                shiny_1.9.1                  
 [59] KEGGREST_1.42.0               evaluate_1.0.1               
 [61] lambda.r_1.2.4                futile.logger_1.4.3          
 [63] xml2_1.3.6                    Biostrings_2.70.3            
 [65] pillar_1.10.0                 BiocManager_1.30.25          
 [67] filelock_1.0.3                MatrixGenerics_1.14.0        
 [69] renv_1.0.11                   generics_0.1.3               
 [71] vroom_1.6.5                   RCurl_1.98-1.16              
 [73] BiocVersion_3.18.1            hms_1.1.3                    
 [75] ggplot2_3.5.1                 munsell_0.5.1                
 [77] scales_1.3.0                  xtable_1.8-4                 
 [79] glue_1.8.0                    tools_4.3.0                  
 [81] BiocIO_1.12.0                 locfit_1.5-9.10              
 [83] GenomicAlignments_1.38.2      XML_3.99-0.17                
 [85] grid_4.3.0                    colorspace_2.1-1             
 [87] AnnotationDbi_1.64.1          GenomeInfoDbData_1.2.11      
 [89] restfulr_0.0.15               cli_3.6.3                    
 [91] rappdirs_0.3.3                futile.options_1.0.1         
 [93] S4Arrays_1.2.1                gtable_0.3.6                 
 [95] R.methodsS3_1.8.2             DESeq2_1.42.1                
 [97] digest_0.6.37                 SparseArray_1.2.4            
 [99] rjson_0.2.23                  memoise_2.0.1                
[101] htmltools_0.5.8.1             R.oo_1.27.0                  
[103] lifecycle_1.0.4               httr_1.4.7                   
[105] mime_0.12                     bit64_4.5.2  
```