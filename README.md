# OpDEA v2.0.0
OpDEA v2.0.0 is an updated version of our previous published R shiny application OpDEA (https://github.com/PennHui2016/OpDEA).
It now can be used for:
1. presenting results of our benchmarking of proteomics data differential expression analysis workflow;
2. guiding users to select optimal workflows for analyzing their proteomics data;
3. For analyzing the user's proteomics data; 
4. providing links for downloading datasets used for benchmarking and for downloading our results;
5. extracting different quantification views (we define different expression matrices obatined by 
   different quantification techniques as different views, e.g., quantified by directLFQ or by MaxLFQ) 
   from quantification reports of a given quantification platform e.g., FragPipe or DIA-NN;
6. conducting differential expression analysis with our newly designed multi-view learning-based tool M-VIDIA;
7. conducting patient diagnosis when both training and testing proteomics data are provided;
8. conducting cell clustering if single cell proteomics data is provided.
   
## webserver and standalone toolkit (recommend)
We prepared a freely-accessible webserver for helping users to use OpDEA without installation of the package, 
see http://www.ai4pro.tech:3838/. The webserver requires uploading expression data, so we recommend using our
standalone toolkit (available at: https://zenodo.org/records/13911554, no R package needs to be installed and no need to upload data,
just decompress it and use it) instead. If you still hope to try our R package, please see following installation instructions.

## Requirements for installation of OpDEA v2.0.0
You should install the following packages with the same versionor higher:

R base: R-4.3.1;
R packages: shiny 1.9.1; shinydashboard 0.7.2; threejs 0.3.3; DT 0.33; ggplot2 3.5.1; reshape2 1.4.4; ggpubr 0.6.0; ggsci 3.2.0; readxl 1.4.3; ggalluvial 0.12.5; golem 0.5.1;
iq 1.9.12; dplyr 1.1.4; aggregation 1.0.1; stringr 1.5.1; tidyverse 2.0.0; matrixStats 1.4.1; readr 2.1.5; rrcovNA 0.5.2; BiocManager 1.30.25; NormalyzerDE 1.20.0; limma 3.58.1;
ROTS 1.30.0; MSnbase 2.28.1; edgeR 4.0.16; proDA 1.16.0; DEqMS 1.20.0; plgem 1.70.0; DEP 1.24.0; MSstats 4.10.1; samr 3.0.0; mice 3.16.0; missForest 1.5; SeqKnn 1.0.1; 
GMSimpute 0.0.1.0 (see https://github.com/wangshisheng/NAguideR); umap 0.2.10.0; UpSetR 1.4.0; 

The whole R session of my R environment are as follows:

```
sessionInfo()
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 22631)

Matrix products: default


locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Asia/Singapore
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] zen4R_0.10           NormalyzerDE_1.20.0  rrcovNA_0.5-2        rrcov_1.7-6         
 [5] robustbase_0.99-4-1  SeqKnn_1.0.1         GMSimpute_0.0.1.1    MSnbase_2.28.1      
 [9] ProtGenerics_1.34.0  S4Vectors_0.40.2     mzR_2.36.0           Rcpp_1.0.13         
[13] Biobase_2.62.0       BiocGenerics_0.48.1  limma_3.58.1         umap_0.2.10.0       
[17] UpSetR_1.4.0         ggpubr_0.6.0         reshape2_1.4.4       ggrepel_0.9.6       
[21] ggplot2_3.5.1        iq_1.9.12            shinydashboard_0.7.2 shiny_1.9.1         

loaded via a namespace (and not attached):
  [1] splines_4.3.1               later_1.3.2                 norm_1.0-11.1              
  [4] bitops_1.0-9                tibble_3.2.1                preprocessCore_1.64.0      
  [7] XML_3.99-0.17               lifecycle_1.0.4             rstatix_0.7.2              
 [10] doParallel_1.0.17           vroom_1.6.5                 processx_3.8.4             
 [13] lattice_0.22-6              MASS_7.3-60                 backports_1.5.0            
 [16] magrittr_2.0.3              sass_0.4.9                  jquerylib_0.1.4            
 [19] yaml_2.3.10                 remotes_2.5.0               httpuv_1.6.15              
 [22] zip_2.3.1                   askpass_1.2.1               sessioninfo_1.2.2          
 [25] pkgbuild_1.4.4              reticulate_1.39.0           cowplot_1.1.3              
 [28] MsCoreUtils_1.14.1          golem_0.5.1                 abind_1.4-8                
 [31] pkgload_1.4.0               zlibbioc_1.48.2             GenomicRanges_1.54.1       
 [34] purrr_1.0.2                 RCurl_1.98-1.16             atom4R_0.3-3               
 [37] GenomeInfoDbData_1.2.11     OpDEA_0.0.0.9000            IRanges_2.36.0             
 [40] keyring_1.3.2               RSpectra_0.16-2             DelayedArray_0.28.0        
 [43] ncdf4_1.23                  codetools_0.2-20            xml2_1.3.6                 
 [46] DT_0.33                     tidyselect_1.2.1            shape_1.4.6.1              
 [49] farver_2.1.2                rdflib_0.2.9                matrixStats_1.4.1          
 [52] roxygen2_7.3.2              jsonlite_1.8.9              ellipsis_0.3.2             
 [55] Formula_1.2-5               survival_3.7-0              iterators_1.0.14           
 [58] systemfonts_1.1.0           foreach_1.5.2               tools_4.3.1                
 [61] ragg_1.3.3                  glue_1.8.0                  SparseArray_1.2.4          
 [64] gridExtra_2.3               xfun_0.48                   MatrixGenerics_1.14.0      
 [67] usethis_3.0.0               GenomeInfoDb_1.38.8         dplyr_1.1.4                
 [70] withr_3.0.1                 BiocManager_1.30.25         fastmap_1.2.0              
 [73] fansi_1.0.6                 openssl_2.2.2               callr_3.7.6                
 [76] digest_0.6.37               R6_2.5.1                    mime_0.12                  
 [79] textshaping_0.4.0           colorspace_2.1-1            config_0.3.2               
 [82] utf8_1.2.4                  tidyr_1.3.1                 generics_0.1.3             
 [85] httr_1.4.7                  S4Arrays_1.2.1              htmlwidgets_1.6.4          
 [88] pkgconfig_2.0.3             gtable_0.3.5                impute_1.76.0              
 [91] XVector_0.42.0              pcaPP_2.0-5                 htmltools_0.5.8.1          
 [94] carData_3.0-5               profvis_0.4.0               MALDIquant_1.22.3          
 [97] clue_0.3-65                 scales_1.3.0                png_0.1-8                  
[100] attempt_0.3.1               knitr_1.48                  rstudioapi_0.16.0          
[103] tzdb_0.4.0                  curl_5.2.3                  cachem_1.1.0               
[106] stringr_1.5.1               parallel_4.3.1              miniUI_0.1.1.1             
[109] mzID_1.40.0                 vsn_3.70.0                  desc_1.4.3                 
[112] pillar_1.9.0                grid_4.3.1                  vctrs_0.6.5                
[115] pcaMethods_1.94.0           urlchecker_1.0.1            promises_1.3.0             
[118] car_3.1-3                   xtable_1.8-4                cluster_2.1.6              
[121] readr_2.1.5                 mvtnorm_1.3-1               cli_3.6.3                  
[124] compiler_4.3.1              rlang_1.1.4                 crayon_1.5.3               
[127] ggsignif_0.6.4              labeling_0.4.3              ps_1.8.0                   
[130] affy_1.80.0                 plyr_1.8.9                  fs_1.6.4                   
[133] stringi_1.8.4               redland_1.0.17-18           BiocParallel_1.36.0        
[136] assertthat_0.2.1            munsell_0.5.1               devtools_2.4.5             
[139] glmnet_4.1-8                Matrix_1.5-4.1              hms_1.1.3                  
[142] bit64_4.5.2                 statmod_1.5.0               SummarizedExperiment_1.32.0
[145] fontawesome_0.5.2           broom_1.0.7                 memoise_2.0.1              
[148] affyio_1.72.0               bslib_0.8.0                 bit_4.5.0                  
[151] DEoptimR_1.1-3 
```
## Installation

It can be installed via two ways:
1. Install the package "devtools" if you have not installed it before

    ```
    if(!requireNamespace("devtools")){
       install.packages("devtools")
    }
    ```

Then, the package can be installed from github via the following code:

    library(devtools)
    install_github('PennHui2016/OpDEA')
   

2.Or via downloading **"OpDEA_v2.0.0"** from this site, then installed with the following command:
    
    install.packages(pkgs = '~/OpDEA_0.0.0.9000.tar.gz', repos = NULL, type = "source")

3. download the python environment from Zendo: https://zenodo.org/records/13911848
   decompress the zip file to the folder of OpDEA package in the R library, e.g., /R/library/OpDEA/app/www/
   !!! the whole folder of "directlfq/" should be copyed to the /R/library/OpDEA/app/www/ folder, otherwise
   errors will be returned.
   
At last, the shiny app can be launched via:

    OpDEA::run_app()

If success, the page showing the introduction of our OpDEA and M-VIDIA will be presented. It can be used according to the contents in the help page.

## Cite these article
Hui Peng, He Wang, Weijia Kong, Jinyan Li*, Wilson Wen Bin Goh*. (2024). Optimizing differential expression analysis for proteomics data via high-performing rules and ensemble inference. Nat Commun 15, 3922 (2024). https://doi.org/10.1038/s41467-024-47899-w

Hui Peng, Wilson Wen Bin Goh*. (2024). M-VIDIA: A Multi-View representation to Increase modality Depth using Integrative AI.
 
## Contact
Any problems or requesting source codes for reproducing results in our paper please contact

    Hui Peng: hui.peng@ntu.edu.sg or cdph2009@163.com;  Wilson Wen Bin Goh: wilsongoh@ntu.edu.sg

