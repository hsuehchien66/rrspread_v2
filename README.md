RRspread: Quantifying Spread of Bacteria Using Relative Risk Ratio Framework
================

Last update: 7 February 2024


# Required Input Data:

1. Tree file(s) generated from Gubbins (Recombination-corrected)
2. Metadata including following columns:
  - GPSC_PoPUNK2 (Strain)
  - Region (Categorical geo-information)
  - Longitude
  - Latitude
  - Year_ollection
  
------------------------------------------------------------------------

# Overview of functions in the packages:

**Match_tree_meta**: Match tree label to metadata label

**Time_resolved**: generate estimated divergence time between isolates

**RR_matrix**: generate pairwise matrices from metadata file

**Geo_Strain_Calc_RR**: calculate relative risk ratio

**RR_bootstrap**: generate confidence interval through bootstrap resampling

**PltRR**: generate relative risk ratio plot (for only one calculation)

------------------------------------------------------------------------

## Workflow:

1.  using *Match_tree_meta* to generate a table that match the metadata and tree file
2.  using *Time_resolved* to generate bactdating results, including pairwise divergence time between samples in each GPSC
3.  using *bind_multitree* to bind all the trees and combine with metadata (optional)
4.  using *RR_matrix* to generate input matrices for relative risk ratio calculation
5.  using *RR_bootstrap* to resample the isolates and calculate relative risk ratio (by *Geo_Strain_Calc_RR*)
6.  divergence time v.s. distance
7.  rolling divergence time v.s. relative ratio (fixed distance)
8.  rolling distance v.s. relative ratio (fixed divergence time)

Reference: <https://xavierdidelot.github.io/BactDating/articles/yourData.html>

------------------------------------------------------------------------

# Example Script on Israel Pneumococcal Dataset

## Install RRspread package from github
```r
devtools::install_github("hsuehchien66/rrspread_v2")
```



