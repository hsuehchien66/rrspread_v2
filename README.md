RRspread: Quantifying Spread of Bacteria Using Relative Risk Ratio Framework
================

Last update: 7 February 2024

# Install RRspread package from github
```r
devtools::install_github("hsuehchien66/rrspread_v2")
```

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

# Workflow:

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

## Load libraries

```r
library(ape)
library(BactDating)
library(dplyr)
library(geosphere)
library(ggplot2)
library(data.table)
library(cowplot)
library(rrspread)
```

## Load input data
```r
metadata_file_path = "Israel_1174isolates_metadata_input.txt"
metadata = read.csv(metadata_file_path, sep="\t")
metadata = metadata[, c("lane_id", "GPSC_PoPUNK2", "REGION", "Longitude", "Latitude", "Year_collection")]
```

## Load Gubbins tree

The tree should be recombination-corrected by Gubbins.  
GPSC6,8,10,47,55 are selected for Israel migration analysis

```r
## This is the example for GPSC6, run this section for other GPSCs
gpsc6_gubbin_tree_path = "GPSC6_reference_bwa.final_tree.tre"
  gpsc6_gubbin_tree = read.tree(gpsc6_gubbin_tree_path)
```

## Run BactDating (time-consuming)

**Input:**   
1. Gubbins (recombination-corrected tree),   
2. metadata with collection time   
**Output:**   
1. bactdating results   
2. root-to-tip distance vs. collection time   
**Notes:**   
In **match_tree_meta()**, *drop_tip= TRUE* to make sure reference genome is not included in the analysis,     
In **time_resolved()**, *branch_model* can be assigned for evolutionary rates of branches,   
In **time_resolved()**, *n_it* to assign number of iterations in the MCMC process

```r
## This is the example for GPSC6, run this section for other GPSCs
gpsc6_gubbin_tree_plot = plt_gubbin_tree(gpsc6_gubbin_tree)
gpsc6_time_tree_list = match_tree_meta(metadata, gpsc6_gubbin_tree, "GPSC6_reference", drop_tip = TRUE)
gpsc6_bd_res = time_resolved(gpsc6_time_tree_list$match_tree_meta, gpsc6_time_tree_list$tree, "/Users/hc14/Documents/SpneumoIsrael/analysis_results/bd_res_GPSC6", n_it = 3000, branch_model="mixedgamma")

```

## Bind trees and combine trees with metadata

Match the laneid in the tips of the trees and the laneid in the metadata

```r
gpsc08_bd_tree = gpsc8_bd_res$BactDating$tree
gpsc55_bd_tree = gpsc55_bd_res$BactDating$tree
gpsc47_bd_tree = gpsc47_bd_res$BactDating$tree
gpsc06_bd_tree = gpsc6_bd_res$BactDating$tree
gpsc10_bd_tree = gpsc10_bd_res$BactDating$tree

multitrees = bind_multitree(gpsc08_bd_tree, gpsc55_bd_tree, gpsc47_bd_tree, gpsc06_bd_tree, gpsc10_bd_tree)
plot(multitrees, cex=0.3, no.margin = TRUE)

lanes <- multitrees$tip.label

## subset only those present in the trees
selected_gpsc_metadata = subset(metadata, metadata$lane_id %in% lanes)
lanes_gpsc <- subset(lanes, lanes  %in% selected_gpsc_metadata$lane_id)
multitrees_overall <- keep.tip(multitrees, lanes_gpsc)

## Reorder tree based on order of tip labels
order <- multitrees_overall$tip.label
selected_gpsc_metadata <- selected_gpsc_metadata %>%
  dplyr::slice(match(order, lane_id))
```

## Load Metadata of GPSC6,8,10,47,55

```r
gpsc6_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "6"),]
gpsc8_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "8"),]
gpsc10_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "10"),]
gpsc47_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "47"),]
gpsc55_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "55"),]
```

## Generate matrices for relative risk calculations for each GPSC

**Input:**   
1. Metadata,   
2. Bactdating results  
**Output:**  
11 matrices for following relative risk calculations

```r
gpsc6_rrmatrix <- RR_matrix(meta_table = gpsc6_metadata, bactdate_tree = gpsc6_bd_res$BactDating$tree, ref_distance_min = 0)
gpsc8_rrmatrix <- RR_matrix(meta_table = gpsc8_metadata, bactdate_tree = gpsc8_bd_res$BactDating$tree, ref_distance_min = 0)
gpsc10_rrmatrix <- RR_matrix(meta_table = gpsc10_metadata, bactdate_tree = gpsc10_bd_res$BactDating$tree, ref_distance_min = 0)
gpsc47_rrmatrix <- RR_matrix(meta_table = gpsc47_metadata, bactdate_tree = gpsc47_bd_res$BactDating$tree, ref_distance_min = 0)
gpsc55_rrmatrix <- RR_matrix(meta_table = gpsc55_metadata, bactdate_tree = gpsc55_bd_res$BactDating$tree, ref_distance_min = 0)
```

## Distance vs. Divergence time relationship

```r
## This is the example for GPSC6, run this section for other GPSCs
gpsc6_distance_dist <- as.data.frame(to.upper(gpsc6_rrmatrix$distance_num_mat))
gpsc6_divtime_dist <- as.data.frame(to.upper(gpsc6_rrmatrix$divtime_num_mat))
gpsc6_divtime_distance <- cbind(gpsc6_divtime_dist, gpsc6_distance_dist)
colnames(gpsc6_divtime_distance) <- c("divtime", "distance")
```

## Distance vs. Divergence time with confidence interval by sampling the posterior bactdating trees

### Procedures:
Take each GPSC  
1. Sample the posterior bactdating trees (we want to sample posterior trees sequentially through the MCMC trace)  
2. Reconstruct the sampled posterior trees  
3. Calculate mean pairwise evolutionary divergence time from all the sampled trees  
4. Identify pairs within evolutionary divergence time  
5. Sample pairs with replacement (bootstrap)  
6. Extract spatial distance per evolutionary window  
7. Calculate mean spatial  per evolutionary window

### Run Dtime_Dist_CI_posterior_tree for each GPSC
**Inputs:**   
1. Metadata,  
2. BactDating results
```r
## This is the example for GPSC6, run this section for other GPSCs
gpsc6_Dtime_Dist_CI_posterior_tree <- Dtime_Dist_CI_posterior_tree(gpsc6_metadata, gpsc6_bd_res)
```

## Plot Distance CI vs. Divergence (GPSC6)
```{r}
## This is the example for GPSC6, run this section for other GPSCs
gpsc6_dist_div <- ggplot(gpsc6_Dtime_Dist_CI_posterior_tree$window_divtime_dist, aes(x=mid_win, y=distance.mean, ymin=distance.lowerci, ymax=distance.upperci))+
  geom_line()+
  geom_ribbon(fill="pink", alpha=0.5)+
  #xlim(0,60)+
  xlab("Divergence Time \n (years)")+
  ylab("Distance (km)")+
  ylim(c(0,87))+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))
```










