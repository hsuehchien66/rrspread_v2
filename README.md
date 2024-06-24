RRspread: Quantifying Spread of Bacteria Using Relative Risk Ratio Framework
================

Last update: 24 June 2024

# Install RRspread package from github
```r
devtools::install_github("hsuehchien66/rrspread_v2")
```

# Required Input Data:

1. Tree file(s) generated from Gubbins (Recombination-corrected)
2. Metadata including following columns:
  - GPSC_PoPUNK2 (Strain)
  - REGION (Categorical geo-information)
  - Longitude
  - Latitude
  - Year_collection
  
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

Reference:   
 - [Didelot et al (2018). Bayesian inference of ancestral dates on bacterial phylogenetic trees. Nucleic Acids Research, 46:e134.](https://doi.org/10.1093/nar/gky783)  
 - [Belman S et al (2023). Geographic migration and vaccine-induced fitness changes of Streptococcus pneumoniae. bioRxiv](https://pubmed.ncbi.nlm.nih.gov/36711799/)

------------------------------------------------------------------------
# Citations:
If you use this package, please cite:
 - [Belman S et al (2023). Geographic migration and vaccine-induced fitness changes of Streptococcus pneumoniae. bioRxiv](https://pubmed.ncbi.nlm.nih.gov/36711799/)
 - [Cheng H et al (2024). Estimating geographical spread of Streptococcus pneumoniae within Israel using genomic data](https://doi.org/10.1099/mgen.0.001262)

------------------------------------------------------------------------
# Example script on Israel Pneumococcal Dataset:

Please find the Rmarkdown file here: [RRspread_toy_example.rmd](https://github.com/hsuehchien66/rrspread_v2/blob/main/scripts/RRspread_toy_example.Rmd)  
Please find the example dataset here: [RRspread_example_dataset](https://figshare.com/articles/dataset/Input_files_for_RRspread_on_Israel_Pneumococcal_Dataset/25180613)

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
gpsc6_bd_res = time_resolved(gpsc6_time_tree_list$match_tree_meta, gpsc6_time_tree_list$tree, "bd_res_GPSC6", n_it = 3000, branch_model="arc")

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
#### Dtime_Dist_CI_posterior_tree()
Take each GPSC  
1. Sample the posterior bactdating trees (we want to sample posterior trees sequentially through the MCMC trace)  
2. Reconstruct the sampled posterior trees  
3. Calculate mean pairwise evolutionary divergence time from all the sampled trees  
4. Identify pairs within evolutionary divergence time  
5. Sample pairs with replacement (bootstrap)  
6. Extract spatial distance per evolutionary window  
7. Calculate mean spatial per evolutionary window

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

## Relative risk of Spread: Lineage v.s. District
```r
rr_lineage_location_mat = RR_matrix(selected_gpsc_metadata, multitrees)

colyear_mat_ones = matrix(1, nrow=dim(selected_gpsc_metadata)[1], ncol=dim(selected_gpsc_metadata)[1])  ## don't select by collection year
rr = RR_bootstrap(metadataset = selected_gpsc_metadata,
                  des_geo_matrix = rr_lineage_location_mat$des_location_mat,
                  ref_geo_matrix = rr_lineage_location_mat$ref_location_mat,
                  des_strain_matrix = rr_lineage_location_mat$des_lineage_mat,
                  ref_strain_matrix = rr_lineage_location_mat$ref_lineage_mat,
                  colyear_mat_ones, nboot = 100, subsample_min = 70)

## load library: ggplot2
district_lineage_rr = PltRR(rr, "Within District Relative Ratio", "District")

# metadata: metadata
lineage_matrix <- RR_matrix_notree(meta_table = metadata)

colyear_mat_ones = matrix(1, nrow=dim(metadata)[1], ncol=dim(metadata)[1])  ## don't select by collection year
rr = RR_bootstrap(metadataset = metadata,
                  des_geo_matrix = lineage_matrix$des_location_mat,
                  ref_geo_matrix = lineage_matrix$ref_location_mat,
                  des_strain_matrix = lineage_matrix$des_lineage_mat,
                  ref_strain_matrix = lineage_matrix$ref_lineage_mat,
                  colyear_mat_ones, nboot = 100, subsample_min = 70)

district_lineage_rr = PltRR(rr, "Within District Relative Ratio", "District")
```


## Relative risk of Spread: Lineage v.s. Physical Distance
Rolling window from 0-30 to 120-150 km 

```r
des_distance_min1 = c(0,1)
des_distance_min2 = seq(from=5, to=25, by=5)
des_distance_min3 = seq(from=30, to=90, by=10)
des_distance_min = append(append(des_distance_min1, des_distance_min2), des_distance_min3)
des_distance_max1 = c(1,5)
des_distance_max2 = seq(from=10, to=30, by=5)
des_distance_max3 = seq(from=40, to=100, by=10)
des_distance_max = append(append(des_distance_max1, des_distance_max2), des_distance_max3)

des_distance_range = length(des_distance_max)

res_list = c()
res_df = data.frame()

colyear_mat_ones = matrix(1, nrow=dim(metadata)[1], ncol=dim(metadata)[1])
diag(colyear_mat_ones) <- NA

RRplotLists = list()
RRdistanceLists = list()
RRdframLists = list()

for(distances in 1:des_distance_range){
  distance_name = paste0("Distance: ",des_distance_min[distances],"~",des_distance_max[distances],"km")
  dist_lineage_mat = RR_matrix_notree(meta_table = metadata, 
                               des_distance_min = des_distance_min[distances],
                               des_distance_max = des_distance_max[distances], 
                               ref_distance_min = 50, ref_distance_max = 300)
  dist_lineage_mat_rr = RR_bootstrap(metadataset = metadata,
                                     des_geo_matrix = dist_lineage_mat$des_distance_mat,
                                     ref_geo_matrix = dist_lineage_mat$ref_distance_mat,
                                     des_strain_matrix = dist_lineage_mat$des_lineage_mat,
                                     ref_strain_matrix = dist_lineage_mat$ref_lineage_mat,
                                     colyear_matrix = colyear_mat_ones,
                                     nboot = 100,
                                     subsample_min = 100)
  rr_plot = PltRR(dist_lineage_mat_rr, distance_name, distance_name)
  RRplotLists[[distances]] = print(rr_plot)
  
  distances_col = data.frame(distance_min=des_distance_min[distances], distance_max=des_distance_max[distances], distance_mid = (des_distance_min[distances]+des_distance_max[distances])/2)
  new_row = cbind(dist_lineage_mat_rr, distances_col)
  res_df = rbind(res_df, new_row)
  
  assign(distance_name, rr_plot)
  res_list=c(res_list, distance_name)
}

```

## At lineage level,  [Lineage vs. District] & [Lineage vs. Physical Distance] 
```r
dist_district_lineage_rr_table <- rbind(res_df, rr, fill=TRUE)
for(i in 1:nrow(dist_district_lineage_rr_table)){
  dist_district_lineage_rr_table[i, "distance_window"] = paste0(dist_district_lineage_rr_table[i, "distance_min"],"-" ,dist_district_lineage_rr_table[i, "distance_max"])
}
dist_district_lineage_rr_table[nrow(dist_district_lineage_rr_table), "distance_window"] = "within \ndistrict"
dist_district_lineage_rr_table$distance_window = as.factor(dist_district_lineage_rr_table$distance_window)

dist_district_lineage_rr_table$distance_window <- ordered(dist_district_lineage_rr_table$distance_window, levels=c("0-1", "1-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100", "within \ndistrict"))
rolling_distance_RR_plot_50km <- ggplot(dist_district_lineage_rr_table, aes(x=distance_window, y=RR_mean))+
  geom_pointrange(data=dist_district_lineage_rr_table, mapping=aes(x=as.factor(distance_window), y=RR_mean, ymin = RR_lower_ci, ymax = RR_upper_ci), alpha=0.9) +
  geom_hline(yintercept=1, linetype="dashed", color = "red") + ##dashed line at 1
  #scale_y_log10()+
  ylab("Risk Ratio")+
  xlab("Distance window (km)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_rect(xmin = 14.5, xmax = 15.5, ymin = -1,  ymax = 5,
            fill = 'blue', alpha = 0.01) +
  geom_rect(xmin = 0, xmax = 15.5, ymin = -1, ymax = 5,
            fill = 'red', alpha = 0.01)

rolling_distance_RR_plot_50km
```
![RR vs. Rolling Distance | at Lineage level](https://github.com/hsuehchien66/rrspread_v2/blob/main/figures/rolling_distance_RR_allref_main.png)


## Rolling divergence vs. relative ratio | Within district vs. Between districts

```r
des_divtime_min1 = seq(from=0, to=8, by=2)
des_divtime_min2 = seq(from=10, to=15, by=5)
des_divtime_min = c(des_divtime_min1, des_divtime_min2, 20)
des_divtime_max1 = seq(from=2, to=10, by=2)
des_divtime_max2 = seq(from=15, to=20, by=5)
des_divtime_max = c(des_divtime_max1, des_divtime_max2, 30)
des_divtime_range = length(des_divtime_max)

res_list = c()
res_df_divtime = data.frame()

colyear_mat_ones = matrix(1, nrow=dim(selected_gpsc_metadata)[1], ncol=dim(selected_gpsc_metadata)[1])
diag(colyear_mat_ones) <- NA

RRplotLists = list()
RRdivtimeLists = list()
RRdframLists = list()

for(divtimes in 1:des_divtime_range){
  divtime_name = paste0("Divergence: ",des_divtime_min[divtimes],"-",des_divtime_max[divtimes],"yrs")
  loc_divtime_mat = RR_matrix(meta_table = selected_gpsc_metadata, bactdate_tree = multitrees,
                              des_divtime_min = des_divtime_min[divtimes],
                              des_divtime_max = des_divtime_max[divtimes])
  loc_divtime_mat_rr = RR_bootstrap(metadataset = selected_gpsc_metadata,
                                    des_geo_matrix = loc_divtime_mat$des_location_mat,
                                    ref_geo_matrix = loc_divtime_mat$ref_location_mat,
                                    des_strain_matrix = loc_divtime_mat$des_divtime_mat,
                                    ref_strain_matrix = loc_divtime_mat$ref_divtime_mat,
                                    colyear_matrix = colyear_mat_ones,
                                    nboot = 100,
                                    subsample_min = 70)
  rr_plot = PltRR(loc_divtime_mat_rr, divtime_name, divtime_name)
  RRplotLists[[divtimes]] = print(rr_plot)
  
  divtimes_col = data.frame(divtime_min=des_divtime_min[divtimes], divtime_max=des_divtime_max[divtimes], divtime_mid = (des_divtime_min[divtimes]+des_divtime_max[divtimes])/2)
  new_row = cbind(loc_divtime_mat_rr, divtimes_col)
  res_df_divtime = rbind(res_df_divtime, new_row)
  
  assign(divtime_name, rr_plot)
  res_list=c(res_list, divtime_name)
}

rolling_distance_divtime <- PltRR_rolling_divtime(res_df_divtime,"rolling divergence")+
  scale_x_discrete(labels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10-15", "15-20", "20-30"))+
  #scale_y_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), labels = c("0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "4.0", "4.5", "5.0"))+
  scale_y_log10(limits=c(0.1,7))+
  ggtitle("By District")+
  theme(axis.text = element_text(size=18),title = element_text(size=20))
```

## Rolling divergence vs. relative ratio | Close area (0-30 km) vs. Distant area

```r
des_divtime_min1 = seq(from=0, to=8, by=2)
des_divtime_min2 = seq(from=10, to=15, by=5)
des_divtime_min = c(des_divtime_min1, des_divtime_min2, 20)
des_divtime_max1 = seq(from=2, to=10, by=2)
des_divtime_max2 = seq(from=15, to=20, by=5)
des_divtime_max = c(des_divtime_max1, des_divtime_max2, 30)
des_divtime_range = length(des_divtime_max)

res_list = c()
res_df = data.frame()

colyear_mat_ones = matrix(1, nrow=dim(selected_gpsc_metadata)[1], ncol=dim(selected_gpsc_metadata)[1])
diag(colyear_mat_ones) <- NA

RRplotLists = list()
RRdivtimeLists = list()
RRdframLists = list()

for(divtimes in 1:des_divtime_range){
  divtime_name = paste0("Divergence: ",des_divtime_min[divtimes],"~",des_divtime_max[divtimes],"yrs")
  loc_divtime_mat = RR_matrix(meta_table = selected_gpsc_metadata, bactdate_tree = multitrees,
                              des_divtime_min = des_divtime_min[divtimes],
                              des_divtime_max = des_divtime_max[divtimes],
                              des_distance_min = 0,
                              des_distance_max = 30, 
                              ref_distance_min = 150,
                              ref_distance_max = 300)
  loc_divtime_mat_rr = RR_bootstrap(metadataset = selected_gpsc_metadata,
                                    des_geo_matrix = loc_divtime_mat$des_distance_mat,
                                    ref_geo_matrix = loc_divtime_mat$ref_distance_mat,
                                    des_strain_matrix = loc_divtime_mat$des_divtime_mat,
                                    ref_strain_matrix = loc_divtime_mat$ref_divtime_mat,
                                    colyear_matrix = colyear_mat_ones,
                                    nboot = 100,
                                    subsample_min = 70)
  rr_plot = PltRR(loc_divtime_mat_rr, divtime_name, divtime_name)
  RRplotLists[[divtimes]] = print(rr_plot)
  
  divtimes_col = data.frame(divtime_min=des_divtime_min[divtimes], divtime_max=des_divtime_max[divtimes], divtime_mid = (des_divtime_min[divtimes]+des_divtime_max[divtimes])/2)
  new_row = cbind(loc_divtime_mat_rr, divtimes_col)
  res_df = rbind(res_df, new_row)
  
  assign(divtime_name, rr_plot)
  res_list=c(res_list, divtime_name)
}


rolling_distance_divtime_30km_150kmRef <- PltRR_rolling_divtime(res_df,"rolling divergence")+
  scale_x_discrete(labels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10-15", "15-20", "20-30"))+
  #scale_y_continuous(breaks = c(0.5, 1, 1.5, 2, 3, 4, 5, 6, 7), labels = c("0.5", "1.0", "1.5", "2.0", "3.0",  "4.0", "5.0", "6.0", "7.0"))+
  scale_y_log10(limits=c(0.1,60))+
  ggtitle("By Distance")+
  theme(axis.text = element_text(size=18),title = element_text(size=20))
```

![RR vs. Rolling Divergence](https://github.com/hsuehchien66/rrspread_v2/blob/main/figures/figure4.png)

