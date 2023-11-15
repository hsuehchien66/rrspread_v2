library(ape)
library(BactDating)
library(dplyr)
library(geosphere)
library(ggplot2)
library(data.table)
library(cowplot)


metadata_file_path = "/Users/hc14/Documents/SpneumoIsrael/Israel_raw_data/Israel_GPSC_city_region.csv"
metadata = read.csv(metadata_file_path)
### generate lane_id files for all 1,174 isolates
lane_id_1174 = metadata[,"lane_id"]
write.table(lane_id_1174, "Israel_1174_isolates_laneid.txt", row.names = F, col.names = F, quote = F)


## GPSC8, GPSC55, GPSC47, GPSC6, GPSC10
#read tree files
{
  gpsc8_gubbin_tree_path = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/05_topGPSC_gubbin_output_tree/GPSC8_reference_bwa.final_tree.tre"
  gpsc8_gubbin_tree = read.tree(gpsc8_gubbin_tree_path)

  gpsc55_gubbin_tree_path = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/05_topGPSC_gubbin_output_tree/GPSC55_reference_bwa.final_tree.tre"
  gpsc55_gubbin_tree = read.tree(gpsc55_gubbin_tree_path)

  gpsc47_gubbin_tree_path = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/05_topGPSC_gubbin_output_tree/GPSC47_reference_bwa.final_tree.tre"
  gpsc47_gubbin_tree = read.tree(gpsc47_gubbin_tree_path)

  gpsc6_gubbin_tree_path = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/05_topGPSC_gubbin_output_tree/GPSC6_reference_bwa.final_tree.tre"
  gpsc6_gubbin_tree = read.tree(gpsc6_gubbin_tree_path)

  gpsc10_gubbin_tree_path = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/05_topGPSC_gubbin_output_tree/GPSC10_reference_bwa.final_tree.tre"
  gpsc10_gubbin_tree = read.tree(gpsc10_gubbin_tree_path)
}

#generate bactdating results
{
  gpsc8_gubbin_tree_plot = plt_gubbin_tree(gpsc8_gubbin_tree)
  gpsc8_time_tree_list = match_tree_meta(metadata, gpsc8_gubbin_tree, "GPSC8_reference", drop_tip = TRUE)
  gpsc8_bd_res = time_resolved(gpsc8_time_tree_list$match_tree_meta, gpsc8_time_tree_list$tree, "/Users/hc14/Documents/SpneumoIsrael/analysis_results/bd_res_GPSC8", time="Year_collection", n_it = 300000, branch_model="mixedgamma")

  gpsc55_gubbin_tree_plot = plt_gubbin_tree(gpsc55_gubbin_tree)
  gpsc55_time_tree_list = match_tree_meta(metadata, gpsc55_gubbin_tree, "GPSC55_reference", drop_tip = TRUE)
  gpsc55_bd_res = time_resolved(gpsc55_time_tree_list$match_tree_meta, gpsc55_time_tree_list$tree, "/Users/hc14/Documents/SpneumoIsrael/analysis_results/bd_res_GPSC55", n_it = 300000, branch_model="mixedgamma")

  gpsc47_gubbin_tree_plot = plt_gubbin_tree(gpsc47_gubbin_tree)
  gpsc47_time_tree_list = match_tree_meta(metadata, gpsc47_gubbin_tree, "GPSC47_reference", drop_tip = TRUE)
  gpsc47_bd_res = time_resolved(gpsc47_time_tree_list$match_tree_meta, gpsc47_time_tree_list$tree, "/Users/hc14/Documents/SpneumoIsrael/analysis_results/bd_res_GPSC47", n_it = 300000, branch_model="mixedgamma")

  gpsc6_gubbin_tree_plot = plt_gubbin_tree(gpsc6_gubbin_tree)
  gpsc6_time_tree_list = match_tree_meta(metadata, gpsc6_gubbin_tree, "GPSC6_reference", drop_tip = TRUE)
  gpsc6_bd_res = time_resolved(gpsc6_time_tree_list$match_tree_meta, gpsc6_time_tree_list$tree, "/Users/hc14/Documents/SpneumoIsrael/analysis_results/bd_res_GPSC6", n_it = 300000, branch_model="mixedgamma")

  gpsc10_gubbin_tree_plot = plt_gubbin_tree(gpsc10_gubbin_tree)
  gpsc10_time_tree_list = match_tree_meta(metadata, gpsc10_gubbin_tree, "GPSC10_reference", drop_tip = TRUE)
  gpsc10_bd_res = time_resolved(gpsc10_time_tree_list$match_tree_meta, gpsc10_time_tree_list$tree, "/Users/hc14/Documents/SpneumoIsrael/analysis_results/bd_res_GPSC10", n_it = 300000, branch_model="mixedgamma")
}

### check mcmc convergence
mcmc_gpsc6=as.mcmc.resBactDating(gpsc6_bd_res$BactDating) #check whether converge
coda::effectiveSize(mcmc_gpsc6)
mcmc_gpsc8=as.mcmc.resBactDating(gpsc8_bd_res$BactDating) #check whether converge
coda::effectiveSize(mcmc_gpsc8)
mcmc_gpsc10=as.mcmc.resBactDating(gpsc10_bd_res$BactDating) #check whether converge
coda::effectiveSize(mcmc_gpsc10)
mcmc_gpsc47=as.mcmc.resBactDating(gpsc47_bd_res$BactDating) #check whether converge
coda::effectiveSize(mcmc_gpsc47)
mcmc_gpsc55=as.mcmc.resBactDating(gpsc55_bd_res$BactDating) #check whether converge
coda::effectiveSize(mcmc_gpsc55)

##bactdate tree
{
  gpsc08_bd_tree = gpsc8_bd_res$BactDating$tree
  gpsc55_bd_tree = gpsc55_bd_res$BactDating$tree
  gpsc47_bd_tree = gpsc47_bd_res$BactDating$tree
  gpsc06_bd_tree = gpsc6_bd_res$BactDating$tree
  gpsc10_bd_tree = gpsc10_bd_res$BactDating$tree
}


#bind trees and combined with metadata
{
  multitrees = bind_multitree(gpsc08_bd_tree, gpsc55_bd_tree, gpsc47_bd_tree, gpsc06_bd_tree, gpsc10_bd_tree)
  plot(multitrees, cex=0.3, no.margin = TRUE)
  write.tree(multitrees, file = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/GPSC8_55_47_6_10_Israel.nwk")
  lanes <- multitrees$tip.label

  ## subset only those present in the trees
  ####metadata = read.csv(metadata_file)
  selected_gpsc_metadata = subset(metadata, metadata$lane_id %in% lanes)
  lanes_gpsc <- subset(lanes, lanes  %in% selected_gpsc_metadata$lane_id)
  multitrees_overall <- keep.tip(multitrees, lanes_gpsc)

  #Reorder tree based on order of tip labels
  order <- multitrees_overall$tip.label
  selected_gpsc_metadata <- selected_gpsc_metadata %>%
    dplyr::slice(match(order, lane_id))
}
write.table(selected_gpsc_metadata, file = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/GPSC8_55_47_6_10_Israel_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(selected_gpsc_metadata, file = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/GPSC8_55_47_6_10_Israel_metadata.csv", quote = FALSE, row.names = FALSE)

gpsc6_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "6"),]
gpsc8_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "8"),]
gpsc10_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "10"),]
gpsc47_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "47"),]
gpsc55_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "55"),]

gpsc6_rrmatrix <- RR_matrix(meta_table = gpsc6_metadata, bactdate_tree = gpsc6_bd_res$BactDating$tree)
gpsc8_rrmatrix <- RR_matrix(meta_table = gpsc8_metadata, bactdate_tree = gpsc8_bd_res$BactDating$tree)
gpsc10_rrmatrix <- RR_matrix(meta_table = gpsc10_metadata, bactdate_tree = gpsc10_bd_res$BactDating$tree)
gpsc47_rrmatrix <- RR_matrix(meta_table = gpsc47_metadata, bactdate_tree = gpsc47_bd_res$BactDating$tree)
gpsc55_rrmatrix <- RR_matrix(meta_table = gpsc55_metadata, bactdate_tree = gpsc55_bd_res$BactDating$tree)

gpsc6_distance_dist <- as.data.frame(to.upper(gpsc6_rrmatrix$distance_num_mat))
gpsc6_divtime_dist <- as.data.frame(to.upper(gpsc6_rrmatrix$divtime_num_mat))
gpsc6_divtime_distance <- cbind(gpsc6_divtime_dist, gpsc6_distance_dist)
colnames(gpsc6_divtime_distance) <- c("divtime", "distance")
png("/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC6_distance_divtime_dots.png")
plot(gpsc6_divtime_distance$divtime, gpsc6_divtime_distance$distance)
dev.off()

gpsc8_distance_dist <- as.data.frame(to.upper(gpsc8_rrmatrix$distance_num_mat))
gpsc8_divtime_dist <- as.data.frame(to.upper(gpsc8_rrmatrix$divtime_num_mat))
gpsc8_divtime_distance <- cbind(gpsc8_divtime_dist, gpsc8_distance_dist)
colnames(gpsc8_divtime_distance) <- c("divtime", "distance")
png("/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC8_distance_divtime_dots.png")
plot(gpsc8_divtime_distance$divtime, gpsc8_divtime_distance$distance)
dev.off()

gpsc10_distance_dist <- as.data.frame(to.upper(gpsc10_rrmatrix$distance_num_mat))
gpsc10_divtime_dist <- as.data.frame(to.upper(gpsc10_rrmatrix$divtime_num_mat))
gpsc10_divtime_distance <- cbind(gpsc10_divtime_dist, gpsc10_distance_dist)
colnames(gpsc10_divtime_distance) <- c("divtime", "distance")
png("/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC10_distance_divtime_dots.png")
plot(gpsc10_divtime_distance$divtime, gpsc10_divtime_distance$distance)
dev.off()

gpsc47_distance_dist <- as.data.frame(to.upper(gpsc47_rrmatrix$distance_num_mat))
gpsc47_divtime_dist <- as.data.frame(to.upper(gpsc47_rrmatrix$divtime_num_mat))
gpsc47_divtime_distance <- cbind(gpsc47_divtime_dist, gpsc47_distance_dist)
colnames(gpsc47_divtime_distance) <- c("divtime", "distance")
png("/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC47_distance_divtime_dots.png")
plot(gpsc47_divtime_distance$divtime, gpsc47_divtime_distance$distance)
dev.off()

gpsc55_distance_dist <- as.data.frame(to.upper(gpsc55_rrmatrix$distance_num_mat))
gpsc55_divtime_dist <- as.data.frame(to.upper(gpsc55_rrmatrix$divtime_num_mat))
gpsc55_divtime_distance <- cbind(gpsc55_divtime_dist, gpsc55_distance_dist)
colnames(gpsc55_divtime_distance) <- c("divtime", "distance")
png("/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC55_distance_divtime_dots.png")
plot(gpsc55_divtime_distance$divtime, gpsc55_divtime_distance$distance)
dev.off()

gpsc_all_divtime_distance <- rbind(gpsc6_divtime_distance, gpsc8_divtime_distance, gpsc10_divtime_distance,
                                   gpsc47_divtime_distance, gpsc55_divtime_distance)
png("/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC_all_divtime_dots.png")
plot(gpsc_all_divtime_distance$divtime, gpsc_all_divtime_distance$distance)
dev.off()

#GPSC6
{
  nboot=20
  x_divtime_max = seq(from=5, to=45, by=1)
  x_divtime_min = rep(0, length(x_divtime_max))
  #x_divtime_min = seq(from=0, to=40, by=1)

  boot_select <- data.frame()
  for (win in 1:length(x_divtime_max)) {
    boot.mean <- c()
    for (n in 1:nboot) {
      a = which(gpsc6_divtime_distance$divtime > x_divtime_min[win] &
                  gpsc6_divtime_distance$divtime < x_divtime_max[win])        #select a certain region (index)
      if (length(a)!=0){
        b = sample(x=a, size=100, replace = T)  #subsample without replacement (index)
        resample.mean <- mean(gpsc6_divtime_distance[b, ]$distance)
        boot.mean <- c(boot.mean, resample.mean)
      }
    }
    distance.med <- quantile(boot.mean, probs = c(0.5), na.rm=T) #calculate confidence interval
    distance.mean <- mean(boot.mean)
    distance.upperci <- quantile(boot.mean, probs = c(0.975), na.rm=T)
    distance.lowerci <- quantile(boot.mean, probs = c(0.025), na.rm=T)
    mid_win <- mean(c(x_divtime_min[win], x_divtime_max[win]))
    window_divtime_distance <- cbind(x_divtime_min[win], x_divtime_max[win], mid_win, distance.med, distance.mean, distance.upperci, distance.lowerci)
    boot_select <- rbind(boot_select, window_divtime_distance)
  }

  gpsc6_dist_div <- ggplot(boot_select, aes(x=mid_win, y=distance.mean, ymin=distance.lowerci, ymax=distance.upperci))+
    geom_line()+
    geom_ribbon(fill="pink", alpha=0.5)+
    #xlim(0,60)+
    xlab("Divergence Time \n (years)")+
    ylab("Distance (km)")+
    ylim(c(0,87))+
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

  #png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/gpsc6_distance_divtime_CI.png")
  png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/gpsc6_distance_divtime_rollingCI.png")
  gpsc6_dist_div
  dev.off()
}

#GPSC8
{
  nboot=20
  x_divtime_max = seq(from=5, to=80, by=1)
  #x_divtime_min = rep(0, length(x_divtime_max))
  x_divtime_min = seq(from=0, to=75, by=1)

  boot_select <- data.frame()
  for (win in 1:length(x_divtime_max)) {
    boot.mean <- c()
    for (n in 1:nboot) {
      a = which(gpsc8_divtime_distance$divtime > x_divtime_min[win] &
                  gpsc8_divtime_distance$divtime < x_divtime_max[win])        #select a certain region (index)
      if (length(a)!=0){
        b = sample(x=a, size=100, replace = T)  #subsample without replacement (index)
        resample.mean <- mean(gpsc8_divtime_distance[b, ]$distance)
        boot.mean <- c(boot.mean, resample.mean)
      }
    }
    distance.med <- quantile(boot.mean, probs = c(0.5), na.rm=T) #calculate confidence interval
    distance.mean <- mean(boot.mean)
    distance.upperci <- quantile(boot.mean, probs = c(0.975), na.rm=T)
    distance.lowerci <- quantile(boot.mean, probs = c(0.025), na.rm=T)
    mid_win <- mean(c(x_divtime_min[win], x_divtime_max[win]))
    window_divtime_distance <- cbind(x_divtime_min[win], x_divtime_max[win], mid_win, distance.med, distance.mean, distance.upperci, distance.lowerci)
    boot_select <- rbind(boot_select, window_divtime_distance)
  }
  boot_select <- na.omit(boot_select)

  gpsc8_dist_div <- ggplot(boot_select, aes(x=mid_win, y=distance.mean, ymin=distance.lowerci, ymax=distance.upperci))+
    geom_line()+
    geom_ribbon(fill="pink", alpha=0.5)+
    #xlim(0,60)+
    xlab("Divergence Time \n (years)")+
    ylab("Distance (km)")+
    #ylim(c(0,87))+
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

  png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/gpsc8_distance_divtime_CI.png")
  #png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/gpsc8_distance_divtime_rollingCI.png")
  gpsc8_dist_div
  dev.off()
}

#GPSC10
{
  nboot=20
  x_divtime_max = seq(from=10, to=370, by=1)
  #x_divtime_min = rep(0, length(x_divtime_max))
  x_divtime_min = seq(from=0, to=360, by=1)

  boot_select <- data.frame()
  for (win in 1:length(x_divtime_max)) {
    boot.mean <- c()
    for (n in 1:nboot) {
      a = which(gpsc10_divtime_distance$divtime > x_divtime_min[win] &
                  gpsc10_divtime_distance$divtime < x_divtime_max[win])        #select a certain region (index)
      if (length(a)!=0){
        b = sample(x=a, size=100, replace = T)  #subsample without replacement (index)
        resample.mean <- mean(gpsc10_divtime_distance[b, ]$distance)
        boot.mean <- c(boot.mean, resample.mean)
      }
    }
    distance.med <- quantile(boot.mean, probs = c(0.5), na.rm=T) #calculate confidence interval
    distance.mean <- mean(boot.mean)
    distance.upperci <- quantile(boot.mean, probs = c(0.975), na.rm=T)
    distance.lowerci <- quantile(boot.mean, probs = c(0.025), na.rm=T)
    mid_win <- mean(c(x_divtime_min[win], x_divtime_max[win]))
    window_divtime_distance <- cbind(x_divtime_min[win], x_divtime_max[win], mid_win, distance.med, distance.mean, distance.upperci, distance.lowerci)
    boot_select <- rbind(boot_select, window_divtime_distance)
  }
  boot_select <- na.omit(boot_select)

  gpsc10_dist_div <- ggplot(boot_select, aes(x=mid_win, y=distance.mean, ymin=distance.lowerci, ymax=distance.upperci))+
    geom_line()+
    geom_ribbon(fill="pink", alpha=0.5)+
    #xlim(0,60)+
    xlab("Divergence Time \n (years)")+
    ylab("Distance (km)")+
    #ylim(c(0,87))+
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

  png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/gpsc10_distance_divtime_CI.png")
  #png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/gpsc10_distance_divtime_rollingCI.png")
  gpsc10_dist_div
  dev.off()
}

#GPSC47
{
  nboot=20
  x_divtime_max = seq(from=70, to=80, by=1)
  #x_divtime_min = rep(0, length(x_divtime_max))
  x_divtime_min = seq(from=60, to=70, by=1)

  boot_select <- data.frame()
  for (win in 1:length(x_divtime_max)) {
    boot.mean <- c()
    for (n in 1:nboot) {
      a = which(gpsc47_divtime_distance$divtime > x_divtime_min[win] &
                  gpsc47_divtime_distance$divtime < x_divtime_max[win])        #select a certain region (index)
      if (length(a)!=0){
        b = sample(x=a, size=100, replace = T)  #subsample without replacement (index)
        resample.mean <- mean(gpsc47_divtime_distance[b, ]$distance)
        boot.mean <- c(boot.mean, resample.mean)
      }
    }
    distance.med <- quantile(boot.mean, probs = c(0.5), na.rm=T) #calculate confidence interval
    distance.mean <- mean(boot.mean)
    distance.upperci <- quantile(boot.mean, probs = c(0.975), na.rm=T)
    distance.lowerci <- quantile(boot.mean, probs = c(0.025), na.rm=T)
    mid_win <- mean(c(x_divtime_min[win], x_divtime_max[win]))
    window_divtime_distance <- cbind(x_divtime_min[win], x_divtime_max[win], mid_win, distance.med, distance.mean, distance.upperci, distance.lowerci)
    boot_select <- rbind(boot_select, window_divtime_distance)
  }
  boot_select <- na.omit(boot_select)

  gpsc47_dist_div <- ggplot(boot_select, aes(x=mid_win, y=distance.mean, ymin=distance.lowerci, ymax=distance.upperci))+
    geom_line()+
    geom_ribbon(fill="pink", alpha=0.5)+
    #xlim(0,60)+
    xlab("Divergence Time \n (years)")+
    ylab("Distance (km)")+
    #ylim(c(0,87))+
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

  png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/gpsc47_distance_divtime_60_80_CI.png")
  #png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/gpsc47_distance_divtime_rollingCI.png")
  gpsc47_dist_div
  dev.off()
}

#GPSC55
{
nboot=20
x_divtime_max = seq(from=5, to=25, by=1)
x_divtime_min = rep(0, length(x_divtime_max))
#x_divtime_min = seq(from=0, to=20, by=1)

boot_select <- data.frame()
for (win in 1:length(x_divtime_max)) {
  boot.mean <- c()
  for (n in 1:nboot) {
    a = which(gpsc55_divtime_distance$divtime > x_divtime_min[win] &
                gpsc55_divtime_distance$divtime < x_divtime_max[win])        #select a certain region (index)
    if (length(a)!=0){
      b = sample(x=a, size=100, replace = T)  #subsample without replacement (index)
      resample.mean <- mean(gpsc55_divtime_distance[b, ]$distance)
      boot.mean <- c(boot.mean, resample.mean)
    }
  }
  distance.med <- quantile(boot.mean, probs = c(0.5), na.rm=T) #calculate confidence interval
  distance.mean <- mean(boot.mean)
  distance.upperci <- quantile(boot.mean, probs = c(0.975), na.rm=T)
  distance.lowerci <- quantile(boot.mean, probs = c(0.025), na.rm=T)
  mid_win <- mean(c(x_divtime_min[win], x_divtime_max[win]))
  window_divtime_distance <- cbind(x_divtime_min[win], x_divtime_max[win], mid_win, distance.med, distance.mean, distance.upperci, distance.lowerci)
  boot_select <- rbind(boot_select, window_divtime_distance)
}

gpsc55_dist_div <- ggplot(boot_select, aes(x=mid_win, y=distance.mean, ymin=distance.lowerci, ymax=distance.upperci))+
  geom_line()+
  geom_ribbon(fill="pink", alpha=0.5)+
  #xlim(0,60)+
  xlab("Divergence Time \n (years)")+
  ylab("Distance (km)")+
  ylim(c(0,165))+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

#png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC55_distance_divtime_CI.png")
png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC55_distance_divtime_rollingCI.png")
gpsc55_dist_div
dev.off()
}

#GPSCall
{
  gpsc_all_divtime_distance <- rbind(gpsc6_divtime_distance, gpsc8_divtime_distance, gpsc10_divtime_distance,
                                     gpsc47_divtime_distance, gpsc55_divtime_distance)


  nboot=20
  x_divtime_max = seq(from=10, to=90, by=1)
  #x_divtime_min = rep(0, length(x_divtime_max))
  x_divtime_min = seq(from=0, to=80, by=1)

  boot_select <- data.frame()
  for (win in 1:length(x_divtime_max)) {
    boot.mean <- c()
    for (n in 1:nboot) {
      a = which(gpsc_all_divtime_distance$divtime > x_divtime_min[win] &
                  gpsc_all_divtime_distance$divtime < x_divtime_max[win])        #select a certain region (index)
      if (length(a)!=0){
        b = sample(x=a, size=100, replace = T)  #subsample without replacement (index)
        resample.mean <- mean(gpsc_all_divtime_distance[b, ]$distance)
        boot.mean <- c(boot.mean, resample.mean)
      }
    }
    distance.med <- quantile(boot.mean, probs = c(0.5), na.rm=T) #calculate confidence interval
    distance.mean <- mean(boot.mean)
    distance.upperci <- quantile(boot.mean, probs = c(0.975), na.rm=T)
    distance.lowerci <- quantile(boot.mean, probs = c(0.025), na.rm=T)
    mid_win <- mean(c(x_divtime_min[win], x_divtime_max[win]))
    window_divtime_distance <- cbind(x_divtime_min[win], x_divtime_max[win], mid_win, distance.med, distance.mean, distance.upperci, distance.lowerci)
    boot_select <- rbind(boot_select, window_divtime_distance)
  }
  boot_select <- na.omit(boot_select)

  gpsc_all_dist_div <- ggplot(boot_select, aes(x=mid_win, y=distance.mean, ymin=distance.lowerci, ymax=distance.upperci))+
    geom_line()+
    geom_ribbon(fill="pink", alpha=0.5)+
    #xlim(0,60)+
    xlab("Divergence Time \n (years)")+
    ylab("Distance (km)")+
    #ylim(c(0,87))+
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

  png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/gpsc_all_distance_divtime_0_90_CI.png")
  #png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/gpsc_all_distance_divtime_0_90_cumCI.png")
  gpsc_all_dist_div
  dev.off()
}

## lineage v.s. location
## load library: geosphere
rr_lineage_location_mat = RR_matrix(selected_gpsc_metadata, multitrees)

## visualise distance and divergence time distribution
distance_dist <- as.data.frame(to.upper(rr_lineage_location_mat$distance_num_mat))
colnames(distance_dist) = "distance"
divtime_dist <- as.data.frame(to.upper(rr_lineage_location_mat$divtime_num_mat))
colnames(divtime_dist) = "divtime"

distance_dist_p <- distance_dist %>%
  ggplot( aes(x=distance)) +
  geom_histogram( binwidth=10, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Bin size = 10 km") +
  theme(
    plot.title = element_text(size=15)
  )

divtime_dist_p <- divtime_dist %>%
  ggplot( aes(x=divtime)) +
  geom_histogram( binwidth=5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Bin size = 5 years") +
  theme(plot.title = element_text(size=15))+
  theme_bw()

### divergence time vs. distance
div_dist = cbind(distance_dist, divtime_dist)
attach(div_dist)
div_dist_sorted = div_dist[order(divtime),]

dist_div_p = ggplot(div_dist_sorted, aes(x=divtime, y=distance))+
  geom_point()+
  #xlim(0,100)+
  xlab("Divergence Time \n (years)")+
  ylab("Mean Distance (km)")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))


cum_mean_dist = rep(NA, nrow(div_dist_sorted))
for (i in 1:nrow(div_dist_sorted)){
  cum_mean_dist[i] = mean(div_dist_sorted[1:i,"distance"])
}

cum_mean_dist_div = cbind(div_dist_sorted, cum_mean_dist)
cum_mean_dist_div_p = ggplot(cum_mean_dist_div, aes(x=divtime, y=cum_mean_dist))+
  geom_point()+
  #xlim(0,60)+
  xlab("Cumulative Divergence Time \n (years)")+
  ylab("Mean Distance (km)")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

#x_divtime_min = seq(from=0, to=130, by=10)
x_divtime_max = c(0.25, 0.3,0.35,0.4, 0.45, 0.5,0.75, 1,1.25 ,1.5 , 1.75,2,2.5,3,3.5, 4, 6, 10, 20, 30, 40, 50, 60, 80, 90, 100)
x_divtime_max = seq(from=0.25, to=40, by=0.25)

x_divtime_min = rep(0, length(x_divtime_max))
#x_divtime_min = c(0, 0.5, 1 , 1.5, 2, 4, 6, 10, 20, 30, 40, 50, 60, 80, 90)
nboot = 200
res_dist_divtime = data.frame()
for (tim in 1:length(x_divtime_min)) {
  n = which(div_dist_sorted$divtime < x_divtime_max[tim])[length(which(div_dist_sorted$divtime <  x_divtime_max[tim]))]
  resample <- sample(1:n, nboot, replace = TRUE)
  ci.dist.mean <- mean(div_dist_sorted$distance[resample])
  print(ci.dist.med)
  ci.dist.upperci <- quantile(div_dist_sorted$distance[resample], probs = c(0.9))
  ci.dist.lowerci <- quantile(div_dist_sorted$distance[resample], probs = c(0.1))
  res_dist_divtime <- rbind(res_dist_divtime,as.data.frame(cbind(x_divtime_min[tim], x_divtime_max[tim], ci.dist.mean, ci.dist.lowerci, ci.dist.upperci)))
}
## calculate mean first and bootstrap
colnames(res_dist_divtime)[c(1,2)] <- c("divtime_min", "divtime_max")

############################################################
### Within district / between district: lineage
## load library: data.table
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


############################################################
### Short distance / Long distance: Lineage
rr_lineage_distance_0_20km_mat = RR_matrix(selected_gpsc_metadata, multitrees,
                                    des_distance_min = 0, des_distance_max = 20)
rr_lineage_distance_0_20km = RR_bootstrap(metadataset = selected_gpsc_metadata,
                                   des_geo_matrix = rr_lineage_distance_0_20km_mat$des_distance_mat,
                                   ref_geo_matrix = rr_lineage_distance_0_20km_mat$ref_distance_mat,
                                   des_strain_matrix = rr_lineage_distance_0_20km_mat$des_lineage_mat,
                                   ref_strain_matrix = rr_lineage_distance_0_20km_mat$ref_lineage_mat,
                                   colyear_matrix = colyear_mat_ones, nboot = 300, subsample_min = 70)

distance_0_20km_lineage_rr = PltRR(rr_lineage_distance_0_20km, "within 20 km Relative Ratio", "0-20 km")
#
rr_lineage_distance_30_40km_mat = RR_matrix(selected_gpsc_metadata, multitrees,
                                    des_distance_min = 30, des_distance_max = 40)
rr_lineage_distance_30_40km = RR_bootstrap(metadataset = selected_gpsc_metadata,
                                   des_geo_matrix = rr_lineage_distance_30_40km_mat$des_distance_mat,
                                   ref_geo_matrix = rr_lineage_distance_30_40km_mat$ref_distance_mat,
                                   des_strain_matrix = rr_lineage_distance_30_40km_mat$des_lineage_mat,
                                   ref_strain_matrix = rr_lineage_distance_30_40km_mat$ref_lineage_mat,
                                   colyear_matrix = colyear_mat_ones, nboot = 1000, subsample_min = 70)

distance_30_40km_lineage_rr = PltRR(rr_lineage_distance_30_40km, "30-40 km Relative Ratio", "30-40")
#
rr_lineage_distance_40_50km_mat = RR_matrix(selected_gpsc_metadata, multitrees,
                                            des_distance_min = 40, des_distance_max = 50)
rr_lineage_distance_40_50km = RR_bootstrap(metadataset = selected_gpsc_metadata,
                                           des_geo_matrix = rr_lineage_distance_40_50km_mat$des_distance_mat,
                                           ref_geo_matrix = rr_lineage_distance_40_50km_mat$ref_distance_mat,
                                           des_strain_matrix = rr_lineage_distance_40_50km_mat$des_lineage_mat,
                                           ref_strain_matrix = rr_lineage_distance_40_50km_mat$ref_lineage_mat,
                                           colyear_matrix = colyear_mat_ones, nboot = 1000, subsample_min = 70)

distance_40_50km_lineage_rr = PltRR(rr_lineage_distance_40_50km, "40-50 km Relative Ratio", "40-50")
#
rr_lineage_distance_80_150km_mat = RR_matrix(selected_gpsc_metadata, multitrees,
                                    des_distance_min = 80, des_distance_max = 150)
rr_lineage_distance_80_150km = RR_bootstrap(metadataset = selected_gpsc_metadata,
                                   des_geo_matrix = rr_lineage_distance_80_150km_mat$des_distance_mat,
                                   ref_geo_matrix = rr_lineage_distance_80_150km_mat$ref_distance_mat,
                                   des_strain_matrix = rr_lineage_distance_80_150km_mat$des_lineage_mat,
                                   ref_strain_matrix = rr_lineage_distance_80_150km_mat$ref_lineage_mat,
                                   colyear_matrix = colyear_mat_ones, nboot = 1000, subsample_min = 70)

distance_80_150km_lineage_rr = PltRR(rr_lineage_distance_80_150km, "80-150 km Relative Ratio", "80-150")
#

dist_lineage_rr_table <- rbind(rr_lineage_distance_0_20km,rr_lineage_distance_30_40km, rr_lineage_distance_40_50km, rr_lineage_distance_80_150km, rr)
dist_lineage_rr_table[,"distance_mid"] = c("0-20", "30-40","40-50", ">80", "within district")
dist_lineage_rr_table$distance_mid <- factor(dist_lineage_rr_table$distance_mid, levels = c("0-20", "30-40","40-50", ">80", "within district"))
rolling_distance <- PltRR_rolling(dist_lineage_rr_table, "rolling distance")


p <- ggplot(dist_lineage_rr_table, aes(x=as.factor(distance_mid), y=RR_mean))+
  geom_pointrange(data=dist_lineage_rr_table, mapping=aes(x=as.factor(distance_mid), y=RR_mean, ymin = RR_lower_ci, ymax = RR_upper_ci), alpha=0.9) +
  geom_hline(yintercept=1, linetype="dashed", color = "red") + ##dashed line at 1
  scale_y_log10()+
  ylab("Risk Ratio")+
  xlab("Distance window (km)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_rect(xmin = 4.5, xmax = 6, ymin = -0.5, ymax = 1.5,
            fill = 'blue', alpha = 0.05) +
  geom_rect(xmin = 0, xmax = 4.5, ymin = -0.5, ymax = 1.5,
            fill = 'red', alpha = 0.05)


############################################################
### rolling distance vs. relative ratio, lineage
des_distance_min = seq(from=0, to=75, by=15)
des_distance_max = seq(from=20, to=95, by=15)
des_distance_range = length(des_distance_max)

res_list = c()
res_df = data.frame()

colyear_mat_ones = matrix(1, nrow=dim(selected_gpsc_metadata)[1], ncol=dim(selected_gpsc_metadata)[1])
diag(colyear_mat_ones) <- NA

RRplotLists = list()
RRdistanceLists = list()
RRdframLists = list()

for(distances in 1:des_distance_range){
  distance_name = paste0("Distance: ",des_distance_min[distances],"~",des_distance_max[distances],"km")
  dist_lineage_mat = RR_matrix(meta_table = selected_gpsc_metadata, bactdate_tree = multitrees,
                               des_distance_min = des_distance_min[distances],
                               des_distance_max = des_distance_max[distances])
  dist_lineage_mat_rr = RR_bootstrap(metadataset = selected_gpsc_metadata,
                                     des_geo_matrix = dist_lineage_mat$des_distance_mat,
                                     ref_geo_matrix = dist_lineage_mat$ref_distance_mat,
                                     des_strain_matrix = dist_lineage_mat$des_lineage_mat,
                                     ref_strain_matrix = dist_lineage_mat$ref_lineage_mat,
                                     colyear_matrix = colyear_mat_ones,
                                     nboot = 200,
                                     subsample_min = 100)
  rr_plot = PltRR(dist_lineage_mat_rr, distance_name, distance_name)
  RRplotLists[[distances]] = print(rr_plot)

  distances_col = data.frame(distance_min=des_distance_min[distances], distance_max=des_distance_max[distances], distance_mid = (des_distance_min[distances]+des_distance_max[distances])/2)
  new_row = cbind(dist_lineage_mat_rr, distances_col)
  res_df = rbind(res_df, new_row)

  assign(distance_name, rr_plot)
  res_list=c(res_list, distance_name)
}


### bind rolling distance and within district plot
dist_district_lineage_rr_table <- rbind(res_df, rr, fill=TRUE)
for(i in 1:nrow(dist_district_lineage_rr_table)){
  dist_district_lineage_rr_table[i, "distance_window"] = paste0(dist_district_lineage_rr_table[i, "distance_min"],"-" ,dist_district_lineage_rr_table[i, "distance_max"])
}
dist_district_lineage_rr_table[nrow(dist_district_lineage_rr_table), "distance_window"] = "within \ndistrict"
dist_district_lineage_rr_table$distance_window = as.factor(dist_district_lineage_rr_table$distance_window)


png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/rolling_distance_RR.png", res=1200, width = 5000, height = 3000)
rolling_distance_RR_plot <- ggplot(dist_district_lineage_rr_table, aes(x=as.factor(distance_window), y=RR_mean))+
  geom_pointrange(data=dist_district_lineage_rr_table, mapping=aes(x=as.factor(distance_window), y=RR_mean, ymin = RR_lower_ci, ymax = RR_upper_ci), alpha=0.9) +
  geom_hline(yintercept=1, linetype="dashed", color = "red") + ##dashed line at 1
  scale_y_log10()+
  ylab("Risk Ratio")+
  xlab("Distance window (km)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_rect(xmin = 6.5, xmax = 9, ymin = -0.5, ymax = 1.5,
            fill = 'blue', alpha = 0.05) +
  geom_rect(xmin = 0, xmax = 6.5, ymin = -0.5, ymax = 1.5,
            fill = 'red', alpha = 0.05)

rolling_distance_RR_plot
dev.off()


### cumulative distance vs. relative ratio, lineage
des_distance_min = rep(0,5)
des_distance_max = seq(from=20, to=100, by=20)
des_distance_range = length(des_distance_max)

res_list = c()
res_df = data.frame()

colyear_mat_ones = matrix(1, nrow=dim(selected_gpsc_metadata)[1], ncol=dim(selected_gpsc_metadata)[1])
diag(colyear_mat_ones) <- NA

RRplotLists = list()
RRdistanceLists = list()
RRdframLists = list()

for(distances in 1:des_distance_range){
  distance_name = paste0("Distance: ",des_distance_min[distances],"~",des_distance_max[distances],"km")
  dist_lineage_mat = RR_matrix(meta_table = selected_gpsc_metadata, bactdate_tree = multitrees,
                               des_distance_min = des_distance_min[distances],
                               des_distance_max = des_distance_max[distances])
  dist_lineage_mat_rr = RR_bootstrap(metadataset = selected_gpsc_metadata,
                                     des_geo_matrix = dist_lineage_mat$des_distance_mat,
                                     ref_geo_matrix = dist_lineage_mat$ref_distance_mat,
                                     des_strain_matrix = dist_lineage_mat$des_lineage_mat,
                                     ref_strain_matrix = dist_lineage_mat$ref_lineage_mat,
                                     colyear_matrix = colyear_mat_ones,
                                     nboot = 50,
                                     subsample_min = 70)
  rr_plot = PltRR(dist_lineage_mat_rr, distance_name, distance_name)
  RRplotLists[[distances]] = print(rr_plot)

  distances_col = data.frame(distance_min=des_distance_min[distances], distance_max=des_distance_max[distances], distance_mid = (des_distance_min[distances]+des_distance_max[distances])/2)
  new_row = cbind(dist_lineage_mat_rr, distances_col)
  res_df = rbind(res_df, new_row)

  assign(distance_name, rr_plot)
  res_list=c(res_list, distance_name)
}

PltRR_rolling(res_df,"rolling distance")

###
### rolling divergence vs. relative ratio, location
des_divtime_min = seq(from=0, to=70, by=10)
des_divtime_max = seq(from=10, to=80, by=10)
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
                                    nboot = 20,
                                    subsample_min = 70)
  rr_plot = PltRR(loc_divtime_mat_rr, divtime_name, divtime_name)
  RRplotLists[[divtimes]] = print(rr_plot)

  divtimes_col = data.frame(divtime_min=des_divtime_min[divtimes], divtime_max=des_divtime_max[divtimes], divtime_mid = (des_divtime_min[divtimes]+des_divtime_max[divtimes])/2)
  new_row = cbind(loc_divtime_mat_rr, divtimes_col)
  res_df_divtime = rbind(res_df_divtime, new_row)

  assign(divtime_name, rr_plot)
  res_list=c(res_list, divtime_name)
}

rolling_distance_divtime <- PltRR_rolling_divtime(res_df_divtime,"rolling divergence")


############################################################
### rolling divergence vs. relative ratio | fixed distance (0-30 km)

des_divtime_min = seq(from=0, to=70, by=10)
des_divtime_max = seq(from=10, to=80, by=10)
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
                              des_distance_max = 30)
  loc_divtime_mat_rr = RR_bootstrap(metadataset = selected_gpsc_metadata,
                                    des_geo_matrix = loc_divtime_mat$des_distance_mat,
                                    ref_geo_matrix = loc_divtime_mat$ref_distance_mat,
                                    des_strain_matrix = loc_divtime_mat$des_divtime_mat,
                                    ref_strain_matrix = loc_divtime_mat$ref_divtime_mat,
                                    colyear_matrix = colyear_mat_ones,
                                    nboot = 50,
                                    subsample_min = 70)
  rr_plot = PltRR(loc_divtime_mat_rr, divtime_name, divtime_name)
  RRplotLists[[divtimes]] = print(rr_plot)

  divtimes_col = data.frame(divtime_min=des_divtime_min[divtimes], divtime_max=des_divtime_max[divtimes], divtime_mid = (des_divtime_min[divtimes]+des_divtime_max[divtimes])/2)
  new_row = cbind(loc_divtime_mat_rr, divtimes_col)
  res_df = rbind(res_df, new_row)

  assign(divtime_name, rr_plot)
  res_list=c(res_list, divtime_name)
}

rolling_distance_divtime_30km <- PltRR_rolling_divtime(res_df,"rolling divergence")
rolling_distance_divtime_30km

rolling_divergence_rr <- plot_grid(rolling_distance_divtime, NULL,rolling_distance_divtime_30km, rel_widths=c(1,0,1),rel_heights=c(1,0,1), labels=c("A","", "B"), scale=c(1,0,1), label_size = 12, align = "hv", nrow = 1)
png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/Israel_rolling_divergence_rr.png",res=600, width = 6000, height = 3000)
rolling_divergence_rr
dev.off()

### confidence interal for divergence time
### function to sample divergence time CI
Dtime_CI <- function(metadata, nboot, bd_res){
  boot_array <- array(NA, dim=c(nrow(metadata), nrow(metadata), nboot))
  for(bt in 1:nboot){
    ##########sample from the confidence intervals of the nodes##################################################
    nod.times<-bd_res$BactDating$CI
    ####sample node date from across CIs and calculate time to tips
    node.date<-rep(NA,nrow(nod.times))
    for(k in 1:nrow(nod.times)){
      node.date[k]<-runif(1,min=nod.times[k,][1],max=nod.times[k,][2])
    }
    nod.dist<-abs(outer(node.date,node.date,"-"))

    new.edge.length<-rep(NA,length(node.date))
    edges<-bd_res$BactDating$tree$edge
    new.edge.length<-mapply(function(i,j) nod.dist[i,j], i=edges[,1],j=edges[,2])
    ####need edge, Nnode, tip.label
    samp.tree<-bd_res$BactDating$tree
    samp.tree$edge.length<-new.edge.length
    res_tree<-samp.tree

    rr_lineage_location_mat_boot = RR_matrix(metadata, res_tree)

    boot_array[,,bt] <- rr_lineage_location_mat_boot$divtime_num_mat
  }

  boot.med = apply(boot_array, c(1,2), FUN=function(iboot) quantile(iboot, probs = c(0.5), na.rm=T))
  boot.lower.ci = apply(boot_array, c(1,2), FUN=function(iboot) quantile(iboot, probs = c(0.025), na.rm=T))
  boot.upper.ci = apply(boot_array, c(1,2), FUN=function(iboot) quantile(iboot, probs = c(0.975), na.rm=T))

  divtime_dist_med <- as.data.frame(to.upper(boot.med))
  colnames(divtime_dist_med) = "divtime_med"

  divtime_dist_lowerci <- as.data.frame(to.upper(boot.lower.ci))
  colnames(divtime_dist_lowerci) = "divtime_lowerci"

  divtime_dist_upperci <- as.data.frame(to.upper(boot.upper.ci))
  colnames(divtime_dist_upperci) = "divtime_upperci"

  distance_dist <- as.data.frame(to.upper(rr_lineage_location_mat_boot$distance_num_mat))
  colnames(distance_dist) = "distance"

  div_dist_gpsc = cbind(distance_dist, divtime_dist_med, divtime_dist_lowerci, divtime_dist_upperci)
  attach(div_dist_gpsc)
  div_dist_sorted_gpsc = div_dist_gpsc[order(divtime_med),]

  cum_mean_dist_gpsc = rep(NA, nrow(div_dist_sorted_gpsc))
  for (i in 1:nrow(div_dist_sorted_gpsc)){
    cum_mean_dist_gpsc[i] = mean(div_dist_sorted_gpsc[1:i,"distance"])
  }
  cum_mean_dist_div_gpsc = cbind(div_dist_sorted_gpsc, cum_mean_dist_gpsc)

  return(cum_mean_dist_div_gpsc)
}

### metadata
gpsc6_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "6"),]
gpsc8_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "8"),]
gpsc10_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "10"),]
gpsc47_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "47"),]
gpsc55_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "55"),]
### sampling confidence interval of node dates from bactdating tree
gpsc8_gubbin_tree_plot = plt_gubbin_tree(gpsc8_gubbin_tree)
gpsc8_time_tree_list = match_tree_meta(metadata, gpsc8_gubbin_tree, "GPSC8_reference", drop_tip = TRUE)
gpsc8_bd_res = time_resolved(gpsc8_time_tree_list$match_tree_meta, gpsc8_time_tree_list$tree, "/Users/hc14/Documents/SpneumoIsrael/analysis_results/bd_res_GPSC8", time="Year_collection")

gpsc55_gubbin_tree_plot = plt_gubbin_tree(gpsc55_gubbin_tree)
gpsc55_time_tree_list = match_tree_meta(metadata, gpsc55_gubbin_tree, "GPSC55_reference", drop_tip = TRUE)
gpsc55_bd_res = time_resolved(gpsc55_time_tree_list$match_tree_meta, gpsc55_time_tree_list$tree, "/Users/hc14/Documents/SpneumoIsrael/analysis_results/bd_res_GPSC55")

gpsc47_gubbin_tree_plot = plt_gubbin_tree(gpsc47_gubbin_tree)
gpsc47_time_tree_list = match_tree_meta(metadata, gpsc47_gubbin_tree, "GPSC47_reference", drop_tip = TRUE)
gpsc47_bd_res = time_resolved(gpsc47_time_tree_list$match_tree_meta, gpsc47_time_tree_list$tree, "/Users/hc14/Documents/SpneumoIsrael/analysis_results/bd_res_GPSC47")

gpsc6_gubbin_tree_plot = plt_gubbin_tree(gpsc6_gubbin_tree)
gpsc6_time_tree_list = match_tree_meta(metadata, gpsc6_gubbin_tree, "GPSC6_reference", drop_tip = TRUE)
gpsc6_bd_res = time_resolved(gpsc6_time_tree_list$match_tree_meta, gpsc6_time_tree_list$tree, "/Users/hc14/Documents/SpneumoIsrael/analysis_results/bd_res_GPSC6")

gpsc10_gubbin_tree_plot = plt_gubbin_tree(gpsc10_gubbin_tree)
gpsc10_time_tree_list = match_tree_meta(metadata, gpsc10_gubbin_tree, "GPSC10_reference", drop_tip = TRUE)
gpsc10_bd_res = time_resolved(gpsc10_time_tree_list$match_tree_meta, gpsc10_time_tree_list$tree, "/Users/hc14/Documents/SpneumoIsrael/analysis_results/bd_res_GPSC10")

mcmc_gpsc8=as.mcmc.resBactDating(gpsc8_bd_res$BactDating) #check whether converge
coda::effectiveSize(mcmc_gpsc8)

mcmc_gpsc55=as.mcmc.resBactDating(gpsc55_bd_res$BactDating)
coda::effectiveSize(mcmc_gpsc55)

mcmc_gpsc47=as.mcmc.resBactDating(gpsc47_bd_res$BactDating)
coda::effectiveSize(mcmc_gpsc47)

mcmc_gpsc6=as.mcmc.resBactDating(gpsc6_bd_res$BactDating)
coda::effectiveSize(mcmc_gpsc6)

mcmc_gpsc10=as.mcmc.resBactDating(gpsc10_bd_res$BactDating)
coda::effectiveSize(mcmc_gpsc10)


### sampling
nboot = 20

gpsc6_dist_divtime_ci <- Dtime_CI(metadata=gpsc6_metadata, nboot=20, bd_res=gpsc6_bd_res)
gpsc8_dist_divtime_ci <- Dtime_CI(metadata=gpsc8_metadata, nboot=20, bd_res=gpsc8_bd_res)
gpsc10_dist_divtime_ci <- Dtime_CI(metadata=gpsc10_metadata, nboot=20, bd_res=gpsc10_bd_res)
gpsc47_dist_divtime_ci <- Dtime_CI(metadata=gpsc47_metadata, nboot=20, bd_res=gpsc47_bd_res)
gpsc55_dist_divtime_ci <- Dtime_CI(metadata=gpsc55_metadata, nboot=20, bd_res=gpsc55_bd_res)

gpsc6_8_10_47_55_dist_divtime_ci <- rbind(gpsc6_dist_divtime_ci,
                                          gpsc8_dist_divtime_ci,
                                          gpsc10_dist_divtime_ci,
                                          gpsc47_dist_divtime_ci,
                                          gpsc55_dist_divtime_ci)
### plotting cumulative divergence time vs. distance
#GPSC6
cum_mean_dist_div_gpsc6 <- ggplot(gpsc6_dist_divtime_ci, aes(x=divtime_med, y=cum_mean_dist_gpsc, xmin=divtime_lowerci, xmax=divtime_upperci))+
  geom_line()+
  geom_ribbon(fill="blue", alpha=0.7)+
  #xlim(0,60)+
  xlab("Divergence Time \n (years)")+
  ylab("Cumulative Mean Distance (km)")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC6_distance_divtime_CI.png")
cum_mean_dist_div_gpsc6
dev.off()

#GPSC8
cum_mean_dist_div_gpsc8 <- ggplot(gpsc8_dist_divtime_ci, aes(x=divtime_med, y=cum_mean_dist_gpsc, xmin=divtime_lowerci, xmax=divtime_upperci))+
  geom_line()+
  geom_ribbon(fill="blue", alpha=0.7)+
  #xlim(0,60)+
  xlab("Divergence Time \n (years)")+
  ylab("Cumulative Mean Distance (km)")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC8_distance_divtime_CI.png")
cum_mean_dist_div_gpsc8
dev.off()

#GPSC10
cum_mean_dist_div_gpsc10 <- ggplot(gpsc10_dist_divtime_ci, aes(x=divtime_med, y=cum_mean_dist_gpsc, xmin=divtime_lowerci, xmax=divtime_upperci))+
  geom_line()+
  geom_ribbon(fill="blue", alpha=0.7)+
  #xlim(0,60)+
  xlab("Divergence Time \n (years)")+
  ylab("Cumulative Mean Distance (km)")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC10_distance_divtime_CI.png")
cum_mean_dist_div_gpsc10
dev.off()

#GPSC47
cum_mean_dist_div_gpsc47 <- ggplot(gpsc47_dist_divtime_ci, aes(x=divtime_med, y=cum_mean_dist_gpsc, xmin=divtime_lowerci, xmax=divtime_upperci))+
  geom_line()+
  geom_ribbon(fill="blue", alpha=0.7)+
  #xlim(0,60)+
  xlab("Divergence Time \n (years)")+
  ylab("Cumulative Mean Distance (km)")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC47_distance_divtime_CI.png")
cum_mean_dist_div_gpsc47
dev.off()

#GPSC55
cum_mean_dist_div_gpsc55 <- ggplot(gpsc55_dist_divtime_ci, aes(x=divtime_med, y=cum_mean_dist_gpsc, xmin=divtime_lowerci, xmax=divtime_upperci))+
  geom_line()+
  geom_ribbon(fill="blue", alpha=0.7)+
  #xlim(0,60)+
  xlab("Divergence Time \n (years)")+
  ylab("Cumulative Mean Distance (km)")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC55_distance_divtime_CI.png")
cum_mean_dist_div_gpsc55
dev.off()

#GPSC6+8+10+47+55
cum_mean_dist_div_gpsc_all <- ggplot(gpsc6_8_10_47_55_dist_divtime_ci, aes(x=divtime_med, y=cum_mean_dist_gpsc, xmin=divtime_lowerci, xmax=divtime_upperci))+
  geom_line()+
  geom_ribbon(fill="blue", alpha=0.7)+
  xlab("Divergence Time \n (years)")+
  xlim(0,160)+
  ylab("Cumulative Mean Distance (km)")+
  ylim(0,80)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSCall_distance_divtime_CI.png", res=1200, width = 5000, height = 5000)
cum_mean_dist_div_gpsc_all
dev.off()

save.image("/Users/hc14/Documents/SpneumoIsrael/analysis_results/IsraelRR.RData")
load("/Users/hc14/Documents/SpneumoIsrael/analysis_results/IsraelRR.RData")
####
#### Divtime Distance CI
####
'''
Take each GPSC
1. Sample bactdating tree (internal node)
2. Calculate pairwise evolutionary divergence time
3. Identify pairs within evolutionary divergence time
4. Sample pairs with replacement (bootstrap)
5. Extract spatial distance per evolutionary window
6. Calculate mean spatial  per evolutionary window
'''

Dtime_Dist_CI <- function(metadata, bd_res, bdtree_nboot=20, divtime_nboot=20, upper_window=50)
{
  gpsc_rrmatrix <- RR_matrix(meta_table = metadata, bactdate_tree = bd_res$BactDating$tree)
  nboot = bdtree_nboot
  boot_array <- array(NA, dim=c(nrow(metadata), nrow(metadata), nboot))
  for(bt in 1:nboot){
    ##########sample from the confidence intervals of the nodes##################################################
    nod.times<-bd_res$BactDating$CI #181 nodes + tips
    ####sample node date from across CIs and calculate time to tips
    node.date<-rep(NA,nrow(nod.times))

    for(k in 1:nrow(nod.times)){
      lowerci = nod.times[k,][1]
      upperci = nod.times[k,][2]
      midci = mean(lowerci, upperci)
      node.date[k]<-rtnorm(n=1,mean= midci, lower = lowerci, upper=upperci)
    }

    nod.dist <- outer(node.date,node.date,"-")

    new.edge.length<-rep(NA,length(node.date))
    edges<-bd_res$BactDating$tree$edge
    new.edge.length<-mapply(function(i,j) nod.dist[i,j], i=edges[,2],j=edges[,1])
    ####need edge, Nnode, tip.label
    samp.tree<-bd_res$BactDating$tree
    samp.tree$edge.length<-new.edge.length
    res_tree<-samp.tree

    rr_lineage_location_mat_boot = RR_matrix(metadata, res_tree)
    boot_array[,,bt] <- rr_lineage_location_mat_boot$divtime_num_mat
  }

  for (i in 1:nrow(metadata)) {
    gpsc_divtime_boot_mean=apply(boot_array, c(1,2), FUN=function(iboot) mean(iboot))
  }
  gpsc_distance_dist <- as.data.frame(to.upper(gpsc_rrmatrix$distance_num_mat))
  gpsc_divtime_dist <- as.data.frame(to.upper(gpsc_divtime_boot_mean))
  gpsc_dist_divtime_boot <- cbind(gpsc_distance_dist, gpsc_divtime_dist)
  colnames(gpsc_dist_divtime_boot) <- c("distance", "divtime")

  nboot = divtime_nboot
  x_divtime_max = seq(from=5, to=upper_window, by=1)
  x_divtime_min = rep(0, length(x_divtime_max))
  #x_divtime_min = seq(from=0, to=20, by=1)

  boot_select <- data.frame()
  for (win in 1:length(x_divtime_max)) {
    boot.mean <- c()
    for (n in 1:nboot) {
      a = which(gpsc_dist_divtime_boot$divtime > x_divtime_min[win] &
                  gpsc_dist_divtime_boot$divtime < x_divtime_max[win])        #select a certain region (index)
      if (length(a)!=0){
        b = sample(x=a, size=100, replace = T)  #subsample with replacement
        resample.mean <- mean(gpsc_dist_divtime_boot[b, ]$distance)
        boot.mean <- c(boot.mean, resample.mean)
      }
    }
    distance.med <- quantile(boot.mean, probs = c(0.5), na.rm=T) #calculate confidence interval
    distance.mean <- mean(boot.mean)
    distance.upperci <- quantile(boot.mean, probs = c(0.975), na.rm=T)
    distance.lowerci <- quantile(boot.mean, probs = c(0.025), na.rm=T)
    mid_win <- mean(c(x_divtime_min[win], x_divtime_max[win]))
    window_divtime_distance <- cbind(x_divtime_min[win], x_divtime_max[win], mid_win, distance.med, distance.mean, distance.upperci, distance.lowerci)
    boot_select <- rbind(boot_select, window_divtime_distance)
  }

  res_bootstrap_divtime_dist = list("window_divtime_dist" = boot_select, "divtime_dist_table" = gpsc_dist_divtime_boot)
  return(res_bootstrap_divtime_dist)
}

gpsc6_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "6"),] # 91 samples
gpsc6_rrmatrix <- RR_matrix(meta_table = gpsc6_metadata, bactdate_tree = gpsc6_bd_res$BactDating$tree)

gpsc8_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "8"),]
gpsc8_rrmatrix <- RR_matrix(meta_table = gpsc8_metadata, bactdate_tree = gpsc8_bd_res$BactDating$tree)

gpsc10_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "10"),]
gpsc10_rrmatrix <- RR_matrix(meta_table = gpsc10_metadata, bactdate_tree = gpsc10_bd_res$BactDating$tree)

gpsc47_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "47"),]
gpsc47_rrmatrix <- RR_matrix(meta_table = gpsc47_metadata, bactdate_tree = gpsc47_bd_res$BactDating$tree)

gpsc55_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$GPSC_PoPUNK2 == "55"),]
gpsc55_rrmatrix <- RR_matrix(meta_table = gpsc55_metadata, bactdate_tree = gpsc55_bd_res$BactDating$tree)

gpsc6_boot_divtime_dist_ci <- Dtime_Dist_CI(metadata = gpsc6_metadata, bd_res = gpsc6_bd_res, bdtree_nboot=100, divtime_nboot=100, upper_window = 50)
gpsc6_dist_div <- ggplot(gpsc6_boot_divtime_dist_ci$window_divtime_dist, aes(x=mid_win, y=distance.mean, ymin=distance.lowerci, ymax=distance.upperci))+
  geom_line()+
  geom_ribbon(fill="pink", alpha=0.5)+
  #xlim(0,60)+
  xlab("Divergence Time \n (years)")+
  ylab("Distance (km)")+
  ylim(c(0,165))+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC6_distance_divtime_cumCI.png")
gpsc6_dist_div
dev.off()

gpsc8_boot_divtime_dist_ci <- Dtime_Dist_CI(metadata = gpsc8_metadata, bd_res = gpsc8_bd_res, upper_window = 80)
gpsc8_dist_div <- ggplot(gpsc8_boot_divtime_dist_ci$window_divtime_dist, aes(x=mid_win, y=distance.mean, ymin=distance.lowerci, ymax=distance.upperci))+
  geom_line()+
  geom_ribbon(fill="pink", alpha=0.5)+
  #xlim(0,60)+
  xlab("Divergence Time \n (years)")+
  ylab("Distance (km)")+
  ylim(c(0,165))+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC8_distance_divtime_cumCI.png")
gpsc8_dist_div
dev.off()

gpsc10_boot_divtime_dist_ci <- Dtime_Dist_CI(metadata = gpsc10_metadata, bd_res = gpsc10_bd_res, upper_window = 365)
gpsc10_dist_div <- ggplot(gpsc10_boot_divtime_dist_ci$window_divtime_dist, aes(x=mid_win, y=distance.mean, ymin=distance.lowerci, ymax=distance.upperci))+
  geom_line()+
  geom_ribbon(fill="pink", alpha=0.5)+
  #xlim(0,60)+
  xlab("Divergence Time \n (years)")+
  ylab("Distance (km)")+
  ylim(c(0,165))+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC10_distance_divtime_cumCI.png")
gpsc10_dist_div
dev.off()

gpsc47_boot_divtime_dist_ci <- Dtime_Dist_CI(metadata = gpsc47_metadata, bd_res = gpsc47_bd_res, upper_window = 200)
gpsc47_dist_div <- ggplot(gpsc47_boot_divtime_dist_ci$window_divtime_dist, aes(x=mid_win, y=distance.mean, ymin=distance.lowerci, ymax=distance.upperci))+
  geom_line()+
  geom_ribbon(fill="pink", alpha=0.5)+
  #xlim(0,60)+
  xlab("Divergence Time \n (years)")+
  ylab("Distance (km)")+
  ylim(c(0,165))+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC47_distance_divtime_cumCI.png")
gpsc47_dist_div
dev.off()

gpsc55_boot_divtime_dist_ci <- Dtime_Dist_CI(metadata = gpsc55_metadata, bd_res = gpsc55_bd_res, upper_window = 25)
gpsc55_dist_div <- ggplot(gpsc55_boot_divtime_dist_ci$window_divtime_dist, aes(x=mid_win, y=distance.mean, ymin=distance.lowerci, ymax=distance.upperci))+
  geom_line()+
  geom_ribbon(fill="pink", alpha=0.5)+
  #xlim(0,60)+
  xlab("Divergence Time \n (years)")+
  ylab("Distance (km)")+
  ylim(c(0,165))+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC55_distance_divtime_cumCI.png")
gpsc55_dist_div
dev.off()

GPSCall_treeboot_divtime_dist <- rbind(gpsc6_boot_divtime_dist_ci$divtime_dist_table,
                                       gpsc8_boot_divtime_dist_ci$divtime_dist_table,
                                       gpsc10_boot_divtime_dist_ci$divtime_dist_table,
                                       gpsc47_boot_divtime_dist_ci$divtime_dist_table,
                                       gpsc55_boot_divtime_dist_ci$divtime_dist_table)

nboot = 100
x_divtime_max = seq(from=5, to=365, by=1)
x_divtime_min = rep(0, length(x_divtime_max))
#x_divtime_min = seq(from=0, to=20, by=1)

boot_select <- data.frame()
for (win in 1:length(x_divtime_max)) {
  boot.mean <- c()
  for (n in 1:nboot) {
    a = which(GPSCall_treeboot_divtime_dist$divtime > x_divtime_min[win] &
                GPSCall_treeboot_divtime_dist$divtime < x_divtime_max[win])        #select a certain region (index)
    if (length(a)!=0){
      b = sample(x=a, size=100, replace = T)  #subsample without replacement (index)
      resample.mean <- mean(GPSCall_treeboot_divtime_dist[b, ]$distance)
      boot.mean <- c(boot.mean, resample.mean)
    }
  }
  distance.med <- quantile(boot.mean, probs = c(0.5), na.rm=T) #calculate confidence interval
  distance.mean <- mean(boot.mean)
  distance.upperci <- quantile(boot.mean, probs = c(0.975), na.rm=T)
  distance.lowerci <- quantile(boot.mean, probs = c(0.025), na.rm=T)
  mid_win <- mean(c(x_divtime_min[win], x_divtime_max[win]))
  window_divtime_distance <- cbind(x_divtime_min[win], x_divtime_max[win], mid_win, distance.med, distance.mean, distance.upperci, distance.lowerci)
  boot_select <- rbind(boot_select, window_divtime_distance)
}

GPSCall_dist_div <- ggplot(boot_select, aes(x=mid_win, y=distance.mean, ymin=distance.lowerci, ymax=distance.upperci))+
  geom_line()+
  geom_ribbon(fill="pink", alpha=0.5)+
  #xlim(0,20)+
  xlab("Divergence Time \n (years)")+
  ylab("Distance (km)")+
  ylim(c(0,165))+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSCall_distance_divtime_rollingCI.png")
GPSCall_dist_div
dev.off()

Dtime_Dist_CI_posterior_tree <- function(metadata, bd_res, Ntreestosample=100, divtime_nboot=20, upper_window=50)
{
  ################################################################
  ## Sample trees from BactDating posterior
  ## Noemie Lefrancq
  ## Last update 21/09/2023
  ################################################################

  ## Load Bactdating result

  #resbd <- readRDS('bd_GPSC8')
  Ntips = length(bd_res$BactDating$tree$tip.label)
  Ntreestosample = 100
  ## ID of tree to sample
  tree_ID = round(seq(1, nrow(bd_res$BactDating$record), length.out = Ntreestosample))
  ################################################################
  ## Reconstruct one tree from the trace
  ################################################################
  ## Withdraw node heights

  posterior_sampling <- array(NA, dim=c(Ntips, Ntips, Ntreestosample))
  for(samp in 1:Ntreestosample){
    tip_dates = bd_res$BactDating$record[tree_ID[samp],1:Ntips]
    node_dates = bd_res$BactDating$record[tree_ID[samp],(Ntips+1):(2*Ntips-1)]
    tipandnodes_dates = c(tip_dates, node_dates)

    ## Take BD final tree for template
    reconstructed_tree = bd_res$BactDating$tree
    ## Replace edge lengths, based on the tip and node dates that were extracted above
    reconstructed_tree$edge.length = tipandnodes_dates[reconstructed_tree$edge[,2]] - tipandnodes_dates[reconstructed_tree$edge[,1]]
    ## pairwise distance
    reconstructed_pairwise_dist = cophenetic.phylo(reconstructed_tree)
    posterior_sampling[,,samp] <- reconstructed_pairwise_dist
  }

  for (i in 1:Ntips) {
    gpsc_divtime_posterior_mean <- apply(posterior_sampling, c(1,2), FUN=function(isamp) mean(isamp))
    gpsc_divtime_posterior_upperci <- apply(posterior_sampling, c(1,2), FUN = function(isamp) quantile(isamp, probs = 0.975, na.rm = T))
    gpsc_divtime_posterior_lowerci <- apply(posterior_sampling, c(1,2), FUN = function(isamp) quantile(isamp, probs = 0.075, na.rm = T))
  }

  gpsc_rrmatrix <- RR_matrix(meta_table = metadata, bactdate_tree = bd_res$BactDating$tree)
  gpsc_distance_dist <- as.data.frame(to.upper(gpsc_rrmatrix$distance_num_mat))
  gpsc_divtime_dist <- as.data.frame(to.upper(gpsc_divtime_posterior_mean))
  gpsc_dist_divtime_boot <- cbind(gpsc_distance_dist, gpsc_divtime_dist)
  colnames(gpsc_dist_divtime_boot) <- c("distance", "divtime")

  nboot = divtime_nboot
  x_divtime_max = seq(from=5, to=upper_window, by=1)
  x_divtime_min = rep(0, length(x_divtime_max))
  #x_divtime_min = seq(from=0, to=20, by=1)

  boot_select <- data.frame()
  for (win in 1:length(x_divtime_max)) {
    boot.mean <- c()
    for (n in 1:nboot) {
      a = which(gpsc_dist_divtime_boot$divtime > x_divtime_min[win] &
                  gpsc_dist_divtime_boot$divtime < x_divtime_max[win])        #select a certain region (index)
      if (length(a)!=0){
        b = sample(x=a, size=100, replace = T)  #subsample with replacement
        resample.mean <- mean(gpsc_dist_divtime_boot[b, ]$distance)
        boot.mean <- c(boot.mean, resample.mean)
      }
    }
    distance.med <- quantile(boot.mean, probs = c(0.5), na.rm=T) #calculate confidence interval
    distance.mean <- mean(boot.mean)
    distance.upperci <- quantile(boot.mean, probs = c(0.975), na.rm=T)
    distance.lowerci <- quantile(boot.mean, probs = c(0.025), na.rm=T)
    mid_win <- mean(c(x_divtime_min[win], x_divtime_max[win]))
    window_divtime_distance <- cbind(x_divtime_min[win], x_divtime_max[win], mid_win, distance.med, distance.mean, distance.upperci, distance.lowerci)
    boot_select <- rbind(boot_select, window_divtime_distance)
  }

  res_bootstrap_divtime_dist = list("window_divtime_dist" = boot_select, "divtime_dist_table" = gpsc_dist_divtime_boot)
  return(res_bootstrap_divtime_dist)
}

gpsc6_Dtime_Dist_CI_posterior_tree <- Dtime_Dist_CI_posterior_tree(gpsc6_metadata, gpsc6_bd_res)
gpsc8_Dtime_Dist_CI_posterior_tree
gpsc10_Dtime_Dist_CI_posterior_tree
gpsc47_Dtime_Dist_CI_posterior_tree
gpsc55_Dtime_Dist_CI_posterior_tree


