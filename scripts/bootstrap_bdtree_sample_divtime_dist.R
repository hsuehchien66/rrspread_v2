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
        b = sample(x=a, size=100, replace = T)  #subsample without replacement (index)
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
