### function to sample divergence time CI
Dtime_CI <- function(metadata, nboot, bd_res){
  boot_array <- array(NA, dim=c(nrow(metadata), nrow(metadata), nboot))
  for(bt in 1:nboot){
    ##########sample from the confidence intervals of the nodes##################################################
    nod.times<-bd_res$BactDating$CI
    ####sample node date from across CIs and calculate time to tips
    node.date<-rep(NA,nrow(nod.times))
    for(k in 1:nrow(nod.times)){
      node.date[k]<-runif(1,min=nod.times[k,][1],max=nod.times[k,][2])  ##rnorm
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

  res_list = list("sampled_tree" = res_tree,
                  "cum_mean_dist_div_gpsc" = cum_mean_dist_div_gpsc)
  return(res_list)
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

gpsc6_dist_divtime_ci <- Dtime_CI(metadata=gpsc6_metadata, nboot=1000, bd_res=gpsc6_bd_res)
gpsc8_dist_divtime_ci <- Dtime_CI(metadata=gpsc8_metadata, nboot=1000, bd_res=gpsc8_bd_res)
gpsc10_dist_divtime_ci <- Dtime_CI(metadata=gpsc10_metadata, nboot=1000, bd_res=gpsc10_bd_res)
gpsc47_dist_divtime_ci <- Dtime_CI(metadata=gpsc47_metadata, nboot=1000, bd_res=gpsc47_bd_res)
gpsc55_dist_divtime_ci <- Dtime_CI(metadata=gpsc55_metadata, nboot=1000, bd_res=gpsc55_bd_res)

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

ggplot(gpsc6_8_10_47_55_dist_divtime_ci, aes(x=divtime_med, y=cum_mean_dist_gpsc))+
  geom_line()+
  geom_ribbon(fill="blue", alpha=0.7,aes(x=divtime_med, ymin=divtime_lowerci, ymax=divtime_upperci))+
  xlab("Divergence Time \n (years)")+
  xlim(0,160)+
  ylab("Cumulative Mean Distance (km)")+
  ylim(0,80)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(family = "Helvetica"))


png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSCall_distance_divtime_CI.png", res=1200, width = 5000, height = 5000)
cum_mean_dist_div_gpsc_all
dev.off()

write.csv(x = gpsc6_8_10_47_55_dist_divtime_ci,file = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSCall_distance_divtime_CI.csv", quote = FALSE, row.names = FALSE, col.names = TRUE)


png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC55_bdtree.png", res=1200, width = 5000, height = 5000)
plot(gpsc55_bd_res$BactDating,show.tip.label = F, "treeCI", title)
axisPhylo(backward = F, side = 1)
dev.off()

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC6_bdtree.png", res=1200, width = 5000, height = 5000)
plot(gpsc6_bd_res$BactDating,show.tip.label = F, "treeCI")
axisPhylo(backward = F, side = 1)
dev.off()

p1 <- ggdraw() +
  draw_image("/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC6_bdtree.png", scale = 1)+
  theme(plot.margin = unit(c(-2,-2,-2,-2), "cm"))
p2 <- ggdraw() +
  draw_image("/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC6_distance_divtime_cumCI.png", scale = 1)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

time_resolved_phylo_distance <- plot_grid(p1, NULL,p2, rel_widths=c(1,0,1),rel_heights=c(1,0,1), labels=c("A","", "B"), scale=c(1.2,0,1), label_size = 12, align = "hv", nrow = 1)

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/Israel_time_resolved_phylo_distance_boots.png",res=600, width = 6000, height = 3000)
time_resolved_phylo_distance
dev.off()

gpsc6_8_10_47_55_dist_divtime_ci_sorted = gpsc6_8_10_47_55_dist_divtime_ci[order(gpsc6_8_10_47_55_dist_divtime_ci$divtime_lowerci), ]
boot_dist_divtime_df = data.frame()
x_divtime_min = seq(from=0, to=140, by=10)
x_divtime_max = seq(from=10, to=150, by=10)
nboot = 100
res_dist_divtime = data.frame(matrix(NA, nrow = length(x_divtime_max), ncol = 5))
colnames(res_dist_divtime) = c("min_divtime", "max_divtime", "med_dist", "lowerci_dist", "upperci_dist")
for(win in 1:length(x_divtime_max)){
  window_accepted_divtime = data.frame()
  for(bt in 1:nboot){
    while (pairi in 1:nrow(gpsc6_8_10_47_55_dist_divtime_ci)) {
      if (gpsc6_8_10_47_55_dist_divtime_ci[pairi,"divtime_lowerci"] > x_divtime_max[win]){

      }else{
        boot_divtime <- runif(1,min=gpsc6_8_10_47_55_dist_divtime_ci[pairi,"divtime_lowerci"],
                              max=gpsc6_8_10_47_55_dist_divtime_ci[pairi,"divtime_upperci"])
        if(boot_divtime > x_divtime_min[win] & boot_divtime < x_divtime_max[win]){
          tmp_accepted_divtime <- cbind(boot_divtime, x_divtime_min[win], x_divtime_max[win], gpsc6_8_10_47_55_dist_divtime_ci[pairi,"distance"])
          window_accepted_divtime <- rbind(window_accepted_divtime, tmp_accepted_divtime)
        }
      }
    }
  }
  boot.dist.med = quantile(window_accepted_divtime[,4], probs = c(0.5), na.rm=T)
  boot.dist.lower.ci = quantile(window_accepted_divtime[,4], probs = c(0.025), na.rm=T)
  boot.dist.upper.ci = quantile(window_accepted_divtime[,4], probs = c(0.975), na.rm=T)
  window_dist_boot = cbind(x_divtime_min[win], x_divtime_max[win], boot.dist.med, boot.dist.lower.ci, boot.dist.upper.ci)
  window_dist_boot = as.data.frame(window_dist_boot)
  colnames(window_dist_boot) = c("min_divtime", "max_divtime", "med_dist", "lowerci_dist", "upperci_dist")
  res_dist_divtime[win,] <- window_dist_boot
  print(window_dist_boot)
  print(res_dist_divtime)
}


png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC8_bdtree.png", res=1200, width = 5000, height = 5000)
plot(gpsc8_bd_res$BactDating,show.tip.label = F, "treeCI")
axisPhylo(backward = F, side = 1)
dev.off()

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC10_bdtree.png", res=1200, width = 5000, height = 5000)
plot(gpsc10_bd_res$BactDating,show.tip.label = F, "treeCI")
axisPhylo(backward = F, side = 1)
dev.off()

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC47_bdtree.png", res=1200, width = 5000, height = 5000)
plot(gpsc47_bd_res$BactDating,show.tip.label = F, "treeCI")
axisPhylo(backward = F, side = 1)
dev.off()

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC55_bdtree.png", res=1200, width = 5000, height = 5000)
plot(gpsc55_bd_res$BactDating,show.tip.label = F, "treeCI")
axisPhylo(backward = F, side = 1)
dev.off()

GPSC8_bdtree <- ggdraw() +
  draw_image("/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC8_bdtree.png", scale = 1)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
GPSC10_bdtree <- ggdraw() +
  draw_image("/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC10_bdtree.png", scale = 1)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
GPSC47_bdtree <- ggdraw() +
  draw_image("/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC47_bdtree.png", scale = 1)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
GPSC55_bdtree <- ggdraw() +
  draw_image("/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/GPSC55_bdtree.png", scale = 1)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))


gpscs_bd_tree <- plot_grid(GPSC8_bdtree, GPSC10_bdtree, GPSC47_bdtree, GPSC55_bdtree,  labels=c("A", "B", "C", "D"), scale=c(1,1,1,1), label_size = 12, align = "hv", nrow = 2)

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/gpscs_bd_tree.png",res=600, width = 6000, height = 6000)
gpscs_bd_tree
dev.off()


### metadata
pwDistancePlot <- function(metadata, district_name){
  dist <- gen_dist_mat(metadata = metadata, min_dist = 0, max_dist = 40)

  pw_dist <- as.data.frame(to.upper(dist$pairwise_distance_matrix))
  colnames(pw_dist) = "distance"

  pw_dist_p <- pw_dist %>%
    ggplot( aes(x=distance)) +
    geom_histogram( binwidth=5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    ggtitle(district_name) +
    xlab("Distance (km)")+
    ylab("Pair Count")+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())

  return(pw_dist_p)
}

centre_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$REGION == "Center District"),]
haifa_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$REGION == "Haifa District"),]
jerusalem_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$REGION == "Jerusalem District"),]
north_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$REGION == "North District"),]
south_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$REGION == "South District"),]
telaviv_metadata <- selected_gpsc_metadata[which(selected_gpsc_metadata$REGION == "Tel Aviv District"),]

centre_pwdist_p <- pwDistancePlot(centre_metadata, "Central District")
haifa_pwdist_p <- pwDistancePlot(haifa_metadata, "Haifa District")
jerusalem_pwdist_p <- pwDistancePlot(jerusalem_metadata, "Jerusalem District")
north_pwdist_p <- pwDistancePlot(north_metadata, "North District")
south_pwdist_p <- pwDistancePlot(south_metadata, "South District")
telaviv_pwdist_p <- pwDistancePlot(telaviv_metadata, "Tel Aviv District")

district_pw_dist <- plot_grid(centre_distance_dist_p, haifa_pwdist_p,
                              jerusalem_pwdist_p, north_pwdist_p,
                              south_pwdist_p, telaviv_pwdist_p,
                              labels = "AUTO", label_size = 12, align = "hv", nrow = 3)

png(filename = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/figures/Israel_district_pwdistance.png",res=600, width = 4000, height = 6000)
district_pw_dist
dev.off()

### roottotip check
GPSC6_rtip <- ggdraw() +
  draw_image("/Users/hc14/Documents/SpneumoIsrael/analysis_results/06_bactdating_ouptput/GPSC6_Israel_rtd_tree.pdf", scale = 1)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))





x_divtime_min = seq(from=0, to=130, by=10)
x_divtime_max = seq(from=30, to=160, by=10)
res_dist_divtime = data.frame(matrix(NA, nrow = length(x_divtime_max), ncol = 5))
colnames(res_dist_divtime) = c("min_divtime", "max_divtime", "med_dist", "lowerci_dist", "upperci_dist")
for (i in 1:length(x_divtime_max)){
  temp.window.divtime <- gpsc6_8_10_47_55_dist_divtime_ci[which(gpsc6_8_10_47_55_dist_divtime_ci$divtime_med > x_divtime_min[i] & gpsc6_8_10_47_55_dist_divtime_ci$divtime_med < x_divtime_max[i]),  ]
  ci.dist.med <- quantile(temp.window.divtime$distance, probs = c(0.5))
  ci.dist.upperci <- quantile(temp.window.divtime$distance, probs = c(0.9))
  ci.dist.lowerci <- quantile(temp.window.divtime$distance, probs = c(0.1))
  res_dist_divtime[i, ] <- as.data.frame(cbind(x_divtime_min[i], x_divtime_max[i], ci.dist.med, ci.dist.lowerci, ci.dist.upperci))
}

for (i in 1:nrow(res_dist_divtime)){
  res_dist_divtime[i, "mid_divtime"] <- mean(c(res_dist_divtime[i, "min_divtime"], res_dist_divtime[i, "max_divtime"]))
}



ggplot(res_dist_divtime, aes(x=mid_divtime, y=med_dist, ymin=lowerci_dist, ymax=upperci_dist))+
  geom_line()+
  geom_ribbon(alpha=0.7)+
  xlab("Divergence Time \n (years)")+
  xlim(0,160)+
  ylab("Distance (km)")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

cum_mean_dist_div_gpsc_all <- ggplot(gpsc6_8_10_47_55_dist_divtime_ci, aes(x=divtime_med, y=cum_mean_dist_gpsc, xmin=divtime_lowerci, xmax=divtime_upperci))+
  geom_line()+
  geom_ribbon(fill="blue", alpha=0.7)+
  xlab("Divergence Time \n (years)")+
  xlim(0,160)+
  ylab("Cumulative Mean Distance (km)")+
  ylim(0,80)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

save(gpsc6_metadata, gpsc8_metadata, gpsc10_metadata, gpsc47_metadata, gpsc55_metadata,
     gpsc6_bd_res, gpsc8_bd_res, gpsc10_bd_res, gpsc47_bd_res, gpsc55_bd_res,
     file = "/Users/hc14/Documents/SpneumoIsrael/analysis_results/metadata_bdtree.Rdata")

nboot = 10
x_divtime_min = seq(from=0, to=100, by=10)
x_divtime_max = seq(from=30, to=130, by=10)

res_dist_divtime = data.frame(matrix(NA, nrow = length(x_divtime_max), ncol = 5))
colnames(res_dist_divtime) = c("min_divtime", "max_divtime", "med_dist", "lowerci_dist", "upperci_dist")
nboot = 10
for(win in 1:length(x_divtime_max)){
  window_accepted_divtime = data.frame()
  for(bt in 1:nboot){
    for(pairi in 1:nrow(gpsc6_8_10_47_55_dist_divtime_ci)){
      if (gpsc6_8_10_47_55_dist_divtime_ci[pairi,"divtime_lowerci"] > x_divtime_min[win] | gpsc6_8_10_47_55_dist_divtime_ci[pairi,"divtime_upperci"] < x_divtime_max[win]){
        #print(gpsc6_8_10_47_55_dist_divtime_ci[pairi,"divtime_lowerci"])
        #print(x_divtime_max[win])
        boot_divtime <- runif(1,min=gpsc6_8_10_47_55_dist_divtime_ci[pairi,"divtime_lowerci"],
                              max=gpsc6_8_10_47_55_dist_divtime_ci[pairi,"divtime_upperci"])
        if (boot_divtime > x_divtime_min[win] & boot_divtime < x_divtime_max[win]){
          tmp_accepted_divtime <- cbind(boot_divtime, x_divtime_min[win], x_divtime_max[win], gpsc6_8_10_47_55_dist_divtime_ci[pairi,"distance"])
          window_accepted_divtime <- rbind(window_accepted_divtime, tmp_accepted_divtime)
        }
      }
    }
  }
  boot.dist.med = quantile(window_accepted_divtime[,4], probs = c(0.5), na.rm=T)
  boot.dist.lower.ci = quantile(window_accepted_divtime[,4], probs = c(0.025), na.rm=T)
  boot.dist.upper.ci = quantile(window_accepted_divtime[,4], probs = c(0.975), na.rm=T)
  window_dist_boot = cbind(x_divtime_min[win], x_divtime_max[win], boot.dist.med, boot.dist.lower.ci, boot.dist.upper.ci)
  window_dist_boot = as.data.frame(window_dist_boot)
  colnames(window_dist_boot) = c("min_divtime", "max_divtime", "med_dist", "lowerci_dist", "upperci_dist")
  res_dist_divtime[win,] <- window_dist_boot
  print(bt)
}

