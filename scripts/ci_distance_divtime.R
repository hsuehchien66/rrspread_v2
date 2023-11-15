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

  return(cum_mean_dist_div_gpsc)
}


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

### boostrapping distance within divtime windows
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
        break
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
