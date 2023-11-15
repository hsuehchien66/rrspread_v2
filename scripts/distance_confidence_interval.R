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


