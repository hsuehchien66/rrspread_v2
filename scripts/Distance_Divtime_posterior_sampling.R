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
