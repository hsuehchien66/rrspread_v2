### using nod.quants
metadata = gpsc6_metadata
nboot = 20
bactdate_tree = gpsc6_bd_res$BactDating$tree
all.date <- allDates(gpsc6_bd_res$BactDating$tree) ## i internal nodes + (i+1) leaves = 2*i +1 = 161
node.date <- nodeDates(gpsc6_bd_res$BactDating$tree) ## i internal nodes, i=80

matched_tree_meta_table = match_tree_meta(Metadata = metadata, Tree = bactdate_tree,
                                          drop_tip_name = "", drop_tip = FALSE)
matched_table = matched_tree_meta_table$match_tree_meta
lane.names <- matched_table$lane_id
vector_coltime = matched_table$Year_collection
vector_coltime = vector_coltime + runif(length(vector_coltime), min = -0.01, max = 0.01)

coltime_mat = abs(outer(vector_coltime, vector_coltime, "-"))
diag(coltime_mat) <- NA
coltime_num_mat = coltime_mat

nod.times <- gpsc6_bd_res$BactDating$CI
all.date <- allDates(gpsc6_bd_res$BactDating$tree)
nod.times.ci <- cbind(nod.times, all.date)
colnames(nod.times.ci) <- c("lowerci", "upperci", "mean")

nboot=10
nod.quants<-seq(0,1,0.0275)
res_dist_divtime = data.frame(matrix(NA, nrow = length(x_divtime_max), ncol = 5))
colnames(res_dist_divtime) = c("min_divtime", "max_divtime", "med_dist", "lowerci_dist", "upperci_dist")
x_divtime_max = seq(from=10, to=150, by=10)
x_divtime_min = rep(0, length(x_divtime_max))
window_data <- vector("list", length(x_divtime_max))

for (win in 1:length(x_divtime_max)) {
  window_accepted_divtime = data.frame()
  for (j in 1:nboot) {
      node.date<-rep(NA,nrow(nod.times.ci))
      for(k in 1:nrow(nod.times)){
        node.date[k]<-quantile(nod.times[k,],probs=nod.quants[j])
        # node.date[k]<-runif(1,min=nod.times[k,][1],max=nod.times[k,][2])
      }

      nod.dist<-abs(outer(node.date,node.date,"-")) ## why not use (coltime1 + coltime2 - 2*internal nodes)/2
      new.edge.length<-rep(NA,length(node.date))
      edges<-gpsc6_bd_res$BactDating$tree$edge ## dimension of nod.dist and edges don't match
      new.edge.length<-mapply(function(i,j) nod.dist[i,j], i=edges[,1],j=edges[,2])

      samp.tree<-gpsc6_bd_res$BactDating$tree
      samp.tree$edge.length<-new.edge.length
      res_tree<-samp.tree

      rr_lineage_location_mat_boot = RR_matrix(gpsc6_metadata, res_tree)
      pw_divtime <- rr_lineage_location_mat_boot$divtime_num_mat
      boot_divtime <- rr_lineage_location_mat_boot$divtime_num_mat

      divtime_dist <- as.data.frame(to.upper(boot_divtime))
      colnames(divtime_dist) = "divtime"

      distance_dist <- as.data.frame(to.upper(rr_lineage_location_mat_boot$distance_num_mat))
      colnames(distance_dist) = "distance"

      div_dist_gpsc = cbind(distance_dist, divtime_dist)
      attach(div_dist_gpsc)
      div_dist_sorted_gpsc = div_dist_gpsc[order(divtime),]

      for (pairi in 1:nrow(div_dist_sorted_gpsc)) {
        if (div_dist_sorted_gpsc[pairi,"divtime"] > x_divtime_max[win]){
          break ## divtime out of window
        }else{
          tmp_accepted_divtime <- cbind(div_dist_sorted_gpsc[pairi,"divtime"], x_divtime_min[win], x_divtime_max[win], div_dist_sorted_gpsc[pairi,"distance"])
          # mean geographic distance
          window_accepted_divtime <- rbind(window_accepted_divtime, tmp_accepted_divtime)
        }
      }
  }
  colnames(window_accepted_divtime) <- c("divtime", "divtime_window_min", "divtime_window_max", "distance")
  window_data[[win]] <- window_accepted_divtime
  boot.dist.med = quantile(window_accepted_divtime[,4], probs = c(0.5), na.rm=T)
  boot.dist.lower.ci = quantile(window_accepted_divtime[,4], probs = c(0.025), na.rm=T)
  boot.dist.upper.ci = quantile(window_accepted_divtime[,4], probs = c(0.975), na.rm=T)
  window_dist_boot = cbind(x_divtime_min[win], x_divtime_max[win], boot.dist.med, boot.dist.lower.ci, boot.dist.upper.ci)
  window_dist_boot = as.data.frame(window_dist_boot)
  colnames(window_dist_boot) = c("min_divtime", "max_divtime", "med_dist", "lowerci_dist", "upperci_dist")
  res_dist_divtime[win,] <- window_dist_boot
}







