# euclidean_distance between two vectors
euclidean_distance <- function(vector1, vector2) { sqrt( sum( (vector1 - vector2) ^ 2)) }

# returns a vector of the top three signatures' rows identified in a vector of
# mut sigs intensities; will return tied intensities in low --> high row order
# and others low --> high intensity
top_three_sigs <- function(mut_matrix_weights) {
  intesities <- order(mut_matrix_weights)
  top_three <- intesities[(length(intesities) - 2) : length(intesities)]
  return(top_three)
}

# Takes one vector of mut sig intensities and calculates 
# 1) cosine similary MEAN and VARIANCE
# 2) euclidean distance MEAN and VARIANCE
# variance is NA for 1 sample comparisons
calculate_all_stats <- function(mut_sub_matrix, mut_full_df) {
  distance_vals <- sapply(1:ncol(mut_sub_matrix), function(colnum){
    cosine_val <- MutationalPatterns::cos_sim(mut_sub_matrix[,colnum], mut_full_df[,1])
    euclid_val <- euclidean_distance(mut_sub_matrix[,colnum], mut_full_df[,1])

    both_stats <- c(cosine_val, euclid_val)
    return(both_stats)
  })

  mean_cos_val <- mean(distance_vals[1,])
  mean_eu_val <- mean(distance_vals[2,])
  
  variance_cos_val <- var(distance_vals[1,])
  variance_eu_val <- var(distance_vals[2,])

  all_stats <- t(matrix( c(mean_cos_val, mean_eu_val, variance_cos_val, variance_eu_val),
                      dimnames = list(c("mean cos", "mean euclid",
                                        "cos variance", "euclid variance"), 
                                      c("values"))))
  return(all_stats)
}

# generates a data frame that contains the cosine mean/var and euclidean mean/var
# values at each of five, ten, twenty-five, fifty, and hundred variant subsamples
create_plot_input <- function(software, mut_full_df) {
  graphable_all_muts <- lapply(list(five_subs, ten_subs, twenty_five_subs, fifty_subs, hundred_subs),
         function(subsample_df){
           software_output <- subsample_to_mut_sigs(software, subsample_df)
           stats_matrix <- calculate_all_stats(software_output, mut_full_df)
           return(stats_matrix)
         })

  graphable_all_muts <- do.call(rbind.data.frame, graphable_all_muts)
  rownames(graphable_all_muts) <- (paste(c(5, 10, 25, 50, 100),"vars"))
  colnames(graphable_all_muts) <- c("mean cos", "mean euclid", "cos variance", "euclid variance")
  return(graphable_all_muts)
}

# generate graph inputs for all three software
qp_plot_input <- create_plot_input("qp", qp_mut_df)
ds_plot_input <- create_plot_input("ds", ds_mut_df)
sl_plot_input <- create_plot_input("sl", sl_mut_df)
