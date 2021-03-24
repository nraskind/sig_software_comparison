library(ggplot2)
library(gridExtra)

ggbar_plot_sigs <- function(mut_sigs_df){
  gg_sig_plot <- ggplot(mut_sigs_df, aes(x=c(1:nrow(mut_sigs_df)),
                                      y=mut_sigs_df[1:nrow(mut_sigs_df), 1])) +
    geom_bar(stat="identity") + ggtitle(colnames(mut_sigs_df)) + xlab("Signature") + ylab("Frequency") +
    scale_x_continuous(n.breaks = nrow(mut_sigs_df)) +
    scale_y_continuous(n.breaks = 10) 
  
  return(gg_sig_plot)
}

plot_stats <- function(plot_input, method="cosine", software="") {
  if (method == "cosine") {
    mean_vals <- as.numeric(plot_input[,1])
    variance_vals <- as.numeric(plot_input[,3])
    measuring <- paste(method, "similarity")
  }
  if (method == "euclidean") {
    mean_vals <- as.numeric(plot_input[,2])
    variance_vals <- as.numeric(plot_input[,4])
    measuring <- paste(method, "distance")
  }
  
  ggplot(plot_input, aes( x=c(5, 10, 25, 50, 100) )) +
    geom_line(aes(y=mean_vals), colour="black") +
    ggtitle(paste(software, measuring)) + xlab("number of variants") + ylab(measuring) +
    geom_errorbar( aes(ymin=(mean_vals - variance_vals),
                       ymax=(mean_vals + variance_vals)))
}