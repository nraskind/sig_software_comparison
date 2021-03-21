# subsamples a trinucleotide mutation matrix
# returns a mutation matrix with numvars mutations
subsample_mut_matrix <- function(tri_mut_matrix, numvars=25, replacement=TRUE) {
  # creates a vector of the 96 trinucleotides and their frequencies for sampling
  samplevec <- sapply(rownames(tri_mut_matrix), function(trinucleotide){
    rep(trinucleotide, tri_mut_matrix[trinucleotide, 1])
  })
  samplevec <- unlist(samplevec)
  names(samplevec) <- NULL

  # samples numvars mutations with/without replacement  
  sampled_mutations <- as.data.frame(table(sample(samplevec, numvars, replace=replacement)))

  # recreates the original trinucleotide matrix, this time with numvars mutations
  sub_sigs_input <- sapply(rownames(tri_mut_matrix), function(trinucleotide) {
    if(trinucleotide %in% as.character(sampled_mutations[,1])) {
      row <- which(as.character(sampled_mutations[,1]) == trinucleotide) 
      return(sampled_mutations[row, 2])
    } else {
       return(0)
    }
  })

  # matrix with sample name = column name
  sub_sigs_input <- as.matrix(sub_sigs_input)
  colnames(sub_sigs_input) <- colnames(tri_mut_matrix)
  
  return(sub_sigs_input)
}
