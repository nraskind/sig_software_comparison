# siglasso
library("devtools")
library("siglasso")

#deconstructSigs
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)

# SignatureEstimation (qp)
library(SignatureEstimation)

# generatng trinucleotide context
library(MutationalPatterns)

# reference genome
bsg38 <- BSgenome.Hsapiens.UCSC.hg38

# 5% cutoff for mutational signature intensity - used for qp output
trim_low_intesity_sigs <- function(mut_sigs_matrix) {
  trimmed_sigs <- sapply(1:nrow(mut_sigs_matrix), function(signature){
    if (mut_sigs_matrix[signature, 1] < .05) {
      return(0)
    } else {
      return(mut_sigs_matrix[signature, 1])
    }
  })
  trimmed_sigs <- data.frame(trimmed_sigs,
                             row.names = rownames(mut_sigs_matrix))
  colnames(trimmed_sigs) <- colnames(mut_sigs_matrix)
  return(trimmed_sigs)
}

# creates a trinucleotide context matrix from (a) vcf file
create_mut_matrix <- function(vcf_files_input, reference_genome=bsg38) {
  vcf_list <- MutationalPatterns::read_vcfs_as_granges(
    vcf_files=vcf_files_input,
    sample_names=basename(vcf_files_input),
    genome=reference_genome
  )

  mut_matrix <- MutationalPatterns::mut_matrix(vcf_list, ref_genome=reference_genome, extension=1)
  return(mut_matrix)
}
mut_matrix <- create_mut_matrix(vcf_file_input)

# creates a data frame of mutational signatures for the given software
create_signature_df <- function(software, tri_mut_matrix) {
  if (missing(tri_mut_matrix)) {
    stop("Please enter a 96 x 1 trinucleotide matrix")
  }
  
  if (software == "sl") {
    sl_mut_sigs <- siglasso::siglasso(tri_mut_matrix,
                                      plot=FALSE,
                                      signature=t(signatures.cosmic))
    sl_mut_sigs <- as.data.frame(sl_mut_sigs)
    colnames(sl_mut_sigs) <- paste("sl", colnames(sl_mut_sigs), sep=".")
    return (sl_mut_sigs)
  } else if (software == "ds") {
    tumor_ref = as.data.frame(t(tri_mut_matrix))
    sample_id = rownames(tumor_ref)

    signature_list <- deconstructSigs::whichSignatures(
                                      tumor.ref = tumor_ref, 
                                      signatures.ref = signatures.cosmic, 
                                      sample.id = sample_id, 
                                      contexts.needed = TRUE,
                                      tri.counts.method = 'default')
    
    ds_mut_sigs <- as.data.frame(t(signature_list[["weights"]][]))
    colnames(ds_mut_sigs) <- paste("ds", colnames(ds_mut_sigs), sep=".")
    return(ds_mut_sigs)
  } else if (software == "qp") {
    qp_mut_sigs <- as.data.frame( SignatureEstimation::decomposeQP(m = tri_mut_matrix,
                                                                   P = t(signatures.cosmic)),
                                  row.names = rownames(signatures.cosmic))
    colnames(qp_mut_sigs) <- paste("qp", colnames(tri_mut_matrix), sep=".")
    qp_mut_sigs <- trim_low_intesity_sigs(qp_mut_sigs)
    return(qp_mut_sigs)
  } else {
    stop("Please enter one of: \"sl\", \"ds\", or \"qp\" to be used in generating the mutational signatures")
  }
}

# "ground truth" dataframes
ds_mut_df <- create_signature_df("ds", tri_mut_matrix=mut_matrix)
sl_mut_df <- create_signature_df("sl", tri_mut_matrix=mut_matrix)
qp_mut_df <- create_signature_df("qp", tri_mut_matrix=mut_matrix)

# create subsamples with 5, 10, 25, 50, and 100 variants
# simulation_repeats (25) simulations/columns each
create_subsample_dfs <- function(numvars, mut_matrix, simulation_repeats=25) {
  num_subs <- as.data.frame(sapply(1:simulation_repeats, function(x){
    subsample_mut_matrix(mut_matrix, numvars)
  }))
  colnames(num_subs) <- rep(colnames(mut_matrix), ncol(num_subs))
  rownames(num_subs) <- rownames(mut_matrix)
  return(num_subs)
}

five_subs <- create_subsample_dfs(5, mut_matrix)
ten_subs <- create_subsample_dfs(10, mut_matrix)
twenty_five_subs <- create_subsample_dfs(25, mut_matrix)
fifty_subs <- create_subsample_dfs(50, mut_matrix)
hundred_subs <- create_subsample_dfs(100, mut_matrix)

# generates a data frame of mutational signatures for every column of
# the above subsample df's with the specified software
software_subsample_dfs <- function(software, subsample_df) {
  if (software == "sl") {
    # right now siglasso is the only one that can take multiple samples at once
    mut_sigs_output <- create_signature_df("sl", subsample_df)
  }else {
   mut_sigs_output <- as.data.frame(sapply(1:ncol(subsample_df), function(x){
     create_signature_df(software, subsample_df[,x,drop=FALSE])
    }))
  }
  
  return(mut_sigs_output)
}
