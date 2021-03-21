library(VariantAnnotation)

temp_dir <- tempdir()

# filters a file, creates a .bgz, and .bgz.tbi, and stores them in temp_dir
filter_one_sample <- function(vcf_input_file){
  
  filter_params <- function(vcf_input_file){
    ## Getting rid of unknown chromosomes
    notUn <- !(grepl(pattern="chrUn", x=rowRanges(vcf_input_file)))
    
    tlod <- info(vcf_input_file)$TLOD > 12
    nlod <- info(vcf_input_file)$NLOD > 10
    
    pv2 <- info(vcf_input_file)$PV2 < .05
    pv2 <- !is.na(pv2)
    
    ## filter by AD
    ad <- geno(vcf_input_file)$AD
    
    # percentage filter
    ad_pct <- sapply(ad[,1], function(allele){
      # .05 or 1/200 cutoff
      (allele[2] / sum(allele[1], allele[2])) > .005
    })
    
    # constant filter
    ## checks if reference and variant allele frequencies are greater than 6
    ad_const <- min(as(ad[,1], "List")) > 6.0
    
    tlod & nlod & pv2 & notUn & ad_const & ad_pct
  }
  filter_input <- FilterRules(list(isSNV=isSNV, filters=filter_params))
  
  filter_destination <- paste(temp_dir, basename(vcf_input_file), sep="/")
  filter_destination_zipped <- paste(filter_destination, "bgz", sep=".")
  filter_destination_indexed <- paste(filter_destination_zipped, "tbi", sep=".")
  
  if(file.exists(filter_destination)) {
    file.remove(filter_destination)
  }
  
  if(file.exists(filter_destination_zipped)){
    file.remove(filter_destination_zipped)
  }
  
  if(file.exists(filter_destination_indexed)){
    file.remove(filter_destination_indexed)
  }
  
  #creates a bgzip and .tbi file
  filterVcf(file=vcf_input_file, 
            genome="hg38",
            destination=filter_destination,
            filters=filter_input,
            verbose=FALSE,
            index=TRUE)
  
  # creates a vcf file to accompany the bgzip and tabix files for siglasso
  # system(paste("bgzip -cd ", filter_destination_zipped, " > ", filter_destination, sep=""))
  
  ## return file path to the bgzipped vcf file
  return(filter_destination_zipped)
}