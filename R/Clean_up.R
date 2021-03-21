# temporary file clean up
temp_files <- list.files(temp_dir, full.names = T, pattern = "vcf")
file.remove(temp_files)
