#escribir bedGraphs para datos con zScore

my_trackPath <- "tracks_minus0", 

writeTrackFiles(object = fcf, 
                assay = "zScore", 
                folder = my_trackPath,
                format = "bedGraph", 
                removeZeros = FALSE)

#list bedGraphs
PATH_TO_DATA <- paste0(metadata(fcf)$projectPath, my_trackPath)

my_samples <- list.files(path = PATH_TO_DATA, pattern = ".bedGraph")


lapply(X = my_samples, FUN = function(i){
  my_sample  <- i
  my_path    <- paste0(PATH_TO_DATA, i)
  my_tbl     <- data.table::fread(input = my_path)
  my_logical <- assay(fcf, "peaks")[,my_sample]
  data.table::fwrite(x = testing[my_logical], sep = "\t", 
                     file = paste0(PATH_TO_DATA, "clean_", my_sample, ".clean.bedGraph")
  )
})

my_sample = "testdata_hf_zt6_1"
my_path   = paste0("DataPPARG/tracks_minus0/zScore_", my_sample, ".bedGraph")
testing <-data.table::fread(input = my_path)
testing



data.table::fwrite(x = testing[my_logical], sep = "\t", 
                   file = paste0("DataPPARG/tracks_minus0/zScore_Clean_", my_sample, ".clean.bedGraph")
                   )


my_sample = "testdata_hf_zt6_1"
my_path   = paste0("DataPPARG/tracks/counts_", my_sample, ".bedGraph")
testing <-data.table::fread(input = my_path)
testing
