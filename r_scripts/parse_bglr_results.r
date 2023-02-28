# INPUT CONFIGURATION MANAGEMENT ------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 1){
  #loading the parameters
  source(args[1])
} else {
  #this is the default configuration, used for development and debug
  writeLines('Using default config')
  
  #this dataframe should be always present in config files, and declared
  #as follows
  config = NULL
  config = rbind(config, data.frame(
    base_folder = '~/Documents/freeclimb/g_x_e',
    X = 'data/MD_2019/markers.csv',
    result_folder = 'Analysis/BGLR',
    outdir = 'Analysis/BGLR/MD_2019',
    traits = 'MD_2019,MD_2020',
    prefix = "GxE_kinship_incomplete",
    ntraits = 3,
    burnIn = 100,
    thin = 5, ## default thin in BGLR is 5
    force_overwrite = FALSE
  ))
  
}

library("stringr")
library("data.table")

traits = strsplit(config$traits, ",")[[1]]

res = data.frame("trait"=NULL, "test_corr"=NULL, "train_corr"=NULL)
for (trt in traits) {
  
  trt = str_trim(trt, side = "both")
  fpath = file.path(config$base_folder, config$result_folder, trt)
  filenames = list.files(fpath, pattern="*.RData", full.names=TRUE)
  fname = filenames[grepl(pattern = config$prefix, filenames)][1]
  print(paste("Loading data from",fname))
  load(fname)
  
  y <- to_save$y
  fmGRM <- to_save$fmGRM
  tst = to_save$test_rows
  
  #5# Assesment of correlation in TRN and TST data sets
  test_corr = cor(fmGRM$yHat[tst],y$value[tst])
  train_corr = cor(fmGRM$yHat[-tst],y$value[-tst]) #TRN
  
  temp = data.frame("trait"=trt, "test_corr"=test_corr, "train_corr"=train_corr)
  res = rbind.data.frame(res,temp)
  
  rm(y,fmGRM,tst)
}

print("DONE!!")

