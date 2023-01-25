
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
    outdir = 'Analysis/BGLR',
    prefix = "GxE_",
    ntraits = 4,
    burnIn = 500,
    thin = 5, ## default thin in BGLR is 5
    force_overwrite = FALSE
  ))
  
}

#### READ FILES
writeLines(" - reading results files")

varU_main = scan(file.path(config$base_folder, config$outdir, paste(config$prefix, 'ETA_main_varB.dat', sep="")))
varU_main = varU_main[-c(1:(config$burnIn/config$thin))]

for (i in 1:config$ntraits) {
  
  fname = paste(config$prefix, 'ETA_int', i, '_varB.dat', sep="")
  temp = scan(file.path(config$base_folder, config$outdir, fname))
  assign(paste("varU_int", i, sep = ""), temp)
}

varE = read.table(file.path(config$base_folder, config$outdir,paste(config$prefix, 'varE.dat', sep="")),header=FALSE)[-c(1:200),]  # 1000, from 201 to 1200 ?

## PROCESS FILES
## hic sunt leones
varU1=varU_main+varU_int1
varU2=varU_main+varU_int2
varU3=varU_main+varU_int3
varU4=varU_main+varU_int4

h2_1=varU1/(varU1+varE[,1])
h2_2=varU2/(varU2+varE[,2])
h2_3=varU3/(varU3+varE[,3])
h2_4=varU4/(varU4+varE[,4])

m <- cbind.data.frame(varU1, varU2, varU3, varU4)
# sqrt(m)
COR=varU_main/sqrt(varU1*varU2*varU3*varU4)
mean(h2_1)
mean(h2_2)
mean(COR)
