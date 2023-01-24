######################################################################################################################

# GxE using marker-by-environment interactions

# (1) As a random regression on markers

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
    y = 'data/phenotypes.csv',
    X = 'data/markers.csv',
    nIter = 3000,
    burnIn = 500,
    outdir = 'Analysis/BGLR',
    force_overwrite = FALSE
  ))
  
}

# SETUP -------------------------------------------------------------------
library("BGLR")
library("knitr")
library("tidyverse")
library("data.table")


# READ THE DATA -------------------------------------------------------------------
## phenotypes (same phenotypes, multiple envs (maybe also multiple years))
fname = file.path(config$base_folder, config$y) 
Y = fread(fname, header = TRUE) # grain yield evaluated in 4 different environments (599 samples, 4 columns)
## markers
fname = file.path(config$base_folder, config$X) 
X = fread(fname, header = TRUE) ## 599 samples, 1279 markers (DArT markers --> 0/1)

ntraits = ncol(Y) 

print("raw correlations between phenotypes")
print(kable(round(cor(Y),2))) # correlation between yields in the 4 environments (round with 2 decimals)

writeLines(" - scaling marker data")
X = scale(X)/sqrt(ncol(X))   # scaled genotypes (599 samples x 1279 markers) ## why /sqrt(ncols) ??


y <- Y |> gather(key = "trait")

## make a zero matrix the same size of the marker data
X0 = matrix(nrow=nrow(X),ncol=ncol(X),0) # a matrix full of zeros (599 x 1279)

## make "p" copies of X  <-- "p" phenotypes (same pheno in "p" envs, same sample IDs: this is why it is duplicated)
X_main <- matrix(rep(t(X), ntraits), ncol=ncol(X),byrow=TRUE) 
# X_main = rbind(X,X) # (1198 x 1279) 

for (i in 1:ntraits) {
  
  n = ntraits
  print(paste(i,n))
  temp1 = matrix(rep(t(X0), i-1), ncol=ncol(X0),byrow=TRUE)
  temp2 = matrix(rep(t(X0), n-i), ncol=ncol(X0),byrow=TRUE)
  temp = rbind(temp1, X, temp2)
  assign(paste("X", i, sep = "_"), temp)
}

print("Writing out example data")
dir.create(file.path(config$base_folder, config$outdir), recursive = TRUE, showWarnings = FALSE)

outpath = file.path(config$base_folder, config$outdir, "GxE_")
fm = BGLR(y=y$value,ETA=list(             
  main=list(X=X_main,model='BRR'),
  int1=list(X=X_1,model='BRR'),
  int2=list(X=X_2,model='BRR'),
  int3=list(X=X_3,model='BRR'),
  int4=list(X=X_4,model='BRR')
  ),
  nIter=config$nIter, burnIn=config$burnIn, saveAt=outpath, groups=rep(1:ntraits,each=nrow(X))
)

## hic sunt leones
varU_main=scan(file.path(config$base_folder, config$outdir,'GxE_ETA_main_varB.dat'))[-c(1:200)] #1000 ?
varU_int1=scan(file.path(config$base_folder, config$outdir,'GxE_ETA_int1_varB.dat'))[-c(1:200)] #1000 ?
varU_int2=scan(file.path(config$base_folder, config$outdir,'GxE_ETA_int2_varB.dat'))[-c(1:200)] #1000 ?
varU_int3=scan(file.path(config$base_folder, config$outdir,'GxE_ETA_int3_varB.dat'))[-c(1:200)] #1000 ?
varU_int4=scan(file.path(config$base_folder, config$outdir,'GxE_ETA_int4_varB.dat'))[-c(1:200)] #1000 ?

varE=read.table(file.path(config$base_folder, config$outdir,'GxE_varE.dat'),header=FALSE)[-c(1:200),]  # 1000, from 201 to 1200 ?
varU1=varU_main+varU_int1
varU2=varU_main+varU_int2
varU3=varU_main+varU_int3
varU4=varU_main+varU_int4

h2_1=varU1/(varU1+varE[,1])
h2_2=varU2/(varU2+varE[,2])
h2_3=varU3/(varU3+varE[,3])
h2_4=varU4/(varU4+varE[,4])

m <- cbind.data.frame(varU1, varU2, varU3, varU4)
sqrt(m)
COR=varU_main/sqrt(varU1*varU2*varU3*varU4)
mean(h2_1)
mean(h2_2)
mean(COR)







