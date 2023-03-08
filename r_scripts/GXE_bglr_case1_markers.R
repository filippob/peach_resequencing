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
    base_folder = '/home/freeclimb/GxE',
    y = 'data/BD_2019/phenotypes.csv',
    X = 'data/BD_2019/markers.csv',
    trait = 'BD',
    year = 2019,
    nIter = 2000,
    burnIn = 200,
    thin = 5, ## default value in BGLR is 5
    outdir = 'Analysis/BGLR',
    prefix = "GxE_mrk_reg_",
    subsample = 1000, ## n. of SNPs to subsample randomly
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
X <- as.data.frame(X)
vec <- sample(1:ncol(X), config$subsample)
X <- X[,vec]
row.names(Y) <- paste("s",1:nrow(Y),sep="")
row.names(X) <- paste("s",1:nrow(Y),sep="")

ntraits = ncol(Y) 

print("raw correlations between phenotypes")
print(kable(round(cor(Y, use = "complete.obs"),2))) # correlation between yields in the 4 environments (round with 2 decimals)

writeLines(" - scaling marker data")
X = scale(X, center = TRUE, scale = FALSE)/sqrt(ncol(X))   # scaled genotypes (599 samples x 1279 markers) ## why /sqrt(ncols) ??

writeLines(" - preparing long vector of phenotypes")
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

################################################
print("Running the BGLR model - marker matrix")
experiment = paste(config$trait, config$year, sep="_")
dir.create(file.path(config$base_folder, config$outdir, experiment), recursive = TRUE, showWarnings = FALSE)


#############################################
## !! MANUAL EDITING HERE IF NTRAITS > 4!! ##
#############################################
if (ntraits == 2) {
  
  LP = list(main=list(X=X_main,model='BRR'), 
            int1=list(X=X_1,model='BRR'),
            int2=list(X=X_2,model='BRR')
  )
} else if (ntraits == 3) {
  
  LP = list(main=list(X=X_main,model='BRR'), 
            int1=list(X=X_1,model='BRR'),
            int2=list(X=X_2,model='BRR'),
            int3=list(X=X_3,model='BRR')
  )
} else if (ntraits == 4) {
  
  LP = list(main=list(X=X_main,model='BRR'), 
            int1=list(X=X_1,model='BRR'),
            int2=list(X=X_2,model='BRR'),
            int3=list(X=X_3,model='BRR'),
            int4=list(X=X_4,model='BRR')
  )
}
###############################


outpath = file.path(config$base_folder, config$outdir, experiment, config$prefix)
fm = BGLR(y=y$value,ETA=LP,
  nIter=config$nIter, burnIn=config$burnIn, thin = config$thin,
  saveAt=outpath, groups=rep(1:ntraits,each=nrow(X))
)

print("Writing out results")
fname = paste(config$prefix, "BRR_res.RData", sep="")
save(fm, file = file.path(config$base_folder, config$outdir, experiment, fname))

print("DONE!")

## stub of code to plot marker effects (not p-values)
# bHat = fm$ETA[[1]]$b ## the index [[1]] refers to the common genetic effect from the model
# plot(bHat^2, ylab='Estimated Squared-Marker Effect',
#      type='o',cex=.5,col=4,main='Marker Effects')
