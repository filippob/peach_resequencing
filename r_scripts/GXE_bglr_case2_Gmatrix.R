######################################################################################################################

# GxE using marker-by-environment interactions


# (2) Using genomic relationships

#A model equivalent to the one presented above can be implemented using G-matrices (or factorizations of it) 
# with off-diagonal blocks zeroed out for interactions

######################################################################################################################

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
    evd_threshold = 1e-4,
    nIter = 3000,
    burnIn = 500,
    thin = 5, ## default value in BGLR is 5
    outdir = 'Analysis/BGLR/kinship',
    prefix = "GxE_kin_",
    force_overwrite = FALSE
  ))
  
}

# SETUP -------------------------------------------------------------------
library("BGLR")
library("knitr")
library("BGData")
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

####################
writeLines(" - preparing long vector of phenotypes")
y <- Y |> gather(key = "trait")

writeLines(" - calculating matrix of genetic kinship")
G = getG(as.matrix(X),center=TRUE,scaleG=TRUE,scale=TRUE)
# G=getG(wheat.X,center=TRUE,scaleG=TRUE,scale=TRUE, nCores = getOption("mc.cores", 1L))
# If center = FALSE, scale = FALSE and scaleG = FALSE, getG produces the same outcome than tcrossprod

writeLines(" - extracting eigenvalues and eigenvectors from the kinship matrix")
EVD = eigen(G)
PC = EVD$vectors[,EVD$values > config$evd_threshold]   # 599 x 598 ? SELECT THRESHOLD TO RETAIN PC's
for(i in 1:ncol(PC)){ PC[,i]=EVD$vectors[,i]*sqrt(EVD$values[i]) }  #? eigenvectors multiplied by the sqrt of the corresponding eigenvalues

writeLines(" - preparing design matrices for multiple envs")
## make "p" copies of PC <-- "p" phenotypes (same pheno in "p" envs, same sample IDs: this is why it is duplicated)
XMain = matrix(rep(t(PC), ntraits), ncol=ncol(PC),byrow=TRUE) 


X0 = matrix(nrow=nrow(X),ncol=ncol(PC),0) # a matrix full of zeros
for (i in 1:ntraits) {
  
  n = ntraits
  print(paste(i,n))
  temp1 = matrix(rep(t(X0), i-1), ncol=ncol(X0),byrow=TRUE)
  temp2 = matrix(rep(t(X0), n-i), ncol=ncol(X0),byrow=TRUE)
  temp = rbind(temp1, PC, temp2)
  assign(paste("X", i, sep = "_"), temp)
}

###########################
# BGLR MODEL -------------------------------------------------------------------
print("Running the BGLR model - kinship matrix")
dir.create(file.path(config$base_folder, config$outdir), recursive = TRUE, showWarnings = FALSE)

outpath = file.path(config$base_folder, config$outdir, config$prefix)
LP = list(main=list(X=XMain,model='BRR'), 
          int1=list(X=X_1,model='BRR'),
          int2=list(X=X_2,model='BRR'),
          int3=list(X=X_3,model='BRR'),
          int4=list(X=X_4,model='BRR')
      )

fmGRM = BGLR(y=y$value,ETA=LP,
             nIter=config$nIter, burnIn=config$burnIn, thin = config$thin,
             saveAt=outpath, groups=rep(1:ntraits,each=nrow(X))
             )

print("Writing out results")
fname = paste(config$prefix, "BRR_res.RData", sep="")
save(fm, file = file.path(config$base_folder, config$outdir, fname))

print("DONE!")



