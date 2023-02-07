######################################################################################################################

# GxE using marker-by-environment interactions

# (3) Incomplete, unbalanced designs

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
    base_folder = '/home/freeclimb/GxE',
    y = 'data/phenotypes.csv',
    X = 'data/markers.csv',
    evd_threshold = 1e-3,
    nIter = 2000,
    burnIn = 100,
    thin = 5, ## default value in BGLR is 5
    outdir = 'Analysis/BGLR/kinship',
    prefix = "GxE_kin_",
    subsample = 1000, ## n. of SNPs to subsample randomly
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
X <- as.data.frame(X)
vec <- sample(1:ncol(X), config$subsample)
X <- X[,vec]
row.names(Y) <- paste("s",1:nrow(Y),sep="")
row.names(X) <- paste("s",1:nrow(Y),sep="")

ntraits = ncol(Y) 

y <- Y |> rownames_to_column('id') |> gather(key = "trait", value = "value", -id)
# y <- y |> group_by(trait) |> sample_n(sample(nrow(Y),1))
y <- filter(y, !is.na(value))
v <- group_by(y, trait) |> summarise(N=n()) |> pull(N)

print(kable(v))

env = NULL
for (i in 1:length(v)) {
  print(i)
  env = c(env, rep(i,v[i]))
}

### KINSHIP MATRIX
writeLines(" - calculating matrix of genetic kinship")
G = getG(as.matrix(X),center=TRUE,scaleG=TRUE,scale=TRUE)

writeLines(" - extracting eigenvalues and eigenvectors from the kinship matrix")
EVD = eigen(G)
PC = EVD$vectors[,EVD$values > config$evd_threshold]   # 599 x 598 ? SELECT THRESHOLD TO RETAIN PC's
for(i in 1:ncol(PC)){ PC[,i]=EVD$vectors[,i]*sqrt(EVD$values[i]) }  #? eigenvectors multiplied by the sqrt of the corresponding eigenvalues

writeLines(" - preparing design matrices for multiple envs")
## make "p" copies of PC <-- "p" phenotypes (same pheno in "p" envs, same sample IDs: this is why it is duplicated)
row.names(PC) <- row.names(X)
X0 = PC[y$id,] # Matrix for main effects 

# now interactions
for (i in 1:ntraits) {
  n = ntraits
  print(paste(i,n))
  assign(paste("X", i, sep = ""), X0)
}

###############################
## !! MANUAL EDITING HERE !! ##
###############################
for(i in 1:nrow(X0)){
  X1[i,]<-(env[i]==1)*X0[i,]
  X2[i,]<-(env[i]==2)*X0[i,]	
  X3[i,]<-(env[i]==3)*X0[i,]
  X4[i,]<-(env[i]==4)*X0[i,]
}

LP = list(main=list(X=X0,model='BRR'), 
          int1=list(X=X1,model='BRR'),
          int2=list(X=X2,model='BRR'),
          int3=list(X=X3,model='BRR'),
          int4=list(X=X4,model='BRR')
)
###############################

print("Running the BGLR model - kinship matrix")
dir.create(file.path(config$base_folder, config$outdir), recursive = TRUE, showWarnings = FALSE)
outpath = file.path(config$base_folder, config$outdir, config$prefix)
fmGRM = BGLR(y=y$value,ETA=LP,
             nIter=config$nIter, burnIn=config$burnIn, thin = config$thin,
             saveAt=outpath, groups=env)

print("Writing out results")
fname = paste(config$prefix, "BRR_res.RData", sep="")
save(fmGRM, file = file.path(config$base_folder, config$outdir, fname))

print("DONE!")

# # Method: using SNPs explicitly (no kinship: SNP-BLUP instead of GBLUP)
# X0 = X[y$id,] # Matrix for main effects
# dim(X0)
# stopifnot(all(gsub("\\..*$","",row.names(X0))==y$id))
# 
# # now interactions
# for (i in 1:ntraits) {
#   n = ntraits
#   print(paste(i,n))
#   assign(paste("X", i, sep = ""), X0)
# }
# 
# ###############################
# ## !! MANUAL EDITING HERE !! ##
# ###############################
# for(i in 1:nrow(X0)) {
#   
#   X1[i,]<-(env[i]==1)*X0[i,]
#   X2[i,]<-(env[i]==2)*X0[i,]	
#   X3[i,]<-(env[i]==3)*X0[i,]	
#   X4[i,]<-(env[i]==4)*X0[i,]	
# }
# 
# ###############################
# 
# fm=BGLR(y=y$value,ETA=list(list(X=X0,model='BRR'),
#                      list(X=X1,model='BRR'),
#                      list(X=X2,model='BRR')),
#         groups=env)
# 
# fm$varE
# fm$ETA[[1]]$varB
# fm$ETA[[2]]$varB
# fm$ETA[[3]]$varB
