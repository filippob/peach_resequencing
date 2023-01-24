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
    outdir = 'Analysis/BGLR',
    force_overwrite = FALSE
  ))
  
}

# SETUP -------------------------------------------------------------------
library("BGLR")
library("knitr")
library("data.table")


# READ THE DATA -------------------------------------------------------------------
## phenotypes (same phenotypes, multiple envs (maybe also multiple years))
fname = file.path(config$base_folder, config$y) 
Y = fread(fname, header = TRUE) # grain yield evaluated in 4 different environments (599 samples, 4 columns)
## markers
fname = file.path(config$base_folder, config$X) 
X = fread(fname, header = TRUE) ## 599 samples, 1279 markers (DArT markers --> 0/1)


print("raw correlations between phenotypes")
print(kable(round(cor(Y),2))) # correlation between yields in the 4 environments (round with 2 decimals)

writeLines(" - scaling marker data")
X = scale(X)/sqrt(ncol(X))   # scaled genotypes (599 samples x 1279 markers) ## why /sqrt(ncols) ??


y2=Y[,2]  # 599 yield in env 2
y3=Y[,3]  # 599 yield in env 3
y=c(y2,y3) # yields in env 2 and 3 together (1198)

X0=matrix(nrow=nrow(X),ncol=ncol(X),0) # a matrix full of zeros (599 x 1279)

X_main=rbind(X,X) # (1198 x 1279) 2 copies of X  ?
X_1=rbind(X,X0)   # (1198 x 1279) X + X0  ?
X_2=rbind(X0,X)   # (1198 x 1279) X0 + X  ?

fm=BGLR( y=y,ETA=list(             
  main=list(X=X_main,model='BRR'),
  int1=list(X=X_1,model='BRR'),
  int2=list(X=X_2,model='BRR')
),
nIter=6000,burnIn=1000,saveAt='GxE_',groups=rep(1:2,each=nrow(X))
)

varU_main=scan('GxE_ETA_main_varB.dat')[-c(1:200)] #1000 ?
varU_int1=scan('GxE_ETA_int1_varB.dat')[-c(1:200)] #1000 ?
varU_int2=scan('GxE_ETA_int2_varB.dat')[-c(1:200)] #1000 ?
varE=read.table('GxE_varE.dat',header=FALSE)[-c(1:200),]  # 1000, from 201 to 1200 ?
varU1=varU_main+varU_int1
varU2=varU_main+varU_int2
h2_1=varU1/(varU1+varE[,1])
h2_2=varU2/(varU2+varE[,2])
COR=varU_main/sqrt(varU1*varU2)
mean(h2_1)
mean(h2_2)
mean(COR)







