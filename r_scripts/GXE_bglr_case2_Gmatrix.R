######################################################################################################################

# GxE using marker-by-environment interactions


# (2) Using genomic relationships

#A model equivalent to the one presented above can be implemented using G-matrices (or factorizations of it) 
# with off-diagonal blocks zeroed out for interactions

library(BGLR)
data(wheat)
Y=wheat.Y # grain yield evaluated in 4 different environments
round(cor(Y),2)
y2=Y[,2]
y3=Y[,3]
y=c(y2,y3)

library(BGData)

G=getG(wheat.X,center=TRUE,scaleG=TRUE,scale=TRUE)
#Error in mclapply(X = seq_len(nChunks), FUN = chunkApply, mc.cores = nCores) : 
#'mc.cores' > 1 is not supported on Windows
G=getG(wheat.X,center=TRUE,scaleG=TRUE,scale=TRUE, nCores = getOption("mc.cores", 1L))
# If center = FALSE, scale = FALSE and scaleG = FALSE, getG produces the same outcome than tcrossprod

EVD=eigen(G)
PC=EVD$vectors[,EVD$values>1e-5]   # 599 x 598 ?
for(i in 1:ncol(PC)){ PC[,i]=EVD$vectors[,i]*sqrt(EVD$values[i]) }  #?

XMain=rbind(PC,PC)
X0=matrix(nrow=nrow(X),ncol=ncol(PC),0) # a matrix full of zeros
X1=rbind(PC,X0)
X2=rbind(X0,PC)

LP=list(main=list(X=XMain,model='BRR'), int1=list(X=X1,model='BRR'),int2=list(X=X2,model='BRR'))
fmGRM=BGLR(y=y,ETA=LP,nIter=12000,burnIn=2000,saveAt='GRM_',groups=rep(1:2,each=nrow(X)))

plot(fm$yHat,fmGRM$yHat)

rbind(fm$varE,fmGRM$varE)



