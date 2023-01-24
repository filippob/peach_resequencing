######################################################################################################################

# GxE using marker-by-environment interactions


# (3) Incomplete, unbalanced designs

library(BGLR)
data(wheat)
Y=wheat.Y
X=scale(wheat.X)/sqrt(ncol(wheat.X))
rownames(X)<-1:nrow(X)# use your IDs, don't need to be 1:n
rownames(Y)<-1:nrow(Y) # Y has data from 4 env (each col=1 env) 
n1=150
n2=305
y=c(sample(Y[,1],size=n1),sample(Y[,2],size=n2)) #unbalanced data from 2 env
env=c(rep(1,n1),rep(2,n2))

# Method 1: using SNPs explicitly

X0=X[names(y),] # Matrix for main effects
stopifnot(all(rownames(X0)==names(y)))

# now interactions
X1=X0
X2=X0
for(i in 1:nrow(X0)){
  X1[i,]<-(env[i]==1)*X0[i,]
  X2[i,]<-(env[i]==2)*X0[i,]	
}
fm=BGLR(y=y,ETA=list(list(X=X0,model='BRR'),
                     list(X=X1,model='BRR'),
                     list(X=X2,model='BRR'))
        ,groups=env)

fm$varE
fm$ETA[[1]]$varB
fm$ETA[[2]]$varB
fm$ETA[[3]]$varB


