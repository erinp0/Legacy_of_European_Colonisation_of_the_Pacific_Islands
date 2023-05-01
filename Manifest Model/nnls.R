#########################################
## BEGIN PROGRAM:

options(scipen=999)
library(nnls)

## General NNLS functions

admix.nnls<-function(X,Y){
    ## Do an nnls admixture analysis (X ~ Y) where X is a vector
##    if(is(X,"numeric")) stop("X must be numeric in admix.nnls")
    ourmix=getoverallfit(Y,X)$x
    ourmix=ourmix/sum(ourmix)
    ourmix
}

admix.nnls.all<-function(X,Y,verbose=TRUE){
    ## Do nnls admixture on each ROW of the MATRIX x
    ##if(is(X,"matrix")) stop("X must be a matrix in admix.nnls.all")
    
    ret<-t(sapply(1:dim(X)[1],function(i){
        if(verbose) print(paste("Processing ind",i,"of",dim(X)[1]))
        admix.nnls(X[i,],Y)
    }))
    rownames(ret)<-rownames(X)
    ret
}


#######################
## MAIN NNLS FUNCTIONS, from GLOBETROTTER paper

getfit=function(predmat,fitdata,restrict=1){
    ## Gets the nnls fit with the restriction that a given row is not used.
    ## This addresses the sum to one constraint.
    ## If we succeed in getting a valid return, it is also guarenteed to be the
    ## best fitting. Also, we are guaranteed to get this from some restriction
    ## choice
  temp=predmat[-restrict,,drop=FALSE]
  for(i in 1:nrow(temp)) temp[i,]=temp[i,]-predmat[restrict,]

  fitdata2=fitdata-predmat[restrict,]

  v=nnls(t(temp),fitdata2)
  
  x=v$x
  newx=1:nrow(predmat)
  newx[!((1:nrow(predmat))==restrict)]=x
  newx[restrict]=1-sum(x)
  v$x=newx
  names(v$x)=rownames(predmat)

  return(v)
}

getoverallfit=function(predmat,fitdata){
  restrict=1
  rep=1
  i=1
  while(rep==1){
    q=getfit(predmat,fitdata,restrict=i)

    if(q$x[i]>0) rep=0
    i=i+1
  }

  return(q)
}

##############################################
## COMMAND TO RUN:

#ourmix=getoverallfit(donor.mat,recipient.vec)$x
#ourmix=ourmix[ourmix>0]
#ourmix=ourmix/sum(ourmix)

## END PROGRAM
###############################################
