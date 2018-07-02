# Copyright (c) 2018 - 2020 Chong Ma
#
# This file contains the utility functions for the WAP R package. The WAP R package
# is written for the weight adjusted p-value approach, which enables to increase
# the power while controlling the familywise error rate for discovering the true 
# significant gene assosiated with the response variable. 

#' Projection Matrix
#' 
#' Calculate the projection of matrix
#' 
#' @param X matrix with dimensions \eqn{n \times p}.
#' @param ... additional matrix or vectors that have same column 
#'            dimension with X.
#'            
#' @return A diagonal matrix.
#' 
#' @examples 
#' 
#' \dontrun{
#'   X=matrix(10,rnorm(100))
#'   projM(X)
#'   
#'   # add additional vector Y
#'   Y=rep(1,10)
#'   projM(X,Y) 
#' }
#' 
#' @export
projM=function(X,...){
  X=cbind(X,...)
  
  if(!is.null(X)){
    u=svd(X)$u
    return(u%*%t(u))
  }else{
    return(0)
  }
}

#' Weighted p-value 
#' 
#' Calcuate adjusted p-value by using weighted p-value approach
#' 
#' @param Y A vector of numerical values. response variable.
#' @param X A design matrix. 
#' @param Z A matrix for intermediate variables.
#' @param V A matrix for external variables. Default is NULL.
#' @param subset A vector of indices for a sub-sample of the data or a subsample size. 
#'               subset is TRUE for method ``Lw". Default is NULL.
#' @param Raw.Weight Whether or not return the raw weight. Default is FALSE.
#' @param method weight methods. Must choose one from ``mean",``max", and ``Lw". 
#' @param ncores number of cpus. See \code{\link[parallel]{makeCluster}}.
#' @param center Whether or not centerize the data. Default is FALSE.
#' @param scale Whether or not scale the data. Default is FALSE. 
#' 
#' @return A dataframe consists of raw p-values, weight, and adjusted p-values. When the Raw.Weight
#'         is TRUE, a list will be returned with the result dataframe and the raw weight. 
#' @examples 
#' 
#' \dontrun{
#' X=matrix(rnorm(500),50)
#' beta1=0.5;beta2=0.1
#' Z1=X[,1]*beta1+rnorm(50)
#' Y=Z1*beta2+rnorm(50)
#' Z=cbind(Z1,matrix(rnorm(50000),50))
#' test=wpvalue(Y,X,Z,method="max",ncores=4)
#' } 
#'
#' @export
wpvalue<-function(Y,X,Z,V=NULL,subset=NULL,Raw.Weight=FALSE,
                  method=c("max","mean","Lw"),
                  ncores=1,center=TRUE,scale=FALSE){
  
  # whether or not scale or center data
  N=length(Y)
  Y=scale(Y,center,scale)
  X=scale(X,center,scale)
  Z=scale(Z,center,scale)
  if(!is.null(V)) V=scale(V,center,scale)
  
  p=ifelse(is.null(V),0,ncol(V)) ## ndim of V
  pvalue=weight=adj.pvalue=vector(mode = "numeric",length=ncol(X))
  raw.weight=vector(mode = "list",length = ncol(X))
  i=j=NULL
  
  # error checking for method
  if(N<p+3) stop("The design matrix is degenerate because of N<p!")
  if(all(method %in% c("max","mean","Lw")) && length(method)==3) method="max"
  if(length(method) > 1 || !(method %in% c("max","mean","Lw"))){
    stop("Must choose one method from max, mean or Lw!")
  } 
  
  cl<-parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  if(is.null(subset)){
    cnst=rep(1,N)
    for(i in 1:ncol(X)){
      F_XY=(N-p-2)*t(Y)%*%(projM(cnst,V,X[,i])-projM(cnst,V))%*%Y/
        (t(Y)%*%(diag(N)-projM(cnst,V,X[,i]))%*%Y)
      
      raw.weight[[i]]=foreach::foreach(j=1:ncol(Z),.combine = rbind,
                                       .multicombine = TRUE,
                                       .errorhandling = "pass",
                                       .export = c("projM")) %dopar%{
                                         F_XZ=(N-p-3)*t(Z[,j])%*%(projM(cnst,V,Y,X[,i])-projM(cnst,V,Y))%*%Z[,j]/
                                           (t(Z[,j])%*%(diag(N)-projM(cnst,V,Y,X[,i]))%*%Z[,j])
                                         F_ZY=(N-p-3)*t(Y)%*%(projM(cnst,V,X[,i],Z[,j])-projM(cnst,V,X[,i]))%*%Y/
                                           (t(Y)%*%(diag(N)-projM(cnst,V,X[,i],Z[,j]))%*%Y)
                                         
                                         as.vector(c(F_XZ,F_ZY))
                                       }
      pvalue[i]=stats::pf(F_XY,1,N-p-2,lower.tail = FALSE)
    }
  }else{
    # one subset to calculate p-value, and
    # the other to calculate weight
    if(length(subset)>1){
      Nl=length(subset);Nt=N-Nl
      Yl=Y[subset];Yt=Y[-subset]
      Xl=X[subset,];Xt=X[-subset,]
      Zl=Z[subset,];Zt=Z[-subset,]
      if(!is.null(V)){
        Vl=V[subset,];Vt[-subset,]
      }
    }else{
      Nl=subset;Nt=N-Nl
      subsam=sample(1:N,Nl)
      Yl=Y[subsam];Yt=Y[-subsam]
      Xl=X[subsam,];Xt=X[-subsam,]
      Zl=Z[subsam,];Zt=Z[-subsam,]
      if(is.null(V)){
        Vl=Vt=NULL
      }else{
        Vl=V[subsam,];Vt[-subsam,]
      }
    }
    
    if(Nl<p+3) stop("The design matrix is degenerate because of N<p!")
    cnstl=rep(1,Nl);cnstt=rep(1,Nt)
    for(i in 1:ncol(X)){
      F_XY=(Nl-p-2)*t(Yl)%*%(projM(cnstl,Vl,Xl[,i])-projM(cnstl,Vl))%*%Yl/
        (t(Yl)%*%(diag(Nl)-projM(cnstl,Vl,Xl[,i]))%*%Yl)
      
      raw.weight[[i]]=foreach::foreach(j=1:ncol(Zt),.combine = rbind,
                                       .multicombine = TRUE,
                                       .errorhandling = "pass",
                                       .export = c("projM")) %dopar%{
                                         F_XZ=(Nt-p-3)*t(Zt[,j])%*%(projM(cnstt,Vt,Yt,Xt[,i])-projM(cnstt,Vt,Yt))%*%Zt[,j]/
                                           (t(Zt[,j])%*%(diag(Nt)-projM(cnstt,Vt,Yt,Xt[,i]))%*%Zt[,j])
                                         F_ZY=(Nt-p-3)*t(Yt)%*%(projM(cnstt,Vt,Xt[,i],Zt[,j])-projM(cnstt,Vt,Xt[,i]))%*%Yt/
                                           (t(Yt)%*%(diag(Nt)-projM(cnstt,Vt,Xt[,i],Zt[,j]))%*%Yt)
                                         
                                         as.vector(c(F_XZ,F_ZY))
                                       }
      pvalue[i]=stats::pf(F_XY,1,Nt-p-2,lower.tail = FALSE)
    }
  }
  parallel::stopCluster(cl)
  
  if(method %in% "max"){
    tmp.weight=sapply(raw.weight,function(t) max(t[,1]*t[,2],na.rm = TRUE))
    weight=tmp.weight/mean(tmp.weight)
  }
  
  if(method %in% "mean"){
    tmp.weight=sapply(raw.weight,function(t) mean(t[,1]*t[,2],na.rm = TRUE))
    weight=tmp.weight/mean(tmp.weight)
  }
  
  if(method %in% "Lw"){
    tmp.weight=sapply(raw.weight,function(t) max(t[,1],na.rm = TRUE))
    weight=tmp.weight/mean(tmp.weight)
  }
  
  adj.pvalue=sapply(pvalue/weight,function(t) min(t,1))
  
  if(Raw.Weight){
    return(list(wapTab=data.frame(p.value=pvalue,
                                  weight=weight,
                                  adj.pvalue=adj.pvalue),
                raw.weight=raw.weight))
  }else{
    return(wapTabe=data.frame(p.value=pvalue,
                              weight=weight,
                              adj.pvalue=adj.pvalue))
  }
}


