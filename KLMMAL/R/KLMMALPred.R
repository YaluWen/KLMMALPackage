#' @title Genetic Risk Prediction Using a Multi-Kernel Linear Mixed Model with Adaptive Lasso
#'
#' @description
#' \code{KLMMALPred} analyses the data prepared by \code{ReadKLMM} and predict the outcome values based on Multi-Kernel Linear Mixed Model with Adaptive Lasso method
#' @details
#' This function is used to predict phenotype using the Multi-Kernel Linear Mixed Model with Adaptive Lasso method. 
#'
#' @param Data The input data prepared by \code{ReadKLMM}. 
#' @param method Newton-Raphson(NR) or Expectation-Maximization(EM). Default is \code{NR}. 
#' @param lambdaall The vector of tuning parameters, and default is \code{NULL}.
#' @param stdK A bool indicating whether the GSMs should be standardized, and default is \code{TRUE}.  
#' @param tol The tolerance level, and default is \code{1e-6}.
#' @param tolK The tolerance level of SVD decomposition, and default is \code{1e-6}.
#' @param stepmax Maximum number of steps for the algorithms to run, and default is \code{300}.
#' @return A list that contains two elements (MLE and Pred). The Pred element stores the predicted values for the individuals in the testing samples (i.e. not in the \code{train.index}),
#' The MLE element contains \code{paras}: parameter estimates for each genomic region, \code{sigma2}: variance estimates,
#' \code{beta}: parameter estimates for demographic variables (i.e. \code{X}), \code{lambda}: the selected tuning parameters,
#' \code{allzero}: indicating whether a gBLUP model has been fitted (if \code{TRUE}, then a gBLUP model is fit.
#' @examples
#' ## We recommend to use "NR" unless your sample size is very small ##
#' ## EM can take a long time if the number of regions is large or the sample size is relatively moderate # 
#'  
#' ## Eg1: read both genotype and phenotype files from the disk, calculate the GSMs using linear kernels, and get the predicted values ##
#' trainID=read.table(system.file("extdata", "TrainID.txt", package = "KLMMAL"));
#' annotation=read.table(system.file("extdata", "Annotation.txt", package = "KLMMAL"),header=T);
#' Data=ReadKLMM(bed=system.file("extdata", "testing.bed", package = "KLMMAL"),trainID=trainID[,1], annotation = annotation, phenofile=system.file("extdata", "pheno.txt", package = "KLMMAL"), kernels = c("linear"), AllRegions = 0);
#' Result1=KLMMALPred(Data,method="NR");
#' print(cor(Result1$Pred,Data$Y[-Data$train.index]))
#' Result2=KLMMALPred(Data,method="EM", lambdaall=c(0,0.1,0.2)); # don't run, can take a very long time #
#' 
#' ## Eg2: read both genotype and phenotype files from the disk, and calculate the GSMs using both linear and poly2 kernels ##
#' trainID=read.table(system.file("extdata", "TrainID.txt", package = "KLMMAL"));
#' annotation=read.table(system.file("extdata", "Annotation.txt", package = "KLMMAL"),header=T);
#' Data=ReadKLMM(bed=system.file("extdata", "testing.bed", package = "KLMMAL"),trainID=trainID[,1], annotation = annotation, phenofile=system.file("extdata", "pheno.txt", package = "KLMMAL"), kernels = c("linear","poly2"), AllRegions = 0);
#' Result2=KLMMALPred(Data,method="NR")
#' print(cor(Result2$Pred,Data$Y[-Data$train.index]))
#' 
#' ## Eg3: use pre-calculated GSMs
#' load(system.file("extdata", "Data.rda", package = "KLMMAL"))
#' Result3=KLMMALPred(Data,method="NR")
#' print(cor(Result3$Pred,Data$Y[-Data$train.index]))
#' 
#' @export
KLMMALPred<-function(Data,method="NR", lambdaall=NULL, stdK=TRUE,tol=1e-6,tolK=1e-6,stepmax=300)
{
  Y=Data$Y;
  K=Data$K;
  train.index=Data$train.index;
  X=NULL
  if(sum(names(Data) == "X")!=0) X=Data$X;
  reduce=FALSE
  reduce=TRUE;speedup=TRUE;epsilon=tol;minstraint=0.01;
  if(method=="NR")
  {
    Result=LMEPredictWrap(Y=Y,K=K,X=X,stdK=stdK,tol=tol, stepmax=stepmax,minstraint=minstraint,epsilon=epsilon,train.index=train.index,reduce=reduce,speedup=speedup)
  }
  if(method=="EM")
  {
    Kinship=Data$K;
    minvar=minstraint;
    FitsLinear=list();BICs=NULL;
    KinshipTrain=list(); for(i in 1:length(Kinship)) KinshipTrain[[i]]=Kinship[[i]][train.index,train.index]
    system.time(remlint<-LMEMIXEDQUICK(MyData =  list("kinship"=KinshipTrain,"pheno"=Y[train.index]), kinIn = 1,covVarRead = 0,remliter = 200))
    betaint_dint=c(remlint$beta,remlint$delta[-1]/max(1,remlint$delta[1]))
    weight=betaint_dint[-(1:length(remlint$beta))]; weight[weight==0]=minvar; weight=(1/weight);
    
    if(is.null(lambdaall))
    {
      print("lambda is not supplied: Use NR to find lambda range")
      Result=LMEPredictWrap(Y=Y,K=K,X=X,stdK=stdK,tol=tol, stepmax=stepmax,minstraint=minstraint,epsilon=epsilon,train.index=train.index,reduce=reduce,speedup=speedup)
      lambdaall=Result$MLE$fit$lambda;
    }
    for(i in 1:length(lambdaall))
    {
      lambda=lambdaall[i];
      cat("PLEStep: Obtaining MLE! Lambda=",lambda,"\n")
      FitsLinear[[i]]=EMAlgorithm(X=X,Y=Y[train.index],K=KinshipTrain, tol=tol, tolK=tolK, maxstep=stepmax,lambda=lambda,betaint_dint=betaint_dint,minvar=minstraint,memsave=TRUE)
      paras=FitsLinear[[i]]$sigmaKStep
      eigen0=eigenMLE(KinshipTrain,paras);
      beta=FitsLinear[[i]]$betaStep; beta=matrix(beta,ncol=1)
      sigma2=FitsLinear[[i]]$sigma2Step[1,1];
      if(is.null(X)) X1=matrix(1,nrow=length(train.index));
      if(!is.null(X)) {X1=cbind(1,X); X1=X1[train.index,]};
      N=length(Y[train.index])
      BICs=rbind(BICs,log(N)*sum(paras!=0)-2*logLikeP(sigma2,paras,K=KinshipTrain,beta,Y=Y[train.index],eigenV0=eigen0,X=X1,lambda=lambda,para_int=1/weight))
    }
      
    lambdaselect=(lambdaall)[which(BICs==min(BICs))]
    ResultEM=list();
    ResultEM$fit=FitsLinear;
    ResultEM$BICs=BICs;
    ResultEM$select=list();
    ResultEM$select$paras=FitsLinear[[which(BICs==min(BICs))]]$sigmaKStep;
    ResultEM$select$sigma2=FitsLinear[[which(BICs==min(BICs))]]$sigma2Step;
    ResultEM$select$beta=FitsLinear[[which(BICs==min(BICs))]]$betaStep;
    ResultEM$select$lambda=lambdaselect;
    ResultEM$select$allzero=FALSE;
    if(sum(ResultEM$select$paras!=0)==0) ResultEM$select$allzero=TRUE;
    paras=list();paras$beta=matrix(ResultEM$select$beta,ncol=1); paras$paras=ResultEM$select$paras; paras$allzero=ResultEM$select$allzero;
    ResultEM$Pred=predictLME(paras,YTrain=Y[train.index],K=Kinship,X=X,train.index=train.index)
    
    Result=list();
    Result$MLE=ResultEM$select;
    Result$Pred=ResultEM$Pred;
  }
  Result
}  
