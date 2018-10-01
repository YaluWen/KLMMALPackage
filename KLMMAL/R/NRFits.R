MyGBLUP<-function(Amat,Y,trainingID, X=NULL)
{
  YFit=Y;
  YFit[-trainingID]=NA;
  ans <- rrBLUP::mixed.solve(YFit,K=Amat, X=X)
  Ypred=ans$u;
  Ypred
}

NRMLE<-function(Y,K,X=NULL,stdK=TRUE,tol=1e-6, stepmax=300,minstraint=0.01,epsilon=1e-6,speedup=TRUE)
{
  NUM0=NULL;
  N=length(Y);
  NK=length(K);
  ## standardized Kinships ##
  if(stdK)
  {
    for(i in 1:length(K))
    {
      tt1=psych::tr(K[[i]])/N
      tt=K[[i]]/tt1;
      K[[i]]=tt
    }
  }
    
  
  if(!is.null(X)) X=cbind(1,X)
  if(is.null(X)) X=matrix(1,length(Y))
  
  dsnew=dsold=rep(0.5/NK+0.01,NK);
  if(sum(dsnew)>1) dsnew=dsnew/sum(dsnew)*0.99;
  dsold=dsnew;

  betanew=betaold=solve(t(X) %*% X) %*% t(X) %*% Y
  
  sigmanull=sigma2new=sigma2old=sqrt(t(Y-X%*% betanew) %*% (Y-X %*% betanew ) /N)[1,1]
  
  loglikenull=logLike(sigma2old,dsold,K,betaold,Y,X=X)
  
  loglikeold=loglikenew=loglikenull;
  loglikediff=100
  step=1;
  parasold=parasnew=c(dsnew); heri=rep(0,length(dsnew));
  num0=rep(0,length(parasold));
  alphastep=1;rerun=FALSE;
  
  while(step < stepmax & abs(loglikediff) > tol)
  {
    if((step %% 100 ==1) & step>1) cat("Estimating MLEs Step=",step,"\n");
   # cat(step,"-:",abs(loglikediff),"-",num0,"\n");
   # cat(parasnew,"\n");
 
    if(!rerun)
    {
      loglikeold=loglikenew;
      ds=parasold=parasnew;
      betaold=betanew;
      if(step>1) {res=resnew;solveV0K=solveV0Knew;sigma2=sigma2new}
      if(step==1)
      {
        res=Y-X %*% betaold;
        if(speedup) {
          K1=list();
          paraindex0=which(parasold!=0);
          for(i in 1:length(paraindex0)) K1[[i]]=K[[paraindex0[i]]];
          solveV0K=sovleV0Kfun(K1,parasold[paraindex0]);
          rm(paraindex0,K1);
        }
        if(!speedup) solveV0K=sovleV0Kfun(K,ds);
        sigma2=sqrt(t(res) %*% solveV0K$invV0 %*% res/N);
      }
      if(!speedup) derivs<-derivsFunFisher(solveV0K,res,paras=parasold,N=N,sigma2=sigma2)
      if(speedup) derivs<-derivsFunFisher(solveV0K,res,paras=parasold[parasold!=0],N=N,sigma2=sigma2)
      firstderi=derivs$firstderi;
      secondderi=derivs$secondderi;
      stepsize=solve(secondderi) %*% firstderi;
      if(speedup){stepsize1=rep(0,length(parasold)); stepsize1[parasold!=0]=stepsize; stepsize=stepsize1;rm(stepsize1);}
    }
   
    # adding checks following MultiBLUP #
    maxstepsize=max(abs(stepsize));
    
    # if the stepsize is too large, then reduce the moves #
   # if(maxstepsize>.1)	
  #  {
  #    print("reduce Step")
  #    stepsize=stepsize*0.5/maxstepsize;
  #  };
    
    if(!rerun) alphastep=1;
    parasnewtemp=parasold + alphastep*stepsize; 
    parasnewtemp[parasnewtemp<0]=0;
    eigen0=eigenMLE(K,parasnewtemp);
    if(!speedup) solveV0Knew=sovleV0Kfun(K,parasnewtemp);
    if(speedup)
    {
      #cat("AAA 2","\n");
      K1=list();
      paraindex0=which(parasnewtemp!=0);
      if(length(paraindex0)>0)
      {
       for(i in 1:length(paraindex0)) K1[[i]]=K[[paraindex0[i]]];
       solveV0Knew=sovleV0Kfun(K1,parasnewtemp[paraindex0]);
       rm(K1);
      }
      if(length(paraindex0)<1) solveV0Knew=sovleV0Kfun(K,parasnewtemp);
      rm(paraindex0)
    }
    betanew=solve(t(X) %*% solveV0Knew$invV0 %*% X) %*% t(X)  %*% (solveV0Knew$invV0  %*% Y);
    resnew=Y-X %*% betanew;
    sigma2new=sqrt(t(resnew) %*% solveV0Knew$invV0 %*% resnew /N)
    loglikenew=logLike(sigma2=sigma2new,ds=parasnewtemp,K=K,beta=betanew,Y=Y,eigenV0=eigen0,X=X)
    loglikediff=loglikenew-loglikeold;
    #print(loglikediff) 
    #print(loglikenew)
    #print(loglikeold)
    #print("start") 

    # last move was very bad, return towards previous state
    if(is.na(loglikediff)) save.image("Error.RData");
    if(is.na(loglikenull)) save.image("Error.RData");
    if(loglikediff< (-.5*abs(loglikenull)))
    {
      #print("The last move was very bad and reduce step size!")
      alphastep=alphastep/2;
      if(alphastep<epsilon) {mlerun=0;rflag=3;} #Can move step stuck at the poor location 
      rerun=TRUE;
    }  
    if(loglikediff >= (-.5*abs(loglikenull)))
    {
      rerun=FALSE;

      tot_heri=1+sum(parasnewtemp);
      index1=num0>=3 # set zeros 
      stepsize[index1]=-parasold[index1]/alphastep; #num0[index1]=num0[index1]+1; # set zeros
      index=(parasnewtemp/tot_heri)<minstraint & (num0<4); # set to mins
      if(sum(index)>0)
      {
        stepsize[index]=(minstraint-parasold[index])/alphastep
        num0[index]=num0[index]+1;
      }
      index2=((parasnewtemp/tot_heri)>=minstraint) & (parasold!=0) & num0<4 
      if(sum(index2)>0) num0[index2]=0;
      
      NUM0=rbind(NUM0,num0)
      
      # check if stuck at one position, 4 consecutive 01010101 or 012012012012 #
      if(nrow(NUM0)>12)
      {
        tt1=NUM0[(nrow(NUM0)-8+1):nrow(NUM0),]; 
        tt2=NUM0[(nrow(NUM0)-12+1):nrow(NUM0),]; 
        tt1=apply((tt1),2,paste,collapse="")
        tt2=apply((tt2),2,paste,collapse="");
        #print(tt1)
        if(sum(tt1=="01010101")>0 | sum(tt1=="10101010")>0) {num0[tt1=="01010101"|tt1=="10101010"]=4}
        if(sum(tt2=="012012012012")>0 | sum(tt2=="120120120120")>0 | sum(tt2=="201201201201")>0) {num0[tt2=="012012012012"|tt2=="120120120120"|tt2=="201201201201"]=4;};
      }
      
      # update deltasnew #
      parasnew=parasold + alphastep*stepsize; 
      heri=parasnew/(sum(parasnew)+1)
      
      if(sum(parasnew!=parasnewtemp)!=0)
      {
        # Calculate the heritability #
        #cat("AAA 1","\n");
        eigen0=eigenMLE(K,parasnew);
        if(!speedup) solveV0Knew=sovleV0Kfun(K,parasnew);
        if(speedup)
        {
          K1=list();
          paraindex0=which(parasnew!=0);
          if(length(paraindex0)>0)
          {  
            for(i in 1:length(paraindex0)) K1[[i]]=K[[paraindex0[i]]];
            solveV0Knew=sovleV0Kfun(K1,parasnew[paraindex0]);
            rm(K1)
          }
          if(length(paraindex0)<1) solveV0Knew=sovleV0Kfun(K,parasnew);
          rm(paraindex0);
        }
        betanew=solve(t(X) %*% solveV0Knew$invV0 %*% X) %*% t(X)  %*% (solveV0Knew$invV0  %*% Y);
        resnew=Y-X %*% betanew;
        sigma2new=sqrt(t(resnew) %*% solveV0Knew$invV0 %*% resnew /N)
        loglikenew=logLike(sigma2=sigma2new,ds=parasnew,K=K,beta=betanew,Y=Y,eigenV0=eigen0, X=X)
        loglikediff=loglikenew-loglikeold;
      }
      step=step+1;
    }
  }
  result=list();
  result$conv=0; if(step<maxstepsize) result$conv=1;
  result$beta=betanew;
  result$sigma2=sigma2new;
  result$ds=parasnew;
  result$heri=heri;
  result;
}

eigenMLE<-function(K,ds)
{
  V=diag((nrow(K[[1]])));
  for(i in 1:length(ds)) V=V+ds[i]*K[[i]];
  return(eigen(V));
}

logLike<-function(sigma2,ds,K,beta,Y,eigenV0=NULL,X=X)
{
  res=Y-X %*% beta;
  if(is.null(eigenV0))
  {
    V=diag(length(Y))*sigma2;
    for(i in 1:length(ds)) V=V+ds[i]*K[[i]]*sigma2
    bb=eigen(V);
    inv=bb$vectors %*% diag(1/bb$values) %*% solve(bb$vectors)
    logdetV=sum(log(abs(bb$values)));
  }
  if(!is.null(eigenV0))
  {
    inv=eigenV0$vectors %*% diag(1/eigenV0$values/sigma2) %*% solve(eigenV0$vectors)
    logdetV=sum(log(abs(eigenV0$values*sigma2)));
  }
  -(logdetV+t(res) %*% inv %*% res)/2
}

sovleV0Kfun<-function(K,ds)
{
  V=diag((nrow(K[[1]])));
  for(i in 1:length(ds)) V=V+ds[i]*K[[i]];
  #print(ds)
  invV0=solve(V);
  invV0K=list();
  for(i in 1:length(ds)) invV0K[[i]]=invV0 %*% K[[i]];
  result=list();
  result$invV0=invV0;
  result$invV0K=invV0K;
  result
}

derivsFunFisher<-function(solveV0K,res,paras,N,sigma2)
{
  ds=paras[1:length(paras)]
  tmp=res %*% (t(res) %*% solveV0K$invV0);
  firstderi=rep(0,length(ds));
  secondderi=matrix(0,length(ds),length(ds));
  diag(secondderi)=1;
  select=which(paras!=0);
  if(length(select)>0)
  {
    for(i in 1:length(select))
    {
      tmp1=sum(diag(solveV0K$invV0K[[select[i]]]));
      firstderi[select[i]]=-1/2*tmp1+1/2/sigma2*sum(tmp * (solveV0K$invV0K[[select[i]]]));
    }
    for(i in 1:length(select))
      for(j in i:length(select))
      {
        secondderi[select[i],select[j]]=secondderi[select[j],select[i]]=1/2*sum(solveV0K$invV0K[[select[i]]] * t(solveV0K$invV0K[[select[j]]]))
      }
  }
  secondderi
  result=list();
  result$firstderi=firstderi;
  result$secondderi=secondderi;
  result
}


logLikeP<-function(sigma2,ds,K,beta,Y,eigenV0=NULL,X=X,lambda=0,para_int)
{
  res=Y-X %*% beta;
  if(is.null(eigenV0))
  {
    V=diag(length(Y))*sigma2;
    for(i in 1:length(ds)) V=V+ds[i]*K[[i]]*sigma2
    bb=eigen(V);
    inv=bb$vectors %*% diag(1/bb$values) %*% solve(bb$vectors)
    logdetV=sum(log(abs(bb$values)));
  }
  if(!is.null(eigenV0))
  {
    inv=eigenV0$vectors %*% diag(1/eigenV0$values/sigma2) %*% solve(eigenV0$vectors)
    logdetV=sum(log(abs(eigenV0$values*sigma2)));
  }
  -(logdetV+t(res) %*% inv %*% res)/2-lambda*sum(abs(ds)/para_int);
}

NRMLE_L1App<-function(Y,K,X=NULL,stdK=TRUE,tol=1e-6, stepmax=300,minstraint=0.01,epsilon=1e-6,para_int=NULL,reduce=FALSE,speedup=TRUE)
{
  
  N=length(Y);
  NK=length(K);
  ## standardized Kinships ##
  if(stdK)
  {
    for(i in 1:length(K))
    {
      tt1=psych::tr(K[[i]])/N
      tt=K[[i]]/tt1;
      K[[i]]=tt
    }
  }
  if(is.null(para_int)) {print("PLEStep: Obtaining MLE!"); result=NRMLE(Y,K,X,stdK,tol,stepmax,minstraint,epsilon,speedup=speedup);  para_int=result$ds; }
  
  if(length(K)>500)
  {
    print("PLEStep: Will use the reduced version by setting reduce=TRUE!")
    reduce=TRUE;
  }
  if(!is.null(X)) X=cbind(1,X)
  if(is.null(X)) X=matrix(1,length(Y))
  
  if(sum(result$ds!=0)<5) reduce =FALSE; 
  if(reduce)
  {
    Ksave=K;para_int_save=para_int;
    select=which(para_int!=0);
    if(length(select)==0) stop("None of the genetic are predictive! Use linear regression should be enough!")
    K=list();
    for(i in 1:length(select))
    {
      K[[i]]=Ksave[[select[i]]]
    }
    para_int=para_int[select];
  }
  eigen0=eigenMLE(K,para_int);
  solveV0K=sovleV0Kfun(K,para_int);
  beta=solve(t(X) %*% solveV0K$invV0 %*% X) %*% t(X)  %*% (solveV0K$invV0  %*% Y);
  res=Y-X %*% beta;
  sigma2=sqrt(t(res) %*% solveV0K$invV0 %*% res/N);
  derivs<-derivsFunFisher(solveV0K=solveV0K,res=res,paras=para_int,N=N,sigma2=sigma2)
  firstderi=derivs$firstderi;
  secondderi=derivs$secondderi; # Fisher at MLE #
  
  ## One step approaximation ##
  print("PLEStep: One step approaximation!"); 
  weight=para_int; weight[weight==0]=minstraint; weight=(1/weight);
  
  l=KFAS::ldl(secondderi,tol=1e-10);
  d=matrix(0,nrow=nrow(l),ncol=ncol(l))
  diag(d)=sqrt(diag(l));
  diag(l)=1;
  A=l %*% d;
  A=t(A);
  
  tt=diag(nrow(A)); diag(tt)=1/weight;
  Xstar=A %*% tt;
  Ystar=t(solve(A)) %*% (firstderi+secondderi %*% para_int);
  fit<-glmnet::glmnet(x=Xstar,y=Ystar,lower.limits = 0, standardize = F, intercept = F,lambda = seq(from=2, to=0, by=-0.1));
 
  BICs=NULL; COEFS=as.matrix(coef(fit))[-1,]; scale=1
  while(sum(COEFS[,1])!=0)
  {
    scale=scale*2;
    fit<-glmnet::glmnet(x=Xstar,y=Ystar,lower.limits = 0, standardize = F, intercept = F,lambda = scale*seq(from=2, to=0, by=-0.1));
    COEFS=as.matrix(coef(fit))[-1,];
  }
  select=which(apply(COEFS,2,sum)!=0)
  COEFS=tt %*% COEFS
  if(length(select)>0)
  {
    for(i in 1:length(select))
    {
      #print(select[i])
      paras=COEFS[,select[i]];
      eigen0=eigenMLE(K,paras);
      solveV0K=sovleV0Kfun(K,paras);
      beta=solve(t(X) %*% solveV0K$invV0 %*% X) %*% t(X)  %*% (solveV0K$invV0  %*% Y);
      res=Y-X %*% beta;
      sigma2=sqrt(t(res) %*% solveV0K$invV0 %*% res/N);
      BICs=rbind(BICs,
                 log(N)*sum(paras!=0)-2*logLikeP(sigma2,paras,K,beta,Y,eigenV0=eigen0,X=X,lambda=fit$lambda[select[i]],para_int=1/weight))
      
    }
    
    lambdaselect=(fit$lambda[select])[which(BICs==min(BICs))]
    
    ## tuning lambda ##;
    fit<-glmnet::glmnet(x=Xstar,y=Ystar,lower.limits = 0, standardize = F, intercept = F,lambda = seq(from=lambdaselect+0.1*scale, to=max(0,round(lambdaselect-0.1*scale,6)), by=-0.01*scale));
    BICs=NULL; COEFS=as.matrix(coef(fit))[-1,]; 
    select=which(apply(COEFS,2,sum)!=0)
    COEFS=tt %*% COEFS
    for(i in 1:length(select))
    {
      #print(select[i])
      paras=COEFS[,select[i]];
      eigen0=eigenMLE(K,paras);
      solveV0K=sovleV0Kfun(K,paras);
      beta=solve(t(X) %*% solveV0K$invV0 %*% X) %*% t(X)  %*% (solveV0K$invV0  %*% Y);
      res=Y-X %*% beta;
      sigma2=sqrt(t(res) %*% solveV0K$invV0 %*% res/N);
      BICs=rbind(BICs,
                 log(N)*sum(paras!=0)-2*logLikeP(sigma2,paras,K,beta,Y,eigenV0=eigen0,X=X,lambda=fit$lambda[select[i]],para_int=1/weight))
      
    }
    
    lambdaselect=(fit$lambda[select])[which(BICs==min(BICs))]
    paras=COEFS[,fit$lambda==lambdaselect];
    eigen0=eigenMLE(K,paras);
    solveV0K=sovleV0Kfun(K,paras);
    beta=solve(t(X) %*% solveV0K$invV0 %*% X) %*% t(X)  %*% (solveV0K$invV0  %*% Y);
    res=Y-X %*% beta;
    sigma2=sqrt(t(res) %*% solveV0K$invV0 %*% res/N);
    
    result=list();
    result$fit=fit;
    result$BICs=BICs;
    result$COEFs=COEFS;
    result$select=list();
    result$select$paras=paras;
    if(reduce)
    {
      #print(para_int_save)
      paratmp=rep(0,length(para_int_save));
      paratmp[para_int_save!=0]=paras;
      result$select$paras=paratmp;
      result$select$reduce=which(para_int_save!=0);
      para_int=para_int_save;
    }
    result$select$sigma2=sigma2;
    result$select$beta=beta;
    result$select$lambda=lambdaselect;
    result$para_int=para_int;
    result$select$allzero=FALSE;
   # result;
  }
  if(length(select)==0)
  {
    paras=COEFS[,1];
    eigen0=eigenMLE(K,paras);
    solveV0K=sovleV0Kfun(K,paras);
    beta=solve(t(X) %*% solveV0K$invV0 %*% X) %*% t(X)  %*% (solveV0K$invV0  %*% Y);
    res=Y-X %*% beta;
    sigma2=sqrt(t(res) %*% solveV0K$invV0 %*% res/N);
    
    result=list();
    result$fit=fit;
    result$BICs=BICs;
    result$COEFs=COEFS;
    result$select=list();
    result$select$paras=paras;
    result$select$sigma2=sigma2;
    result$select$beta=beta;
    result$select$lambda=0;
    result$para_int=para_int;
    result$select$allzero=TRUE;
  }
  result;
}


### prediction ###
predictLME<-function(paras,YTrain,K,X,train.index)
{
  if(!paras$allzero)
  {
    if(!is.null(X)) X=cbind(1,X);
    if(is.null(X)) X=matrix(1,nrow(K[[1]]))
    beta=paras$beta;
    V=diag(nrow(K[[1]]));
    for(i in 1:length(K)) V=V+paras$paras[i]*K[[i]];
    V12=V[-train.index,train.index];
    V22=V[train.index,train.index];
    Ypred=X[-train.index,] %*% beta+V12 %*% solve(V22) %*% (YTrain-X[train.index,] %*% beta)
  }
  if(paras$allzero)
  {
    ## Go back to gBLUP ##
    YgBLUP=rep(NA,nrow(K[[1]]));
    YgBLUP[train.index]=YTrain;
    Amat=K[[1]]
    if(length(K)>1) for(i in 2:length(K)) Amat=Amat+K[[i]];
    Ypred=MyGBLUP(Amat=Amat,Y=YgBLUP,trainingID=train.index, X=X)
    Ypred=Ypred[-train.index]
  }
  Ypred
}

LMEPredictWrap<-function(Y,K,X=NULL,stdK=TRUE,tol=1e-6, stepmax=300,minstraint=0.01,epsilon=1e-6,train.index=NULL,reduce=FALSE,speedup=TRUE)
{
  NK=length(K); 
  print("Main: Preparing the data ........")
  ## standardized Kinships ##
  if(stdK)
  {
    for(i in 1:length(K))
    {
      tt1=psych::tr(K[[i]])/nrow(K[[1]])
      tt=K[[i]]/tt1;
      K[[i]]=tt
    }
  }
  if(!is.null(train.index))
  {
    KTrain=list();
    for(i in 1:length(K))
    {
      KTrain[[i]]=K[[i]][train.index,train.index];
    }
    YTrain=Y[train.index];
    if(!is.null(X)) XTrain=X[train.index,]
    if(is.null(X)) XTrain=X;
  }
  if(is.null(train.index))
  {
    print("The index for training data is not provided, and thus will not do prediction!")
    KTrain=K; YTrain=Y;XTrain=X;
  }
  N=length(YTrain);
  print("Main: Obtaining the penalised estimates .......")
  
  Result1=NRMLE_L1App(Y=YTrain,K=KTrain,X=XTrain,stdK=FALSE,tol=tol, stepmax=stepmax,minstraint=minstraint,epsilon=epsilon,para_int=NULL,reduce=reduce,speedup = speedup)
  Result2=list();
  if(!is.null(train.index))
  {
    print("Main: Predicting new observations  .......")
    Result2=predictLME(paras=Result1$select,Y=YTrain,K=K,X=X,train.index=train.index)
  }
  
  print("Main: All fun is yours!")
  
  Result=list();
  Result$MLE=Result1$select;
  Result$Pred=Result2;
  Result
}
