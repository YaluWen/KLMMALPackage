
predict_data_size <- function(numeric_size, number_type = "numeric") {
  if(number_type == "integer") {
    byte_per_number = 4
  } else if(number_type == "numeric") {
    byte_per_number = 8 #[ 8 bytes por numero]
  } else {
    stop(sprintf("Unknown number_type: %s", number_type))
  }
  estimate_size_in_bytes = (numeric_size * byte_per_number)
  class(estimate_size_in_bytes) = "object_size";
  estimate_size_in_bytes
}


EMAlgorithm<-function(X=NULL,Y,K, tol=1e-6, tolK=1e-8, maxstep=1000,lambda,betaint_dint=NULL,minvar=0.01,memsave=FALSE)
{
  
  bb=predict_data_size(length(Y)^2*length(K)*5, "numeric") # 1.5 Gb
  if(bb>8*1024*1024*1024) 
  {
    message("It is likely the big matrices will use all the memory, and thus forced to use memory save function")
    memsave=TRUE;
  }
  if(length(Y)>1000 & length(K))
  step=1;
  
  if(!is.null(X)) {nbeta=1+ncol(X);X=cbind(1,X);}
  if(is.null(X)) {nbeta=1;X=matrix(1,nrow=length(Y))}
  
  KmatrixList<-KmatrixCalBigAll(K,tol=tolK);
  
  if(is.null(betaint_dint)) 
  {
    system.time(remlint<-LMEMIXEDQUICK(MyData =  list("kinship"=K,"pheno"=Y), kinIn = 1,covVarRead = 0,remliter = 100))
    betaint_dint=c(remlint$beta,remlint$delta[-1]/max(1,remlint$delta[1]))
  }
  
  
  betaint_dint[betaint_dint==0]=minvar;
  paras=betaint_dint;
  paras=round(paras,-log10(tol));
  betaStep=paras[1:nbeta]
  sigmaKStep=paras[(nbeta+1):length(paras)]
  dint=1/sigmaKStep
  

  valuesold=0; valuesnew=100; step=1;
  
  while(abs(1-valuesold/valuesnew)>tol & step< maxstep)
  {
    #print(sigmaKStep);
    #print(step);
    #print(abs(1-valuesold/valuesnew),digits=15)
    valuesold=valuesnew;
    
    invV=diag(length(Y))
    for(i in 1:length(K)) invV=invV+sigmaKStep[i]^2*K[[i]];
    invV=solve(invV);
    res=Y-X %*% matrix(betaStep,ncol=1);
    sigma2Step=sqrt(t(res) %*% invV %*% res/length(Y));

    ### Create A matrix, it is of the dimension n*M ###
    A=matrix(0,nrow=length(Y),ncol=length(K));
    D1=matrix(0,nrow=length(K),ncol=length(K));
    D2=matrix(0,nrow=length(K),ncol=length(K));
    
    invVres= invV %*% res
    select=which(sigmaKStep!=0);notselect=which(sigmaKStep==0);
    if(memsave)
    {
      if(length(select)>0)
      {
        for(i in 1:length(select))
        {
          if(sigmaKStep[select[i]]!=0)
          {
            A[,select[i]]=K[[select[i]]] %*% invVres * sigmaKStep[select[i]]
          }
        }
        for(i in 1:length(select))
          for(j in i:length(select))
          {
            #cat(i,"-",j)
            D2[select[j],select[i]]=D2[select[i],select[j]]=sigmaKStep[select[i]] * sigmaKStep[select[j]] * (t(invVres) %*% K[[select[i]]] )%*% (K[[select[j]]] %*% invVres)
            
            if(i==j)
            {
              CovVars=(diag(length(Y))-(sigmaKStep[select[i]]*t(KmatrixList$Kmatrix[[select[i]]])) %*%invV %*% (KmatrixList$Kmatrix[[select[i]]] * sigmaKStep[select[i]]))*sigma2Step[1,1]
            }
            if(i!=j)
            {
              #cov1=cbind(KmatrixList$Kmatrix[[select[i]]]* sigmaKStep[select[i]],KmatrixList$Kmatrix[[select[j]]]* sigmaKStep[select[j]])  
              #CovVars=(diag(2*length(Y))-t(cov1) %*% invV %*% cov1)*sigma2Step[1,1]
              #CovVars=CovVars[1:length(Y),(length(Y)+1):(2*length(Y))]
              CovVars=-t(KmatrixList$Kmatrix[[select[i]]]* sigmaKStep[select[i]]) %*% invV %*% KmatrixList$Kmatrix[[select[j]]]* sigmaKStep[select[j]] *sigma2Step[1,1]
            }
            Bij=t(KmatrixList$Kmatrix[[select[i]]]) %*% KmatrixList$Kmatrix[[select[j]]]
            D1[select[j],select[i]]=D1[select[i],select[j]]=sum(matrixcalc::hadamard.prod(Bij,CovVars));
          }
      }
      
      if(length(notselect)>0)
      {
        for(i in 1:length(notselect))
        {
          D1[notselect[i],notselect[i]]=sum(matrixcalc::hadamard.prod(KmatrixList$Kmatrix[[notselect[i]]],KmatrixList$Kmatrix[[notselect[i]]]))*sigma2Step[1,1]
        }
      }
      
    }
   
    if(!memsave)
    {
      ## Calculate something that will be used again and again ##
      if(step==1)
      {
        KmatrixAllSigma=KmatrixAll=NULL;
        for(j in 1:length(KmatrixList$Kmatrix))
        {
          KmatrixAll=cbind(KmatrixAll,KmatrixList$Kmatrix[[j]])  
        }
      }
        
      KmatrixAllSigma=t(t(KmatrixAll)*rep(sigmaKStep,each=length(Y)))

      
      if(length(select)>0)
      {
        
        for(i in 1:length(select))
        {
          if(sigmaKStep[select[i]]!=0)
          {
            A[,select[i]]=K[[select[i]]] %*% invVres * sigmaKStep[select[i]]
          }
        }
        selectStart=(select-1)*length(Y)+1;        selectEnd=(select)*length(Y);              selectall=NULL;
        for(i in 1:length(selectStart))selectall=c(selectall,selectStart[i]:selectEnd[i])
        
        for(i in 1:length(select))
          for(j in i:length(select))
          {
            D2[select[j],select[i]]=D2[select[i],select[j]]=sigmaKStep[select[i]] * sigmaKStep[select[j]] * (t(invVres) %*% K[[select[i]]] )%*% (K[[select[j]]] %*% invVres)
          }
        for(i in 1:length(select))
        {
         # print(i);
          selectcurrent=selectall[((i-1)*length(Y)+1):length(selectall)]
          KmatrixAllSigmaSelect=KmatrixAllSigma[,selectcurrent];
          KmatrixAllSelect=KmatrixAll[,selectcurrent];
          
          CovVars=-t(KmatrixList$Kmatrix[[select[i]]]* sigmaKStep[select[i]]) %*% invV %*% KmatrixAllSigmaSelect *sigma2Step[1,1]
          CovVars[,1:length(Y)]=diag(length(Y))*sigma2Step[1,1]+CovVars[,1:length(Y)]
          Bij=t(KmatrixList$Kmatrix[[select[i]]]) %*% KmatrixAllSelect
          D1[select[i],select[i:length(select)]]=apply(matrix(apply(matrixcalc::hadamard.prod(Bij,CovVars),2,sum),nrow=length(Y)),2,sum)
        } 
      }
      
      if(length(notselect)>0)
      {
        for(i in 1:length(notselect))
        {
          D1[notselect[i],notselect[i]]=sum(matrixcalc::hadamard.prod(KmatrixList$Kmatrix[[notselect[i]]],KmatrixList$Kmatrix[[notselect[i]]]))*sigma2Step[1,1]
        }
      }
      D1[lower.tri(D1)] <- t(D1)[lower.tri(D1)]
    }
    D=D1+D2;
    
    tX_X=t(X) %*% X;tX_A=t(X) %*% A;
    Dmat1=cbind(tX_X,tX_A);
    Dmat2=cbind(t(tX_A),D);
    Dmat=rbind(Dmat1,Dmat2)
    
   # print("Dmat");
  #  print(Dmat); 
    dvec1=t(Y) %*% cbind(X,A);
    dvec2=-lambda * c(rep(0,length(betaStep)),1/dint)
    dvec=dvec1+dvec2
    
    Amat1=diag(length(sigmaKStep))
    Amat2=matrix(0,nrow=nrow(Amat1),ncol=length(betaStep))
    Amat=cbind(Amat2,Amat1);
#    print(Amat);
    bvec=rep(0,length(sigmaKStep))
    solution=quadprog::solve.QP(Dmat,dvec,t(Amat),bvec=bvec)
    
    paras=round(solution$solution,-log10(tol));
    betaStep=paras[1:nbeta]
    sigmaKStep=paras[(nbeta+1):length(paras)]
    
    valuesnew=solution$value
    step=step+1;
    
  }
  
  invV=diag(length(Y))
  for(i in 1:length(K)) invV=invV+sigmaKStep[i]^2*K[[i]];
  invV=solve(invV);
  res=Y-X %*% matrix(betaStep,ncol=1);
  sigma2Step=sqrt(t(res) %*% invV %*% res/length(Y));
  
  
  result=list(); result$converge=0;
  if(abs(1-valuesold/valuesnew)<=tol) result$converge=1;
  result$sigmaKStep=sigmaKStep
  result$betaStep=betaStep
  result$sigma2Step=sigma2Step
  result;
 ## checking conditional variance ##
  #bb=(diag(length(Y))-(sigmaKStep[i]*t(KmatrixList$Kmatrix[[i]])) %*%invV %*% (KmatrixList$Kmatrix[[i]] * sigmaKStep[i]))*sigma2Step[1,1]
  #gg=t(KmatrixList$Kmatrix[[i]]) %*% invV %*% res * sigmaKStep[i]
  ## covariance ##
  #cov1=cbind(KmatrixList$Kmatrix[[1]]* sigmaKStep[1],KmatrixList$Kmatrix[[2]]* sigmaKStep[2])
  #bb=(diag(1000)-t(cov1) %*% invV %*% cov1)*sigma2Step[1,1]
}

KmatrixCalBigAll<-function(K, tol = 1e-16,memsave=TRUE)
{
  # K is the proposed variance covariance structure, it is a list #
  KmatrixList=list();KmatrixList$zeroindex=NULL;KmatrixList$Kmatrix=list();
  for(i in 1:length(K))
  {
    l=round(KFAS::ldl(K[[i]], tol = tol),-log10(tol))
    d <- diag(sqrt(diag(l)))
    diag(l) <- 1
    tmp=l %*% d
    tmp=round(tmp,-log10(tol))
    KmatrixList$Kmatrix[[i]]=tmp
  }
  
  KmatrixList
}

