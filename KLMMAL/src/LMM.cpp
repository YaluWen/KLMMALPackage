// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// This functin is used to get the initial estimates for REML, mles. //
// [[Rcpp::export]]
Rcpp::List  LMEMIXEDQUICK(Rcpp::List MyData, bool kinIn=1, bool covVarRead=0,double epsilon=1e-6, int remliter=30) {
  //set random variable seed
  srand((unsigned)time(NULL));
  double minstraint=0.0001;
  /* Declare variables */
  Rcpp::StringVector filenameall;
  arma::mat Covariate;
  arma::vec Response;
  int num_fixed =0;
  
  arma::mat kinship;
  arma::vec kintraces;
  int num_kins=0;
  int num_samples_use;
  arma::cube KINSHIP(1,1,1);
  
  /* Setting Phenotypes */
  {
    printf("Main: Preparing the data-Obtaining Phenotypes...... \n");
    num_samples_use=Rcpp::as<arma::vec>(MyData["pheno"]).size();
    Response=Rcpp::as<arma::vec>(MyData["pheno"]);
  }
  
  std::cout<<"The sample size for this analysis is "<<num_samples_use<<std::endl;
  
  /* Getting Covariates */;
  printf("Main: Preparing the data- Obtaining covariates......\n");
  if(covVarRead)
  {
    Covariate=Rcpp::as<arma::mat>(MyData["cov"]);
    num_fixed = Covariate.n_cols;
  }
  std::cout<<"There are  "<<num_fixed<<" covariates considered in the analysis"<<std::endl;
  
  
  /* Geting Kinships */
  {
    //printf("--- Obtaining Kinships --- \n");
    if(kinIn) 
    {
      // Directly passed kinship to C //
      num_kins=Rcpp::as<Rcpp::List>(MyData["kinship"]).length();
      KINSHIP.resize(num_samples_use,num_samples_use,num_kins);
      kintraces.resize(num_kins);
      std::cout<<"There are "<<num_kins <<" kinship regions considered!" <<std::endl;
      for(int k=0;k<num_kins;k++)
      {
        KINSHIP.slice(k) = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(MyData["kinship"])[k]);
        kintraces[k]=arma::trace(KINSHIP.slice(k))/num_samples_use;
      } 
    }
    else
    {
      // Reading kinship from disk //
      filenameall=Rcpp::as<Rcpp::StringVector>(MyData["kinshipfiles"]);
      num_kins=filenameall.size();
      std::cout<<"There are "<<num_kins <<"kinship regions considered!" <<std::endl;
      kintraces.resize(num_kins);
      for(int k=0; k<num_kins; k++)
      {
        std::string filename=Rcpp::as< std::string >(filenameall(k));
        kinship.load(filename.c_str());  
        kintraces[k]=arma::trace(kinship)/num_samples_use;
      }
    }  
    
  }
  
  printf("Main: Obtaining the penalised estimates......\n");
  printf("PLEStep: Obtaining MLE! \n");
  
  /* Declear variables used for getting MLEs */
  bool mlerun=1;
  arma::vec deltas(num_kins+1,arma::fill::ones);
  arma::vec deltasnew(num_kins+1,arma::fill::ones);
  arma::vec pstepdirection(num_kins+1,arma::fill::ones); //direction of updates
  
  arma::mat Vmatrix(num_samples_use,num_samples_use,arma::fill::eye); //Vmatrix=deltas[0]*I+deltas[1]*K1+deltas[2]*K2+...
  arma::mat VmatrixInv(num_samples_use,num_samples_use);
  arma::vec beta(num_fixed+1,arma::fill::ones); // For covariates
  double loglike;
  double loglikeold;
  double loglikenull;
  double diffloglike;
  
  int count=0;
  
  arma::mat Pmatrix; //Pmatrix=solve(Vmatrix)-solve(Vmatrix) %*% X %*% solve(t(X) %*% solve(V) %*% X) %*% t(X) %*% solve(V)
  arma::mat XTinvVX; //XTinvVX=t(X) %*% solve(V) %*% X;
  
  arma::vec BI(num_kins+1); //first derivatives
  arma::vec Aactive; //first derivatives active sets
  
  arma::mat AI(num_kins+1,num_kins+1); //average information, second derivatives
  
  arma::mat tmp;
  arma::mat kinshipi(num_samples_use,num_samples_use,arma::fill::eye);
  
  arma::vec heritabiliy(num_kins+1,arma::fill::ones);
  arma::vec scales(num_kins+1,arma::fill::ones);
  arma::vec fixed(num_kins+1,arma::fill::zeros);
  
  double alphastep=1.0;
  double varnull=0.01;
  double varmin=0.01;
  int nfree=num_samples_use-num_fixed;
  
  
  int rflag=0; 
  int nlost=0;
  
  arma::vec eigvalV;
  arma::mat eigvecV;
  arma::mat eigvecinvV;
  
  arma::vec eigvalXTinvVX;
  arma::mat eigvecXTinvVX;
  arma::mat eigvecinvXTinvVX;
  
  arma::mat invXTinvVX;
  
  // Getting fixed effects //
  if(!covVarRead) Covariate=arma::mat(num_samples_use,1,arma::fill::ones);
  else{
    Covariate=arma::join_rows(arma::mat(num_samples_use,1,arma::fill::ones),Covariate);
  }
  
  /* Calculation is based on MultiBLUP*/;
  /* Using Newton Raspson To Update Parameters */
  // printf("--- Obtaining Estimates:  Initilization \n");
  // Initilization //
  {
    beta=arma::inv(Covariate.t() * Covariate) * Covariate.t() * Response;
    deltas(0)=arma::as_scalar(arma::cov(Response-Covariate*beta))*num_samples_use;
    varnull=deltas(0)/nfree*1.0;
    varmin=0.0001*varnull; //need likenull
    
    //set scales
    for(int k=0;k<num_kins;k++){scales(1+k)=kintraces(k);}
    
    //initialize heritabilities, variances, fixed and set count to zero
    for(int k=0;k<num_kins;k++){heritabiliy(1+k)=.5/(num_kins)+0.01;}
    heritabiliy(0)=1-0.5-0.01*num_kins;
    // heritabiliy(0)=max(heritabiliy(0),0.2);
    for(int k=0;k<1+num_kins;k++){deltas(k)=heritabiliy(k)*varnull/scales(k);}
    
    //calculate the loglikelihood under null (all kins=0)
    Vmatrix=Vmatrix*varnull;   
    
    arma::eig_sym(eigvalV, eigvecV, Vmatrix);
    eigvecinvV=arma::inv(eigvecV);
    VmatrixInv=eigvecV * arma::diagmat(1.0/eigvalV) * eigvecinvV;
    
    //VmatrixInv=arma::inv(Vmatrix);
    
    
    XTinvVX=Covariate.t()*VmatrixInv*Covariate;
    
    arma::eig_sym(eigvalXTinvVX, eigvecXTinvVX, XTinvVX);
    eigvecinvXTinvVX=arma::inv(eigvecXTinvVX);
    invXTinvVX=eigvecXTinvVX * arma::diagmat(1.0/eigvalXTinvVX) * eigvecinvXTinvVX;
    
    //Pmatrix=VmatrixInv-VmatrixInv*Covariate*arma::inv(XTinvVX)*Covariate.t()*VmatrixInv;
    Pmatrix=VmatrixInv-VmatrixInv*Covariate*invXTinvVX*Covariate.t()*VmatrixInv;
    
    double val;     double val1;
    val=arma::sum(log(arma::abs(eigvalV)));
    val1=arma::sum(log(arma::abs(eigvalXTinvVX)));
    //double sign;
    
    //arma::log_det(val, sign, Vmatrix);
    //arma::log_det(val1, sign, XTinvVX);
    
    loglikenull=nfree*log(2*M_PI)+val+val1+arma::as_scalar(Response.t() * Pmatrix.t() * Response);
    //loglikenull=nfree*log(2*M_PI)+log(fabs(detV))+log(fabs(arma::det(XTinvVX)))+arma::as_scalar(Response.t() * Pmatrix.t() * Response);
    
    loglikenull=-loglikenull/2.0;
    
  }
  loglike=loglikeold=loglikenull;
  
  //printf("--- Obtaining Estimates:  Estimation \n");
  
  while(mlerun)
  {
    /* compute Vmatrix and related :shortcut=0, get V=(delta0 I + delta1 K1 + delta2 K2 + ...)^-1 directly, currently does not support shortcut */
    {
      Vmatrix=Vmatrix.eye();
      Vmatrix=Vmatrix*deltas[0];
      // Adding kinship to it //
      for(int k=0;k<num_kins;k++)
      {
        //read kinships (and get traces direct)
        if(kinIn){Vmatrix=Vmatrix+deltas[k+1]*KINSHIP.slice(k);}
        else{
          std::string filename=Rcpp::as< std::string >(filenameall(k));
          kinship.load(filename.c_str());  
          Vmatrix=Vmatrix+deltas[k+1]*kinship;
        }
      } 
      
      //VmatrixInv=arma::inv(Vmatrix);
      arma::eig_sym(eigvalV, eigvecV, Vmatrix);
      eigvecinvV=arma::inv(eigvecV);
      VmatrixInv=eigvecV * arma::diagmat(1.0/eigvalV) * eigvecinvV;
    }
    
    /* Compuate Pmatrix :Getting P and XTinvVX*/
    {
      XTinvVX=Covariate.t()*VmatrixInv*Covariate;
      
      arma::eig_sym(eigvalXTinvVX, eigvecXTinvVX, XTinvVX);
      eigvecinvXTinvVX=arma::inv(eigvecXTinvVX);
      invXTinvVX=eigvecXTinvVX * arma::diagmat(1.0/eigvalXTinvVX) * eigvecinvXTinvVX;
      
      //Pmatrix=VmatrixInv-VmatrixInv*Covariate*arma::inv(XTinvVX)*Covariate.t()*VmatrixInv;
      Pmatrix=VmatrixInv-VmatrixInv*Covariate*invXTinvVX*Covariate.t()*VmatrixInv;
    }
    /* Compuate likelihood */
    {
      if(count==0){
        loglike=loglikeold=loglikenull;       diffloglike=0;
      }
      else 
      {
        loglikeold=loglike;
        double val;     double val1;
        //double sign;
        val=arma::sum(log(arma::abs(eigvalV)));
        val1=arma::sum(log(arma::abs(eigvalXTinvVX)));
        
        //arma::log_det(val, sign, Vmatrix);
        //arma::log_det(val1, sign, XTinvVX);
        
        //loglike=num_samples_use*log(2*M_PI)+log(fabs(detV))+log(fabs(arma::det(XTinvVX)))+arma::as_scalar(Response.t() * Pmatrix.t() * Response); 
        loglike=nfree*log(2*M_PI)+val+val1+arma::as_scalar(Response.t() * Pmatrix.t() * Response); 
        
        loglike=loglike/(-2.0);
        diffloglike=loglike-loglikeold;
      }
    }
    
    /* Get first derivaties (BI*) and second derivatives (AI), using average information */
    /* For those that has been set to zero, first derivatives =0, and second, along diag=1, off=0 */;
    //time_t my_time = time(NULL);
    
    //printf("Derivatives New : %s", ctime(&my_time));
    
    arma::mat PY= Pmatrix * Response;
    arma::cube YTPK(1, num_samples_use, num_kins+1);
    for(int i=0; i<num_kins+1; i++)
    {
      if(fixed[i]<3)
      {
        //std::cout<<"Enter i"<<i<<std::endl;
        if(i==0) YTPK.slice(i)=Response.t() * Pmatrix;
        else{
          if(kinIn) YTPK.slice(i)= Response.t() * Pmatrix * KINSHIP.slice(i-1);
          else{
            std::string filename=Rcpp::as< std::string >(filenameall(i-1));
            kinship.load(filename.c_str()); 
            YTPK.slice(i)=Response.t() * Pmatrix * kinship;
          }
        }
      }
    }
    
    // ctime() used to give the present time
    //time_t my_time1 = time(NULL);
    
    //printf("Derivatives New First: %s", ctime(&my_time1));

    for(int i=0; i<num_kins+1; i++)
    {
      if(fixed[i]<3)// the ith is not zero
      {
        //arma::mat YTPK=Response.t() * PK.slice(i);
        double PK=0.0;
        if(i!=0)
        {
          for(int kk=0; kk<Pmatrix.n_rows;kk++)
          {
            if(kinIn) PK=PK+arma::as_scalar(Pmatrix.row(kk)*KINSHIP.slice(i-1).col(kk));
            else{
              std::string filename=Rcpp::as< std::string >(filenameall(i-1));
              kinship.load(filename.c_str()); 
              PK=PK+arma::as_scalar(Pmatrix.row(kk)*kinship.col(kk));;
            }
          }
        }
        else PK=arma::trace(Pmatrix);
        

        BI(i)=-0.5*(PK-arma::as_scalar(YTPK.slice(i) * PY));
        //BI(i)=-0.5*(arma::trace(tmp)-arma::as_scalar(Response.t() * Pmatrix * tmp *Response));
        
        for(int j=0; j<num_kins+1; j++)
        {
          if(fixed[j]<3) // The jth is not zero 
          {
            if(kinIn){
              if(j==0) AI(i,j)=-0.5*(arma::as_scalar(YTPK.slice(i) *  Pmatrix  * PY ));
              else AI(i,j)=-0.5*(arma::as_scalar((YTPK.slice(i) *  Pmatrix) * (KINSHIP.slice(j-1) * PY )));
            }
            else
            {
              if(j==0) AI(i,j)=-0.5*(arma::as_scalar(YTPK.slice(i) *  Pmatrix  * PY ));
              else {
                std::string filename=Rcpp::as< std::string >(filenameall(j-1));
                kinship.load(filename.c_str()); 
                AI(i,j)=-0.5*(arma::as_scalar((YTPK.slice(i) *  Pmatrix) * (kinship * PY )));
              }
              
            }
          }
          else AI(i,j)=0;
        }
      }
      else //The ith is zero
      {
        BI(i)=0;
        for(int j=0; j<num_kins+1; j++) AI(i,j)=0;
        AI(i,i)=1;
      }
    }
   // time_t my_time2 = time(NULL);
    
   // printf("Derivatives New Second: %s", ctime(&my_time2));
  

    YTPK.reset();
    
    /* Adding checking here */
    if(diffloglike<-.5*fabs(loglikenull))//last move was very bad, return towards previous state
    {
      loglike=loglikeold;
      //printf("The last move was very poor - will return towards previous state\n");
      deltas=deltas-alphastep*pstepdirection; //return to the previous values
      alphastep=alphastep/2;//reduce stepsize to be half
      deltasnew=deltas+alphastep*pstepdirection; //get new delats
      if(alphastep<epsilon) {mlerun=0;rflag=3;} //Can move step stuck at the poor location //
    }  
    else
    {
      alphastep=1;
      /* compute direction pK:pstepdirection using MultiBLUP */
      {

        pstepdirection= - arma::inv( AI ) * BI;

        double max=0;
        for(int i=0;i<pstepdirection.size();i++)
        {
          
          if(fixed[i]<2)
          {
            max=std::max(max,fabs(pstepdirection[i]));
          }
        }
        
        if(max>.1*varnull)	//then reduce moves
        {
          for(int k=0;k<pstepdirection.size();k++){pstepdirection[k]*=.1*varnull/max;}
        };
        
        deltasnew=deltas+alphastep*pstepdirection;
        double sum1=0.0;for(int kk=0;kk<1+num_kins;kk++){sum1+=scales[kk]*deltasnew[kk];}
        
        for(int k=0;k<pstepdirection.size();k++)
        {
          if(fixed[k]<3)	// does not meet the constraint criteria to set to be zero
          {
            double tmp=scales[k]*deltasnew[k]/sum1    ; // heritability of the new values //
            if(tmp<minstraint){pstepdirection[k]=varmin-deltas[k];fixed[k]++;} else{fixed[k]=0;}
            if(fixed[k]==3){pstepdirection[k]=-deltas[k];} //seting the parameters to be zero
          }
          else{fixed[k]++;pstepdirection[k]=-deltas[k];}	//will already have pstepdirection=0
        }
      }
      
      
      
      /* update deltasnew */
      deltasnew=deltas+alphastep*pstepdirection;
      
      /* Calculate the heritability */
      double sum=0.0;for(int k=0;k<1+num_kins;k++){sum+=scales[k]*deltas[k];}
      for(int k=0;k<1+num_kins;k++){heritabiliy[k]=scales[k]*deltas[k]/sum;}
      
    }
    
    /* Print out results */
    /*
    if(count==0)
    {
      printf("Inter\t");
      for(int k=0;k<num_kins;k++){printf("K_%d\t",k+1);}
      printf("loglikelihood\tdiffloglike\tTarget\t\Num_constraint \n");
      printf("Start\t");
      for(int k=0;k<num_kins;k++){printf("%.4f\t", heritabiliy[1+k]);}
      nlost=0; 
      printf("%.4f\tNA\t\t%.4f\t%d \n", loglike, sqrt(epsilon), nlost );
      
    }
    else{
      printf("Inter%d\t",count);
      for(int k=0;k<num_kins;k++){printf("%.4f\t", heritabiliy[1+k]);}
      nlost=0; for(int k=0;k<num_kins+1;k++) {if(fixed[k]>2 )nlost=nlost+1;}
      printf("%.4f\t%.6f\t%.4f\t%d \n", loglike, diffloglike, sqrt(epsilon), nlost );
    }
    */
    
    
    
    
    
    
    if(count>0 & fabs(diffloglike)<epsilon) {rflag=0; mlerun=0;}//likelihood converge
    //if(count>0 & arma::as_scalar(BI.t() * BI) < epsilon) {rflag=1; mlerun=0;}//first derivaties
    
    if(count>remliter) {
      mlerun=0;  rflag=2;
    }
    
    count++;
    if(mlerun) deltas=deltasnew;
    //mlerun=0;
    
  }
  
  /* Print out the final Result */;
  printf("Final%d\t",count);
  for(int k=0;k<num_kins;k++){printf("%.4f\t", heritabiliy[1+k]);}
  nlost=0; for(int k=0;k<num_kins;k++) {if(fixed[k]>2 )nlost=nlost+1;}
  printf("%.4f\t%.6f\t%.4f\t%d \n", loglike, diffloglike, sqrt(epsilon), nlost );
  
  
  
  //std::cout<<rflag<<"--"<<alphastep<<std::endl;
  if(rflag==0)
    printf("Likelihood converges \n");
  if(rflag==1)
    printf("First derivatives converges \n");
  if(rflag==2)
    printf("Warning, REML did not converge after %d iterations; consider using \"--reml-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance", count);
  if(rflag==3)
    printf("Warning, steps are too small %.4f", alphastep);
  
  
  /* Calculating parameters for the fixed effects */
  Vmatrix.eye();
  Vmatrix=Vmatrix*deltas[0];
  // Adding kinship to it //
  for(int k=0;k<num_kins;k++)
  {
    //read kinships (and get traces direct)
    if(kinIn)kinship = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(MyData["kinship"])[k]);
    else{
      std::string filename=Rcpp::as< std::string >(filenameall(k));
      kinship.load(filename.c_str());  
    }
    Vmatrix=Vmatrix+deltas[k+1]*kinship;
  } 
  VmatrixInv=arma::inv(Vmatrix);
  beta=arma::inv(Covariate.t() * VmatrixInv * Covariate)*Covariate.t()*VmatrixInv * Response;
  
  double sum=0.0;for(int k=0;k<1+num_kins;k++){sum+=scales[k]*deltas[k];}
  for(int k=0;k<1+num_kins;k++){heritabiliy[k]=scales[k]*deltas[k]/sum; }
  
  return Rcpp::List::create(Rcpp::Named("beta")=beta,
                            Rcpp::Named("delta")=deltas,
                            Rcpp::Named("her")=heritabiliy);
  
}



