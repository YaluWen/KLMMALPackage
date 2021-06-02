#' @title Genetic Risk Prediction Using a Multi-Kernel Linear Mixed Model with Adaptive Lasso
#'
#' @description
#' \code{ReadKLMM} prepares the data for analyses
#' @details
#' This function is used to pre-process the data used by the Multi-Kernel Linear Mixed Model with Adaptive Lasso method. 
#'
#' @param bed The name of the file containing the packed binary SNP genotype data. It should have the extension .bed; if it doesn't, then this extension will be appended
#' @param bim The file containing the SNP descriptions
#' @param fam The file containing subject (and, possibly, family) identifiers. This is basically a tab-delimited "pedfile"
#' @param na.strings Strings in .bam and .fam files to be recoded as NA
#' @param phenofile The file containing phenotype information. If provided, will ignore the affected values from the .fam file. The file should contain two columns: pheno: phenotypic values. ID: subject ID for each indidivuals. 
#' @param trainID Subject ID for the training individuals. 
#' @param annotation A data frame providing information to cut the genomes. It has at least four columns: gene, chr, start, end. The gene column must be unique. The other three columns contain information about chromosome, start and end of each regomic region. If Kinship is not provided, then annotation must be provided. 
#' @param Y A vector of phenotypes with each name being subject ID. If provided, will ignore \code{phenofile} and \code{fam}. 
#' @param X A matrix of demographic variables, should be of the same order as Y (i.e. \code{rownames(X)=names(Y)}). The intercept column is not needed. 
#' @param Kinship A list of the pre-calculated genomic similarity matrix, with each element being a genomic region. The row and column names for each element should be the subjectID. If provided, will ignore \code{bed} and \code{annotation}.
#' @param kernels A vector of strings only taking two values: linear (linear kernel) and poly2 (polynomial kernel with 2 df).
#' @param AllRegion A bool to indicate whether GSM should be calculated from the entire genome. Default=0. 
#' @return A list that contains Y, Kinship and train.index, all of which are needed for KLMMAL algorithm(i.e. the function \code{KLMMALPred}). 
#' @examples
#' ## Using the data with pre-calculated GSMs ##
#' load(system.file("extdata", "Data.rda", package = "KLMMAL"))
#' 
#' ## read both genotype and phenotype files from the disk, and calculate the GSMs using linear kernels ##
#' trainID=read.table(system.file("extdata", "TrainID.txt", package = "KLMMAL"));
#' annotation=read.table(system.file("extdata", "Annotation.txt", package = "KLMMAL"),header=T);
#' Data=ReadKLMM(bed=system.file("extdata", "testing.bed", package = "KLMMAL"),trainID=trainID[,1], annotation = annotation, phenofile=system.file("extdata", "pheno.txt", package = "KLMMAL"), kernels = c("linear"), AllRegions = 0);
#' 
#' ## read genotype from the disk and calculate the GSMs using linear and poly2 kernels##
#' pheno=read.table(system.file("extdata", "pheno.txt", package = "KLMMAL"),header=T);
#' Y=pheno$pheno;names(Y)=pheno$ID;
#' trainID=read.table(system.file("extdata", "TrainID.txt", package = "KLMMAL"));
#' annotation=read.table(system.file("extdata", "Annotation.txt", package = "KLMMAL"),header=T);
#' Data=ReadKLMM(bed=system.file("extdata", "testing.bed", package = "KLMMAL"),trainID=trainID[,1], annotation = annotation, Y = Y, kernels = c("linear", "poly2"), AllRegions = 0);
#' @export
ReadKLMM<-function(bed,bim,fam,  na.strings = c("0", "-9"), phenofile, trainID, annotation=NULL, Y=NULL, X=NULL, Kinship=NULL, kernels=c("linear","poly2"), AllRegions=0)
{
  ## import phenotypes ##
  if(is.null(Y) & missing (phenofile) & missing(bed)) stop("Error: You must provide phenotypes for the training data!")
  if(is.null(Y))
  {
    if(missing (phenofile)){
      warning(paste("Warning: Phenotype files are not provided and thus get phenotypes from the",fam, ".fam"))
      genodata=snpStats::read.plink(bed);
      Y=genodata$fam$affected;
      names(Y)=paste(genodata$fam$FID, genodata$fam$IID, sep="_")
    }
    if(!missing (phenofile)){
      pheno=read.table(phenofile,header = T) # pheno, and ID #
      if(sum(colnames(pheno) %in% c("pheno", "ID"))!=2) stop("Error: Phenotype files have contain pheno and ID columns")
      Y=pheno$pheno;
      names(Y)=pheno$ID;
    }
  }
  if(is.null(names(Y))) stop("Error: There is no IDs for phenotype data!")
  trainID=unique(trainID);
  if(sum(trainID %in% names(Y))!=length(trainID)) stop("Error: The data does not have all the training samples. Check!")
  train.index=which(names(Y) %in% trainID)
  if(sum(!is.na(Y[train.index])<5)) stop("Error: Less than 5 subjects for the training samples. The sample size is too small. ")
  
  ## Get GSMs ##
  if(missing(bed) & is.null(Kinship)) stop("Error: You must provide either the genotype files (in plink format) or the GSMs");
  if(!is.null(Kinship))
  {   
    # Kinships have been provided #
    n=nrow(Kinship[[1]]);
    if(length(Kinship)>1){
      for(i in 2:length(Kinship)) if(nrow(Kinship[[i]])!=n) stop("The GSMs for each region do not have the same dimension!")
    }
  }
  if(is.null(Kinship))
  {
    if(sum(kernels %in% c("linear","poly2"))==0) stop("Error: Currently only support linear and polynomial with degree 2 kernels, other kernels are not supported. You can supply the GSMs using your preferred kernels instead.")
    if(is.null(annotation)) stop("Error: You must provide the start and end of each genomic regions")
    if(nrow(annotation)<1 & !AllRegions) stop("Error: Must have at least one genomic regions")
    if(sum(colnames(annotation) %in% c("gene","chr","start", "end"))!=4) stop("Error: the annotation must be a dataframe with at least four columns: gene, chr, start, end (gene/region name(unique), chromosome, start and end of the genomic regions)")
    if(length(unique(annotation$gene))!=nrow(annotation)) stop("Error: gene/region names must be unique")
    
    genodata=snpStats::read.plink(bed);
    Geno=as(genodata$genotypes,"numeric")
    
    Kinship=list();IncludeRegions=NULL;index=1;
    print("Calculating Kinshps!")
    if(AllRegions) {Kinship[[1]]=rrBLUP::A.mat(Geno-1); IncludeRegions=rbind(IncludeRegions,"chr1:chr22:all:all")}
    
    for(i in 1:nrow(annotation))
    {
      if(i %% 20 ==1) cat("Calculating Kinships for the ",i, "th region... \n")
      inc=(genodata$map$chromosome==annotation$chr[i]) & (genodata$map$position <= annotation$end[i]) & (genodata$map$position >= annotation$start[i]);
      if(sum(inc)>0)
      {
        genenamestmp=paste(format(annotation$gene[i]), sum(inc), sep=":")
        IncludeRegions=rbind(IncludeRegions,genenamestmp)
        
        geno=Geno[,inc];
        nas=is.na(geno);
        # impute the missings with the average values #
        if(!is.null(dim(geno))){
          if(ncol(geno)>1){
            tmpimpute=apply(geno,2,mean,na.rm=T);
            tmpimputegs=matrix(rep(tmpimpute,nrow(geno)),nrow=nrow(geno),ncol=ncol(geno),byrow=T)
            geno[nas]=tmpimputegs[nas];
           }
          if(ncol(geno)==1){
            tmpimpute=mean(geno,na.rm=T);
            geno[nas]=tmpimpute;
            geno=matrix(geno,ncol=1)
          }
        }
        if(is.null(dim(geno))){
          tmpimpute=mean(geno,na.rm=T);
          geno[nas]=tmpimpute;
          geno=matrix(geno,ncol=1)
        }
                       
        tmp3=geno %*% t(geno)/ncol(geno); #tmp3=tmp/tmp2;
        add=TRUE;
        if(length(Kinship)>0)
        {
          for(kk in 1:length(Kinship))
          {
            if(all.equal(Kinship[[kk]],tmp3)==TRUE) {add=FALSE;break;}
          }
        }
        
        if(add)
        {
          Kinship[[index]]=tmp3;names(Kinship)[index]=genenamestmp;
          index=index+1;
        }
        if(!add)
        {
          names(Kinship)[length(Kinship)]=paste(names(Kinship)[length(Kinship)],genenamestmp,sep=";")
        }
      }
    }
    
    if(sum(kernels=="linear")>0 & sum(kernels=="poly2")>0) 
    {
      nK=length(Kinship); 
      for(j in 1:nK) 
        Kinship[[j+nK]]=Kinship[[j]] * Kinship[[j]]; 
      names(Kinship)=c(paste(names(Kinship)[1:nK],"_linear",sep=""), paste(names(Kinship)[1:nK],"_poly2",sep=""))
    }
    if(sum(kernels=="linear")>0 & sum(kernels=="poly2")==0) names(Kinship)= paste(names(Kinship)[1:length(Kinship)],"_linear",sep="")
    if(sum(kernels=="linear")==0 & sum(kernels=="poly2")>=0) 
    {
      nK=length(Kinship); 
      for(j in 1:nK) 
        Kinship[[j]]=Kinship[[j]] * Kinship[[j]]; 
      names(Kinship)= paste(names(Kinship)[1:nK],"_poly2",sep="")
    }
    for(i in 1:length(Kinship)) colnames(Kinship[[i]])=rownames(Kinship[[i]])=paste(genodata$fam$pedigree, genodata$fam$member, sep="_")
    
}  

  
  ## Compare the consistencies b/w GSMs and Ys ##
  if(length(Y)!=nrow(Kinship[[1]])) stop("Error: the number of individuals in phenotype data and the number of individuals in the genotype data does not match!")
  if(sum(!(names(Y) %in% colnames(Kinship[[1]])))!=0) stop("Error: The individuals in the phenotype files and the individuals in the genotypes files do not match")
  if(sum(!(colnames(Kinship[[1]]) %in% names(Y)))!=0) stop("Error: The individuals in the phenotype files and the individuals in the genotypes files do not match")
  if(sum(names(Y)!=colnames(Kinship[[1]]))!=0) {
    warning("The individual orders in the phenotype data and genotype data do not match, rearrange the genotype data!")
    sn=rank(names(Y));
    or=order(colnames(Kinship[[1]]));
    for(i in 1:length(Kinship)) Kinship[[i]]=(Kinship[[i]][or,or])[sn,sn]
    if(sum(names(Y)!=colnames(Kinship[[1]]))!=0) stop("Error: check IDs!")
  }
  
  ## Compare the consistencies between X and Ys ##
  if(!is.null(X))
  {
    if(rownames(X)!=names(Y)) stop("X and Y should be of the same order")
  }
  
  Data=list();
  Data$Y=Y;
  Data$train.index=train.index;
  Data$K=Kinship;
  Data$X=X;
  Data
}
