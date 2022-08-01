#' SNP marker cleaning and QC
#'
#' @description Uses a SNP marker matrix and sequentially filters for missing data, imputes remaining missing NA SNP alleles per marker, filters for minor allele frequency and 'prunes' markers based on marker correlation.
#'
#' @param marker.mat NxM marker matrix with row names and column names. Markers coded as 0,1,2 for hom, het,hom.
#' @param impute Logical whether to run imputation step. Default is TRUE
#' @param cutoff correlation threshold for pruning. Numeric between 0 and 1. Default is 0.8
#' @param min.maf Minimum minor allele frequency to keep. Numeric between 0 and 0.5. Default is 0.05
#' @param max.NAs Maximum proportion of NA missing data Numeric between 0 and 1. Default is 0.1
#' @param verbose Logical whether to print progress. Default is TRUE.
#' @param n.cores Number of cores that can be used to run imputation in parallel. The Default detectcores() from Parallel will use available multiple cores on a Windows OS but n.cores should be defined when running on a HPC or other OS.
#' @return  A NxM marker matrix.
#'
#' @examples
#' library(MTensGS)
#' data(TGdata)
#'
#' QCd.markers<-marker.QC(marker.mat = TGgenos)
#'
#' #Check size of new marker dataset
#' dim(TGgenos)
#' dim(QCd.markers)
#'
#'
#' @export
marker.QC<-function(marker.mat,#NxM marker matrix with row names and column names. Markers coded as 0,1,2 for hom, het,hom.
                    impute=T, #Logical whether to run imputation step. Default is TRUE
                    cutoff=0.8, #correlation threshold for pruning. Numeric between 0 and 1. Default is 0.8
                    min.maf=0.05, #Minimum minor allele frequency to keep. Numeric between 0 and 0.5. Default is 0.05
                    max.NAs=0.1, #Maximum proportion of NA missing data Numeric between 0 and 1. Default is 0.1
                    n.cores=NULL, #
                    verbose=T) #Logical whether to print progress. Default is TRUE.
{
  if(is.null(n.cores)){
    n.cores<-parallel::detectCores()
  }

  if(verbose==T){if(min.maf>0.2){print("High MAF!")}}
  if(verbose==T){if(max.NAs>0.3){print("High max NAs!")}}
  if(verbose==T){print("Removing high NA markers")}
  marker.mat[marker.mat==1]<-NA #Hets to NAs
  missing<-(apply(marker.mat,2,FUN = function(x) sum(is.na(x))))/nrow(marker.mat) #calculate NA freq
  marker.mat<-marker.mat[,!missing>max.NAs] #remove markers with >10% NAs
  is.mono<-apply(marker.mat,2,FUN = function(x) sum(na.omit(x=="0"))==0 | sum(na.omit(x=="2"))==0) #remove monomorphics
  marker.mat<-marker.mat[,!is.mono]

  if(impute){#impute missing markers
    require(foreach)
    if(verbose==T){print("Imputing missing data")}
    rest.of.data<-marker.mat
    median.impute<-function(x){
      x[is.na(x)]<-median(na.omit(x))
      return(x)
    }
    rest.of.data<-cbind(apply(rest.of.data,2,FUN = median.impute))
    file.remove("Imputing log.txt")
    if(verbose==T){print("Imputing log is output in directory")}
    cl <- parallel::makeCluster(n.cores,outfile="Imputing log.txt")
    doParallel::registerDoParallel(cl)
    imp.start<-Sys.time()
    imputed<-foreach::foreach(i=1:ncol(marker.mat),.combine = cbind,.multicombine = TRUE,.packages='randomForest',.verbose = F) %dopar% {
      c.sub<-apply(marker.mat,2,function(x) cor(x,marker.mat[,i],use = "pairwise.complete.obs"))
      rest.of.data.sub<-rest.of.data[,names(abs(c.sub)[order(abs(c.sub),decreasing=T)[1:100]])]
      rest.of.data.sub<-rest.of.data.sub[,!colnames(rest.of.data.sub)==colnames(marker.mat)[i]]
      isna<-is.na(marker.mat[,i])
      imp.marker<-marker.mat[,i]
      if(sum(isna)>0){
        rfmod<-randomForest::randomForest(x = rest.of.data.sub[!isna,],y = as.factor(marker.mat[!isna,i]),ntree = 20,importance = F)
        imp.marker[isna]<-as.numeric(as.character(predict(rfmod,newdata = rest.of.data.sub[isna,])))
      }
      if(isTRUE(!round((1:ncol(marker.mat)/ncol(marker.mat))*100)[i]==round((1:ncol(marker.mat)/ncol(marker.mat))*100)[i-1])){
        if(verbose==T){print(paste("|",round((1:ncol(marker.mat)/ncol(marker.mat))*100)[i],"%",sep=""))
          elapsed<-Sys.time()-imp.start
          percent.done<-((i/ncol(marker.mat))*100)
          print("Estimated time left:")
          print((100-percent.done)*(elapsed/percent.done))
        }
      }
      return(imp.marker)
    }
    parallel::stopCluster(cl)
    doParallel::stopImplicitCluster()
  }else{
    imputed<-marker.mat
  }
  dimnames(imputed)<-dimnames(marker.mat)

  if(verbose==T){print("Removing low MAF markers")}
  maf<-apply(imputed,2,function(x) min(sum(na.omit(x==2)),sum(na.omit(x==0)))/nrow(imputed)) #Calculate MAF
  all.markers.mafrm<-imputed[,!maf<min.maf] #Remove MAF<0.05 markers
  if(cutoff<1){
    if(verbose==T){print("Pruning:")}
    to.remove<-c()
    for(i in 1:ncol(all.markers.mafrm)){
      if(!colnames(all.markers.mafrm)[i]%in%to.remove){
        pruned<-all.markers.mafrm[,!colnames(all.markers.mafrm)%in%to.remove]
        tocheck<-colnames(all.markers.mafrm)[(i+1):ncol(all.markers.mafrm)]
        tocheck<-tocheck[tocheck%in%colnames(pruned)]
        if (length(tocheck)>1) {
          c.sub<-apply(pruned[,tocheck],2,function(x) cor(x,pruned[,colnames(all.markers.mafrm)[i]],use = "pairwise.complete.obs"))
        }
        linked<-names(c.sub)[abs(c.sub)>cutoff]
        to.remove<-c(to.remove,linked[!linked==colnames(all.markers.mafrm)[i]])
      }
      if(isTRUE(!round((1:ncol(all.markers.mafrm)/ncol(all.markers.mafrm))*100)[i]==round((1:ncol(all.markers.mafrm)/ncol(all.markers.mafrm))*100)[i-1])){
        if(verbose==T){cat("|",sep="")}
      }
    }
    if(verbose==T){cat(":)",sep="")}
    out<-pruned
  }else{
    out<-all.markers.mafrm
    rownames(out)<-rownames(marker.mat)
  }
  return(out)
}




#' Cross validate multi-trait ensembles
#'
#' @description Performs cross validation for several base genomic prediction models as well as multiple genetic model and multi-trait ensemble prediction models. Minimum input requirements are SNP genotype data and multi-trait phenotype data base prediction models should run well on
#' default parameters.
#'
#' @param Genos NxM marker \emph{matrix} with row names and column names as output from marker.QC(). Markers coded as 0,2 for each SNP allele. Inbred lines
#' are assumed. Must have no missing (NA) values. Row names are Genotype ID names and column are SNP marker names.
#' @param Phenos Data frame of N x p trait data. Row names are Genotype ID matching row names of Genos. Column names are trait names. Missing data is allowed.
#' @param base.models \emph{vector} of base model names to use. Options include:
#' c("rrblup","EGBLUP","RKHS","LASSO","RF","GBM","XGB")
#'
#'-\strong{rrblup} Ridge regression/GBLUP implemented in the rrblup package using a linear GRM.
#'
#'-\strong{EGBLUP} The Extended GBLUP incorporating a linear and epistatic GRM kernel implemented in the BGLR package.
#'
#'-\strong{RKHS} Replicating Kernel Hilbert Space fitting multiple Gaussian kernels using the kernel averaging method for several bandwidth
#'     parameters. Implemented in the BGLR package.
#'
#'-\strong{LASSO} Linear Shrinkage model using the LASSO penalty to shrink most SNP effects to 0. Implemented in the glmnet package. Each model is
#'     cross validated with several values of lambda.
#'
#'-\strong{RF} Random forest decision tree ensemble learning method implemented in the randomForest package.  Forests of 300 trees are grown as defualt but
#'     can be increased using rf.ntrees if needed.
#'
#'-\strong{GBM} Gradient Boosting Machine learning models implemented in the gbm package. May take a long time... 5 fold cross validation is performed on each model to
#'     determine optimum tree iteration to use. cross validations can be run in parallel on multiple cores by setting n.cores if not defined by default.
#'

#' @param n.test.folds \emph{integer} number of test folds to run. Default if 10.
#' @param n.valid.folds \emph{integer} number of validation folds to run within each test fold. Default if 8.
#' @param n.CV.rounds \emph{integer} number of rounds of k-fold cross validation for test folds to perform. Default is 2
#' @param n.valid.rounds \emph{integer} number of rounds of k-fold cross validation for internal validation folds to perform. Default is 3
#' @param CV.groups \emph{list} of vectors containing Genotype ID names. Cross validation can also be performed across predefined subsets of genotypes.
#' @param n.cores Number of cores that can be used to run cross validation of RF, GBM models. The Default detectcores() from Parallel will use available multiple
#' cores on a Windows OS but n.cores should be defined when running on a HPC or other OS.
#' @param rf.ntrees \emph{integer} Number of trees to use for each random forest model. Default is 400.
#' @param gbm.ntrees \emph{integer} Max number of trees to use for each GBM model. Default is 300.
#' @param verbose Level of progress to print throughout. 0 = none; 1 = some; 2 = all.
#' @param bglr.saveAt file path location to save temporary data files for 'BGLR' model functions. Using a local folder may speed up computing time. Default is the current working directory.
#'
#' @return Returns a list object of:
#'
#'-\strong{All out-of-fold predictions} for all genotypes for base prediction models, genetic model ensemble and multi-trait ensemble models for each trait.
#'
#'-\strong{Phenotype data} all observed multi-trait input data.
#'
#'-\strong{Prediction accuracies} list of correlation coefficients (r) for each model for each trait.
#'
#'-\strong{Base models} \emph{vector} of base prediction models used.
#'
#'-\strong{Traits} \emph{vector} of trait names used.
#'
#'-\strong{CV fold allocations} A two level list of length = n.CV.rounds, each containing the genotype ID names allocated to each test fold.
#'
#' @details Ensemble predictions for multiple genetic models and multiple traits are trained on out-of-validation fold predictions using several model approaches.
#' Including a range of genetic models that capture a diversity of genetic architectures may improve accuracy of ensembles. Unbalanced phenotyping across
#' genotypes may increase accuracy of multi-trait ensembles when several related but partially phenotyped traits are included.
#'
#' @examples
#' library(MTensGS)
#' data(TGdata)
#'
#' QCd.markers<-marker.QC(marker.mat = TGgenos)
#'
#' MTens.CV.results<-CV.MTens(Genos = QCd.markers
#'                            ,Phenos = TGphenos
#'                            ,base.models = c("rrblup","LASSO")
#'                            ,n.CV.rounds = 1,n.valid.rounds=1)
#'
#' @export
#'
CV.MTens<-function(Genos,Phenos,
                   base.models=c("rrblup","EGBLUP","RKHS","LASSO","RF","GBM"),
                   n.test.folds=10,
                   n.valid.folds=8,
                   n.CV.rounds=2,
                   n.valid.rounds=3,
                   CV.groups=NULL,
                   n.cores=NULL,
                   rf.ntrees=400,
                   gbm.ntrees=300,
                   bglr.saveAt=NULL,
                   verbose=2){
  library(foreach)
  run.start.time<-Sys.time()
  if(verbose>0){
    if(length(rownames(Genos)[!rownames(Genos)%in%rownames(Phenos)])>0){
      print("Lines in Genos but not in Phenos:")
      print(rownames(Genos)[!rownames(Genos)%in%rownames(Phenos)])}
    if(length(rownames(Phenos)[!rownames(Phenos)%in%rownames(Genos)])>0){
      print("Lines in Phenos but not in Genos:")
      print(rownames(Phenos)[!rownames(Phenos)%in%rownames(Genos)])}
    invalid.base.models<-base.models[!base.models %in% c("rrblup","EGBLUP","RKHS","LASSO","RF","GBM")]
    if (length(invalid.base.models)>0) {
      print("Invalid base model names!  :")
      print(invalid.base.models)
      print(c("Options include:",c("rrblup","EGBLUP","RKHS","LASSO","RF","GBM")))
    }
  }

  if(is.null(n.cores)){
    n.cores<-parallel::detectCores()
  }

  if(is.null(bglr.saveAt)){
    bglr.saveAt<-tempdir()
  }
  Genos<-as.matrix(Genos[rownames(Genos)%in%rownames(Phenos),])
  Phenos<-as.data.frame(Phenos[rownames(Phenos)%in%rownames(Genos),])

  if(!is.null(CV.groups)){
    n.CV.rounds<-1
    n.test.folds<-length(CV.groups)
  }

  all.CV.rounds<-list()
  traits<-colnames(Phenos)
  for(r in 1:n.CV.rounds){
    if(verbose>0){print(paste("Starting test CV round",r))}

    if(is.null(CV.groups)){
      test.folds<-split(sample(1:nrow(Phenos)),1:n.test.folds)#make folds
    }else{
      test.folds<-lapply(CV.groups,function(x) which(rownames(Phenos)%in%x) )
    }

    all.out.of.fold.preds<-list()
    for (i in 1:n.test.folds) {

      all.trait.out.of.test.fold.preds<-list()
      all.trait.out.of.valid.fold.preds<-list()

      valid.folds<-list()
      for (vr in 1:n.valid.rounds) {
        valid.folds[[vr]]<-split(sample(c(1:nrow(Phenos))[!1:nrow(Phenos)%in%test.folds[[i]]]),1:n.valid.folds)#make folds
      }
      traits<-colnames(Phenos)
      for (t in 1:length(traits)) {
        if(verbose>0){print(paste("Starting cross validation training for",traits[t]))}
        #CV each ST model for ensemble base.models----
        all.STmod.out.of.valid.fold.preds<-array(NA,dim = c(nrow(Phenos[-test.folds[[i]],]),length(base.models),n.valid.rounds),
                                                 dimnames = list(rownames(Phenos[-test.folds[[i]],]),base.models,1:n.valid.rounds))
        valid.start<-Sys.time()
        validation.grid<-expand.grid(1:n.valid.folds,1:n.valid.rounds)
        file.remove("Validation folds log.txt")
        cl <- parallel::makeCluster(n.cores,outfile="Validation folds log.txt")
        doParallel::registerDoParallel(cl)
        validation.out<-foreach::foreach(v=1:nrow(validation.grid),.combine = list,.errorhandling="stop",.verbose = F,.multicombine = T) %dopar% {
          valid.names<-rownames(Phenos)[valid.folds[[validation.grid[v,2]]][[validation.grid[v,1]]]]
          valid.fold.preds<-matrix(data=NA,nrow=length(valid.names),ncol=length(base.models),dimnames = list(valid.names,base.models))
          if(verbose>0){print(paste("========Starting",traits[t],"validation fold",validation.grid[v,1],"validation round",validation.grid[v,2]," test fold",i," CV round",r,"=========="))}
          #RRBLUP----
          if("rrblup"%in%base.models){
            start<-Sys.time()
            Phenos.na<-Phenos
            Phenos.na[test.folds[[i]],]<-NA
            Phenos.na[valid.names,]<-NA
            Phenos.na<-cbind.data.frame(rownames(Phenos.na),Phenos.na)
            colnames(Phenos.na)[1]<-"line"
            A<-rrBLUP::A.mat(Genos)
            A<-A[rownames(Phenos)[-test.folds[[i]]],rownames(Phenos)[-test.folds[[i]]]]
            gblup.model<-rrBLUP::kin.blup(data=Phenos.na,pheno = traits[t],geno = "line",K = A)
            valid.fold.preds[valid.names,"rrblup"]<-gblup.model$pred[rownames(valid.fold.preds)][valid.names]
            if(verbose>1){print("running rrBLUP")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #EGBLUP----
          if("EGBLUP"%in%base.models){
            start<-Sys.time()
            Phenos.na<-Phenos[-test.folds[[i]],]
            Phenos.na<-Phenos.na[!is.na(Phenos.na[,t])|rownames(Phenos.na)%in%valid.names,]#Remove lines with na but not in validation fold
            Phenos.na[valid.names,]<-NA
            A <- rrBLUP::A.mat(Genos[rownames(Phenos.na),])
            H<-matrixcalc::hadamard.prod(A,A)
            Phenos.na<-Phenos.na[rownames(A),]
            ETA<-list(list(K=A,model="RKHS"),list(K=H,model="RKHS"))
            fm<-BGLR::BGLR(y =Phenos.na[,t],ETA = ETA,nIter=3000,burnIn=500,df0=5,S0=2,verbose = F,saveAt = bglr.saveAt)
            fm.preds<-fm$yHat
            names(fm.preds)<-rownames(A)
            valid.fold.preds[valid.names,"EGBLUP"]<-fm.preds[rownames(valid.fold.preds)][valid.names]
            if(verbose>1){print("running EGBLUP")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #RF----
          if("RF"%in%base.models){
            start<-Sys.time()
            X<-as.matrix(Genos[rownames(Phenos)[-test.folds[[i]]],])[!is.na(Phenos[-test.folds[[i]],t]),]
            X<-X[!rownames(X)%in%valid.names,]
            rf.fit<-ranger::ranger(x=X,y=Phenos[rownames(X),t],num.tree=rf.ntrees,max.depth=NULL
                                   ,min.node.size=15,mtry = round(ncol(X)/3),num.threads = 1)
            valid.fold.preds[valid.names,"RF"]<-predict(rf.fit,data = as.matrix(Genos[valid.names,]))$predictions
            if(verbose>1){print("running RF")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #GBM----
          if("GBM"%in%base.models){
            start<-Sys.time()
            X<-as.matrix(Genos[rownames(Phenos)[-test.folds[[i]]],])[!is.na(Phenos[-test.folds[[i]],t]),]
            X<-X[!rownames(X)%in%valid.names,]
            X<-X[sample(1:nrow(X)),]
            Y<-Phenos[rownames(X),t]
            gbmdata<-cbind.data.frame(Y,X)

            gbm.tune<-expand.grid(depth=c(2:8),
                                  bag.fraction=c(0.4,0.5,0.6,0.7,0.8,0.9),
                                  minobs=round(seq(2,nrow(X)*0.05,length.out=5)))
            gbm.tune<-gbm.tune[sample(1:nrow(gbm.tune)),]
            gbm.tune.error<-c()
            if(verbose>1){print("GBM param tune")}
            for (h in 1:10){
              if(verbose>1){cat("|",sep="")}
              gbm.fit<-gbm::gbm(Y~.,data = gbmdata,distribution = "gaussian",train.fraction = 0.8
                                ,interaction.depth = gbm.tune$depth[h],bag.fraction =gbm.tune$bag.fraction[h],
                                n.minobsinnode = gbm.tune$minobs[h],n.trees = 50,shrinkage = 0.2,verbose = F,n.cores=1)
              gbm.tune.error[h]<-min(gbm.fit$valid.error)
            }
            for (h in 1:5){
              if(verbose>1){cat("|",sep="")}
              tune.model<-randomForest::randomForest(x=gbm.tune[which(!is.na(gbm.tune.error)),],y=gbm.tune.error[!is.na(gbm.tune.error)],ntree=500,importance=F)
              full.error.preds<-predict(tune.model,newdata = gbm.tune)
              full.error.preds[which(!is.na(gbm.tune.error))]<-NA
              next.index<-which.min(full.error.preds)
              gbm.fit<-gbm::gbm(Y~.,data = gbmdata,distribution = "gaussian",train.fraction = 0.8
                                ,interaction.depth = gbm.tune$depth[next.index],bag.fraction =gbm.tune$bag.fraction[next.index],
                                n.minobsinnode = gbm.tune$minobs[next.index],n.trees = 50,shrinkage = 0.2,verbose = F,n.cores=1)
              gbm.tune.error[next.index]<-min(gbm.fit$valid.error)
            }
            opt.params<-gbm.tune[which.min(gbm.tune.error),]

            gbm.fit<-gbm::gbm(Y~.,data = gbmdata,cv.folds = 5,distribution = "gaussian",interaction.depth = opt.params$depth,n.minobsinnode = opt.params$minobs,
                              n.trees = gbm.ntrees,shrinkage = 0.03,verbose = F,n.cores=1)
            opt.tree<-gbm::gbm.perf(gbm.fit,plot.it = F,method = "cv")
            if(verbose>1){print(paste("GBM opt run tree no. =",opt.tree))}
            valid.fold.preds[valid.names,"GBM"]<-gbm::predict.gbm(gbm.fit,n.trees = opt.tree,newdata = as.data.frame(Genos[valid.names,]))
            if(verbose>1){print("running GBM")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #LASSO----
          if("LASSO"%in%base.models){
            start<-Sys.time()
            X<-as.matrix(Genos[rownames(Phenos)[-test.folds[[i]]],])[!is.na(Phenos[-test.folds[[i]],t]),]
            X<-X[!rownames(X)%in%valid.names,]
            cv.fit<-glmnet::cv.glmnet(x=X,y=Phenos[rownames(X),t],nfolds = 8,lambda = seq(0.000001,0.9,length.out=50)^4)
            lasso.fit<-glmnet::glmnet(x=X,y=Phenos[rownames(X),t],lambda = cv.fit$lambda.min,alpha = 1)
            valid.fold.preds[valid.names,"LASSO"]<-predict(lasso.fit,newx = as.matrix(Genos[valid.names,]))
            if(verbose>1){print("running LASSO")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #RKHS----
          if("RKHS"%in%base.models){
            start<-Sys.time()
            Phenos.na<-Phenos[-test.folds[[i]],]
            Phenos.na<-Phenos.na[!is.na(Phenos.na[,t])|rownames(Phenos.na)%in%valid.names,]#Remove lines with na but not in validation fold
            Phenos.na[valid.names,]<-NA
            Genos.sub<-Genos[rownames(Phenos)[-test.folds[[i]]],]
            D<-as.matrix(dist(Genos.sub,method="euclidean"))^2 #Compute Gaussian kernel
            D<-D/mean(D)
            h<-0.5*c(1/5,1,5)
            ETA<-list(list(K=exp(-h[1]*D),model='RKHS'),
                      list(K=exp(-h[2]*D),model='RKHS'),
                      list(K=exp(-h[3]*D),model='RKHS'))

            Phenos.na<-Phenos.na[rownames(D),]
            fm<-BGLR::BGLR(y =Phenos.na[,t],ETA = ETA,nIter=3000,burnIn=500,df0=5,S0=2,verbose = F,saveAt = bglr.saveAt)
            fm.preds<-fm$yHat
            names(fm.preds)<-rownames(D)
            valid.fold.preds[valid.names,"RKHS"]<-fm.preds[rownames(valid.fold.preds)]
            if(verbose>1){print("running RKHS")}
            if(verbose>1){print(Sys.time()-start)}
          }
          return(valid.fold.preds)
        }
        parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
        if(verbose>0){print(paste("Validation for rounds for",traits[t]))}
        if(verbose>0){print(Sys.time()-valid.start)}


        for (v in 1:length(validation.out)) {
          all.STmod.out.of.valid.fold.preds[rownames(validation.out[[v]]),,validation.grid[v,2]]<-validation.out[[v]]
        }
        all.STmod.out.of.valid.fold.preds<-apply(all.STmod.out.of.valid.fold.preds,1:2,FUN = function(x) mean(na.omit(x))) #Average accross validaton rounds

        base.model.valid.rs<-apply(all.STmod.out.of.valid.fold.preds,2,
                                   function(x) cor(x,Phenos[rownames(all.STmod.out.of.valid.fold.preds),t],use="pairwise.complete.obs"))

        if(verbose>1){print("Out-of-validation fold prediction accuracies:")}
        if(verbose>1){print(round(base.model.valid.rs,2))}

        all.trait.out.of.valid.fold.preds[[t]]<-all.STmod.out.of.valid.fold.preds

        #Run each ST Gmodel on full training data----
        all.STmod.out.of.test.fold.preds<-matrix(NA,nrow = nrow(Phenos[test.folds[[i]],]),ncol=length(base.models),
                                                 dimnames = list(rownames(Phenos[test.folds[[i]],]),base.models))
        {
          if(verbose>0){print(paste("========Starting",traits[t]," test fold",i," CV round",r,"====Full traing set fitting=========="))}

          #RRBLUP----
          if("rrblup"%in%base.models){
            start<-Sys.time()
            Phenos.na<-Phenos
            Phenos.na[test.folds[[i]],]<-NA
            Phenos.na<-cbind.data.frame(rownames(Phenos.na),Phenos.na)
            colnames(Phenos.na)[1]<-"line"
            A <-rrBLUP::A.mat(Genos)
            gblup.model<-rrBLUP::kin.blup(data=Phenos.na,pheno = traits[t],geno = "line",K = A)
            all.STmod.out.of.test.fold.preds[,"rrblup"]<-gblup.model$pred[rownames(all.STmod.out.of.test.fold.preds)]
            if(verbose>1){print("running rrBLUP")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #EGBLUP----
          if("EGBLUP"%in%base.models){
            start<-Sys.time()
            Phenos.na<-Phenos
            Phenos.na<-Phenos.na[!is.na(Phenos.na[,t]),] #Remove NA rows
            Phenos.na[rownames(Phenos)[test.folds[[i]]],]<-NA #Make test folds genos NA
            A <- rrBLUP::A.mat(Genos)
            H<-matrixcalc::hadamard.prod(A,A)
            Phenos.na<-Phenos.na[rownames(A),]
            ETA<-list(list(K=A,model="RKHS"),list(K=H,model="RKHS"))
            fm<-BGLR::BGLR(y =Phenos.na[,t],ETA = ETA,nIter=5000,burnIn=1000,df0=5,S0=2,verbose = F,saveAt = bglr.saveAt)
            fm.preds<-fm$yHat
            names(fm.preds)<-rownames(A)
            all.STmod.out.of.test.fold.preds[,"EGBLUP"]<-fm.preds[rownames(all.STmod.out.of.test.fold.preds)]
            if(verbose>1){print("running EGBLUP")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #RF----
          if("RF"%in%base.models){
            start<-Sys.time()
            X<-as.matrix(Genos[rownames(Phenos)[-test.folds[[i]]],])[!is.na(Phenos[-test.folds[[i]],t]),]
            rf.fit<-ranger::ranger(x=X,y=Phenos[rownames(X),t],num.tree=rf.ntrees,max.depth=NULL
                                   ,min.node.size=15,mtry = round(ncol(X)/3),num.threads = n.cores)
            all.STmod.out.of.test.fold.preds[,"RF"]<-predict(rf.fit,data = as.matrix(Genos[rownames(all.STmod.out.of.test.fold.preds),]))$predictions
            if(verbose>1){print("running RF")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #GBM----
          if("GBM"%in%base.models){
            start<-Sys.time()
            X<-as.matrix(Genos[rownames(Phenos)[-test.folds[[i]]],])[!is.na(Phenos[-test.folds[[i]],t]),]
            X<-X[sample(1:nrow(X)),]
            Y<-Phenos[rownames(X),t]
            gbmdata<-cbind.data.frame(Y,X)
            gbm.tune<-expand.grid(depth=c(2:6),bag.fraction=c(0.4,0.5,0.6,0.7,0.8),minobs=round(seq(2,nrow(X)*0.05,length.out=5)))
            gbm.tune<-gbm.tune[sample(1:nrow(gbm.tune)),]
            gbm.tune.error<-c()
            if(verbose>1){print("GBM param tune")}
            for (h in 1:10){
              if(verbose>1){cat("|",sep="")}
              gbm.fit<-gbm::gbm(Y~.,data = gbmdata,distribution = "gaussian",train.fraction = 0.8
                                ,interaction.depth = gbm.tune$depth[h],bag.fraction =gbm.tune$bag.fraction[h],
                                n.minobsinnode = gbm.tune$minobs[h],n.trees = 50,shrinkage = 0.2,verbose = F,n.cores=1)
              gbm.tune.error[h]<-min(gbm.fit$valid.error)
            }
            for (h in 1:5){
              if(verbose>1){cat("|",sep="")}
              tune.model<-randomForest::randomForest(x=gbm.tune[which(!is.na(gbm.tune.error)),],y=gbm.tune.error[!is.na(gbm.tune.error)],ntree=500,importance=F)
              full.error.preds<-predict(tune.model,newdata = gbm.tune)
              full.error.preds[which(!is.na(gbm.tune.error))]<-NA
              next.index<-which.min(full.error.preds)
              gbm.fit<-gbm::gbm(Y~.,data = gbmdata,distribution = "gaussian",train.fraction = 0.8
                                ,interaction.depth = gbm.tune$depth[next.index],bag.fraction =gbm.tune$bag.fraction[next.index],
                                n.minobsinnode = gbm.tune$minobs[next.index],n.trees = 50,shrinkage = 0.2,verbose = F,n.cores=1)
              gbm.tune.error[next.index]<-min(gbm.fit$valid.error)
            }
            opt.params<-gbm.tune[which.min(gbm.tune.error),]
            gbm.cores<-min(n.cores,5)
            gbm.fit<-gbm::gbm(Y~.,data = gbmdata,cv.folds = 5,distribution = "gaussian",interaction.depth = opt.params$depth,n.minobsinnode = opt.params$minobs,
                              n.trees = gbm.ntrees,shrinkage = 0.03,verbose = F,n.cores = gbm.cores)
            opt.tree<-gbm::gbm.perf(gbm.fit,plot.it = F,method = "cv")
            doParallel::stopImplicitCluster()
            if(verbose>1){print(paste("GBM opt run tree no. =",opt.tree))}
            all.STmod.out.of.test.fold.preds[,"GBM"]<-gbm::predict.gbm(gbm.fit,n.trees = opt.tree,newdata = as.data.frame(Genos[rownames(all.STmod.out.of.test.fold.preds),]))
            if(verbose>1){print("running GBM")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #LASSO----
          if("LASSO"%in%base.models){
            start<-Sys.time()
            X<-as.matrix(Genos[rownames(Phenos)[-test.folds[[i]]],])[!is.na(Phenos[-test.folds[[i]],t]),]
            cv.fit<-glmnet::cv.glmnet(x=X,y=Phenos[rownames(X),t],nfolds = 8,lambda = seq(0.000001,0.9,length.out=50)^4)
            lasso.fit<-glmnet::glmnet(x=X,y=Phenos[rownames(X),t],lambda = cv.fit$lambda.min,alpha = 1)
            all.STmod.out.of.test.fold.preds[,"LASSO"]<-predict(lasso.fit,newx = as.matrix(Genos[rownames(all.STmod.out.of.test.fold.preds),]))
            if(verbose>1){print("running LASSO")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #RKHS----
          if("RKHS"%in%base.models){
            start<-Sys.time()
            Phenos.na<-Phenos
            Phenos.na<-Phenos.na[!is.na(Phenos.na[,t]),] #Remove NA rows
            Phenos.na[rownames(Phenos)[test.folds[[i]]],]<-NA #Make test folds genos NA
            D<-as.matrix(dist(Genos[rownames(Phenos.na),],method="euclidean"))^2 #Compute Gaussian kernel
            D<-D/mean(D)
            h<-0.5*c(1/5,1,5)
            ETA<-list(list(K=exp(-h[1]*D),model='RKHS'),
                      list(K=exp(-h[2]*D),model='RKHS'),
                      list(K=exp(-h[3]*D),model='RKHS'))
            Phenos.na<-Phenos.na[rownames(D),]
            fm<-BGLR::BGLR(y =Phenos.na[,t],ETA = ETA,nIter=5000,burnIn=1000,df0=5,S0=2,verbose = F,saveAt = bglr.saveAt)
            fm.preds<-fm$yHat
            names(fm.preds)<-rownames(D)
            all.STmod.out.of.test.fold.preds[,"RKHS"]<-fm.preds[rownames(all.STmod.out.of.test.fold.preds)]
            if(verbose>1){print("running RKHS")}
            if(verbose>1){print(Sys.time()-start)}
          }
        }
        if(verbose>1){print("Correlations among out-of-test-fold predictions:")}
        if(verbose>1){print(round(cor(all.STmod.out.of.test.fold.preds),2))}

        all.trait.out.of.test.fold.preds[[t]]<-all.STmod.out.of.test.fold.preds
      } #End of traits loop


      #MT model ensembles----
      if(verbose>0){print(paste("=============MT ensemble training for test fold",i," CV round",r,"================"))}
      X<-matrix(unlist(all.trait.out.of.valid.fold.preds),
                nrow = nrow(all.trait.out.of.valid.fold.preds[[1]]),
                ncol=length(base.models)*length(traits),
                byrow = F,dimnames = list(rownames(all.trait.out.of.valid.fold.preds[[1]]),
                                          paste(expand.grid(base.models,traits)[,1],
                                                expand.grid(base.models,traits)[,2])[1:(length(all.trait.out.of.valid.fold.preds)*length(base.models))]))

      newX<-matrix(unlist(all.trait.out.of.test.fold.preds),
                   nrow = nrow(all.trait.out.of.test.fold.preds[[1]]),
                   ncol=length(base.models)*length(traits),
                   byrow = F,dimnames = list(rownames(all.trait.out.of.test.fold.preds[[1]]),
                                             paste(expand.grid(base.models,traits)[,1],
                                                   expand.grid(base.models,traits)[,2])[1:(length(all.trait.out.of.test.fold.preds)*length(base.models))]))

      for (t in 1:length(traits)) {
        if(verbose>1){print(paste("GA optimisation for",traits[t]))}
        Y<-Phenos[rownames(X),t]
        Y<-Phenos[rownames(X),t]
        fitness <- function(params) {
          xsubsub<-X[!is.na(Y),]
          ysubsub<-Y[!is.na(Y)]
          rfcvs<-split(sample(c(1:nrow(xsubsub))),1:8)
          cv.preds<-c()
          for(cv in 1:length(rfcvs)){
            rfmod<-ranger::ranger(x=xsubsub[-rfcvs[[cv]],],y=ysubsub[-rfcvs[[cv]]],num.tree=300,mtry = floor(ncol(xsubsub)*params[1])
                                  ,min.node.size = floor(params[2]),oob.error = F)
            cv.preds[rfcvs[[cv]]]<-predict(rfmod,data=xsubsub[rfcvs[[cv]],])$predictions
          }
          cor(cv.preds,ysubsub)
        }
        GA <- GA::ga("real-valued",lower = c(0.05,2),upper = c(0.95,floor(nrow(X)*0.2)),
                     fitness = fitness, maxiter = 20,keepBest = T,popSize = 10,elitism = 5,
                     parallel = n.cores,monitor = F)
        opt.params<- GA@solution
        rf.fit<-ranger::ranger(x=X[!is.na(Y),],y=Y[!is.na(Y)],num.tree=400,max.depth=NULL
                               ,min.node.size=floor(opt.params[1,2]),mtry = floor(ncol(X)*opt.params[1,1]),oob.error = F)
        MTens.prediction<-predict(rf.fit,data=newX)$predictions
        names(MTens.prediction)<-row.names(newX)
        all.trait.out.of.test.fold.preds[[t]]<-cbind(all.trait.out.of.test.fold.preds[[t]],"MT ensemble"=MTens.prediction)

        if(verbose>1){print(paste("Test fold",i,"prediction accuracy:"))}
        if(verbose>1){print(round(apply(all.trait.out.of.test.fold.preds[[t]],2,function(x) cor(x,Phenos[names(x),t],use="pairwise.complete.obs")),3))}
      } #End of traits loop
      all.out.of.fold.preds[[i]]<-all.trait.out.of.test.fold.preds
      if(verbose>1){print(paste("===========FINISHED test fold",i," CV round",r,"==========="))}
    }
    all.CV.rounds[[r]]<-all.out.of.fold.preds

    if(verbose>1){print(paste("===========FINISHED CV round",r,"==========="))}
  }

  fold.names<-lapply(all.CV.rounds,function(r)
    c(lapply(r,function(i) rownames(i[[1]]))))
  names(fold.names)<-paste("Round",1:length(fold.names))
  for (r in 1:length(fold.names)) {
    names(fold.names[[r]])<-paste("CV fold",1:n.test.folds)
  }

  base.preds<-array(NA,dim=c(length(rownames(Phenos)),length(base.models)+1,length(traits),length(all.CV.rounds))
                    ,dimnames=list(rownames(Phenos),c(base.models,"MT.ensemble"),traits,1:length(all.CV.rounds)))
  for(r in 1:length(all.CV.rounds)){
    for(t in 1:length(traits)){
      trait.mat<-as.matrix(data.table::rbindlist(lapply(all.CV.rounds[[r]],function(x) data.frame(x[[t]]))))
      rownames(trait.mat)<-unlist(fold.names[[r]])
      base.preds[rownames(trait.mat),,t,r]<- trait.mat
    }
  }

  all.traits.pred.r<-list()
  for(t in 1:length(traits)){
    base.rs<-matrix(NA,nrow=length(all.CV.rounds),ncol=length(base.models)+1,
                    dimnames = list(1:length(all.CV.rounds),c(base.models,"MT.ensemble")))

    for(i in 1:ncol(base.rs)){
      for(r in 1:length(all.CV.rounds)) {
        base.rs[r,i]<-cor(Phenos[,t],base.preds[,i,t,r], use="pairwise.complete.obs")
      }
    }

    all.traits.pred.r[[t]]<-base.rs
  }
  names(all.traits.pred.r)<-traits

  all.out<-list("All out-of-fold predictions"=base.preds,
                "phenotype data"=Phenos,
                "Prediction accuracies"=all.traits.pred.r,
                "Base models"=base.models,
                "Traits"=traits,
                "CV fold allocations"=fold.names)
  if(verbose>1){print("FINISHED!!!")}
  if(verbose>1){print(Sys.time()-run.start.time)}
  return(all.out)}








#' Train and predict multi-trait ensembles
#'
#' @description Trains several base genomic prediction models as well as multiple genetic model and multi-trait ensemble prediction models. Minimum input requirements are SNP genotype data and multi-trait phenotype data base prediction models should run well on default parameters.
#'
#' @param Genos.train NxM marker \emph{matrix} for genotypes in the training set with row names and column names as output from marker.QC(). Markers coded as 0,2 for each SNP allele. Inbred lines are assumed. Must have no missing (NA)
#' values. Row names are Genotype ID names and column are SNP marker names.
#' @param Phenos Data frame of N x p trait data for genotypes in the training set. Row names are Genotype ID matching row names of Genos. Column names are trait names. Missing data is allowed.
#' @param New.Genos NxM marker \emph{matrix} for new genotypes to be predicted with row names and column names as output from marker.QC(). Markers coded as 0,2 for each SNP allele. Inbred lines are assumed. Must have no missing (NA)
#' values. Row names are Genotype ID names and column are SNP marker names.
#' @param base.models \emph{vector} of base model names to use. Options include:
#' c("rrblup","EGBLUP","RKHS","LASSO","RF","GBM","XGB")
#'
#'-\strong{rrblup} Ridge regression/GBLUP implemented in the rrblup package using a linear GRM.
#'
#'-\strong{EGBLUP} The Extended GBLUP incorporating a linear and epistatic GRM kernel implemented in the BGLR package.
#'
#'-\strong{RKHS} Replicating Kernel Hilbert Space fitting multiple Gaussian kernels using the kernel averaging method for several bandwidth parameters. Implemented in the BGLR package.
#'
#'-\strong{LASSO} Linear Shrinkage model using the LASSO penalty to shrink most SNP effects to 0. Implemented in the glmnet package. Each model is cross validated with several values of lambda.
#'
#'-\strong{RF} Random forest decision tree ensemble learning method implemented in the randomForest package. Forests of 300 trees are grown as defualt but can be increased using rf.ntrees if needed.
#'
#'-\strong{GBM} Gradient Boosting Machine learning models implemented in the gbm package. May take a long time... 5 fold cross validation is performed on each model to
#'     determine optimum tree iteration to use. cross validations can be run in parallel on multiple cores by setting n.cores if not defined by default.
#'
#'-\strong{XGB} Extreme Gradient Boosting machine learning implemented in the xgboost package. A hyperparameter grid search is performed accross different tree depth and shrinkage
#'     parameters and each model is 5-fold cross validated  to determine optimum tree iteration. May take a long time...
#'
#' @param n.valid.folds (integer) number of validation folds to run for each trait and base model. Default is 8.
#' @param n.CV.rounds \emph{integer} number of rounds of k-fold cross validation for test folds to perform. Default is 2
#' @param n.cores Number of cores that can be used to run cross validation of RF, GBM and XGB models. The Default detectcores() from Parallel will use available multiple
#' cores on a Windows OS but n.cores should be defined when running on a HPC or other OS.
#' @param rf.ntrees (integer) Number of trees to use for each Random forest model
#' @param gbm.ntrees (integer) Max number of trees to use for each GBM model
#' @param verbose Level of progress to print throughout. 0 = none; 1 = some; 2 = all.
#' @param bglr.saveAt file path location to save temporary data files for 'BGLR' model functions. Using a local folder may speed up computing time. Default is the current working directory.
#'
#'
#' @return Returns a list object of:
#'
#'-\strong{All new predictions} \emph{array} Predictions for all new genotypes for base prediction models. Array with dimensions of genotypes,  base models and traits.
#'
#'-\strong{Traits} \emph{vector} of trait names used.
#'
#'-\strong{Base models} \emph{vector} of base prediction models used.
#'
#'@details Ensemble predictions for for genotypes in New.Genos for multiple genetic models and multiple traits are made after trained on out-of-validation fold predictions using several model
#' approaches in the Genos.train and Phenos objects. Including a range of genetic models that capture a diversity of genetic architectures may improve accuracy of ensembles but running all available base models takes a long time. Unbalanced phenotyping across
#' genotypes may increase accuracy of multi-trait ensembles when several related but partially phenotyped traits are included.
#'
#' @examples
#' library(MTensGS)
#' data(TGdata)
#'
#' QCd.markers<-marker.QC(marker.mat = TGgenos)
#'
#' #Split the dataset in two
#' old.genos<-QCd.markers[1:300,]
#' new.genos<-QCd.markers[301:376,]
#' old.phens<-TGphenos[1:300,]
#'
#' new_preds<-Train.Predict.MTens(Genos.train = old.genos,
#'                                New.Genos = new.genos,
#'                                Phenos = old.phens,
#'                                n.valid.rounds = 1,
#'                                n.valid.folds = 6,
#'                                base.models = c("rrblup","LASSO"))
#'
#'
#'
#' @export
Train.Predict.MTens<-function(Genos.train,New.Genos,Phenos,
                              n.valid.folds=8,
                              n.valid.rounds=3,
                              rf.ntrees=400,
                              gbm.ntrees=300,
                              n.cores=NULL,
                              bglr.saveAt=NULL,
                              base.models=c("rrblup","EGBLUP","RKHS","LASSO","RF","GBM"),
                              verbose=2)
{
  run.start.time<-Sys.time()
  if(length(rownames(Genos.train)[!rownames(Genos.train)%in%rownames(Phenos)])>0){
    print("Lines in Genos but not in Phenos:")
    print(rownames(Genos.train)[!rownames(Genos.train)%in%rownames(Phenos)])}
  if(length(rownames(Phenos)[!rownames(Phenos)%in%rownames(Genos.train)])>0){
    print("Lines in Phenos but not in Genos:")
    print(rownames(Phenos)[!rownames(Phenos)%in%rownames(Genos.train)])}

  if(is.null(n.cores)){
    n.cores<-parallel::detectCores()
  }
  if(is.null(bglr.saveAt)){
    bglr.saveAt<-tempdir()
  }

  Genos.train<-as.matrix(Genos.train[rownames(Genos.train)%in%rownames(Phenos),])
  Phenos<-as.data.frame(Phenos[rownames(Phenos)%in%rownames(Genos.train),])

  traits<-colnames(Phenos)
  combined.Genos<-rbind(Genos.train,New.Genos)

  valid.folds<-list()
  for (vr in 1:n.valid.rounds) {
    valid.folds[[vr]]<-split(sample(c(1:nrow(Phenos))),1:n.valid.folds)#make folds
  }

  all.trait.out.of.test.fold.preds<-list()
  all.trait.out.of.valid.base.model.preds<-list()
  for (t in 1:length(traits)) {
    if(verbose>1){print(paste("Starting",traits[t]))}
    #CV each ST model for ensemble base.models----
    all.STmod.out.of.valid.fold.preds<-array(NA,dim = c(nrow(Phenos),length(base.models),n.valid.rounds),
                                             dimnames = list(rownames(Phenos),base.models,1:n.valid.rounds))
    valid.start<-Sys.time()
    validation.grid<-expand.grid(1:n.valid.folds,1:n.valid.rounds)
    file.remove("Validation folds log.txt")
    cl <- parallel::makeCluster(n.cores,outfile="Validation folds log.txt")
    doParallel::registerDoParallel(cl)
    validation.out<-foreach::foreach(v=1:nrow(validation.grid),.combine = list,.errorhandling="stop",.verbose = F,.multicombine = T) %dopar% {
      print(paste("========Starting",traits[t],"validation fold",v,"=========="))
      valid.names<-rownames(Phenos)[valid.folds[[validation.grid[v,2]]][[validation.grid[v,1]]]]
      valid.fold.preds<-matrix(data=NA,nrow=length(valid.names),ncol=length(base.models),dimnames = list(valid.names,base.models))

      #RRBLUP----
      if("rrblup"%in%base.models){
        start<-Sys.time()
        Phenos.na<-Phenos
        Phenos.na[valid.names,]<-NA
        Phenos.na<-cbind.data.frame(rownames(Phenos.na),Phenos.na)
        colnames(Phenos.na)[1]<-"line"

        A <-rrBLUP::A.mat(Genos.train)
        A<-A[rownames(Phenos),rownames(Phenos)]
        gblup.model<-rrBLUP::kin.blup(data=Phenos.na,pheno = traits[t],geno = "line",K = A)
        valid.fold.preds[valid.names,"rrblup"]<-gblup.model$pred[rownames(valid.fold.preds)][valid.names]
        if(verbose>0){print("running rrBLUP")}
        if(verbose>0){print(Sys.time()-start)}
      }

      #EGBLUP----
      if("EGBLUP"%in%base.models){
        start<-Sys.time()
        Phenos.na<-Phenos
        Phenos.na<-Phenos.na[!is.na(Phenos.na[,t])|rownames(Phenos.na)%in%valid.names,]#Remove lines with na but not in validation fold
        Phenos.na[valid.names,]<-NA
        A <- rrBLUP::A.mat(Genos.train[rownames(Phenos.na),])
        H<-matrixcalc::hadamard.prod(A,A)
        Phenos.na<-Phenos.na[rownames(A),]
        ETA<-list(list(K=A,model="RKHS"),list(K=H,model="RKHS"))
        fm<-BGLR::BGLR(y =Phenos.na[,t],ETA = ETA,nIter=3000,burnIn=500,df0=5,S0=2,verbose = F,saveAt = bglr.saveAt)
        fm.preds<-fm$yHat
        names(fm.preds)<-rownames(A)
        valid.fold.preds[valid.names,"EGBLUP"]<-fm.preds[rownames(valid.fold.preds)][valid.names]
        if(verbose>0){print("running EGBLUP")}
        if(verbose>0){print(Sys.time()-start)}
      }

      #RF----
      if("RF"%in%base.models){
        start<-Sys.time()
        X<-as.matrix(Genos.train[rownames(Phenos),])[!is.na(Phenos[,t]),]
        X<-X[!rownames(X)%in%valid.names,]
        rf.fit<-ranger::ranger(x=X,y=Phenos[rownames(X),t],num.tree=rf.ntrees,max.depth=NULL
                               ,min.node.size=15,mtry = round(ncol(X)/3),num.threads = 1)
        valid.fold.preds[valid.names,"RF"]<-predict(rf.fit,data = as.matrix(Genos.train[valid.names,]))$predictions
        if(verbose>0){print("running RF")}
        if(verbose>0){print(Sys.time()-start)}
      }

      #GBM----
      if("GBM"%in%base.models){
        start<-Sys.time()
        X<-as.matrix(Genos.train[rownames(Phenos),])[!is.na(Phenos[,t]),]
        X<-X[!rownames(X)%in%valid.names,]
        Y<-Phenos[rownames(X),t]
        gbmdata<-cbind.data.frame(Y,X)
        gbm.tune<-expand.grid(depth=c(2:8),
                              bag.fraction=c(0.4,0.5,0.6,0.7,0.8,0.9),
                              minobs=round(seq(2,nrow(X)*0.05,length.out=5)))
        gbm.tune<-gbm.tune[sample(1:nrow(gbm.tune)),]
        gbm.tune.error<-c()
        if(verbose>1){print("GBM param tune")}
        for (h in 1:10){
          if(verbose>1){cat("|",sep="")}
          gbm.fit<-gbm::gbm(Y~.,data = gbmdata,distribution = "gaussian",train.fraction = 0.8
                            ,interaction.depth = gbm.tune$depth[h],bag.fraction =gbm.tune$bag.fraction[h],
                            n.minobsinnode = gbm.tune$minobs[h],n.trees = 50,shrinkage = 0.2,verbose = F,n.cores=1)
          gbm.tune.error[h]<-min(gbm.fit$valid.error)
        }
        for (h in 1:5){
          if(verbose>1){cat("|",sep="")}
          tune.model<-randomForest::randomForest(x=gbm.tune[which(!is.na(gbm.tune.error)),],y=gbm.tune.error[!is.na(gbm.tune.error)],ntree=500,importance=F)
          full.error.preds<-predict(tune.model,newdata = gbm.tune)
          full.error.preds[which(!is.na(gbm.tune.error))]<-NA
          next.index<-which.min(full.error.preds)
          gbm.fit<-gbm::gbm(Y~.,data = gbmdata,distribution = "gaussian",train.fraction = 0.8
                            ,interaction.depth = gbm.tune$depth[next.index],bag.fraction =gbm.tune$bag.fraction[next.index],
                            n.minobsinnode = gbm.tune$minobs[next.index],n.trees = 50,shrinkage = 0.2,verbose = F,n.cores=1)
          gbm.tune.error[next.index]<-min(gbm.fit$valid.error)
        }
        opt.params<-gbm.tune[which.min(gbm.tune.error),]

        gbm.fit<-gbm::gbm(Y~.,data = gbmdata,cv.folds = 5,distribution = "gaussian",interaction.depth = opt.params$depth,n.minobsinnode = opt.params$minobs,
                          n.trees = gbm.ntrees,shrinkage = 0.03,verbose = F,n.cores=1)
        opt.tree<-gbm::gbm.perf(gbm.fit,plot.it = F,method = "cv")
        valid.fold.preds[valid.names,"GBM"]<-gbm::predict.gbm(gbm.fit,n.trees = opt.tree,newdata = as.data.frame(Genos.train[valid.names,]))
        if(verbose>0){print("running GBM")}
        if(verbose>0){print(Sys.time()-start)}
      }

      #LASSO----
      if("LASSO"%in%base.models){
        start<-Sys.time()
        X<-as.matrix(Genos.train[rownames(Phenos),])[!is.na(Phenos[,t]),]
        X<-X[!rownames(X)%in%valid.names,]
        cv.fit<-glmnet::cv.glmnet(x=X,y=Phenos[rownames(X),t],nfolds = 5,lambda = seq(0.000001,0.9,length.out=50)^4)
        lasso.fit<-glmnet::glmnet(x=X,y=Phenos[rownames(X),t],lambda = cv.fit$lambda.min,alpha = 1)
        valid.fold.preds[valid.names,"LASSO"]<-predict(lasso.fit,newx = as.matrix(Genos.train[valid.names,]))
        if(verbose>0){print("running LASSO")}
        if(verbose>0){print(Sys.time()-start)}
      }

      #RKHS----
      if("RKHS"%in%base.models){
        start<-Sys.time()
        Phenos.na<-Phenos
        Phenos.na<-Phenos.na[!is.na(Phenos.na[,t])|rownames(Phenos.na)%in%valid.names,]#Remove lines with na but not in validation fold
        Phenos.na[valid.names,]<-NA
        Genos.sub<-Genos.train[rownames(Phenos.na),]
        D<-as.matrix(dist(Genos.sub,method="euclidean"))^2 #Compute Gaussian kernel
        D<-D/mean(D)
        h<-0.5*c(1/5,1,5)
        ETA<-list(list(K=exp(-h[1]*D),model='RKHS'),
                  list(K=exp(-h[2]*D),model='RKHS'),
                  list(K=exp(-h[3]*D),model='RKHS'))

        Phenos.na<-Phenos.na[rownames(D),]
        fm<-BGLR::BGLR(y =Phenos.na[,t],ETA = ETA,nIter=3000,burnIn=500,df0=5,S0=2,verbose = F,saveAt = bglr.saveAt)
        fm.preds<-fm$yHat
        names(fm.preds)<-rownames(D)
        valid.fold.preds[valid.names,"RKHS"]<-fm.preds[rownames(valid.fold.preds)][valid.names]
        if(verbose>0){print("running RKHS")}
        if(verbose>0){print(Sys.time()-start)}
      }
      return(valid.fold.preds)
    }
    parallel::stopCluster(cl)
    doParallel::stopImplicitCluster()
    if(verbose>0){print(paste("Validation for rounds for",traits[t]))}
    if(verbose>0){print(Sys.time()-valid.start)}

    for (v in 1:length(validation.out)) {
      all.STmod.out.of.valid.fold.preds[rownames(validation.out[[v]]),,validation.grid[v,2]]<-validation.out[[v]]
    }
    all.STmod.out.of.valid.fold.preds<-apply(all.STmod.out.of.valid.fold.preds,1:2,FUN = function(x) mean(na.omit(x))) #Average accross validaton rounds

    all.trait.out.of.valid.base.model.preds[[t]]<-all.STmod.out.of.valid.fold.preds

    if(verbose>0){print("Out-of-validation fold prediction accuracies:")}
    if(verbose>0){print(apply(all.STmod.out.of.valid.fold.preds,2,function(x) cor(x,Phenos[,t],use="pairwise.complete.obs")))}


    all.STmod.out.of.test.fold.preds<-matrix(NA,nrow = nrow(New.Genos),ncol=length(base.models),
                                             dimnames = list(rownames(New.Genos),base.models))

    {
      if(verbose>0){print(paste("========Starting",traits[t],"====Full training set fitting=========="))}

      #RRBLUP----
      if("rrblup"%in%base.models){
        start<-Sys.time()
        Phenos.na<-Phenos
        Phenos.na<-cbind.data.frame(rownames(Phenos.na),Phenos.na)
        colnames(Phenos.na)[1]<-"line"
        A <-rrBLUP::A.mat(combined.Genos)
        Phenos.na<-Phenos.na[rownames(A),]
        gblup.model<-rrBLUP::kin.blup(data=Phenos.na,pheno = traits[t],geno = "line",K = A)
        all.STmod.out.of.test.fold.preds[,"rrblup"]<-gblup.model$pred[rownames(all.STmod.out.of.test.fold.preds)]
        if(verbose>0){print("running rrBLUP")}
        if(verbose>0){print(Sys.time()-start)}
      }

      #EGBLUP----
      if("EGBLUP"%in%base.models){
        start<-Sys.time()
        Phenos.na<-Phenos[!is.na(Phenos[,t]),]
        A <- rrBLUP::A.mat(combined.Genos[unique(c(rownames(Phenos.na),rownames(New.Genos))),])
        H<-matrixcalc::hadamard.prod(A,A)
        Phenos.na<-Phenos.na[rownames(A),]
        ETA<-list(list(K=A,model="RKHS"),list(K=H,model="RKHS"))
        fm<-BGLR::BGLR(y =Phenos.na[,t],ETA = ETA,nIter=5000,burnIn=1000,df0=5,S0=2,verbose = F)
        fm.preds<-fm$yHat
        names(fm.preds)<-rownames(A)
        all.STmod.out.of.test.fold.preds[,"EGBLUP"]<-fm.preds[rownames(all.STmod.out.of.test.fold.preds)]
        if(verbose>0){print("running EGBLUP")}
        if(verbose>0){print(Sys.time()-start)}
      }

      #RF----
      if("RF"%in%base.models){
        start<-Sys.time()
        X<-as.matrix(Genos.train[rownames(Phenos),])[!is.na(Phenos[,t]),]
        rf.fit<-ranger::ranger(x=X,y=Phenos[rownames(X),t],num.tree=rf.ntrees,max.depth=NULL
                               ,min.node.size=15,mtry = round(ncol(X)/3),num.threads = 1)
        all.STmod.out.of.test.fold.preds[,"RF"]<-predict(rf.fit,data = as.matrix(New.Genos))$predictions
        if(verbose>0){print("running RF")}
        if(verbose>0){print(Sys.time()-start)}
      }

      #GBM----
      if("GBM"%in%base.models){
        start<-Sys.time()
        X<-as.matrix(Genos.train[rownames(Phenos),])[!is.na(Phenos[,t]),]
        Y<-Phenos[rownames(X),t]
        gbmdata<-cbind.data.frame(Y,X)
        gbm.tune<-expand.grid(depth=c(2:8),
                              bag.fraction=c(0.4,0.5,0.6,0.7,0.8,0.9),
                              minobs=round(seq(2,nrow(X)*0.05,length.out=5)))
        gbm.tune<-gbm.tune[sample(1:nrow(gbm.tune)),]
        gbm.tune.error<-c()
        if(verbose>1){print("GBM param tune")}
        for (h in 1:10){
          if(verbose>1){cat("|",sep="")}
          gbm.fit<-gbm::gbm(Y~.,data = gbmdata,distribution = "gaussian",train.fraction = 0.8
                            ,interaction.depth = gbm.tune$depth[h],bag.fraction =gbm.tune$bag.fraction[h],
                            n.minobsinnode = gbm.tune$minobs[h],n.trees = 50,shrinkage = 0.2,verbose = F,n.cores=1)
          gbm.tune.error[h]<-min(gbm.fit$valid.error)
        }
        for (h in 1:5){
          if(verbose>1){cat("|",sep="")}
          tune.model<-randomForest::randomForest(x=gbm.tune[which(!is.na(gbm.tune.error)),],y=gbm.tune.error[!is.na(gbm.tune.error)],ntree=500,importance=F)
          full.error.preds<-predict(tune.model,newdata = gbm.tune)
          full.error.preds[which(!is.na(gbm.tune.error))]<-NA
          next.index<-which.min(full.error.preds)
          gbm.fit<-gbm::gbm(Y~.,data = gbmdata,distribution = "gaussian",train.fraction = 0.8
                            ,interaction.depth = gbm.tune$depth[next.index],bag.fraction =gbm.tune$bag.fraction[next.index],
                            n.minobsinnode = gbm.tune$minobs[next.index],n.trees = 50,shrinkage = 0.2,verbose = F,n.cores=1)
          gbm.tune.error[next.index]<-min(gbm.fit$valid.error)
        }
        opt.params<-gbm.tune[which.min(gbm.tune.error),]

        gbm.fit<-gbm::gbm(Y~.,data = gbmdata,cv.folds = 5,distribution = "gaussian",interaction.depth = opt.params$depth,n.minobsinnode = opt.params$minobs,
                          n.trees = gbm.ntrees,shrinkage = 0.03,verbose = F,n.cores=n.cores)
        opt.tree<-gbm::gbm.perf(gbm.fit,plot.it = F,method = "cv")
        all.STmod.out.of.test.fold.preds[,"GBM"]<-gbm::predict.gbm(gbm.fit,n.trees = opt.tree,newdata = as.data.frame(New.Genos[rownames(all.STmod.out.of.test.fold.preds),]))
        if(verbose>0){print("running GBM")}
        if(verbose>0){print(Sys.time()-start)}
      }

      #LASSO----
      if("LASSO"%in%base.models){
        start<-Sys.time()
        X<-as.matrix(Genos.train[rownames(Phenos),])[!is.na(Phenos[,t]),]
        cv.fit<-glmnet::cv.glmnet(x=X,y=Phenos[rownames(X),t],nfolds = 8,lambda = seq(0.000001,0.9,length.out=50)^4)
        lasso.fit<-glmnet::glmnet(x=X,y=Phenos[rownames(X),t],lambda = cv.fit$lambda.min,alpha = 1)
        all.STmod.out.of.test.fold.preds[,"LASSO"]<-predict(lasso.fit,newx = as.matrix(New.Genos[rownames(all.STmod.out.of.test.fold.preds),]))
        if(verbose>0){print("running LASSO")}
        if(verbose>0){print(Sys.time()-start)}
      }

      #RKHS----
      if("RKHS"%in%base.models){
        start<-Sys.time()
        Phenos.na<-Phenos[!is.na(Phenos[,t]),]
        D<-as.matrix(dist(combined.Genos[unique(c(rownames(Phenos.na),rownames(New.Genos))),],method="euclidean"))^2 #Compute Gaussian kernel
        D<-D/mean(D)
        Phenos.na<-Phenos.na[rownames(D),]
        h<-0.5*c(1/5,1,5)
        ETA<-list(list(K=exp(-h[1]*D),model='RKHS'),
                  list(K=exp(-h[2]*D),model='RKHS'),
                  list(K=exp(-h[3]*D),model='RKHS'))
        fm<-BGLR::BGLR(y =Phenos.na[,t],ETA = ETA,nIter=5000,burnIn=1000,df0=5,S0=2,verbose = F)
        fm.preds<-fm$yHat
        names(fm.preds)<-rownames(D)
        all.STmod.out.of.test.fold.preds[,"RKHS"]<-fm.preds[rownames(all.STmod.out.of.test.fold.preds)]
        if(verbose>0){print("running RKHS")}
        if(verbose>0){print(Sys.time()-start)}
      }
    }
    all.trait.out.of.test.fold.preds[[t]]<-all.STmod.out.of.test.fold.preds
  }
  #MT model ensembles----
  if(verbose>0){print(paste("=============MT ensemble training================"))}
  X<-matrix(unlist(all.trait.out.of.valid.base.model.preds),
            nrow = nrow(all.trait.out.of.valid.base.model.preds[[1]]),
            ncol=length(base.models)*length(traits),
            byrow = F,dimnames = list(rownames(all.trait.out.of.valid.base.model.preds[[1]]),
                                      paste(expand.grid(base.models,traits)[,1],
                                            expand.grid(base.models,traits)[,2])[1:(length(all.trait.out.of.valid.base.model.preds)*length(base.models))]))

  newX<-matrix(unlist(all.trait.out.of.test.fold.preds),
               nrow = nrow(all.trait.out.of.test.fold.preds[[1]]),
               ncol=length(base.models)*length(traits),
               byrow = F,dimnames = list(rownames(all.trait.out.of.test.fold.preds[[1]]),
                                         paste(expand.grid(base.models,traits)[,1],
                                               expand.grid(base.models,traits)[,2])[1:(length(all.trait.out.of.test.fold.preds)*length(base.models))]))

  for (t in 1:length(traits)) {
    if(verbose>1){print(paste("GA optimisation for",traits[t]))}
    Y<-Phenos[rownames(X),t]
    fitness <- function(params) {
      xsubsub<-X[!is.na(Y),]
      ysubsub<-Y[!is.na(Y)]
      rfcvs<-split(sample(c(1:nrow(xsubsub))),1:8)
      cv.preds<-c()
      for(cv in 1:length(rfcvs)){
        rfmod<-ranger::ranger(x=xsubsub[-rfcvs[[cv]],],y=ysubsub[-rfcvs[[cv]]],num.tree=300,mtry = floor(ncol(xsubsub)*params[1])
                              ,min.node.size = floor(params[2]),oob.error = F)
        cv.preds[rfcvs[[cv]]]<-predict(rfmod,data=xsubsub[rfcvs[[cv]],])$predictions
      }
      cor(cv.preds,ysubsub)
    }
    GA <- GA::ga("real-valued",lower = c(0.05,2),upper = c(0.95,floor(nrow(X)*0.2)),
                 fitness = fitness, maxiter = 20,keepBest = T,popSize = 10,elitism = 5,
                 parallel = n.cores,monitor = F)
    opt.params<- GA@solution
    rf.fit<-ranger::ranger(x=X[!is.na(Y),],y=Y[!is.na(Y)],num.tree=400,max.depth=NULL
                           ,min.node.size=floor(opt.params[1,2]),mtry = floor(ncol(X)*opt.params[1,1]),oob.error = F)
    MTens.prediction<-predict(rf.fit,data=newX)$predictions
    names(MTens.prediction)<-row.names(newX)
    all.trait.out.of.test.fold.preds[[t]]<-cbind(all.trait.out.of.test.fold.preds[[t]],"MT ensemble"=MTens.prediction)
  } #End of traits loop
  base.preds<-array(NA,dim=c(length(rownames(New.Genos)),length(base.models)+1,length(traits))
                    ,dimnames=list(rownames(New.Genos),c(base.models,"MT.ensemble"),traits))
  for(t in 1:length(traits)){
    base.preds[rownames(all.trait.out.of.test.fold.preds[[t]]),,t]<-all.trait.out.of.test.fold.preds[[t]]
  }

  all.out<-list("All new predictions"=base.preds,
                "Base models"=base.models,
                "Traits"=traits)
  if(verbose>1){print("FINISHED!!!")}
  if(verbose>1){print(Sys.time()-run.start.time)}

  return(all.out)
}





#' Plot multi-trait Principle component selection indices
#'
#' @description Plot a PCA biplot visualisation of multiple predicted traits to aid selection. Positive and negative trade offs among traits can be considered and weightings applied to subsets of desired traits.
#'
#' @param X Object from Train.Predict.MTens()
#' @param SI.traits \emph{vector} of trait names included in X to visualise. Default is all traits.
#' @param weights \emph{vector} of weightings for each trait in SI.traits. Default is equal weightings (all = 1).
#' @param model Prediction model type to use. Options include. Default is the MT.ensemble
#' @param col1,col2 Colours to use for trait and genotype ID names in biplot.
#'
#'@examples
#' library(MTensGS)
#' data(TGdata)
#'
#' QCd.markers<-marker.QC(marker.mat = TGgenos)
#'
#' #Split the dataset in two
#' old.genos<-QCd.markers[1:300,]
#' new.genos<-QCd.markers[301:376,]
#' old.phens<-TGphenos[1:300,]
#'
#' new_preds<-Train.Predict.MTens(Genos.train = old.genos,
#'                                New.Genos = new.genos,
#'                                Phenos = old.phens,
#'                                n.valid.rounds = 1,
#'                                n.valid.folds = 6,
#'                                base.models = c("rrblup","LASSO"))
#'
#' sel.ind<-MT.PCSI.plot(new_preds,SI.traits = c("YLD_GBR_2011","YLD_GBR_2010","YLD_DEU_2010","PROT")
#'                       ,weights = c(1,1,1,2),cex1 = 0.8)
#'
#' @export
MT.PCSI.plot<-function(X,SI.traits=NULL,
                        weights=NULL,
                        model="MT.ensemble",
                        col1=NULL,
                        col2=NULL,
                        cex1=0.7,
                        cex2=1){
  if(is.null(col1)){
    col1<-viridis::magma(10,alpha = 0.7)[2]
  }
  if(is.null(col2)){
    col2<-viridis::magma(10,alpha = 0.7)[7]
  }
  if(is.null(SI.traits)){
    SI.traits<-colnames(preds)
  }
  if(is.null(weights)){
    weights<-rep(1,length(SI.traits))
  }
  preds<-X$`All new predictions`[,model,SI.traits]
  if (!length(SI.traits)==length(weights)) {
    print("SI.traits and weights vectors are not the same length...")
  }
  index<-apply(scale(preds),1,function(x) weighted.mean(x,w = weights))
  weights.mat<-matrix(rep(weights,nrow(preds)),nrow=nrow(preds),ncol = length(SI.traits),byrow = T)
  scaled<-scale(preds)*weights.mat
  pc<-prcomp(scaled)
  par(mar=c(4,4,1,1))
  plot(pc$x[,1],pc$x[,2],pch=16,col=viridis::mako(nrow(preds),alpha = 0.7,direction = -1)[rank(index)],
       xlim=rep(max(abs(range(pc$x[,1]))),2)*c(-1,1),ylim=rep(max(abs(range(pc$x[,2]))),2)*c(-1,1),
       xlab="Dim 1",ylab="Dim 2")
  abline(h=0,lty=3,col="grey")
  abline(v=0,lty=3,col="grey")
  text(pc$x[,1],pc$x[,2],labels = rownames(pc$x),col=col1,cex=cex1,pos=3)
  plot.window(xlim = rep(max(abs(range(pc$rotation[,1]))),2)*c(-1,1),
              ylim =rep(max(abs(range(pc$rotation[,2]))),2)*c(-1,1)
              ,asp = 1)
  text(pc$rotation[,1],pc$rotation[,2],labels = SI.traits,col=col2,cex=cex2)
  arrows(x0 = 0,y0 = 0,pc$rotation[,1]*0.9,pc$rotation[,2]*0.9,length = 0,col=col2)

  return(index)
  }


