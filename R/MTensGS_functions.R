#' SNP marker cleaning and QC
#'
#' @description Uses a SNP marker matrix and sequentially filters for missing data, imputes
#' remaining missing NA SNP alleles per marker,
#' filters for minor allele frequency and 'prunes' markers based on marker correlation.
#'
#' @param marker.mat NxM marker matrix with row names and column names. Markers coded as 0,1,2
#' for hom, het,hom.
#' @param impute Logical whether to run imputation step. Default is TRUE
#' @param cutoff correlation threshold for pruning. Numeric between 0 and 1. Default is 0.8
#' @param min.maf Minimum minor allele frequency to keep. Numeric between 0 and 0.5. Default
#'  is 0.05
#' @param max.NAs Maximum proportion of NA missing data Numeric between 0 and 1. Default is
#' 0.1
#' @param verbose Logical whether to print progress. Default is TRUE.
#'
#' @return  A NxM marker matrix.
#' @export
#'
marker.QC<-function(marker.mat,#NxM marker matrix with row names and column names. Markers coded as 0,1,2 for hom, het,hom.
                    impute=T, #Logical whether to run imputation step. Default is TRUE
                    cutoff=0.8, #correlation threshold for pruning. Numeric between 0 and 1. Default is 0.8
                    min.maf=0.05, #Minimum minor allele frequency to keep. Numeric between 0 and 0.5. Default is 0.05
                    max.NAs=0.1, #Maximum proportion of NA missing data Numeric between 0 and 1. Default is 0.1
                    verbose=T) #Logical whether to print progress. Default is TRUE.
{
  if(verbose==T){if(min.maf>0.2){print("High MAF!")}}
  if(verbose==T){if(max.NAs>0.3){print("High max NAs!")}}
  if(verbose==T){print("Removing high NA markers")}
  marker.mat[marker.mat==1]<-NA #Hets to NAs
  missing<-(apply(marker.mat,2,FUN = function(x) sum(is.na(x))))/nrow(marker.mat) #calculate NA freq
  marker.mat<-marker.mat[,!missing>max.NAs] #remove markers with >10% NAs

  if(verbose==T){print("Calculating marker correlations....")}
  c<-cor(marker.mat,use="pairwise.complete.obs")#Calculate marker cors
  if(impute==T){
    if(verbose==T){print("Imputing missing data")}
    imputed<-marker.mat
    rest.of.data<-marker.mat
    median.impute<-function(x){
      x[is.na(x)]<-median(na.omit(x))
      return(x)
    }
    rest.of.data<-cbind(apply(rest.of.data,2,FUN = median.impute))
    for(i in 1:ncol(marker.mat)){
      rest.of.data.sub<-rest.of.data[,names(abs(c[,i])[order(abs(c[,i]),decreasing=T)[1:100]])]
      rest.of.data.sub<-rest.of.data.sub[,!colnames(rest.of.data.sub)==colnames(marker.mat)[i]]

      isna<-is.na(marker.mat[,i])
      if(sum(isna)>0){
        rfmod<-randomForest::randomForest(x = rest.of.data.sub[!isna,],y = as.factor(marker.mat[!isna,i]),ntree = 20,importance = T)
        imputed[isna,i]<-as.numeric(as.character(predict(rfmod,newdata = rest.of.data.sub[isna,])))
      }
      if(isTRUE(!round((1:ncol(marker.mat)/ncol(marker.mat))*100)[i]==round((1:ncol(marker.mat)/ncol(marker.mat))*100)[i-1])){
        if(verbose==T){cat("|",sep="")}
      }
    }
  }
  if(verbose==T){print("Removing low MAF markers")}
  maf<-apply(imputed,2,function(x) min(sum(na.omit(x==2)),sum(na.omit(x==0)))/nrow(imputed)) #Calculate MAF
  all.markers.mafrm<-imputed[,!maf<min.maf] #Remove MAF<0.05 markers
  if(cutoff<1){
    if(verbose==T){print("Pruning:")}
    to.remove<-c()
    c<-c[colnames(all.markers.mafrm),colnames(all.markers.mafrm)]
    for(i in 1:ncol(all.markers.mafrm)){
      if(!colnames(c)[i]%in%to.remove){
        linked<-rownames(c)[!rownames(c)%in%to.remove][abs(c[!rownames(c)%in%to.remove,i])>cutoff]
        to.remove<-c(to.remove,linked[!linked==colnames(c)[i]])
      }
      if(isTRUE(!round((1:ncol(c)/ncol(c))*100)[i]==round((1:ncol(c)/ncol(c))*100)[i-1])){
        if(verbose==T){cat("|",sep="")}
      }
    }
    to.remove<-unique(to.remove)
    pruned<-all.markers.mafrm[,!colnames(all.markers.mafrm)%in%to.remove]
    out<-pruned
  }else{
    out<-all.markers.mafrm
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
#'-\strong{XGB} Extreme Gradient Boosting machine learning implemented in the xgboost package. A hyperparameter grid search is performed accross different tree depth and shrinkage
#'     parameters and each model is 5-fold cross validated  to determine optimum tree iteration. May take a long time...
#'
#' @param n.test.folds \emph{integer} number of test folds to run. Default if 10.
#' @param n.valid.folds \emph{integer} number of validation folds to run within each test fold. Default if 8.
#' @param n.fold.rounds \emph{integer} number of rounds of k-fold cross validation to perform. Default is 2
#' @param CV.groups \emph{list} of vectors containing Genotype ID names. Cross validation can also be performed across predefined subsets of genotypes.
#' @param n.cores Number of cores that can be used to run cross validation of RF, GBM and XGB models. The Default detectcores() from Parallel will use available multiple
#' cores on a Windows OS but n.cores should be defined when running on a HPC or other OS.
#' @param rf.ntrees \emph{integer} Number of trees to use for each random forest model
#' @param gbm.ntrees \emph{integer} Max number of trees to use for each GBM model
#' @param verbose Level of progress to print throughout. 0 = none; 1 = some; 2 = all.
#'
#' @return Returns a list object of:
#'
#'-\strong{All out-of-fold predictions} for all genotypes for base prediction models, genetic model ensemble and multi-trait ensemble models for each trait.
#'
#'-\strong{Phenotype data} all observed multi-trait input data.
#'
#'-\strong{Prediction accuracies} list of correlation coeficient (r) for each model for each trait.
#'
#'-\strong{Base models} \emph{vector} of base prediction models used.
#'
#'-\strong{Traits} \emph{vector} of trait names used.
#'
#'-\strong{CV fold allocations} A two level list of length = n.fold.rounds, each containing the genotype ID names allocated to each test fold.
#'
#' @details Ensemble predictions for multiple genetic models and multiple traits are trained on out-of-validation fold predictions using several model approaches.
#' Including a range of genetic models that capture a diversity of genetic architectures may improve accuracy of ensembles. Unbalanced phenotyping across
#' genotypes may increase accuracy of multi-trait ensembles when several related but partially phenotyped traits are included.
#'
#'
#' @export
#'
CV.MTens<-function(Genos,Phenos,
                   base.models=c("rrblup","EGBLUP","RKHS","LASSO","RF","GBM","XGB"),
                   n.test.folds=10,
                   n.valid.folds=8,
                   n.fold.rounds=2,
                   CV.groups=NULL,
                   n.cores=detectCores(),
                   rf.ntrees=300,
                   gbm.ntrees=200,
                   verbose=2){
  run.start.time<-Sys.time()
  if(verbose>0){
    if(length(rownames(Genos)[!rownames(Genos)%in%rownames(Phenos)])>0){
      print("Lines in Genos but not in Phenos:")
      print(rownames(Genos)[!rownames(Genos)%in%rownames(Phenos)])}
    if(length(rownames(Phenos)[!rownames(Phenos)%in%rownames(Genos)])>0){
      print("Lines in Phenos but not in Genos:")
      print(rownames(Phenos)[!rownames(Phenos)%in%rownames(Genos)])}
  }
  Genos<-as.matrix(Genos[rownames(Genos)%in%rownames(Phenos),])
  Phenos<-as.data.frame(Phenos[rownames(Phenos)%in%rownames(Genos),])

  if(!is.null(CV.groups)){
    n.fold.rounds<-1
    n.test.folds<-length(CV.groups)
  }

  all.CV.rounds<-list()
  traits<-colnames(Phenos)
  for(r in 1:n.fold.rounds){
    if(verbose>0){print(paste("Starting test CV round",r))}

    if(is.null(CV.groups)){
      test.folds<-split(sample(1:nrow(Phenos)),1:n.test.folds)#make folds
    }else{
      test.folds<-lapply(CV.groups,function(x) which(rownames(all.trait.blups)%in%x) )
    }

    #INSERIES---------------------------
    all.out.of.fold.preds<-list()
    for (i in 1:n.test.folds) {

      all.trait.in.test.fold.Gmodens.preds<-matrix(NA,nrow = nrow(Phenos[-test.folds[[i]],]),length(traits),
                                                   dimnames = list(rownames(Phenos[-test.folds[[i]],]),paste("GMod_ens",traits)))
      all.trait.out.of.test.fold.Gmodens.preds<-matrix(NA,nrow = nrow(Phenos[test.folds[[i]],]),length(traits),
                                                       dimnames = list(rownames(Phenos[test.folds[[i]],]),paste("GMod_ens",traits)))

      all.trait.out.of.test.fold.preds<-list()
      all.trait.in.test.fold.preds<-list()
      all.trait.Gmod.ens<-list()

      valid.folds<-split(sample(c(1:nrow(Phenos))[!1:nrow(Phenos)%in%test.folds[[i]]]),1:n.valid.folds)#make folds

      for (t in 1:length(traits)) {
        if(verbose>0){print(paste("Starting",traits[t]))}
        #CV each ST model for ensemble base.models----
        all.STmod.out.of.valid.fold.preds<-matrix(NA,nrow = nrow(Phenos[-test.folds[[i]],]),ncol=length(base.models),
                                                  dimnames = list(rownames(Phenos[-test.folds[[i]],]),base.models))
        for (v in 1:length(valid.folds)) {
          if(verbose>0){print(paste("========Starting",traits[t],"validation fold",v," test fold",i," CV round",r,"=========="))}
          valid.names<-rownames(Phenos)[valid.folds[[v]]]
          #RRBLUP----
          if("rrblup"%in%base.models){
            start<-Sys.time()
            Phenos.na<-Phenos
            Phenos.na[test.folds[[i]],]<-NA
            Phenos.na[valid.folds[[v]],]<-NA
            Phenos.na<-cbind.data.frame(rownames(Phenos.na),Phenos.na)
            colnames(Phenos.na)[1]<-"line"

            A<-rrBLUP::A.mat(Genos)
            A<-A[rownames(Phenos)[-test.folds[[i]]],rownames(Phenos)[-test.folds[[i]]]]
            gblup.model<-rrBLUP::kin.blup(data=Phenos.na,pheno = traits[t],geno = "line",K = A)
            all.STmod.out.of.valid.fold.preds[valid.names,"rrblup"]<-gblup.model$pred[rownames(all.STmod.out.of.valid.fold.preds)][valid.names]
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
            fm<-BGLR::BGLR(y =Phenos.na[,t],ETA = ETA,nIter=3000,burnIn=500,df0=5,S0=2,verbose = F)
            fm.preds<-fm$yHat
            names(fm.preds)<-rownames(A)
            all.STmod.out.of.valid.fold.preds[valid.names,"EGBLUP"]<-fm.preds[rownames(all.STmod.out.of.valid.fold.preds)][valid.names]
            if(verbose>1){print("running EGBLUP")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #RF----
          if("RF"%in%base.models){
            start<-Sys.time()
            X<-as.matrix(Genos[rownames(Phenos)[-test.folds[[i]]],])[!is.na(Phenos[-test.folds[[i]],t]),]
            X<-X[!rownames(X)%in%valid.names,]
            cl <- parallel::makeCluster(n.cores)
            doParallel::registerDoParallel(cl)
            rf.fit <- foreach::foreach(ntree=rep(round(rf.ntrees/n.cores),n.cores),.combine=combine,.multicombine=TRUE,.packages='randomForest') %dopar%{
              randomForest::randomForest(x=X,y=Phenos[rownames(X),t],ntree=ntree,importance=F)
            }
            parallel::stopCluster(cl)
            all.STmod.out.of.valid.fold.preds[valid.names,"RF"]<-predict(rf.fit,newdata = as.matrix(Genos[valid.names,]))
            if(verbose>1){print("running RF")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #GBM----
          if("GBM"%in%base.models){
            start<-Sys.time()
            X<-as.matrix(Genos[rownames(Phenos)[-test.folds[[i]]],])[!is.na(Phenos[-test.folds[[i]],t]),]
            X<-X[!rownames(X)%in%valid.names,]
            Y<-Phenos[rownames(X),t]
            gbmdata<-cbind.data.frame(Y,X)
            gbm.fit<-gbm::gbm(Y~.,data = gbmdata,cv.folds = 5,distribution = "gaussian",interaction.depth = 3,n.minobsinnode = 5,
                         n.trees = gbm.ntrees,shrinkage = 0.1,verbose = F,n.cores = n.cores)
            doParallel::stopImplicitCluster()
            opt.tree<-gbm::gbm.perf(gbm.fit,plot.it = F,method = "cv")
            all.STmod.out.of.valid.fold.preds[valid.names,"GBM"]<-gbm::predict.gbm(gbm.fit,n.trees = opt.tree,newdata = as.data.frame(Genos[valid.names,]))
            if(verbose>1){print("running GBM")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #XGB----
          if("XGB"%in%base.models){
            start<-Sys.time()
            X<-as.matrix(Genos[rownames(Phenos)[-test.folds[[i]]],])[!is.na(Phenos[-test.folds[[i]],t]),]
            X<-X[!rownames(X)%in%valid.names,]
            Y<-Phenos[rownames(X),t]
            gbmdata<-as.matrix(X)
            xgb.tune<-expand.grid(depth=c(2:5),shrink=c(0.05,0.1,0.2))
            xgb.tune<-xgb.tune[sample(1:nrow(xgb.tune))[1:(nrow(xgb.tune)*0.5)],]
            xgb.tune.error<-c()
            for (h in 1:nrow(xgb.tune)) {
              if(verbose>1){cat("|",sep="")}
              xgb.cv<-xgboost::xgb.cv(data = gbmdata,label = Y,max_depth = xgb.tune$depth[h], eta =xgb.tune$shrink[h], nrounds = 500,early_stopping_rounds = 5,
                             nfold=5, objective = "reg:squarederror",verbose=F)
              xgb.tune.error[h]<-min(xgb.cv$evaluation_log$test_rmse_mean)
              if(verbose>1){cat(xgb.cv$best_iteration)}
            }
            xgb.cv<-xgboost::xgb.cv(data = gbmdata,label = Y,max_depth = xgb.tune$depth[which.min(xgb.tune.error)], eta =xgb.tune$shrink[which.min(xgb.tune.error)], nrounds = 500,early_stopping_rounds = 5,
                           nfold=5, objective = "reg:squarederror",verbose=F)
            xgb.model<-xgboost::xgboost(data = gbmdata,label = Y,max_depth = xgb.tune$depth[which.min(xgb.tune.error)],
                               eta = xgb.tune$shrink[which.min(xgb.tune.error)], nrounds = xgb.cv$best_iteration,
                               objective = "reg:squarederror",verbose=F)
            all.STmod.out.of.valid.fold.preds[valid.names,"XGB"]<-predict(xgb.model,newdata = Genos[valid.names,])
            if(verbose>1){print("running XGB")}
            if(verbose>1){print(Sys.time()-start)}
          }


          #LASSO----
          if("LASSO"%in%base.models){
            start<-Sys.time()
            X<-as.matrix(Genos[rownames(Phenos)[-test.folds[[i]]],])[!is.na(Phenos[-test.folds[[i]],t]),]
            X<-X[!rownames(X)%in%valid.names,]
            cv.fit<-glmnet::cv.glmnet(x=X,y=Phenos[rownames(X),t],nfolds = 6)
            lasso.fit<-glmnet::glmnet(x=X,y=Phenos[rownames(X),t],lambda = cv.fit$lambda.min,alpha = 1)
            all.STmod.out.of.valid.fold.preds[valid.names,"LASSO"]<-predict(lasso.fit,newx = as.matrix(Genos[valid.names,]))
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
            fm<-BGLR::BGLR(y =Phenos.na[,t],ETA = ETA,nIter=3000,burnIn=500,df0=5,S0=2,verbose = F)
            fm.preds<-fm$yHat
            names(fm.preds)<-rownames(D)
            all.STmod.out.of.valid.fold.preds[valid.names,"RKHS"]<-fm.preds[rownames(all.STmod.out.of.valid.fold.preds)][valid.names]
            if(verbose>1){print("running RKHS")}
            if(verbose>1){print(Sys.time()-start)}
          }
        }
        if(verbose>1){print("Out-of-validation fold prediction accuracies:")}
        if(verbose>1){print(apply(all.STmod.out.of.valid.fold.preds,2,function(x) cor(x,Phenos[rownames(all.STmod.out.of.valid.fold.preds),t],use="pairwise.complete.obs")))}

        #train Gmod ensemble model----
        if(verbose>0){print("Fitting Gmod ensemble base.models")}
        #RF.Ensemble---
        Y<-Phenos[-test.folds[[i]],t]
        rf.Gmod.ens.fit<-randomForest::randomForest(x=all.STmod.out.of.valid.fold.preds[!is.na(Y),],y=Y[!is.na(Y)],
                                      ntree=500,importance=F)
        #GBM.Ensemble----
        gbmdata<-cbind.data.frame(Y,all.STmod.out.of.valid.fold.preds)[!is.na(Y),]

        gbm.tune<-expand.grid(depth=c(2:5),min.obs=c(2,4,6),shrink=c(0.1,0.05,0.01,0.005))
        gbm.tune<-gbm.tune[sample(1:nrow(gbm.tune))[1:(nrow(gbm.tune)*0.5)],]
        gbm.tune.error<-c()
        for (h in 1:nrow(gbm.tune)) {
          if(verbose>1){cat("|",sep="")}
          gbm.fit<-gbm::gbm(gbmdata$Y~.,data = gbmdata,cv.folds = 4,
                       interaction.depth = gbm.tune$depth[h],
                       n.minobsinnode = gbm.tune$min.obs[h],
                       n.trees = 500,n.cores = n.cores,
                       shrinkage =gbm.tune$shrink[h],verbose = F,distribution = "gaussian")
          doParallel::stopImplicitCluster()
          gbm.tune.error[h]<-min(gbm.fit$cv.error)
        }
        if(verbose>1){print(":)")}
        if(verbose>1){print(gbm.tune[which.min(gbm.tune.error),])}
        GBM.Gmod.ens.fit<-gbm::gbm(gbmdata$Y~.,data = gbmdata,cv.folds = 8,
                              interaction.depth = gbm.tune$depth[which.min(gbm.tune.error)],
                              n.minobsinnode = gbm.tune$min.obs[which.min(gbm.tune.error)],
                              n.trees = 500,n.cores = n.cores,
                              shrinkage =gbm.tune$shrink[which.min(gbm.tune.error)],verbose = F,distribution = "gaussian")
        doParallel::stopImplicitCluster()
        GBM.opt.tree<-gbm::gbm.perf(GBM.Gmod.ens.fit,plot.it = F,method = "cv")

        #LASSO Ensemble----
        cv.fit.ens<-glmnet::cv.glmnet(x=all.STmod.out.of.valid.fold.preds[!is.na(Y),],y=Y[!is.na(Y)],nfolds = 8)
        lasso.fit.ens<-glmnet::glmnet(x=all.STmod.out.of.valid.fold.preds[!is.na(Y),],y=Y[!is.na(Y)],lambda = cv.fit.ens$lambda.min,alpha = 1)

        #PCR Ensemble----
        PCA<-prcomp(all.STmod.out.of.valid.fold.preds)
        pcr.lmod<-lm(Phenos[-test.folds[[i]],t]~.,data=as.data.frame(PCA$x))
        beta.Z <- as.matrix(pcr.lmod$coefficients[-1])
        V <- as.matrix(PCA$rotation)
        pcr.beta.X <- V %*% beta.Z


        #Run each ST Gmodel on full training data----
        all.STmod.out.of.test.fold.preds<-matrix(NA,nrow = nrow(Phenos[test.folds[[i]],]),ncol=length(base.models),
                                                 dimnames = list(rownames(Phenos[test.folds[[i]],]),base.models))
        all.STmod.in.test.fold.preds<-matrix(NA,nrow = nrow(Phenos[-test.folds[[i]],]),ncol=length(base.models),
                                             dimnames = list(rownames(Phenos[-test.folds[[i]],]),base.models))

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
            all.STmod.in.test.fold.preds[,"rrblup"]<-gblup.model$pred[rownames(all.STmod.in.test.fold.preds)]
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
            fm<-BGLR::BGLR(y =Phenos.na[,t],ETA = ETA,nIter=5000,burnIn=1000,df0=5,S0=2,verbose = F)
            fm.preds<-fm$yHat
            names(fm.preds)<-rownames(A)
            all.STmod.out.of.test.fold.preds[,"EGBLUP"]<-fm.preds[rownames(all.STmod.out.of.test.fold.preds)]
            all.STmod.in.test.fold.preds[,"EGBLUP"]<-fm.preds[rownames(all.STmod.in.test.fold.preds)]
            if(verbose>1){print("running EGBLUP")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #RF----
          if("RF"%in%base.models){
            start<-Sys.time()
            X<-as.matrix(Genos[rownames(Phenos)[-test.folds[[i]]],])[!is.na(Phenos[-test.folds[[i]],t]),]
            cl <- Parallel::makeCluster(n.cores)
            doParallel::registerDoParallel(cl)
            rf.fit <-foreach::foreach(ntree=rep(round(rf.ntrees/n.cores),n.cores),.combine=combine,.multicombine=TRUE,.packages='randomForest') %dopar%{
            randomForest::randomForest(x=X,y=Phenos[rownames(X),t],ntree=ntree,importance=F)
            }
            stopCluster(cl)
            all.STmod.out.of.test.fold.preds[,"RF"]<-predict(rf.fit,newdata = as.matrix(Genos[rownames(all.STmod.out.of.test.fold.preds),]))
            all.STmod.in.test.fold.preds[,"RF"]<-predict(rf.fit,newdata = as.matrix(Genos[rownames(all.STmod.in.test.fold.preds),]))
            if(verbose>1){print("running RF")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #GBM----
          if("GBM"%in%base.models){
            start<-Sys.time()
            X<-as.matrix(Genos[rownames(Phenos)[-test.folds[[i]]],])[!is.na(Phenos[-test.folds[[i]],t]),]
            Y<-Phenos[rownames(X),t]
            gbmdata<-cbind.data.frame(Y,X)
            gbm.fit<-gbm::gbm(Y~.,data = gbmdata,cv.folds = 5,distribution = "gaussian",n.cores = n.cores,
                         interaction.depth = 3,n.minobsinnode = 5,n.trees = gbm.ntrees,shrinkage = 0.025,verbose = F)
            doParallel::stopImplicitCluster()
            opt.tree<-gbm::gbm.perf(gbm.fit,plot.it = F,method = "cv")
            all.STmod.out.of.test.fold.preds[,"GBM"]<-gbm::predict.gbm(gbm.fit,n.trees = opt.tree,newdata = as.data.frame(Genos[rownames(all.STmod.out.of.test.fold.preds),]))
            all.STmod.in.test.fold.preds[,"GBM"]<-gbm::predict.gbm(gbm.fit,n.trees = opt.tree,newdata = as.data.frame(Genos[rownames(all.STmod.in.test.fold.preds),]))
            if(verbose>1){print("running GBM")}
            if(verbose>1){print(Sys.time()-start)}
          }


          #XGB----
          if("XGB"%in%base.models){
            start<-Sys.time()
            X<-as.matrix(Genos[rownames(Phenos)[-test.folds[[i]]],])[!is.na(Phenos[-test.folds[[i]],t]),]
            Y<-Phenos[rownames(X),t]
            gbmdata<-as.matrix(X)
            xgb.tune<-expand.grid(depth=c(2:5),shrink=c(0.05,0.1,0.2))
            xgb.tune<-xgb.tune[sample(1:nrow(xgb.tune))[1:(nrow(xgb.tune)*0.5)],]
            xgb.tune.error<-c()
            for (h in 1:nrow(xgb.tune)) {
              if(verbose>1){cat("|",sep="")}
              xgb.cv<-xgboost::xgb.cv(data = gbmdata,label = Y,max_depth = xgb.tune$depth[h], eta =xgb.tune$shrink[h], nrounds = 500,early_stopping_rounds = 5,
                             nfold=5, objective = "reg:squarederror",verbose=F)
              xgb.tune.error[h]<-min(xgb.cv$evaluation_log$test_rmse_mean)
            }
            xgb.cv<-xgboost::xgb.cv(data = gbmdata,label = Y,max_depth = xgb.tune$depth[which.min(xgb.tune.error)], eta =xgb.tune$shrink[which.min(xgb.tune.error)], nrounds = 500,early_stopping_rounds = 5,
                           nfold=5,nthread = n.cores, objective = "reg:squarederror",verbose=F)
            xgb.model<-xgboost::xgboost(data = gbmdata,label = Y,max_depth = xgb.tune$depth[which.min(xgb.tune.error)],
                               eta = xgb.tune$shrink[which.min(xgb.tune.error)], nrounds = xgb.cv$best_iteration,
                               objective = "reg:squarederror",verbose=F)
            all.STmod.out.of.test.fold.preds[,"XGB"]<-predict(xgb.model,newdata = Genos[rownames(all.STmod.out.of.test.fold.preds),])
            all.STmod.in.test.fold.preds[,"XGB"]<-predict(xgb.model,newdata = Genos[rownames(all.STmod.in.test.fold.preds),])
            if(verbose>1){print("running XGB")}
            if(verbose>1){print(Sys.time()-start)}
          }

          #LASSO----
          if("LASSO"%in%base.models){
            start<-Sys.time()
            X<-as.matrix(Genos[rownames(Phenos)[-test.folds[[i]]],])[!is.na(Phenos[-test.folds[[i]],t]),]
            cv.fit<-glmnet::cv.glmnet(x=X,y=Phenos[rownames(X),t],nfolds = 8)
            lasso.fit<-glmnet::glmnet(x=X,y=Phenos[rownames(X),t],lambda = cv.fit$lambda.min,alpha = 1)
            all.STmod.out.of.test.fold.preds[,"LASSO"]<-predict(lasso.fit,newx = as.matrix(Genos[rownames(all.STmod.out.of.test.fold.preds),]))
            all.STmod.in.test.fold.preds[,"LASSO"]<-predict(lasso.fit,newx = as.matrix(Genos[rownames(all.STmod.in.test.fold.preds),]))
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
            fm<-BGLR::BGLR(y =Phenos.na[,t],ETA = ETA,nIter=5000,burnIn=1000,df0=5,S0=2,verbose = F)
            fm.preds<-fm$yHat
            names(fm.preds)<-rownames(D)
            all.STmod.out.of.test.fold.preds[,"RKHS"]<-fm.preds[rownames(all.STmod.out.of.test.fold.preds)]
            all.STmod.in.test.fold.preds[,"RKHS"]<-fm.preds[rownames(all.STmod.in.test.fold.preds)]
            if(verbose>1){print("running RKHS")}
            if(verbose>1){print(Sys.time()-start)}
          }
        }


        #Make Gmod ensemble predictions----
        all.g.mod.ens.out.of.test.fold.preds<-matrix(NA,nrow = nrow(Phenos[test.folds[[i]],]),ncol=4,
                                                     dimnames = list(rownames(Phenos[test.folds[[i]],]),c("RF Ens","GBM Ens","LASSO Ens","PCR Ens")))
        all.g.mod.ens.out.of.test.fold.preds[,"RF Ens"]<-predict(rf.Gmod.ens.fit,newdata = all.STmod.out.of.test.fold.preds)
        all.g.mod.ens.out.of.test.fold.preds[,"GBM Ens"]<-gbm::predict.gbm(GBM.Gmod.ens.fit,newdata = as.data.frame(all.STmod.out.of.test.fold.preds),
                                                                      n.trees = GBM.opt.tree)
        all.g.mod.ens.out.of.test.fold.preds[,"LASSO Ens"]<-predict(lasso.fit.ens,newx = as.matrix(all.STmod.out.of.test.fold.preds))
        all.g.mod.ens.out.of.test.fold.preds[,"PCR Ens"]<-as.matrix(all.STmod.out.of.test.fold.preds)%*%pcr.beta.X

        all.trait.out.of.test.fold.preds[[t]]<-all.STmod.out.of.test.fold.preds
        all.trait.in.test.fold.preds[[t]]<-all.STmod.out.of.valid.fold.preds
        all.trait.Gmod.ens[[t]]<-all.g.mod.ens.out.of.test.fold.preds
        all.trait.out.of.test.fold.Gmodens.preds[,t]<-all.g.mod.ens.out.of.test.fold.preds[,"GBM Ens"]
        all.trait.in.test.fold.Gmodens.preds[,t]<-gbm::predict.gbm(GBM.Gmod.ens.fit,newdata = as.data.frame(all.STmod.in.test.fold.preds),
                                                              n.trees = GBM.opt.tree)


      } #End of traits loop


      #MT model ensembles----
      all.trait.MT.out.of.test.fold.preds<-list()
      if(verbose>0){print(paste("=============MT ensemble training for test fold",i," CV round",r,"================"))}
      for (t in 1:length(traits)) {
        if(verbose>1){print(traits[t])}
        MT.Ens<-matrix(NA,nrow = nrow(Phenos[test.folds[[i]],]),ncol=3,
                       dimnames = list(rownames(Phenos[test.folds[[i]],]),c("RF MT Ens","GBM MT Ens","LASSO MT Ens")))

        X<-matrix(unlist(all.trait.in.test.fold.preds),
                  nrow = nrow(all.trait.in.test.fold.preds[[1]]),
                  ncol=ncol(all.trait.in.test.fold.preds[[1]])*length(all.trait.in.test.fold.preds),
                  byrow = F,dimnames = list(rownames(all.trait.in.test.fold.preds[[1]]),
                                            paste(expand.grid(base.models,traits)[,1],
                                                  expand.grid(base.models,traits)[,2])[1:(length(all.trait.in.test.fold.preds)*length(base.models))]))

        newX<-matrix(unlist(all.trait.out.of.test.fold.preds),
                     nrow = nrow(all.trait.out.of.test.fold.preds[[1]]),
                     ncol=ncol(all.trait.out.of.test.fold.preds[[1]])*length(all.trait.out.of.test.fold.preds),
                     byrow = F,dimnames = list(rownames(all.trait.out.of.test.fold.preds[[1]]),
                                               paste(expand.grid(base.models,traits)[,1],
                                                     expand.grid(base.models,traits)[,2])[1:(length(all.trait.out.of.test.fold.preds)*length(base.models))]))
        Y<-Phenos[rownames(X),t]

        {#RF Ensemble----
          rf.fit<-randomForest::randomForest(x=X[!is.na(Y),],y=Y[!is.na(Y)],
                               ntree=300,maxnodes=100,importance=F)
          MT.Ens[,"RF MT Ens"]<-predict(rf.fit,newdata = newX)
        }

        {#GBM Ensemble----
          gbmdata<-cbind.data.frame(Y,X)

          gbm.tune<-expand.grid(depth=c(2:5),min.obs=c(2,4,6),shrink=c(0.1,0.05,0.01,0.005))
          gbm.tune<-gbm.tune[sample(1:nrow(gbm.tune))[1:(nrow(gbm.tune)*0.5)],]
          gbm.tune.error<-c()
          for (h in 1:nrow(gbm.tune)) {
            if(verbose>1){cat("|",sep="")}
            gbm.fit<-gbm::gbm(Y[!is.na(Y)]~.,data = gbmdata[!is.na(Y),],cv.folds = 4,
                         interaction.depth = gbm.tune$depth[h],
                         n.minobsinnode = gbm.tune$min.obs[h],
                         n.trees = 500,n.cores = n.cores,
                         shrinkage =gbm.tune$shrink[h],verbose = F,distribution = "gaussian")
            doParallel::stopImplicitCluster()
            gbm.tune.error[h]<-min(gbm.fit$cv.error)
          }
          if(verbose>1){cat(":)",sep="")}
          if(verbose>1){print(gbm.tune[which.min(gbm.tune.error),])}
          gbm.fit<-gbm::gbm(Y[!is.na(Y)]~.,data = gbmdata[!is.na(Y),],cv.folds = 8,
                       interaction.depth = gbm.tune$depth[which.min(gbm.tune.error)],
                       n.minobsinnode = gbm.tune$min.obs[which.min(gbm.tune.error)],
                       n.trees = 500,n.cores = n.cores,
                       shrinkage =gbm.tune$shrink[which.min(gbm.tune.error)],verbose = F,distribution = "gaussian")
          doParallel::stopImplicitCluster()
          opt.tree<-gbm.perf(gbm.fit,plot.it = F,method = "cv")

          MT.Ens[,"GBM MT Ens"]<-gbm::predict.gbm(gbm.fit,newdata = as.data.frame(newX),n.trees = opt.tree)
        }

        {#LASSO Ensemble----
          cv.fit.ens<-glmnet::cv.glmnet(x=X[!is.na(Y),],y=Y[!is.na(Y)],nfolds = 8)
          lasso.fit.ens<-glmnet::glmnet(x=X[!is.na(Y),],y=Y[!is.na(Y)],lambda = cv.fit.ens$lambda.min,alpha = 1)
          MT.Ens[,"LASSO MT Ens"]<-predict(lasso.fit.ens,newx = newX)
        }
        all.trait.MT.out.of.test.fold.preds[[t]]<-MT.Ens
      }

      fold.preds<-list("all.ST.Gmod.preds"=all.trait.out.of.test.fold.preds,
                       "all.ST.Ens.preds"=all.trait.Gmod.ens,
                       "all.MT.ens.preds"=all.trait.MT.out.of.test.fold.preds)
      all.out.of.fold.preds[[i]]<-fold.preds
      if(verbose>1){print(paste("===========FINISHED test fold",i," CV round",r,"==========="))}
    }

    all.CV.rounds[[r]]<-all.out.of.fold.preds

    if(verbose>1){print(paste("===========FINISHED CV round",r,"==========="))}
  }

  fold.names<-lapply(all.CV.rounds,function(r)
    c(lapply(r,function(i) rownames(i$all.ST.Ens.preds[[1]]))))
  names(fold.names)<-paste("Round",1:length(fold.names))
  for (r in 1:length(fold.names)) {
    names(fold.names[[r]])<-paste("CV fold",1:n.test.folds)
  }

  base.preds<-array(NA,dim=c(length(rownames(Phenos)),length(base.models),length(traits),length(all.CV.rounds)))
  for(r in 1:length(all.CV.rounds)){
    for(t in 1:length(traits)){
      trait.mat<-as.matrix(rbindlist(lapply(all.CV.rounds[[r]],function(x) data.frame(x$all.ST.Gmod.preds[[t]]))))
      rownames(trait.mat)<-unlist(fold.names[[r]])
      base.preds[,,t,r]<- trait.mat[rownames(Phenos),]
      dimnames(base.preds)<-list(rownames(Phenos),base.models,traits,1:length(all.CV.rounds))
    }
  }

  gmod.ens.preds<-array(NA,dim=c(length(rownames(Phenos)),ncol(all.CV.rounds[[1]][[1]]$all.ST.Ens.preds[[1]]),length(traits),length(all.CV.rounds)))
  for(r in 1:length(all.CV.rounds)){
    for(t in 1:length(traits)){
      trait.mat<-as.matrix(rbindlist(lapply(all.CV.rounds[[r]],function(x) data.frame(x$all.ST.Ens.preds[[t]]))))
      rownames(trait.mat)<-unlist(fold.names[[r]])
      gmod.ens.preds[,,t,r]<- trait.mat[rownames(Phenos),]
      dimnames(gmod.ens.preds)<-list(rownames(Phenos),colnames(all.CV.rounds[[1]][[1]]$all.ST.Ens.preds[[1]]),traits,1:length(all.CV.rounds))
    }
  }

  MT.ens.preds<-array(NA,dim=c(length(rownames(Phenos)),ncol(all.CV.rounds[[1]][[1]]$all.MT.ens.preds[[1]]),length(traits),length(all.CV.rounds)))
  for(r in 1:length(all.CV.rounds)){
    for(t in 1:length(traits)){
      trait.mat<-as.matrix(rbindlist(lapply(all.CV.rounds[[r]],function(x) data.frame(x$all.MT.ens.preds[[t]]))))
      rownames(trait.mat)<-unlist(fold.names[[r]])
      MT.ens.preds[,,t,r]<- trait.mat[rownames(Phenos),]
      dimnames(MT.ens.preds)<-list(rownames(Phenos),colnames(all.CV.rounds[[1]][[1]]$all.MT.ens.preds[[1]]),traits,1:length(all.CV.rounds))
    }
  }

  all.out.of.fold.preds<-list("Base models"=base.preds,"Gmod ensemble models"=gmod.ens.preds,"MT ensemble models"=MT.ens.preds)

  all.traits.pred.r<-list()
  for(t in 1:length(traits)){
    base.rs<-matrix(NA,nrow=length(all.CV.rounds),ncol=length(base.models),
                    dimnames = list(1:length(all.CV.rounds),base.models))
    Gmod.ens.rs<-matrix(NA,nrow=length(all.CV.rounds),ncol=ncol(all.CV.rounds[[1]][[1]]$all.ST.Ens.preds[[1]]),
                        dimnames = list(1:length(all.CV.rounds),colnames(all.CV.rounds[[1]][[1]]$all.ST.Ens.preds[[1]])))
    MT.ens.rs<-matrix(NA,nrow=length(all.CV.rounds),ncol=ncol(all.CV.rounds[[1]][[1]]$all.MT.ens.preds[[1]]),
                      dimnames = list(1:length(all.CV.rounds),colnames(all.CV.rounds[[1]][[1]]$all.MT.ens.preds[[1]])))

    for(i in 1:ncol(base.rs)){
      for(r in 1:length(all.CV.rounds)) {
        base.rs[r,i]<-cor(Phenos[,t],all.out.of.fold.preds$`Base models`[,i,t,r], use="pairwise.complete.obs")
      }
    }

    for(i in 1:ncol(all.CV.rounds[[1]][[1]]$all.ST.Ens.preds[[1]])){
      for(r in 1:length(all.CV.rounds)) {
        Gmod.ens.rs[r,i]<-cor(Phenos[,t],all.out.of.fold.preds$`Gmod ensemble models`[,i,t,r], use="pairwise.complete.obs")
      }
    }

    for(i in 1:ncol(all.CV.rounds[[1]][[1]]$all.MT.ens.preds[[1]])){
      for(r in 1:length(all.CV.rounds)) {
        MT.ens.rs[r,i]<-cor(Phenos[,t],all.out.of.fold.preds$`MT ensemble models`[,i,t,r], use="pairwise.complete.obs")
      }
    }

    all.traits.pred.r[[t]]<-list("Base models"=base.rs,"Gmod ensemble models"=Gmod.ens.rs,"MT ensemble models"=MT.ens.rs)
  }
  names(all.traits.pred.r)<-traits

  all.out<-list("All out-of-fold predictions"=all.out.of.fold.preds,
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
#' @param n.valid.folds (integer) number of validation folds to run within each test fold. Default if 8.
#' @param CV.groups (list) of vectors containing Genotype ID names. Cross validation can also be performed across predefined subsets of genotypes.
#' @param n.cores Number of cores that can be used to run cross validation of RF, GBM and XGB models. The Default detectcores() from Parallel will use available multiple
#' cores on a Windows OS but n.cores should be defined when running on a HPC or other OS.
#' @param rf.ntrees (integer) Number of trees to use for each random forest model
#' @param gbm.ntrees (integer) Max number of trees to use for each GBM model
#' @param verbose Level of progress to print throughout. 0 = none; 1 = some; 2 = all.
#'
#' @return Returns a list object of:
#'
#'-\strong{All out-of-fold predictions} for all genotypes for base prediction models, genetic model ensemble and multi-trait ensemble models for each trait.
#'
#'-\strong{Phenotype data} all observed multi-trait input data.
#'
#'-\strong{Prediction accuracies} list of correlation coeficient (r) for each model for each trait.
#'
#'-\strong{Base models} \emph{vector} of base prediction models used.
#'
#'-\strong{Traits} \emph{vector} of trait names used.
#'
#'-\strong{CV fold allocations} A two level list of length = n.fold.rounds, each containing the genotype ID names allocated to each test fold.
#'
#'#' @details Ensemble predictions for multiple genetic models and multiple traits are trained on out-of-validation fold predictions using several model approaches.
#' Including a range of genetic models that capture a diversity of genetic architectures may improve accuracy of ensembles. Unbalanced phenotyping across
#' genotypes may increase accuracy of multi-trait ensembles when several related but partially phenotyped traits are included.
#'
#' @export

Train.Predict.MTens<-function(Genos.train,New.Genos,Phenos,
                              n.valid.folds=8,
                              rf.ntrees=300,
                              gbm.ntrees=200,
                              n.cores=detectCores()-1,
                              base.models=c("rrblup","EGBLUP","RKHS","LASSO","RF","GBM","XGB")
){
  run.start.time<-Sys.time()
  if(length(rownames(Genos.train)[!rownames(Genos.train)%in%rownames(Phenos)])>0){
    print("Lines in Genos but not in Phenos:")
    print(rownames(Genos.train)[!rownames(Genos.train)%in%rownames(Phenos)])}
  if(length(rownames(Phenos)[!rownames(Phenos)%in%rownames(Genos.train)])>0){
    print("Lines in Phenos but not in Genos:")
    print(rownames(Phenos)[!rownames(Phenos)%in%rownames(Genos.train)])}

  Genos.train<-as.matrix(Genos.train[rownames(Genos.train)%in%rownames(Phenos),])
  Phenos<-as.data.frame(Phenos[rownames(Phenos)%in%rownames(Genos.train),])

  traits<-colnames(Phenos)
  combined.Genos<-rbind(Genos.train,New.Genos)

  valid.folds<-split(sample(c(1:nrow(Phenos))),1:n.valid.folds)#make folds


  all.trait.base.model.preds<-list()
  all.trait.out.of.valid.base.model.preds<-list()
  all.trait.Gmod.ens.preds<-list()
  for (t in 1:length(traits)) {
    print(paste("Starting",traits[t]))
    #CV each ST model for ensemble base.models----
    all.STmod.out.of.valid.fold.preds<-matrix(NA,nrow = nrow(Phenos),ncol=length(base.models),
                                              dimnames = list(rownames(Phenos),base.models))
    for (v in 1:length(valid.folds)) {
      print(paste("========Starting",traits[t],"validation fold",v,"=========="))
      valid.names<-rownames(Phenos)[valid.folds[[v]]]
      #RRBLUP----
      if("rrblup"%in%base.models){
        start<-Sys.time()
        Phenos.na<-Phenos
        Phenos.na[valid.folds[[v]],]<-NA
        Phenos.na<-cbind.data.frame(rownames(Phenos.na),Phenos.na)
        colnames(Phenos.na)[1]<-"line"

        A <- A.mat(Genos.train)
        A<-A[rownames(Phenos),rownames(Phenos)]
        gblup.model<-kin.blup(data=Phenos.na,pheno = traits[t],geno = "line",K = A)
        all.STmod.out.of.valid.fold.preds[valid.names,"rrblup"]<-gblup.model$pred[rownames(all.STmod.out.of.valid.fold.preds)][valid.names]
        print("running rrBLUP")
        print(Sys.time()-start)
      }

      #EGBLUP----
      if("EGBLUP"%in%base.models){
        start<-Sys.time()
        Phenos.na<-Phenos
        Phenos.na<-Phenos.na[!is.na(Phenos.na[,t])|rownames(Phenos.na)%in%valid.names,]#Remove lines with na but not in validation fold
        Phenos.na[valid.names,]<-NA
        A <- A.mat(Genos.train[rownames(Phenos.na),])
        H<-hadamard.prod(A,A)
        Phenos.na<-Phenos.na[rownames(A),]
        ETA<-list(list(K=A,model="RKHS"),list(K=H,model="RKHS"))
        fm<-BGLR(y =Phenos.na[,t],ETA = ETA,nIter=3000,burnIn=500,df0=5,S0=2,verbose = F)
        fm.preds<-fm$yHat
        names(fm.preds)<-rownames(A)
        all.STmod.out.of.valid.fold.preds[valid.names,"EGBLUP"]<-fm.preds[rownames(all.STmod.out.of.valid.fold.preds)][valid.names]
        print("running EGBLUP")
        print(Sys.time()-start)
      }

      #RF----
      if("RF"%in%base.models){
        start<-Sys.time()
        X<-as.matrix(Genos.train[rownames(Phenos),])[!is.na(Phenos[,t]),]
        X<-X[!rownames(X)%in%valid.names,]
        cl <- makeCluster(n.cores)
        registerDoParallel(cl)
        rf.fit <- foreach(ntree=rep(round(rf.ntrees/n.cores),n.cores),.combine=combine,.multicombine=TRUE,.packages='randomForest') %dopar%{
          randomForest(x=X,y=Phenos[rownames(X),t],ntree=ntree,importance=F)
        }
        stopCluster(cl)
        all.STmod.out.of.valid.fold.preds[valid.names,"RF"]<-predict(rf.fit,newdata = as.matrix(Genos.train[valid.names,]))
        print("running RF")
        print(Sys.time()-start)
      }

      #GBM----
      if("GBM"%in%base.models){
        start<-Sys.time()
        X<-as.matrix(Genos.train[rownames(Phenos),])[!is.na(Phenos[,t]),]
        X<-X[!rownames(X)%in%valid.names,]
        Y<-Phenos[rownames(X),t]
        gbmdata<-cbind.data.frame(Y,X)
        gbm.fit<-gbm(Y~.,data = gbmdata,cv.folds = 5,distribution = "gaussian",interaction.depth = 3,n.minobsinnode = 5,
                     n.trees = gbm.ntrees,shrinkage = 0.1,verbose = F,n.cores = n.cores)
        opt.tree<-gbm.perf(gbm.fit,plot.it = F,method = "cv")
        all.STmod.out.of.valid.fold.preds[valid.names,"GBM"]<-predict.gbm(gbm.fit,n.trees = opt.tree,newdata = as.data.frame(Genos.train[valid.names,]))
        print("running GBM")
        print(Sys.time()-start)
      }

      #XGB----
      if("XGB"%in%base.models){
        start<-Sys.time()
        X<-as.matrix(Genos.train[rownames(Phenos),])[!is.na(Phenos[,t]),]
        X<-X[!rownames(X)%in%valid.names,]
        Y<-Phenos[rownames(X),t]
        gbmdata<-as.matrix(X)
        xgb.tune<-expand.grid(depth=c(2:5),shrink=c(0.05,0.1,0.2))
        xgb.tune<-xgb.tune[sample(1:nrow(xgb.tune))[1:(nrow(xgb.tune)*0.7)],]
        xgb.tune.error<-c()
        for (h in 1:nrow(xgb.tune)) {
          cat("|",sep="")
          xgb.cv<-xgb.cv(data = gbmdata,label = Y,max_depth = xgb.tune$depth[h], eta =xgb.tune$shrink[h], nrounds = 500,early_stopping_rounds = 5,
                         nfold=5,nthread = n.cores, objective = "reg:squarederror",verbose=F)
          xgb.tune.error[h]<-min(xgb.cv$evaluation_log$test_rmse_mean)
          cat(xgb.cv$best_iteration,min(xgb.cv$evaluation_log$test_rmse_mean))
        }
        xgb.cv<-xgb.cv(data = gbmdata,label = Y,max_depth = xgb.tune$depth[which.min(xgb.tune.error)], eta =xgb.tune$shrink[which.min(xgb.tune.error)],
                       nrounds = 500,early_stopping_rounds = 5,nfold=5,nthread = n.cores, objective = "reg:squarederror",verbose=F)
        xgb.model<-xgboost(data = gbmdata,label = Y,max_depth = xgb.tune$depth[which.min(xgb.tune.error)],
                           eta = xgb.tune$shrink[which.min(xgb.tune.error)], nrounds = xgb.cv$best_iteration,
                           nthread = n.cores, objective = "reg:squarederror",verbose=F)
        all.STmod.out.of.valid.fold.preds[valid.names,"XGB"]<-predict(xgb.model,newdata = as.matrix(Genos.train[valid.names,]))
        print("running XGB")
        print(Sys.time()-start)
      }


      #LASSO----
      if("LASSO"%in%base.models){
        start<-Sys.time()
        X<-as.matrix(Genos.train[rownames(Phenos),])[!is.na(Phenos[,t]),]
        X<-X[!rownames(X)%in%valid.names,]
        cv.fit<-cv.glmnet(x=X,y=Phenos[rownames(X),t],nfolds = 5)
        lasso.fit<-glmnet(x=X,y=Phenos[rownames(X),t],lambda = cv.fit$lambda.min,alpha = 1)
        all.STmod.out.of.valid.fold.preds[valid.names,"LASSO"]<-predict(lasso.fit,newx = as.matrix(Genos.train[valid.names,]))
        print("running LASSO")
        print(Sys.time()-start)
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
        fm<-BGLR(y =Phenos.na[,t],ETA = ETA,nIter=3000,burnIn=500,df0=5,S0=2,verbose = F)
        fm.preds<-fm$yHat
        names(fm.preds)<-rownames(D)
        all.STmod.out.of.valid.fold.preds[valid.names,"RKHS"]<-fm.preds[rownames(all.STmod.out.of.valid.fold.preds)][valid.names]
        print("running RKHS")
        print(Sys.time()-start)
      }
    }

    all.trait.out.of.valid.base.model.preds[[t]]<-all.STmod.out.of.valid.fold.preds

    print("Out-of-validation fold prediction accuracies:")
    print(apply(all.STmod.out.of.valid.fold.preds,2,function(x) cor(x,Phenos[,t],use="pairwise.complete.obs")))

    #train Gmod ensemble model----
    print("Fitting Gmod ensemble base.models")
    #RF.Ensemble---
    Y<-Phenos[,t]
    rf.Gmod.ens.fit<-randomForest(x=all.STmod.out.of.valid.fold.preds[!is.na(Y),],y=Y[!is.na(Y)],
                                  ntree=500,importance=F)
    #GBM.Ensemble----
    gbmdata<-cbind.data.frame(Y,all.STmod.out.of.valid.fold.preds)[!is.na(Y),]

    gbm.tune<-expand.grid(depth=c(2:5),min.obs=c(2,4,6),shrink=c(0.1,0.05,0.01,0.005))
    gbm.tune<-gbm.tune[sample(1:nrow(gbm.tune))[1:(nrow(gbm.tune)*0.5)],]
    gbm.tune.error<-c()
    for (h in 1:nrow(gbm.tune)) {
      cat("|",sep="")
      gbm.fit<-gbm(gbmdata$Y~.,data = gbmdata,cv.folds = 4,
                   interaction.depth = gbm.tune$depth[h],
                   n.minobsinnode = gbm.tune$min.obs[h],
                   n.trees = 500,n.cores = n.cores,
                   shrinkage =gbm.tune$shrink[h],verbose = F,distribution = "gaussian")
      stopImplicitCluster()
      gbm.tune.error[h]<-min(gbm.fit$cv.error)
      cat(which.min(gbm.fit$cv.error))
    }
    print(":)")
    print(gbm.tune[which.min(gbm.tune.error),])
    GBM.Gmod.ens.fit<-gbm(gbmdata$Y~.,data = gbmdata,cv.folds = 8,
                          interaction.depth = gbm.tune$depth[which.min(gbm.tune.error)],
                          n.minobsinnode = gbm.tune$min.obs[which.min(gbm.tune.error)],
                          n.trees = 500,n.cores = n.cores,
                          shrinkage =gbm.tune$shrink[which.min(gbm.tune.error)],verbose = F,distribution = "gaussian")
    stopImplicitCluster()
    GBM.opt.tree<-gbm.perf(GBM.Gmod.ens.fit,plot.it = F,method = "cv")

    #LASSO Ensemble----
    cv.fit.ens<-cv.glmnet(x=all.STmod.out.of.valid.fold.preds[!is.na(Y),],y=Y[!is.na(Y)],nfolds = 8)
    lasso.fit.ens<-glmnet(x=all.STmod.out.of.valid.fold.preds[!is.na(Y),],y=Y[!is.na(Y)],lambda = cv.fit.ens$lambda.min,alpha = 1)

    #PCR Ensemble----
    PCA<-prcomp(all.STmod.out.of.valid.fold.preds)
    pcr.lmod<-lm(Phenos[,t]~.,data=as.data.frame(PCA$x))
    beta.Z <- as.matrix(pcr.lmod$coefficients[-1])
    V <- as.matrix(PCA$rotation)
    pcr.beta.X <- V %*% beta.Z


    #Run each ST Gmodel on full training data----
    all.new.data.base.preds<-matrix(NA,nrow = nrow(New.Genos),ncol=length(base.models),
                                    dimnames = list(rownames(New.Genos),base.models))

    {
      print(paste("========Starting",traits[t],"====Full training set fitting=========="))

      #RRBLUP----
      if("rrblup"%in%base.models){
        start<-Sys.time()
        Phenos.na<-Phenos
        Phenos.na<-cbind.data.frame(rownames(Phenos.na),Phenos.na)
        colnames(Phenos.na)[1]<-"line"
        A <- A.mat(combined.Genos)
        Phenos.na<-Phenos.na[rownames(A),]
        gblup.model<-kin.blup(data=Phenos.na,pheno = traits[t],geno = "line",K = A)
        all.new.data.base.preds[,"rrblup"]<-gblup.model$pred[rownames(all.new.data.base.preds)]
        print("running rrBLUP")
        print(Sys.time()-start)
      }

      #EGBLUP----
      if("EGBLUP"%in%base.models){
        start<-Sys.time()
        Phenos.na<-Phenos[!is.na(Phenos[,t]),]
        A <- A.mat(combined.Genos[unique(c(rownames(Phenos.na),rownames(New.Genos))),])
        H<-hadamard.prod(A,A)
        Phenos.na<-Phenos.na[rownames(A),]
        ETA<-list(list(K=A,model="RKHS"),list(K=H,model="RKHS"))
        fm<-BGLR(y =Phenos.na[,t],ETA = ETA,nIter=5000,burnIn=1000,df0=5,S0=2,verbose = F)
        fm.preds<-fm$yHat
        names(fm.preds)<-rownames(A)
        all.new.data.base.preds[,"EGBLUP"]<-fm.preds[rownames(all.new.data.base.preds)]
        print("running EGBLUP")
        print(Sys.time()-start)
      }

      #RF----
      if("RF"%in%base.models){
        start<-Sys.time()
        X<-as.matrix(Genos.train[rownames(Phenos),])[!is.na(Phenos[,t]),]
        cl <- makeCluster(n.cores)
        registerDoParallel(cl)
        rf.fit <- foreach(ntree=rep(round(rf.ntrees/n.cores),n.cores),.combine=combine,.multicombine=TRUE,.packages='randomForest') %dopar%{
          randomForest(x=X,y=Phenos[rownames(X),t],ntree=ntree,importance=F)
        }
        stopCluster(cl)
        all.new.data.base.preds[,"RF"]<-predict(rf.fit,newdata = as.matrix(New.Genos[rownames(all.new.data.base.preds),]))
        print("running RF")
        print(Sys.time()-start)
      }

      #GBM----
      if("GBM"%in%base.models){
        start<-Sys.time()
        X<-as.matrix(Genos.train[rownames(Phenos),])[!is.na(Phenos[,t]),]
        Y<-Phenos[rownames(X),t]
        gbmdata<-cbind.data.frame(Y,X)
        gbm.fit<-gbm(Y~.,data = gbmdata,cv.folds = 5,distribution = "gaussian",
                     interaction.depth = 3,n.minobsinnode = 5,n.trees = gbm.ntrees,shrinkage = 0.025,verbose = F,n.cores = n.cores)
        opt.tree<-gbm.perf(gbm.fit,plot.it = F,method = "cv")

        all.new.data.base.preds[,"GBM"]<-predict.gbm(gbm.fit,n.trees = opt.tree,newdata = as.data.frame(New.Genos[rownames(all.new.data.base.preds),]))
        print("running GBM")
        print(Sys.time()-start)
      }


      #XGB----
      if("XGB"%in%base.models){
        start<-Sys.time()
        X<-as.matrix(Genos.train[rownames(Phenos),])[!is.na(Phenos[,t]),]
        Y<-Phenos[rownames(X),t]
        gbmdata<-as.matrix(X)
        xgb.tune<-expand.grid(depth=c(2:5),shrink=c(0.05,0.1,0.2))
        xgb.tune<-xgb.tune[sample(1:nrow(xgb.tune))[1:(nrow(xgb.tune)*0.7)],]
        xgb.tune.error<-c()
        for (h in 1:nrow(xgb.tune)) {
          cat("|",sep="")
          xgb.cv<-xgb.cv(data = gbmdata,label = Y,max_depth = xgb.tune$depth[h], eta =xgb.tune$shrink[h], nrounds = 500,early_stopping_rounds = 5,
                         nfold=5,nthread = n.cores, objective = "reg:squarederror",verbose=F)
          xgb.tune.error[h]<-min(xgb.cv$evaluation_log$test_rmse_mean)
        }
        xgb.cv<-xgb.cv(data = gbmdata,label = Y,max_depth = xgb.tune$depth[which.min(xgb.tune.error)], eta =xgb.tune$shrink[[which.min(xgb.tune.error)]], nrounds = 500,early_stopping_rounds = 5,
                       nfold=5,nthread = n.cores, objective = "reg:squarederror",verbose=F)
        xgb.model<-xgboost(data = gbmdata,label = Y,max_depth = xgb.tune$depth[which.min(xgb.tune.error)],
                           eta = xgb.tune$shrink[which.min(xgb.tune.error)], nrounds = xgb.cv$best_iteration,
                           nthread = n.cores, objective = "reg:squarederror",verbose=F)
        all.new.data.base.preds[,"XGB"]<-predict(xgb.model,newdata = as.matrix(New.Genos[rownames(all.new.data.base.preds),]))
        print("running XGB")
        print(Sys.time()-start)
      }

      #LASSO----
      if("LASSO"%in%base.models){
        start<-Sys.time()
        X<-as.matrix(Genos.train[rownames(Phenos),])[!is.na(Phenos[,t]),]
        cv.fit<-cv.glmnet(x=X,y=Phenos[rownames(X),t],nfolds = 8)
        lasso.fit<-glmnet(x=X,y=Phenos[rownames(X),t],lambda = cv.fit$lambda.min,alpha = 1)
        all.new.data.base.preds[,"LASSO"]<-predict(lasso.fit,newx = as.matrix(New.Genos[rownames(all.new.data.base.preds),]))
        print("running LASSO")
        print(Sys.time()-start)
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
        fm<-BGLR(y =Phenos.na[,t],ETA = ETA,nIter=5000,burnIn=1000,df0=5,S0=2,verbose = F)
        fm.preds<-fm$yHat
        names(fm.preds)<-rownames(D)
        all.new.data.base.preds[,"RKHS"]<-fm.preds[rownames(all.new.data.base.preds)]
        print("running RKHS")
        print(Sys.time()-start)
      }
    }
    all.trait.base.model.preds[[t]]<-all.new.data.base.preds

    #Make Gmod ensemble predictions----
    all.new.data.gmod.ens.preds<-matrix(NA,nrow = nrow(New.Genos),ncol=4,
                                        dimnames = list(rownames(New.Genos),c("RF Ens","GBM Ens","LASSO Ens","PCR Ens")))
    all.new.data.gmod.ens.preds[,"RF Ens"]<-predict(rf.Gmod.ens.fit,newdata = all.new.data.base.preds)
    all.new.data.gmod.ens.preds[,"GBM Ens"]<-predict.gbm(GBM.Gmod.ens.fit,newdata = as.data.frame(all.new.data.base.preds),
                                                         n.trees = GBM.opt.tree)
    all.new.data.gmod.ens.preds[,"LASSO Ens"]<-predict(lasso.fit.ens,newx = as.matrix(all.new.data.base.preds))
    all.new.data.gmod.ens.preds[,"PCR Ens"]<-as.matrix(all.new.data.base.preds)%*%pcr.beta.X

    all.trait.Gmod.ens.preds[[t]]<-all.new.data.gmod.ens.preds

  } #End of traits loop


  #MT model ensembles----
  all.trait.MT.new.data.preds<-list()
  print(paste("=============MT ensemble training=============="))
  for (t in 1:length(traits)) {
    print(traits[t])
    MT.Ens<-matrix(NA,nrow = nrow(New.Genos),ncol=3,
                   dimnames = list(rownames(New.Genos),c("RF MT Ens","GBM MT Ens","LASSO MT Ens")))

    X<-matrix(unlist(all.trait.out.of.valid.base.model.preds),
              nrow = nrow(all.trait.out.of.valid.base.model.preds[[1]]),
              ncol=ncol(all.trait.out.of.valid.base.model.preds[[1]])*length(all.trait.out.of.valid.base.model.preds),
              byrow = F,dimnames = list(rownames(all.trait.out.of.valid.base.model.preds[[1]]),
                                        paste(expand.grid(base.models,traits)[,1],
                                              expand.grid(base.models,traits)[,2])[1:(length(all.trait.out.of.valid.base.model.preds)*length(base.models))]))

    newX<-matrix(unlist(all.trait.base.model.preds),
                 nrow = nrow(all.trait.base.model.preds[[1]]),
                 ncol=ncol(all.trait.base.model.preds[[1]])*length(all.trait.base.model.preds),
                 byrow = F,dimnames = list(rownames(all.trait.base.model.preds[[1]]),
                                           paste(expand.grid(base.models,traits)[,1],
                                                 expand.grid(base.models,traits)[,2])[1:(length(all.trait.base.model.preds)*length(base.models))]))
    Y<-Phenos[rownames(X),t]

    {#RF Ensemble----
      rf.fit<-randomForest(x=X[!is.na(Y),],y=Y[!is.na(Y)],
                           ntree=300,maxnodes=100,importance=F)
      MT.Ens[,"RF MT Ens"]<-predict(rf.fit,newdata = newX)
    }

    {#GBM Ensemble----
      gbmdata<-cbind.data.frame(Y,X)

      gbm.tune<-expand.grid(depth=c(2:5),min.obs=c(2,4,6),shrink=c(0.1,0.05,0.01,0.005))
      gbm.tune<-gbm.tune[sample(1:nrow(gbm.tune))[1:(nrow(gbm.tune)*0.5)],]
      gbm.tune.error<-c()
      for (h in 1:nrow(gbm.tune)) {
        cat("|",sep="")
        gbm.fit<-gbm(Y[!is.na(Y)]~.,data = gbmdata[!is.na(Y),],cv.folds = 4,
                     interaction.depth = gbm.tune$depth[h],
                     n.minobsinnode = gbm.tune$min.obs[h],
                     n.trees = 500,n.cores = n.cores,
                     shrinkage =gbm.tune$shrink[h],verbose = F,distribution = "gaussian")
        stopImplicitCluster()
        gbm.tune.error[h]<-min(gbm.fit$cv.error)
      }
      cat(":)",sep="")
      print(gbm.tune[which.min(gbm.tune.error),])
      gbm.fit<-gbm(Y[!is.na(Y)]~.,data = gbmdata[!is.na(Y),],cv.folds = 8,
                   interaction.depth = gbm.tune$depth[which.min(gbm.tune.error)],
                   n.minobsinnode = gbm.tune$min.obs[which.min(gbm.tune.error)],
                   n.trees = 500,n.cores = n.cores,
                   shrinkage =gbm.tune$shrink[which.min(gbm.tune.error)],verbose = F,distribution = "gaussian")
      stopImplicitCluster()
      opt.tree<-gbm.perf(gbm.fit,plot.it = F,method = "cv")

      MT.Ens[,"GBM MT Ens"]<-predict.gbm(gbm.fit,newdata = as.data.frame(newX),n.trees = opt.tree)
    }

    {#LASSO Ensemble----
      cv.fit.ens<-cv.glmnet(x=X[!is.na(Y),],y=Y[!is.na(Y)],nfolds = 8)
      lasso.fit.ens<-glmnet(x=X[!is.na(Y),],y=Y[!is.na(Y)],lambda = cv.fit.ens$lambda.min,alpha = 1)
      MT.Ens[,"LASSO MT Ens"]<-predict(lasso.fit.ens,newx = newX)
    }
    all.trait.MT.new.data.preds[[t]]<-MT.Ens
  }
  names(all.trait.MT.new.data.preds)<-traits
  names(all.trait.base.model.preds)<-traits
  names(all.trait.Gmod.ens.preds)<-traits

  out<-list("Base model preds"=all.trait.base.model.preds,
            "ST Ensemble preds"=all.trait.Gmod.ens.preds,
            "MT Ensemble preds"=all.trait.MT.new.data.preds,
            "Traits"=traits,
            "Base models"=base.models)

  print(paste("===========FINISHED========================"))
  print(Sys.time()-run.start.time)
  return(out)
}




















#' Plot multi-trait Principle component selection indices

plot.MT.PC.SI<-function(X,
                        SI.traits=NULL,
                        weights=NULL,
                        model="MT Ensemble preds",
                        sub.model="RF MT Ens",
                        col1=rgb(0.1,0.2,0.8,0.7,maxColorValue = 1),
                        col2=rgb(1,0,0,0.6,maxColorValue = 1)){
  preds<-X[[which(names(X)==model)]]
  preds<-matrix(unlist(lapply(preds,function(x) x[,colnames(x)==sub.model])),nrow = nrow(preds[[1]]),ncol = length(preds),
                dimnames = list(rownames(preds[[1]]),names(preds)))
  if(is.null(SI.traits)){
    SI.traits<-colnames(preds)
  }
  if(is.null(weights)){
    weights<-rep(1,length(SI.traits))
  }
  if (!length(SI.traits)==length(weights)) {
    print("SI.traits and weights vectors are not the same length...")
  }
  weights.mat<-matrix(rep(weights,nrow(preds)),nrow=nrow(preds),ncol = length(SI.traits),byrow = T)
  scaled<-scale(preds[,SI.traits])*weights.mat
  biplot(prcomp(scaled),cex=c(0.7,1),col =c(col1,col2))
}


