# MTensGS
## Multi-trait and multi-model ensemble Genomic Selection

### Installation
 To install from Github using devtools:
 
```
devtools::install_github("NickFlagleaf/MTensGS")
```
 
The package has an example dataset of wheat genotypes and phenotypes from:
```
library(MTensGS)
data(TGdata)
```
 
### Marker proccessing
The `marker.QC()` function can first be used to QC SNP marker data by imputing missing data, pruning based on LD and skimming based on a minor allele frequency and missing data threshold. The SNP markers should be coded as 0,1,2 for hom, het and hom calls respectively.

```
QCd.markers<-marker.QC(marker.mat = TGgenos)
```
With the default pruning threshold of 0.8, this should reduce the marker dataset from 2535 down to 867 markers.


### Cross validation of base and ensemble models
The `CV.MTens()` function can be then used to cross validate the prediction accuracy of a selection of possible base models and ensemble models. Although this will take a long time for large datasets, a quick example can be run with the QCd TG example dataset with three base models and just one round of cross validation:

```
MTens.CV.results<-CV.MTens(Genos = QCd.markers
                           ,Phenos = TGphenos
                           ,base.models = c("rrblup","LASSO","RF")
                           ,n.CV.rounds = 1,n.valid.rounds=1)
```

### Predicting into new data
Once accuracy of the ensembles have been confirmed through cross validation with existing data, the `Train.Predict.MTens()` function can be used to deploy the models with new marker data for genotypes that have not been phenotyped yet. In this example the first 300 lines are assumed to be the training set and the remaining 76 lines are the new lines to predict.

```
#Split the dataset in two
old.genos<-QCd.markers[1:300,]
new.genos<-QCd.markers[301:376,]
old.phens<-TGphenos[1:300,]

new_preds<-Train.Predict.MTens(Genos.train = old.genos,
                               New.Genos = new.genos,
                               Phenos = old.phens,
                               n.valid.rounds = 1,
                               base.models = c("rrblup","LASSO","RF"))
```

### Plotting and setting a multi-trait selection index
Once predictions for new lines have been made, the `MT.PCSI.plot()` function can be used to plot a PCA visualisation of the of the multi-trait predictions and set a selection index based on weightings for each trait. In this example protein (PROT) is given double the weighting in comparison to the other three yield traits.

```
sel.ind<-MT.PCSI.plot(new_preds,SI.traits = c("YLD_GBR_2011","YLD_GBR_2010","YLD_DEU_2010","PROT")
                      ,weights = c(1,1,1,2),cex1 = 0.8)
```
