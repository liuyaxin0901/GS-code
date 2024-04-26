rm(list=ls())
source("#File storage path")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

myGD=myGD
myGM=myGM
myCV=myCV #gender,age,pca
myY=myY #year

nrep=100
nfold=5
set.seed(980901)
WG.dif.y=matrix(NA,nrep,#Phenotypic column number 
	)
WMA.dif.y=matrix(NA,nrep,#Phenotypic column number
 )
#选择表型
for(k in 3:8){

  #去掉NA
  Y=myY[complete.cases(myY[,k]),c(1,k)]
  pred.Y=matrix(NA,nrow(Y),nrep)
  pred.Y.MA=matrix(NA,nrow(Y),nrep)
#-------------------------------------------------------------prediction
  gblupw=NULL 
  mablupw=NULL

  for(rep in 1:nrep){
  #running gBLUP
  sets=sample(cut(1:nrow(Y),nfold,labels=F),nrow(Y))  
  pred.Y.rep=NULL
  pred.Y.MA.rep=NULL

    for(i in 1:nfold){
      training=sets!=i 
        #------------------↓↓↓↓↓↓↓↓↓↓-----------------gBLUP
       GBLUP<- GAPIT(
         Y=Y[training,],
         GD=myGD,
         GM=myGM,
         PCA.total=0,
         CV=myCV,
         model="gBLUP",
         SNP.test=FALSE,
         file.output=F
        )

      inf.y=GBLUP$Pred[GBLUP$Pred[,3]==2,]
      inf.y=merge(Y[!training,],inf.y,by.x=colnames(Y)[1],by.y=colnames(inf.y)[1])
      pred.Y=cor(Y[!training,2],inf.y[,9])
      pred.Y.rep=append(pred.Y.rep,pred.Y)
      #------------------↓↓↓↓↓↓↓↓↓↓-----------------MABLUP
        MABLUP<- GAPIT(
        Y=Y[training,],
        GD=myGD,
        GM=myGM,
        PCA.total=0,
        CV=myCV,
        buspred=T,
        lmpred=F,
        model="BLINK",
        #SNP.test=FALSE,
        file.output=F
      ) 
      inf.y.MA=MABLUP$Pred[MABLUP$Pred[,3]==2,]
      inf.y.MA=merge(Y[!training,],inf.y.MA,by.x=colnames(Y)[1],by.y=colnames(inf.y.MA)[1])
      pred.Y.MA=cor(Y[!training,2],inf.y.MA[,9])
      pred.Y.MA.rep=append(pred.Y.MA.rep,pred.Y.MA)

    }# end of i

      gblupw=mean(pred.Y.rep)
      WG.dif.y[rep,(k-1)]=gblupw

      mablupw=mean(pred.Y.MA.rep)
      WMA.dif.y[rep,(k-1)]=mablupw

      write.csv(WG.dif.y,"weight_gblup.result.csv",quote=F)
      write.csv(WMA.dif.y,"weight_MABLUP.result.csv",quote=F)
         
  } #end of rep
  
} #end of choose Y

rm(list=ls()) 

#install.packages("BGLR")
#library(BGLR)
nrep=100
nfold=5
set.seed(980901)
lasso.result=matrix(NA,nrep,#
	)
bayesB.result=matrix(NA,nrep,#
	)
for(k in 3:8){

  #去掉NA
  Y=myY[complete.cases(myY[,k]),c(1,k)]
  GD=myGD[complete.cases(myY[,k]),]
  CV=myCV[complete.cases(myY[,k]),]
#-------------------------------------------------------------prediction
  lasso=NULL 
  bayesB=NULL

  for(rep in 1:nrep){
    #running gBLUP
    sets=sample(cut(1:nrow(Y),nfold,labels=F),nrow(Y)) 
   
    pred.Y.lasso.rep=NULL
    pred.Y.B.rep=NULL
    for(i in 1:nfold){
      training=sets!=i #选择测试和预测
      bayes.y=Y
      bayes.y[!training,2]=NA
       #--------------------------------------------------------------------------------------------------Bayes lasso
      ETA1=list(list(X=data.frame(GD[,-1],CV[,-1]),model="BL"))
      LASSO<- BGLR(y=bayes.y[,2],ETA=ETA1, nIter=100000, burnIn=60000 )
      pred.Y.lasso=cor(as.numeric(Y[!training,2]),as.numeric(LASSO$yHat[!training]))
      pred.Y.lasso.rep=append(pred.Y.lasso,pred.Y.lasso.rep)   
       #--------------------------------------------------------------------------------------------------Bayes B
      ETA2<-list(list(X=data.frame(GD[,-1],CV[,-1]),model='BayesB'))
      BayesB<-BGLR(y=bayes.y[,2],ETA=ETA2, nIter=100000, burnIn=60000)
      pred.Y.B=cor(as.numeric(Y[!training,2]),as.numeric(BayesB$yHat[!training]))
      pred.Y.B.rep=append(pred.Y.B, pred.Y.B.rep)
    }# end of i

    lasso=append(lasso,mean(pred.Y.lasso.rep))
    bayesB=append(bayesB,mean(pred.Y.B.rep))    
      
  } #end of rep
  lasso.result[,(k-2)]=lasso
  bayesB.result[,(k-2)]=bayesB
  write.csv(lasso.result,"lasso.result.csv",quote=F)
  write.csv(bayesB.result,"bayeaB.result.csv",quote=F) 

} #end of choose Y
#############################################################################KRR
rm(list=ls()) 
#install.packages("KRMM")
#library(KRMM)
nrep=100
nfold=5
set.seed(980901)
KRR.result=matrix(NA,nrep,#
	)

for(k in 3:8){   
	Y=myY[complete.cases(myY[,k]),c(1,k)]
	GD=myGD[complete.cases(myY[,k]),]
	krr=NULL
	for(rep in 1:nrep){
		print(rep)    
		sets=sample(cut(1:nrow(Y),nfold,labels=F),nrow(Y))
		pred.Y.KRR.rep=NULL
		for(i in 1:nfold){
		training=sets!=i
		#-----------------------------------------------------------------KRR
		KRR = Kernel_Ridge_MM(Y_train=Y[training,-1],
                      		Matrix_covariates_train=as.matrix(GD[training,-1]), 
                      		method="RKHS", kernel = "Gaussian")
		pre_krr = Predict_kernel_Ridge_MM(KRR,
                                   		 Matrix_covariates_target=as.matrix(GD[!training,-1]))
		pred.Y.KRR=cor(as.numeric(pre_krr),Y[!training,-1])
		pred.Y.KRR.rep=append(pred.Y.KRR,pred.Y.KRR.rep)
	
     }# end of i
    krr=mean(pred.Y.KRR.rep)
    KRR.result[rep,(k-2)]=krr
    write.csv(KRR.result,"KRR.result.csv",quote=F)
	} #end of rep
   
} #end of choose Y

#############################################################################SVM
rm(list=ls()) 
#install.packages("e1071")
#install.packages("caret")
#library(e1071)
#library(caret)
nrep=100
nfold=5
set.seed(980901)
SVM.result=matrix(NA,nrep,#
	)
for(k in 3:8){
	Y=myY[complete.cases(myY[,k]),c(1,k)]
	GD=myGD[complete.cases(myY[,k]),]
	svm=NULL
	for(rep in 1:nrep){
		sets=sample(cut(1:nrow(Y),nfold,labels=F),nrow(Y))
		pred.Y.SVM.rep=NULL
		for(i in 1:nfold){
			training=sets!=i
			x_train<-GD[training,-1];
			y_train<-Y[training,-1];
			x_test<-GD[!training,-1];
			y_test<-Y[!training,-1]
	       #------------------------------------------------------------------------------------SVM
			SVM<-svm(x_train, y_train, scale=FALSE, kernel='radial') 
			pre_svm <- predict(SVM, x_test);
			pred.Y.SVM=cor(as.numeric(pre_svm),y_test)
			pred.Y.SVM.rep=append(pred.Y.SVM,pred.Y.SVM.rep)

    }# end of i
    svm=mean(pred.Y.SVM.rep)
    SVM.result[rep,(k-2)]=svm
    write.csv(SVM.result,"SVM.result.csv",quote=F)

	} #end of rep
   
} #end of choose Y


########################################################### gblupstar 5k
rm(list=ls())
myGD=myGD
myGM=myGM
myCV=myCV #gender,pca
myY=myY #year old
blup=blup #genetic correlation
rawCV=myCV

cc=cor(myY[,2:7],use="pairwise.complete.obs")
cc[cc==1]=0
cc[is.na(cc)]=0

set.seed(980901)
ll=1
nrep=100
nfold=5
WG.dif.y=matrix(NA,nrep,#
	)

for(v in 1:ncol(cc)){
  max_pheno=which.max(cc[,v])+1
  max_data=blup[,max_pheno]
  myCV=cbind(rawCV,max_data
  
  #选择表型
  for(k in 3:7){
  #去掉NA
  Y=myY[complete.cases(myY[,k]),c(1,k)]
  pred.Y=matrix(NA,nrow(Y),nrep)
  #-------------------------------------------------------------blup
  gblupw=NULL 
  for(rep in 1:nrep){
    #running gBLUP
    sets=sample(cut(1:nrow(Y),nfold,labels=F),nrow(Y))  
    pred.Y.rep=NULL
    
    for(i in 1:nfold){
      training=sets!=i 
      #------------------↓↓↓↓↓↓↓↓↓↓-----------------gBLUP
      GBLUP<- GAPIT(
        Y=Y[training,],
        GD=myGD,
        GM=myGM,
        PCA.total=0,
        CV=myCV,
        model="gBLUP",
        SNP.test=FALSE,
        file.output=F
      )
      
      inf.y=GBLUP$Pred[GBLUP$Pred[,3]==2,]
      inf.y=merge(Y[!training,],inf.y,by.x=colnames(Y)[1],by.y=colnames(inf.y)[1])
      pred.Y=cor(Y[!training,2],inf.y[,9])
      pred.Y.rep=append(pred.Y.rep,pred.Y)
      
    }# end of i
    
    gblupw=mean(pred.Y.rep)
    WG.dif.y[rep,ll]=gblupw
    
    
    
  } #end of rep
  ll=ll+1
  
  } #end of choose Y
  
}

write.csv(WG.dif.y,"MtGBLUP_2-6.result.csv",quote=F)