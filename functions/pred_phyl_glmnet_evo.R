pred_phyl_glmnet_evo=function(tree, traits, opt.thres=FALSE, opt.lam=FALSE,lam=0.05,prob=0.5, mc=TRUE){
 
   ##The function to predict binary traits of AOA. The output is:
  ##1. A list with the predicted values of the NAs of the trait table. If no NAs, predictions are not done.
  ##2. The results of the leave one out cross validation for each gene. Output:Predictions, accuracy, sensitivity (TPR), specificity (TNR)
  #and tuning parameters (lambda and threshold). See below 
  #Two parameters can be optimized:
    #1. Lambda. The penalization parameter of the elastic net regression. If opt.lam=TRUE, otherwise lam=0.05
    #2. Threshold for binary classification of probabilities. If opt.thres=TRUE, otherwise prob=0.5
  ##Requires packages MPSEM, MASS, parallel, PROC, caret, and glmnet
  
  
##THREE AUXILIARY STEPS FOR THE MAIN FUNCTION####

  #1.Names of tree tips and rownames of traits must be same order
spmatch<-match(tree$tip.label, rownames(traits))
traits2<-traits[spmatch,drop=FALSE,]
#If species are not the same in tree and table, stop.
if(identical(rownames(traits2), tree$tip.label)!=TRUE) {
  stop("Species in tree and in trait table are different")
}

#2.The predictions are only made on the columns with NAs. I select them. 
traits3<-traits2[,colSums(is.na(traits2))>0, drop=FALSE]#The predict function will be applied to this object

#3. Create the list with the names from which select for the test dataset.
#I did not find a way that lapply would keep the rownames
labels=rownames(traits2)


##THE TWO INTERNAL FUNCTIONS########

#1. FUNCTION TO PERFORM CROSS VALIDATION. WITH AND WITHOUT LAMBDA AND THRESHOLD OPTIMIZATION####

cross_validation_opt=function(var1,tree,labels, opt.thres, prob, lam, opt.lam){#this function will be later run in lapply. var1 is each column of the trait table
  
  #In case of traits with NAs, need to filter the tree to delete the tips for which I do not have trait info  
  train_data<-var1[!is.na(var1)]#I take only the rows with trait values
  train_data<-as.data.frame(train_data)
  train_names<-labels[!is.na(var1)]#Select the row names from the original rownames
  rownames(train_data)<-train_names#assignn them to the train dataset
  
  #Prune the tree 
  spmatch_val<-match(tree$tip.label, train_names)#select the tips that are not in the train dataset
  tree_val<-drop.tip(tree, tree$tip.label[is.na(spmatch_val)])#delete the tips that are not in the train dataset
  spmatch_val<-match(tree_val$tip.label, train_names)#Match tips and names of the train dataset
  data_val<-train_data[spmatch_val,drop=FALSE,]#Reorder tips and train dataset names
  
  #Check that it is correct  
  if(identical(rownames(data_val), tree_val$tip.label)!=TRUE){
    stop("species in tree and trait table in the cross validation differ")
  }
  
  #Optimizing the lambda factor for the later leave on out predictions
  #Obtaining the range of lambdas for the whole dataset (it is almost the same as that after deleting one row, I use it)
  opt.pgraph<-Phylo2DirectedGraph(tree_val)#I convert the tree into a graph object for the PEM
  
  train_PEM_val<-PEM.fitSimple(#I obtain the PEMs, calculating the a parameter
    y = data_val[,1],
    x = NULL,
    w = opt.pgraph) 
  
  train_eigen_val<-train_PEM_val$u#I extract the eigenvectors
  
  el_glm_val_temp<-glmnet(#I build the elastic net regression with the eigenvectors as variables and the gene as response
    x= train_eigen_val,
    y = data_val[,1],
    family = "binomial",
    alpha=0.5)   #alpha 0.5 combines the lasso and ridge shrinkage. 
  
  set_lambdas<-el_glm_val_temp$lambda#I obtain all the values of lambda for the iteration
  if(opt.lam==TRUE){s=set_lambdas}else{s=lam} #I will make the predictions for all lambdas or only the user-defined
  cross_val<-matrix(ncol=length(s), nrow=nrow(data_val))
  rownames(cross_val)<-rownames(data_val)

  for (i in 1:nrow(data_val)){#This is the LOOCV
    
    #1.Get graph locations for the target species and the graph info for the other species
    target_graph_val<-getGraphLocations(tree_val,
                                        targets = rownames(data_val)[i])
    
    #2. Calculate the PEM object, calculating the a parameter. 
    
    train_PEM_val<-PEM.fitSimple(
      y = data_val[-i,1],
      x = NULL,
      w = target_graph_val$x) 
    
    train_eigen_val<-train_PEM_val$u#Extract the igenvectors
    
    #3. Obtain the loadings of the target species projected on the training PEM
    target_scores_val<-Locations2PEMscores(object=train_PEM_val,
                                           gsc=target_graph_val)$scores
    #4. Build the elastic net regression
    
    el_glm_val<-glmnet(
      x= train_eigen_val,
      y = data_val[-i,1],
      family = "binomial",
      alpha=0.5) 
    
    #5. Make predictions on the target
    pred_val<-predict(el_glm_val, target_scores_val, s=s, type="response")
    cross_val[i,]<-pred_val
  }
  
  base_accuracy<-c(0)
  for (i in 1:ncol(cross_val)){
    #i=1
    if (opt.thres==TRUE){#Optimizing the threshold value
      roc_score<-suppressMessages(roc(response=data_val[,1], predictor=cross_val[,i],direction="<" ))#I run the ROC curve with the observed and predicted values
      best_thres<-coords(roc_score, "best", best.method="youden")#I select the threshold optimizing for the youden index Y=(sensitivity+specificy)-1
      prob<-best_thres[1,1]#I obtain the optimized threshold
      compar<-cbind(data_val,ifelse(cross_val[,i]>prob,1,0),cross_val[,i])#Merge observed, classified predicted and predicted probability
      colnames(compar)=c("obs", "pred", "pred_exact")#Rename columns
      av=mean(compar[,1]==compar[,2])#Calculate accuracy of predictions. Compare predicted vs observed
      if (av>base_accuracy){
        base_accuracy<-av
        opt_lam<-av#Save accuracy on the vector
        opt_pred<-compar#Save predicted values
        sensit<-best_thres[1,"sensitivity"]#Save sensitivity
        spec<-best_thres[1,"specificity"]#Save specificity
        thresh<-prob
        s_opt<-s[i]}#Save optimized threshold
      
    }else{
      prob=prob#Delete this line?#Sets the probability to the one specified by user in the main function
      roc_score<-suppressMessages(roc(response=data_val[,1], predictor=cross_val[,i], direction="<"))#Run the ROC curve for sensitivity and specificity values
      compar<-cbind(data_val,ifelse(cross_val[,i]>prob,1,0),cross_val[,i])#Merge predicted and observed
      colnames(compar)=c("obs", "pred", "pred_exact")#Rename columns
      av=mean(compar[,1]==compar[,2])#Calculate accuracy of predictions. Compare predicted vs observed
      if (av>base_accuracy){
        base_accuracy<-av
        opt_lam<-av#Save accuracy on the vector
        opt_pred<-compar#Save predicted values
        sensit<-coords(roc_score, x=prob, input="threshold", ret="sensitivity")[1,1]#Save sensitivity
        spec<-coords(roc_score, x=prob, input="threshold", ret="specificity")[1,1]#Save specificity
        thresh<-prob
        s_opt<-s[i]}
    }
    
    
  }#End of loop for choosing optimum lambda and parameter
  
  
  #Now I summarize the results of the LOOCV
  
  list("cross_predictions"=opt_pred,"cross_mean"=opt_lam, "opt_lambda"=s_opt, "threshold"=thresh, "sensitivity"=sensit, "specificity"=spec)
}


#2. FUNCTION FOR THE PREDICTIONS WITH AND WITHOUT OPTIMIZED LAMBDA AND THRESHOLD####
predict_bin_opt=function(nam,traits_sel,tree, labels, opt.thres, s, thres, lam, prob){
  #non NA values to fit the model
  var1<-traits_sel[,nam, drop=FALSE]
  train_data<-var1[!is.na(var1),,drop=FALSE]
  
  #The names of the genomes with NAs. we will get the locations of them
  test_names<-labels[is.na(var1)] 
  
  #1.Get graph locations for the target species and the graph info for the other species
  target_graph<-getGraphLocations(tree,
                                  targets = test_names)
  
  #2. Calculate the PEM object. Calculate the a parameter 
  
  #train_data<-as.data.frame(train_data)
  #rownames(train_data)<-labels[!is.na(var1)]
  #print(train_data)
  
  train_PEM<-PEM.fitSimple(
    y = train_data[,1],
    x = NULL,
    w = target_graph$x) 
  #We extract the loadings on the evectors to fit the model
  train_eigen<-train_PEM$u
  
  #3. Obtain the loadings of the target species projected on the training PEM
  target_scores<-Locations2PEMscores(object=train_PEM,
                                     gsc=target_graph)$scores
  #print(target_scores)
  
  #4. Build the elastic net regression
  #parameters
  #alpha = 0.5
  #s=lam
  
  #4.Building the glmnet
  el_glm<-glmnet(
    x= train_eigen,
    y = train_data[,1],
    family = "binomial",
    alpha=0.5) 
  
  #5. Make predictions on the target
  #Obtain predicted probabilities
  if(opt.lam==TRUE){
    s=s[nam]
  }else{s=lam}
  #print(s)
  pred<-predict(el_glm, target_scores, s=s, type="response")
  #Convert probabilities to 1/0 and merge
  if(opt.thres==TRUE){
    prob=thres[nam]
  }else{prob=prob}#Delete the else
  #print(prob)
  pred_class<-ifelse(pred>prob,1,0)
  list("pred_class"=pred_class, "pred_prob"=pred, "lambda"=s[[1]], "threshold"=prob[[1]])
}


##THE 2 DIFFERENT OPTIONS OF THE FUNCTION:WITH AND WITHOUT MULTICORE PROCESSING####

  if(mc==TRUE){#Multicore
    validation_opt<-mclapply(traits2,mc.cores=detectCores(),cross_validation_opt, tree=tree,labels=labels, opt.thres=opt.thres, prob=prob, lam=lam, opt.lam=opt.lam)#I run the optimization
    s=sapply(validation_opt, function(x) x[["opt_lambda"]])
    thres=sapply(validation_opt, function(x) x[["threshold"]])
    traits_list<-colnames(traits3)
    names(traits_list)<-traits_list
    predictions<-mclapply(traits_list,mc.cores=detectCores(),predict_bin_opt, traits_sel=traits3,tree=tree,labels=labels,opt.thres=opt.thres, thres=thres, s=s, lam=lam, prob=prob)
  }else{#No multicore
    validation_opt<-lapply(traits2,cross_validation_opt, tree=tree,labels=labels, opt.thres=opt.thres, prob=prob, lam=lam, opt.lam=opt.lam)#I run the optimization
    s=sapply(validation_opt, function(x) x[["opt_lambda"]])
    thres=sapply(validation_opt, function(x) x[["threshold"]])
    traits_list<-colnames(traits3)
    names(traits_list)<-traits_list
    predictions<-lapply(traits_list,predict_bin_opt, traits_sel=traits3,tree=tree,labels=labels,opt.thres=opt.thres, thres=thres, s=s, lam=lam, prob=prob)
  }
  #Cleaning the output of the whole function
  mean_crossval<-lapply(validation_opt, function(x) x[["cross_mean"]])
  predicted_val<-lapply(validation_opt, function(x) x[["cross_predictions"]])
  opt_lambda<-lapply(validation_opt, function(x) x[["opt_lambda"]])
  threshold_fin<-lapply(validation_opt, function(x) x[["threshold"]])
  sens_fin<-lapply(validation_opt, function(x) x[["sensitivity"]])
  specs_fin<-lapply(validation_opt, function(x) x[["specificity"]])
  list("predictions"=predictions, "CV_predictions"=predicted_val, "CV_accuracy"=mean_crossval,"CV_opt_lambda"=opt_lambda, "CV_threshold"=threshold_fin, "CV_sensitivity"=sens_fin, "CV_specificity"=specs_fin)

}#End of function



