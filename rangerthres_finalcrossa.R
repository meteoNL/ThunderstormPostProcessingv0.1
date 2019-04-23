rm(list=ls()[! ls() %in% c("LT_i","VT_i","LT_val","VT_val","varindex_shell")])
# load libraries:
# if they need to be installed that can be done by: install.packages("dplyr")
library(dplyr)
library(tidyr)
library(ranger)
library(verification)
library(ggplot2)
library(tree)
library(SpecsVerification)
library(parallel)

print("rangerthres active for")
print("LT_i equals:")
print(LT_i)
print("------")
print("VT_i equals:")
print(VT_i)

#hyperparameters settings
numbtree = 500
m= 2
min_length = 1
node_size = 9#50

# import the data frame and generate training and testing set
minpredictant = 1.5 #1.5 discharges
setwd("/usr/people/groote")
ObsPV = read.csv(file = "full_final00z_dataset2.csv")
colnames(ObsPV)[c(3,4,5,6,7,10,11)] <- c("Year","validdate","Month","Day","validtime2","region","leadtime_count")
ObsPV <- filter(ObsPV, (is.na(Ndischarge)+is.na(Bradbury_min))==0)
ObsPV$leadtime_count = ObsPV$leadtime_count
ObsPV$validdate = ObsPV$validdate+ObsPV$Year*10000
nmem = 25 #number of members for ensemble CRPS
p=0.25

#valid time, regions, lead time values and apply selection for one combination of VT and LT
years = c(as.numeric(unique(ObsPV$Year)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[LT_i]
#VT = unique(ObsPV$validtime)[VT_i]
VT = paste0("Init_00z_finallcross_wconf",length(varindex_shell))
regions = c(unique(ObsPV$region))
#ObsPV=cbind(ObsPV, selector = as.numeric(as.character(VT)==as.character(ObsPV$validtime))*as.numeric(as.character(LT)==as.character(ObsPV$leadtime_count)))
ObsPV=cbind(ObsPV, selector = as.numeric(as.character(LT)==as.character(ObsPV$leadtime_count)))
climset <- filter(ObsPV, selector == 1 & Ndischarge > minpredictant)

# initial predictor set containing all predictors; predictand and evaluation thresholds
orig_varindex = varindex_shell
predictant_ind = 627
thres_eval = seq(25,400,25)
climset[predictant_ind]=climset[predictant_ind]^p

### Above this point, the settings for a run have been defined!! #####
# declare some names
remove_variable = c()
overall_scores = data.frame()
overall_scores_quan = data.frame()
importances_dataset = data.frame()

### Functions defined: ####
qrf_procedure_thres <- function(train_set, test_set, predictant_index, varindexset, m_hyp, ntree_hyp, node_size_hyp, min_length_setting, wval, yval){
  #this function does the full QRF fitting procedure for threshold predictions, elimination and predictions with test data

  #declare some names to save results
  overall_scores_local = data.frame()
  overall_quan = data.frame()
  importances_save = data.frame()
  result = list()

  #loop for any number of predictors from 2 to total available number
  while (length(train_set[varindexset])>min_length_setting){
    #predictors still in; predictand
    pot_preds <- names(train_set[varindexset])
    predictant <- names(train_set[predictant_index])

    #assign the variables and formula for the fit and then do a fit
    formul = train_set[predictant]~.
    assign("train_set", train_set, .GlobalEnv)
    assign("predictant", predictant, .GlobalEnv)
    fit1=ranger(formul, data = data.frame(train_set[c(pot_preds)]), num.trees = ntree_hyp, mtry = min(m_hyp, length(varindexset)), min.node.size = node_size_hyp, importance = "permutation", quantreg = TRUE)    #plot(qrf_fit)

    #remove worst variable, based on permutation
    remove_variable = varindexset[fit1$variable.importance == min(fit1$variable.importance)][1]
    varindexset = varindexset[varindexset != remove_variable]
   # print("removed variable:")
  #  print(remove_variable)

    #evaluate continuous ranked probability score
    qrf_pred_quan = predict(fit1, data.frame(test_set[, pot_preds]), type = "quantiles", quantiles = seq(0.5,nmem-0.5)/nmem)$predictions
    qrf_pred_quan = data.frame(qrf_pred_quan, occurence = test_set[predictant_index], region = test_set["region"], npred = length(varindexset), mtry = m_hyp, effmtry = min(m_hyp, length(varindexset)), min_n_size = node_size_hyp, w = wval, y = yval)
    ######

    qrf_pred = data.frame()
    #do evaluation with thresholds
    for(threshold in thres_eval){
      #predict discrete CDF
      fit1_pred_thres = predict(fit1, data = test_set, predict.all = TRUE)

      #get observations
      obs = as.numeric(test_set[predictant_index]^(1/p) > threshold)

      #calculate probabilities from discrete CDF and make data frame with model settigns and regions in it
      probs = rowMeans(fit1_pred_thres$predictions^(1/p)>threshold)
      qrf_pred = rbind(qrf_pred, data.frame(probability = probs, observed = obs, region = test_set["region"], thres = threshold, npred = (length(varindexset)+1), mtry = m_hyp, effmtry = min(m_hyp, (length(varindexset)+1)), min_n_size = node_size_hyp, w = wval, y = yval, validdate = test_set$validdate))
    }

    #remember the five most important predictors of the fit
    test_importance_save = data.frame(npred = (length(varindexset)+1), mtry = m_hyp, effmtry = min(m_hyp, (length(varindexset)+1)), min_n_size = node_size_hyp, w = wval, y = yval, importances = sort(fit1$variable.importance)[1:length(fit1$variable.importance)])

    #remove assigned variables
    rm("predictant", "train_set")

    #add predictions (probabilities & ensemble members for CRPS) for this number of predictors to other predictions for any number of predictors and the most important predictors
    overall_scores_local = rbind(overall_scores_local, qrf_pred)
    overall_quan = rbind(overall_quan, qrf_pred_quan)
    importances_save = rbind(importances_save, test_importance_save)
  }

  #return results for all numbers of predictors
  result$overall_scores_local = overall_scores_local
  result$overall_quan = overall_quan
  result$importances_save = importances_save
  return(result)
}

make_subset <- function(i){
  verifdate = verifdates[i]
  subset = predictiondataframe[predictiondataframe$validdate == verifdate,]
  return(subset)

}
do_bootstrap <- function(i){
  datelist = unique(predictiondataframe$validdate)
  set.seed(i+950)
  randdates = round(0.5+length(unique(predictiondataframe$validdate))*runif(length(unique(predictiondataframe$validdate))))
  verifdates = datelist[randdates]
  verifydataframe = data.frame()
  #assign("verifydataframe", verifydataframe, .GlobalEnv)
  assign("verifdates", verifdates, .GlobalEnv)
  verifydataframe = do.call(rbind,lapply(seq(1,length(verifdates)),make_subset))
  skillscore = brier(verifydataframe$observed, verifydataframe$probability, rep(mean(predictiondataframe$observed),length(predictiondataframe$observed)), bins = FALSE)$ss
  #print(skillscore)
  return(skillscore)
}
block_bootstrapping <- function(predictiondataframe,ntimes,interval){
  skillscores = data.frame()
  assign("predictiondataframe", predictiondataframe, .GlobalEnv)
  print(predictiondataframe$thres)
  clusterExport(cl = cls, c("make_subset", "do_bootstrap","predictiondataframe","brier"))
  skillscores = unlist(parLapply(cl = cls, X = seq(1,ntimes), fun = do_bootstrap))
  #print(skillscores)
  funcresult = list()
  funcresult$up = quantile(skillscores, (0.5+interval/2))
  funcresult$low = quantile(skillscores, (0.5-interval/2))
  rm("predictiondataframe")
  return(funcresult)

}
# -----------------------------------------------------------------------------
### End of the functions ####

q = 1

#do procedure for 9-fold cross validation: select two years of data from data frame
for (y in years){

  #select first 2 years out of three
  train_y = filter(climset, Year != y)
  test_y = filter(climset, Year == y)
  varindex = orig_varindex

  result = qrf_procedure_thres(train_y, test_y, predictant_ind, varindex, m, numbtree, node_size, min_length, 1, y)
  overall_scores = rbind(overall_scores, result$overall_scores_local)
  overall_scores_quan = rbind(overall_scores_quan, result$overall_quan)
  importances_dataset = rbind(importances_dataset, result$importances_save)


}
#this table will contain for all hyperparameters (mtry, node_size, npredictors), per region, per threshold and per testsubset
setwd("/usr/people/groote/ThunderstormPostProcessingv1/rangerres")
#write.csv(overall_scores, file = "overall scores3.csv") #this file is very very large

#calculate Brier Skill Score
cls <- makeCluster(getOption("cl.cores", 8))
qrf_ss <- overall_scores %>% group_by(npred, mtry, min_n_size, thres) %>% summarise(bs = brier(obs = observed, pred = probability, bins = FALSE)$ss, bs_conf = list(block_bootstrapping(data.frame(observed, probability, validdate),1000,0.95)))
stopCluster(cls)
qrf_bs <- overall_scores %>% group_by(npred, mtry, min_n_size, thres) %>% summarise(bs = brier(obs = observed, pred = probability, bins = FALSE)$bs)
overall_scores_quan = data.frame(overall_scores_quan, newcol = overall_scores_quan$npred*10000+overall_scores_quan$mtry*100+overall_scores_quan$min_n_size)
qrf_ss = data.frame(npred = qrf_ss$npred, thres = qrf_ss$thres, mtry = qrf_ss$mtry, min_n_size = qrf_ss$min_n_size, bss = qrf_ss$bs, bss_upper=unlist(qrf_ss$bs_conf)[seq(1,2*length(qrf_ss$bs_conf),2)],bss_lower = unlist(qrf_ss$bs_conf)[seq(2,2*length(qrf_ss$bs_conf),2)])
#calculate CRPS score
qrf_crps = data.frame()
for(val in unique(overall_scores_quan$newcol)){
  subset = filter(overall_scores_quan, newcol == val)
  newscore = data.frame(crps = mean(EnsCrps(as.matrix(subset[1:nmem]),as.numeric(unlist(subset[nmem+1])))), npred = round(val/10000), mtry = round(val%%10000/100), min_n_size=val%%100)
  qrf_crps = rbind(qrf_crps, newscore)
}
refcrps = mean(EnsCrps(t(matrix(rep(t(climset[predictant_ind]),dim(climset[predictant_ind])[1]), nrow=dim(climset[predictant_ind])[1])),as.matrix(climset[predictant_ind])))
sksc=1-(qrf_crps$crps)/refcrps
qrf_crps = cbind(qrf_crps, skillscore = sksc)

#write all scores and list with important predictors for all of the models to a file
write.csv(qrf_ss, file=paste0("qrf_thresholds_ss_imp_newrandom_",VT,"_LT_",LT,".csv"))
write.csv(qrf_bs, file=paste0("qrf_thresholds_bs_imp_newrandom_",VT,"_LT_",LT,".csv"))
write.csv(qrf_crps, file=paste0("qrf_thresholds_crps_imp_newrandom_",VT,"_LT_",LT,".csv"))
write.csv(importances_dataset, file=paste0("imp_",VT,"_LT_",LT,".csv"))
