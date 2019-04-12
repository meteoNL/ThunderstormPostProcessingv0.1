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

print("rangertrying active for")
print("LT_i equals:")
print(LT_i)
print("------")
print("VT_i equals:")
print(VT_i)

#hyperparameters settings
numbtree = 500
m = 2
min_length = 1
node_size = 15#50

# import the data frame and generate training and testing set
minpredictant = 1.5 #1.5 discharges
setwd("/usr/people/groote")
ObsPV = read.csv(file = "full_final00z_dataset2.csv")
colnames(ObsPV)[c(3,4,5,6,7,10,11)] <- c("Year","validdate","Month","Day","validtime2","region","leadtime_count")
ObsPV$leadtime_count = ObsPV$leadtime_count
ObsPV$validdate = ObsPV$validdate+ObsPV$Year*10000
ObsPV <- filter(ObsPV, region < 25, Ndischarge > -100) #last one to prevent issues for missing values
#ObsPV <- filter(ObsPV, region != 1)

#valid time, regions, lead time values and apply selection for one combination of VT and LT
years = c(as.numeric(unique(ObsPV$Year)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[LT_i]
#VT = unique(ObsPV$validtime)[VT_i]
VT = paste0("Init_00z_few_finalcross",length(varindex_shell))
regions = c(unique(ObsPV$region))
#ObsPV=cbind(ObsPV, selector = as.numeric(as.character(VT)==as.character(ObsPV$validtime))*as.numeric(as.character(LT)==as.character(ObsPV$leadtime_count)))
ObsPV=cbind(ObsPV, selector = as.numeric(as.character(LT)==as.character(ObsPV$leadtime_count)))
climset <- filter(ObsPV, selector == 1)

print(dim(climset))

# set indices of the predictors and predictand
occurence = as.numeric(climset[625]>minpredictant)
climset = data.frame(climset, occurence)
orig_varindex = varindex_shell
predictant_ind = length(climset)

### Above this point, the settings for a run have been defined!! #####
# declare some names
remove_variable = c()
overall_scores = data.frame()
importances_df = data.frame()

### Functions defined: ####
qrf_procedure <- function(train_set, test_set, predictant_index, varindexset, m_hyp, ntree_hyp, node_size_hyp, min_length_setting, wval, yval){
  #this function does the full QRF fitting procedure for thunderstorm occurence (yes/no), elimination and predictions with test data

  #declare some names to save results
  overall_scores_local = data.frame()
  importances_save = data.frame()
  result=list()

  #loop for any number of predictors from 2 to total available number
  while (length(train_set[varindexset])>min_length_setting){
    #predictors still in; predictand
    pot_preds <- names(train_set[varindexset])
    predictant <- names(train_set[predictant_index])

    #assign the variables and formula for the fit and do a fit
    formul = train_set[predictant]~.
    assign("train_set", train_set, .GlobalEnv)
    assign("predictant", predictant, .GlobalEnv)
    fit1=ranger(formul, data = data.frame(train_set[c(pot_preds)]), num.trees = ntree_hyp, mtry = min(m_hyp, (length(varindexset))), min.node.size = node_size_hyp, importance = "permutation", quantreg = TRUE)    #plot(qrf_fit)

    #remove worst variable, based on permutation
    remove_variable = varindexset[fit1$variable.importance == min(fit1$variable.importance)][1]
    varindexset = varindexset[varindexset != remove_variable]
  #  print("removed variable:")
  #  print(remove_variable)

    ######
    #predict with fit the probabilities and remember the five most important predictors of the fit
    fit1_pred = predict(fit1, data = test_set)
    test_importance_save = data.frame(npred = (length(varindexset)+1), mtry = m_hyp, effmtry = min(m_hyp, (length(varindexset))), min_n_size = node_size_hyp, w = wval, y = yval, importances = sort(fit1$variable.importance)[(length(fit1$variable.importance)-min(5,length(fit1$variable.importance))):length(fit1$variable.importance)])
    #make data frame with predictions
    qrf_pred = data.frame(prob = fit1_pred$predictions, occurence = test_set[predictant_index], region = test_set["region"], npred = (length(varindexset)), mtry = m_hyp, effmtry = min(m_hyp, (length(varindexset))), min_n_size = node_size_hyp, w = wval, y = yval)

    #remove globals and collect the data
    rm("predictant", "train_set")
    overall_scores_local = rbind(overall_scores_local, qrf_pred)
    importances_save = rbind(importances_save, test_importance_save)
  }

  #return results for all numbers of predictors
  result$overall_scores_local = overall_scores_local
  result$importances_save = importances_save
  return(result)
}

# -----------------------------------------------------------------------------
### End of the functions ####

q = 1
#do procedure for 9-fold cross validation: select two years of data from data frame
for (y in years){
  train_y = filter(climset, Year != y)
  test_y = filter(climset, Year == y)

  varindex = orig_varindex

  #do the fits procedure
  result = qrf_procedure(train_y, test_y, predictant_ind, varindex, m, numbtree, node_size, min_length, 1, y)

  #save the predictions and important predictors
  overall_scores = rbind(overall_scores, result$overall_scores_local)
  importances_df = rbind(importances_df, result$importances_save)
  write.csv(importances_df, file=paste0("importances_newrandom_",VT,"_LT_",LT,".csv"))
}

#this table will contain all Brier Scores for all hyperparameters (mtry, node_size, npredictors), per region, per threshold and per testsubset
#calculate Brier (Skill) Scores for 9-fold
setwd("/usr/people/groote/ThunderstormPostProcessingv1/rangerres")
qrf_ss <- overall_scores %>% group_by(region, npred, mtry, min_n_size) %>% summarise(bs = brier(obs = occurence, pred = prob, bins = FALSE)$ss)
qrf_bs <- overall_scores %>% group_by(region, npred, mtry, min_n_size) %>% summarise(bs = brier(obs = occurence, pred = prob, bins = FALSE)$bs)

#write them to files
write.csv(qrf_ss, file=paste0("qrf_ss_imp_newrandom_",VT,"_LT_",LT,".csv"))
write.csv(qrf_bs, file=paste0("qrf_bs_imp_newrandom_",VT,"_LT_",LT,".csv"))
write.csv(importances_df, file=paste0("importances_newrandom_",VT,"_LT_",LT,".csv"))

