rm(list=ls(all=TRUE))
# load libraries:
# if they need to be installed that can be done by: install.packages("dplyr")
library(dplyr)
library(tidyr)
library(ranger)
library(verification)
library(ggplot2)
library(tree)
library(SpecsVerification)

#hyperparameters settings
numbtree = 500
m_settings = c(2, 6, 10)
min_length = 1
node_size_settings = c(3, 9, 15)#50

# import the data frame and generate training and testing set
minpredictant = 1.5 #1.5 discharges
setwd("/usr/people/groote")
ObsPV = read.csv(file = "Thunderstorm_radar_merged.csv")
numbsubset = 3
p=0.25
nmem = 10 #number of members for ensemble CRPS

#valid time, regions, lead time values and apply selection for one combination of VT and LT
years = c(as.numeric(unique(ObsPV$Year.x)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[1]
VT = unique(ObsPV$validtime.x)[2]
regions = c(unique(ObsPV$region.x))
climset <- filter(ObsPV, validtime.x == VT & leadtime_count == LT & Ndischarge > minpredictant)

# initial predictor set containing all predictors; predictand and evaluation thresholds
orig_varindex = c(seq(18,37),seq(39,79),seq(81,101))
predictant_ind = 113
climset[predictant_ind]=climset[predictant_ind]^p
thres_eval = seq(1.8,3.0,0.05)^4

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
    print("removed variable:")
    print(remove_variable)

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
      probs = rowMeans(fit1_pred_thres$predictions^(1/p) > threshold)
      qrf_pred = rbind(qrf_pred, data.frame(probability = probs, observed = obs, region = test_set["region"], thres = threshold, npred = (length(varindexset)+1), mtry = m_hyp, effmtry = min(m_hyp, (length(varindexset)+1)), min_n_size = node_size_hyp, w = wval, y = yval))
    }

    #remember the five most important predictors of the fit
    test_importance_save = data.frame(npred = (length(varindexset)+1), mtry = m_hyp, effmtry = min(m_hyp, (length(varindexset)+1)), min_n_size = node_size_hyp, w = wval, y = yval, importances = sort(fit1$variable.importance)[(length(fit1$variable.importance)-min(5,length(fit1$variable.importance))):length(fit1$variable.importance)])

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

# -----------------------------------------------------------------------------
### End of the functions ####

q = 1

#do procedure for 9-fold cross validation: select two years of data from data frame
for (y in years){
  #select first 2 years out of three
  train_y = filter(climset, Year != y)
  test_y = filter(climset, Year == y)
  set.seed(15+seq(years)[q]) #for reproducability purposes

  #then make three random subsets based on date
  testdf = data.frame(validdate = unique(train_y$validdate), subset = round(runif(unique(train_y$validdate))*numbsubset+0.5))
  train_sub <- left_join(train_y, testdf, by = names(testdf)[1])
  for (w in seq(numbsubset)){
    print(y)
    print(w)

    #declare testing and training datasets
    train_q = filter(train_sub, subset != w)
    test_q = filter(train_sub, subset == w)
    #print(train_q)
    #do a loop over the hyperparameter settings
    for (node_size in node_size_settings){
      #print("node size equals: ")
      #print(node_size)
      for (m in m_settings){
        #print("m equals: ")
        #print(m)

        #reset to all predictors
        varindex = orig_varindex

        #apply QRF fitting procedures
        result = qrf_procedure_thres(train_q, test_q, predictant_ind, varindex, m, numbtree, node_size, min_length, w, y)

        #save predictions for ENS CRPS and Brier Skill Score calculation
        overall_scores = rbind(overall_scores, result$overall_scores_local)
        overall_scores_quan = rbind(overall_scores_quan, result$overall_quan)
        importances_dataset = rbind(importances_dataset, result$importances_save)
      }
    }
  }
  q = q + 1
}
#this table will contain all brier scores for all hyperparameters (mtry, node_size, npredictors), per region, per threshold and per testsubset
setwd("/usr/people/groote/ThunderstormPostProcessingv1/rangerres")
#write.csv(overall_scores, file = "overall scores3.csv") #this file is very very large

#calculate Brier Skill Score
qrf_ss <- overall_scores %>% group_by(npred, mtry, min_n_size, thres) %>% summarise(bs = brier(obs = observed, pred = probability, bins = FALSE)$ss)
qrf_bs <- overall_scores %>% group_by(npred, mtry, min_n_size, thres) %>% summarise(bs = brier(obs = observed, pred = probability, bins = FALSE)$bs)
overall_scores_quan = data.frame(overall_scores_quan, newcol = overall_scores_quan$npred*10000+overall_scores_quan$mtry*100+overall_scores_quan$min_n_size)

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

#-----------------------------------------------------------------
## Testing the functions

set.seed(712)
x1 = rnorm(1000,0,5)
x2 = rnorm(1000,0,5)
x3 = rnorm(1000,0,5)
pert = rnorm(1000,0,5)
yvalues = x1*25+x2*100+pert
test_df = data.frame(x1, x2, x3, yvalues)
library(devtools)
library(testthat)
usethis::use_testthat()
form = test_df["yvalues"] ~.
qrf_fit_test <- ranger(form, data = data.frame(test_df[seq(1,3)]), num.trees=250, mtry = 1, min.node.size = 5, quantreg = TRUE, importance ="permutation")
importance_table = qrf_fit_test$variable.importance
test_that("Predictor ranking of dataset",{
  expect_gt(importance_table[2],importance_table[1])
  expect_gt(importance_table[1],importance_table[3])
})
test_that("Dataset complete?",{
  expect_equal(climset %>% arrange(Year, Month, Day), rbind(train_y, test_y) %>% arrange(Year, Month, Day))
  expect_equal(rbind(train_q,test_q) %>% arrange(Year, Month, Day, region),train_sub %>% arrange(Year, Month, Day, region))
})

test_that("Random subset numbers",{
  expect_equal(min(testdf$subset),1)
  expect_equal(max(testdf$subset),3)
})
