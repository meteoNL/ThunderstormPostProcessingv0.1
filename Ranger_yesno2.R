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

#hyperparameters
numbtree = 250
m_settings = c(2, 6, 10)
min_length = 2
node_size_settings = c(3, 9, 15)#50

# import the data frame and generate training and testing set
minpredictant = 1.5 #1.5 discharges
setwd("/usr/people/groote")
ObsPV = read.csv(file = "Thunderstorm_radar_merged.csv")
#ObsPV <- filter(ObsPV, region != 1)
numbsubset = 3
years = c(as.numeric(unique(ObsPV$Year.x)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[1]
VT = unique(ObsPV$validtime.x)[2]
regions = c(unique(ObsPV$region.x))
climset <- filter(ObsPV, validtime.x == VT & leadtime_count == LT)
occurence = as.numeric(climset[107]>minpredictant)
climset = data.frame(climset, occurence)
orig_varindex = c(seq(18,101))
predictant_ind = 114

remove_variable = c()
overall_scores = data.frame()
importances_df = data.frame()
qrf_procedure <- function(train_set, test_set, predictant_index, varindexset, m_hyp, ntree_hyp, node_size_hyp, min_length_setting, wval, yval){
  overall_scores_local = data.frame()
  importances_save = data.frame()
  result=list()
  while (length(train_set[varindexset])>min_length_setting){
    #predictors still in
    pot_preds <- names(train_set[varindexset])
    predictant <- names(train_set[predictant_index])

    #assign the variables and formula for the fit and do a fit
    formul = train_set[predictant]~.
    assign("train_set", train_set, .GlobalEnv)
    assign("predictant", predictant, .GlobalEnv)
    fit1=ranger(formul, data = data.frame(train_set[c(pot_preds)]), num.trees = ntree_hyp, mtry = min(m_hyp, length(varindexset)), min.node.size = node_size_hyp, importance = "impurity", quantreg = TRUE)    #plot(qrf_fit)

    #remove variable
    remove_variable = varindexset[fit1$variable.importance == min(fit1$variable.importance)][1]
    varindexset = varindexset[varindexset != remove_variable]
    print("removed variable:")
    print(remove_variable)

    #evaluate continuous ranked probability score
    #qrf_pred_quan = predict(qrf_fit, data.frame(test_set[, pot_preds]), what = seq(0.05,0.95,0.1))

    ######
    fit1_pred = predict(fit1, data = test_set)
    test_importance_save = data.frame(npred = length(varindexset), mtry = m_hyp, effmtry = min(m_hyp, length(varindexset)), min_n_size = node_size_hyp, w = wval, y = yval, importances = sort(fit1$variable.importance)[(length(fit1$variable.importance)-min(5,length(fit1$variable.importance))):length(fit1$variable.importance)])
    print(test_importance_save)
    qrf_pred = data.frame(prob = fit1_pred$predictions, occurence = test_set[predictant_index], region = test_set["region"], npred = length(varindexset), mtry = m_hyp, effmtry = min(m_hyp, length(varindexset)), min_n_size = node_size_hyp, w = wval, y = yval)

    #remove globals and collect the data
    rm("predictant", "train_set")
    result$overall_scores_local = rbind(overall_scores_local, qrf_pred)
    result$importances_save = rbind(importances_save, test_importance_save)
  }
  return(result)
}

q = 1
for (y in years){
  train_y = filter(climset, Year != y)
  test_y = filter(climset, Year == y)
  set.seed(seq(years)[q])
  randomsubset = round(runif(nrow(train_y))*numbsubset+0.5)
  train_sub = cbind(train_y, subset = randomsubset)
  for (w in seq(numbsubset)){
    print(y)
    print(w)
    train_q = filter(train_sub, subset != w)
    test_q = filter(train_sub, subset == w)
    #print(train_q)
    for (node_size in node_size_settings){
      print("node size equals: ")
      print(node_size)
      for (m in m_settings){
        print("m equals: ")
        print(m)
        varindex = orig_varindex
        result = qrf_procedure(train_q, test_q, predictant_ind, varindex, m, numbtree, node_size, min_length, w, y)
        overall_scores = rbind(overall_scores, result$overall_scores_local)
        importances_df = rbind(importances_df, result$importances_save)
      }
    }
  }
  q = q + 1
}
#this table contains all brier scores for all hyperparameters (mtry, node_size, npredictors), per region, per threshold and per testsubset
setwd("/usr/people/groote/ThunderstormPostProcessingv1/rangerres")
qrf_ss <- overall_scores %>% group_by(region, npred, mtry, min_n_size) %>% summarise(bs = brier(obs = occurence, pred = prob, bins = FALSE)$ss)
qrf_bs <- overall_scores %>% group_by(region, npred, mtry, min_n_size) %>% summarise(bs = brier(obs = occurence, pred = prob, bins = FALSE)$bs)
write.csv(qrf_ss, file="qrf_ss_imp.csv")
write.csv(qrf_bs, file="qrf_bs_imp.csv")
write.csv(importance_df, file="importances.csv")

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
qrf_fit_test <- ranger(form, data = data.frame(test_df[seq(1,3)]), num.trees=250, mtry = 1, min.node.size = 5, quantreg = TRUE, importance ="impurity")
importance_table = qrf_fit_test$variable.importance
test_that("Predictor ranking of dataset",{
  expect_gt(importance_table[2],importance_table[1])
  expect_gt(importance_table[1],importance_table[3])
})
test_that("Dataset complete?",{
  expect_equal(climset %>% arrange(Year, Month, Day), rbind(train_y, test_y) %>% arrange(Year, Month, Day))
  expect_equal(rbind(train_q,test_q) %>% arrange(Year, Month, Day),train_sub %>% arrange(Year, Month, Day))
})

test_that("Random subset numbers",{
  expect_equal(min(randomsubset),1)
  expect_equal(max(randomsubset),3)
})

