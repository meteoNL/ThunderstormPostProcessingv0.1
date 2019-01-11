rm(list=ls(all=TRUE))
# load libraries:
# if they need to be installed that can be done by: install.packages("dplyr")
library(dplyr)
library(tidyr)
library(quantregForest)
library(verification)
library(ggplot2)
library(tree)
library(SpecsVerification)

#hyperparameters
numbtree = 250
m_settings = c(2, 6, 18)
min_length = 2
node_size_settings = c(3, 9, 27)#50

# import the data frame and generate training and testing set
minpdis = 1.5
ObsPV = readRDS(file = "/usr/people/whan/Research/Whanetal_HarmoniePr_2017/data/ObsPV.rds")
#ObsPV <- filter(ObsPV, region != 1)
numbsubset = 3
years = c(as.numeric(unique(ObsPV$Year)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[1]
VT = c(unique(ObsPV$validtime))[2]
regions = c(unique(ObsPV$region))
th = exp(seq(0,5)/3)
climset = filter(ObsPV, Ndischarges > mindis & validtime == VT & leadtime_count == LT)
orig_varindex = seq(16, length(climset))
predictant_ind = 6

remove_variable = c()
overall_scores = data.frame()

qrf_procedure <- function(train_set, test_set, predictant_index, varindexset, m_hyp, ntree_hyp, node_size_hyp, min_length_setting, wval, yval){
  overall_scores_local = data.frame()
  while (length(train_set[varindexset])>min_length_setting){
    #predictors still in
    pot_preds <- names(train_set[varindexset])
    #do a fit
    qrf_fit <- quantregForest(x = data.frame(train_set[varindexset]), y = unlist(train_set[predictant_index]), ntree = ntree_hyp, mtry = m_hyp, nodesize = node_size_hyp)
    #plot(qrf_fit)

    #remove variable
    remove_variable = varindexset[qrf_fit$importance == min(qrf_fit$importance)][1]
    varindexset = varindexset[varindexset != remove_variable]

    #make prediction for thresholds
    qrf_pred <- predict(qrf_fit, data.frame(test_set[, pot_preds]), what = ecdf)
    qrf_probs <- do.call(rbind, lapply(th, function(t) {
      data.frame(th = t, prob = unlist(lapply(qrf_pred, function(nr) 1-nr(t))), binobs = unlist(test_set[predictant_index]) > t, reg = test_set["region"])
    }))

    #evaluate continuous ranked probability score
    qrf_pred_quan = predict(qrf_fit, data.frame(test_set[, pot_preds]), what = seq(0.05,0.95,0.1))
    Crps_samples = EnsCrps(qrf_pred_quan, unlist(test_set[predictant_index]))
    df_with_crpsscore = data.frame(test_set, crpsscore = Crps_samples)
    CRPS_score_mean <- df_with_crpsscore %>% group_by(region) %>% summarise(mean(crpsscore))
    print(CRPS_score_mean)

    #evaluate brier score and add to data frame
    qrf_bs <- qrf_probs %>% group_by(th, region) %>% summarise(bs = brier(obs = binobs, pred = prob, bins = FALSE)$bs)
    print(length(varindexset))
    print(qrf_bs)
    brierscore = data.frame(matrix(qrf_bs$bs, ncol = length(th)), region = matrix(qrf_bs$region, ncol = length(th)), crpsscore = CRPS_score_mean[2], testsubset = wval, nodesize = node_size_hyp, mtry = m_hyp, test_year = yval, npredictors = length(varindexset))
    overall_scores_local = rbind(overall_scores_local, brierscore)
    print(q)
  }
  return(overall_scores_local)
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
    print(train_q)
    for (node_size in node_size_settings){
      print("node size equals: ")
      print(node_size)
      for (m in m_settings){
        print("m equals: ")
        print(m)
        varindex = orig_varindex
        brierscore = qrf_procedure(train_q, test_q, predictant_ind, varindex, m, numbtree, node_size, min_length, w, y)
        overall_scores = rbind(overall_scores, brierscore)

      }
    }
  }
  q = q + 1
}
#this table contains all brier scores for all hyperparameters (mtry, node_size, npredictors), per region, per threshold and per testsubset
write.csv(overall_scores, file = "overall scores.csv")
#plot(data.frame(overall_scores$X1, overall_scores$nodesize))
#plot(data.frame(overall_scores$X3, overall_scores$nodesize))
#plot(data.frame(overall_scores$X5, overall_scores$nodesize))

quantiles_frame = data.frame()
for (m in m_settings){
  for (node_size in node_size_settings){
    for (threshold in seq(th)){
      for (npredictorsval in unique(c(overall_scores$npredictors))){
        print("threshold index equals:")
        print(threshold)
        print("m equals:")
        print(m)
        print("node size equals:")
        print(node_size)
        print("number of predictors equals:")
        print(npredictorsval)
        subset <- filter(overall_scores, mtry == m & nodesize == node_size & npredictors == npredictorsval)
        retrieve_q = subset[threshold]
        retrieve_q = unlist(c(retrieve_q))
        values = quantile(retrieve_q, c(0.1, 0.25, 0.5, 0.75, 0.9))
        print(data.frame(values, m, th[threshold], node_size))
        quantiles_frame = rbind(data.frame(q0.10 = values[1], q0.25 = values[2], q0.50 = values[3], q0.75 = values[4], q0.90 = values[5], m, th[threshold], node_size, npredictorsval), quantiles_frame)
      }
    }
  }
}

quantiles_frame_crps = data.frame()
for (m in m_settings){
  for (node_size in node_size_settings){
    for (npredictorsval in unique(c(overall_scores$npredictors))){
      print("m equals:")
      print(m)
      print("node size equals:")
      print(node_size)
      print("number of predictors equals:")
      print(npredictorsval)
      subset <- filter(overall_scores, mtry == m & nodesize == node_size & npredictors == npredictorsval)
      retrieve_q = subset$mean.crpsscore.
      retrieve_q = unlist(c(retrieve_q))
      values = c(mean(retrieve_q), min(retrieve_q), max(retrieve_q))
      print(data.frame(values, m, node_size))
      quantiles_frame_crps = rbind(data.frame(mean = values[1], min = values[2], max = values[3], m, node_size, npredictorsval), quantiles_frame_crps)
    }
  }
}

#this dataframe contains all brier scores per hyperparameter setting and per threshold as quantiles
write.csv(quantiles_frame, file = "quantiles frame.csv")
write.csv(quantiles_frame_crps, file = "quantiles frame with crps.csv")

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
qrf_fit_test <- quantregForest(x = data.frame(test_df[seq(1,3)]), y = unlist(test_df["yvalues"]), ntree=numbtree, mtry = 1, nodesize = 5)
importance_table = qrf_fit_test$importance
test_that("Predictor ranking of dataset",{
  expect_gt(importance_table[2],importance_table[1])
  expect_gt(importance_table[1],importance_table[3])
})
test_that("Dataset complete?",{
  expect_equal(rbind(test_y,train_y),climset)
  expect_equal(rbind(train_q,test_q),train_sub)
})

test_that("Random subset numbers",{
  expect_equal(min(randomsubset),1)
  expect_equal(max(randomsubset),3)
})

