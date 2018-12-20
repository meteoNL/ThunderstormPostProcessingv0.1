rm(list=ls(all=TRUE))
# load libraries:
# if they need to be installed that can be done by: install.packages("dplyr")
library(dplyr)
library(tidyr)
library(quantregForest)
library(verification)
library(ggplot2)
library(tree)

#hyperparameters
numbtree = 250
m_settings = c(2, 6, 18)
min_length = 2
node_size_settings = c(9, 27, 81)#50

# import the data frame and generate training and testing set
ObsPV = readRDS(file = "/usr/people/whan/Research/Whanetal_HarmoniePr_2017/data/ObsPV.rds")
#ObsPV <- filter(ObsPV, region != 1)
numbsubset = 3
years = c(as.numeric(unique(ObsPV$Year)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[1]
VT = c(unique(ObsPV$validtime))[2]
regions = c(unique(ObsPV$region))
th = exp(seq(0,5)/3)
climset = filter(ObsPV, validtime == VT & leadtime_count == LT)
orig_varindex = seq(16,25)#seq(16, length(climset))
predictant_ind = 6

#verify_qrf_per_reg <- function(fit, test_set, predictant, reg_set, npredictors, numbertrees, mval, nsize){
# for(reg in reg_set){
#    observed_subset = filter(test_set, region = reg)
#    observed = observed_subset[predictant]
#    qrf_pred <- predict(qrf_fit, data.frame(test_q[, pot_preds]), what = ecdf)
#    brier_score = brier(as.numeric(qrf_pred > thres), )
#  }
#}
#varindex = seq(16,length(climset)) #indices of potential predictors

# ## in this data frames we have:
# region = 2:12
# Year = 2010, 2011, 2013
# Month = 04:10
# Day = 1:31
# validtime = "VT_0612" "VT_1218" "VT_1800" "VT_0006"
# radarmax = the radar observation
# validdate
# T0 = initialization time
# leadtime = how many hours lead time
# leadtime_count 1:7
# init = initialization of the model
# vt_start = start of the verification time
# vt_end = end of the verification time
# LT_start = start of the lead time
# LT_end = end of the lead time
# Columns 17-99 are the potential predictor variables. maximum and minimum values of each predictor.
remove_variable = c()
overall_scores = data.frame()


q = 1
for (y in years){
  train_y = filter(climset, Year != y)
  test_y = filter(climset, Year == y)
  set.seed(seq(years)[q])
  randomsubset = round(runif(nrow(train_y))*numbsubset+0.5)
  train_sub = cbind(train_y,subset = randomsubset)
  j = 1
  for (w in seq(numbsubset)){
    print(y)
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
        while (length(train_q[varindex])>min_length){
          #predictors still in
          pot_preds <- names(train_q[varindex])
          #do a fit
          qrf_fit <- quantregForest(x = data.frame(train_q[varindex]), y = unlist(train_q[predictant_ind]), ntree=numbtree, mtry = m, nodesize = node_size)
          plot(qrf_fit)

          #remove variable
          remove_variable = varindex[qrf_fit$importance == min(qrf_fit$importance)][1]
          varindex = varindex[varindex != remove_variable]

          #make prediction for thresholds
          qrf_pred <- predict(qrf_fit, data.frame(test_q[, pot_preds]), what = ecdf)
          qrf_probs <- do.call(rbind, lapply(th, function(t) {
            data.frame(th = t, prob = unlist(lapply(qrf_pred, function(nr) 1-nr(t))), binobs = unlist(test_q[predictant_ind]) > t, reg = test_q["region"])
          }))

          #evaluate
          qrf_bs <- qrf_probs %>% group_by(th, region) %>% summarise(bs = brier(obs = binobs, pred = prob, bins = FALSE)$bs)
          print(length(varindex))
          print(qrf_bs)
          brierscore = data.frame(matrix(qrf_bs$bs, ncol = length(th)), region = matrix(qrf_bs$region, ncol = length(th)), testsubset = w, nodesize = node_size, mtry = m, test_year = y, npredictors = length(varindex))
          overall_scores = rbind(overall_scores, brierscore)
          j = j+1
          print(q)
          print(j)
        }
      }
    }
  }
  q = q + 1
}
write.csv(overall_scores, file = "overall scores.csv")
plot(data.frame(overall_scores$X1, overall_scores$nodesize))
plot(data.frame(overall_scores$X3, overall_scores$nodesize))
plot(data.frame(overall_scores$X5, overall_scores$nodesize))

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
write.csv(quantiles_frame, file = "quantiles frame.csv")

set.seed(712)
x1 = rnorm(1000,0,5)
x2 = rnorm(1000,0,5)
x3 = rnorm(1000,0,5)
pert = rnorm(1000,0,5)


#-----------------------------------------------------------------
## Testing the functions
library(devtools)
library(testthat)
usethis::use_testthat()

# fit a QRF that predicts radar data from a set of potential predictors
#pot_preds <- names(train[varindex])
#qrf_fit <- quantregForest(x = data.frame(train[, pot_preds]), y = unlist(train[, "radarmax"]),ntree=numbtree, mtry = m)

# predict the ECDF for each day in the test data set
#qrf_pred <- predict(qrf_fit, data.frame(test[, pot_preds]), what = ecdf)
#plot(qrf_fit)
# make probabilitys of exceeding a set of thresholds:
#th = seq(1, 15, 2)
#qrf_probs <- do.call(rbind, lapply(th, function(t) {
#  data.frame(th = t, prob = unlist(lapply(qrf_pred, function(nr) 1-nr(t))), binobs = unlist(test[, "radarmax"]) > t)
#}))
# calculate the BS for each thrshold
#qrf_bs <- qrf_probs %>% group_by(th) %>% summarise(bs = brier(obs = binobs, pred = prob, bins = FALSE)$bs)

# make a plot
#print(ggplot(data = qrf_bs) +
#       geom_point(aes(x = th, y = bs)) +
#        theme_bw() + ylim(0, 0.2) +
#         xlab("threshold"))
# dev.off()

#for(i in seq(1,length(th))){
# plot(data.frame(overall_scores[1], overall_scores[i+1]), type = "o", main = paste("Threshold in mm/h: ", th[i]))
#}

#show splits
#alltrees = data.frame(getTree(qrf_fit,1))
#for(i in seq(2,numbtree)){
#  fulltree = data.frame(getTree(qrf_fit, i))
#  alltrees = rbind(alltrees, fulltree)
#}
#varnames = names(qrf_fit$forest$xlevels)
#for (i in seq(1,length(varnames))){
#  variable_splits <- filter(alltrees,split.var == i)
#  plot(variable_splits$split.point, variable_splits$prediction, main = varnames[i])
#}
#fullt1 = getTree(qrf_fit, 1)
#t1 = data.frame(fullt1)
#tree1 <- tree(data = t1)
#plot(tree1)

