########################################
#This is without any transformations!!!#
#For now, sample climatology is used!!!#
########################################



#### BUILD IN SOME TESTING!!!!!

### BUILD IN RPS SCORE


#rm(list=ls(all=TRUE))
### LOAD EXTERNAL CODE
library(dplyr)
library(MASS)
library(verification)
library(devtools)
usethis::use_testthat()

### SET GENERAL CONDITIONS FOR THE MODEL
#read dataset
ObsPV = readRDS(file = "/usr/people/whan/Research/Whanetal_HarmoniePr_2017/data/ObsPV.rds")
years = c(as.numeric(unique(ObsPV$Year)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[1]
VT = c(unique(ObsPV$validtime))[2]
regions = c(unique(ObsPV$region))

#set default available variables: predictant and predictor
numsubset = 3 #number of subsets for hyperparametersetting
ind_predictant = 6
varindex=seq(16,99)
pot_preds=names(ObsPV[varindex])
ndec = 7
maxsteps = 6

fit_logitModels_and_predict <- function(train_set, test_set, predictant = ind_predictant, pot_pred_indices = varindex,
                         thres = 3.00, maxstepsAIC = 10, print_conf = FALSE){

  #do the reference fit and compute observations
  nullfit = glm(unlist(train_set[predictant] > thres) ~ 1, data = train_set[pot_pred_indices], family=binomial)
  observed = as.numeric(test_set[predictant] > thres)

  #results storage
  modellist = list()
  briers = c()
  results = list()

  #stepwise LR model fitting
  for (nstepsAIC in seq(maxstepsAIC)){
    logitMod = stepAIC(nullfit, scope = list(upper = lm(unlist(train_set[predictant] > thres) ~ .,
                                                        data=train_set[pot_pred_indices]),
                                                        lower = ~ 1), trace = 0, steps=nstepsAIC)
    if(print_conf == TRUE){
      conf=exp(confint.default(logitMod))
      print(conf)
    }

    #verification
    values <- predict(logitMod, newdata = test_set)
    probability = exp(values)/(1+exp(values))
    last_brier = brier(observed, probability, bins = FALSE)$bs #alternative with bins - - verification_set$bs
    verification_set <- verify(observed, probability, frcst.type = "prob", obs.type = "binary", title = )
    briers = append(briers,round(last_brier,ndec))


    #visualize skill
    reliability.plot(verification_set, titl = paste("Number of predictors = ", (nstepsAIC), " - Brier score = ", round(last_brier,ndec)), legend.names = paste(logitMod$formula[3]))

    #make lists of model, probabilities and brier scores including skill score
    modellist = append(modellist,list(logitMod))
  }

  #returned values
  results$model = modellist
  results$briers = briers
  results$nullfit = list(nullfit)
  return(results)
}

#create memory for models and their evaluation
brierdataframe = data.frame()
models = list()
nullfits = list()
for(y in years){

  #create random subsets to train on
  train = filter(ObsPV, Year != y & validtime == VT & leadtime_count == LT)
  randomsubset = round(runif(nrow(train))*numsubset+0.5)
  train_sub = cbind(train,subset = randomsubset)

  for(j in seq(numsubset)){
    #check approximately equal length of random subsets by printing relative length
    relweight_subset = sum(train_sub$subset[train_sub$subset == j])/nrow(train_sub)/j
    print(relweight_subset)
    print(j)

    #select training and testing dataset
    train_j = filter(train_sub, subset != j)[seq(length(ObsPV))]
    test_j = filter(train_sub, subset == j)[seq(length(ObsPV))]

    #find a fit
    result = fit_logitModels_and_predict(train = train_j, test = test_j, predictant = ind_predictant, pot_pred_indices = varindex, thres = 3.00, maxstepsAIC = maxsteps, print_conf = FALSE)

    #put results in dataframes and vectors
    brierdataframe = rbind(brierdataframe, data.frame(seq(maxsteps), rep(y, maxsteps), result$briers))
    models = append(models, result$model)
    nullfits = append(nullfits, result$nullfit)
  }
}

scores = matrix(c(brierdataframe$result.briers), nrow = maxsteps, ncol = (numsubset*length(years)))
meanscores = rowMeans(scores)
matplot(scores, type = "o", xlab = "Number of predictors (-)", ylab = "Brier score (-)")
plot(meanscores, xlab = "Number of predictors (-)", ylab = "Mean brier score (-)")

for(y in years){
  train_fin = filter(ObsPV, Year != y & validtime == VT & leadtime_count == LT)
  test_fin = filter(ObsPV, Year == y & validtime == VT & leadtime_count == LT)
  result = fit_logitModels_and_predict(train_set = train_fin, test_set = test_fin,
                                       predictant = ind_predictant, pot_pred_indices = varindex,
                                       thres = 3.00, maxstepsAIC = maxsteps, print_conf = FALSE)
  brierdataframe = rbind(brierdataframe, data.frame(seq(maxsteps), rep(y, maxsteps), result$briers))
  models = append(models, result$model)
  nullfits = append(nullfits, result$nullfit)
}

test_that("Test dataset complete?", {
  expect_equal(filter(ObsPV, validtime == VT & leadtime_count == LT), rbind(train_fin, test_fin))
  })
