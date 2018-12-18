########################################
#This is without any transformations!!!#
#For now, sample climatology is used!!!#
########################################


#######


#THIS IS STILL USING COMPLEMENTARY PROBABILITIES P_y=0 = (1-P_y=1)

#######

### LOAD EXTERNAL CODE
library(dplyr)
library(MASS)
library(verification)
library(arm)
library(crch)
library(profvis)
#rm(list=ls(all=TRUE))
### SET GENERAL CONDITIONS FOR THE MODEL
#set thresholds and hyperparameter; determine test dataset and training dataset
p = 0.1 #power transformation to linearize thresholds
maxvars = 4
numsubset = 3 #number of subsets for hyperparameter selection
thres = (0.3*exp(seq(0,2)/2)) #precipitation thresholds used for "g(q)"
thres_eval = 15.0 * thres[1] #precipitation threshold for evaluation
minpredictant = 0.1 #minimum precipitation sum considered as precipitation
ObsPV = readRDS(file = "/usr/people/whan/Research/Whanetal_HarmoniePr_2017/data/ObsPV.rds")
years = c(as.numeric(unique(ObsPV$Year)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[1]
VT = c(unique(ObsPV$validtime))[2]
regions = c(unique(ObsPV$region))#[1:2]
climset = filter(ObsPV, radarmax > minpredictant & validtime == VT & leadtime_count == LT)

#do transformations for thresholds
climset$radarmax = climset$radarmax^p
thres = thres^p
thres_eval = thres_eval^p
ndec = 4 #number of decimals usedi when appending scores to list of scores

#set available variables & predictant
varindex=seq(16,99)
pot_preds=names(climset[varindex])
ind_predictant = 6

##
#transformations Richardson numbers (otherwise this script for now protests; more transformations will be done in different script)
cnst_RItrans = 1e-5
RIindex=c(seq(34,35),seq(76,77))
climset[RIindex] = log(cnst_RItrans+climset[RIindex])
plot(climset[RIindex])
##

fit_test_all_pot_pred <- function(train_set, predictant, pot_pred_indices, train_thresholds, used_preds = c(0)){
  #selects the best predictor by forward fittng, trying potential predictors and selecting the one with lowest AIC
  AICscores = list()
  for(i in names(train_set[pot_pred_indices])){
    model = hxlr(reformulate(termlabels = names(data.frame(train_set[i], train_set[used_preds])), response = as.name(names(train_set[predictant]))), data=data.frame(train_set[predictant], train_set[i], train_set[used_preds]), thresholds = train_thresholds)
    AICscores = append(AICscores,AIC(model))
  }
  added = pot_pred_indices[unlist(AICscores[seq(length(pot_pred_indices))])==min(unlist(AICscores))]
  return(added)
}

verify_ELRmodel_per_reg <- function(test_set, model, reg_set, predictant, test_threshold, i, reliabilityplot = FALSE){
  briers = c()
  for(reg in reg_set){

    #select subset
    region_subset = filter(test_set, region == reg)
    observed = data.frame(region_subset[predictant])

    #predict with model and verify, calculate bs
    values <- predict(model, newdata = region_subset, type = "cumprob", thresholds = test_threshold)
    verification_set <- verify(as.numeric(observed > test_threshold), (1-values), frcst.type = "prob", obs.type = "binary", title = "")
    eval = brier(as.numeric(observed > test_threshold), (1-values), bins = FALSE)
    brierval = eval$bs
    brierbase = eval$bs.baseline
    briers = append(briers, c(y, i, reg, test_threshold, brierval, brierbase))

    #if required plot reliability plot
    if (reliabilityplot == TRUE){
      reliability.plot(verification_set, titl = paste(paste(names(model$start)[seq(3,length(model$start))],collapse=" + "), " - Brier score = ", round(verification_set$bs,ndec)))
    }
  }
  return(briers)
}


fit_extended_logitModels <- function(train_set, test_set, predictant = ind_predictant, pot_pred_indices = varindex,
                                     train_thresholds = thres, test_thresholds = thres, maxnumbervars = maxvars){

  #initialize vector of chosen variables and brier scores
  variables = c()
  briers = c()
  results = list()
  n_variables = c()

  if(maxnumbervars > length(pot_pred_indices)){
    message("You can not add more variables to the model than the number of available predictors")
    return("failed")
  }

  #select the best predictor and add to first model
  added = fit_test_all_pot_pred(train_set, predictant, pot_pred_indices, train_thresholds)
  variables = c(variables, added)

  #add this model to the model list
  firstmodel = hxlr(reformulate(termlabels = names(data.frame(train_set[variables])), response = as.name(names(train_set[predictant]))), data=data.frame(train_set[predictant], train_set[variables]), thresholds = train_thresholds)
  modellist = list(firstmodel)
  n_variables = append(n_variables, length(variables))

  ### ITERATION, ADDING VARIABLES
  while(length(variables) < maxnumbervars){
    #update potential predictors remaining
    remaining_indices = pot_pred_indices[!pot_pred_indices%in%variables]

    #get best new predictor variable and add to variables
    added = fit_test_all_pot_pred(train_set, predictant, remaining_indices, train_thresholds, used_preds = variables)
    variables = c(variables, added)

    #add model based on these variables to model list
    bestmodel = hxlr(reformulate(termlabels = names(data.frame(train_set[variables])), response = as.name(names(train_set[predictant]))), data=data.frame(train_set[predictant], train_set[variables]), thresholds = train_thresholds)
    modellist = append(modellist, list(bestmodel))
    n_variables = append(n_variables, length(variables))
  }

  #select model and apply verification result with verification model per region function
  for (i in seq(modellist)){
    model_ver = modellist[[i]]
    n_variables_i = n_variables[i]

    for(j in seq(test_thresholds)){
      test_threshold = test_thresholds[j]
      briers = append(briers, verify_ELRmodel_per_reg(test_set, model_ver, regions, predictant, test_threshold, n_variables_i, reliabilityplot = TRUE))

    }
  }

  #add to results list
  results$models = modellist
  results$brierscores = briers
  return(results)
}

# -----------------------------------------------------------------------------

brierdataframe = data.frame()
models = list()

q = 1
for(y in years){
  train_fin = filter(climset, Year != y)
  set.seed(seq(years)[q]) #for reproducability purposes
  randomsubset = round(runif(nrow(train_fin))*numsubset+0.5)
  train_sub = cbind(train_fin, subset = randomsubset)
  for(j in seq(numsubset)){
    #check approximately equal length of random subsets by printing relative length
    relweight_subset = sum(train_sub$subset[train_sub$subset == j])/nrow(train_sub)/j
    print(relweight_subset)
    print(j)

    #select training and testing dataset
    train_j = filter(train_sub, subset != j)[seq(length(climset))]
    test_j = filter(train_sub, subset == j)[seq(length(climset))]

    #find a fit
    result = fit_extended_logitModels(train_j, test_j, predictant = ind_predictant, pot_pred_indices = varindex,
                                      train_thresholds = thres, test_thresholds = thres_eval, maxnumbervars = maxvars)

    #put results in dataframes and vectors
    brierdataframe = rbind(brierdataframe, data.frame(test_year = (result$briers)[seq(1,length(result$briers),6)],
                                                      npredictors = (result$briers)[seq(2,length(result$briers),6)],
                                                      region = (result$briers)[seq(3,length(result$briers),6)],
                                                      test_threshold = (result$briers)[seq(4,length(result$briers),6)],
                                                      brier_score = (result$briers)[seq(5,length(result$briers),6)],
                                                      brier_base = (result$briers)[seq(6,length(result$briers),6)]))

    models = append(models, result$models)
  }
  q = q+1 #update seed set
}

for(y in years){
  test_fin = filter(climset, Year == y)
  result = fit_extended_logitModels(train_fin, test_fin, predictant = ind_predictant, pot_pred_indices = varindex,
                                    train_thresholds = thres, test_thresholds = thres_eval, maxnumbervars = maxvars)
  brierdataframe = rbind(brierdataframe, data.frame(test_year = (result$briers)[seq(1,length(result$briers),6)],
                                                    npredictors = (result$briers)[seq(2,length(result$briers),6)],
                                                    region = (result$briers)[seq(3,length(result$briers),6)],
                                                    test_threshold = (result$briers)[seq(4,length(result$briers),6)],
                                                    brier_score = (result$briers)[seq(5,length(result$briers),6)],
                                                    brier_base = (result$briers)[seq(6,length(result$briers),6)]))
  models = append(models, result$models)
}

plot(brierdataframe$npredictors, brierdataframe$brier_score)

#-----------------------------------------------------------------
## Testing the functions
library(devtools)
library(testthat)
usethis::use_testthat()
test_that("Test dataset and potential predictors complete?", {
  expect_equal(climset, rbind(train_fin, test_fin))
  expect_equal(length(pot_preds), length(varindex))
})
test_that("Transformation with power p - check", {
  expect_equal(c(climset[ind_predictant]^(1/p)),c(filter(ObsPV, radarmax > minpredictant & validtime == VT &
                                                  leadtime_count == LT)[ind_predictant]))
  expect_gte(min(climset[ind_predictant]^(1/p)),minpredictant)
})

y = years[1]
regions = c(1)
testthat_df = data.frame(a=(seq(10,20)+2*rnorm(11)),b=seq(20,40,2),d=rnorm(11), region = rep(1,11))
thresholds_testthat = c(quantile(testthat_df$a,0.25)[[1]],quantile(testthat_df$a,0.95)[[1]])
model_testthat = fit_extended_logitModels(train_set = testthat_df, test_set = testthat_df, predictant = 1, pot_pred_indices = c(2,3), train_thresholds = thresholds_testthat, test_thresholds = thresholds_testthat, maxnumbervars  = 1)$models
test_that("Testing function fit_test_all_pot_pred",{
  expect_equal(fit_test_all_pot_pred(train_j, 6, 28, thres, used_preds = 30), 28)
  expect_error(fit_test_all_pot_pred(train_j, 6, 28, thres, used_preds = 28))
  expect_error(fit_test_all_pot_pred(train_j, 28, 28, thres, used_preds = 30))
  expect_error(fit_test_all_pot_pred(train_j, 450, 28, thres, used_preds = 30))
})

test_that("Testing function verify_ELRmodel_per_reg: Brier scores of probabilities between 0.0 and 0.5 (0.5 implies a highly disappointing model as
          typically high probabilities of 50-100% will then be associated with no occurence and the reverse, 0-50% with occurences",{
  expect_gt(verify_ELRmodel_per_reg(testthat_df, model_testthat[[1]], regions, 1, thresholds_testthat[1], 1, reliabilityplot = FALSE)[5],0)
  expect_gt(verify_ELRmodel_per_reg(testthat_df, model_testthat[[1]], regions, 1, thresholds_testthat[1], 1, reliabilityplot = FALSE)[6],0)
  expect_lt(verify_ELRmodel_per_reg(testthat_df, model_testthat[[1]], regions, 1, thresholds_testthat[1], 1, reliabilityplot = FALSE)[5],0.5)
  expect_lt(verify_ELRmodel_per_reg(testthat_df, model_testthat[[1]], regions, 1, thresholds_testthat[1], 1, reliabilityplot = FALSE)[6],0.5)
  expect_error(verify_ELRmodel_per_reg(train_j, model_testthat[[1]], regions, 1, thresholds_testthat[1], 1, reliabilityplot = FALSE)[6])
})

test_that("Testing function fit_extended_logitModels compared to verification",{
  expect_equal(class(model_testthat[[1]]),"hxlr")
  expect_equal(fit_extended_logitModels(testthat_df, testthat_df, predictant = 1, pot_pred_indices = c(2,3), train_thresholds = thresholds_testthat, test_thresholds = thresholds_testthat, maxnumbervars = 1)$brierscores[1:6],verify_ELRmodel_per_reg(testthat_df, model_testthat[[1]], regions, 1, thresholds_testthat[1], 1, reliabilityplot = FALSE))
  expect_error(fit_extended_logitModels(testthat_df, testthat_df, predictant = 1, pot_pred_indices = c(2,7), train_thresholds = thresholds_testthat, test_thresholds = thresholds_testthat, maxnumbervars = 1))
  expect_error(fit_extended_logitModels(testthat_df, train_j, predictant = 1, pot_pred_indices = c(2,3), train_thresholds = thresholds_testthat, test_thresholds = thresholds_testthat, maxnumbervars = 1))
  expect_equal(fit_extended_logitModels(testthat_df, testthat_df, predictant = 1, pot_pred_indices = c(2,3), train_thresholds = thresholds_testthat, test_thresholds = thresholds_testthat, maxnumbervars = 3),"failed")
})

set.seed(121)
x1 = rnorm(500,0,5)
x2 = rnorm(500,0,5)
x3 = rnorm(500,0,5)
x4 = rnorm(500,0,5)
pert = rnorm(500,0,10)
y1 = x1*100 + x4*10 + pert
testthat_df2 = data.frame(y1, x1, x2, x3, x4, region = rep(1,500))
thresholds_testthat2 = c(quantile(testthat_df2$y1,0.25)[[1]],quantile(testthat_df2$y1,0.75)[[1]])
model_testthat2 = fit_extended_logitModels(train_set = testthat_df2, test_set = testthat_df2, predictant = 1, pot_pred_indices = seq(2,5), train_thresholds = thresholds_testthat2, test_thresholds = thresholds_testthat2, maxnumbervars  = 1)$models
model_testthat22 = fit_extended_logitModels(train_set = testthat_df2, test_set = testthat_df2, predictant = 1, pot_pred_indices = seq(2,5), train_thresholds = thresholds_testthat2, test_thresholds = thresholds_testthat2, maxnumbervars  = 2)$models

test_that("Testing the predictor choice among 4 predictors",{
  expect_equal(fit_test_all_pot_pred(testthat_df2, 1, seq(2,5), thresholds_testthat2), 2)
  expect_equal(fit_test_all_pot_pred(testthat_df2, 1, seq(3,5), thresholds_testthat2, used_preds = 2), 5)
  expect_equal(fit_extended_logitModels(testthat_df2, testthat_df2, predictant = 1, pot_pred_indices = c(2,5), train_thresholds = thresholds_testthat2, maxnumbervars = 2),fit_extended_logitModels(testthat_df2, testthat_df2, predictant = 1, pot_pred_indices = seq(2,5), train_thresholds = thresholds_testthat2, maxnumbervars = 2))
  expect_equal(model_testthat2[[1]],model_testthat22[[1]])
})

modelres_manually = exp(0.07722407+0.02944015*thresholds_testthat2[1]-2.94777089*testthat_df2$x1)
modelres2_manually = exp((-1.0189620+0.2509291*thresholds_testthat2[1]-25.5650782*x1-2.8128970*x4))
test_that("Verify function test with the second testing dataframe",{
  expect_equal(verify_ELRmodel_per_reg(testthat_df2, model_testthat22[[1]], regions, 1, thresholds_testthat2[1], 1, reliabilityplot = FALSE)[5],brier(as.numeric(testthat_df2$y1 > thresholds_testthat2[1]),1/(1+modelres_manually), bins = FALSE)$bs)
  expect_equal(verify_ELRmodel_per_reg(testthat_df2, model_testthat22[[1]], regions, 1, thresholds_testthat2[1], 1, reliabilityplot = FALSE)[6],brier(as.numeric(testthat_df2$y1 > thresholds_testthat2[1]),1/(1+modelres_manually), bins = FALSE)$bs.baseline)
  expect_equal(verify_ELRmodel_per_reg(testthat_df2, model_testthat22[[2]], regions, 1, thresholds_testthat2[1], 1, reliabilityplot = FALSE)[5],brier(as.numeric(testthat_df2$y1 > thresholds_testthat2[1]),1/(1+modelres2_manually), bins = FALSE)$bs)
  expect_equal(verify_ELRmodel_per_reg(testthat_df2, model_testthat22[[2]], regions, 1, thresholds_testthat2[1], 1, reliabilityplot = FALSE)[6],brier(as.numeric(testthat_df2$y1 > thresholds_testthat2[1]),1/(1+modelres2_manually), bins = FALSE)$bs.baseline)
})
