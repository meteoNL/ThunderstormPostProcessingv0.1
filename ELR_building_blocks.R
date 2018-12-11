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

rm(list=ls(all=TRUE))
### SET GENERAL CONDITIONS FOR THE MODEL
#set thresholds and hyperparameter; determine test dataset and training dataset
p = 0.001 #power transformation to linearize thresholds
maxvars = 4
numsubset = 3 #number of subsets for hyperparametersetting
thres = (0.1*exp(seq(0,9)/3))
thres_eval = 10.0 * thres
ObsPV = readRDS(file = "/usr/people/whan/Research/Whanetal_HarmoniePr_2017/data/ObsPV.rds")
years = c(as.numeric(unique(ObsPV$Year)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[1]
VT = c(unique(ObsPV$validtime))[2]
regions = c(unique(ObsPV$region))#[1:2]
climset = filter(ObsPV, radarmax > 0 & validtime == VT & leadtime_count == LT)

#do transformations for threshold
climset$radarmax = climset$radarmax^p
thres = thres^p
thres_eval = thres_eval^p
ndec = 4 #number of decimals usedi when appending scores to list of scores

#set available variables
varindex=seq(16,99)
pot_preds=names(climset[varindex])
ind_predictant = 6

##
#transformations Richardson numbers
cnst_RItrans = 1e-5
RIindex=c(seq(34,35),seq(76,77))
climset[RIindex] = log(cnst_RItrans+climset[RIindex])
plot(climset[RIindex])
##

fit_test_all_pot_pred <- function(train_set, predictant, pot_pred_indices, train_thresholds, used_preds = c(0)){
  #selects the best predictor by forward fittng, trying all remaining possible predictors and selecting the one with lowest AIC
  AICscores = list()
  for(i in names(train_set[pot_pred_indices])){
    model = hxlr(reformulate(termlabels = names(data.frame(train_set[i], train_set[used_preds])), response = as.name(names(train_set[predictant]))), data=data.frame(train_set[predictant], train_set[i], train_set[used_preds]), thresholds = train_thresholds)
    AICscores = append(AICscores,AIC(model))
  }
  added = pot_pred_indices[unlist(AICscores[seq(length(pot_pred_indices))])==min(unlist(AICscores))]
  # print(AICscores)
  return(added)
}

verify_model_per_reg <- function(test_set, model, reg_set, predictant, test_thresholds, i, reliabilityplot = FALSE){
  briers = c()
  for(reg in reg_set){

    #select subset
    region_subset = filter(test_set, region == reg)

    #get observations for this region
    observed = data.frame(matrix(rep(as.matrix(region_subset[predictant]),length(test_thresholds)),ncol=length(test_thresholds)))

    #May make one more separate function from what is below

    #predict with model and verify, calculate bs
    values <- predict(model, newdata = region_subset, type = "cumprob", thresholds = test_thresholds)
    verification_set <- verify(as.numeric(observed < test_thresholds), values, frcst.type = "prob", obs.type = "binary", title = "")
    brierval = brier(as.numeric(observed < test_thresholds), values, bins = FALSE)$bs
    briers = append(briers, c(i, reg, y, brierval))

    #if required plot reliability plot
    if (reliabilityplot == TRUE){
      reliability.plot(verification_set, titl = paste(paste(names(model_ver$start)[seq(3,length(model_ver$start))],collapse=" + "), " - Brier score = ", round(verification_set$bs,ndec)))
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

  #select the best predictor and add to first model
  added = fit_test_all_pot_pred(train_set, predictant, pot_pred_indices, train_thresholds)
  variables = c(variables, added)

  #add this model to the model list
  firstmodel = hxlr(reformulate(termlabels = names(data.frame(train_set[variables])), response = as.name(names(train_set[predictant]))), data=data.frame(train_set[predictant], train_set[variables]), thresholds = train_thresholds)
  modellist = list(firstmodel)

  ### ITERATION, ADDING VARIABLES
  while(length(variables) < maxvars){
    #update potential predictors remaining
    remaining_indices = pot_pred_indices[!pot_pred_indices%in%variables]

    #get best new predictor variable and add to variables
    added = fit_test_all_pot_pred(train_set, predictant, remaining_indices, train_thresholds, used_preds = variables)
    variables = c(variables, added)

    #add model based on these variables to model list
    bestmodel = hxlr(reformulate(termlabels = names(data.frame(train_set[variables])), response = as.name(names(train_set[predictant]))), data=data.frame(train_set[predictant], train_set[variables]), thresholds = train_thresholds)
    modellist = append(modellist, list(bestmodel))
  }

  #select model and apply verification result with verification model per region function
  for (i in seq(modellist)){
    model_ver = modellist[[i]]
    briers = append(briers, verify_model_per_reg(test_set, model_ver, regions, predictant, test_thresholds, length(variables), reliabilityplot = FALSE))
  }

  #add to results list
  results$models = modellist
  results$brierscores = briers
  return(results)
}

# -----------------------------------------------------------------------------

brierdataframe = data.frame()
models = list()

for(y in years){
  train_fin = filter(climset, Year != y)
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
                                      train_thresholds = thres, test_thresholds = thres, maxnumbervars = maxvars)

    #put results in dataframes and vectors
    brierdataframe = data.frame(result$brierscores) #rbind(brierdataframe, data.frame(test_year = (result$briers)[seq(1,length(result$briers),4)],
    #                                npredictors = (result$briers)[seq(2,length(result$briers),4)],
    #                                region = (result$briers)[seq(3,length(result$briers),4)],
    #                                brier_score = (result$briers)[seq(4,length(result$briers),4)]))

    models = append(models, result$models)
  }
}

for(y in years){
  test_fin = filter(climset, Year == y)
  result = fit_extended_logitModels(train_fin, test_fin, predictant = ind_predictant, pot_pred_indices = varindex,
                                    train_thresholds = thres, test_thresholds = thres, maxnumbervars = maxvars)
  brierdataframe = data.frame(result$brierscores)#rbind(brierdataframe, data.frame(test_year = (result$briers)[seq(1,length(result$briers),4)],
  #                                 npredictors = (result$briers)[seq(2,length(result$briers),4)],
  #                                region = (result$briers)[seq(3,length(result$briers),4)],
  #                               brier_score = (result$briers)[seq(4,length(result$briers),4)]))
  models = append(models, result$models)
}
