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
rm(list=ls(all=TRUE))
### SET GENERAL CONDITIONS FOR THE MODEL
#set thresholds and hyperparameter; determine test dataset and training dataset
p = 0.25 #power transformation to linearize thresholds
maxvars = 3
numsubset = 3 #number of subsets for hyperparameter selection
thres = c(4,6,10,15,25,36,50,70,100)
thres_eval = 20 #precipitation threshold for evaluation
minpredictant = 1.5 #minimum precipitation sum considered as precipitation
setwd("/usr/people/groote/")
ObsPV = read.csv(file = "Thunderstorm_radar_merged.csv")
setwd("/usr/people/groote/ThunderstormPostProcessingv1/")
years = c(as.numeric(unique(ObsPV$Year)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[1]
VT = unique(ObsPV$validtime.x)[2]
regions = c(unique(ObsPV$region))#[1:2]

#change radarmax into other name if necessary:
climset = filter(ObsPV, Ndischarge > minpredictant & validtime.x == VT & leadtime_count == LT)

#do transformations for thresholds
climset$Dischargerate = climset$Dischargerate^p
thres = thres^p
thres_eval = thres_eval^p
ndec = 4 #number of decimals usedi when appending scores to list of scores

#set available variables & predictant
#varindex=c(seq(18,71),seq(73,76),seq(80,82),seq(85,101))
varindex=c(seq(18,27),seq(30,37),seq(39,42),seq(44,57),seq(59,67),seq(70,71),seq(73,76),seq(81,82),seq(86,99))
pot_preds=names(climset[varindex])
ind_predictant = 113

##
#transformations Richardson numbers (otherwise this script for now protests; more transformations will be done in different script)
cnst_RItrans = 1e-5
RIindex=c(seq(36,37),seq(78,79))
climset[RIindex] = log(cnst_RItrans+climset[RIindex])
plot(climset[RIindex])
##

fit_test_all_pot_pred <- function(train_set, predictant, pot_pred_indices, train_thresholds, used_preds = c(0), used_sc_preds = c(0)){
  #selects the best predictor by forward fittng, trying potential predictors and selecting the one with lowest AIC
  AICscores = list()
  for(i in names(train_set[pot_pred_indices])){
    model = hxlr(reformulate(termlabels = paste0(reformulate(termlabels = names(data.frame(train_set[i], train_set[used_preds]))),"|",
                                                 reformulate(termlabels = names(data.frame(train_set[used_sc_preds]))))[2],
                                                 response = as.name(names(train_set[predictant]))),
                             data=data.frame(train_set[predictant], train_set[i], train_set[used_preds], train_set[used_sc_preds]),
                             thresholds = train_thresholds)
    AICscores = append(AICscores,AIC(model))
  }
  added = pot_pred_indices[unlist(AICscores[seq(length(pot_pred_indices))])==min(unlist(AICscores))]
  return(added)
}

fit_test_all_pot_pred_no_sc <- function(train_set, predictant, pot_pred_indices, train_thresholds, used_preds = c(0)){
  #selects the best predictor by forward fittng, trying potential predictors and selecting the one with lowest AIC
  AICscores = list()
  for(i in names(train_set[pot_pred_indices])){
    model = hxlr(reformulate(termlabels = paste0(reformulate(termlabels = names(data.frame(train_set[i], train_set[used_preds]))))[2],
                             response = as.name(names(train_set[predictant]))),
                 data=data.frame(train_set[predictant], train_set[i], train_set[used_preds]),
                 thresholds = train_thresholds)
    AICscores = append(AICscores,AIC(model))
  }
  added = pot_pred_indices[unlist(AICscores[seq(length(pot_pred_indices))])==min(unlist(AICscores))]
  return(added)
}
fit_test_all_pot_pred_sc <- function(train_set, predictant, pot_pred_indices, train_thresholds, used_preds = c(0), used_sc_preds = c(0)){
  #selects the best predictor by forward fittng, trying potential predictors and selecting the one with lowest AIC
  AICscores = list()
  for(i in names(train_set[pot_pred_indices])){
    model = hxlr(reformulate(termlabels = paste0(reformulate(termlabels = names(data.frame(train_set[used_preds]))),"|",
                                                 reformulate(termlabels = names(data.frame(train_set[i],train_set[used_sc_preds]))))[2],
                             response = as.name(names(train_set[predictant]))),
                 data=data.frame(train_set[predictant], train_set[i], train_set[used_preds], train_set[used_sc_preds]),
                 thresholds = train_thresholds)
    AICscores = append(AICscores,AIC(model))
  }
  added = pot_pred_indices[unlist(AICscores[seq(length(pot_pred_indices))])==min(unlist(AICscores))]
  return(added)
}
new_verify_ELR <- function(test_set, model, predictant, test_threshold, npred, reliabilityplot = FALSE){
  observed = data.frame(observation = as.numeric(test_set[predictant]>test_threshold))

  #predict with model and verify, calculate bs
  values <- as.numeric(1-predict(model, newdata = test_set, type = "cumprob", test_threshold))
  result = data.frame(observation = observed, probability = values, npredictors = npred, threshold = test_threshold)
  return(result)
}

fit_extended_logitModels <- function(train_set, test_set, predictant = ind_predictant, pot_pred_indices = varindex,
                                     train_thresholds = thres, test_thresholds = thres, maxnumbervars = maxvars){

  #initialize vector of chosen variables and brier scores
  variables = c()
  sc_variables = c()
  verified_set = data.frame()
  results = list()
  n_variables = c()

  if(maxnumbervars > length(pot_pred_indices)){
    message("You can not add more variables to the model than the number of available predictors")
    return("failed")
  }

  #select the best predictor and add to first model
  added = fit_test_all_pot_pred_no_sc(train_set, predictant, pot_pred_indices, train_thresholds)
  variables = c(variables, added)
  added = fit_test_all_pot_pred_sc(train_set, predictant, pot_pred_indices, train_thresholds, variables)
  sc_variables = c(sc_variables, added)

  #add this model to the model list
  firstmodel = hxlr(reformulate(termlabels = paste0(reformulate(termlabels = names(data.frame(train_set[variables]))),
                                                    "|",reformulate(termlabels = names(data.frame(train_set[sc_variables]))))[2],
                                response = as.name(names(train_set[predictant]))),
                    data=data.frame(train_set[predictant], train_set[variables], train_set[sc_variables]), thresholds = train_thresholds)
  modellist = list(firstmodel)
  n_variables = append(n_variables, length(variables))

  ### ITERATION, ADDING VARIABLES
  while(length(variables) < maxnumbervars){
    #update potential predictors remaining
    remaining_indices = pot_pred_indices[!pot_pred_indices%in%variables]

    #get best new predictor variable and add to variables
    added = fit_test_all_pot_pred(train_set, predictant, remaining_indices, train_thresholds, used_preds = variables, used_sc_preds = sc_variables)
    variables = c(variables, added)

# Switch following line on if multiple scale predictors have to be used!! Then there are always as many location as scale predictors. Otherwise, there is one scale predictor.
 #   added = fit_test_all_pot_pred_sc(train_set, predictant, pot_pred_indices, train_thresholds, used_preds = variables, used_sc_preds = sc_variables)
#    sc_variables = c(sc_variables, added)

    #add model based on these variables to model list
    bestmodel = hxlr(reformulate(termlabels = paste0(reformulate(names(data.frame(train_set[variables]))),"|",
                                                     reformulate(names(data.frame(train_set[sc_variables]))))[2],
                                 response = as.name(names(train_set[predictant]))),
                     data=data.frame(train_set[predictant], train_set[variables], train_set[sc_variables]), thresholds = train_thresholds)
    modellist = append(modellist, list(bestmodel))
    n_variables = append(n_variables, length(variables))
  }

  #select model and apply verification result with verification model per region
  for (i in seq(modellist)){
    model_ver = modellist[[i]]
    n_variables_i = n_variables[i]

    for(j in seq(length(test_thresholds))){
      test_threshold = test_thresholds[j]
      verified_set = rbind(verified_set, new_verify_ELR(test_set, model_ver, predictant, test_threshold, n_variables_i, reliabilityplot = FALSE))

    }
  }

  #add to results list
  results$models = modellist
  results$verification = verified_set
  return(results)
}

# -----------------------------------------------------------------------------

brierdataframe = data.frame()
models = list()

q = 1
for(y in years){
  train_fin = filter(climset, Year != y)
  set.seed(15+seq(years)[q]) #for reproducability purposes
  testdf = data.frame(validdate = unique(train_fin$validdate), subset = round(runif(unique(train_fin$validdate))*numsubset+0.5))
  train_sub <- left_join(train_fin, testdf, by = names(testdf)[1])
  #train_sub = cbind(train_fin, subset = randomsubset)
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
   # print(result$briers)

    #put results in dataframes and vectors
    brierdataframe = rbind(brierdataframe, result$verification)

    models = append(models, result$models)
  }
  q = q+1 #update seed set
}
models2 = list()
brierdataframe2 = data.frame()
for(y in years){
  test_fin = filter(climset, Year == y)
  train_fin = filter(climset, Year != y)
  result = fit_extended_logitModels(train_fin, test_fin, predictant = ind_predictant, pot_pred_indices = varindex,
                                    train_thresholds = thres, test_thresholds = thres_eval, maxnumbervars = maxvars)
  brierdataframe2 = rbind(brierdataframe2, data.frame(result$verification))
  models2 = append(models2, result$models)
}

thescores=data.frame()
for(npred in unique(brierdataframe2$npredictors)){
  for(thres in unique(brierdataframe2$threshold)){
    subset = filter(brierdataframe, npredictors == npred, threshold == thres)
    subset2 = filter(brierdataframe2, npredictors == npred, threshold == thres)
    score = brier(subset$observation, subset$probability, bins = FALSE)$ss
    score2 = brier(subset2$observation, subset2$probability, bins = FALSE)$ss
    thescores = rbind(thescores, data.frame(numpredictors = npred, threshold = thres, ss_9fold = score, ss_years = score2))
  }
}
plot(thescores$numpredictors, thescores$ss_years)
print(thescores)
write.csv(thescores, "HELR_scores2.csv")
