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
p = 0.5 #power transformation to linearize thresholds
maxvars = 4
numsubset = 3 #number of subsets for hyperparametersetting
thres = (0.3*exp(seq(0,2)/2))
thres_eval = 10.0 * thres
ObsPV = readRDS(file = "/usr/people/whan/Research/Whanetal_HarmoniePr_2017/data/ObsPV.rds")
years = c(as.numeric(unique(ObsPV$Year)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[1]
VT = c(unique(ObsPV$validtime))[2]
regions = c(unique(ObsPV$region))
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

fit_extended_logitModels <- function(train_set, test_set, predictant = ind_predictant,
                                     pot_pred_indices = varindex, train_thresholds = thres,
                                     test_thresholds = thres, maxnumbervars = maxvars){

  #initialize vector of chosen variables, AICscores and brier scores
  variables = c()
  briers = c()
  AICscores = list()
  results = list()

  #fit a model for each potential predictor and evaluate AIC score
  for(i in names(train_set[pot_pred_indices])){
    model = hxlr(reformulate(termlabels = names(data.frame(train_set[i])), response = as.name(names(train_set[predictant]))),
                 data=data.frame(train_set[predictant], train_set[i]), thresholds = train_thresholds)
    AICscores = append(AICscores,AIC(model))
  }

  #select the best predictor and add to first model
  added = pot_pred_indices[unlist(AICscores[seq(length(pot_pred_indices))])==min(unlist(AICscores))]
  variables = c(variables, added)
  firstmodel = hxlr(reformulate(termlabels = names(data.frame(train_set[variables])), response = as.name(names(train_set[predictant]))),
                    data=data.frame(train_set[predictant], train_set[variables]), thresholds = train_thresholds)
  modellist = list(firstmodel)

  ### ITERATION, ADDING VARIABLES
  while(length(variables) < maxvars){
    AICscores = list()
    for(i in names(train_set[pot_pred_indices[!pot_pred_indices%in%variables]])){
      model = hxlr(reformulate(termlabels = names(data.frame(train_set[variables], train_set[i])), response = as.name(names(train_set[predictant]))),
                   data=data.frame(train_set[predictant], train_set[variables], train_set[i]), thresholds = train_thresholds)
      AICscores = append(AICscores, AIC(model))
    }
    print(max(unlist(AICscores)))
    print(min(unlist(AICscores)))

    #update variable selection
    added = pot_pred_indices[!pot_pred_indices%in%variables][unlist(AICscores[seq(length(pot_pred_indices[!pot_pred_indices%in%variables]))])
                                                             ==min(unlist(AICscores))]
    variables = c(variables, added)
    bestmodel = hxlr(reformulate(termlabels = names(data.frame(train_set[variables])), response = as.name(names(train_set[predictant]))),
                     data=data.frame(train_set[predictant], train_set[variables]), thresholds = train_thresholds)
    modellist = append(modellist, list(bestmodel))
  }

  for (i in seq(modellist)){
    for(reg in regions){
      region_subset = filter(test_set, region == reg)
      observed = data.frame(matrix(rep(as.matrix(region_subset[predictant]),length(test_thresholds)),ncol=length(test_thresholds)))
      model_ver = modellist[[i]]
      values <- predict(model_ver, newdata = region_subset, type = "cumprob", thresholds = test_thresholds)
      verification_set <- verify(as.numeric(observed < test_thresholds), values, frcst.type = "prob", obs.type = "binary", title = "")
      eval = brier(as.numeric(observed < test_thresholds), values, bins = FALSE)
      brierval = eval$bs
      brierbase = eval$bs.baseline
      briers = append(briers, c(y, i, reg, brierval,brierbase))
      reliability.plot(verification_set, titl = paste(paste(names(model_ver$start)[seq(3,length(model_ver$start))],collapse=" + "),
                                                      " - Brier score = ", round(verification_set$bs,ndec)))
    }
  }
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
  set.seed(seq(numsubset)[q])
  randomsubset = round(runif(nrow(train_fin))*numsubset+0.5)
  train_sub = cbind(train_fin, subset = randomsubset)
  for(j in seq(numsubset)){
    #check approximately equal length of random subsets by printing relative length
    relweight_subset = sum(train_sub$subset[train_sub$subset == j])/nrow(train_sub)/j
    print(relweight_subset)

    #select training and testing dataset
    train_j = filter(train_sub, subset != j)[seq(length(climset))]
    test_j = filter(train_sub, subset == j)[seq(length(climset))]

    #find a fit
    result = fit_extended_logitModels(train_j, test_j, predictant = ind_predictant, pot_pred_indices = varindex,
                                      train_thresholds = thres, test_thresholds = thres_eval, maxnumbervars = maxvars)

    #put results in dataframes and vectors
    brierdataframe = rbind(brierdataframe, data.frame(test_year = (result$brierscores)[seq(1,length(result$brierscores),5)],
                                                      npredictors = (result$brierscores)[seq(2,length(result$brierscores),5)],
                                                      region = (result$brierscores)[seq(3,length(result$brierscores),5)],
                                                      brier_score = (result$brierscores)[seq(4,length(result$brierscores),5)],
                                                      brier_base = (result$brierscores)[seq(5,length(result$brierscores),5)]))
    models = append(models, result$models)
  }
  q = q+1 #update seed of random subset
}

for(y in years){
  test_fin = filter(climset, Year == y)
  result = fit_extended_logitModels(train_fin, test_fin, predictant = ind_predictant, pot_pred_indices = varindex,
                                    train_thresholds = thres, test_thresholds = thres_eval, maxnumbervars = maxvars)
  brierdataframe = rbind(brierdataframe, data.frame(test_year = (result$brierscores)[seq(1,length(result$brierscores),5)],
                                                    npredictors = (result$brierscores)[seq(2,length(result$brierscores),5)],
                                                    region = (result$brierscores)[seq(3,length(result$brierscores),5)],
                                                    brier_score = (result$brierscores)[seq(4,length(result$brierscores),5)],
                                                    brier_base = (result$brierscores)[seq(5,length(result$brierscores),5)]))
  models = append(models, result$models)
}
