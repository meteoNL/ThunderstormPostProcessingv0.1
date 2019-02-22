########################################
#This is without any transformations!!!# (except RI)
########################################

#######

### LOAD EXTERNAL CODE
library(dplyr)
library(MASS)
library(Hmisc)
library(verification)
library(SpecsVerification)
library(arm)
library(crch)
library(profvis)
rm(list=ls(all=TRUE))
### SET GENERAL CONDITIONS FOR THE MODEL
#set thresholds and hyperparameter; determine test dataset and training dataset
p = 0.25 #power transformation to linearize thresholds
maxvars = 4
nmembers = 10 #number of members/quantiles to calculate CRPS
numsubset = 3 #number of subsets for hyperparameter selection
percmin = 50
percmax = 95
percint = 5
thres_eval = seq(1.8,3.0,0.05)^4 #discharge threshold for evaluation
minpredictant = 1.5 #minimum number of discharges sum considered as thunderstorm occurence

#read data, years, VT, LT and regions. Change VT, LT and regions for subset fitting.
setwd("/usr/people/groote/")
ObsPV = read.csv(file = "Thunderstorm_radar_merged.csv")
setwd("/usr/people/groote/ThunderstormPostProcessingv1/")
years = c(as.numeric(unique(ObsPV$Year)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[1]
VT = unique(ObsPV$validtime.x)[2]
regions = c(unique(ObsPV$region))#[1:2]

#selection from dataset
climset = filter(ObsPV, Ndischarge > minpredictant & validtime.x == VT & leadtime_count == LT)

#do transformations for thresholds
climset$Dischargerate = climset$Dischargerate^p
thres = quantile(climset$Dischargerate,seq(percmin,percmax,percint)/100)
thres_eval = thres_eval^p
ndec = 4 #number of decimals usedi when appending scores to list of scores

#set available variables & predictant
#varindex=c(seq(18,71),seq(73,76),seq(80,82),seq(85,101))
varindex=c(seq(18,27),seq(30,37),seq(39,42),seq(44,57),seq(59,67),seq(70,71),seq(73,76),seq(81,82),seq(86,99))
pot_preds=names(climset[varindex])
ind_predictant = 113

### Above this point, the settings for a run have been defined!! #####
### Functions defined: ####
fit_test_all_pot_pred <- function(train_set, predictant, pot_pred_indices, train_thresholds, used_preds = c(0), used_sc_preds = c(0)){
  #selects the best predictor by forward fittng, trying potential predictors and selecting the one with lowest AIC
  # for general case here, with both scale and location predictors
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
  if(length(added)>1){
    added=added[1]
  }
  return(added)
}

fit_test_all_pot_pred_no_sc <- function(train_set, predictant, pot_pred_indices, train_thresholds, used_preds = c(0)){
  #selects the best predictor by forward fittng, trying potential predictors and selecting the one with lowest AIC
  #here for models without scale parameter variables selected
  AICscores = list()
  for(i in names(train_set[pot_pred_indices])){
    model = hxlr(reformulate(termlabels = paste0(reformulate(termlabels = names(data.frame(train_set[i], train_set[used_preds]))))[2],
                             response = as.name(names(train_set[predictant]))),
                 data=data.frame(train_set[predictant], train_set[i], train_set[used_preds]),
                 thresholds = train_thresholds)
    AICscores = append(AICscores,AIC(model))
  }
  added = pot_pred_indices[unlist(AICscores[seq(length(pot_pred_indices))])==min(unlist(AICscores))]
  if(length(added)>1){
    added=added[1]
    print(added)
  }
  return(added)
}
fit_test_all_pot_pred_sc <- function(train_set, predictant, pot_pred_indices, train_thresholds, used_preds = c(0), used_sc_preds = c(0)){
  #selects the best predictor by forward fittng, trying potential predictors and selecting the one with lowest AIC
  # here for finding best scale predictor from predictor set
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
  if(length(added)>1){
    added=added[1]
    print(added)
  }
  return(added)
}
new_verify_ELR <- function(test_set, model, predictant, test_threshold, npred, reliabilityplot = FALSE){
  #get observations
  observed = data.frame(observation = as.numeric(test_set[predictant]>test_threshold))

  #predict with model to verify model; do predictions
  values <- as.numeric(1-predict(model, newdata = test_set, type = "cumprob", test_threshold))
  result = data.frame(observation = observed, probability = values, npredictors = npred, threshold = test_threshold)
  return(result)
}

getEnsCRPS <- function(model, test_set, members, observations){
  quans = (seq(members)-0.5)/members # these are the quantiles to be computed with given number of members
  logtrans = log(1/quans - 1)
  logtrans = matrix(rep(logtrans, each = dim(test_set)[1]), ncol = members) #repetition of transformed quantiles in matrix

  #apply the XLR model
  if(length(model$coefficients$scale)!=0){
    logtrans = logtrans * exp(matrix(rep(as.matrix(test_set[names(model$coefficients$scale)])%*%model$coefficients$scale, members),ncol = members))
  }
  modelpred = logtrans +matrix(rep(-model$coefficients$intercept[1]+as.matrix(test_set[names(model$coefficients$location)])%*%model$coefficients$location, members), ncol = members)
  modelpred = modelpred/model$coefficients$intercept[2]

  #calculate ensemble CRPS
  result = modelpred/model$coefficients$intercept[2]
  result = cbind(result, npredictors = length(model$coefficients$location)+length(model$coefficients$scale), observations)
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
  crpsscores = data.frame()

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

  #select model (with accompanying number of variables) and put its predictions for verification in a data frame
  for (i in seq(modellist)){
    model_ver = modellist[[i]]
    n_variables_i = n_variables[i]

    for(j in seq(length(test_thresholds))){
      test_threshold = test_thresholds[j]
      verified_set = rbind(verified_set, new_verify_ELR(test_set, model_ver, predictant, test_threshold, n_variables_i, reliabilityplot = FALSE))
    }
    #calculate model Ensemble CRPS-score
    crpsscores = rbind(crpsscores, getEnsCRPS(model_ver, test_set,nmembers, as.numeric(unlist(test_set[predictant]))))
  }

  #add results (model, its predictions and its CRPS-score) to results list
  results$models = modellist
  results$crpsscores = crpsscores
  results$verification = verified_set
  return(results)
}

# -----------------------------------------------------------------------------
### End of the functions ####

brierdataframe = data.frame()
models = list()
crpsscorelist = data.frame()

q = 1
#do procedure for 9-fold cross validation: select two years of data from data frame
for(y in years){
  train_fin = filter(climset, Year != y)
  set.seed(15+seq(years)[q]) #for reproducability purposes

  #make three data sets for about 2/3 vs. 1/3 of the days within 2 years data set
  testdf = data.frame(validdate = unique(train_fin$validdate), subset = round(runif(unique(train_fin$validdate))*numsubset+0.5))
  train_sub <- left_join(train_fin, testdf, by = names(testdf)[1])
  #train_sub = cbind(train_fin, subset = randomsubset)
  for(j in seq(numsubset)){
    #check approximately equal length of random subsets by printing relative length
    relweight_subset = sum(train_sub$subset[train_sub$subset == j])/nrow(train_sub)/j
    print(relweight_subset)
    print(j)

    ##finally select training and testing dataset
    train_j = filter(train_sub, subset != j)[seq(length(climset))]
    test_j = filter(train_sub, subset == j)[seq(length(climset))]

    #find fitting models with the function
    result = fit_extended_logitModels(train_j, test_j, predictant = ind_predictant, pot_pred_indices = varindex,
                                      train_thresholds = thres, test_thresholds = thres_eval, maxnumbervars = maxvars)
   # print(result$briers)

    #put results (predictions, models and CRPS-score) in dataframes and vectors
    brierdataframe = rbind(brierdataframe, result$verification)
    crpsscorelist = rbind(crpsscorelist, result$crpsscores)
    models = append(models, result$models)
  }
  q = q+1 #update seed set
}

# final cross validation section

models2=list()
crpsscorelist2 = data.frame()
brierdataframe2 = data.frame()
for(y in years){
  #select two out of three years each time
  test_fin = filter(climset, Year == y)
  train_fin = filter(climset, Year != y)

  #apply function to find models with given training sets
  result = fit_extended_logitModels(train_fin, test_fin, predictant = ind_predictant, pot_pred_indices = varindex,
                                    train_thresholds = thres, test_thresholds = thres_eval, maxnumbervars = maxvars)
  #calculate results: models, verification and CRPS-scores
  brierdataframe2 = rbind(brierdataframe2, data.frame(result$verification))
  models2 = append(models2, result$models)
  crpsscorelist2 = rbind(crpsscorelist2, result$crpsscores)
}
#make reliability plot
nr=max(brierdataframe2$npred)
for(thresh in thres_eval){
  subsetting = filter(brierdataframe2, threshold == thresh)
  plot.new()
  barplotdata = data.frame()
  for(pred in unique(brierdataframe2$npredictors)){
    subset = filter(subsetting, npredictors == pred)

    #get data for verification plot
    verified = verify(subset$obs, subset$prob)
    barplotdata = rbind(barplotdata, verified[[10]])
    names(barplotdata) = verified[[8]]

    #and then plot
    par(mar=c(4,4,4,28))
    if(pred == 1){
      plot(data.frame(verified[[8]],verified[[9]]), xlim = c(0, 1), ylim = c(0, 1),
           legend.names = pred,
           col = rainbow(nr)[pred], type = "o", lwd = 2, xlab = "Forecasted probability (-)", ylab = "Observed relative frequency (-)",
           main = paste0("Reliability plot thresholds at interval ",round(thresh^4)," dis./5 min."))

    }else{
      lines(verified[[8]],verified[[9]], legend.names = pred, col = rainbow(nr)[pred], type = "o", lwd = 2, main = "lineplot")#, col = c(1-0.1*pred,1,1))
    }
    lines(c(0,1),c(0,1))
    legend(0,1, legend = seq(nr), col = rainbow(nr), lty = 1, lwd = 2)
  }
  barplotdata[is.na(barplotdata)]<-0
  text(0.14,1.02,"No. of predictors:")
  par(mar=c(4,28,4,4))
  subplot(plot(data.frame(rep(verified[[8]],each=nr),unlist(barplotdata)),
               ylim = c(0,1),xlim=c(0,1),xlab = "Forecasted probability (-)",
               ylab = "Relative forecasting frequency (-)",
               col = rainbow(nr)[seq(nr)], main = "Forecasts issued for each no. of pred.", lwd=2),
          x = c(0.0,1.0),y = c(0.0,1.0), size = c(5,5))
  text(0.8,1.02,"No. of predictors:")
  legend(0.75,1, legend = seq(nr), col = rainbow(nr), lty = 0, pch=1, lwd = 2)
}
#calculate Brier Skill Scores for each number of predictors with threshold for both 9-fold and final cross validation
#for final cross validation, a 2 is added in each name (see above)
thescores=data.frame()
crpsdata=data.frame()
for(npred in unique(brierdataframe$npredictors)){
  for(thresh in unique(brierdataframe$threshold)){
    subset = filter(brierdataframe, npredictors == npred, threshold == thresh)
    subset2 = filter(brierdataframe2, npredictors == npred, threshold == thresh)
    score = brier(subset$observation, subset$probability, bins = FALSE)$ss
    score2 = brier(subset2$observation, subset2$probability, bins = FALSE)$ss
    thescores = rbind(thescores, data.frame(numpredictors = npred, threshold = thresh^4, ss_9fold = score, ss_years = score2))
    #crpsdata = mean(EnsCrps(subset4[-length(subset4)],subset4[length(subset4)]))
  }
  subset3 = filter(crpsscorelist, npredictors == npred)
  subset4 = filter(crpsscorelist2, npredictors == npred)
  colval = length(subset4)
  crpsdata = rbind(crpsdata, data.frame(npredictors = npred, score_9fold = mean(EnsCrps(ens = as.matrix(subset3[-c(colval-1,colval)]),obs = as.matrix(subset3[colval]))), score_years = mean(EnsCrps(ens = as.matrix(subset4[-c(colval-1,colval)]), obs = as.matrix(subset4[colval])))))
}

refcrps = mean(EnsCrps(t(matrix(rep(t(climset[ind_predictant]),dim(climset[ind_predictant])[1]), nrow=dim(climset[ind_predictant])[1])),as.matrix(climset[ind_predictant])))
ss_years=1-(crpsdata$score_years)/refcrps
crpsdata = cbind(crpsdata, SS = ss_years)
print(crpsdata)
#plot, print and write to CSV
plot(thescores$threshold,thescores$ss_years, xlab = "Lightning intensity threshold (dis./5 min)", ylab="BSS", main = "Verification score as function of threshold", col=thescores$numpredictors, legend.names=thescores$numpredictors)
legend(75,max(thescores$ss_years-0.02),legend=seq(nr),col=unique(thescores$numpredictors),pch=1)
setwd("/usr/people/groote/ThunderstormPostProcessingv1/HELRres")
write.csv(thescores,paste0("ELR_scores_",VT,"_LT_",LT,"_npred_",length(varindex),".csv"))
#crpsdata=data.frame(npred=rep(seq(1,maxvars),numsubset),crps=matrix(as.character(crpsscorelist2)))
name = paste0("HELR_CRPSscores_",VT,"_LT_",LT,"_npred_",length(varindex),".csv")
write.csv(crpsdata,name)
saveRDS(models2, paste0("Models_",VT,"_LT_",LT,"_npred_",length(varindex)))
#-----------------------------------------------------------------
