########################################
#This is without any transformations!!!# (except RI)
########################################

#######

rm(list=ls()[! ls() %in% c("LT_i","VT_i","LT_val","VT_val","varindex_shell")])
### LOAD EXTERNAL CODE
library(dplyr)
library(MASS)
library(Hmisc)
library(verification)
library(SpecsVerification)
library(arm)
library(crch)
library(profvis)
library(parallel)

writedir = "/usr/people/groote/ThunderstormPostProcessingv1/ELRres"

print("ELR_blocks4 active for")
print("LT_i equals:")
print(LT_i)
print("------")
print("VT_i equals:")
#print(VT_i)

### SET GENERAL CONDITIONS FOR THE MODEL
#set thresholds and hyperparameter; determine test dataset and training dataset
p = 0.25 #power transformation to linearize thresholds
maxvars = 4
nmembers = 25 #number of ensemble members to calculate CRPS
numsubset = 3 #number of subsets for hyperparameter selection
#thresholds used for training; percentiles
percmin = 50
percmax = 95
percint = 5
thres_eval = seq(25,400,25) #discharge threshold for evaluation
minpredictant = 1.5 #minimum sum of discharges considered as thunderstorm case

#read data, years, VT, LT and regions. Change VT, LT and regions for subset fitting.
setwd("/usr/people/groote/")
ObsPV = read.csv(file = "full_final00z_dataset2.csv")
#ObsPV$Dischargerate[ObsPV$Year==2015]=ObsPV$Dischargerate_flits[ObsPV$Year==2015]
#ObsPV$Ndischarge[ObsPV$Year==2015]=ObsPV$Ndischarge_flits[ObsPV$Year==2015]

colnames(ObsPV)[c(3,4,5,6,7,10,11)] <- c("Year","validdate","Month","Day","validtime2","region","leadtime_count")
ObsPV$leadtime_count = ObsPV$leadtime_count
ObsPV$validdate = ObsPV$validdate+ObsPV$Year*10000
setwd("/usr/people/groote/ThunderstormPostProcessingv1/")
years = c(as.numeric(unique(ObsPV$Year)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[LT_i]
#VT = unique(ObsPV$validtime)[VT_i]
VT = paste0("Init_00z_wconf_othercross_flits",length(varindex_shell))
regions = c(unique(ObsPV$region))#[1:2]

###
ObsPV$MSLP_mean = ObsPV$MSLP_mean - 100000
ObsPV$MSLP_max = ObsPV$MSLP_max - 100000
###

#selection from dataset
#ObsPV=cbind(ObsPV, selector = as.numeric(as.character(VT)==as.character(ObsPV$validtime))*as.numeric(as.character(LT)==as.character(ObsPV$leadtime_count)))
ObsPV=cbind(ObsPV, selector = as.numeric(as.character(LT)==as.character(ObsPV$leadtime_count)))
climset = filter(ObsPV, Ndischarge > minpredictant & selector == 1)
print(dim(climset))
#do transformations for thresholds
climset$Dischargerate = climset$Dischargerate^p
thres = unique(quantile(climset$Dischargerate,seq(percmin,percmax,percint)/100))
thres_eval = thres_eval^p
ndec = 4 #number of decimals usedi when appending scores to list of scores

#set available variables & predictant
varindex=varindex_shell
pot_preds=names(climset[varindex])
ind_predictant = 627

### Above this point, the settings for a run have been defined!! #####
### Functions defined: ####
fit_test_all_pot_pred <- function(train_set, predictant, pot_pred_indices, train_thresholds, used_preds = c(0)){
  #this function selects the best predictor by forward fittng; trying potential predictors and selecting the one with lowest AIC
  AICscores = list()
  for(i in names(train_set[pot_pred_indices])){
   # print(names(train_set[i]))
   # print(names(train_set[used_preds]))
    model = hxlr(reformulate(termlabels = names(data.frame(train_set[i], train_set[used_preds])),
                             response = as.name(names(train_set[predictant]))),
                 data=data.frame(train_set[predictant], train_set[i], train_set[used_preds]),
                 thresholds = train_thresholds)
    AICscores = append(AICscores,AIC(model))
  }
  added = pot_pred_indices[unlist(AICscores[seq(length(pot_pred_indices))])==min(unlist(AICscores))]
  if(length(added)>1){
    added=added[1]
  }
  return(added)
}

new_verify_ELR <- function(test_set, model, predictant, test_threshold, npred, reliabilityplot = FALSE){

  #get observations
  observed = data.frame(observation = as.numeric(test_set[predictant]>test_threshold))

  #predict with model to do the verification
  values <- as.numeric(1-predict(model, newdata = test_set, type = "cumprob", test_threshold))
  result = data.frame(observation = observed, probability = values, npredictors = npred, threshold = test_threshold, validdate = test_set$validdate)
  return(result)
}

getEnsCRPS <- function(model, test_set, members, observations){
  #this function calculates the ensemble CRPS score, making use of quantiles and inverting the ELR model that it gets
  quans = (seq(members)-0.5)/members # these are the quantiles to be computed with given number of members
  logtrans = log(1/quans - 1)
  logtrans = matrix(rep(logtrans, each = dim(test_set)[1]), ncol = members) #repetition of transformed quantiles in matrix

  #apply the XLR model
  if(length(model$coefficients$scale)!=0){
    logtrans = logtrans * exp(matrix(rep(as.matrix(test_set[names(model$coefficients$scale)])%*%model$coefficients$scale, members),ncol = members))
  }
  modelpred = logtrans +matrix(rep(-model$coefficients$intercept[1]+as.matrix(test_set[names(model$coefficients$location)])%*%model$coefficients$location, members), ncol = members)
  result = modelpred/model$coefficients$intercept[2]
  result = cbind(result, npredictors = length(model$coefficients$location)+length(model$coefficients$scale), observations)
  return(result)
}

fit_extended_logitModels <- function(train_set, test_set, predictant = ind_predictant, pot_pred_indices = varindex,
                                     train_thresholds = thres, test_thresholds = thres, maxnumbervars = maxvars){

  #initialize vector of chosen variables and brier scores
  variables = c()
  verified_set = data.frame()
  results = list()
  n_variables = c()
  crpsscores = data.frame()

  if(maxnumbervars > length(pot_pred_indices)){
    message("You can not add more variables to the model than the number of available predictors")
    return("failed")
  }

  #select the best predictor and add to first model in list of variables
  added = fit_test_all_pot_pred(train_set, predictant, pot_pred_indices, train_thresholds)
  variables = c(variables, added)

  #add this first model with best predictor variable to the model listh
  firstmodel = hxlr(reformulate(termlabels = names(data.frame(train_set[variables])), response = as.name(names(train_set[predictant]))), data=data.frame(train_set[predictant], train_set[variables]), thresholds = train_thresholds)
  modellist = list(firstmodel)
  n_variables = append(n_variables, length(variables))

  ### ITERATION, ADDING VARIABLES
  while(length(variables) < maxnumbervars){
    #update potential predictors remaining
    remaining_indices = pot_pred_indices[!pot_pred_indices%in%variables]

    #get best new predictor variable and add to variables and save the index
    added = fit_test_all_pot_pred(train_set, predictant, remaining_indices, train_thresholds, used_preds = variables)
    variables = c(variables, added)

    #add model based on these variables to model list, including number of variables in the model
    bestmodel = hxlr(reformulate(termlabels = names(data.frame(train_set[variables])), response = as.name(names(train_set[predictant]))), data=data.frame(train_set[predictant], train_set[variables]), thresholds = train_thresholds)
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
make_subset <- function(i){
  verifdate = verifdates[i]
  subset = predictiondataframe[predictiondataframe$validdate == verifdate,]
  return(subset)

}
do_bootstrap <- function(i){
  datelist = unique(predictiondataframe$validdate)
  set.seed(i+950)
  randdates = round(0.5+length(unique(predictiondataframe$validdate))*runif(length(unique(predictiondataframe$validdate))))
  verifdates = datelist[randdates]
  verifydataframe = data.frame()
  #assign("verifydataframe", verifydataframe, .GlobalEnv)
  assign("verifdates", verifdates, .GlobalEnv)
  verifydataframe = do.call(rbind,lapply(seq(1,length(verifdates)),make_subset))
  skillscore = brier(verifydataframe$observation, verifydataframe$probability, rep(mean(predictiondataframe$observation),length(predictiondataframe$observation)), bins = FALSE)$ss
  #print(skillscore)
  return(skillscore)
}
block_bootstrapping <- function(predictiondataframe,ntimes,interval){
  skillscores = data.frame()
  assign("predictiondataframe", predictiondataframe, .GlobalEnv)
  print(predictiondataframe$threshold)
  clusterExport(cl = cls, c("make_subset", "do_bootstrap","predictiondataframe","brier"))
  skillscores = unlist(parLapply(cl = cls, X = seq(1,ntimes), fun = do_bootstrap))
  #print(skillscores)
  funcresult = list()
  funcresult$up = quantile(skillscores, (0.5+interval/2))
  funcresult$low = quantile(skillscores, (0.5-interval/2))
  rm("predictiondataframe")
  return(funcresult)

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
  set.seed(18+seq(years)[q]) #for reproducability purposes

  #make three data sets for about 2/3 vs. 1/3 of the days within 2 years data set
  testdf = data.frame(validdate = unique(train_fin$validdate), subset = round(runif(unique(train_fin$validdate))*numsubset+0.5))
  train_sub <- left_join(train_fin, testdf, by = names(testdf)[1])
  #train_sub = cbind(train_fin, subset = randomsubset)
  for(j in seq(numsubset)){
    #check approximately equal length of random subsets by printing relative length
    relweight_subset = sum(train_sub$subset[train_sub$subset == j])/nrow(train_sub)/j
    print(relweight_subset)
    print(j)

    #finally select training and testing dataset
    train_j = filter(train_sub, subset != j)[seq(length(climset))]
    test_j = filter(train_sub, subset == j)[seq(length(climset))]

    #find fitting models with the function
    result = fit_extended_logitModels(train_j, test_j, predictant = ind_predictant, pot_pred_indices = varindex,
                                      train_thresholds = thres, test_thresholds = thres_eval, maxnumbervars = maxvars)
    # print(result$briers)

    #put results (predictions, models and CRPS-score) in dataframes and vectors
    brierdataframe = rbind(brierdataframe, result$verification)
    models = append(models, result$models)
    crpsscorelist = rbind(crpsscorelist, result$crpsscores)
  }
  q = q+1 #update seed set for the next year in 9-fold
}

# final cross validation section

models2 = list()
crpsscorelist2 = data.frame()
brierdataframe2 = data.frame()
set.seed(25) #for reproducability purposes
#make three data sets for about 2/3 vs. 1/3 of the days within 2 years data set
testdf = data.frame(validdate = unique(climset$validdate), subset = round(runif(unique(climset$validdate))*numsubset+0.5))
train_sub <- left_join(climset, testdf, by = names(testdf)[1])

for(j in seq(numsubset)){

  #select two out of three years each time
  test_fin = filter(train_sub, subset == j)
  train_fin = filter(train_sub, subset != j)

  #apply function to find models with given training sets
  result = fit_extended_logitModels(train_fin, test_fin, predictant = ind_predictant, pot_pred_indices = varindex,
                                    train_thresholds = thres, test_thresholds = thres_eval, maxnumbervars = maxvars)

  #calculate results: models, verification and CRPS-scores
  brierdataframe2 = rbind(brierdataframe2, data.frame(result$verification))
  crpsscorelist2 = rbind(crpsscorelist2, result$crpsscores)
  models2 = append(models2, result$models)
}
#calculate Brier Skill Scores for each number of predictors with threshold for both 9-fold and final cross validation
#for final cross validation, a 2 is added in each name (see above)
thescores=data.frame()
crpsdata=data.frame()
cls <- makeCluster(getOption("cl.cores", 8))
for(npred in unique(brierdataframe$npredictors)){
  for(thresh in unique(brierdataframe$threshold)){
    subset = filter(brierdataframe, npredictors == npred, threshold == thresh)
    subset2 = filter(brierdataframe2, npredictors == npred, threshold == thresh)
    score = brier(subset$observation, subset$probability, bins = FALSE)$ss
    conf = block_bootstrapping(subset, 1000, 0.95)
    score2 = brier(subset2$observation, subset2$probability, bins = FALSE)$ss
    conf2 = block_bootstrapping(subset2, 1000, 0.95)
    thescores = rbind(thescores, data.frame(numpredictors = npred, threshold = thresh^4, ss_9fold = score,ss_9fold_up = conf$up, ss_9fold_low = conf$low, ss_years = score2, ss_years_up = conf2$up, ss_years_low = conf2$low))
    #crpsdata = mean(EnsCrps(subset4[-length(subset4)],subset4[length(subset4)]))
  }
  subset3 = filter(crpsscorelist, npredictors == npred)
  subset4 = filter(crpsscorelist2, npredictors == npred)
  colval = length(subset4)
  crpsdata = rbind(crpsdata, data.frame(npredictors = npred, score_9fold = mean(EnsCrps(ens = as.matrix(subset3[-c(colval-1,colval)]),obs = as.matrix(subset3[colval]))), score_years = mean(EnsCrps(ens = as.matrix(subset4[-c(colval-1,colval)]), obs = as.matrix(subset4[colval])))))
}
stopCluster(cls)
#make reliability plot
setwd(writedir)
nr=max(brierdataframe2$npred)
for(thresh in thres_eval){
  subsetting = filter(brierdataframe2, threshold == thresh)
  png(file=paste0("Relplot_",VT,"_LT_",LT,"_npred_",length(varindex),"_thres_",thresh,".png"), width=720, height = 400)
  plot.new()
  barplotdata = data.frame()
  for(pred in unique(brierdataframe2$npredictors)){
    subset = filter(subsetting, npredictors == pred)
    SkSc = filter(thescores, threshold == thresh^4)$ss_years

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
  text(0.7,1.02,"No. of predictors:")
  legend(0.45,1, legend = paste0(seq(nr), " with BSS ",round(SkSc[seq(nr)],3)), col = rainbow(nr), lty = 0, pch=1, lwd = 2)
  dev.off()
  print(length(dev.list()))
}

refcrps = mean(EnsCrps(t(matrix(rep(t(climset[ind_predictant]),dim(climset[ind_predictant])[1]), nrow=dim(climset[ind_predictant])[1])),as.matrix(climset[ind_predictant])))
ss_years=1-(crpsdata$score_years)/refcrps
crpsdata = cbind(crpsdata, SS = ss_years)
print(crpsdata)
#plot, print and write to CSV
png(file=paste0("Npred_",VT,"_LT_",LT,"_npred_",length(varindex),".png"), width=720, height = 400)
plot(thescores$threshold,thescores$ss_9fold, xlab = "Lightning intensity threshold (dis./5 min) - 9-fold cross validation", ylab="BSS", main = "Verification score as function of threshold", col=thescores$numpredictors, legend.names=thescores$numpredictors, ylim = c(-0.3,0.6))
legend(72,0.5,legend=seq(nr),col=unique(thescores$numpredictors),pch=1)
text(70,0.55,"No. predictors:")
dev.off()
print(length(dev.list()))
png(file=paste0("SSyears_",VT,"_LT_",LT,"_npred_",length(varindex),".png"), width=720, height = 400)
plot(thescores$threshold,thescores$ss_years, xlab = "Lightning intensity threshold (dis./5 min) - final verification", ylab="BSS", main = "Verification score as function of threshold", col=thescores$numpredictors, legend.names=thescores$numpredictors, ylim = c(-0.3,0.6))
legend(72,0.5,legend=seq(nr),col=unique(thescores$numpredictors),pch=1)
text(70,0.55,"No. predictors:")
dev.off()
write.csv(thescores,paste0("ELR_scores_",VT,"_LT_",LT,"_npred_",length(varindex),".csv"))
#crpsdata=data.frame(npred=rep(seq(1,maxvars),numsubset),crps=matrix(as.character(crpsscorelist2)))
name = paste0("ELR_CRPSscores_",VT,"_LT_",LT,"_npred_",length(varindex),".csv")
write.csv(crpsdata,name)
saveRDS(models2, paste0("Models_",VT,"_LT_",LT,"_npred_",length(varindex)))

