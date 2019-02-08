########################################
#This is without any transformations!!!#
########################################

rm(list=ls(all=TRUE))
### LOAD EXTERNAL CODE
library(dplyr)
library(MASS)
library(verification)
library(devtools)
library(testthat)
usethis::use_testthat()

### SET GENERAL CONDITIONS FOR THE MODEL
#read dataset
setwd("/usr/people/groote/")
ObsPV = read.csv(file = "Thunderstorm_radar_merged.csv")
setwd("/usr/people/groote/ThunderstormPostProcessingv1/")
years = c(as.numeric(unique(ObsPV$Year)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[1]
VT = unique(ObsPV$validtime)[2]
regions = c(unique(ObsPV$region))
threshold = 1.50

#ObsPV = cbind(ObsPV, log(0.01+ObsPV$LidStrength_0mabovegnd_6hrlymin/ObsPV$PrecipitableWater_0mabovegnd_6hrlymax))


#set default available variables: predictant and predictor
numsubset = 3 #number of subsets for hyperparametersetting
ind_predictant = 107
varindex=c(seq(18,101))
pot_preds=names(ObsPV[varindex])
ndec = 4
maxsteps = 6

verify_LRmodel_per_reg <- function(train_set, test_set, model, predictant, nstepsAIC, thres){
  predictions = data.frame()

  if (length(train_set) != length(test_set)){
    message("Training and test sets have to have the same number of columns")
    return("Failed")
  }
  observed = as.numeric(test_set[predictant] > thres)
  values <- predict(model, newdata = test_set)
  probability = exp(values)/(1+exp(values))
  #last_brier = brier(observed, probability, bins = FALSE)$bs #alternative with bins - - verification_set$bs
  predictions = data.frame(prob = probability, obs = observed, npred = nstepsAIC, region = test_set$region)#append(briers, c(y, nstepsAIC, reg, round(last_brier, ndec)))
  return(predictions)

}

fit_logitModels_and_predict <- function(train_set, test_set, region_set = regions, predictant = ind_predictant, pot_pred_indices = varindex,
                                        thres = threshold, maxstepsAIC = 10, print_conf = FALSE){

  if(maxstepsAIC > length(pot_pred_indices)){
    message("You can not add more variables to the model than the number of available predictors")
    return("Failed")
  }

  #do the reference fit
  nullfit = glm(unlist(train_set[predictant] > thres) ~ 1, data = train_set[pot_pred_indices], family=binomial)

  #results storage
  modellist = list()
  briers = data.frame()
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

    verification_result = verify_LRmodel_per_reg(train_set, test_set, logitMod, predictant, nstepsAIC, thres)

    #make lists of model, probabilities and brier scores including skill score
    modellist = append(modellist,list(logitMod))
    briers = rbind(briers, verification_result)
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
q = 1
for(y in years){

  #create random subsets to train on
  train_fin = filter(ObsPV, Year.x != y & validtime.x == VT & leadtime_count == LT)
  set.seed(15+seq(years)[q]) #for reproducability purposes
  testdf = data.frame(validdate = unique(train_fin$validdate), subset = round(runif(unique(train_fin$validdate))*numsubset+0.5))
  train_sub <- left_join(train_fin, testdf, by = names(testdf)[1])

  for(j in seq(numsubset)){
    #check approximately equal length of random subsets by printing relative length
    relweight_subset = sum(train_sub$subset[train_sub$subset == j])/nrow(train_sub)/j
    print(relweight_subset)
    print(j)

    #select training and testing dataset
    train_j = filter(train_sub, subset != j)[seq(length(ObsPV))]
    test_j = filter(train_sub, subset == j)[seq(length(ObsPV))]

    #find a fit
    result = fit_logitModels_and_predict(train = train_j, test = test_j, region_set = regions, predictant = ind_predictant, pot_pred_indices = varindex, thres = threshold, maxstepsAIC = maxsteps, print_conf = FALSE)

    #put results in dataframes and vectors
    brierdataframe = rbind(brierdataframe, result$briers)
    models = append(models, result$model)
    nullfits = append(nullfits, result$nullfit)
  }
  q = q + 1
}

brierdataframe2 = data.frame()
for(y in years){
  train_fin = filter(ObsPV, Year.x != y & validtime.x == VT & leadtime_count == LT)
  test_fin = filter(ObsPV, Year.x == y & validtime.x == VT & leadtime_count == LT)
  result = fit_logitModels_and_predict(train_set = train_fin, test_set = test_fin,
                                       predictant = ind_predictant, pot_pred_indices = varindex,
                                       thres = 3.00, maxstepsAIC = maxsteps, print_conf = FALSE)
  brierdataframe2 = rbind(brierdataframe2, result$briers)
  models = append(models, result$model)
  nullfits = append(nullfits, result$nullfit)
}

LR_ss <- brierdataframe %>% group_by(npred, region) %>% summarise(bs = brier(obs = obs, pred = prob, bins = FALSE)$ss)
LR_bs <- brierdataframe %>% group_by(npred, region) %>% summarise(bs = brier(obs = obs, pred = prob, bins = FALSE)$bs)
LR_ss2 <- brierdataframe2 %>% group_by(npred, region) %>% summarise(bs = brier(obs = obs, pred = prob, bins = FALSE)$ss)
LR_bs2 <- brierdataframe2 %>% group_by(npred, region) %>% summarise(bs = brier(obs = obs, pred = prob, bins = FALSE)$bs)

nr=max(brierdataframe2$npred)
for(reg in regions){
  plot.new()
  for(pred in unique(brierdataframe2$npred)){
    subset = filter(brierdataframe2, npred == pred, region == reg)
    if(pred == 1){
      plot(data.frame(verify(subset$obs, subset$prob)[[8]],verify(subset$obs, subset$prob)[[9]]), xlim = c(0, 1), ylim = c(0, 1),
           legend.names = pred,
           col = rainbow(nr)[pred], type = "o", lwd = 2, xlab = "Forecasted probability", ylab = "Observed relative frequency",
           main = paste0("Reliability plot 0/1, region = ", reg))

    }else{
      lines(verify(subset$obs,subset$prob)[[8]],verify(subset$obs,subset$prob)[[9]], legend.names = pred, col = rainbow(nr)[pred], type = "o", lwd = 2)#, col = c(1-0.1*pred,1,1))
    }
    abline(0,1)
    legend(0,1, legend = seq(nr), col = rainbow(nr), lty = 1, lwd = 2)
  }
  text(0.05,1.02,"No. of predictors:")
}

# --------------------------------------------------
test_that("Test dataset complete?", {
  expect_equal(filter(ObsPV, validtime.x == VT & leadtime_count == LT) %>% arrange(Year, Month, Day), rbind(train_fin, test_fin) %>% arrange(Year, Month, Day))
})

test_that("Test resulting brier data frame for obvious errors",{
  expect_equal(unique(brierdataframe$region), regions)
  expect_lte(min(brierdataframe$npred),maxsteps)
})
set.seed(228)
x1 = rnorm(500,0,5)
x2 = rnorm(500,0,5)
x3 = rnorm(500,0,5)
x4 = rnorm(500,0,5)
pert = rnorm(500,0,10)
y1 = x1*100+x4*10+pert
testthat_dfLR = data.frame(y1,x1,x2,x3,x4, region = 1)
test_model = fit_logitModels_and_predict(testthat_dfLR, testthat_dfLR, region_set = c(1), 1, pot_pred_indices = c(2,5),
                                         thres = threshold, maxstepsAIC = 2, print_conf = FALSE)$model[[1]]
test_that("Function fit_logitModels_and_predict",{
  expect_equal(fit_logitModels_and_predict(testthat_dfLR, testthat_dfLR, region_set = c(1), 1, pot_pred_indices = seq(2,5),
                                           thres = threshold, maxstepsAIC = 2, print_conf = FALSE)[[1]][[1]][1:24],
               fit_logitModels_and_predict(testthat_dfLR, testthat_dfLR, region_set = c(1), 1, pot_pred_indices = c(2,5),
                                           thres = threshold, maxstepsAIC = 2, print_conf = FALSE)[[1]][[1]][1:24]) ## at index 25, the data frame used for fitting is in the result; these should be different
  expect_equal(fit_logitModels_and_predict(testthat_dfLR, testthat_dfLR, region_set = c(1), 1, pot_pred_indices = seq(2,5),
                                           thres = threshold, maxstepsAIC = 2, print_conf = FALSE)[[1]][[2]][1:24],
               fit_logitModels_and_predict(testthat_dfLR, testthat_dfLR, region_set = c(1), 1, pot_pred_indices = c(2,5),
                                           thres = threshold, maxstepsAIC = 2, print_conf = FALSE)[[1]][[2]][1:24]) ## at index 25, the data frame used for fitting is in the result; these should be different
  expect_equal(fit_logitModels_and_predict(testthat_dfLR, testthat_dfLR, region_set = c(1), 1, pot_pred_indices = c(2,5),
                                           thres = threshold, maxstepsAIC = 3, print_conf = FALSE),"Failed") ## at index 25, the data frame used for fitting is in the result; these should be different
  expect_error(fit_logitModels_and_predict(testthat_dfLR, testthat_dfLR, region_set = c(1), 1, pot_pred_indices = c(2,190),
                                           thres = threshold, maxstepsAIC = 2, print_conf = FALSE))
  expect_error(fit_logitModels_and_predict(testthat_dfLR, testthat_dfLR, region_set = c(3), 1, pot_pred_indices = c(2,190),
                                           thres = threshold, maxstepsAIC = 2, print_conf = FALSE))
  expect_error(fit_logitModels_and_predict(train_j[18,20], testthat_dfLR, region_set = c(3), 25, pot_pred_indices = c(2,5),
                                           thres = threshold, maxstepsAIC = 2, print_conf = FALSE))
})
