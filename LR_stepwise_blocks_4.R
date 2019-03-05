########################################
#This is without any transformations!!!#
########################################

rm(list=ls()[! ls() %in% c("LT_i","VT_i","LT_val","VT_val","varindex_shell")])
### LOAD EXTERNAL CODE
library(dplyr)
library(MASS)
library(Hmisc)
library(verification)
library(devtools)
library(testthat)
usethis::use_testthat()

print("LR_stepwise_blocks active for")
print("LT_i equals:")
print(LT_i)
print("------")
print("VT_i equals:")
print(VT_i)

### SET GENERAL CONDITIONS FOR THE MODEL
#read dataset
setwd("/usr/people/groote/")
ObsPV = read.csv(file = "ECMWF_merged4.csv")
setwd("/usr/people/groote/ThunderstormPostProcessingv1/")

# get years, regions, VT and LT in data set and set values for selection/subset to fit
years = c(as.numeric(unique(ObsPV$Year)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[LT_i]
VT = unique(ObsPV$validtime)[VT_i]
regions = c(unique(ObsPV$region))
threshold = 1.50 #number of discharges threshold for thunderstorm case

#set default available variables: predictant and predictor indices and number subsets, max number of predictors etc
numsubset = 3 #number of subsets for hyperparametersetting
ind_predictant = 107
varindex=varindex_shell
pot_preds=names(ObsPV[varindex])
ndec = 4
maxsteps = 6

### Above this point, the settings for a run have been defined!! #####
### Functions defined: ####
verify_LRmodel_per_reg <- function(train_set, test_set, model, predictant, nstepsAIC, thres){
  # this function predicts model outcomes for each model setting and discriminates by region
  predictions = data.frame()

  if (length(train_set) != length(test_set)){
    message("Training and test sets have to have the same number of columns")
    return("Failed")
  }

  # read observations
  observed = as.numeric(test_set[predictant] > thres)

  # do predictions with the model and convert it to appropriate probability of lightning
  values <- predict(model, newdata = test_set)
  probability = exp(values)/(1+exp(values))

  # make data frame with predicted thunderstorm probability, observations, region and number of predictors (all used for verification)
  predictions = data.frame(prob = probability, obs = observed, npred = nstepsAIC, region = test_set$region)
  return(predictions)
}

fit_logitModels_and_predict <- function(train_set, test_set, region_set = regions, predictant = ind_predictant, pot_pred_indices = varindex,
                                        thres = threshold, maxstepsAIC = 10, print_conf = FALSE){
  # this function does the fitting procedure
  if(maxstepsAIC > length(pot_pred_indices)){
    message("You can not add more variables to the model than the number of available predictors")
    return("Failed")
  }

  #do the reference fit
  nullfit = glm(unlist(train_set[predictant] > thres) ~ 1, data = train_set[pot_pred_indices], family=binomial)

  #results storage space declared
  modellist = list()
  briers = data.frame()
  results = list()

  #stepwise LR model fitting procedure
  for (nstepsAIC in seq(maxstepsAIC)){
    logitMod = stepAIC(nullfit, scope = list(upper = lm(unlist(train_set[predictant] > thres) ~ .,
                                                        data=train_set[pot_pred_indices]),
                                             lower = ~ 1), trace = 0, steps=nstepsAIC)

    #optionally print confidence interval (which is based on independence assumptions between predictors that won't be realistic)
    if(print_conf == TRUE){
      conf=exp(confint.default(logitMod))
      print(conf)
    }

    #do predictions with the model and get outcomes
    verification_result = verify_LRmodel_per_reg(train_set, test_set, logitMod, predictant, nstepsAIC, thres)

    #make lists of model and predicted probabilities which is later used for the scores
    modellist = append(modellist,list(logitMod))
    briers = rbind(briers, verification_result)
  }

  #returned values
  results$model = modellist
  results$briers = briers
  results$nullfit = list(nullfit)
  return(results)
}
#-------------------------------------
### End of the functions ####

#create memory for models and their evaluation
brierdataframe = data.frame()
models = list()
nullfits = list()
q = 1
ObsPV=cbind(ObsPV, selector = as.numeric(as.character(VT)==as.character(ObsPV$validtime))*as.numeric(as.character(LT)==as.character(ObsPV$leadtime_count)))

#do procedure for 9-fold cross validation: select two years of data from data frame
for(y in years){

  #create random subsets to train on
  train_fin = filter(ObsPV, Year != y & selector == 1)
  set.seed(15+seq(years)[q]) #for reproducability purposes

  #make three data sets for about 2/3 vs. 1/3 of the days within 2 years data set
  testdf = data.frame(validdate = unique(train_fin$validdate), subset = round(runif(unique(train_fin$validdate))*numsubset+0.5))
  train_sub <- left_join(train_fin, testdf, by = names(testdf)[1])

  for(j in seq(numsubset)){
    #check approximately equal length of random subsets by printing relative length
    relweight_subset = sum(train_sub$subset[train_sub$subset == j])/nrow(train_sub)/j
    print(relweight_subset)
    print(j)

    #select training and testing dataset for each model
    train_j = filter(train_sub, subset != j)[seq(length(ObsPV))]
    test_j = filter(train_sub, subset == j)[seq(length(ObsPV))]

    #do the fits
    result = fit_logitModels_and_predict(train = train_j, test = test_j, region_set = regions, predictant = ind_predictant, pot_pred_indices = varindex, thres = threshold, maxstepsAIC = maxsteps, print_conf = FALSE)

    #calculate results: models, verification/predictions
    brierdataframe = rbind(brierdataframe, result$briers)
    models = append(models, result$model)
    nullfits = append(nullfits, result$nullfit)
  }
  q = q + 1
}

brierdataframe2 = data.frame()
# final cross validation section

for(y in years){
  #select two of the three years of data
  train_fin = filter(ObsPV, Year != y & selector == 1)
  test_fin = filter(ObsPV, Year == y & selector == 1)

  #do the fitting procedure
  result = fit_logitModels_and_predict(train_set = train_fin, test_set = test_fin,
                                       predictant = ind_predictant, pot_pred_indices = varindex,
                                       thres = 3.00, maxstepsAIC = maxsteps, print_conf = FALSE)
  #calculate results: models, verification/predictions
  brierdataframe2 = rbind(brierdataframe2, result$briers)
  models = append(models, result$model)
  nullfits = append(nullfits, result$nullfit)
}

#calculate brier scores and skill scores
LR_ss <- brierdataframe %>% group_by(npred, region) %>% summarise(bs = brier(obs = obs, pred = prob, bins = FALSE)$ss)
LR_bs <- brierdataframe %>% group_by(npred, region) %>% summarise(bs = brier(obs = obs, pred = prob, bins = FALSE)$bs)
LR_ss2 <- brierdataframe2 %>% group_by(npred, region) %>% summarise(bs = brier(obs = obs, pred = prob, bins = FALSE)$ss)
LR_bs2 <- brierdataframe2 %>% group_by(npred, region) %>% summarise(bs = brier(obs = obs, pred = prob, bins = FALSE)$bs)

setwd("/usr/people/groote/ThunderstormPostProcessingv1/LRres")
#make reliability plot
print(length(dev.list()))
nr=max(brierdataframe2$npred)
for(reg in regions){
  png(file=paste0("Relplot_",VT,"_LT_",LT,"_npred_",length(varindex),"_reg_",reg,".png"), width=720, height = 400)
  plot.new()
  barplotdata = data.frame()
  for(pred in unique(brierdataframe2$npred)){
    subset = filter(brierdataframe2, npred == pred, region == reg)

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
           main = paste0("Reliability plot 0/1, region = ", reg))

    }else{
      lines(verified[[8]],verified[[9]], legend.names = pred, col = rainbow(nr)[pred], type = "o", lwd = 2, main = "lineplot")#, col = c(1-0.1*pred,1,1))
    }
    lines(c(0,1),c(0,1))
    legend(0,1, legend = seq(nr), col = rainbow(nr), lty = 1, lwd = 2)
  }
  barplotdata[is.na(barplotdata)]<-0
  text(0.15,1.02,"No. of predictors:")
  par(mar=c(4,28,4,4))
  subplot(plot(data.frame(rep(verified[[8]],each=nr),unlist(barplotdata)),
               ylim = c(0,1),xlim=c(0,1),xlab = "Forecasted probability (-)",
               ylab = "Relative forecasting frequency (-)",
               col = rainbow(nr)[seq(nr)], main = "Forecasts issued for each no. of pred.", lwd=2),
          x = c(0.0,1.0),y = c(0.0,1.0), size = c(5,5))
  text(0.8,1.02,"No. of predictors:")
  legend(0.75,1, legend = seq(nr), col = rainbow(nr), lty = 0, pch=1, lwd = 2)
  #dev.copy(pdf,paste0("Relplot_",VT,"_LT_",LT,"_npred_",length(varindex),"_reg_",reg,".pdf"))
  dev.off()
  print(length(dev.list()))
}
png(file=paste0("Npred_",VT,"_LT_",LT,"_npred_",length(varindex),".png"), width=720, height = 400)
plot(LR_ss$npred,LR_ss$bs,col=LR_ss$region, xlab="Number of predictors (-)", ylab = "BSS",xlim=c(0,nr+2),ylim=c(-0.3,0.6), legend.names = LR_ss$region)
legend(nr+0.5,0.5, legend = seq(nr), col = rainbow(nr), lty = 0, pch=1, lwd = 2)
dev.off()
print(length(dev.list()))
print("dimensions of dataset - train, test")
print(dim(train_fin))
print(dim(test_fin))

write.csv(LR_ss,paste0("LR_scores_9fold_",VT,"_LT_",LT,"npred_",length(varindex),".csv"))
write.csv(LR_ss2,paste0("LR_scores_fin_",VT,"_LT_",LT,"npred_",length(varindex),".csv"))
submodels = list()
for(i in seq(9*maxsteps+1,12*maxsteps)){
  submodels = append(submodels, list(models[[i]]$coefficients))
}
saveRDS(submodels, paste0("Models_",VT,"_LT_",LT,"_npred_",length(varindex)))

# --------------------------------------------------
test_that("Test dataset complete?", {
  expect_equal(filter(ObsPV, selector == 1) %>% arrange(Year, Month, Day), rbind(train_fin, test_fin) %>% arrange(Year, Month, Day))
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
