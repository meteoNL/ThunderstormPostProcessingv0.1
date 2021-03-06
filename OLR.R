########################################
#This model type was used to assess the linearity of transformed threshold with power (0.25)

### LOAD EXTERNAL CODE
library(dplyr)
library(MASS)
library(verification)
library(ordinal)
library(arm)
library(crch)
library(profvis)
#rm(list=ls(all=TRUE))
### SET GENERAL CONDITIONS FOR THE MODEL
#set thresholds and hyperparameter; determine test dataset and training dataset
p = 0.25 #power transformation to linearize thresholds
maxvars = 4
numsubset = 3 #number of subsets for hyperparameter selection
thres = seq(5,70,5)
minpredictant = 1.5 #minimum precipitation sum considered as precipitation
ObsPV = read.csv(file = "Thunderstorm_radar_merged.csv")
years = c(as.numeric(unique(ObsPV$Year)))
LT = c(as.numeric(unique(ObsPV$leadtime_count)))[1]
VT = unique(ObsPV$validtime.x)[2]
regions = c(unique(ObsPV$region))#[1:2]

#change radarmax into other name if necessary:
climset = filter(ObsPV, Ndischarge > minpredictant & validtime.x == VT & leadtime_count == LT)
predictant = 113

#selects the best predictor by forward fittng, trying potential predictors and selecting the one with lowest AIC
observed = seq(length(climset))*0
for(i in seq(length(thres))){
  threshold = thres[i]
  observed = observed+as.numeric(climset[predictant]>threshold)
  print(observed)
}

model = clm(as.factor(observed)~ .,
               data=data.frame(climset[32]), link = "logit", start = c((thres^p),-0.05))
plot(model$coefficients[1:(length(model$coefficients)-1)], thres^p)
