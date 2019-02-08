library(dplyr)

## READING ECMWF DATA ###
setwd("/usr/people/groote/ThunderstormPostProcessingv1")
setwd("/usr/people/groote/ThunderstormPostProcessingv1/ECWMF_Data2")
xreg = c("W","M","E")
yreg = c("XN","MN","MS","XS")
fperiods = seq(9,27,3)
dataset = data.frame()
for (period in fperiods){
  region = 1
  for(yregion in yreg){
    for (xregion in xreg){
      name = sprintf("ECMO12_FP0%02d_%1s_%2s_1013.dat",period,xregion,yregion)
      print(name)
      #print(name)
      first=read.csv(name, sep = "\t", header = TRUE)
      year = round((first[1])/10000,0)
      month = round(round(first[1]-10000*round(first[1]/10000,0),0)/100,0)
      day = round(first[1]-100*round(first[1]/100),0)
      df= data.frame(Region = region, Year = year, Month = month, Day = day, VT = sprintf("VT_%02d%02d", (9+period)%%24, (15+period)%%24), leadtimeH = sprintf("%02d-%02d", (period-9), (period-3)), first)
      colnames(df) = c("region","Year","Month","Day", "VT", "leadtimeH", names(first))
      dataset = rbind(df,dataset)
      region = region + 1
    }
  }
}
fperiods = seq(30,39,3)
dataset2 = data.frame()
for (period in fperiods){
  region = 1
  for(yregion in yreg){
    for (xregion in xreg){
      name = sprintf("ECMO12_FP0%02d_%1s_%2s_1013.dat",period,xregion,yregion)
      print(name)
      #print(name)
      first=read.csv(name, sep = "\t", header = TRUE)
      year = round((first[1])/10000,0)
      month = round(round(first[1]-10000*round(first[1]/10000,0),0)/100,0)
      day = round(first[1]-100*round(first[1]/100),0)
      df= data.frame(Region = region, Year = year, Month = month, Day = day, VT = sprintf("VT_%02d%02d", (9+period)%%24, (15+period)%%24), leadtimeH = sprintf("%02d-%02d", (period-9), (period-3)), first)
      colnames(df) = c("region","Year","Month","Day", "VT", "leadtimeH", names(first))
      dataset2 = rbind(df,dataset2)
      region = region + 1
    }
  }
}

## GENERATING IDS FOR COMBINING DATAFRAMES ####
dataset2 = data.frame(dataset2, id = paste("id", dataset2$region, dataset2$Year, dataset2$Month, dataset2$Day, dataset2$VT, dataset2$leadtimeH))
dataset = data.frame(dataset, id = paste("id", dataset$region, dataset$Year, dataset$Month, dataset$Day, dataset$VT, dataset$leadtimeH))

## REPLACE MISSINGS WITH NA ####
dataset2[dataset2 == 9999] <- NA
dataset2[dataset2 == 99999] <- NA
dataset[dataset == 9999] <- NA
dataset[dataset == 99999] <- NA
setwd("/usr/people/groote")

### USE OLD DATA FRAME AND COMBINE WITH IDS ####
Test = read.csv("Thunderstorm_radar_merged.csv")
Test$id = paste(Test$id, Test$leadtime)
Test$id = as.character(Test$id)
dataset2$id = as.character(dataset2$id)
dataset$id = as.character(dataset$id)
dataset = data.frame(dataset, leadtime = dataset$leadtimeH)
dataset2 = data.frame(dataset2, leadtime = dataset2$leadtimeH)
new_Obs <- left_join(Test, dataset, by = "id")
new_Obs <- left_join(new_Obs, dataset2, by = "id")

#### GENERATE SOME EXTRA PREDICTORS ####
PW1 = log(new_Obs$PrecipitableWater_0mabovegnd_6hrlymax)
PW2 = new_Obs$PrecipitableWater_0mabovegnd_6hrlymax/4.
newpredictors = data.frame(
  Boyden_PW1 = (new_Obs$Boyden_0mabovegnd_6hrlymax-85)*PW1,
  Boyden_PW2 = (new_Obs$Boyden_0mabovegnd_6hrlymax-85)*PW2,
  Bradbury_PW1 = (new_Obs$Bradbury_0mabovegnd_6hrlymin-14)*PW1,
  Bradbury_PW2 = (new_Obs$Bradbury_0mabovegnd_6hrlymin-14)*PW2,
  Edward_PW1 = (new_Obs$theataW_925mb_6hrlymax-new_Obs$theataW_500mb_6hrlymin)*PW1,
  Edward_PW2 = (new_Obs$theataW_925mb_6hrlymax-new_Obs$theataW_500mb_6hrlymin)*PW2,
  LidS_PW = log(0.01+new_Obs$LidStrength_0mabovegnd_6hrlymin*PW2),
  Lifted_PW1 = PW1/(new_Obs$Lifted_0mabovegnd_6hrlymin+6),
  Lifted_PW2 = PW2/(new_Obs$Lifted_0mabovegnd_6hrlymin+6),
  binary_1 = as.numeric((new_Obs$Boyden_0mabovegnd_6hrlymax > 93) *(new_Obs$PrecipitableWater_0mabovegnd_6hrlymax > 10) *(new_Obs$LidStrength_0mabovegnd_6hrlymin < 3) ),
  binary_2 = as.numeric((new_Obs$Boyden_0mabovegnd_6hrlymax > 93) *(new_Obs$PrecipitableWater_0mabovegnd_6hrlymax > 10) *(new_Obs$LidStrength_0mabovegnd_6hrlymin < 3)*(new_Obs$modJefferson_0mabovegnd_6hrlymax > 20)*(new_Obs$TotalTotals_2e_0mabovegnd_6hrlymax > 20) ),
  cos_WDIR_500_max = cos(new_Obs$WDIR_500mb_6hrlymax),
  sin_WDIR_500_max = sin(new_Obs$WDIR_500mb_6hrlymax),
  cos_WDIR_850_max = cos(new_Obs$WDIR_850mb_6hrlymax),
  sin_WDIR_850_max = sin(new_Obs$WDIR_850mb_6hrlymax),
  cos_WDIR_500_min = cos(new_Obs$WDIR_500mb_6hrlymin),
  sin_WDIR_500_min = sin(new_Obs$WDIR_500mb_6hrlymin),
  cos_WDIR_850_min = cos(new_Obs$WDIR_850mb_6hrlymin),
  sin_WDIR_850_min = sin(new_Obs$WDIR_850mb_6hrlymin))

new_Obs = cbind(new_Obs, newpredictors)

## WRITE NEW DATASET #### 
write.csv(new_Obs, "ECMWF_merged3.csv")