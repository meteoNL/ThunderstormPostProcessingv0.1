library(dplyr)

## READING ECMWF DATA ###
setwd("/usr/people/groote/ThunderstormPostProcessingv1")
setwd("/usr/people/groote/ThunderstormPostProcessingv1/ECWMF_Data2")

# get all the regions and lead times
xreg = c("W","M","E")
yreg = c("XN","MN","MS","XS")
fperiods = seq(9,27,3)
dataset = data.frame()

# loop over ECMWF data available
for (period in fperiods){
  region = 1
  for(yregion in yreg){
    for (xregion in xreg){
      name = sprintf("ECMO12_FP0%02d_%1s_%2s_1013.dat",period,xregion,yregion)
      print(name)
      first=read.csv(name, sep = "\t", header = TRUE)
      
      # read year, month, day from file
      year = round((first[1])/10000,0)
      month = round(round(first[1]-10000*round(first[1]/10000,0),0)/100,0)
      day = round(first[1]-100*round(first[1]/100),0)
      df= data.frame(Region = region, Year = year, Month = month, Day = day, VT = sprintf("VT_%02d%02d", (9+period)%%24, (15+period)%%24), leadtimeH = sprintf("%02d-%02d", (period-9), (period-3)), first)
      
      #adjust column names; leadtimeH is harmonie lead time belonging to the ECMWF file (6 hours less +/- 3 hours from the central point in interval)
      colnames(df) = c("region","Year","Month","Day", "VT", "leadtimeH", names(first))
      dataset = rbind(df,dataset)
      
      #go to next region
      region = region + 1
    }
  }
}

# same as above, but for different file names at different lead times
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
# id gives unique characteristics of the record: combining region, date, lead time and valid time
dataset2 = data.frame(dataset2, id = paste("id", dataset2$region, dataset2$Year, dataset2$Month, dataset2$Day, dataset2$VT, dataset2$leadtimeH))
dataset = data.frame(dataset, id = paste("id", dataset$region, dataset$Year, dataset$Month, dataset$Day, dataset$VT, dataset$leadtimeH))

## REPLACE MISSINGS WITH NA ####
dataset2[dataset2 == 9999] <- NA
dataset2[dataset2 == 99999] <- NA
dataset[dataset == 9999] <- NA
dataset[dataset == 99999] <- NA
setwd("/usr/people/groote")

### USE OLD DATA FRAME AND COMBINE WITH IDS ####
# read old dataframe from file and get similar ids as are above for the ECMWF files
Test = read.csv("Thunderstorm_radar_merged.csv")
Test$id = paste(Test$id, Test$leadtime)
Test$id = as.character(Test$id)
dataset2$id = as.character(dataset2$id)
dataset$id = as.character(dataset$id)
dataset = data.frame(dataset, leadtime = dataset$leadtimeH)
dataset2 = data.frame(dataset2, leadtime = dataset2$leadtimeH)
#and merge datasets
new_Obs <- left_join(Test, dataset, by = "id")
new_Obs <- left_join(new_Obs, dataset2, by = "id")

#### GENERATE SOME EXTRA PREDICTORS ####
# use PW
PW1 = log(new_Obs$PrecipitableWater_0mabovegnd_6hrlymax)
PW2 = new_Obs$PrecipitableWater_0mabovegnd_6hrlymax
#combine CAPE And CIN
CAPE_CIN_add = new_Obs$Surface_CAPE_35mabovegnd_6hrlymax+new_Obs$Surface_ConvInhib_0mabovegnd_6hrlymax
CAPE_CIN_add_sfc = new_Obs$Surface_CAPE_0mabovegnd_6hrlymax+new_Obs$Surface_ConvInhib_0mabovegnd_6hrlymax
#combine harmonie precip
hprecipsqrtmax = new_Obs$Rain_0mabovegnd_6hrlymax
#some stepwise predictors
stepwise_first = (new_Obs$modJefferson_0mabovegnd_6hrlymax > 32.5) + (PW2 > 44)
conditional_first = (new_Obs$modJefferson_0mabovegnd_6hrlymax > 32.5) * (PW2 > 44)
#transform Ri number
cnst_RItrans = 1e-5
RIindex=c(seq(36,37),seq(78,79))
new_Obs[RIindex] = log(cnst_RItrans+new_Obs[RIindex])

# some predictors and combined predictors in new dataframe
newpredictors = data.frame(
  Boyden_PW1 = (new_Obs$Boyden_0mabovegnd_6hrlymax-85)*PW1,
  Boyden_PW2 = (new_Obs$Boyden_0mabovegnd_6hrlymax-85)*PW2,
  Boyden_hprecipmax = (new_Obs$Boyden_0mabovegnd_6hrlymax-85)*hprecipsqrtmax,
  Bradbury_PW1 = (new_Obs$Bradbury_0mabovegnd_6hrlymin-14)*PW1,
  Bradbury_PW2 = (new_Obs$Bradbury_0mabovegnd_6hrlymin-14)*PW2,
  Bradbury_hprecipmax = (new_Obs$Bradbury_0mabovegnd_6hrlymin-14)*hprecipsqrtmax,
  Edward_PW1 = (new_Obs$theataW_925mb_6hrlymax-new_Obs$theataW_500mb_6hrlymin)*PW1,
  Edward_PW2 = (new_Obs$theataW_925mb_6hrlymax-new_Obs$theataW_500mb_6hrlymin)*PW2,
  Edward_hprecipmax = (new_Obs$theataW_925mb_6hrlymax-new_Obs$theataW_500mb_6hrlymin)*hprecipsqrtmax,
  LidS_PW = log(0.01+new_Obs$LidStrength_0mabovegnd_6hrlymin*PW2),
  LidS_hprecipmax = log(0.01+new_Obs$LidStrength_0mabovegnd_6hrlymin*hprecipsqrtmax),
  Lifted_PW1 = PW1/(new_Obs$Lifted_0mabovegnd_6hrlymin+6),
  Lifted_PW2 = PW2/(new_Obs$Lifted_0mabovegnd_6hrlymin+6),
  Lifted_hprecipmax = hprecipsqrtmax/(new_Obs$Lifted_0mabovegnd_6hrlymin+6),
  binary_1 = as.numeric((new_Obs$Boyden_0mabovegnd_6hrlymax > 93) *(new_Obs$PrecipitableWater_0mabovegnd_6hrlymax > 10) *(new_Obs$LidStrength_0mabovegnd_6hrlymin < 3) ),
  binary_2 = as.numeric((new_Obs$Boyden_0mabovegnd_6hrlymax > 93) *(new_Obs$PrecipitableWater_0mabovegnd_6hrlymax > 10) *(new_Obs$LidStrength_0mabovegnd_6hrlymin < 3)*(new_Obs$modJefferson_0mabovegnd_6hrlymax > 20)*(new_Obs$TotalTotals_2e_0mabovegnd_6hrlymax > 20) ),
  cos_WDIR_500_max = cos(new_Obs$WDIR_500mb_6hrlymax),
  sin_WDIR_500_max = sin(new_Obs$WDIR_500mb_6hrlymax),
  cos_WDIR_850_max = cos(new_Obs$WDIR_850mb_6hrlymax),
  sin_WDIR_850_max = sin(new_Obs$WDIR_850mb_6hrlymax),
  cos_WDIR_500_min = cos(new_Obs$WDIR_500mb_6hrlymin),
  sin_WDIR_500_min = sin(new_Obs$WDIR_500mb_6hrlymin),
  cos_WDIR_850_min = cos(new_Obs$WDIR_850mb_6hrlymin),
  sin_WDIR_850_min = sin(new_Obs$WDIR_850mb_6hrlymin),
  sfccape_pow_0.2max = new_Obs$Surface_CAPE_0mabovegnd_6hrlymax^0.2,
  mucape_pow_0.2max = new_Obs$Surface_CAPE_35mabovegnd_6hrlymax^0.2,
  sfccape_pow_0.2min = new_Obs$Surface_CAPE_0mabovegnd_6hrlymin^0.2,
  mucape_pow_0.2min = new_Obs$Surface_CAPE_35mabovegnd_6hrlymin^0.2,
  CAPE_CIN_addmax = CAPE_CIN_add,
  CAPE_CIN_addmax_sfc = CAPE_CIN_add_sfc,
  CAPE_CIN_pow_0.5max <- ifelse(CAPE_CIN_add<0,0,CAPE_CIN_add^0.5),
  CAPE_CIN_pow_0.5sfc <- ifelse(CAPE_CIN_add_sfc<0,0,CAPE_CIN_add_sfc^0.5),
  SWEAT_sfc_max_pow0.2 = new_Obs$SWEAT_0mabovegnd_6hrlymax^0.2,
  abs_helicity_pow0.1_max = abs(new_Obs$Helicity_0mabovegnd_6hrlymax)^0.1,
  abs_helicity_pow0.1_min = abs(new_Obs$Helicity_0mabovegnd_6hrlymin)^0.1,
  mod_Jeff_pow0.2_max_0replace <- ifelse(new_Obs$modJefferson_0mabovegnd_6hrlymax<0,0,new_Obs$modJefferson_0mabovegnd_6hrlymax^0.2),
  mod_Jeff_pow0.2_min_0replace <- ifelse(new_Obs$modJefferson_0mabovegnd_6hrlymin<0,0,new_Obs$modJefferson_0mabovegnd_6hrlymin^0.2),
  logPWmax = PW1,
  logPWmin = log(new_Obs$PrecipitableWater_0mabovegnd_6hrlymin),
  Theta_w850_max_pow0.2 = new_Obs$theataW_850mb_6hrlymax^0.2,
  Theta_w850_min_pow0.2 = new_Obs$theataW_850mb_6hrlymin^0.2,
  ff300.x_pow0.2 = new_Obs$ff300.x^0.2,
  ff300.y_pow0.2 = new_Obs$ff300.y^0.2,
  ta5010.x_pow2 = new_Obs$ta5010.x^2,
  ta5010.y_pow2 = new_Obs$ta5010.y^2,
  stepwise_Inst_PW = stepwise_first,
  stepwise_Inst_PW_rain = stepwise_first + (new_Obs$Rain_0mabovegnd_6hrlymax > 5), 
  conditional_Inst_PW = conditional_first,
  conditional_Inst_PW_rain = conditional_first * (new_Obs$Rain_0mabovegnd_6hrlymax > 5),
  new_Obs[RIindex]
  )

# merge with old dataframe
new_Obs = cbind(new_Obs, newpredictors)
names(new_Obs)[241:244]=paste0(names(new_Obs)[241:244],"_transformedlog+cnst")

###### REMOVE WRONG VALUES ######
new_Obs <- filter(new_Obs, LCL_0mabovegnd_6hrlymax < 12500)
new_Obs = new_Obs[2:length(new_Obs)]
colnames(new_Obs)[1:10]=c("id","region","Year","Month","Day","validtime","radarmax","validdate","T0","leadtime")

## WRITE NEW DATASET TO FILE #### 
write.csv(new_Obs, "ECMWF_merged4.csv")
