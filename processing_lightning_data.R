first=read.csv("processed_lightning_dataset/kldn_nil_2010_E_MN", sep = " ", header = FALSE)

#define available regions matrix (KOUW regions) and years
xreg = c("W","M","E")
yreg = c("XN","MN","MS","XS")
years = c(2015,2016,2017)
dataset = data.frame()
for (y in years){
  region = 1
  for(yregion in yreg){
    for (xregion in xreg){
      name = paste("processed_lightning_dataset/regions00061218/kldn_nil_",y,"_",xregion,"_",yregion,sep ="")
      #print(name)
      #read relevant columns from file: year, month day, hour, number of discharges
      first=read.csv(name, sep = " ", header = FALSE)
      year = round(first[1]/1000000,0)
      month = round((first[1]-1000000*round(first[1]/1000000))/10000,0)
      day = round(round(first[1]-10000*round(first[1]/10000,0),0)/100,0)
      hour = round(first[1]-100*round(first[1]/100),0)
      ndischarges=first[2]
      
      #make dataframe and manipulate its column names
      df= data.frame(region, year, month, day, sprintf("VT_%02d%02d",unlist(c((hour-3)%%24)),unlist(c((hour+3)%%24))), ndischarges)
      colnames(df) = c("region","Year","Month","Day","validtime","Ndischarge")
      dataset = rbind(df,dataset)
      
      #update region
      region = region + 1
    }
  }
}

# write csv of dataset with number of discharges and using an id in data frame to merge data frames later
dataset = data.frame(id = paste("id",dataset$region,dataset$Year,dataset$Month,dataset$Day,dataset$validtime),dataset)
write.csv(dataset, "Ndischarges_dataset.csv")

#apply the procedure once more for discharge rates
maxdataset=data.frame()
for (y in years){
  region = 1
  for(yregion in yreg){
    for (xregion in xreg){
      name = paste("processed_lightning_dataset/regions00061218/kldn_nil_m5mi2_",y,"_",xregion,"_",yregion,sep ="")
      #print(name)
      first=read.csv(name, sep = " ", header = FALSE)
      year = round(first[1]/1000000,0)
      month = round((first[1]-1000000*round(first[1]/1000000))/10000,0)
      day = round(round(first[1]-10000*round(first[1]/10000,0),0)/100,0)
      hour = round(first[1]-100*round(first[1]/100),0)
      hour_re = 
      ndischarges=first[2]
      df= data.frame(region, year, month, day, sprintf("VT_%02d%02d",unlist(c((hour-3)%%24)),unlist(c((hour+3)%%24))), ndischarges)
      colnames(df) = c("region","Year","Month","Day","validtime","Dischargerate")
      maxdataset = rbind(df,maxdataset)
      region = region + 1
    }
  }
}

# replace missing/wrong values with NA
maxdataset[maxdataset == 99999] <- NA

# write csv of dataset with number of discharges and using an id in data frame to merge data frames later
maxdataset = data.frame(id = paste("id",dataset$region,dataset$Year,dataset$Month,dataset$Day,dataset$validtime),maxdataset)
write.csv(maxdataset, "Maxdischarges_dataset.csv")
library(dplyr)

#merge with predictors data frame using id
ObsPV = read.csv("Final_predictors_00zset.csv")
year = round(ObsPV$Init_date/10000)
dates=sort(unique(ObsPV$Init_date%%10000))
dates =c(dates, dates[length(dates)]+1)
matrix_one = matrix(rep(dates, each = nrow(ObsPV)),ncol=length(dates))
matrix_two = matrix(rep(seq(1,length(dates)),each = nrow(ObsPV)),ncol=length(dates))
matrix_three = matrix(rep((ObsPV$Init_date%%10000),length(dates)),ncol=length(dates))
matrix_four = as.numeric(matrix_one==matrix_three)*matrix_two
new=ObsPV$Init_date*0
for(i in seq(1, nrow(ObsPV))){
  new[i]=max(matrix_four[i,])
}
# change this line for the other runs!!! 
#new = new+as.numeric(ObsPV$Forecast_hours_mean>23.999)
print("change the forecast_hours_mean_threshold!!!!!!!")
print(somewhere)
valid_dates = dates[new]
valid_time = sprintf("VT_%02d%02d",ObsPV$Forecast_hours_min%%24,ObsPV$Forecast_hours_max%%24)
ObsPV = cbind(Year = year, Valid_date = valid_dates, Month = round(valid_dates/100), Day = valid_dates%%100, validtime = valid_time, ObsPV)
ObsPV = data.frame(id = paste("id",ObsPV$Region,ObsPV$Year,as.numeric(ObsPV$Month),as.numeric(ObsPV$Day),ObsPV$validtime),ObsPV)
ObsPV$id = as.character(ObsPV$id)
dataset$id = as.character(dataset$id)
maxdataset$id = as.character(maxdataset$id)
new_Obs <- left_join(ObsPV, dataset, by = "id")
new_Obs <- left_join(new_Obs, maxdataset, by = "id")

#write the result to a new file
write.csv(new_Obs, "full_final00z_dataset.csv")
