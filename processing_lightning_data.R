first=read.csv("processed_lightning_dataset/kldn_nil2_2010_E_MN", sep = " ", header = FALSE)
xreg = c("W","M","E")
yreg = c("XN","MN","MS","XS")
years = c(2010,2011,2013)
dataset = data.frame()
for (y in years){
  region = 1
  for(yregion in yreg){
    for (xregion in xreg){
      name = paste("processed_lightning_dataset/kldn_nil2_",y,"_",xregion,"_",yregion,sep ="")
      #print(name)
      first=read.csv(name, sep = " ", header = FALSE)
      year = round(first[1]/1000000,0)
      month = round((first[1]-1000000*round(first[1]/1000000))/10000,0)
      day = round(round(first[1]-10000*round(first[1]/10000,0),0)/100,0)
      hour = round(first[1]-100*round(first[1]/100),0)
      ndischarges=first[2]
      df= data.frame(region, year, month, day, sprintf("VT_%02d%02d",unlist(c(hour-3)),unlist(c((hour+3)%%24))), ndischarges)
      colnames(df) = c("region","Year","Month","Day","validtime","Ndischarge")
      dataset = rbind(df,dataset)
      region = region + 1
    }
  }
}
dataset = data.frame(id = paste("id",dataset$region,dataset$Year,dataset$Month,dataset$Day,dataset$validtime),dataset)
write.csv(dataset, "Ndischarges_dataset.csv")
maxdataset=data.frame()
for (y in years){
  region = 1
  for(yregion in yreg){
    for (xregion in xreg){
      name = paste("processed_lightning_dataset/kldn_nil2_m5mi2_",y,"_",xregion,"_",yregion,sep ="")
      #print(name)
      first=read.csv(name, sep = " ", header = FALSE)
      year = round(first[1]/1000000,0)
      month = round((first[1]-1000000*round(first[1]/1000000))/10000,0)
      day = round(round(first[1]-10000*round(first[1]/10000,0),0)/100,0)
      hour = round(first[1]-100*round(first[1]/100),0)
      hour_re = 
      ndischarges=first[2]
      df= data.frame(region, year, month, day, sprintf("VT_%02d%02d",unlist(c(hour-3)),unlist(c((hour+3)%%24))), ndischarges)
      colnames(df) = c("region","Year","Month","Day","validtime","Dischargerate")
      maxdataset = rbind(df,maxdataset)
      region = region + 1
    }
  }
}
maxdataset[maxdataset == 99999] <- NA
maxdataset = data.frame(id = paste("id",dataset$region,dataset$Year,dataset$Month,dataset$Day,dataset$validtime),maxdataset)
write.csv(maxdataset, "Maxdischarges_dataset.csv")
library(dplyr)
ObsPV = readRDS(file = "/usr/people/whan/Research/Whanetal_HarmoniePr_2017/data/ObsPV.rds")
ObsPV = data.frame(id = paste("id",ObsPV$region,ObsPV$Year,as.numeric(ObsPV$Month),as.numeric(ObsPV$Day),ObsPV$validtime),ObsPV)
ObsPV$id = as.character(ObsPV$id)
dataset$id = as.character(dataset$id)
maxdataset$id = as.character(maxdataset$id)
new_Obs <- left_join(ObsPV, dataset, by = "id")
new_Obs <- left_join(new_Obs, maxdataset, by = "id")
write.csv(new_Obs, "Thunderstorm_radar_merged.csv")