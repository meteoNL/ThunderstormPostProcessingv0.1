library(parallel)
#library(foreach)
cls <- makeCluster(getOption("cl.cores", 8))
#Rlib <- "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3/"
#.libPaths(c(.libPaths(), Rlib))
#setwd("/net/pc150388/nobackup_1/users/stat/para/INDECS/HA40REFC/GWI/")
#setwd("/net/pc150388/nobackup_1/users/stat/para/INDECS/HA40REFC/THUNFC/")

# load the mepsr library from Kiri's directory
library(meteogrid, lib.loc = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3/")
library(Rgrib2, lib.loc = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3/")
library(mepsr, lib.loc = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3/")
library(raster, lib.loc = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3/")
library(dplyr, lib.loc = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3/")
library(shiny, lib.loc = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3/")
library(magrittr, lib.loc = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3/")
# get KOUW raster:
setwd("/usr/people/groote")
#ha40_file <- "testgrib" #path to the file you want (each file contains one date, one leadtime and all variables)

# define generic domain information:
#h40_grib_info <- Gopen(filename = ha40_file)
datelist = c(rep(4,16),rep(5,31),rep(6,30),rep(7,31),rep(8,31),rep(9,30),rep(10,15))*100+c(seq(15,30),seq(1,31),seq(1,30),seq(1,31),seq(1,31),seq(1,30),seq(1,15))
datelist = rep(datelist,3)+c(rep(2015,length(datelist))*10000,rep(2016,length(datelist))*10000,rep(2017,length(datelist))*10000)
rasterKOUW = readRDS("harmonie_mask.rds")

#define column names for the data frame with GRIB data
column_names <- c("Forecast_hours","Region","PrecipWater","Conv_cloud_cov","MSLP","Cloud_base","Graupel_col",
                  "Lightning_H40","Cloud_top","Cloud_ice","Cloud_water","Rain","Snow","Vertical_velocity",
                  "DPT_500","DPT_600","DPT_700","DPT_850","WDIR_500","WDIR_850","WSPD_500","WSPD_850",
                  "U_700","V_700","ThetaW_500","ThetaW_850","ThetaW_925","ThetaWs_500",
                  "Boyden","Bradbury","Showalter","Rackliff","Jefferson","ModJefferson","K_index","Fateev",
                  "Totals_Totals","Vertical_Totals","Cross_Totals","SWEAT","sin_dd_dif_500_850","TQ",
                  "LiftedIndex","SBCAPE","MUCAPE","SBCIN","MUCIN","LFC","LNB","Shear","Ri_nr_lvl1","Ri_nr_lvl2",
                  "Storm_Travel","Helicity","Lid_Strength","PrecipWater_acc","Rain_acc","Cloud_top2","Graupel_col2","Moisture_convergence")

rotate <- function(x) t(apply(x, 2, rev))

# collect all data from the two gribfiles in THUNFC & GWI
readarun <- function(Lts){
  library(meteogrid, lib.loc = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3/")
  library(Rgrib2, lib.loc = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3/")
  library(mepsr, lib.loc = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3/")
  library(raster, lib.loc = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3/")
  library(dplyr, lib.loc = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3/")
  library(shiny, lib.loc = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3/")
  library(magrittr, lib.loc = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3/")
  date=datum
  max_level = 14000
#  message(date)
  no_cloud_value = -10
  cloud_base_classifications_lower = c(0,2500, 5500)
  cloud_base_classifications_upper = c(2500, 5500,15000)
  overall_frame = data.frame()
  ind_thunfc=c(1,11,seq(17,25),2)
  for(L in Lts*6){
    print(1)
    setwd("/net/pc150388/nobackup_1/users/stat/para/INDECS/HA40REFC/")
    all_data_run = data.frame()
    for(i in seq(0,6)){
      i=i+L
      #        print(i)
      #        # read the files for different forecasting hours & put them on the KOUW grid for THUNFC
      dataset = data.frame(as.vector(rasterKOUW))
      fn1 = sprintf("THUNFC/HA40_N25_THUNFC_%08d0600_0%02d00_GB",date,i)
      # 
      #        #domain data; are assumed to be always the same in both files
      domain_data <- meteogrid::DomainExtent(Gdec(Ghandle(fn1, 1)))
      if(nrow(Gopen(filename = fn1))==25){
        for(j in ind_thunfc){
          geofield_data <- Gdec(Ghandle(fn1, j))
          tmp <- raster(t(rotate(t(geofield_data))), xmn=domain_data$x0, xmx=domain_data$x1, ymn=domain_data$y0, ymx=domain_data$y1)
          column = resample(tmp, rasterKOUW, method = "ngb")
          #            # print(column)
          # 
          #            #combine with dataset
          dataset = cbind(dataset, as.vector(column))
        }
      } else{
        dataset = cbind(dataset, data.frame(matrix(NA, ncol = 12, nrow = nrow(dataset))))
      }
      #        # read the files for different forecasting hours & put them on the KOUW grid for GWI
      fn2 = sprintf("GWI/HA40_GWI_%08d060000_%02d00_GB",date,i)
      if(nrow(Gopen(filename = fn2))==46){
        for(j in seq(1,46)){
          #print(getwd())
          geofield_data <- Gdec(Ghandle(fn2, j))
          tmp <- raster(t(rotate(t(geofield_data))), xmn=domain_data$x0, xmx=domain_data$x1, ymn=domain_data$y0, ymx=domain_data$y1)
          column = resample(tmp, rasterKOUW, method = "ngb")
          if(j==35){
            column[column>max_level]<-NA
          }
          if(j==34){
            column[column>max_level]<-NA
          }
          
          #combine with dataset
          dataset = cbind(dataset, as.vector(column))
        }
      } else{
        dataset = cbind(dataset, data.frame(matrix(NA, ncol = 46, nrow = nrow(dataset))))
      }
      
      #        #add forecasting hour and declare column names
      dataset = cbind(FT = i, dataset)
      colnames(dataset) <- column_names
      #   
      #        #add extra predictors on KOUW grid
      extra_predictors = data.frame(Bradbury_base925 = dataset$ThetaW_925-dataset$ThetaW_500,ThetaW_925_850_diff = dataset$ThetaW_925-dataset$ThetaW_850, U_500 = -sin(dataset$WDIR_500/360*2*pi)*dataset$WSPD_500, V_500 = -cos(dataset$WDIR_500/360*2*pi)*dataset$WSPD_500, U_850 = -sin(dataset$WDIR_850/360*2*pi)*dataset$WSPD_850, V_850 = -cos(dataset$WDIR_850/360*2*pi)*dataset$WSPD_850,
                                    MUCAPE_graupel_pows = (dataset$MUCAPE^0.25)*(dataset$Graupel_col^0.1), MUCAPE_snow_pows = (dataset$MUCAPE^0.25)*(dataset$Snow^(1.0/6)), MUCAPE_water_pows = (dataset$MUCAPE^0.25)*(dataset$Cloud_water^0.10),ifelse((dataset$SBCAPE+dataset$SBCIN)<0,0,(dataset$SBCAPE+dataset$SBCIN)^0.5),ifelse((dataset$MUCAPE+dataset$MUCIN)<0,0,(dataset$MUCAPE+dataset$MUCIN)^0.5))
      extra_predictors = cbind(extra_predictors, data.frame(Bulk_shear_850_500 = ((extra_predictors$U_850-extra_predictors$U_500)^2 + (extra_predictors$V_850-extra_predictors$V_500)^2)^0.5,
                                                            Bulk_shear_850_700 = ((extra_predictors$U_850-dataset$U_700)^2 + (extra_predictors$V_850-dataset$V_700)^2)^0.5,
                                                            Cloud_layers_depth = dataset$Cloud_top-dataset$Cloud_base))
      dataset = cbind(dataset, extra_predictors)
      #   
      #        #maybe: here apply filter: if LFC/LNB/cloud_base/cloud_top > 12.500 m, replace by another value
      #   
      #        #combine with previous hours
      all_data_run = rbind(all_data_run, dataset)
    }
    #   
    #      # mean KOUW grid cell pressure tendency and Boyden tendency
    pressure_data <- all_data_run %>% group_by(Region, Forecast_hours) %>% summarise_at(c("MSLP","Boyden"),mean)
    p_of_t = t(matrix(pressure_data$MSLP,ncol=length(unique(pressure_data$Region))))
    Boyden_of_t = t(matrix(pressure_data$Boyden,ncol=length(unique(pressure_data$Region))))
    dpdt = p_of_t[,ncol(p_of_t)]-p_of_t[,1]
    dBoydendt = Boyden_of_t[,ncol(Boyden_of_t)]-Boyden_of_t[,1]
    mean_cc <- all_data_run %>% group_by(Region) %>% summarise(1-sum(is.na(Cloud_base))/length(Cloud_base))
    
    quantilesperregion = data.frame()
    for(regionsel in unique(all_data_run$Region)){
      #        #get all quantiles and subsequently mean, per region filter
      test <- filter(all_data_run, Region == regionsel)
      test2 <- sapply(test, quantile, probs = c(0.0,0.02,0.1,0.5,0.9,0.98,1.00), na.rm = TRUE)
      test3 <- sapply(test, mean, na.rm = TRUE)
      test2 = rbind(test3,test2)
      
      #        #ad to previous quantiles
      quantilesperregion = rbind(quantilesperregion, data.frame(Region = regionsel, t(as.vector(test2))))
    }
    quantilesperregion = cbind(quantilesperregion, dp_dt = dpdt, dBoyden_dt = dBoydendt, mean_cloud_cover = mean_cc[2])
    
    colnames(quantilesperregion) <- c("Region",paste0(rep(names(all_data_run[1:length(all_data_run)]),each=8),rep(c("_mean","_min","_q0.02","_q0.10","_median","_q0.90","_q0.98","_max"),length(all_data_run))),
                                      "dp_dt","dBoyden_dt","mean_cloud_cover")
    # cloud base categorial data frame
    cloud_base_df = data.frame()
    cloud_base_vector = c()
    for(pred in t(matrix(c(quantilesperregion$Cloud_base_q0.10, quantilesperregion$Cloud_base_median, quantilesperregion$Cloud_base_q0.90),ncol=3,nrow=13))){
      no1 = as.numeric(pred<cloud_base_classifications_upper)
      no2 = no1 * as.numeric(1-(pred<cloud_base_classifications_lower))
      no3 = c(no2, 0)
      if(max(no2, na.rm = TRUE)==-Inf){
        no3 = c(0,0,0,1)
      }
      cloud_base_vector = c(cloud_base_vector, no3)
    }
    cloud_base_matrix = t(matrix(cloud_base_vector, ncol = length(unique(quantilesperregion$Region))))
    colnrs_base = c(seq(1,3),seq(5,7),seq(9,12))
    #      #print(data.frame(cloud_base_matrix)[colnrs_base])
    quantilesperregion = cbind(quantilesperregion, data.frame(cloud_base_matrix)[colnrs_base])
    #      #give names to the potential predictors
    colnames(quantilesperregion) <- c(names(quantilesperregion[1:(length(quantilesperregion)-10)]),
                                      "Low_q0.1","Middle_q0.1","High_q0.1",
                                      "Low_median","Middle_median","High_median",
                                      "Low_q0.9","Middle_q0.9","High_q0.9","No_clouds")
    #   
    #      # replace no cloud cases with negative values for cloud base, top and depth to distinguish
    Delta_MSLP_maxmin = quantilesperregion$MSLP_max - quantilesperregion$MSLP_min
    Delta_MSLP_quan = quantilesperregion$MSLP_q0.98 - quantilesperregion$MSLP_q0.02
    quantilesperregion = cbind(quantilesperregion, Delta_MSLP_maxmin, Delta_MSLP_quan)
    quantilesperregion[quantilesperregion$No_clouds==TRUE,c(seq(42,49),seq(66,73),seq(458,465),seq(586,593))]<- no_cloud_value
    quantilesperregion = quantilesperregion[c(seq(1,529),537,545,553,561,seq(567,569),seq(570,length(quantilesperregion)))]
    print(names(quantilesperregion))
    quantilesperregion = cbind(Init_date = date, quantilesperregion)
    #      #print(quantilesperregion[seq(1,40)])
    overall_frame = rbind(overall_frame, quantilesperregion)
    #      #print(overall_frame[seq(1,40)])
    #    write.csv(overall_frame,"/usr/people/groote/Dataset_00z_predictorspar.csv")
  }
  return(overall_frame)
  
}
readgrib <- function(i){
  assign("datelist", datelist, .GlobalEnv)
  print(datelist[i])
  #    print(datelist)
  #    print(dim(rasterKOUW))
  datum=datelist[i]
  assign("datum",datum, envir = .GlobalEnv)
  set_of_date=data.frame()
  clusterExport(cl = cls, c("column_names", "readgrib","rotate","rasterKOUW","readarun","datum"))
  result = parLapply(cl=cls, X = 0.5:6.5, fun = readarun)
  set_of_date <- Reduce(rbind, result)
  #result_df <- Reduce(rbind, result)
  
  #    print(date)
  
  return(set_of_date)
}  

all_data = data.frame()
for(i in seq(518,length(datelist))){
  print(i)
  all_data = rbind(readgrib(i),all_data)
  write.csv(all_data,"Dataset_06z_predictors_lastdays.csv")
  
}
# testpar <- function(x){
#   # assign("datelist", datelist, .GlobalEnv)
#   # print(overall_frame)
#   # print(as.vector(x))
#   # print(datelist[x])
#   # print("---")
#   return(as.POSIXct(as.character(x), format = "%Y%m%d"))
# }
# clusterExport(cl = cls, c("datelist","overall_frame","column_names", "readgrib","rotate","rasterKOUW"))
#clusterApply(cl = cls,x = 1:40, fun = testpar)
# for(p in seq(0,8)){
#   print(p)
#   sequence = seq((p*60+10),((p+1)*60+9))
#   result=parLapply(cl=cls, X = sequence, fun = readgrib)
#   result_df <- Reduce(rbind, result)
#   name = sprintf("/usr/people/groote/Dataset_00z_predictorspar%01d.csv",p)
#   print(name)
#   #write.csv(result_df,name)
# }

#result = parLapply(cl=cls, X = 549:552, fun = readgrib)
#result_df <- Reduce(rbind, result)
#name = sprintf("/usr/people/groote/Dataset_00z_predictorspar%01d.csv",9)
#result = parLapply(cl=cls, X = 1:8, fun = readgrib)
#result_df <- Reduce(rbind, result)
#name = sprintf("/usr/people/groote/Dataset_00z_predictorspar%01d.csv",10)
stopCluster(cls)