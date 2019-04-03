library(dplyr)
new_Obs = read.csv("Correct_00z_dataset.csv")
new_Obs <- filter(new_Obs, Region < 25) #just to remove NA from regions

#### GENERATE SOME EXTRA PREDICTORS ####
# use PW
PW1 = log(new_Obs$PrecipWater_max)
PW2 = new_Obs$PrecipWater_max
#combine CAPE And CIN
#CAPE_CIN_add_max = new_Obs$MUCAPE_max+new_Obs$MUCIN.or.Unkown_min
#CAPE_CIN_add_meanCIN = new_Obs$MUCAPE_max+new_Obs$MUCIN.or.Unkown_mean
#CAPE_CIN_add_sfcmax = new_Obs$SBCAPE_max+new_Obs$SBCIN_min
#combine harmonie precip
hprecipcubertmax = new_Obs$Rain_acc_max^(1./3)
#some stepwise predictors
stepwise_first = (new_Obs$ModJefferson_max > 32.5) + (PW2 > 44)
conditional_first = (new_Obs$ModJefferson_max > 32.5) * (PW2 > 44)
#transform Ri number
cnst_RItrans = 1e-5
RIindex=c(seq(397,412))
new_Obs[RIindex] = log(cnst_RItrans+new_Obs[RIindex])
# some predictors and combined predictors in new dataframe
newpredictors = data.frame(
  WSPD_500_Delta = new_Obs$WSPD_500_max-new_Obs$WSPD_500_min,
  WSPD_500_Delta2 = new_Obs$WSPD_500_q0.98-new_Obs$WSPD_500_q0.02,
  WSPD_850_Delta = new_Obs$WSPD_850_max-new_Obs$WSPD_850_min,
  WSPD_850_Delta2 = new_Obs$WSPD_850_q0.98-new_Obs$WSPD_850_q0.02,
  Moisture_convergence_q0.98_q0.02 = new_Obs$Moisture_convergence_q0.98 - new_Obs$Moisture_convergence_q0.02,
  Moisture_convergence_q0.90_q0.10 = new_Obs$Moisture_convergence_q0.90 - new_Obs$Moisture_convergence_q0.10,
  U_500_Delta = new_Obs$U_500_max - new_Obs$U_500_min,
  V_500_Delta = new_Obs$V_500_max - new_Obs$V_500_min,
  U_700_Delta = new_Obs$U_700_max - new_Obs$U_700_min,
  V_700_Delta = new_Obs$V_700_max - new_Obs$V_700_min,
  U_850_Delta = new_Obs$U_850_max - new_Obs$U_850_min,
  V_850_Delta = new_Obs$V_850_max - new_Obs$V_850_min,
  Cloud_layers_depth_Delta = new_Obs$Cloud_layers_depth_max-new_Obs$Cloud_layers_depth_min,
  Bulk_shear_850_700_Delta = new_Obs$Bulk_shear_850_700_max -  new_Obs$Bulk_shear_850_700_min,
  Bulk_shear_850_500_Delta = new_Obs$Bulk_shear_850_500_max -  new_Obs$Bulk_shear_850_500_min,
  Boyden_PW1 = (new_Obs$Boyden_max-85)*PW1,
  Boyden_PW2 = (new_Obs$Boyden_max-85)*PW2,
  Boyden_hprecipmax = (new_Obs$Boyden_max-85)*hprecipcubertmax,
  Bradbury_PW1 = (new_Obs$Bradbury_min-14)*PW1,
  Bradbury_PW2 = (new_Obs$Bradbury_min-14)*PW2,
  Bradbury_hprecipmax = (new_Obs$Bradbury_min-14)*hprecipcubertmax,
  Edward_PW1 = (new_Obs$ThetaW_925_max-new_Obs$ThetaW_500_min)*PW1,
  Edward_PW2 = (new_Obs$ThetaW_925_max-new_Obs$ThetaW_500_min)*PW2,
  Edward_hprecipmax = (new_Obs$ThetaW_925_max-new_Obs$ThetaW_500_min)*hprecipcubertmax,
  LidS_PW = log(0.01+new_Obs$Lid_Strength_min*PW2),
  LidS_hprecipmax = log(0.01+new_Obs$Lid_Strength_min*hprecipcubertmax),
  Lifted_PW1 = PW1/(new_Obs$LiftedIndex_min+6),
  Lifted_PW2 = PW2/(new_Obs$LiftedIndex_min+6),
  Lifted_hprecipmax = hprecipcubertmax/(new_Obs$LiftedIndex_min+6),
  binary_1 = as.numeric((new_Obs$Boyden_max > 93) *(new_Obs$PrecipWater_max > 10) *(new_Obs$Lid_Strength_min < 3) ),
#  binary_2 = as.numeric((new_Obs$Boyden_max > 93) *(new_Obs$PrecipWater_max > 10) *(new_Obs$Lid_Strength_min < 3)*(new_Obs$ModJefferson_ymax > 20)*(new_Obs$TotalTotals_2e_0mabovegnd_6hrlymax > 20) ),
 # cos_WDIR_500_max = cos(new_Obs$WDIR_500mb_6hrlymax),
#  sin_WDIR_500_max = sin(new_Obs$WDIR_500mb_6hrlymax),
#  cos_WDIR_850_max = cos(new_Obs$WDIR_850mb_6hrlymax),
#  sin_WDIR_850_max = sin(new_Obs$WDIR_850mb_6hrlymax),
#  cos_WDIR_500_min = cos(new_Obs$WDIR_500mb_6hrlymin),
#  sin_WDIR_500_min = sin(new_Obs$WDIR_500mb_6hrlymin),
#  cos_WDIR_850_min = cos(new_Obs$WDIR_850mb_6hrlymin),
#  sin_WDIR_850_min = sin(new_Obs$WDIR_850mb_6hrlymin),
  sfccape_pow_0.2max = new_Obs$SBCAPE_max^0.2,
  mucape_pow_0.2max = new_Obs$MUCAPE_max^0.2,
#  CAPE_CIN_addmax = CAPE_CIN_add_max,
#  CAPE_CIN_addsfcmax = CAPE_CIN_add_sfcmax,
#  CAPE_CIN_addmeanCIN = CAPE_CIN_add_meanCIN,
#  CAPE_CIN_pow_0.5max <- ifelse(CAPE_CIN_add_max<0,0,CAPE_CIN_add_max^0.5),
#  CAPE_CIN_pow_0.5sfc <- ifelse(CAPE_CIN_add_sfcmax<0,0,CAPE_CIN_add_sfcmax^0.5),
#  CAPE_CIN_pow_0.5meanCIN <- ifelse(CAPE_CIN_add_meanCIN<0,0,CAPE_CIN_add_meanCIN^0.5),
  SWEAT_max_pow0.2 = new_Obs$SWEAT_max^0.2,
  abs_helicity_pow0.1_max = abs(new_Obs$Helicity_max)^0.1,
  abs_helicity_pow0.1_min = abs(new_Obs$Helicity_min)^0.1,
  mod_Jeff_pow0.2_max_0replace <- ifelse(new_Obs$ModJefferson_max<0,0,new_Obs$ModJefferson_max^0.2),
  mod_Jeff_pow0.2_min_0replace <- ifelse(new_Obs$ModJefferson_min<0,0,new_Obs$ModJefferson_min^0.2),
  logPWmax = PW1,
  logPWmin = log(new_Obs$PrecipWater_min),
  coast = as.numeric(new_Obs$Region %in% c(seq(1,5),7)),
#  inland = as.numeric(new_Obs$Region %in% c(6,seq(8,12))),
  #ff300.x_pow0.2 = new_Obs$ff300.x^0.2,
  #ff300.y_pow0.2 = new_Obs$ff300.y^0.2,
  #ta5010.x_pow2 = new_Obs$ta5010.x^2,
  #ta5010.y_pow2 = new_Obs$ta5010.y^2,
  stepwise_Inst_PW = stepwise_first,
  stepwise_Inst_PW_rain = stepwise_first + (new_Obs$Rain_acc_max > 5),
  conditional_Inst_PW = conditional_first,
  conditional_Inst_PW_rain = conditional_first * (new_Obs$Rain_acc_max > 5),
  new_Obs[RIindex]
#####
# ADD COAST REGION PARTS
)
graupelindex = seq(53,60)
cloudiceindex = seq(77,84)
cloudwaterindex = seq(85,92)
rainindex = seq(93,100)
snowindex = seq(101,108)
moistureconverindex = seq(469,476)
new_Obs[graupelindex] = new_Obs[graupelindex]^0.20
new_Obs[cloudiceindex] = new_Obs[cloudiceindex]^0.50
new_Obs[cloudwaterindex] = new_Obs[cloudwaterindex]^0.20
new_Obs[rainindex] = new_Obs[rainindex]^(1.00/3)
new_Obs[snowindex] = new_Obs[snowindex]^(1.00/3)
new_Obs[moistureconverindex] = abs(new_Obs[moistureconverindex])^(1./7)
colnames(new_Obs)[c(graupelindex, cloudiceindex,cloudwaterindex,rainindex,snowindex,moistureconverindex)] <- paste0("Trans_",
                                                                                                                    colnames(new_Obs[c(graupelindex, cloudiceindex,cloudwaterindex,rainindex,snowindex,moistureconverindex)]))


new_Obs = cbind(new_Obs, newpredictors)
selected_predictors = c(seq(3,6),12,13,seq(21,length(new_Obs)))
new_Obs = new_Obs[selected_predictors]

# merge with old dataframe
#names(new_Obs)[241:244]=paste0(names(new_Obs)[241:244],"_transformedlog+cnst")

###### REMOVE WRONG VALUES ######
#new_Obs <- filter(new_Obs, LCL_0mabovegnd_6hrlymax < 12500)
#new_Obs = new_Obs[2:length(new_Obs)]
#colnames(new_Obs)[1:10]=c("id","region","Year","Month","Day","validtime","radarmax","validdate","T0","leadtime")

## WRITE NEW DATASET TO FILE ####
write.csv(new_Obs, "full_final00z_dataset.csv")
