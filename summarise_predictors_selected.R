library(dplyr)
VTs = c("VT_0006","VT_0612","VT_1218","VT_1800")
LTs = seq(1,7)
modnrELR = c(4,8,12)
modnrLR = c(6,12,18)

records = data.frame()
for(VT in VTs){
  for(LT in LTs){
    for(i in seq(1,3)){
      modnrELR_i = modnrELR[i]
      modnrLR_i = modnrLR[i]
      fn = paste0("Models_",VT,"_LT_",LT,"_npred_111")
      setwd("/usr/people/groote/ThunderstormPostProcessingv1/ELRres")
      predictorsELR = data.frame(t(names(readRDS(fn)[[modnrELR_i]]$coefficients$location[1:4])))
      setwd("/usr/people/groote/ThunderstormPostProcessingv1/LRres")
      predictorsLR = data.frame(t(names(readRDS(fn)[[modnrLR_i]][2:7])))
      add_record = data.frame(Validtime=VT, Leadtime_count=LT, i_th_test_year = i, predictors1=predictorsELR, predictors2=predictorsLR)
      records = rbind(records, add_record)
    }
  }
}

LR_freqtable1 <- records %>% group_by(Validtime, predictors2.X1) %>% summarise(Freq = n())
LR_freqtable2 <- records %>% group_by(Validtime, predictors2.X2) %>% summarise(Freq = n())
LR_freqtable3 <- records %>% group_by(Validtime, predictors2.X3) %>% summarise(Freq = n())

ELR_freqtable1 <- records %>% group_by(Validtime, predictors1.X1) %>% summarise(Freq = n())
ELR_freqtable2 <- records %>% group_by(Validtime, predictors1.X2) %>% summarise(Freq = n())

print.data.frame(LR_freqtable1)
print.data.frame(ELR_freqtable1)
