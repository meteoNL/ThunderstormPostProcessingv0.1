setwd("/usr/people/groote/ThunderstormPostProcessingv1/rangerres")
data = read.csv("ECimp.csv")

#select certain hyperparameter optionally:
mtries = 2
nsize = 3
#data <- filter(data, mtry == mtries, min_n_size == nsize)
maxpredictors = 10
data = filter(data, npred < (maxpredictors + 0.5))
freqtable=prop.table(table(data$varnames))
plot(sort(as.numeric(freqtable),decreasing = TRUE))
freqtable=sort(freqtable[freqtable > 0],decreasing = TRUE)
write.csv(freqtable, file = "EC_case_LT_3_VT_1218.csv")
