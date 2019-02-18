# read data from file
setwd("/usr/people/groote/ThunderstormPostProcessingv1/rangerres")
data = read.csv("ECimp.csv")

# select certain hyperparameter optionally (and if so, switch potentially commented filter on):
mtries = 2
nsize = 3
# data <- filter(data, mtry == mtries, min_n_size == nsize)
maxpredictors = 10
data = filter(data, npred < (maxpredictors + 0.5))

# generate table with relative frequencies of the data table; plot it
freqtable=prop.table(table(data$varnames))
plot(sort(as.numeric(freqtable),decreasing = TRUE))

#sort and write relevant part of the table with predictors that occur
freqtable=sort(freqtable[freqtable > 0],decreasing = TRUE)
write.csv(freqtable, file = "EC_case_LT_3_VT_1218.csv")
