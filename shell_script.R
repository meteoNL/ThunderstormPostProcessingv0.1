LT_val = seq(1,1)
VT_val = seq(1,4)
varindex_shell = c(seq(18,35),seq(39,42),seq(44,71),seq(73,76),seq(81,82),seq(86,101),seq(192,200),seq(203,205),seq(208,217),seq(220,230),seq(237,242))

for(LT_i in LT_val){
  for(VT_i in VT_val){
    setwd("/usr/people/groote/ThunderstormPostProcessingv1/R")
    source("ELR_building_blocks2.R")
    while (length(dev.list())>0){
      dev.off()
    }
    setwd("/usr/people/groote/ThunderstormPostProcessingv1/R")
    source("LR_stepwise_blocks_4.R")
    while (length(dev.list())>0){
      dev.off()
    }
    setwd("/usr/people/groote/ThunderstormPostProcessingv1/R")
    source("rangerthres_newrandom.R")
    setwd("/usr/people/groote/ThunderstormPostProcessingv1/R")
    source("rangertrying_newrandom.R")
  }
}

print("The tasks have been completed")
