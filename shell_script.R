### How to use the scripts??

# first (1), the read_grib_final.R-script should be run to make a dataset with model output derived predictors.
# then (2), the lightning data should be connected to these with the processing_lightning_data.R-script.
# then (3), the new predictors such as BradburyPW's, EdwardPW's and BoydenPW's (and some more) need to be calculated,
#          as well as some applications of transformations with the new_additional_predictors_script.R-script.

# after the three steps above, the dataset to generate models with is ready. This script can be run to apply the four methods.
#          in these four scripts, the variables "VT" (file names of models/verification), "writedir" (writing directory) and
#          "ObsPV" (source data) need to be changed into custom values






# define valid times and lead times
LT_val = seq(1,7)
VT_val = seq(1,4) #in fact NOT used

#define potential predictor set
prednames40 <- c("PrecipWater_max","MSLP_max","Trans_Graupel_col_max","ThetaW_850_max","Bradbury_min","MUCAPE_max","MUCIN_mean","LNB_max","Helicity_min",
                 "Helicity_max","Trans_Moisture_convergence_min","Trans_Moisture_convergence_max","Bulk_shear_850_500_mean","dp_dt","coast","Boyden_max",
                 "ModJefferson_max","K_index_max","Fateev_max","Edward_PW1","Edward_PW2","Trans_Snow_max","DPT_700_max","DPT_850_max","ThetaW_925_max",
                 "Showalter_min","Jefferson_max","LiftedIndex_min","SBCAPE_max","MUCAPE_q0.90","LNB_q0.90","Bradbury_base925_max" ,"MUCAPE_graupel_pows_max",
                 "MUCAPE_snow_pows_max","ifelse..dataset.SBCAPE...dataset.SBCIN....0..0...dataset.SBCAPE..._max","ifelse..dataset.MUCAPE...dataset.MUCIN....0..0...dataset.MUCAPE..._max",
                 "Boyden_PW1","Boyden_PW2","Bradbury_PW1","Bradbury_PW2")
varindexshell = prednames40

#apply all the fitting procedures, or the ones that are not commented
# the dataframe as input for predictors is defined in the four source scripts
# the export names of verification score files annd model fits are also defined in the files. The unique characteristic name of these file names is the variable VT.
for(LT_i in LT_val){
  setwd("/usr/people/groote/ThunderstormPostProcessingv1/R")
  source("LR_stepwise_blocks_4.R") # does the LR fits
  while (length(dev.list())>0){
    dev.off()
  }
  tryCatch({source("ELR_building_blocks4b.R")},
         error = function(err){print(err)
           setwd("/usr/people/groote/ThunderstormPostProcessingv1/R")
            source("ELR_building_blocks3b.R")
        }) # does the ELR fit with two different random seeds (the secodn seed setting if the first is not working)
  while (length(dev.list())>0){
  dev.off()
  }
  setwd("/usr/people/groote/ThunderstormPostProcessingv1/R")
  setwd("/usr/people/groote/ThunderstormPostProcessingv1/R")
  source("rangerthres_finalcrossc.R") #does the threshold predictions for 25-400 dis./5min. with QRF
  setwd("/usr/people/groote/ThunderstormPostProcessingv1/R")
  source("rangertrying_finalcross_reldi.R") #does the thunderstorm occurrence predictions with QRF
}

print("The tasks have been completed")















## Below: find definitions of other predictor sets, not hard coded


#varindex_shell = c(seq(18,35),seq(39,42),seq(44,71),seq(73,76),seq(81,82),seq(86,101),seq(192,200),seq(203,205),seq(208,217),seq(220,230),seq(237,239),c(241,242))
#varindex_shell = c(varindex_shell,85)
#varindex_shell = c(
#  seq(28,35),44,51,seq(52,59),60,seq(65,67),seq(76,83),84,85,seq(88,91),92,seq(96,99),100,seq(104,107),108,seq(112,115),116,117,seq(121,123),124,125,
#  130,131,seq(132,135),seq(137,139),seq(140,143),seq(145,147),164,165,170,171,seq(172,174),178,179,180,181,187,seq(188,191),seq(193,195),196,197,201,
#  203,204,205,207,209,211,212,213,215,219,220,221,227,228,229,231,233,235,seq(236,243),seq(244,247),seq(249,251),seq(252,255),257,259,
#  seq(260,263),seq(265,267),seq(268,275),seq(276,283),seq(284,291),seq(292,299),300,301,303,305,307,seq(308,315),seq(316,323),seq(324,331),
#  seq(332,339),seq(340,343),345,347,348,seq(352,355),356,seq(360,363),364,seq(369,371),seq(372,379),380,381,seq(384,387),388,389,seq(391,395),
#  420,421,423,425,427,429,430,432,434,435,436,443,452,453,455,458,459,seq(477,479),seq(481,483),484,485,
#  487,489,491,492,499,500,501,507,508,509,511,515,517,521,523,524,525,529,531,533,535,537,539,541,545,547,548,550,554,555,seq(556,558),seq(560,562),seq(564,566),seq(568,570),572,573,
#  seq(577,600),seq(603,605),seq(607,617),seq(620,628),seq(633,636),seq(641,643)
#)
#varindex_shell = varindex_shell-17

## start of 228 set
varindex_shell = c(11,seq(15,18),34,42,43,49,50,65,66,74,81,90,91,seq(95,98),106,113,114,115,seq(120,122),seq(128,130),147,148,153,154,156,164,
                   171,173,174,179,180,184,186,187,188,190,192,194,195,196,198,202,203,204,210,211,216,218,seq(219,224),seq(227,230),235,236,238,
                   240,243,246,seq(248,250),251,seq(255,258),259,seq(263,266),274,275,279,281,282,283,284,286,288,290,291,seq(295,298),299,303,312,
                   315,seq(319,322),seq(324,326),337,338,seq(344,346),347,seq(352,354),356,363,371,seq(375,378),403,404,406,408,410,418,419,426,
                   435,441,442,seq(460,462),seq(464,466),467,472,474,475,482,484,492,494,500,514
)
varindex_shell = varindex_shell+4
varindex_shell = c(varindex_shell, seq(519,525),527,529,531,533,535,539,541,549,seq(550,588),seq(591,593),seq(595,601),seq(602,608),623)
## end of 228 set

## start of 90 set
varindex_shell2=c(varindex_shell[c(5,6,7,10,12,13,14,15,20,21,23,27,30,31,35,36,37,43,48,52,55,58,60,66,70,77,80,82,87,88,92,97,102,103,105,110,111,115,116,117,118,122,123,124,125,127,129,134,135,136,140,141,146,149,151,152,153,155,156,157,158,159,160,163,167,170,171,172,173,174,185,187,189,191,193,194,197,198,199,202,203,205,206,208,209,216)])
varindex_shell = c(varindex_shell2, 416, 359, 526, 534, 604)

## end of 90 set

#15 predictors set for 00z run
varindex_shell = c(22,38,54,198,224,350,359,382,416,422,464,470,526,550,604)
#additional 6 for 21 set for 00z run
varindex_shell = c(varindex_shell,222,262,270,278,586,587)
#additional 20 for 40 set for 00z run
varindex_shell = c(varindex_shell, 102,126,134,206,232,254,328,342,348,380,478,519,520,522,525,580,581,583,584)
#varindex_shell = c(22,38,54,198,224,350,359,382,416,422,550)
#15 predictors set for 06z run
varindex_shell = c(22,38,54,206,232,358,367,390,424,430,472,478,534,558,611)
#15 predictors set for 06z run with vertical velocity
varindex_shell = c(22,38,54,103,104,110,206,232,358,367,390,424,430,472,478,534,558,611)
varindex_shell = c(22,38,54,350)
