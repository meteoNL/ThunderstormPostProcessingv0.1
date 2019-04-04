LT_val = seq(1,7)
VT_val = seq(1,4)
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
varindex_shell = c(11,seq(15,18),34,42,43,49,50,65,66,74,81,90,91,seq(95,98),106,113,114,115,seq(120,122),seq(128,130),147,148,153,154,156,164,
                   171,173,174,179,180,184,186,187,188,190,192,194,195,196,198,202,203,204,210,211,216,218,seq(219,224),seq(227,230),235,236,238,
                   240,243,246,seq(248,250),251,seq(255,258),259,seq(263,266),274,275,279,281,282,283,284,286,288,290,291,seq(295,298),299,303,312,
                   315,seq(319,322),seq(324,326),337,338,seq(344,346),347,seq(352,354),356,363,371,seq(375,378),403,404,406,408,410,418,419,426,
                   435,441,442,seq(460,462),seq(464,466),467,472,474,475,482,484,492,494,500,514
)
varindex_shell = varindex_shell+4
varindex_shell = c(varindex_shell, seq(519,525),527,529,531,533,535,539,541,549,seq(550,588),seq(591,593),seq(602,608),623)

varindex_shell2=c(varindex_shell[c(5,6,7,10,12,13,14,15,20,21,23,27,30,31,35,36,37,43,48,52,55,58,60,66,70,77,80,82,87,88,92,97,102,103,105,110,111,115,116,117,118,122,123,124,125,127,129,134,135,136,140,141,146,149,151,152,153,155,156,157,158,159,160,163,167,170,171,172,173,174,185,187,189,191,193,194,197,198,199,202,203,205,206,208,209,216)])
varindex_shell = c(varindex_shell2, 416, 359, 526, 534)
for(LT_i in LT_val){
  setwd("/usr/people/groote/ThunderstormPostProcessingv1/R")
  tryCatch({source("ELR_building_blocks4.R")},
           error = function(err){print(err)
             setwd("/usr/people/groote/ThunderstormPostProcessingv1/R")
             source("ELR_building_blocks3.R")
           })
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

print("The tasks have been completed")
