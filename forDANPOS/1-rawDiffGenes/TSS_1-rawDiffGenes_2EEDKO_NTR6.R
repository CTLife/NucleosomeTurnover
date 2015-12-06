#############################################################################################################################
## Part  6:  Classify DNA regions based on NTR.
#############################################################################################################################





#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################


subdir_1_part6 <- paste(Part6_g,  "/1-one1-diff-NTR", sep = "")
if( ! file.exists(subdir_1_part6) ) { dir.create(subdir_1_part6) }

length(row_Average_one1_EEDheto_NTR2)       
length(row_Average_one1_EEDko_NTR2)         
row_Average_one1_diffNTR <- (row_Average_one1_EEDheto_NTR2 - row_Average_one1_EEDko_NTR2 )
length(row_Average_one1_diffNTR)
summary(row_Average_one1_diffNTR)

one1_NTR_threshold <- 0.05
length( one1_week0_dist[one1_week0_dist < -one1_NTR_threshold] )  ## NTR of one1_week0_EEDheto is less than one1_week0_EEDko
length( one1_week0_dist[one1_week0_dist >  one1_NTR_threshold] )  ## NTR of one1_week0_EEDheto is more than one1_week0_EEDko



## 差异不显著
one1_NTR_diffIndex_1 <- ( abs(row_Average_one1_diffNTR) <= one1_NTR_threshold  )
length(one1_NTR_diffIndex_1[one1_NTR_diffIndex_1])
write.table(x=one1_NTR_diffIndex_1,   file = paste(subdir_1_part6, "/1-one1-NoDiff.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")


## 差异显著  one1_week0_EEDheto is less than one1_week0_EEDko
one1_NTR_diffIndex_2 <- ( row_Average_one1_diffNTR < -one1_NTR_threshold )
length(one1_NTR_diffIndex_2[one1_NTR_diffIndex_2])
write.table(x=one1_NTR_diffIndex_2,   file = paste(subdir_1_part6, "/1-one1-EEDheto-Less-Than-EEDko.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")


## 差异显著  one1_week0_EEDheto is more than one1_week0_EEDko
one1_NTR_diffIndex_3 <- ( row_Average_one1_diffNTR >  one1_NTR_threshold )
length(one1_NTR_diffIndex_3[one1_NTR_diffIndex_3])
write.table(x=one1_NTR_diffIndex_3,   file = paste(subdir_1_part6, "/1-one1-EEDheto-More-Than-EEDko.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")






MyAverageLines_3(vector2=c(colMeans(Average_one1_week0_EEDheto[one1_NTR_diffIndex_1, ]),  colMeans(Average_one1_week0_EEDko[one1_NTR_diffIndex_1, ]), 
                           colMeans(Average_one1_week4_EEDheto[one1_NTR_diffIndex_1, ]),  colMeans(Average_one1_week4_EEDko[one1_NTR_diffIndex_1, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part6,     fileName2="A-one1-NoDiffNTR",  
                 title2="Down-regulated Genes (No DiffNTR)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=-0.3,   Ymax2=1,    height2=3.2,   width2=5.55 , center2="TSS" )


one1_lenNTR_A <- nrow(Average_one1_week0_EEDheto[one1_NTR_diffIndex_1, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_one1_week0_EEDheto[one1_NTR_diffIndex_1, ]),  rowMeans(Average_one1_week0_EEDko[one1_NTR_diffIndex_1, ]), 
                            rowMeans(Average_one1_week4_EEDheto[one1_NTR_diffIndex_1, ]),  rowMeans(Average_one1_week4_EEDko[one1_NTR_diffIndex_1, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", one1_lenNTR_A),   rep("week0_EEDko",  one1_lenNTR_A), 
                                 rep("week4_EEDheto", one1_lenNTR_A),   rep("week4_EEDko",  one1_lenNTR_A) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_1_part6,  fileName2="A-one1-all-NoDiffNTR-BoxViolin",  
                  title2="Down-regulated Genes  (No DiffNTR)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_one1_week0_EEDheto[one1_NTR_diffIndex_1, 201:400]),  rowMeans(Average_one1_week0_EEDko[one1_NTR_diffIndex_1, 201:400]), 
                            rowMeans(Average_one1_week4_EEDheto[one1_NTR_diffIndex_1, 201:400]),  rowMeans(Average_one1_week4_EEDko[one1_NTR_diffIndex_1, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", one1_lenNTR_A),   rep("week0_EEDko",  one1_lenNTR_A), 
                                 rep("week4_EEDheto", one1_lenNTR_A),   rep("week4_EEDko",  one1_lenNTR_A) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_1_part6,  fileName2="A-one1-2kbRegion-NoDiffNTR-BoxViolin",  
                  title2="Down-regulated Genes  (No DiffNTR)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )










MyAverageLines_3(vector2=c(colMeans(Average_one1_week0_EEDheto[one1_NTR_diffIndex_2, ]),  colMeans(Average_one1_week0_EEDko[one1_NTR_diffIndex_2, ]), 
                           colMeans(Average_one1_week4_EEDheto[one1_NTR_diffIndex_2, ]),  colMeans(Average_one1_week4_EEDko[one1_NTR_diffIndex_2, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part6,     fileName2="B-one1-EEDhetoLess",  
                 title2="Down-regulated Genes (EEDheto Less)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=-0.3,   Ymax2=1,     height2=3.2,   width2=5.55 , center2="TSS" )


one1_lenNTR_B <- nrow(Average_one1_week0_EEDheto[one1_NTR_diffIndex_2, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_one1_week0_EEDheto[one1_NTR_diffIndex_2, ]),  rowMeans(Average_one1_week0_EEDko[one1_NTR_diffIndex_2, ]), 
                            rowMeans(Average_one1_week4_EEDheto[one1_NTR_diffIndex_2, ]),  rowMeans(Average_one1_week4_EEDko[one1_NTR_diffIndex_2, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", one1_lenNTR_B),   rep("week0_EEDko",  one1_lenNTR_B), 
                                 rep("week4_EEDheto", one1_lenNTR_B),   rep("week4_EEDko",  one1_lenNTR_B) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_1_part6,  fileName2="B-one1-all-EEDhetoLess-BoxViolin",  
                  title2="Down-regulated Genes  (EEDheto Less)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_one1_week0_EEDheto[one1_NTR_diffIndex_2, 201:400]),  rowMeans(Average_one1_week0_EEDko[one1_NTR_diffIndex_2, 201:400]), 
                            rowMeans(Average_one1_week4_EEDheto[one1_NTR_diffIndex_2, 201:400]),  rowMeans(Average_one1_week4_EEDko[one1_NTR_diffIndex_2, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", one1_lenNTR_B),   rep("week0_EEDko",  one1_lenNTR_B), 
                                 rep("week4_EEDheto", one1_lenNTR_B),   rep("week4_EEDko",  one1_lenNTR_B) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_1_part6,  fileName2="B-one1-2kbRegion-EEDhetoLess-BoxViolin",  
                  title2="Down-regulated Genes  (EEDheto Less)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )












MyAverageLines_3(vector2=c(colMeans(Average_one1_week0_EEDheto[one1_NTR_diffIndex_3, ]),  colMeans(Average_one1_week0_EEDko[one1_NTR_diffIndex_3, ]), 
                           colMeans(Average_one1_week4_EEDheto[one1_NTR_diffIndex_3, ]),  colMeans(Average_one1_week4_EEDko[one1_NTR_diffIndex_3, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part6,     fileName2="C-one1-EEDhetoMore",  
                 title2="Down-regulated Genes (EEDheto More)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=0.1,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


one1_lenNTR_C <- nrow(Average_one1_week0_EEDheto[one1_NTR_diffIndex_3, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_one1_week0_EEDheto[one1_NTR_diffIndex_3, ]),  rowMeans(Average_one1_week0_EEDko[one1_NTR_diffIndex_3, ]), 
                            rowMeans(Average_one1_week4_EEDheto[one1_NTR_diffIndex_3, ]),  rowMeans(Average_one1_week4_EEDko[one1_NTR_diffIndex_3, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", one1_lenNTR_C),   rep("week0_EEDko",  one1_lenNTR_C), 
                                 rep("week4_EEDheto", one1_lenNTR_C),   rep("week4_EEDko",  one1_lenNTR_C) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_1_part6,  fileName2="C-one1-all-EEDhetoMore-BoxViolin",  
                  title2="Down-regulated Genes  (EEDheto More)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_one1_week0_EEDheto[one1_NTR_diffIndex_3, 201:400]),  rowMeans(Average_one1_week0_EEDko[one1_NTR_diffIndex_3, 201:400]), 
                            rowMeans(Average_one1_week4_EEDheto[one1_NTR_diffIndex_3, 201:400]),  rowMeans(Average_one1_week4_EEDko[one1_NTR_diffIndex_3, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", one1_lenNTR_C),   rep("week0_EEDko",  one1_lenNTR_C), 
                                 rep("week4_EEDheto", one1_lenNTR_C),   rep("week4_EEDko",  one1_lenNTR_C) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_1_part6,  fileName2="C-one1-2kbRegion-EEDhetoMore-BoxViolin",  
                  title2="Down-regulated Genes  (EEDheto More)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )













column_NTR2_one1_EEDheto_same   <-  log( colMeans(Average_one1_week0_EEDheto[one1_NTR_diffIndex_1, ]) / colMeans(Average_one1_week4_EEDheto[one1_NTR_diffIndex_1, ]) )/4
column_NTR2_one1_EEDko_same     <-  log( colMeans(Average_one1_week0_EEDko[one1_NTR_diffIndex_1, ])   / colMeans(Average_one1_week4_EEDko[one1_NTR_diffIndex_1, ])   )/4
column_NTR2_one1_EEDheto_less   <-  log( colMeans(Average_one1_week0_EEDheto[one1_NTR_diffIndex_2, ]) / colMeans(Average_one1_week4_EEDheto[one1_NTR_diffIndex_2, ]) )/4
column_NTR2_one1_EEDko_less     <-  log( colMeans(Average_one1_week0_EEDko[one1_NTR_diffIndex_2, ])   / colMeans(Average_one1_week4_EEDko[one1_NTR_diffIndex_2, ])   )/4
column_NTR2_one1_EEDheto_more   <-  log( colMeans(Average_one1_week0_EEDheto[one1_NTR_diffIndex_3, ]) / colMeans(Average_one1_week4_EEDheto[one1_NTR_diffIndex_3, ]) )/4
column_NTR2_one1_EEDko_more     <-  log( colMeans(Average_one1_week0_EEDko[one1_NTR_diffIndex_3, ])   / colMeans(Average_one1_week4_EEDko[one1_NTR_diffIndex_3, ])   )/4

row_NTR2_one1_EEDheto_same   <-  log( row_Average_one1_week0_EEDheto[one1_NTR_diffIndex_1] / row_Average_one1_week4_EEDheto[one1_NTR_diffIndex_1] )/4
row_NTR2_one1_EEDko_same     <-  log( row_Average_one1_week0_EEDko[one1_NTR_diffIndex_1]   / row_Average_one1_week4_EEDko[one1_NTR_diffIndex_1]   )/4
row_NTR2_one1_EEDheto_less   <-  log( row_Average_one1_week0_EEDheto[one1_NTR_diffIndex_2] / row_Average_one1_week4_EEDheto[one1_NTR_diffIndex_2] )/4
row_NTR2_one1_EEDko_less     <-  log( row_Average_one1_week0_EEDko[one1_NTR_diffIndex_2]   / row_Average_one1_week4_EEDko[one1_NTR_diffIndex_2]   )/4
row_NTR2_one1_EEDheto_more   <-  log( row_Average_one1_week0_EEDheto[one1_NTR_diffIndex_3] / row_Average_one1_week4_EEDheto[one1_NTR_diffIndex_3] )/4
row_NTR2_one1_EEDko_more     <-  log( row_Average_one1_week0_EEDko[one1_NTR_diffIndex_3]   / row_Average_one1_week4_EEDko[one1_NTR_diffIndex_3]   )/4


MyAverageLines_3(vector2=c(column_NTR2_one1_EEDheto_same,  column_NTR2_one1_EEDko_same,  column_NTR2_one1_EEDheto_less, column_NTR2_one1_EEDko_less,
                           column_NTR2_one1_EEDheto_more, column_NTR2_one1_EEDko_more ),   
                 numSample2=6,   
                 sampleType2=c( rep("EEDheto_same", numOfColumns1),   rep("EEDko_same", numOfColumns1),  rep("EEDheto_less", numOfColumns1),   rep("EEDko_less", numOfColumns1), 
                                rep("EEDheto_more", numOfColumns1),   rep("EEDko_more", numOfColumns1)  ), 
                 sampleRank2=c( "EEDheto_same",   "EEDko_same",  "EEDheto_less",  "EEDko_less",  "EEDheto_more",  "EEDko_more"  ),  
                 colours2=c("EEDheto_same"="red",   "EEDko_same"="red4",  "EEDheto_less"="blue",  "EEDko_less"="blue4",  "EEDheto_more"="green",  "EEDko_more"="green4"),   
                 path2=subdir_1_part6,     fileName2="1-one1-NTR",  
                 title2="Down-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR2_one1_EEDheto_same,  column_NTR2_one1_EEDko_same ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_same", numOfColumns1),   rep("EEDko_same", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_same",   "EEDko_same"  ),  
                 colours2=c("EEDheto_same"="red",   "EEDko_same"="red4" ),   
                 path2=subdir_1_part6,     fileName2="1-one1-NTR-same",  
                 title2="Down-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR2_one1_EEDheto_less,  column_NTR2_one1_EEDko_less ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_less", numOfColumns1),   rep("EEDko_less", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_less",   "EEDko_less"  ),  
                 colours2=c("EEDheto_less"="blue",   "EEDko_less"="blue4" ),   
                 path2=subdir_1_part6,     fileName2="1-one1-NTR-less",  
                 title2="Down-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR2_one1_EEDheto_more,  column_NTR2_one1_EEDko_more ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_more", numOfColumns1),   rep("EEDko_more", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_more",   "EEDko_more"  ),  
                 colours2=c("EEDheto_more"="blue",   "EEDko_more"="blue4" ),   
                 path2=subdir_1_part6,     fileName2="1-one1-NTR-more",  
                 title2="Down-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyBoxViolinPlot_1(vector2=c(row_NTR2_one1_EEDheto_same,  row_NTR2_one1_EEDko_same,  row_NTR2_one1_EEDheto_less, row_NTR2_one1_EEDko_less,
                            row_NTR2_one1_EEDheto_more,  row_NTR2_one1_EEDko_more  ),  
                  sampleType2=c( rep("EEDheto_same", one1_lenNTR_A),   rep("EEDko_same", one1_lenNTR_A),  rep("EEDheto_less", one1_lenNTR_B),   rep("EEDko_less", one1_lenNTR_B), 
                                 rep("EEDheto_more", one1_lenNTR_C),   rep("EEDko_more", one1_lenNTR_C)  ), 
                  sampleRank2=c( "EEDheto_same",   "EEDko_same",  "EEDheto_less",  "EEDko_less",  "EEDheto_more",  "EEDko_more"  ),  
                  colours2=c("red",  "red4",   "blue",    "blue4", "green", "green4" ),   
                  path2=subdir_1_part6,  fileName2="1-one1-NTR-boxViolin",  
                  title2="Down-regulated Genes",  xLab2="Samples",  yLab2="NTR",   
                  height2=3.88,   width2=4, Ymin2=-0.3, Ymax2=1 )





















########################################################################
subdir_2_part6 <- paste(Part6_g,  "/2-two2-diff-NTR", sep = "")
if( ! file.exists(subdir_2_part6) ) { dir.create(subdir_2_part6) }

length(row_Average_two2_EEDheto_NTR2)       
length(row_Average_two2_EEDko_NTR2)         
row_Average_two2_diffNTR <- (row_Average_two2_EEDheto_NTR2 - row_Average_two2_EEDko_NTR2 )
length(row_Average_two2_diffNTR)
summary(row_Average_two2_diffNTR)

row_Average_two2_diffNTR[is.na(row_Average_two2_diffNTR)] <- 0
row_Average_two2_diffNTR[row_Average_two2_diffNTR>  1]    <- 1
row_Average_two2_diffNTR[row_Average_two2_diffNTR< -1]    <- -1
summary(row_Average_two2_diffNTR)

two2_NTR_threshold <- 0.05
length( two2_week0_dist[two2_week0_dist < -two2_NTR_threshold] )  ## NTR of two2_week0_EEDheto is less than two2_week0_EEDko
length( two2_week0_dist[two2_week0_dist >  two2_NTR_threshold] )  ## NTR of two2_week0_EEDheto is more than two2_week0_EEDko



## 差异不显著
two2_NTR_diffIndex_1 <- ( abs(row_Average_two2_diffNTR) <= two2_NTR_threshold  )
length(two2_NTR_diffIndex_1[two2_NTR_diffIndex_1])
write.table(x=two2_NTR_diffIndex_1,   file = paste(subdir_2_part6, "/1-two2-NoDiff.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")


## 差异显著  two2_week0_EEDheto is less than two2_week0_EEDko
two2_NTR_diffIndex_2 <- ( row_Average_two2_diffNTR < -two2_NTR_threshold )
length(two2_NTR_diffIndex_2[two2_NTR_diffIndex_2])
write.table(x=two2_NTR_diffIndex_2,   file = paste(subdir_2_part6, "/1-two2-EEDheto-Less-Than-EEDko.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")


## 差异显著  two2_week0_EEDheto is more than two2_week0_EEDko
two2_NTR_diffIndex_3 <- ( row_Average_two2_diffNTR >  two2_NTR_threshold )
length(two2_NTR_diffIndex_3[two2_NTR_diffIndex_3])
write.table(x=two2_NTR_diffIndex_3,   file = paste(subdir_2_part6, "/1-two2-EEDheto-More-Than-EEDko.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")




dim(Average_two2_week0_EEDheto[two2_NTR_diffIndex_1, ])

MyAverageLines_3(vector2=c(colMeans(Average_two2_week0_EEDheto[two2_NTR_diffIndex_1, ], na.rm = TRUE),  colMeans(Average_two2_week0_EEDko[two2_NTR_diffIndex_1, ], na.rm = TRUE), 
                           colMeans(Average_two2_week4_EEDheto[two2_NTR_diffIndex_1, ], na.rm = TRUE),  colMeans(Average_two2_week4_EEDko[two2_NTR_diffIndex_1, ], na.rm = TRUE) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_2_part6,     fileName2="A-two2-NoDiffNTR",  
                 title2="Unchanged Genes (No DiffNTR)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=0.1,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


two2_lenNTR_A <- nrow(Average_two2_week0_EEDheto[two2_NTR_diffIndex_1, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_two2_week0_EEDheto[two2_NTR_diffIndex_1, ]),  rowMeans(Average_two2_week0_EEDko[two2_NTR_diffIndex_1, ]), 
                            rowMeans(Average_two2_week4_EEDheto[two2_NTR_diffIndex_1, ]),  rowMeans(Average_two2_week4_EEDko[two2_NTR_diffIndex_1, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", two2_lenNTR_A),   rep("week0_EEDko",  two2_lenNTR_A), 
                                 rep("week4_EEDheto", two2_lenNTR_A),   rep("week4_EEDko",  two2_lenNTR_A) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_2_part6,  fileName2="A-two2-all-NoDiffNTR-BoxViolin",  
                  title2="Unchanged Genes  (No DiffNTR)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_two2_week0_EEDheto[two2_NTR_diffIndex_1, 201:400]),  rowMeans(Average_two2_week0_EEDko[two2_NTR_diffIndex_1, 201:400]), 
                            rowMeans(Average_two2_week4_EEDheto[two2_NTR_diffIndex_1, 201:400]),  rowMeans(Average_two2_week4_EEDko[two2_NTR_diffIndex_1, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", two2_lenNTR_A),   rep("week0_EEDko",  two2_lenNTR_A), 
                                 rep("week4_EEDheto", two2_lenNTR_A),   rep("week4_EEDko",  two2_lenNTR_A) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_2_part6,  fileName2="A-two2-2kbRegion-NoDiffNTR-BoxViolin",  
                  title2="Unchanged Genes  (No DiffNTR)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )










MyAverageLines_3(vector2=c(colMeans(Average_two2_week0_EEDheto[two2_NTR_diffIndex_2, ]),  colMeans(Average_two2_week0_EEDko[two2_NTR_diffIndex_2, ]), 
                           colMeans(Average_two2_week4_EEDheto[two2_NTR_diffIndex_2, ]),  colMeans(Average_two2_week4_EEDko[two2_NTR_diffIndex_2, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_2_part6,     fileName2="B-two2-EEDhetoLess",  
                 title2="Unchanged Genes (EEDheto Less)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NTR_min,   Ymax2=col_NTR_max,    height2=3.2,   width2=5.55 , center2="TSS" )


two2_lenNTR_B <- nrow(Average_two2_week0_EEDheto[two2_NTR_diffIndex_2, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_two2_week0_EEDheto[two2_NTR_diffIndex_2, ]),  rowMeans(Average_two2_week0_EEDko[two2_NTR_diffIndex_2, ]), 
                            rowMeans(Average_two2_week4_EEDheto[two2_NTR_diffIndex_2, ]),  rowMeans(Average_two2_week4_EEDko[two2_NTR_diffIndex_2, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", two2_lenNTR_B),   rep("week0_EEDko",  two2_lenNTR_B), 
                                 rep("week4_EEDheto", two2_lenNTR_B),   rep("week4_EEDko",  two2_lenNTR_B) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_2_part6,  fileName2="B-two2-all-EEDhetoLess-BoxViolin",  
                  title2="Unchanged Genes  (EEDheto Less)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_two2_week0_EEDheto[two2_NTR_diffIndex_2, 201:400]),  rowMeans(Average_two2_week0_EEDko[two2_NTR_diffIndex_2, 201:400]), 
                            rowMeans(Average_two2_week4_EEDheto[two2_NTR_diffIndex_2, 201:400]),  rowMeans(Average_two2_week4_EEDko[two2_NTR_diffIndex_2, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", two2_lenNTR_B),   rep("week0_EEDko",  two2_lenNTR_B), 
                                 rep("week4_EEDheto", two2_lenNTR_B),   rep("week4_EEDko",  two2_lenNTR_B) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_2_part6,  fileName2="B-two2-2kbRegion-EEDhetoLess-BoxViolin",  
                  title2="Unchanged Genes  (EEDheto Less)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )












MyAverageLines_3(vector2=c(colMeans(Average_two2_week0_EEDheto[two2_NTR_diffIndex_3, ]),  colMeans(Average_two2_week0_EEDko[two2_NTR_diffIndex_3, ]), 
                           colMeans(Average_two2_week4_EEDheto[two2_NTR_diffIndex_3, ]),  colMeans(Average_two2_week4_EEDko[two2_NTR_diffIndex_3, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_2_part6,     fileName2="C-two2-EEDhetoMore",  
                 title2="Unchanged Genes (EEDheto More)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NTR_min,   Ymax2=col_NTR_max,    height2=3.2,   width2=5.55 , center2="TSS" )


two2_lenNTR_C <- nrow(Average_two2_week0_EEDheto[two2_NTR_diffIndex_3, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_two2_week0_EEDheto[two2_NTR_diffIndex_3, ]),  rowMeans(Average_two2_week0_EEDko[two2_NTR_diffIndex_3, ]), 
                            rowMeans(Average_two2_week4_EEDheto[two2_NTR_diffIndex_3, ]),  rowMeans(Average_two2_week4_EEDko[two2_NTR_diffIndex_3, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", two2_lenNTR_C),   rep("week0_EEDko",  two2_lenNTR_C), 
                                 rep("week4_EEDheto", two2_lenNTR_C),   rep("week4_EEDko",  two2_lenNTR_C) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_2_part6,  fileName2="C-two2-all-EEDhetoMore-BoxViolin",  
                  title2="Unchanged Genes  (EEDheto More)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_two2_week0_EEDheto[two2_NTR_diffIndex_3, 201:400]),  rowMeans(Average_two2_week0_EEDko[two2_NTR_diffIndex_3, 201:400]), 
                            rowMeans(Average_two2_week4_EEDheto[two2_NTR_diffIndex_3, 201:400]),  rowMeans(Average_two2_week4_EEDko[two2_NTR_diffIndex_3, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", two2_lenNTR_C),   rep("week0_EEDko",  two2_lenNTR_C), 
                                 rep("week4_EEDheto", two2_lenNTR_C),   rep("week4_EEDko",  two2_lenNTR_C) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_2_part6,  fileName2="C-two2-2kbRegion-EEDhetoMore-BoxViolin",  
                  title2="Unchanged Genes  (EEDheto More)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )











column_NTR2_two2_EEDheto_same   <-  log( colMeans(Average_two2_week0_EEDheto[two2_NTR_diffIndex_1, ]) / colMeans(Average_two2_week4_EEDheto[two2_NTR_diffIndex_1, ]) )/4
column_NTR2_two2_EEDko_same     <-  log( colMeans(Average_two2_week0_EEDko[two2_NTR_diffIndex_1, ])   / colMeans(Average_two2_week4_EEDko[two2_NTR_diffIndex_1, ])   )/4
column_NTR2_two2_EEDheto_less   <-  log( colMeans(Average_two2_week0_EEDheto[two2_NTR_diffIndex_2, ]) / colMeans(Average_two2_week4_EEDheto[two2_NTR_diffIndex_2, ]) )/4
column_NTR2_two2_EEDko_less     <-  log( colMeans(Average_two2_week0_EEDko[two2_NTR_diffIndex_2, ])   / colMeans(Average_two2_week4_EEDko[two2_NTR_diffIndex_2, ])   )/4
column_NTR2_two2_EEDheto_more   <-  log( colMeans(Average_two2_week0_EEDheto[two2_NTR_diffIndex_3, ]) / colMeans(Average_two2_week4_EEDheto[two2_NTR_diffIndex_3, ]) )/4
column_NTR2_two2_EEDko_more     <-  log( colMeans(Average_two2_week0_EEDko[two2_NTR_diffIndex_3, ])   / colMeans(Average_two2_week4_EEDko[two2_NTR_diffIndex_3, ])   )/4

row_NTR2_two2_EEDheto_same   <-  log( row_Average_two2_week0_EEDheto[two2_NTR_diffIndex_1] / row_Average_two2_week4_EEDheto[two2_NTR_diffIndex_1] )/4
row_NTR2_two2_EEDko_same     <-  log( row_Average_two2_week0_EEDko[two2_NTR_diffIndex_1]   / row_Average_two2_week4_EEDko[two2_NTR_diffIndex_1]   )/4
row_NTR2_two2_EEDheto_less   <-  log( row_Average_two2_week0_EEDheto[two2_NTR_diffIndex_2] / row_Average_two2_week4_EEDheto[two2_NTR_diffIndex_2] )/4
row_NTR2_two2_EEDko_less     <-  log( row_Average_two2_week0_EEDko[two2_NTR_diffIndex_2]   / row_Average_two2_week4_EEDko[two2_NTR_diffIndex_2]   )/4
row_NTR2_two2_EEDheto_more   <-  log( row_Average_two2_week0_EEDheto[two2_NTR_diffIndex_3] / row_Average_two2_week4_EEDheto[two2_NTR_diffIndex_3] )/4
row_NTR2_two2_EEDko_more     <-  log( row_Average_two2_week0_EEDko[two2_NTR_diffIndex_3]   / row_Average_two2_week4_EEDko[two2_NTR_diffIndex_3]   )/4


MyAverageLines_3(vector2=c(column_NTR2_two2_EEDheto_same,  column_NTR2_two2_EEDko_same,  column_NTR2_two2_EEDheto_less, column_NTR2_two2_EEDko_less,
                           column_NTR2_two2_EEDheto_more, column_NTR2_two2_EEDko_more ),   
                 numSample2=6,   
                 sampleType2=c( rep("EEDheto_same", numOfColumns1),   rep("EEDko_same", numOfColumns1),  rep("EEDheto_less", numOfColumns1),   rep("EEDko_less", numOfColumns1), 
                                rep("EEDheto_more", numOfColumns1),   rep("EEDko_more", numOfColumns1)  ), 
                 sampleRank2=c( "EEDheto_same",   "EEDko_same",  "EEDheto_less",  "EEDko_less",  "EEDheto_more",  "EEDko_more"  ),  
                 colours2=c("EEDheto_same"="red",   "EEDko_same"="red4",  "EEDheto_less"="blue",  "EEDko_less"="blue4",  "EEDheto_more"="green",  "EEDko_more"="green4"),   
                 path2=subdir_2_part6,     fileName2="1-two2-NTR",  
                 title2="Unchanged Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR2_two2_EEDheto_same,  column_NTR2_two2_EEDko_same ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_same", numOfColumns1),   rep("EEDko_same", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_same",   "EEDko_same"  ),  
                 colours2=c("EEDheto_same"="red",   "EEDko_same"="red4" ),   
                 path2=subdir_2_part6,     fileName2="1-two2-NTR-same",  
                 title2="Unchanged Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR2_two2_EEDheto_less,  column_NTR2_two2_EEDko_less ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_less", numOfColumns1),   rep("EEDko_less", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_less",   "EEDko_less"  ),  
                 colours2=c("EEDheto_less"="blue",   "EEDko_less"="blue4" ),   
                 path2=subdir_2_part6,     fileName2="1-two2-NTR-less",  
                 title2="Unchanged Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR2_two2_EEDheto_more,  column_NTR2_two2_EEDko_more ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_more", numOfColumns1),   rep("EEDko_more", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_more",   "EEDko_more"  ),  
                 colours2=c("EEDheto_more"="blue",   "EEDko_more"="blue4" ),   
                 path2=subdir_2_part6,     fileName2="1-two2-NTR-more",  
                 title2="Unchanged Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyBoxViolinPlot_1(vector2=c(row_NTR2_two2_EEDheto_same,  row_NTR2_two2_EEDko_same,  row_NTR2_two2_EEDheto_less, row_NTR2_two2_EEDko_less,
                            row_NTR2_two2_EEDheto_more,  row_NTR2_two2_EEDko_more  ),  
                  sampleType2=c( rep("EEDheto_same", two2_lenNTR_A),   rep("EEDko_same", two2_lenNTR_A),  rep("EEDheto_less", two2_lenNTR_B),   rep("EEDko_less", two2_lenNTR_B), 
                                 rep("EEDheto_more", two2_lenNTR_C),   rep("EEDko_more", two2_lenNTR_C)  ), 
                  sampleRank2=c( "EEDheto_same",   "EEDko_same",  "EEDheto_less",  "EEDko_less",  "EEDheto_more",  "EEDko_more"  ),  
                  colours2=c("red",  "red4",   "blue",    "blue4", "green", "green4" ),   
                  path2=subdir_2_part6,  fileName2="1-two2-NTR-boxViolin",  
                  title2="Unchanged Genes",  xLab2="Samples",  yLab2="NTR",   
                  height2=3.88,   width2=4, Ymin2=-0.3, Ymax2=1 )




























########################################################################
subdir_3_part6 <- paste(Part6_g,  "/3-three3-diff-NTR", sep = "")
if( ! file.exists(subdir_3_part6) ) { dir.create(subdir_3_part6) }

length(row_Average_three3_EEDheto_NTR2)       
length(row_Average_three3_EEDko_NTR2)         
row_Average_three3_diffNTR <- (row_Average_three3_EEDheto_NTR2 - row_Average_three3_EEDko_NTR2 )
length(row_Average_three3_diffNTR)
summary(row_Average_three3_diffNTR)

three3_NTR_threshold <- 0.05
length( three3_week0_dist[three3_week0_dist < -three3_NTR_threshold] )  ## NTR of three3_week0_EEDheto is less than three3_week0_EEDko
length( three3_week0_dist[three3_week0_dist >  three3_NTR_threshold] )  ## NTR of three3_week0_EEDheto is more than three3_week0_EEDko



## 差异不显著
three3_NTR_diffIndex_1 <- ( abs(row_Average_three3_diffNTR) <= three3_NTR_threshold  )
length(three3_NTR_diffIndex_1[three3_NTR_diffIndex_1])
write.table(x=three3_NTR_diffIndex_1,   file = paste(subdir_3_part6, "/1-three3-NoDiff.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")


## 差异显著  three3_week0_EEDheto is less than three3_week0_EEDko
three3_NTR_diffIndex_2 <- ( row_Average_three3_diffNTR < -three3_NTR_threshold )
length(three3_NTR_diffIndex_2[three3_NTR_diffIndex_2])
write.table(x=three3_NTR_diffIndex_2,   file = paste(subdir_3_part6, "/1-three3-EEDheto-Less-Than-EEDko.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")


## 差异显著  three3_week0_EEDheto is more than three3_week0_EEDko
three3_NTR_diffIndex_3 <- ( row_Average_three3_diffNTR >  three3_NTR_threshold )
length(three3_NTR_diffIndex_3[three3_NTR_diffIndex_3])
write.table(x=three3_NTR_diffIndex_3,   file = paste(subdir_3_part6, "/1-three3-EEDheto-More-Than-EEDko.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")






MyAverageLines_3(vector2=c(colMeans(Average_three3_week0_EEDheto[three3_NTR_diffIndex_1, ]),  colMeans(Average_three3_week0_EEDko[three3_NTR_diffIndex_1, ]), 
                           colMeans(Average_three3_week4_EEDheto[three3_NTR_diffIndex_1, ]),  colMeans(Average_three3_week4_EEDko[three3_NTR_diffIndex_1, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_3_part6,     fileName2="A-three3-NoDiffNTR",  
                 title2="Up-regulated Genes (No DiffNTR)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=-0.3,   Ymax2=1,     height2=3.2,   width2=5.55 , center2="TSS" )


three3_lenNTR_A <- nrow(Average_three3_week0_EEDheto[three3_NTR_diffIndex_1, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_three3_week0_EEDheto[three3_NTR_diffIndex_1, ]),  rowMeans(Average_three3_week0_EEDko[three3_NTR_diffIndex_1, ]), 
                            rowMeans(Average_three3_week4_EEDheto[three3_NTR_diffIndex_1, ]),  rowMeans(Average_three3_week4_EEDko[three3_NTR_diffIndex_1, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", three3_lenNTR_A),   rep("week0_EEDko",  three3_lenNTR_A), 
                                 rep("week4_EEDheto", three3_lenNTR_A),   rep("week4_EEDko",  three3_lenNTR_A) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part6,  fileName2="A-three3-all-NoDiffNTR-BoxViolin",  
                  title2="Up-regulated Genes  (No DiffNTR)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_three3_week0_EEDheto[three3_NTR_diffIndex_1, 201:400]),  rowMeans(Average_three3_week0_EEDko[three3_NTR_diffIndex_1, 201:400]), 
                            rowMeans(Average_three3_week4_EEDheto[three3_NTR_diffIndex_1, 201:400]),  rowMeans(Average_three3_week4_EEDko[three3_NTR_diffIndex_1, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", three3_lenNTR_A),   rep("week0_EEDko",  three3_lenNTR_A), 
                                 rep("week4_EEDheto", three3_lenNTR_A),   rep("week4_EEDko",  three3_lenNTR_A) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part6,  fileName2="A-three3-2kbRegion-NoDiffNTR-BoxViolin",  
                  title2="Up-regulated Genes  (No DiffNTR)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )










MyAverageLines_3(vector2=c(colMeans(Average_three3_week0_EEDheto[three3_NTR_diffIndex_2, ]),  colMeans(Average_three3_week0_EEDko[three3_NTR_diffIndex_2, ]), 
                           colMeans(Average_three3_week4_EEDheto[three3_NTR_diffIndex_2, ]),  colMeans(Average_three3_week4_EEDko[three3_NTR_diffIndex_2, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_3_part6,     fileName2="B-three3-EEDhetoLess",  
                 title2="Up-regulated Genes (EEDheto Less)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=-0.3,   Ymax2=1,     height2=3.2,   width2=5.55 , center2="TSS" )


three3_lenNTR_B <- nrow(Average_three3_week0_EEDheto[three3_NTR_diffIndex_2, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_three3_week0_EEDheto[three3_NTR_diffIndex_2, ]),  rowMeans(Average_three3_week0_EEDko[three3_NTR_diffIndex_2, ]), 
                            rowMeans(Average_three3_week4_EEDheto[three3_NTR_diffIndex_2, ]),  rowMeans(Average_three3_week4_EEDko[three3_NTR_diffIndex_2, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", three3_lenNTR_B),   rep("week0_EEDko",  three3_lenNTR_B), 
                                 rep("week4_EEDheto", three3_lenNTR_B),   rep("week4_EEDko",  three3_lenNTR_B) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part6,  fileName2="B-three3-all-EEDhetoLess-BoxViolin",  
                  title2="Up-regulated Genes  (EEDheto Less)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_three3_week0_EEDheto[three3_NTR_diffIndex_2, 201:400]),  rowMeans(Average_three3_week0_EEDko[three3_NTR_diffIndex_2, 201:400]), 
                            rowMeans(Average_three3_week4_EEDheto[three3_NTR_diffIndex_2, 201:400]),  rowMeans(Average_three3_week4_EEDko[three3_NTR_diffIndex_2, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", three3_lenNTR_B),   rep("week0_EEDko",  three3_lenNTR_B), 
                                 rep("week4_EEDheto", three3_lenNTR_B),   rep("week4_EEDko",  three3_lenNTR_B) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part6,  fileName2="B-three3-2kbRegion-EEDhetoLess-BoxViolin",  
                  title2="Up-regulated Genes  (EEDheto Less)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )












MyAverageLines_3(vector2=c(colMeans(Average_three3_week0_EEDheto[three3_NTR_diffIndex_3, ]),  colMeans(Average_three3_week0_EEDko[three3_NTR_diffIndex_3, ]), 
                           colMeans(Average_three3_week4_EEDheto[three3_NTR_diffIndex_3, ]),  colMeans(Average_three3_week4_EEDko[three3_NTR_diffIndex_3, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_3_part6,     fileName2="C-three3-EEDhetoMore",  
                 title2="Up-regulated Genes (EEDheto More)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=-0.3,   Ymax2=1,    height2=3.2,   width2=5.55 , center2="TSS" )


three3_lenNTR_C <- nrow(Average_three3_week0_EEDheto[three3_NTR_diffIndex_3, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_three3_week0_EEDheto[three3_NTR_diffIndex_3, ]),  rowMeans(Average_three3_week0_EEDko[three3_NTR_diffIndex_3, ]), 
                            rowMeans(Average_three3_week4_EEDheto[three3_NTR_diffIndex_3, ]),  rowMeans(Average_three3_week4_EEDko[three3_NTR_diffIndex_3, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", three3_lenNTR_C),   rep("week0_EEDko",  three3_lenNTR_C), 
                                 rep("week4_EEDheto", three3_lenNTR_C),   rep("week4_EEDko",  three3_lenNTR_C) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part6,  fileName2="C-three3-all-EEDhetoMore-BoxViolin",  
                  title2="Up-regulated Genes  (EEDheto More)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_three3_week0_EEDheto[three3_NTR_diffIndex_3, 201:400]),  rowMeans(Average_three3_week0_EEDko[three3_NTR_diffIndex_3, 201:400]), 
                            rowMeans(Average_three3_week4_EEDheto[three3_NTR_diffIndex_3, 201:400]),  rowMeans(Average_three3_week4_EEDko[three3_NTR_diffIndex_3, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", three3_lenNTR_C),   rep("week0_EEDko",  three3_lenNTR_C), 
                                 rep("week4_EEDheto", three3_lenNTR_C),   rep("week4_EEDko",  three3_lenNTR_C) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part6,  fileName2="C-three3-2kbRegion-EEDhetoMore-BoxViolin",  
                  title2="Up-regulated Genes  (EEDheto More)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )










column_NTR2_three3_EEDheto_same   <-  log( colMeans(Average_three3_week0_EEDheto[three3_NTR_diffIndex_1, ]) / colMeans(Average_three3_week4_EEDheto[three3_NTR_diffIndex_1, ]) )/4
column_NTR2_three3_EEDko_same     <-  log( colMeans(Average_three3_week0_EEDko[three3_NTR_diffIndex_1, ])   / colMeans(Average_three3_week4_EEDko[three3_NTR_diffIndex_1, ])   )/4
column_NTR2_three3_EEDheto_less   <-  log( colMeans(Average_three3_week0_EEDheto[three3_NTR_diffIndex_2, ]) / colMeans(Average_three3_week4_EEDheto[three3_NTR_diffIndex_2, ]) )/4
column_NTR2_three3_EEDko_less     <-  log( colMeans(Average_three3_week0_EEDko[three3_NTR_diffIndex_2, ])   / colMeans(Average_three3_week4_EEDko[three3_NTR_diffIndex_2, ])   )/4
column_NTR2_three3_EEDheto_more   <-  log( colMeans(Average_three3_week0_EEDheto[three3_NTR_diffIndex_3, ]) / colMeans(Average_three3_week4_EEDheto[three3_NTR_diffIndex_3, ]) )/4
column_NTR2_three3_EEDko_more     <-  log( colMeans(Average_three3_week0_EEDko[three3_NTR_diffIndex_3, ])   / colMeans(Average_three3_week4_EEDko[three3_NTR_diffIndex_3, ])   )/4

row_NTR2_three3_EEDheto_same   <-  log( row_Average_three3_week0_EEDheto[three3_NTR_diffIndex_1] / row_Average_three3_week4_EEDheto[three3_NTR_diffIndex_1] )/4
row_NTR2_three3_EEDko_same     <-  log( row_Average_three3_week0_EEDko[three3_NTR_diffIndex_1]   / row_Average_three3_week4_EEDko[three3_NTR_diffIndex_1]   )/4
row_NTR2_three3_EEDheto_less   <-  log( row_Average_three3_week0_EEDheto[three3_NTR_diffIndex_2] / row_Average_three3_week4_EEDheto[three3_NTR_diffIndex_2] )/4
row_NTR2_three3_EEDko_less     <-  log( row_Average_three3_week0_EEDko[three3_NTR_diffIndex_2]   / row_Average_three3_week4_EEDko[three3_NTR_diffIndex_2]   )/4
row_NTR2_three3_EEDheto_more   <-  log( row_Average_three3_week0_EEDheto[three3_NTR_diffIndex_3] / row_Average_three3_week4_EEDheto[three3_NTR_diffIndex_3] )/4
row_NTR2_three3_EEDko_more     <-  log( row_Average_three3_week0_EEDko[three3_NTR_diffIndex_3]   / row_Average_three3_week4_EEDko[three3_NTR_diffIndex_3]   )/4


MyAverageLines_3(vector2=c(column_NTR2_three3_EEDheto_same,  column_NTR2_three3_EEDko_same,  column_NTR2_three3_EEDheto_less, column_NTR2_three3_EEDko_less,
                           column_NTR2_three3_EEDheto_more, column_NTR2_three3_EEDko_more ),   
                 numSample2=6,   
                 sampleType2=c( rep("EEDheto_same", numOfColumns1),   rep("EEDko_same", numOfColumns1),  rep("EEDheto_less", numOfColumns1),   rep("EEDko_less", numOfColumns1), 
                                rep("EEDheto_more", numOfColumns1),   rep("EEDko_more", numOfColumns1)  ), 
                 sampleRank2=c( "EEDheto_same",   "EEDko_same",  "EEDheto_less",  "EEDko_less",  "EEDheto_more",  "EEDko_more"  ),  
                 colours2=c("EEDheto_same"="red",   "EEDko_same"="red4",  "EEDheto_less"="blue",  "EEDko_less"="blue4",  "EEDheto_more"="green",  "EEDko_more"="green4"),   
                 path2=subdir_3_part6,     fileName2="1-three3-NTR",  
                 title2="Up-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR2_three3_EEDheto_same,  column_NTR2_three3_EEDko_same ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_same", numOfColumns1),   rep("EEDko_same", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_same",   "EEDko_same"  ),  
                 colours2=c("EEDheto_same"="red",   "EEDko_same"="red4" ),   
                 path2=subdir_3_part6,     fileName2="1-three3-NTR-same",  
                 title2="Up-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR2_three3_EEDheto_less,  column_NTR2_three3_EEDko_less ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_less", numOfColumns1),   rep("EEDko_less", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_less",   "EEDko_less"  ),  
                 colours2=c("EEDheto_less"="blue",   "EEDko_less"="blue4" ),   
                 path2=subdir_3_part6,     fileName2="1-three3-NTR-less",  
                 title2="Up-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR2_three3_EEDheto_more,  column_NTR2_three3_EEDko_more ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_more", numOfColumns1),   rep("EEDko_more", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_more",   "EEDko_more"  ),  
                 colours2=c("EEDheto_more"="blue",   "EEDko_more"="blue4" ),   
                 path2=subdir_3_part6,     fileName2="1-three3-NTR-more",  
                 title2="Up-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyBoxViolinPlot_1(vector2=c(row_NTR2_three3_EEDheto_same,  row_NTR2_three3_EEDko_same,  row_NTR2_three3_EEDheto_less, row_NTR2_three3_EEDko_less,
                            row_NTR2_three3_EEDheto_more,  row_NTR2_three3_EEDko_more  ),  
                  sampleType2=c( rep("EEDheto_same", three3_lenNTR_A),   rep("EEDko_same", three3_lenNTR_A),  rep("EEDheto_less", three3_lenNTR_B),   rep("EEDko_less", three3_lenNTR_B), 
                                 rep("EEDheto_more", three3_lenNTR_C),   rep("EEDko_more", three3_lenNTR_C)  ), 
                  sampleRank2=c( "EEDheto_same",   "EEDko_same",  "EEDheto_less",  "EEDko_less",  "EEDheto_more",  "EEDko_more"  ),  
                  colours2=c("red",  "red4",   "blue",    "blue4", "green", "green4" ),   
                  path2=subdir_3_part6,  fileName2="1-three3-NTR-boxViolin",  
                  title2="Up-regulated Genes",  xLab2="Samples",  yLab2="NTR",   
                  height2=3.88,   width2=4, Ymin2=-0.3, Ymax2=1 )








####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################








