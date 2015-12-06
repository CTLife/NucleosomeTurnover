#############################################################################################################################
## Part  4:  Figures about nucleosome turnover rate (NTR). Compute the correlation between NOL and NTR.
#############################################################################################################################





#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
subdir_1_part4 <- paste(Part4_g,  "/1-NTRfor3classes", sep = "")
if( ! file.exists(subdir_1_part4) ) { dir.create(subdir_1_part4) }



###############  average columns
column_Average_one1_EEDheto_NTR       <-  log( column_Average_one1_week0_EEDheto/column_Average_one1_week4_EEDheto)/4
column_Average_one1_EEDko_NTR         <-  log( column_Average_one1_week0_EEDko/column_Average_one1_week4_EEDko )/4
column_Average_two2_EEDheto_NTR       <-  log( column_Average_two2_week0_EEDheto/column_Average_two2_week4_EEDheto )/4
column_Average_two2_EEDko_NTR         <-  log( column_Average_two2_week0_EEDko/column_Average_two2_week4_EEDko )/4
column_Average_three3_EEDheto_NTR     <-  log( column_Average_three3_week0_EEDheto/column_Average_three3_week4_EEDheto )/4
column_Average_three3_EEDko_NTR       <-  log( column_Average_three3_week0_EEDko/column_Average_three3_week4_EEDko )/4

column_Average_one1_EEDheto_halflife      <- log(2)/column_Average_one1_EEDheto_NTR
column_Average_one1_EEDko_halflife        <- log(2)/column_Average_one1_EEDko_NTR
column_Average_two2_EEDheto_halflife      <- log(2)/column_Average_two2_EEDheto_NTR
column_Average_two2_EEDko_halflife        <- log(2)/column_Average_two2_EEDko_NTR
column_Average_three3_EEDheto_halflife    <- log(2)/column_Average_three3_EEDheto_NTR
column_Average_three3_EEDko_halflife      <- log(2)/column_Average_three3_EEDko_NTR


NTR_min1 <- min(column_Average_one1_EEDheto_NTR, column_Average_one1_EEDko_NTR, column_Average_two2_EEDheto_NTR, column_Average_two2_EEDko_NTR, column_Average_three3_EEDheto_NTR, column_Average_three3_EEDko_NTR)
NTR_max1 <- max(column_Average_one1_EEDheto_NTR, column_Average_one1_EEDko_NTR, column_Average_two2_EEDheto_NTR, column_Average_two2_EEDko_NTR, column_Average_three3_EEDheto_NTR, column_Average_three3_EEDko_NTR)

halflife_min1 <- min(column_Average_one1_EEDheto_halflife, column_Average_one1_EEDko_halflife, column_Average_two2_EEDheto_halflife, column_Average_two2_EEDko_halflife, column_Average_three3_EEDheto_halflife, column_Average_three3_EEDko_halflife)
halflife_max1 <- max(column_Average_one1_EEDheto_halflife, column_Average_one1_EEDko_halflife, column_Average_two2_EEDheto_halflife, column_Average_two2_EEDko_halflife, column_Average_three3_EEDheto_halflife, column_Average_three3_EEDko_halflife)

NTR_min1
NTR_max1
halflife_min1
halflife_max1


MyAverageLines_3(vector2=c(column_Average_one1_EEDheto_NTR, column_Average_one1_EEDko_NTR,      column_Average_two2_EEDheto_NTR, 
                           column_Average_two2_EEDko_NTR,   column_Average_three3_EEDheto_NTR,  column_Average_three3_EEDko_NTR ),   
                numSample2=6,   
                sampleType2=c( rep("Down_EEDheto", numOfColumns1),   rep("Down_EEDko", numOfColumns1),   rep("Unchange_EEDheto", numOfColumns1),  
                               rep("Unchange_EEDko", numOfColumns1), rep("Up_EEDheto", numOfColumns1),   rep("Up_EEDko", numOfColumns1)  ), 
                colours2=c("red", "red4",   "blue", "blue4",  "green", "green4"),   
                sampleRank2=c( "Down_EEDheto",  "Down_EEDko",  "Unchange_EEDheto", 
                               "Unchange_EEDko",  "Up_EEDheto",   "Up_EEDko"  ),     
                path2=subdir_1_part4,      fileName2="1-NTR-6curves",  
                title2="NTR for 3 kinds of genes",    xLab2="Relative distance (kb)",        yLab2="NTR",   
                Ymin2=NTR_min1,   Ymax2=NTR_max1,    height2=3.2,   width2=5.85,    center2="TSS" )




MyAverageLines_3(vector2=c(column_Average_one1_EEDheto_NTR, column_Average_one1_EEDko_NTR ),   
                 numSample2=2,   
                 sampleType2=c( rep("Down_EEDheto", numOfColumns1),  rep("Down_EEDko", numOfColumns1) ), 
                 colours2=c("red", "red4"),   
                 sampleRank2=c( "Down_EEDheto",  "Down_EEDko"  ),     
                 path2=subdir_1_part4,      fileName2="2-down-NTR-2curves",  
                 title2="NTR for Down-regulated genes",    xLab2="Relative distance (kb)",        yLab2="NTR",   
                 Ymin2=NTR_min1,   Ymax2=NTR_max1,    height2=3.2,   width2=5.55, center2="TSS" )


MyAverageLines_3(vector2=c( column_Average_two2_EEDheto_NTR,    column_Average_two2_EEDko_NTR ),   
                 numSample2=2,   
                 sampleType2=c( rep("Unchange_EEDheto", numOfColumns1),  rep("Unchange_EEDko", numOfColumns1) ), 
                 colours2=c(   "blue", "blue4"),   
                 sampleRank2=c(  "Unchange_EEDheto", "Unchange_EEDko"  ),     
                 path2=subdir_1_part4,      fileName2="3-retain-NTR-2curves",  
                 title2="NTR for unchanged genes",    xLab2="Relative distance (kb)",        yLab2="NTR",   
                 Ymin2=NTR_min1,   Ymax2=NTR_max1,    height2=3.2,   width2=5.85, center2="TSS" )


MyAverageLines_3(vector2=c(column_Average_three3_EEDheto_NTR,  column_Average_three3_EEDko_NTR ),   
                 numSample2=2,   
                 sampleType2=c(rep("Up_EEDheto", numOfColumns1),  rep("Up_EEDko", numOfColumns1)  ), 
                 colours2=c("green", "green4"),   
                 sampleRank2=c(  "Up_EEDheto",   "Up_EEDko"  ),     
                 path2=subdir_1_part4,      fileName2="4-up-NTR-2curves",  
                 title2="NTR for up-regulated genes",    xLab2="Relative distance (kb)",        yLab2="NTR",   
                 Ymin2=NTR_min1,   Ymax2=NTR_max1,    height2=3.2,   width2=5.55, center2="TSS" )















###############  average rows
row_Average_one1_EEDheto_NTR2       <-  log( row_Average_one1_week0_EEDheto/row_Average_one1_week4_EEDheto)/4
row_Average_one1_EEDko_NTR2         <-  log( row_Average_one1_week0_EEDko/row_Average_one1_week4_EEDko )/4
row_Average_two2_EEDheto_NTR2       <-  log( row_Average_two2_week0_EEDheto/row_Average_two2_week4_EEDheto )/4
row_Average_two2_EEDko_NTR2         <-  log( row_Average_two2_week0_EEDko/row_Average_two2_week4_EEDko )/4
row_Average_three3_EEDheto_NTR2     <-  log( row_Average_three3_week0_EEDheto/row_Average_three3_week4_EEDheto )/4
row_Average_three3_EEDko_NTR2       <-  log( row_Average_three3_week0_EEDko/row_Average_three3_week4_EEDko )/4


MyBoxViolinPlot_1(vector2=c(row_Average_one1_EEDheto_NTR2,  row_Average_one1_EEDko_NTR2,       row_Average_two2_EEDheto_NTR2,  
                            row_Average_two2_EEDko_NTR2,    row_Average_three3_EEDheto_NTR2,   row_Average_three3_EEDko_NTR2),  
                  sampleType2=c( rep("Down_EEDheto", numOfRows1),       rep("Down_EEDko", numOfRows1),  
                                 rep("Unchange_EEDheto", numOfRows2),   rep("Unchange_EEDko", numOfRows2),
                                 rep("Up_EEDheto", numOfRows3),         rep("Up_EEDko", numOfRows3) ),   
                  sampleRank2=c( "Down_EEDheto",  "Down_EEDko",  "Unchange_EEDheto", 
                                 "Unchange_EEDko",  "Up_EEDheto",   "Up_EEDko"  ),  
                  colours2=c("red", "red4",   "blue", "blue4",  "green", "green4"),   
                  path2=subdir_1_part4,  fileName2="5-NTR-boxViolin",  
                  title2="NTR for 3 kinds of genes",  xLab2="Samples",  yLab2="NTR",   
                  height2=3.88,   width2=4, Ymin2=0, Ymax2=0.6 )


MyBoxViolinPlot_3_s6(vector2=c(row_Average_one1_EEDheto_NTR2,  row_Average_one1_EEDko_NTR2,       row_Average_two2_EEDheto_NTR2,  
                               row_Average_two2_EEDko_NTR2,    row_Average_three3_EEDheto_NTR2,   row_Average_three3_EEDko_NTR2),  
                sampleType2=c( rep("Down_EEDheto", numOfRows1),       rep("Down_EEDko", numOfRows1),  
                                    rep("Unchange_EEDheto", numOfRows2),   rep("Unchange_EEDko", numOfRows2),
                                    rep("Up_EEDheto", numOfRows3),         rep("Up_EEDko", numOfRows3) ),   
                sampleRank2=c( "Down_EEDheto",  "Down_EEDko",  "Unchange_EEDheto", 
                                    "Unchange_EEDko",  "Up_EEDheto",   "Up_EEDko"  ), 
                colours2=c("red", "red4",   "blue", "blue4",  "green", "green4"),   
                path2=subdir_1_part4,  fileName2="7-NTR-6curves",  
                title2="NTR for 3 kinds of genes",  xLab2="Samples",  yLab2="NTR",   
                height2=3.88,   width2=5, Ymin2=0, Ymax2=0.6 )


MyHistogram_7(vector2=c(row_Average_one1_EEDheto_NTR2,  row_Average_one1_EEDko_NTR2,       row_Average_two2_EEDheto_NTR2,  
                        row_Average_two2_EEDko_NTR2,    row_Average_three3_EEDheto_NTR2,   row_Average_three3_EEDko_NTR2),  
              sampleType2=c( rep("Down_EEDheto", numOfRows1),       rep("Down_EEDko", numOfRows1),  
                             rep("Unchange_EEDheto", numOfRows2),   rep("Unchange_EEDko", numOfRows2),
                             rep("Up_EEDheto", numOfRows3),         rep("Up_EEDko", numOfRows3) ), 
              colours2=c("red", "red4",   "blue", "blue4",  "green", "green4"),      
              path2=subdir_1_part4,    
              fileName2="10-NTR-6curves",   
              title2="NTR for 3 kinds of genes",    xLab2="NTR",  
              height2=3.2,  width2=5.6,   xMin2=0,  xMax2=0.6,  yMin2=0,  yMax2=5)


MyHistogram_7(vector2=c(row_Average_one1_EEDheto_NTR2,  row_Average_one1_EEDko_NTR2 ),  
              sampleType2=c( rep("Down_EEDheto", numOfRows1),     rep("Down_EEDko", numOfRows1)  ),  
              colours2=c("red", "red4"  ),      
              path2=subdir_1_part4,    
              fileName2="11-NTR-2curves",   
              title2="NTR for Down-regulated genes",    xLab2="NTR",  
              height2=3.2,  width2=5.6,   xMin2=0,  xMax2=0.6,  yMin2=0,  yMax2=5)

MyHistogram_7(vector2=c( row_Average_two2_EEDheto_NTR2,  row_Average_two2_EEDko_NTR2 ),  
              sampleType2=c(  rep("Unchange_EEDheto", numOfRows2),   rep("Unchange_EEDko", numOfRows2) ),  
              colours2=c(   "blue", "blue4" ),      
              path2=subdir_1_part4,    
              fileName2="12-NTR-2curves",   
              title2="NTR for unchanged genes",    xLab2="NTR",  
              height2=3.2,  width2=5.6,   xMin2=0,  xMax2=0.6,  yMin2=0,  yMax2=5)

MyHistogram_7(vector2=c(   row_Average_three3_EEDheto_NTR2, row_Average_three3_EEDko_NTR2),  
              sampleType2=c( rep("Up_EEDheto", numOfRows3),       rep("Up_EEDko", numOfRows3) ),  
              colours2=c(  "green", "green4"),      
              path2=subdir_1_part4,    
              fileName2="13-NTR-2curves",   
              title2="NTR for up-regulated genes",    xLab2="NTR",  
              height2=3.2,  width2=5.6,   xMin2=0,  xMax2=0.6,  yMin2=0,  yMax2=5)












##########################################################################
subdir_2_part4 <- paste(Part4_g,  "/2-Rows5categaries", sep = "")
if( ! file.exists(subdir_2_part4) ) { dir.create(subdir_2_part4) }



one1_EEDheto_NTR_AA     <-  log( reduceRow1_Average_one1_week0_EEDheto[1, ] / reduceRow1_Average_one1_week4_EEDheto[1, ] ) / 4
one1_EEDko_NTR_AA       <-  log( reduceRow1_Average_one1_week0_EEDko[1, ] / reduceRow1_Average_one1_week4_EEDko[1, ] ) / 4
two2_EEDheto_NTR_AA     <-  log( reduceRow1_Average_two2_week0_EEDheto[1, ] / reduceRow1_Average_two2_week4_EEDheto[1, ] ) / 4
two2_EEDko_NTR_AA       <-  log( reduceRow1_Average_two2_week0_EEDko[1, ] / reduceRow1_Average_two2_week4_EEDko[1, ] ) / 4
three3_EEDheto_NTR_AA   <-  log( reduceRow1_Average_three3_week0_EEDheto[1, ] / reduceRow1_Average_three3_week4_EEDheto[1, ] ) / 4
three3_EEDko_NTR_AA     <-  log( reduceRow1_Average_three3_week0_EEDko[1, ] / reduceRow1_Average_three3_week4_EEDko[1, ] ) / 4

one1_EEDheto_halflife_AA      <- log(2)/one1_EEDheto_NTR_AA
one1_EEDko_halflifea_AA       <- log(2)/one1_EEDko_NTR_AA
two2_EEDheto_halflife_AA      <- log(2)/two2_EEDheto_NTR_AA
two2_EEDko_halflife_AA        <- log(2)/two2_EEDko_NTR_AA
three3_EEDheto_halflife_AA    <- log(2)/three3_EEDheto_NTR_AA
three3_EEDko_halflife_AA      <- log(2)/three3_EEDko_NTR_AA

NTR_min_AA <- min(one1_EEDheto_NTR_AA, one1_EEDko_NTR_AA, two2_EEDheto_NTR_AA, two2_EEDko_NTR_AA, three3_EEDheto_NTR_AA, three3_EEDko_NTR_AA)
NTR_max_AA <- max(one1_EEDheto_NTR_AA, one1_EEDko_NTR_AA, two2_EEDheto_NTR_AA, two2_EEDko_NTR_AA, three3_EEDheto_NTR_AA, three3_EEDko_NTR_AA)
halflife_min_AA <- min(one1_EEDheto_halflife_AA, one1_EEDko_halflifea_AA, two2_EEDheto_halflife_AA, two2_EEDko_halflife_AA, three3_EEDheto_halflife_AA, three3_EEDko_halflife_AA)
halflife_max_AA <- max(one1_EEDheto_halflife_AA, one1_EEDko_halflifea_AA, two2_EEDheto_halflife_AA, two2_EEDko_halflife_AA, three3_EEDheto_halflife_AA, three3_EEDko_halflife_AA)
NTR_min_AA
NTR_max_AA
halflife_min_AA
halflife_max_AA




one1_EEDheto_NTR_BB     <-  log( reduceRow1_Average_one1_week0_EEDheto[2, ] / reduceRow1_Average_one1_week4_EEDheto[2, ] ) / 4
one1_EEDko_NTR_BB       <-  log( reduceRow1_Average_one1_week0_EEDko[2, ] / reduceRow1_Average_one1_week4_EEDko[2, ] ) / 4
two2_EEDheto_NTR_BB     <-  log( reduceRow1_Average_two2_week0_EEDheto[2, ] / reduceRow1_Average_two2_week4_EEDheto[2, ] ) / 4
two2_EEDko_NTR_BB       <-  log( reduceRow1_Average_two2_week0_EEDko[2, ] / reduceRow1_Average_two2_week4_EEDko[2, ] ) / 4
three3_EEDheto_NTR_BB   <-  log( reduceRow1_Average_three3_week0_EEDheto[2, ] / reduceRow1_Average_three3_week4_EEDheto[2, ] ) / 4
three3_EEDko_NTR_BB     <-  log( reduceRow1_Average_three3_week0_EEDko[2, ] / reduceRow1_Average_three3_week4_EEDko[2, ] ) / 4

one1_EEDheto_halflife_BB      <- log(2)/one1_EEDheto_NTR_BB
one1_EEDko_halflifea_BB       <- log(2)/one1_EEDko_NTR_BB
two2_EEDheto_halflife_BB      <- log(2)/two2_EEDheto_NTR_BB
two2_EEDko_halflife_BB        <- log(2)/two2_EEDko_NTR_BB
three3_EEDheto_halflife_BB    <- log(2)/three3_EEDheto_NTR_BB
three3_EEDko_halflife_BB      <- log(2)/three3_EEDko_NTR_BB

NTR_min_BB <- min(one1_EEDheto_NTR_BB, one1_EEDko_NTR_BB, two2_EEDheto_NTR_BB, two2_EEDko_NTR_BB, three3_EEDheto_NTR_BB, three3_EEDko_NTR_BB)
NTR_max_BB <- max(one1_EEDheto_NTR_BB, one1_EEDko_NTR_BB, two2_EEDheto_NTR_BB, two2_EEDko_NTR_BB, three3_EEDheto_NTR_BB, three3_EEDko_NTR_BB)
halflife_min_BB <- min(one1_EEDheto_halflife_BB, one1_EEDko_halflifea_BB, two2_EEDheto_halflife_BB, two2_EEDko_halflife_BB, three3_EEDheto_halflife_BB, three3_EEDko_halflife_BB)
halflife_max_BB <- max(one1_EEDheto_halflife_BB, one1_EEDko_halflifea_BB, two2_EEDheto_halflife_BB, two2_EEDko_halflife_BB, three3_EEDheto_halflife_BB, three3_EEDko_halflife_BB)
NTR_min_BB
NTR_max_BB
halflife_min_BB
halflife_max_BB







one1_EEDheto_NTR_CC     <-  log( reduceRow1_Average_one1_week0_EEDheto[3, ] / reduceRow1_Average_one1_week4_EEDheto[3, ] ) / 4
one1_EEDko_NTR_CC       <-  log( reduceRow1_Average_one1_week0_EEDko[3, ] / reduceRow1_Average_one1_week4_EEDko[3, ] ) / 4
two2_EEDheto_NTR_CC     <-  log( reduceRow1_Average_two2_week0_EEDheto[3, ] / reduceRow1_Average_two2_week4_EEDheto[3, ] ) / 4
two2_EEDko_NTR_CC       <-  log( reduceRow1_Average_two2_week0_EEDko[3, ] / reduceRow1_Average_two2_week4_EEDko[3, ] ) / 4
three3_EEDheto_NTR_CC   <-  log( reduceRow1_Average_three3_week0_EEDheto[3, ] / reduceRow1_Average_three3_week4_EEDheto[3, ] ) / 4
three3_EEDko_NTR_CC     <-  log( reduceRow1_Average_three3_week0_EEDko[3, ] / reduceRow1_Average_three3_week4_EEDko[3, ] ) / 4

one1_EEDheto_halflife_CC      <- log(2)/one1_EEDheto_NTR_CC
one1_EEDko_halflifea_CC       <- log(2)/one1_EEDko_NTR_CC
two2_EEDheto_halflife_CC      <- log(2)/two2_EEDheto_NTR_CC
two2_EEDko_halflife_CC        <- log(2)/two2_EEDko_NTR_CC
three3_EEDheto_halflife_CC    <- log(2)/three3_EEDheto_NTR_CC
three3_EEDko_halflife_CC     <- log(2)/three3_EEDko_NTR_CC

NTR_min_CC <- min(one1_EEDheto_NTR_CC, one1_EEDko_NTR_CC, two2_EEDheto_NTR_CC, two2_EEDko_NTR_CC, three3_EEDheto_NTR_CC, three3_EEDko_NTR_CC)
NTR_max_CC <- max(one1_EEDheto_NTR_CC, one1_EEDko_NTR_CC, two2_EEDheto_NTR_CC, two2_EEDko_NTR_CC, three3_EEDheto_NTR_CC, three3_EEDko_NTR_CC)
halflife_min_CC <- min(one1_EEDheto_halflife_CC, one1_EEDko_halflifea_CC, two2_EEDheto_halflife_CC, two2_EEDko_halflife_CC, three3_EEDheto_halflife_CC, three3_EEDko_halflife_CC)
halflife_max_CC <- max(one1_EEDheto_halflife_CC, one1_EEDko_halflifea_CC, two2_EEDheto_halflife_CC, two2_EEDko_halflife_CC, three3_EEDheto_halflife_CC, three3_EEDko_halflife_CC)
NTR_min_CC
NTR_max_CC
halflife_min_CC
halflife_max_CC









one1_EEDheto_NTR_DD     <-  log( reduceRow1_Average_one1_week0_EEDheto[4, ] / reduceRow1_Average_one1_week4_EEDheto[4, ] ) / 4
one1_EEDko_NTR_DD       <-  log( reduceRow1_Average_one1_week0_EEDko[4, ] / reduceRow1_Average_one1_week4_EEDko[4, ] ) / 4
two2_EEDheto_NTR_DD     <-  log( reduceRow1_Average_two2_week0_EEDheto[4, ] / reduceRow1_Average_two2_week4_EEDheto[4, ] ) / 4
two2_EEDko_NTR_DD       <-  log( reduceRow1_Average_two2_week0_EEDko[4, ] / reduceRow1_Average_two2_week4_EEDko[4, ] ) / 4
three3_EEDheto_NTR_DD   <-  log( reduceRow1_Average_three3_week0_EEDheto[4, ] / reduceRow1_Average_three3_week4_EEDheto[4, ] ) / 4
three3_EEDko_NTR_DD     <-  log( reduceRow1_Average_three3_week0_EEDko[4, ] / reduceRow1_Average_three3_week4_EEDko[4, ] ) / 4

one1_EEDheto_halflife_DD      <- log(2)/one1_EEDheto_NTR_DD
one1_EEDko_halflifea_DD       <- log(2)/one1_EEDko_NTR_DD
two2_EEDheto_halflife_DD      <- log(2)/two2_EEDheto_NTR_DD
two2_EEDko_halflife_DD        <- log(2)/two2_EEDko_NTR_DD
three3_EEDheto_halflife_DD    <- log(2)/three3_EEDheto_NTR_DD
three3_EEDko_halflife_DD     <- log(2)/three3_EEDko_NTR_DD

NTR_min_DD <- min(one1_EEDheto_NTR_DD, one1_EEDko_NTR_DD, two2_EEDheto_NTR_DD, two2_EEDko_NTR_DD, three3_EEDheto_NTR_DD, three3_EEDko_NTR_DD)
NTR_max_DD <- max(one1_EEDheto_NTR_DD, one1_EEDko_NTR_DD, two2_EEDheto_NTR_DD, two2_EEDko_NTR_DD, three3_EEDheto_NTR_DD, three3_EEDko_NTR_DD)
halflife_min_DD <- min(one1_EEDheto_halflife_DD, one1_EEDko_halflifea_DD, two2_EEDheto_halflife_DD, two2_EEDko_halflife_DD, three3_EEDheto_halflife_DD, three3_EEDko_halflife_DD)
halflife_max_DD <- max(one1_EEDheto_halflife_DD, one1_EEDko_halflifea_DD, two2_EEDheto_halflife_DD, two2_EEDko_halflife_DD, three3_EEDheto_halflife_DD, three3_EEDko_halflife_DD)
NTR_min_DD
NTR_max_DD
halflife_min_DD
halflife_max_DD




one1_EEDheto_NTR_EE     <-  log( reduceRow1_Average_one1_week0_EEDheto[5, ] / reduceRow1_Average_one1_week4_EEDheto[5, ] ) / 4
one1_EEDko_NTR_EE       <-  log( reduceRow1_Average_one1_week0_EEDko[5, ] / reduceRow1_Average_one1_week4_EEDko[5, ] ) / 4
two2_EEDheto_NTR_EE     <-  log( reduceRow1_Average_two2_week0_EEDheto[5, ] / reduceRow1_Average_two2_week4_EEDheto[5, ] ) / 4
two2_EEDko_NTR_EE       <-  log( reduceRow1_Average_two2_week0_EEDko[5, ] / reduceRow1_Average_two2_week4_EEDko[5, ] ) / 4
three3_EEDheto_NTR_EE   <-  log( reduceRow1_Average_three3_week0_EEDheto[5, ] / reduceRow1_Average_three3_week4_EEDheto[5, ] ) / 4
three3_EEDko_NTR_EE     <-  log( reduceRow1_Average_three3_week0_EEDko[5, ] / reduceRow1_Average_three3_week4_EEDko[5, ] ) / 4

one1_EEDheto_halflife_EE      <- log(2)/one1_EEDheto_NTR_EE
one1_EEDko_halflifea_EE       <- log(2)/one1_EEDko_NTR_EE
two2_EEDheto_halflife_EE      <- log(2)/two2_EEDheto_NTR_EE
two2_EEDko_halflife_EE        <- log(2)/two2_EEDko_NTR_EE
three3_EEDheto_halflife_EE    <- log(2)/three3_EEDheto_NTR_EE
three3_EEDko_halflife_EE     <- log(2)/three3_EEDko_NTR_EE

NTR_min_EE <- min(one1_EEDheto_NTR_EE, one1_EEDko_NTR_EE, two2_EEDheto_NTR_EE, two2_EEDko_NTR_EE, three3_EEDheto_NTR_EE, three3_EEDko_NTR_EE)
NTR_max_EE <- max(one1_EEDheto_NTR_EE, one1_EEDko_NTR_EE, two2_EEDheto_NTR_EE, two2_EEDko_NTR_EE, three3_EEDheto_NTR_EE, three3_EEDko_NTR_EE)
halflife_min_EE <- min(one1_EEDheto_halflife_EE, one1_EEDko_halflifea_EE, two2_EEDheto_halflife_EE, two2_EEDko_halflife_EE, three3_EEDheto_halflife_EE, three3_EEDko_halflife_EE)
halflife_max_EE <- max(one1_EEDheto_halflife_EE, one1_EEDko_halflifea_EE, two2_EEDheto_halflife_EE, two2_EEDko_halflife_EE, three3_EEDheto_halflife_EE, three3_EEDko_halflife_EE)
NTR_min_EE
NTR_max_EE
halflife_min_EE
halflife_max_EE



NTR_min_GG <- min( NTR_min_EE, NTR_min_DD, NTR_min_CC, NTR_min_BB, NTR_min_AA )
NTR_max_GG <- max( NTR_max_EE, NTR_max_DD, NTR_max_CC, NTR_max_BB, NTR_max_AA )









MyAverageLines_3(vector2=c(one1_EEDheto_NTR_AA, one1_EEDko_NTR_AA, two2_EEDheto_NTR_AA, two2_EEDko_NTR_AA, three3_EEDheto_NTR_AA, three3_EEDko_NTR_AA ),   
                numSample2=6,   
                sampleType2=c( rep("Down_EEDheto", numOfColumns1),   rep("Down_EEDko", numOfColumns1),   rep("Unchange_EEDheto", numOfColumns1),  
                               rep("Unchange_EEDko", numOfColumns1), rep("Up_EEDheto", numOfColumns1),   rep("Up_EEDko", numOfColumns1)  ), 
                colours2=c("red", "red4",   "blue", "blue4",  "green", "green4"),   
                sampleRank2=c( "Down_EEDheto",  "Down_EEDko",  "Unchange_EEDheto", 
                                "Unchange_EEDko",  "Up_EEDheto",   "Up_EEDko"  ),     
                path2=subdir_2_part4,      fileName2="1A-NTR-6curves",  
                title2="NTR for 3 kinds of genes (Lowest)",    xLab2="Relative distance (kb)",        yLab2="NTR",   
                Ymin2=NTR_min_GG,   Ymax2=NTR_max_GG,    height2=3.2,   width2=5.85,    center2="TSS" )




MyAverageLines_3(vector2=c(one1_EEDheto_NTR_BB, one1_EEDko_NTR_BB, two2_EEDheto_NTR_BB, two2_EEDko_NTR_BB, three3_EEDheto_NTR_BB, three3_EEDko_NTR_BB ),   
                 numSample2=6,   
                 sampleType2=c( rep("Down_EEDheto", numOfColumns1),   rep("Down_EEDko", numOfColumns1),   rep("Unchange_EEDheto", numOfColumns1),  
                                rep("Unchange_EEDko", numOfColumns1), rep("Up_EEDheto", numOfColumns1),   rep("Up_EEDko", numOfColumns1)  ), 
                 colours2=c("red", "red4",   "blue", "blue4",  "green", "green4"),   
                 sampleRank2=c( "Down_EEDheto",  "Down_EEDko",  "Unchange_EEDheto", 
                                "Unchange_EEDko",  "Up_EEDheto",   "Up_EEDko"  ),     
                 path2=subdir_2_part4,      fileName2="1B-NTR-6curves",  
                 title2="NTR for 3 kinds of genes (Low)",    xLab2="Relative distance (kb)",        yLab2="NTR",   
                 Ymin2=NTR_min_GG,   Ymax2=NTR_max_GG,    height2=3.2,   width2=5.85,    center2="TSS" )



MyAverageLines_3(vector2=c(one1_EEDheto_NTR_CC, one1_EEDko_NTR_CC, two2_EEDheto_NTR_CC, two2_EEDko_NTR_CC, three3_EEDheto_NTR_CC, three3_EEDko_NTR_CC ),   
                 numSample2=6,   
                 sampleType2=c( rep("Down_EEDheto", numOfColumns1),   rep("Down_EEDko", numOfColumns1),   rep("Unchange_EEDheto", numOfColumns1),  
                                rep("Unchange_EEDko", numOfColumns1), rep("Up_EEDheto", numOfColumns1),   rep("Up_EEDko", numOfColumns1)  ), 
                 colours2=c("red", "red4",   "blue", "blue4",  "green", "green4"),   
                 sampleRank2=c( "Down_EEDheto",  "Down_EEDko",  "Unchange_EEDheto", 
                                "Unchange_EEDko",  "Up_EEDheto",   "Up_EEDko"  ),     
                 path2=subdir_2_part4,      fileName2="1C-NTR-6curves",  
                 title2="NTR for 3 kinds of genes (Medium)",    xLab2="Relative distance (kb)",        yLab2="NTR",   
                 Ymin2=NTR_min_GG,   Ymax2=NTR_max_GG,    height2=3.2,   width2=5.85,    center2="TSS" )





MyAverageLines_3(vector2=c(one1_EEDheto_NTR_DD, one1_EEDko_NTR_DD, two2_EEDheto_NTR_DD, two2_EEDko_NTR_DD, three3_EEDheto_NTR_DD, three3_EEDko_NTR_DD ),   
                 numSample2=6,   
                 sampleType2=c( rep("Down_EEDheto", numOfColumns1),   rep("Down_EEDko", numOfColumns1),   rep("Unchange_EEDheto", numOfColumns1),  
                                rep("Unchange_EEDko", numOfColumns1), rep("Up_EEDheto", numOfColumns1),   rep("Up_EEDko", numOfColumns1)  ), 
                 colours2=c("red", "red4",   "blue", "blue4",  "green", "green4"),   
                 sampleRank2=c( "Down_EEDheto",  "Down_EEDko",  "Unchange_EEDheto", 
                                "Unchange_EEDko",  "Up_EEDheto",   "Up_EEDko"  ),     
                 path2=subdir_2_part4,      fileName2="1D-NTR-6curves",  
                 title2="NTR for 3 kinds of genes (High)",    xLab2="Relative distance (kb)",        yLab2="NTR",   
                 Ymin2=NTR_min_GG,   Ymax2=NTR_max_GG,    height2=3.2,   width2=5.85,    center2="TSS" )






MyAverageLines_3(vector2=c(one1_EEDheto_NTR_EE, one1_EEDko_NTR_EE, two2_EEDheto_NTR_EE, two2_EEDko_NTR_EE, three3_EEDheto_NTR_EE, three3_EEDko_NTR_EE ),   
                 numSample2=6,   
                 sampleType2=c( rep("Down_EEDheto", numOfColumns1),   rep("Down_EEDko", numOfColumns1),   rep("Unchange_EEDheto", numOfColumns1),  
                                rep("Unchange_EEDko", numOfColumns1), rep("Up_EEDheto", numOfColumns1),   rep("Up_EEDko", numOfColumns1)  ), 
                 colours2=c("red", "red4",   "blue", "blue4",  "green", "green4"),   
                 sampleRank2=c( "Down_EEDheto",  "Down_EEDko",  "Unchange_EEDheto", 
                                "Unchange_EEDko",  "Up_EEDheto",   "Up_EEDko"  ),     
                 path2=subdir_2_part4,      fileName2="1E-NTR-6curves",  
                 title2="NTR for 3 kinds of genes (Highest)",    xLab2="Relative distance (kb)",        yLab2="NTR",   
                 Ymin2=NTR_min_GG,   Ymax2=NTR_max_GG,    height2=3.2,   width2=5.85,    center2="TSS" )











##########################################################################
subdir_3_part4 <- paste(Part4_g,  "/3-Rows5categaries-boxplot", sep = "")
if( ! file.exists(subdir_3_part4) ) { dir.create(subdir_3_part4) }

row_Average_one1      <- numClasses( c(1:numOfRows1), binNum1=5)
row_Average_one1_len1 <- length(row_Average_one1[[1]])
row_Average_one1_len2 <- length(row_Average_one1[[2]])
row_Average_one1_len3 <- length(row_Average_one1[[3]])
row_Average_one1_len4 <- length(row_Average_one1[[4]])
row_Average_one1_len5 <- length(row_Average_one1[[5]])

row_Average_two2      <- numClasses( c(1:numOfRows2), binNum1=5)
row_Average_two2_len1 <- length(row_Average_two2[[1]])
row_Average_two2_len2 <- length(row_Average_two2[[2]])
row_Average_two2_len3 <- length(row_Average_two2[[3]])
row_Average_two2_len4 <- length(row_Average_two2[[4]])
row_Average_two2_len5 <- length(row_Average_two2[[5]])

row_Average_three3      <- numClasses( c(1:numOfRows3), binNum1=5)
row_Average_three3_len1 <- length(row_Average_three3[[1]])
row_Average_three3_len2 <- length(row_Average_three3[[2]])
row_Average_three3_len3 <- length(row_Average_three3[[3]])
row_Average_three3_len4 <- length(row_Average_three3[[4]])
row_Average_three3_len5 <- length(row_Average_three3[[5]])




one1_EEDheto_NTR_AA2     <-  log( row_Average_one1_week0_EEDheto[ row_Average_one1[[1]] ] / row_Average_one1_week4_EEDheto[ row_Average_one1[[1]] ] ) / 4
one1_EEDko_NTR_AA2       <-  log( row_Average_one1_week0_EEDko[ row_Average_one1[[1]] ] / row_Average_one1_week4_EEDko[ row_Average_one1[[1]] ] ) / 4
two2_EEDheto_NTR_AA2     <-  log( row_Average_two2_week0_EEDheto[ row_Average_two2[[1]] ] / row_Average_two2_week4_EEDheto[ row_Average_two2[[1]] ] ) / 4
two2_EEDko_NTR_AA2       <-  log( row_Average_two2_week0_EEDko[ row_Average_two2[[1]] ] / row_Average_two2_week4_EEDko[ row_Average_two2[[1]] ] ) / 4
three3_EEDheto_NTR_AA2   <-  log( row_Average_three3_week0_EEDheto[ row_Average_three3[[1]] ] / row_Average_three3_week4_EEDheto[ row_Average_three3[[1]] ] ) / 4
three3_EEDko_NTR_AA2     <-  log( row_Average_three3_week0_EEDko[ row_Average_three3[[1]] ] / row_Average_three3_week4_EEDko[ row_Average_three3[[1]] ] ) / 4

one1_EEDheto_halflife_AA2      <- log(2)/one1_EEDheto_NTR_AA2
one1_EEDko_halflifea_AA2       <- log(2)/one1_EEDko_NTR_AA2
two2_EEDheto_halflife_AA2      <- log(2)/two2_EEDheto_NTR_AA2
two2_EEDko_halflife_AA2        <- log(2)/two2_EEDko_NTR_AA2
three3_EEDheto_halflife_AA2    <- log(2)/three3_EEDheto_NTR_AA2
three3_EEDko_halflife_AA2      <- log(2)/three3_EEDko_NTR_AA2

NTR_min_AA2 <- min(one1_EEDheto_NTR_AA2, one1_EEDko_NTR_AA2, two2_EEDheto_NTR_AA2, two2_EEDko_NTR_AA2, three3_EEDheto_NTR_AA2, three3_EEDko_NTR_AA2)
NTR_max_AA2 <- max(one1_EEDheto_NTR_AA2, one1_EEDko_NTR_AA2, two2_EEDheto_NTR_AA2, two2_EEDko_NTR_AA2, three3_EEDheto_NTR_AA2, three3_EEDko_NTR_AA2)
halflife_min_AA2 <- min(one1_EEDheto_halflife_AA2, one1_EEDko_halflifea_AA2, two2_EEDheto_halflife_AA2, two2_EEDko_halflife_AA2, three3_EEDheto_halflife_AA2, three3_EEDko_halflife_AA2)
halflife_max_AA2 <- max(one1_EEDheto_halflife_AA2, one1_EEDko_halflifea_AA2, two2_EEDheto_halflife_AA2, two2_EEDko_halflife_AA2, three3_EEDheto_halflife_AA2, three3_EEDko_halflife_AA2)
NTR_min_AA2
NTR_max_AA2
halflife_min_AA2
halflife_max_AA2









one1_EEDheto_NTR_BB2     <-  log( row_Average_one1_week0_EEDheto[ row_Average_one1[[2]] ] / row_Average_one1_week4_EEDheto[ row_Average_one1[[2]] ] ) / 4
one1_EEDko_NTR_BB2       <-  log( row_Average_one1_week0_EEDko[ row_Average_one1[[2]] ] / row_Average_one1_week4_EEDko[ row_Average_one1[[2]] ] ) / 4
two2_EEDheto_NTR_BB2     <-  log( row_Average_two2_week0_EEDheto[ row_Average_two2[[2]] ] / row_Average_two2_week4_EEDheto[ row_Average_two2[[2]] ] ) / 4
two2_EEDko_NTR_BB2       <-  log( row_Average_two2_week0_EEDko[ row_Average_two2[[2]] ] / row_Average_two2_week4_EEDko[ row_Average_two2[[2]] ] ) / 4
three3_EEDheto_NTR_BB2   <-  log( row_Average_three3_week0_EEDheto[ row_Average_three3[[2]] ] / row_Average_three3_week4_EEDheto[ row_Average_three3[[2]] ] ) / 4
three3_EEDko_NTR_BB2     <-  log( row_Average_three3_week0_EEDko[ row_Average_three3[[2]] ] / row_Average_three3_week4_EEDko[ row_Average_three3[[2]] ] ) / 4

one1_EEDheto_halflife_BB2      <- log(2)/one1_EEDheto_NTR_BB2
one1_EEDko_halflifea_BB2       <- log(2)/one1_EEDko_NTR_BB2
two2_EEDheto_halflife_BB2      <- log(2)/two2_EEDheto_NTR_BB2
two2_EEDko_halflife_BB2        <- log(2)/two2_EEDko_NTR_BB2
three3_EEDheto_halflife_BB2    <- log(2)/three3_EEDheto_NTR_BB2
three3_EEDko_halflife_BB2      <- log(2)/three3_EEDko_NTR_BB2

NTR_min_BB2 <- min(one1_EEDheto_NTR_BB2, one1_EEDko_NTR_BB2, two2_EEDheto_NTR_BB2, two2_EEDko_NTR_BB2, three3_EEDheto_NTR_BB2, three3_EEDko_NTR_BB2)
NTR_max_BB2 <- max(one1_EEDheto_NTR_BB2, one1_EEDko_NTR_BB2, two2_EEDheto_NTR_BB2, two2_EEDko_NTR_BB2, three3_EEDheto_NTR_BB2, three3_EEDko_NTR_BB2)
halflife_min_BB2 <- min(one1_EEDheto_halflife_BB2, one1_EEDko_halflifea_BB2, two2_EEDheto_halflife_BB2, two2_EEDko_halflife_BB2, three3_EEDheto_halflife_BB2, three3_EEDko_halflife_BB2)
halflife_max_BB2 <- max(one1_EEDheto_halflife_BB2, one1_EEDko_halflifea_BB2, two2_EEDheto_halflife_BB2, two2_EEDko_halflife_BB2, three3_EEDheto_halflife_BB2, three3_EEDko_halflife_BB2)
NTR_min_BB2
NTR_max_BB2
halflife_min_BB2
halflife_max_BB2











one1_EEDheto_NTR_CC2     <-  log( row_Average_one1_week0_EEDheto[ row_Average_one1[[3]] ] / row_Average_one1_week4_EEDheto[ row_Average_one1[[3]] ] ) / 4
one1_EEDko_NTR_CC2       <-  log( row_Average_one1_week0_EEDko[ row_Average_one1[[3]] ] / row_Average_one1_week4_EEDko[ row_Average_one1[[3]] ] ) / 4
two2_EEDheto_NTR_CC2     <-  log( row_Average_two2_week0_EEDheto[ row_Average_two2[[3]] ] / row_Average_two2_week4_EEDheto[ row_Average_two2[[3]] ] ) / 4
two2_EEDko_NTR_CC2       <-  log( row_Average_two2_week0_EEDko[ row_Average_two2[[3]] ] / row_Average_two2_week4_EEDko[ row_Average_two2[[3]] ] ) / 4
three3_EEDheto_NTR_CC2   <-  log( row_Average_three3_week0_EEDheto[ row_Average_three3[[3]] ] / row_Average_three3_week4_EEDheto[ row_Average_three3[[3]] ] ) / 4
three3_EEDko_NTR_CC2     <-  log( row_Average_three3_week0_EEDko[ row_Average_three3[[3]] ] / row_Average_three3_week4_EEDko[ row_Average_three3[[3]] ] ) / 4

one1_EEDheto_halflife_CC2      <- log(2)/one1_EEDheto_NTR_CC2
one1_EEDko_halflifea_CC2       <- log(2)/one1_EEDko_NTR_CC2
two2_EEDheto_halflife_CC2      <- log(2)/two2_EEDheto_NTR_CC2
two2_EEDko_halflife_CC2        <- log(2)/two2_EEDko_NTR_CC2
three3_EEDheto_halflife_CC2    <- log(2)/three3_EEDheto_NTR_CC2
three3_EEDko_halflife_CC2      <- log(2)/three3_EEDko_NTR_CC2

NTR_min_CC2 <- min(one1_EEDheto_NTR_CC2, one1_EEDko_NTR_CC2, two2_EEDheto_NTR_CC2, two2_EEDko_NTR_CC2, three3_EEDheto_NTR_CC2, three3_EEDko_NTR_CC2)
NTR_max_CC2 <- max(one1_EEDheto_NTR_CC2, one1_EEDko_NTR_CC2, two2_EEDheto_NTR_CC2, two2_EEDko_NTR_CC2, three3_EEDheto_NTR_CC2, three3_EEDko_NTR_CC2)
halflife_min_CC2 <- min(one1_EEDheto_halflife_CC2, one1_EEDko_halflifea_CC2, two2_EEDheto_halflife_CC2, two2_EEDko_halflife_CC2, three3_EEDheto_halflife_CC2, three3_EEDko_halflife_CC2)
halflife_max_CC2 <- max(one1_EEDheto_halflife_CC2, one1_EEDko_halflifea_CC2, two2_EEDheto_halflife_CC2, two2_EEDko_halflife_CC2, three3_EEDheto_halflife_CC2, three3_EEDko_halflife_CC2)
NTR_min_CC2
NTR_max_CC2
halflife_min_CC2
halflife_max_CC2











one1_EEDheto_NTR_DD2     <-  log( row_Average_one1_week0_EEDheto[ row_Average_one1[[4]] ] / row_Average_one1_week4_EEDheto[ row_Average_one1[[4]] ] ) / 4
one1_EEDko_NTR_DD2       <-  log( row_Average_one1_week0_EEDko[ row_Average_one1[[4]] ] / row_Average_one1_week4_EEDko[ row_Average_one1[[4]] ] ) / 4
two2_EEDheto_NTR_DD2     <-  log( row_Average_two2_week0_EEDheto[ row_Average_two2[[4]] ] / row_Average_two2_week4_EEDheto[ row_Average_two2[[4]] ] ) / 4
two2_EEDko_NTR_DD2       <-  log( row_Average_two2_week0_EEDko[ row_Average_two2[[4]] ] / row_Average_two2_week4_EEDko[ row_Average_two2[[4]] ] ) / 4
three3_EEDheto_NTR_DD2   <-  log( row_Average_three3_week0_EEDheto[ row_Average_three3[[4]] ] / row_Average_three3_week4_EEDheto[ row_Average_three3[[4]] ] ) / 4
three3_EEDko_NTR_DD2     <-  log( row_Average_three3_week0_EEDko[ row_Average_three3[[4]] ] / row_Average_three3_week4_EEDko[ row_Average_three3[[4]] ] ) / 4

one1_EEDheto_halflife_DD2      <- log(2)/one1_EEDheto_NTR_DD2
one1_EEDko_halflifea_DD2       <- log(2)/one1_EEDko_NTR_DD2
two2_EEDheto_halflife_DD2      <- log(2)/two2_EEDheto_NTR_DD2
two2_EEDko_halflife_DD2        <- log(2)/two2_EEDko_NTR_DD2
three3_EEDheto_halflife_DD2    <- log(2)/three3_EEDheto_NTR_DD2
three3_EEDko_halflife_DD2      <- log(2)/three3_EEDko_NTR_DD2

NTR_min_DD2 <- min(one1_EEDheto_NTR_DD2, one1_EEDko_NTR_DD2, two2_EEDheto_NTR_DD2, two2_EEDko_NTR_DD2, three3_EEDheto_NTR_DD2, three3_EEDko_NTR_DD2)
NTR_max_DD2 <- max(one1_EEDheto_NTR_DD2, one1_EEDko_NTR_DD2, two2_EEDheto_NTR_DD2, two2_EEDko_NTR_DD2, three3_EEDheto_NTR_DD2, three3_EEDko_NTR_DD2)
halflife_min_DD2 <- min(one1_EEDheto_halflife_DD2, one1_EEDko_halflifea_DD2, two2_EEDheto_halflife_DD2, two2_EEDko_halflife_DD2, three3_EEDheto_halflife_DD2, three3_EEDko_halflife_DD2)
halflife_max_DD2 <- max(one1_EEDheto_halflife_DD2, one1_EEDko_halflifea_DD2, two2_EEDheto_halflife_DD2, two2_EEDko_halflife_DD2, three3_EEDheto_halflife_DD2, three3_EEDko_halflife_DD2)
NTR_min_DD2
NTR_max_DD2
halflife_min_DD2
halflife_max_DD2










one1_EEDheto_NTR_EE2     <-  log( row_Average_one1_week0_EEDheto[ row_Average_one1[[5]] ] / row_Average_one1_week4_EEDheto[ row_Average_one1[[5]] ] ) / 4
one1_EEDko_NTR_EE2       <-  log( row_Average_one1_week0_EEDko[ row_Average_one1[[5]] ] / row_Average_one1_week4_EEDko[ row_Average_one1[[5]] ] ) / 4
two2_EEDheto_NTR_EE2     <-  log( row_Average_two2_week0_EEDheto[ row_Average_two2[[5]] ] / row_Average_two2_week4_EEDheto[ row_Average_two2[[5]] ] ) / 4
two2_EEDko_NTR_EE2       <-  log( row_Average_two2_week0_EEDko[ row_Average_two2[[5]] ] / row_Average_two2_week4_EEDko[ row_Average_two2[[5]] ] ) / 4
three3_EEDheto_NTR_EE2   <-  log( row_Average_three3_week0_EEDheto[ row_Average_three3[[5]] ] / row_Average_three3_week4_EEDheto[ row_Average_three3[[5]] ] ) / 4
three3_EEDko_NTR_EE2     <-  log( row_Average_three3_week0_EEDko[ row_Average_three3[[5]] ] / row_Average_three3_week4_EEDko[ row_Average_three3[[5]] ] ) / 4

one1_EEDheto_halflife_EE2      <- log(2)/one1_EEDheto_NTR_EE2
one1_EEDko_halflifea_EE2       <- log(2)/one1_EEDko_NTR_EE2
two2_EEDheto_halflife_EE2      <- log(2)/two2_EEDheto_NTR_EE2
two2_EEDko_halflife_EE2        <- log(2)/two2_EEDko_NTR_EE2
three3_EEDheto_halflife_EE2    <- log(2)/three3_EEDheto_NTR_EE2
three3_EEDko_halflife_EE2      <- log(2)/three3_EEDko_NTR_EE2

NTR_min_EE2 <- min(one1_EEDheto_NTR_EE2, one1_EEDko_NTR_EE2, two2_EEDheto_NTR_EE2, two2_EEDko_NTR_EE2, three3_EEDheto_NTR_EE2, three3_EEDko_NTR_EE2)
NTR_max_EE2 <- max(one1_EEDheto_NTR_EE2, one1_EEDko_NTR_EE2, two2_EEDheto_NTR_EE2, two2_EEDko_NTR_EE2, three3_EEDheto_NTR_EE2, three3_EEDko_NTR_EE2)
halflife_min_EE2 <- min(one1_EEDheto_halflife_EE2, one1_EEDko_halflifea_EE2, two2_EEDheto_halflife_EE2, two2_EEDko_halflife_EE2, three3_EEDheto_halflife_EE2, three3_EEDko_halflife_EE2)
halflife_max_EE2 <- max(one1_EEDheto_halflife_EE2, one1_EEDko_halflifea_EE2, two2_EEDheto_halflife_EE2, two2_EEDko_halflife_EE2, three3_EEDheto_halflife_EE2, three3_EEDko_halflife_EE2)
NTR_min_EE2
NTR_max_EE2
halflife_min_EE2
halflife_max_EE2








MyBoxViolinPlot_3_s6(vector2=c(one1_EEDheto_NTR_AA2, one1_EEDko_NTR_AA2, two2_EEDheto_NTR_AA2, two2_EEDko_NTR_AA2, three3_EEDheto_NTR_AA2, three3_EEDko_NTR_AA2 ),   
                     sampleType2=c( rep("Down_EEDheto", row_Average_one1_len1),         rep("Down_EEDko", row_Average_one1_len1), 
                                    rep("Unchange_EEDheto", row_Average_two2_len1),     rep("Unchange_EEDko", row_Average_two2_len1),
                                    rep("Up_EEDheto", row_Average_three3_len1),         rep("Up_EEDko", row_Average_three3_len1)  ), 
                     sampleRank2=c( "Down_EEDheto",   "Down_EEDko",  "Unchange_EEDheto",  "Unchange_EEDko",   "Up_EEDheto",  "Up_EEDko"),  
                     colours2=c("red",  "red4",   "blue",    "blue4",  "green",  "green4" ),   
                     path2=subdir_3_part4,  fileName2="1-one1-BoxViolin-NTR",  
                     title2="Lowest",  xLab2="Samples",  yLab2="NTR",   
                     height2=3.88,   width2=5, Ymin2=0, Ymax2=0.8 )





MyBoxViolinPlot_3_s6(vector2=c(one1_EEDheto_NTR_BB2, one1_EEDko_NTR_BB2, two2_EEDheto_NTR_BB2, two2_EEDko_NTR_BB2, three3_EEDheto_NTR_BB2, three3_EEDko_NTR_BB2 ),   
                     sampleType2=c( rep("Down_EEDheto", row_Average_one1_len2),         rep("Down_EEDko", row_Average_one1_len2), 
                                    rep("Unchange_EEDheto", row_Average_two2_len2),     rep("Unchange_EEDko", row_Average_two2_len2),
                                    rep("Up_EEDheto", row_Average_three3_len2),         rep("Up_EEDko", row_Average_three3_len2)  ), 
                     sampleRank2=c( "Down_EEDheto",   "Down_EEDko",  "Unchange_EEDheto",  "Unchange_EEDko",   "Up_EEDheto",  "Up_EEDko"),  
                     colours2=c("red",  "red4",   "blue",    "blue4",  "green",  "green4" ),   
                     path2=subdir_3_part4,  fileName2="2-two2-BoxViolin-NTR",  
                     title2="Low",  xLab2="Samples",  yLab2="NTR",   
                     height2=3.88,   width2=5, Ymin2=0, Ymax2=0.8 )



MyBoxViolinPlot_3_s6(vector2=c(one1_EEDheto_NTR_CC2, one1_EEDko_NTR_CC2, two2_EEDheto_NTR_CC2, two2_EEDko_NTR_CC2, three3_EEDheto_NTR_CC2, three3_EEDko_NTR_CC2 ),   
                     sampleType2=c( rep("Down_EEDheto", row_Average_one1_len3),         rep("Down_EEDko", row_Average_one1_len3), 
                                    rep("Unchange_EEDheto", row_Average_two2_len3),     rep("Unchange_EEDko", row_Average_two2_len3),
                                    rep("Up_EEDheto", row_Average_three3_len3),         rep("Up_EEDko", row_Average_three3_len3)  ), 
                     sampleRank2=c( "Down_EEDheto",   "Down_EEDko",  "Unchange_EEDheto",  "Unchange_EEDko",   "Up_EEDheto",  "Up_EEDko"),  
                     colours2=c("red",  "red4",   "blue",    "blue4",  "green",  "green4" ),   
                     path2=subdir_3_part4,  fileName2="3-BoxViolin-NTR",  
                     title2="Medium",  xLab2="Samples",  yLab2="NTR",   
                     height2=3.88,   width2=5, Ymin2=0, Ymax2=0.8 )



MyBoxViolinPlot_3_s6(vector2=c(one1_EEDheto_NTR_DD2, one1_EEDko_NTR_DD2, two2_EEDheto_NTR_DD2, two2_EEDko_NTR_DD2, three3_EEDheto_NTR_DD2, three3_EEDko_NTR_DD2 ),   
                     sampleType2=c( rep("Down_EEDheto", row_Average_one1_len4),         rep("Down_EEDko", row_Average_one1_len4), 
                                    rep("Unchange_EEDheto", row_Average_two2_len4),     rep("Unchange_EEDko", row_Average_two2_len4),
                                    rep("Up_EEDheto", row_Average_three3_len4),         rep("Up_EEDko", row_Average_three3_len4)  ), 
                     sampleRank2=c( "Down_EEDheto",   "Down_EEDko",  "Unchange_EEDheto",  "Unchange_EEDko",   "Up_EEDheto",  "Up_EEDko"),  
                     colours2=c("red",  "red4",   "blue",    "blue4",  "green",  "green4" ),   
                     path2=subdir_3_part4,  fileName2="4-BoxViolin-NTR",  
                     title2="High",  xLab2="Samples",  yLab2="NTR",   
                     height2=3.88,   width2=5, Ymin2=0, Ymax2=0.8 )




MyBoxViolinPlot_3_s6(vector2=c(one1_EEDheto_NTR_EE2, one1_EEDko_NTR_EE2, two2_EEDheto_NTR_EE2, two2_EEDko_NTR_EE2, three3_EEDheto_NTR_EE2, three3_EEDko_NTR_EE2 ),   
                     sampleType2=c( rep("Down_EEDheto", row_Average_one1_len5),         rep("Down_EEDko", row_Average_one1_len5), 
                                    rep("Unchange_EEDheto", row_Average_two2_len5),     rep("Unchange_EEDko", row_Average_two2_len5),
                                    rep("Up_EEDheto", row_Average_three3_len5),         rep("Up_EEDko", row_Average_three3_len5)  ), 
                     sampleRank2=c( "Down_EEDheto",   "Down_EEDko",  "Unchange_EEDheto",  "Unchange_EEDko",   "Up_EEDheto",  "Up_EEDko"),  
                     colours2=c("red",  "red4",   "blue",    "blue4",  "green",  "green4" ),   
                     path2=subdir_3_part4,  fileName2="5-BoxViolin-NTR",  
                     title2="Highest",  xLab2="Samples",  yLab2="NTR",   
                     height2=3.88,   width2=5, Ymin2=0, Ymax2=0.8 )












##########################################################################
subdir_4_part4 <- paste(Part4_g,  "/4-NTR-and-NOL", sep = "")
if( ! file.exists(subdir_4_part4) ) { dir.create(subdir_4_part4) }


row_Average_one1_EEDheto_NTR2       <-  log( row_Average_one1_week0_EEDheto/row_Average_one1_week4_EEDheto)/4
row_Average_one1_EEDko_NTR2         <-  log( row_Average_one1_week0_EEDko/row_Average_one1_week4_EEDko )/4
row_Average_two2_EEDheto_NTR2       <-  log( row_Average_two2_week0_EEDheto/row_Average_two2_week4_EEDheto )/4
row_Average_two2_EEDko_NTR2         <-  log( row_Average_two2_week0_EEDko/row_Average_two2_week4_EEDko )/4
row_Average_three3_EEDheto_NTR2     <-  log( row_Average_three3_week0_EEDheto/row_Average_three3_week4_EEDheto )/4
row_Average_three3_EEDko_NTR2       <-  log( row_Average_three3_week0_EEDko/row_Average_three3_week4_EEDko )/4



MyScatterDiagram_2(vector2X=row_one1_H3_normal_rep1,  
                        vector2Y=row_Average_one1_EEDko_NTR2,  
                        path2=subdir_4_part4,   fileName2="1-one1-EEDko",   
                        xLab2="NOL",   yLab2="NTR",  title2="Down-regulated Genes",   
                        height2=4,  width2=4,  yMin2=0, yMax2=0.6,  xMin2=0, xMax2=2.5, alpha2=1)
        
MyScatterDiagram_2(vector2X=row_one1_H3_normal_rep1,  
                   vector2Y=row_Average_one1_EEDheto_NTR2,  
                   path2=subdir_4_part4,   fileName2="1-one1-EEDheto",   
                   xLab2="NOL",   yLab2="NTR",  title2="Down-regulated Genes",   
                   height2=4,  width2=4,  yMin2=0, yMax2=0.6,  xMin2=0, xMax2=2.5, alpha2=1)

MyCor2Vars_1( vector1=row_one1_H3_normal_rep1, vector2=row_Average_one1_EEDko_NTR2,   file1=paste(subdir_4_part4,  "/1-one1-EEDko.txt",        sep = "") )        
MyCor2Vars_1( vector1=row_one1_H3_normal_rep1, vector2=row_Average_one1_EEDheto_NTR2, file1=paste(subdir_4_part4,  "/1-one1-EEDheto.txt",      sep = "") )        





MyScatterDiagram_2(vector2X=row_two2_H3_normal_rep1,  
                   vector2Y=row_Average_two2_EEDko_NTR2,  
                   path2=subdir_4_part4,   fileName2="2-two2-EEDko",   
                   xLab2="NOL",   yLab2="NTR",  title2="Unchanged Genes",   
                   height2=4,  width2=4,  yMin2=0, yMax2=0.6,  xMin2=0, xMax2=2.5,  alpha2=0.5)

MyScatterDiagram_2(vector2X=row_two2_H3_normal_rep1,  
                   vector2Y=row_Average_two2_EEDheto_NTR2,  
                   path2=subdir_4_part4,   fileName2="2-two2-EEDheto",   
                   xLab2="NOL",   yLab2="NTR",  title2="Unchanged Genes",   
                   height2=4,  width2=4,  yMin2=0, yMax2=0.6,  xMin2=0, xMax2=2.5, alpha2=0.5)

MyCor2Vars_1( vector1=row_two2_H3_normal_rep1, vector2=row_Average_two2_EEDko_NTR2,   file1=paste(subdir_4_part4,  "/2-two2-EEDko.txt",        sep = "") )        
MyCor2Vars_1( vector1=row_two2_H3_normal_rep1, vector2=row_Average_two2_EEDheto_NTR2, file1=paste(subdir_4_part4,  "/2-two2-EEDheto.txt",      sep = "") )        







MyScatterDiagram_2(vector2X=row_three3_H3_normal_rep1,  
                   vector2Y=row_Average_three3_EEDko_NTR2,  
                   path2=subdir_4_part4,   fileName2="3-three3-EEDko",   
                   xLab2="NOL",   yLab2="NTR",  title2="Up-regulated Genes",   
                   height2=4,  width2=4,  yMin2=0, yMax2=0.6,  xMin2=0, xMax2=2.5, alpha2=1)

MyScatterDiagram_2(vector2X=row_three3_H3_normal_rep1,  
                   vector2Y=row_Average_three3_EEDheto_NTR2,  
                   path2=subdir_4_part4,   fileName2="3-three3-EEDheto",   
                   xLab2="NOL",   yLab2="NTR",  title2="Up-regulated Genes",   
                   height2=4,  width2=4,  yMin2=0, yMax2=0.6,  xMin2=0, xMax2=2.5, alpha2=1)

MyCor2Vars_1( vector1=row_three3_H3_normal_rep1, vector2=row_Average_three3_EEDko_NTR2,   file1=paste(subdir_4_part4,  "/3-three3-EEDko.txt",        sep = "") )        
MyCor2Vars_1( vector1=row_three3_H3_normal_rep1, vector2=row_Average_three3_EEDheto_NTR2, file1=paste(subdir_4_part4,  "/3-three3-EEDheto.txt",      sep = "") )        



####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################




