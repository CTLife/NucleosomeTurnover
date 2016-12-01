#############################################################################################################################
## There are ten parts for analyzing the DANPOS2.2.2 processed results of nucleosome turnover ChIP-seq data:
## Part  1:  Load R libraries and define some functions.
## Part  2:  Read in the input matrix data, and reduce the dimensions (rows or columns) of matrix, and save the information of input files.
## Part  3:  Figures about nucleosome occupancy level (NOL).
## Part  4:  Figures about nucleosome turnover rate (NTR) and compute the correlation between NOL and NTR.
## Part  5:  Classify DNA regions based on NOL.
## Part  6:  Classify DNA regions based on NTR.
## Part  7:  Figures about nucleosome contribution value (NCV).  
##                   3D-graph: x=NOL, y=NTR, z=biological function.            (based on Data)
##                   NCV=f(NOL, NTR)=b1/[1+(NOL/k1)^n1] + b2/[1+(k2/NTR)^n2]   (based on Model)
## Part  8:  Figures about a smaller matrix, so we can compute NTR at a single DNA region. Finally, a NTR matrix should be generated.  
## Part  9:  Correaltions between the 3 indexes (NOL, NTR, NCV) and biological functions (GO, KEGG, DEGs, DPs, etc).    
## Part 10:  Other analysis.
#############################################################################################################################


#############################################################################################################################
## Part 3:  Figures about nucleosome occupancy level (NOL).
#############################################################################################################################





#################################################################### Start ##########################################################################################################################################
subdir_1_part3 <- paste(Part3_g,  "/1-AverageColumns-AllRows", sep = "")
if( ! file.exists(subdir_1_part3) ) { dir.create(subdir_1_part3) }

MyAverageLines_1(vector2=c(column_Average_week0, column_Average_week1,  column_Average_week2,  column_Average_week4,  column_Average_week6, column_Average_week8  ),   
             numSample2=6,   
             sampleType2=c( rep("week 0", numOfColumns1),    rep("week 1", numOfColumns1),  rep("week 2", numOfColumns1),  
                            rep("week 4", numOfColumns1),    rep("week 6", numOfColumns1),  rep("week 8",   numOfColumns1)  ), 
             sampleRank2=c( "week 0",  "week 1",   "week 2",  "week 4",  "week 6",   "week 8" ),     
             colours2=c( "week 0"="red2",   "week 1"="cyan2",    "week 2"="green2",  "week 4"="blue2",  "week 6"="purple2",   "week 8"="orange2" ), 
             path2=subdir_1_part3,     fileName2=paste("Part3-1-A-WT-averageCurve",  "5kb", mySample_g, sep = "_") ,  
             title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
             Ymin2=0,   Ymax2=1.25,    height2=3.25,   width2=4.83, center2=myCenter_g  )    ## height=4cm, width=5cm

numOfColumns2 <- 100
MyAverageLines_3(vector2=c(column_Average_week0[200:299], column_Average_week1[200:299],  column_Average_week2[200:299],  column_Average_week4[200:299],  column_Average_week6[200:299], column_Average_week8[200:299]  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week 0", numOfColumns2),    rep("week 1", numOfColumns2),  rep("week 2",   numOfColumns2),  
                                rep("week 4", numOfColumns2),    rep("week 6", numOfColumns2),  rep("week 8",   numOfColumns2)  ), 
                 sampleRank2=c( "week 0",  "week 1",   "week 2",  "week 4",  "week 6",   "week 8" ),     
                 colours2=c( "week 0"="red2",   "week 1"="cyan2",    "week 2"="green2",  "week 4"="blue2",  "week 6"="purple2",   "week 8"="orange2" ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-A-WT-averageCurve",  "1kb", mySample_g, sep = "_") ,  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.25,   width2=4.83, center2=myCenter_g  )    ## 4cm + 5cm

numOfColumns3 <- 200
MyAverageLines_9(vector2=c(column_Average_week0[150:349], column_Average_week1[150:349],  column_Average_week2[150:349],  column_Average_week4[150:349],  column_Average_week6[150:349], column_Average_week8[150:349]  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week 0", numOfColumns3),    rep("week 1", numOfColumns3),  rep("week 2",   numOfColumns3),  
                                rep("week 4", numOfColumns3),    rep("week 6", numOfColumns3),  rep("week 8",   numOfColumns3)  ), 
                 sampleRank2=c( "week 0",  "week 1",   "week 2",  "week 4",  "week 6",   "week 8" ),     
                 colours2=c( "week 0"="red2",   "week 1"="cyan2",    "week 2"="green2",  "week 4"="blue2",  "week 6"="purple2",   "week 8"="orange2" ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-A-WT-averageCurve",  "2kb", mySample_g, sep = "_") ,  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.25,   width2=4.83, center2=myCenter_g  )    ## 4cm + 5cm

numOfColumns4 <- 300
MyAverageLines_13(vector2=c(column_Average_week0[100:399], column_Average_week1[100:399],  column_Average_week2[100:399],  column_Average_week4[100:399],  column_Average_week6[100:399], column_Average_week8[100:399]  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week 0", numOfColumns4),    rep("week 1", numOfColumns4),  rep("week 2",   numOfColumns4),  
                                rep("week 4", numOfColumns4),    rep("week 6", numOfColumns4),  rep("week 8",   numOfColumns4)  ), 
                 sampleRank2=c( "week 0",  "week 1",   "week 2",  "week 4",  "week 6",   "week 8" ),     
                 colours2=c( "week 0"="red2",   "week 1"="cyan2",    "week 2"="green2",  "week 4"="blue2",  "week 6"="purple2",   "week 8"="orange2" ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-A-WT-averageCurve",  "3kb", mySample_g, sep = "_") ,  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.25,   width2=4.83, center2=myCenter_g  )    ## 4cm + 5cm





MyAverageLines_1(vector2=c(column_Average_week0, column_Average_H3  ),   
                 numSample2=2,   
                 sampleType2=c( rep("week 0", numOfColumns1),    rep("H3", numOfColumns1)   ), 
                 sampleRank2=c( "week 0",  "H3"  ),     
                 colours2=c( "week 0"="red2",   "H3"="yellow2"  ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-B-WT-averageCurve-H3-week0",  "5kb", mySample_g, sep = "_") ,  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,     height2=3.25,   width2=4.83, center2=myCenter_g   )

MyAverageLines_3(vector2=c(column_Average_week0[200:299], column_Average_H3[200:299]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("week 0", numOfColumns2),    rep("H3", numOfColumns2)   ), 
                 sampleRank2=c( "week 0",  "H3"  ),     
                 colours2=c( "week 0"="red2",   "H3"="yellow2"  ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-B-WT-averageCurve-H3-week0",  "1kb", mySample_g, sep = "_") ,  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,     height2=3.25,   width2=4.83 , center2=myCenter_g  )

MyAverageLines_9(vector2=c(column_Average_week0[150:349], column_Average_H3[150:349]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("week 0", numOfColumns3),    rep("H3", numOfColumns3)   ), 
                 sampleRank2=c( "week 0",  "H3"  ),     
                 colours2=c( "week 0"="red2",   "H3"="yellow2"  ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-B-WT-averageCurve-H3-week0",  "2kb", mySample_g, sep = "_") ,  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,     height2=3.25,   width2=4.83 , center2=myCenter_g  )

MyAverageLines_13(vector2=c(column_Average_week0[100:399], column_Average_H3[100:399]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("week 0", numOfColumns4),    rep("H3", numOfColumns4)   ), 
                 sampleRank2=c( "week 0",  "H3"  ),     
                 colours2=c( "week 0"="red2",   "H3"="yellow2"  ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-B-WT-averageCurve-H3-week0",  "3kb", mySample_g, sep = "_") ,  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,     height2=3.25,   width2=4.83 , center2=myCenter_g  )





MyAverageLines_1(vector2=c(column_Average_H3,  column_Average_week0, column_Average_week1,  column_Average_week2,  column_Average_week4,  column_Average_week6, column_Average_week8  ),   
                 numSample2=7,   
                 sampleType2=c( rep("H3", numOfColumns1), rep("week 0", numOfColumns1),    rep("week 1", numOfColumns1),  rep("week 2", numOfColumns1),  
                                rep("week 4", numOfColumns1),    rep("week 6", numOfColumns1),  rep("week 8",   numOfColumns1)  ), 
                 sampleRank2=c( "H3", "week 0",  "week 1",   "week 2",  "week 4",  "week 6",   "week 8" ),     
                 colours2=c( "H3"="yellow2",  "week 0"="red2",   "week 1"="cyan2",    "week 2"="green2",  "week 4"="blue2",  "week 6"="purple2",   "week 8"="orange2" ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-C-WT-averageCurve-withH3",  "5kb", mySample_g, sep = "_") ,  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,     height2=3.25,   width2=4.83  , center2=myCenter_g  )

MyAverageLines_3(vector2=c(column_Average_H3[200:299],  column_Average_week0[200:299], column_Average_week1[200:299],  column_Average_week2[200:299],  column_Average_week4[200:299],  column_Average_week6[200:299], column_Average_week8[200:299]  ),   
                 numSample2=7,   
                 sampleType2=c( rep("H3", numOfColumns2), rep("week 0", numOfColumns2),    rep("week 1", numOfColumns2),  rep("week 2", numOfColumns2),  
                                rep("week 4", numOfColumns2),    rep("week 6", numOfColumns2),  rep("week 8",   numOfColumns2)  ), 
                 sampleRank2=c( "H3", "week 0",  "week 1",   "week 2",  "week 4",  "week 6",   "week 8" ),     
                 colours2=c( "H3"="yellow2",  "week 0"="red2",   "week 1"="cyan2",    "week 2"="green2",  "week 4"="blue2",  "week 6"="purple2",   "week 8"="orange2" ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-C-WT-averageCurve-withH3",  "1kb", mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,     height2=3.25,   width2=4.83 , center2=myCenter_g   )

MyAverageLines_9(vector2=c(column_Average_H3[150:349],  column_Average_week0[150:349], column_Average_week1[150:349],  column_Average_week2[150:349],  column_Average_week4[150:349],  column_Average_week6[150:349], column_Average_week8[150:349]  ),   
                 numSample2=7,   
                 sampleType2=c( rep("H3", numOfColumns3), rep("week 0", numOfColumns3),    rep("week 1", numOfColumns3),  rep("week 2", numOfColumns3),  
                                rep("week 4", numOfColumns3),    rep("week 6", numOfColumns3),  rep("week 8",   numOfColumns3)  ), 
                 sampleRank2=c( "H3", "week 0",  "week 1",   "week 2",  "week 4",  "week 6",   "week 8" ),     
                 colours2=c( "H3"="yellow2",  "week 0"="red2",   "week 1"="cyan2",    "week 2"="green2",  "week 4"="blue2",  "week 6"="purple2",   "week 8"="orange2" ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-C-WT-averageCurve-withH3",  "2kb", mySample_g, sep = "_") , 
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,     height2=3.25,   width2=4.83  , center2=myCenter_g  )

MyAverageLines_13(vector2=c(column_Average_H3[100:399],  column_Average_week0[100:399], column_Average_week1[100:399],  column_Average_week2[100:399],  column_Average_week4[100:399],  column_Average_week6[100:399], column_Average_week8[100:399]  ),   
                 numSample2=7,   
                 sampleType2=c( rep("H3", numOfColumns4), rep("week 0", numOfColumns4),    rep("week 1", numOfColumns4),  rep("week 2", numOfColumns4),  
                                rep("week 4", numOfColumns4),    rep("week 6", numOfColumns4),  rep("week 8",   numOfColumns4)  ), 
                 sampleRank2=c( "H3", "week 0",  "week 1",   "week 2",  "week 4",  "week 6",   "week 8" ),     
                 colours2=c( "H3"="yellow2",  "week 0"="red2",   "week 1"="cyan2",    "week 2"="green2",  "week 4"="blue2",  "week 6"="purple2",   "week 8"="orange2" ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-C-WT-averageCurve-withH3",  "3kb", mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,     height2=3.25,   width2=4.83  , center2=myCenter_g  )





MyAverageLines_1(vector2=c(column_Average_week0_EEDheto,  column_Average_week0_EEDko,  column_Average_week4_EEDheto,  column_Average_week4_EEDko  ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-D-EEDko-averageCurve",  "5kb", mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,     height2=3.25,   width2=5.6, center2=myCenter_g   )

MyAverageLines_3(vector2=c(column_Average_week0_EEDheto[200:299],  column_Average_week0_EEDko[200:299],  column_Average_week4_EEDheto[200:299],  column_Average_week4_EEDko[200:299]  ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns2),    rep("week0_EEDko", numOfColumns2),   
                                rep("week4_EEDheto", numOfColumns2),    rep("week4_EEDko", numOfColumns2)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-D-EEDko-averageCurve",  "1kb", mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,     height2=3.25,   width2=5.6 , center2=myCenter_g  )

MyAverageLines_9(vector2=c(column_Average_week0_EEDheto[150:349],  column_Average_week0_EEDko[150:349],  column_Average_week4_EEDheto[150:349],  column_Average_week4_EEDko[150:349]  ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns3),    rep("week0_EEDko", numOfColumns3),   
                                rep("week4_EEDheto", numOfColumns3),    rep("week4_EEDko", numOfColumns3)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_1_part3,      fileName2=paste("Part3-1-D-EEDko-averageCurve",  "2kb", mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,     height2=3.25,   width2=5.6 , center2=myCenter_g  )

MyAverageLines_13(vector2=c(column_Average_week0_EEDheto[100:399],  column_Average_week0_EEDko[100:399],  column_Average_week4_EEDheto[100:399],  column_Average_week4_EEDko[100:399]  ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns4),    rep("week0_EEDko", numOfColumns4),   
                                rep("week4_EEDheto", numOfColumns4),    rep("week4_EEDko", numOfColumns4)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_1_part3,      fileName2=paste("Part3-1-D-EEDko-averageCurve",  "3kb", mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,     height2=3.25,   width2=5.6 , center2=myCenter_g  )





MyAverageLines_1(vector2=c(column_Average_banding,  column_Average_sham  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-E-TAC-averageCurve",  "5kb", mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,     height2=3.25,   width2=4.9 , center2=myCenter_g  )

MyAverageLines_3(vector2=c(column_Average_banding[200:299],  column_Average_sham[200:299]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns2),    rep("sham", numOfColumns2)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-E-TAC-averageCurve",  "1kb", mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,     height2=3.25,   width2=4.9  , center2=myCenter_g  )

MyAverageLines_9(vector2=c(column_Average_banding[150:349],  column_Average_sham[150:349]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns3),    rep("sham", numOfColumns3)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-E-TAC-averageCurve",  "2kb", mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,     height2=3.25,   width2=4.9 , center2=myCenter_g  )

MyAverageLines_13(vector2=c(column_Average_banding[100:399],  column_Average_sham[100:399]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns4),    rep("sham", numOfColumns4)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-E-TAC-averageCurve",  "3kb", mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,     height2=3.25,   width2=4.9  , center2=myCenter_g  )





MyAverageLines_1(vector2=c(column_week0_Rep1,  column_week0_Rep2,  column_week1_Rep1,  column_week1_Rep2,  column_week2_Rep1,  column_week2_Rep2,  
                           column_week4_Rep1,  column_week4_Rep2,  column_week4_Rep3,  column_week6_Rep1,  column_week6_Rep2,  column_week8_Rep1  ),   
                 numSample2=12,   
                 sampleType2=c( rep("week0_Rep1", numOfColumns1),    rep("week0_Rep2", numOfColumns1),    rep("week1_Rep1", numOfColumns1),  rep("week1_Rep2", numOfColumns1),  
                                rep("week2_Rep1", numOfColumns1),    rep("week2_Rep2", numOfColumns1),  
                                rep("week4_Rep1", numOfColumns1),    rep("week4_Rep2", numOfColumns1),    rep("week4_Rep3", numOfColumns1),    
                                rep("week6_Rep1", numOfColumns1),    rep("week6_Rep2", numOfColumns1),    rep("week8_Rep1",   numOfColumns1)  ), 
                 sampleRank2=c( "week0_Rep1",  "week0_Rep2",  "week1_Rep1",   "week1_Rep2",  "week2_Rep1",  "week2_Rep2",  
                                "week4_Rep1",  "week4_Rep2",  "week4_Rep3",   "week6_Rep1",  "week6_Rep2",  "week8_Rep1" ),     
                 colours2=c( "week0_Rep1"="red",       "week0_Rep2"="red4",     "week1_Rep1"="cyan",    "week1_Rep2"="cyan4",    
                             "week2_Rep1"="green",     "week2_Rep2"="green4",  
                             "week4_Rep1"="skyblue",   "week4_Rep2"="blue",    "week4_Rep3"="blue4",  
                             "week6_Rep1"="purple",    "week6_Rep2"="purple4",  "week8_Rep1"="orange2" ), 
                 path2=subdir_1_part3,    fileName2=paste("Part3-1-F-WT-averageCurve-Reps",  "5kb", mySample_g, sep = "_") , 
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,     height2=3.25,   width2=5.3  , center2=myCenter_g  )

MyAverageLines_1(vector2=c(column_week0_Rep1,  column_week0_Rep2,  column_week1_Rep1,  column_week1_Rep2,  column_week2_Rep1,  column_week2_Rep2,  
                           column_week4_Rep2,  column_week4_Rep3,  column_week6_Rep2,  column_week8_Rep1  ),   
                 numSample2=10,   
                 sampleType2=c( rep("week0_Rep1", numOfColumns1),    rep("week0_Rep2", numOfColumns1),    rep("week1_Rep1", numOfColumns1),  rep("week1_Rep2", numOfColumns1),  
                                rep("week2_Rep1", numOfColumns1),    rep("week2_Rep2", numOfColumns1),  
                                rep("week4_Rep2", numOfColumns1),    rep("week4_Rep3", numOfColumns1),    
                                rep("week6_Rep2", numOfColumns1),    rep("week8_Rep1",   numOfColumns1)  ), 
                 sampleRank2=c( "week0_Rep1",  "week0_Rep2",  "week1_Rep1",   "week1_Rep2",  "week2_Rep1",  "week2_Rep2",  
                                "week4_Rep2",  "week4_Rep3",   "week6_Rep2",  "week8_Rep1" ),     
                 colours2=c( "week0_Rep1"="red",       "week0_Rep2"="red4",     "week1_Rep1"="cyan",    "week1_Rep2"="cyan4",    
                             "week2_Rep1"="green",     "week2_Rep2"="green4",  
                             "week4_Rep2"="blue",    "week4_Rep3"="blue4",  
                             "week6_Rep2"="purple4",  "week8_Rep1"="orange2" ), 
                 path2=subdir_1_part3,    fileName2=paste("Part3-1-G-WT-averageCurve-Reps",  "5kb", mySample_g, sep = "_") , 
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,     height2=3.25,   width2=5.3  , center2=myCenter_g  )

MyAverageLines_1(vector2=c(column_week0_EEDheto_Rep1,  column_week0_EEDheto_Rep2,  column_week0_EEDko_Rep1,  
                           column_week4_EEDheto_Rep1,  column_week4_EEDheto_Rep2,  column_week4_EEDko_Rep1,  column_week4_EEDko_Rep2  ),   
                 numSample2=7,   
                 sampleType2=c( rep("week0_EEDheto_Rep1", numOfColumns1),    rep("week0_EEDheto_Rep2", numOfColumns1),    rep("week0_EEDko_Rep1", numOfColumns1),   
                                rep("week4_EEDheto_Rep1", numOfColumns1),    rep("week4_EEDheto_Rep2", numOfColumns1),    rep("week4_EEDko_Rep1", numOfColumns1),  rep("week4_EEDko_Rep2", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto_Rep1",  "week0_EEDheto_Rep2",  "week0_EEDko_Rep1",   
                                "week4_EEDheto_Rep1",  "week4_EEDheto_Rep2",  "week4_EEDko_Rep1",  "week4_EEDko_Rep2"  ),     
                 colours2=c( "week0_EEDheto_Rep1"="red",      "week0_EEDheto_Rep2"="red4",     "week0_EEDko_Rep1"="orange2",   
                             "week4_EEDheto_Rep1"="blue",     "week4_EEDheto_Rep2"="blue4",    "week4_EEDko_Rep1"="green", "week4_EEDko_Rep2"="green4" ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-H-EEDko-averageCurve-Reps",  "5kb", mySample_g, sep = "_") , 
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.25,   width2=6.1  , center2=myCenter_g  )

MyAverageLines_1(vector2=c(column_banding_Rep1,  column_banding_Rep2,  column_sham_Rep1,  column_sham_Rep2  ),   
                 numSample2=4,   
                 sampleType2=c( rep("banding_Rep1", numOfColumns1),  rep("banding_Rep2", numOfColumns1),    rep("sham_Rep1", numOfColumns1),  rep("sham_Rep2", numOfColumns1)  ), 
                 sampleRank2=c( "banding_Rep1",   "banding_Rep2",   "sham_Rep1",  "sham_Rep2"  ),     
                 colours2=c( "banding_Rep1"="blue",  "banding_Rep2"="blue4",  "sham_Rep1"="green",  "sham_Rep2"="green4"  ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-I-TAC-averageCurve-Reps",  "5kb", mySample_g, sep = "_") , 
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.25,   width2=5.4  , center2=myCenter_g  )





MyAverageLines_2(vector2=c(column_Average_H3,  column_Average_week0, column_Average_week1,  column_Average_week2,  column_Average_week4,  column_Average_week6, column_Average_week8  ),  
                 SEM2=c(SEM_Average_H3, SEM_Average_week0, SEM_Average_week1, SEM_Average_week2, SEM_Average_week4, SEM_Average_week6, SEM_Average_week8), 
                 numSample2=7,   
                 sampleType2=c( rep("H3", numOfColumns1), rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8",   numOfColumns1)  ), 
                 sampleRank2=c( "H3", "week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "H3"="yellow2",  "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-J-WT-averageCurve-errorBar",  "5kb", mySample_g, sep = "_") , 
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,     height2=3.25,   width2=4.8  , center2=myCenter_g  )

MyAverageLines_2(vector2=c(column_Average_week0_EEDheto,  column_Average_week0_EEDko,  column_Average_week4_EEDheto,  column_Average_week4_EEDko  ),  
                 SEM2=c(SEM_Average_week0_EEDheto,  SEM_Average_week0_EEDko,  SEM_Average_week4_EEDheto,  SEM_Average_week4_EEDko),
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-K-EEDko-averageCurve-errorBar",  "5kb", mySample_g, sep = "_") , 
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.25,   width2=5.6  , center2=myCenter_g  )

MyAverageLines_2(vector2=c(column_Average_banding,  column_Average_sham  ), 
                 SEM2=c(SEM_Average_banding, SEM_Average_sham),
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_1_part3,     fileName2=paste("Part3-1-L-TAC-averageCurve-errorBar",  "5kb", mySample_g, sep = "_") , 
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.25,   width2=4.9  , center2=myCenter_g  )










## hypothesis test
################################################################################################
subdir_2_part3 <- paste(Part3_g,  "/2-averageColumnsOfAllRows-hypothesisTest", sep = "")
if( ! file.exists(subdir_2_part3) ) { dir.create(subdir_2_part3) }

## Kolmogorov-Smirnov tests
MyHypothesisTest_4(vector1=column_Average_H3,      vector2=column_Average_week0,    file1=paste(subdir_2_part3,  "/Part3-2-A-WT-H3-vs-week0.Kolmogorov-Smirnov.txt",     sep = "")   )   
MyHypothesisTest_4(vector1=column_Average_week0,   vector2=column_Average_week1,    file1=paste(subdir_2_part3,  "/Part3-2-A-WT-week0-vs-week1.Kolmogorov-Smirnov.txt",  sep = "")   )   
MyHypothesisTest_4(vector1=column_Average_week1,   vector2=column_Average_week2,    file1=paste(subdir_2_part3,  "/Part3-2-A-WT-week1-vs-week2.Kolmogorov-Smirnov.txt",  sep = "")   )   
MyHypothesisTest_4(vector1=column_Average_week2,   vector2=column_Average_week4,    file1=paste(subdir_2_part3,  "/Part3-2-A-WT-week2-vs-week4.Kolmogorov-Smirnov.txt",  sep = "")   )   
MyHypothesisTest_4(vector1=column_Average_week4,   vector2=column_Average_week6,    file1=paste(subdir_2_part3,  "/Part3-2-A-WT-week4-vs-week6.Kolmogorov-Smirnov.txt",  sep = "")   )   
MyHypothesisTest_4(vector1=column_Average_week6,   vector2=column_Average_week8,    file1=paste(subdir_2_part3,  "/Part3-2-A-WT-week6-vs-week8.Kolmogorov-Smirnov.txt",  sep = "")   )   
MyHypothesisTest_4(vector1=column_Average_week0_EEDheto,   vector2=column_Average_week0_EEDko,      file1=paste(subdir_2_part3,  "/Part3-2-B-CKO-week0Heto-vs-week0KO.Kolmogorov-Smirnov.txt",  sep = "")   )   
MyHypothesisTest_4(vector1=column_Average_week0_EEDheto,   vector2=column_Average_week4_EEDheto,    file1=paste(subdir_2_part3,  "/Part3-2-B-CKO-week0Heto-vs-week4Heto.Kolmogorov-Smirnov.txt",  sep = "")   )   
MyHypothesisTest_4(vector1=column_Average_week0_EEDko,     vector2=column_Average_week4_EEDko,      file1=paste(subdir_2_part3,  "/Part3-2-B-CKO-week0KO-vs-week4KO.Kolmogorov-Smirnov.txt",  sep = "")   )   
MyHypothesisTest_4(vector1=column_Average_week4_EEDheto,   vector2=column_Average_week4_EEDko,      file1=paste(subdir_2_part3,  "/Part3-2-B-CKO-week4Heto-vs-week4KO.Kolmogorov-Smirnov.txt",  sep = "")   )   
MyHypothesisTest_4(vector1=column_Average_banding,   vector2=column_Average_sham,    file1=paste(subdir_2_part3,  "/Part3-2-C-TAC-banding-vs-sham.Kolmogorov-Smirnov.txt",  sep = "")   )   

## T test and Wilcoxon test  (paired)
MyHypothesisTest_2(vector1=column_Average_H3,      vector2=column_Average_week0,    file1=paste(subdir_2_part3,  "/Part3-2-D-WT-H3-vs-week0.paired.txt",     sep = "")   )   
MyHypothesisTest_2(vector1=column_Average_week0,   vector2=column_Average_week1,    file1=paste(subdir_2_part3,  "/Part3-2-D-WT-week0-vs-week1.paired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=column_Average_week1,   vector2=column_Average_week2,    file1=paste(subdir_2_part3,  "/Part3-2-D-WT-week1-vs-week2.paired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=column_Average_week2,   vector2=column_Average_week4,    file1=paste(subdir_2_part3,  "/Part3-2-D-WT-week2-vs-week4.paired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=column_Average_week4,   vector2=column_Average_week6,    file1=paste(subdir_2_part3,  "/Part3-2-D-WT-week4-vs-week6.paired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=column_Average_week6,   vector2=column_Average_week8,    file1=paste(subdir_2_part3,  "/Part3-2-D-WT-week6-vs-week8.paired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=column_Average_week0_EEDheto,   vector2=column_Average_week0_EEDko,      file1=paste(subdir_2_part3,  "/Part3-2-E-CKO-week0Heto-vs-week0KO.paired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=column_Average_week0_EEDheto,   vector2=column_Average_week4_EEDheto,    file1=paste(subdir_2_part3,  "/Part3-2-E-CKO-week0Heto-vs-week4Heto.paired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=column_Average_week0_EEDko,     vector2=column_Average_week4_EEDko,      file1=paste(subdir_2_part3,  "/Part3-2-E-CKO-week0KO-vs-week4KO.paired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=column_Average_week4_EEDheto,   vector2=column_Average_week4_EEDko,      file1=paste(subdir_2_part3,  "/Part3-2-E-CKO-week4Heto-vs-week4KO.paired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=column_Average_banding,   vector2=column_Average_sham,    file1=paste(subdir_2_part3,  "/Part3-2-F-TAC-banding-vs-sham.paired.txt",  sep = "")   )   

## T test and Wilcoxon test  (unpaired).
MyHypothesisTest_1(vector1=column_Average_H3,      vector2=column_Average_week0,    file1=paste(subdir_2_part3,  "/Part3-2-G-WT-H3-vs-week0.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=column_Average_week0,   vector2=column_Average_week1,    file1=paste(subdir_2_part3,  "/Part3-2-G-WT-week0-vs-week1.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=column_Average_week1,   vector2=column_Average_week2,    file1=paste(subdir_2_part3,  "/Part3-2-G-WT-week1-vs-week2.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=column_Average_week2,   vector2=column_Average_week4,    file1=paste(subdir_2_part3,  "/Part3-2-G-WT-week2-vs-week4.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=column_Average_week4,   vector2=column_Average_week6,    file1=paste(subdir_2_part3,  "/Part3-2-G-WT-week4-vs-week6.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=column_Average_week6,   vector2=column_Average_week8,    file1=paste(subdir_2_part3,  "/Part3-2-G-WT-week6-vs-week8.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=column_Average_week0_EEDheto,   vector2=column_Average_week0_EEDko,      file1=paste(subdir_2_part3,  "/Part3-2-H-CKO-week0Heto-vs-week0KO.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=column_Average_week0_EEDheto,   vector2=column_Average_week4_EEDheto,    file1=paste(subdir_2_part3,  "/Part3-2-H-CKO-week0Heto-vs-week4Heto.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=column_Average_week0_EEDko,     vector2=column_Average_week4_EEDko,      file1=paste(subdir_2_part3,  "/Part3-2-H-CKO-week0KO-vs-week4KO.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=column_Average_week4_EEDheto,   vector2=column_Average_week4_EEDko,      file1=paste(subdir_2_part3,  "/Part3-2-H-CKO-week4Heto-vs-week4KO.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=column_Average_banding,   vector2=column_Average_sham,    file1=paste(subdir_2_part3,  "/Part3-2-I-TAC-banding-vs-sham.unpaired.txt",  sep = "")   )   










## Box-Violin-Plot, ±5kb
###############################################################################  Box-Violin-Plot, ±5kb
subdir_3_part3 <- paste(Part3_g,  "/3-AverageRows-allCols-boxViolin", sep = "")
if( ! file.exists(subdir_3_part3) ) { dir.create(subdir_3_part3) }

sink( file=paste(subdir_3_part3,  "/Part3-3-A-runLog.txt",      sep = "") )

MyBoxViolinPlot_1(vector2=c(row_week0_Rep1,  row_week0_Rep2,  row_week1_Rep1,  row_week1_Rep2,  row_week2_Rep1,  row_week2_Rep2,  
                            row_week4_Rep1,  row_week4_Rep2,  row_week4_Rep3,  row_week6_Rep1,  row_week6_Rep2,  row_week8_Rep1 ),   
                  sampleType2=c( rep("week0_Rep1", numOfRows1),    rep("week0_Rep2", numOfRows1),    rep("week1_Rep1", numOfRows1),  rep("week1_Rep2", numOfRows1),  
                                 rep("week2_Rep1", numOfRows1),    rep("week2_Rep2", numOfRows1),  
                                 rep("week4_Rep1", numOfRows1),    rep("week4_Rep2", numOfRows1),    rep("week4_Rep3", numOfRows1),    
                                 rep("week6_Rep1", numOfRows1),    rep("week6_Rep2", numOfRows1),    rep("week8_Rep1",   numOfRows1)  ), 
                  sampleRank2=c( "week0_Rep1",  "week0_Rep2",  "week1_Rep1",   "week1_Rep2",  "week2_Rep1",  "week2_Rep2",  
                                 "week4_Rep1",  "week4_Rep2",  "week4_Rep3",   "week6_Rep1",  "week6_Rep2",  "week8_Rep1" ),     
                  colours2=c( "week0_Rep1"="red",       "week0_Rep2"="red4",     "week1_Rep1"="cyan",    "week1_Rep2"="cyan4",    
                              "week2_Rep1"="green",     "week2_Rep2"="green4",  
                              "week4_Rep1"="skyblue",   "week4_Rep2"="blue",    "week4_Rep3"="blue4",  
                              "week6_Rep1"="purple",    "week6_Rep2"="purple4",  "week8_Rep1"="orange2" ), 
                  path2=subdir_3_part3,    fileName2=paste("Part3-3-A-WT-12Samples-BoxViolin",  "5kb", mySample_g, sep = "_") , 
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4.27,   width2=7,   Ymin2=0, Ymax2=1.5)   ## width = 1 + 12*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c(row_week0_Rep1,  row_week0_Rep2,  row_week1_Rep1,  row_week1_Rep2,  row_week2_Rep1,  row_week2_Rep2,  
                            row_week4_Rep2,  row_week4_Rep3,  row_week6_Rep2,  row_week8_Rep1 ),   
                  sampleType2=c( rep("week0_Rep1", numOfRows1),    rep("week0_Rep2", numOfRows1),    rep("week1_Rep1", numOfRows1),  rep("week1_Rep2", numOfRows1),  
                                 rep("week2_Rep1", numOfRows1),    rep("week2_Rep2", numOfRows1),  
                                 rep("week4_Rep2", numOfRows1),    rep("week4_Rep3", numOfRows1),    
                                 rep("week6_Rep2", numOfRows1),    rep("week8_Rep1",   numOfRows1)  ), 
                  sampleRank2=c( "week0_Rep1",  "week0_Rep2",  "week1_Rep1",   "week1_Rep2",  "week2_Rep1",  "week2_Rep2",  
                                 "week4_Rep2",  "week4_Rep3",  "week6_Rep2",  "week8_Rep1" ),     
                  colours2=c( "week0_Rep1"="red",       "week0_Rep2"="red4",     "week1_Rep1"="cyan",    "week1_Rep2"="cyan4",    
                              "week2_Rep1"="green",     "week2_Rep2"="green4",  
                              "week4_Rep2"="blue",    "week4_Rep3"="blue4",  
                              "week6_Rep2"="purple4",  "week8_Rep1"="orange2" ), 
                  path2=subdir_3_part3,    fileName2=paste("Part3-3-B-WT-10Samples-BoxViolin",  "5kb", mySample_g, sep = "_") , 
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4.27,   width2=6,   Ymin2=0, Ymax2=1.5)   ## width = 1 + 10*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c(row_Average_week0, row_Average_week1,  row_Average_week2,  row_Average_week4,  row_Average_week6, row_Average_week8  ),   
                  sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                                 rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ), 
                  sampleRank2=c( "week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                  colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                  path2=subdir_3_part3,   fileName2= paste("Part3-3-C-WT-6Samples-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4.0,   width2=4,   Ymin2=0, Ymax2=1.5)    ## width = 1 + 6*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c(row_Average_week0, row_Average_H3  ),   
                  sampleType2=c( rep("week0", numOfRows1),    rep("H3", numOfRows1)   ), 
                  sampleRank2=c( "week0",  "H3"  ),     
                  colours2=c( "week0"="red2",   "H3"="yellow2"  ), 
                  path2=subdir_3_part3,   fileName2= paste("Part3-3-D-WT-2Samples-BoxViolin",  "5kb",   mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4.0,   width2=2,   Ymin2=0, Ymax2=1.5)  ## width = 1 + 2*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c(row_Average_week0, row_Average_H3  ),   
                  sampleType2=c( rep("week0", numOfRows1),    rep("H3", numOfRows1)   ), 
                  sampleRank2=c( "week0",  "H3"  ),     
                  colours2=c( "week0"="red2",   "H3"="yellow2"  ), 
                  path2=subdir_3_part3,   fileName2= paste("Part3-3-E-WT-2Samples-BoxViolin",  "5kb",   mySample_g, sep = "_") , 
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4.0,   width2=2,   Ymin2=0, Ymax2=2.5)   ## width = 1 + 2*0.5, height=5cm 





MyBoxViolinPlot_1(vector2=c(row_week0_EEDheto_Rep1, row_week0_EEDheto_Rep2, row_week0_EEDko_Rep1,   row_week4_EEDheto_Rep1, row_week4_EEDheto_Rep2, row_week4_EEDko_Rep1, row_week4_EEDko_Rep2 ),   
                  sampleType2=c( rep("week0_EEDheto_Rep1", numOfRows1),  rep("week0_EEDheto_Rep2", numOfRows1),    rep("week0_EEDko_Rep1", numOfRows1),  
                                 rep("week4_EEDheto_Rep1", numOfRows1),  rep("week4_EEDheto_Rep2", numOfRows1),    rep("week4_EEDko_Rep1", numOfRows1),  rep("week4_EEDko_Rep2", numOfRows1)  ), 
                  sampleRank2=c( "week0_EEDheto_Rep1",  "week0_EEDheto_Rep2",  "week0_EEDko_Rep1",  "week4_EEDheto_Rep1",  "week4_EEDheto_Rep2",   "week4_EEDko_Rep1",   "week4_EEDko_Rep2" ),  
                  colours2=c( "week0_EEDheto_Rep1"="red",      "week0_EEDheto_Rep2"="red4",     "week0_EEDko_Rep1"="orange2",   
                              "week4_EEDheto_Rep1"="blue",     "week4_EEDheto_Rep2"="blue4",    "week4_EEDko_Rep1"="green", "week4_EEDko_Rep2"="green4" ),  
                  path2=subdir_3_part3,    fileName2= paste("Part3-3-F-EEDko-7Samples-BoxViolin",  "5kb",   mySample_g, sep = "_") , 
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4.65,   width2=4.5,   Ymin2=0, Ymax2=1.0)  ## width = 1 + 7*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c(row_week0_EEDheto_Rep1, row_week0_EEDheto_Rep2, row_week0_EEDko_Rep1,   row_week4_EEDheto_Rep2, row_week4_EEDko_Rep1, row_week4_EEDko_Rep2 ),   
                  sampleType2=c( rep("week0_EEDheto_Rep1", numOfRows1),  rep("week0_EEDheto_Rep2", numOfRows1),    rep("week0_EEDko_Rep1", numOfRows1),  
                                 rep("week4_EEDheto_Rep2", numOfRows1),    rep("week4_EEDko_Rep1", numOfRows1),  rep("week4_EEDko_Rep2", numOfRows1)  ), 
                  sampleRank2=c( "week0_EEDheto_Rep1",  "week0_EEDheto_Rep2",  "week0_EEDko_Rep1",  "week4_EEDheto_Rep2",   "week4_EEDko_Rep1",   "week4_EEDko_Rep2" ),  
                  colours2=c( "week0_EEDheto_Rep1"="red",      "week0_EEDheto_Rep2"="red4",     "week0_EEDko_Rep1"="orange2",   
                              "week4_EEDheto_Rep2"="blue4",    "week4_EEDko_Rep1"="green", "week4_EEDko_Rep2"="green4" ),  
                  path2=subdir_3_part3,    fileName2= paste("Part3-3-G-EEDko-6Samples-BoxViolin",  "5kb",   mySample_g, sep = "_") , 
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4.65,   width2=4.0,   Ymin2=0, Ymax2=1.0)  ## width = 1 + 6*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c(row_Average_week0_EEDheto,  row_Average_week0_EEDko, row_Average_week4_EEDheto,  row_Average_week4_EEDko ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1), rep("week4_EEDheto", numOfRows1),  rep("week4_EEDko", numOfRows1) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ),  
                  path2=subdir_3_part3,  fileName2= paste("Part3-3-H-EEDko-4Samples-BoxViolin",  "5kb",   mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=4.4,   width2=3.0, Ymin2=0, Ymax2=1.0 )  ## width = 1 + 4*0.5, height=5cm 





MyBoxViolinPlot_1(vector2=c(row_banding_Rep1,  row_banding_Rep2,  row_sham_Rep1,  row_sham_Rep2 ),   
                  sampleType2=c( rep("banding_Rep1", numOfRows1),  rep("banding_Rep2", numOfRows1),    rep("sham_Rep1", numOfRows1),  rep("sham_Rep2", numOfRows1)  ), 
                  sampleRank2=c( "banding_Rep1",   "banding_Rep2",   "sham_Rep1",  "sham_Rep2"  ),     
                  colours2=c( "banding_Rep1"="blue",  "banding_Rep2"="blue4",  "sham_Rep1"="green",  "sham_Rep2"="green4"  ), 
                  path2=subdir_3_part3,  fileName2= paste("Part3-3-I-TAC-4Samples-BoxViolin",  "5kb",   mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=4.3,   width2=3.0, Ymin2=0, Ymax2=0.8 )   ## width = 1 + 4*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c(row_Average_banding,  row_Average_sham ),   
                  sampleType2=c( rep("banding", numOfRows1),    rep("sham", numOfRows1)  ), 
                  sampleRank2=c( "banding",   "sham"  ),     
                  colours2=c( "banding"="blue",  "sham"="green4"  ), 
                  path2=subdir_3_part3,  fileName2= paste("Part3-3-J-TAC-2Samples-BoxViolin",  "5kb",   mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=4.07,   width2=2.0, Ymin2=0, Ymax2=0.8 )   ## width = 1 + 2*0.5, height=5cm 


sink()  


















## hypothesis test
################################################################################################
subdir_4_part3 <- paste(Part3_g,  "/4-averageRows-hypothesisTest", sep = "")
if( ! file.exists(subdir_4_part3) ) { dir.create(subdir_4_part3) }

## T test and Wilcoxon test  (unpaired).
MyHypothesisTest_1(vector1=row_Average_H3,      vector2=row_Average_week0,    file1=paste(subdir_4_part3,  "/Part3-4-A-WT-H3-vs-week0.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=row_Average_week0,   vector2=row_Average_week1,    file1=paste(subdir_4_part3,  "/Part3-4-A-WT-week0-vs-week1.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=row_Average_week1,   vector2=row_Average_week2,    file1=paste(subdir_4_part3,  "/Part3-4-A-WT-week1-vs-week2.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=row_Average_week2,   vector2=row_Average_week4,    file1=paste(subdir_4_part3,  "/Part3-4-A-WT-week2-vs-week4.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=row_Average_week4,   vector2=row_Average_week6,    file1=paste(subdir_4_part3,  "/Part3-4-A-WT-week4-vs-week6.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=row_Average_week6,   vector2=row_Average_week8,    file1=paste(subdir_4_part3,  "/Part3-4-A-WT-week6-vs-week8.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=row_Average_week0_EEDheto,   vector2=row_Average_week0_EEDko,      file1=paste(subdir_4_part3,  "/Part3-4-B-CKO-week0Heto-vs-week0KO.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=row_Average_week0_EEDheto,   vector2=row_Average_week4_EEDheto,    file1=paste(subdir_4_part3,  "/Part3-4-B-CKO-week0Heto-vs-week4Heto.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=row_Average_week0_EEDko,     vector2=row_Average_week4_EEDko,      file1=paste(subdir_4_part3,  "/Part3-4-B-CKO-week0KO-vs-week4KO.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=row_Average_week4_EEDheto,   vector2=row_Average_week4_EEDko,      file1=paste(subdir_4_part3,  "/Part3-4-B-CKO-week4Heto-vs-week4KO.unpaired.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=row_Average_banding,   vector2=row_Average_sham,    file1=paste(subdir_4_part3,  "/Part3-4-C-TAC-banding-vs-sham.unpaired.txt",  sep = "")   )   


## T test and Wilcoxon test  (paired).
MyHypothesisTest_2(vector1=row_Average_H3,      vector2=row_Average_week0,    file1=paste(subdir_4_part3,  "/Part3-4-D-WT-H3-vs-week0.unpaired.txt",     sep = "")   )   
MyHypothesisTest_2(vector1=row_Average_week0,   vector2=row_Average_week1,    file1=paste(subdir_4_part3,  "/Part3-4-D-WT-week0-vs-week1.unpaired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=row_Average_week1,   vector2=row_Average_week2,    file1=paste(subdir_4_part3,  "/Part3-4-D-WT-week1-vs-week2.unpaired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=row_Average_week2,   vector2=row_Average_week4,    file1=paste(subdir_4_part3,  "/Part3-4-D-WT-week2-vs-week4.unpaired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=row_Average_week4,   vector2=row_Average_week6,    file1=paste(subdir_4_part3,  "/Part3-4-D-WT-week4-vs-week6.unpaired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=row_Average_week6,   vector2=row_Average_week8,    file1=paste(subdir_4_part3,  "/Part3-4-D-WT-week6-vs-week8.unpaired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=row_Average_week0_EEDheto,   vector2=row_Average_week0_EEDko,      file1=paste(subdir_4_part3,  "/Part3-4-E-CKO-week0Heto-vs-week0KO.unpaired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=row_Average_week0_EEDheto,   vector2=row_Average_week4_EEDheto,    file1=paste(subdir_4_part3,  "/Part3-4-E-CKO-week0Heto-vs-week4Heto.unpaired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=row_Average_week0_EEDko,     vector2=row_Average_week4_EEDko,      file1=paste(subdir_4_part3,  "/Part3-4-E-CKO-week0KO-vs-week4KO.unpaired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=row_Average_week4_EEDheto,   vector2=row_Average_week4_EEDko,      file1=paste(subdir_4_part3,  "/Part3-4-E-CKO-week4Heto-vs-week4KO.unpaired.txt",  sep = "")   )   
MyHypothesisTest_2(vector1=row_Average_banding,   vector2=row_Average_sham,    file1=paste(subdir_4_part3,  "/Part3-4-F-TAC-banding-vs-sham.unpaired.txt",  sep = "")   )   











## Box-Violin-Plot, ±2kb
###############################################################################  Box-Violin-Plot, ±2kb
subdir_5_part3 <- paste(Part3_g,  "/5-AverageRows-2kb-boxViolin", sep = "")
if( ! file.exists(subdir_5_part3) ) { dir.create(subdir_5_part3) }

sink( file=paste(subdir_5_part3,  "/Part3-5-A-runLog.txt",      sep = "") )

MyBoxViolinPlot_1(vector2=c(reduceColumn2_Average_week0, reduceColumn2_Average_week1,  reduceColumn2_Average_week2,  reduceColumn2_Average_week4,  reduceColumn2_Average_week6, reduceColumn2_Average_week8  ),   
                  sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                                 rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ), 
                  sampleRank2=c( "week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                  colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                  path2=subdir_5_part3,   fileName2= paste("Part3-5-A-WT-6Samples-BoxViolin",  "2kb",   mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4,   width2=4,   Ymin2=0, Ymax2=1.5)    ## width = 1 + 6*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c(reduceColumn2_Average_week0, reduceColumn2_Average_H3  ),   
                  sampleType2=c( rep("week0", numOfRows1),    rep("H3", numOfRows1)   ), 
                  sampleRank2=c( "week0",  "H3"  ),     
                  colours2=c( "week0"="red2",   "H3"="yellow2"  ), 
                  path2=subdir_5_part3,   fileName2= paste("Part3-5-B-WT-2Samples-BoxViolin",  "2kb",   mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4,   width2=2,   Ymin2=0, Ymax2=1.5)   ## width = 1 + 2*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c(reduceColumn2_Average_week0, reduceColumn2_Average_H3  ),   
                  sampleType2=c( rep("week0", numOfRows1),    rep("H3", numOfRows1)   ), 
                  sampleRank2=c( "week0",  "H3"  ),     
                  colours2=c( "week0"="red2",   "H3"="yellow2"  ), 
                  path2=subdir_5_part3,   fileName2= paste("Part3-5-C-WT-2Samples-BoxViolin",  "2kb",   mySample_g, sep = "_") , 
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4,   width2=2,   Ymin2=0, Ymax2=2.5)   ## width = 1 + 2*0.5, height=5cm 





MyBoxViolinPlot_1(vector2=c(reduceColumn2_Average_week0_EEDheto,  reduceColumn2_Average_week0_EEDko, reduceColumn2_Average_week4_EEDheto,  reduceColumn2_Average_week4_EEDko ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1), rep("week4_EEDheto", numOfRows1),  rep("week4_EEDko", numOfRows1) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ),  
                  path2=subdir_5_part3,  fileName2= paste("Part3-5-D-EEDko-4Samples-BoxViolin", "2kb",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=4.4,   width2=3.0, Ymin2=0, Ymax2=1.0 ) ## width = 1 + 4*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c(reduceColumn2_Average_banding,  reduceColumn2_Average_sham ),   
                  sampleType2=c( rep("banding", numOfRows1),    rep("sham", numOfRows1)  ), 
                  sampleRank2=c( "banding",   "sham"  ),     
                  colours2=c( "banding"="blue",  "sham"="green4"  ), 
                  path2=subdir_5_part3,  fileName2= paste("Part3-5-E-TAC-2Samples-BoxViolin", "2kb",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=4.2,   width2=2, Ymin2=0, Ymax2=0.8 ) ## width = 1 + 2*0.5, height=5cm 


sink()  
















## Box-Violin-Plot, ±1kb
###############################################################################  Box-Violin-Plot, ±1kb
subdir_6_part3 <- paste(Part3_g,  "/6-AverageRows-1kb-boxViolin", sep = "")
if( ! file.exists(subdir_6_part3) ) { dir.create(subdir_6_part3) }

sink( file=paste(subdir_6_part3,  "/Part3-6-A-runLog.txt",      sep = "") )

MyBoxViolinPlot_1(vector2=c(reduceColumn3_Average_week0, reduceColumn3_Average_week1,  reduceColumn3_Average_week2,  reduceColumn3_Average_week4,  reduceColumn3_Average_week6, reduceColumn3_Average_week8  ),   
                  sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                                 rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ), 
                  sampleRank2=c( "week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                  colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                  path2=subdir_6_part3,   fileName2= paste("Part3-6-A-WT-6Samples-BoxViolin",  "1kb",   mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4,   width2=4,   Ymin2=0, Ymax2=1.5)   ## width = 1 + 6*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c(reduceColumn3_Average_week0, reduceColumn3_Average_H3  ),   
                  sampleType2=c( rep("week0", numOfRows1),    rep("H3", numOfRows1)   ), 
                  sampleRank2=c( "week0",  "H3"  ),     
                  colours2=c( "week0"="red2",   "H3"="yellow2"  ), 
                  path2=subdir_6_part3,   fileName2= paste("Part3-6-B-WT-2Samples-BoxViolin",  "1kb",   mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4,   width2=2,   Ymin2=0, Ymax2=1.5)  ## width = 1 + 2*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c(reduceColumn3_Average_week0, reduceColumn3_Average_H3  ),   
                  sampleType2=c( rep("week0", numOfRows1),    rep("H3", numOfRows1)   ), 
                  sampleRank2=c( "week0",  "H3"  ),     
                  colours2=c( "week0"="red2",   "H3"="yellow2"  ), 
                  path2=subdir_6_part3,   fileName2= paste("Part3-6-C-WT-2Samples-BoxViolin",  "1kb",   mySample_g, sep = "_") , 
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4,   width2=2,   Ymin2=0, Ymax2=2.5)  ## width = 1 + 2*0.5, height=5cm 





MyBoxViolinPlot_1(vector2=c(reduceColumn3_Average_week0_EEDheto,  reduceColumn3_Average_week0_EEDko, reduceColumn3_Average_week4_EEDheto,  reduceColumn3_Average_week4_EEDko ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1), rep("week4_EEDheto", numOfRows1),  rep("week4_EEDko", numOfRows1) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ),  
                  path2=subdir_6_part3,  fileName2= paste("Part3-6-D-EEDko-4Samples-BoxViolin", "1kb",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=4.4,   width2=3.0, Ymin2=0, Ymax2=1.0 )   ## width = 1 + 4*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c(reduceColumn3_Average_banding,  reduceColumn3_Average_sham ),   
                  sampleType2=c( rep("banding", numOfRows1),    rep("sham", numOfRows1)  ), 
                  sampleRank2=c( "banding",   "sham"  ),     
                  colours2=c( "banding"="blue",  "sham"="green4"  ), 
                  path2=subdir_6_part3,  fileName2= paste("Part3-6-E-TAC-2Samples-BoxViolin", "1kb",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=4.2,   width2=2, Ymin2=0, Ymax2=0.8 )   ## width = 1 + 2*0.5, height=5cm 


sink()  













## Box-Violin-Plot, ±500bp
###############################################################################  Box-Violin-Plot, ±500bp
subdir_7_part3 <- paste(Part3_g,  "/7-AverageRows-500bp-boxViolin", sep = "")
if( ! file.exists(subdir_7_part3) ) { dir.create(subdir_7_part3) }

sink( file=paste(subdir_7_part3,  "/Part3-7-A-runLog.txt",      sep = "") )

MyBoxViolinPlot_1(vector2=c(reduceColumn4_Average_week0, reduceColumn4_Average_week1,  reduceColumn4_Average_week2,  reduceColumn4_Average_week4,  reduceColumn4_Average_week6, reduceColumn4_Average_week8  ),   
                  sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                                 rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ), 
                  sampleRank2=c( "week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                  colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                  path2=subdir_7_part3,   fileName2= paste("Part3-7-A-WT-6Samples-BoxViolin",  "500bp",   mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4,   width2=4,   Ymin2=0, Ymax2=1.5)    ## width = 1 + 6*0.5,  height=5cm

MyBoxViolinPlot_1(vector2=c(reduceColumn4_Average_week0, reduceColumn4_Average_H3  ),   
                  sampleType2=c( rep("week0", numOfRows1),    rep("H3", numOfRows1)   ), 
                  sampleRank2=c( "week0",  "H3"  ),     
                  colours2=c( "week0"="red2",   "H3"="yellow2"  ), 
                  path2=subdir_7_part3,   fileName2= paste("Part3-7-B-WT-2Samples-BoxViolin",  "500bp",   mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4,   width2=2,   Ymin2=0, Ymax2=1.5)  ## width = 1 + 2*0.5,  height=5cm

MyBoxViolinPlot_1(vector2=c(reduceColumn4_Average_week0, reduceColumn4_Average_H3  ),   
                  sampleType2=c( rep("week0", numOfRows1),    rep("H3", numOfRows1)   ), 
                  sampleRank2=c( "week0",  "H3"  ),     
                  colours2=c( "week0"="red2",   "H3"="yellow2"  ), 
                  path2=subdir_7_part3,   fileName2= paste("Part3-7-C-WT-2Samples-BoxViolin",  "500bp",   mySample_g, sep = "_") , 
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4,   width2=2,   Ymin2=0, Ymax2=2.5)  ## width = 1 + 2*0.5,  height=5cm





MyBoxViolinPlot_1(vector2=c(reduceColumn4_Average_week0_EEDheto,  reduceColumn4_Average_week0_EEDko, reduceColumn4_Average_week4_EEDheto,  reduceColumn4_Average_week4_EEDko ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1), rep("week4_EEDheto", numOfRows1),  rep("week4_EEDko", numOfRows1) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ),  
                  path2=subdir_7_part3,  fileName2= paste("Part3-7-D-EEDko-4Samples-BoxViolin", "500bp",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=4.4,   width2=3.0, Ymin2=0, Ymax2=1.0)    ## width = 1 + 4*0.5,  height=5cm


MyBoxViolinPlot_1(vector2=c(reduceColumn4_Average_banding,  reduceColumn4_Average_sham ),   
                  sampleType2=c( rep("banding", numOfRows1),    rep("sham", numOfRows1)  ), 
                  sampleRank2=c( "banding",   "sham"  ),     
                  colours2=c( "banding"="blue",  "sham"="green4"  ), 
                  path2=subdir_7_part3,  fileName2= paste("Part3-7-E-TAC-2Samples-BoxViolin", "500bp",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=4.2,   width2=2, Ymin2=0, Ymax2=0.8 )   ## width = 1 + 2*0.5,  height=5cm


sink()  





























###############################################################################
subdir_8_part3 <- paste(Part3_g,  "/8-averageRows-1kb-histogram", sep = "")
if( ! file.exists(subdir_8_part3) ) { dir.create(subdir_8_part3) }

MyHistogram_1(vector2=reduceColumn3_Average_H3,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-A-WT-H3-histogram",  
              title2=myTitle_g,  
              xLab2="H3 signal",   
              height2=3.25,  width2=4.25,   xMin2=0,  xMax2=1.51,  yMin2=0,  yMax2=0.07)  ## height = 4cm,  width = 6cm

MyHistogram_1(vector2=reduceColumn3_Average_week0,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-B-WT-week0-histogram",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=4.25,   xMin2=0,  xMax2=1.51,  yMin2=0,  yMax2=0.07)

MyHistogram_1(vector2=reduceColumn3_Average_week1,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-C-WT-week1-histogram",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=4.25,   xMin2=0,  xMax2=1.51,  yMin2=0,  yMax2=0.07)

MyHistogram_1(vector2=reduceColumn3_Average_week2,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-D-WT-week2-histogram",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=4.25,   xMin2=0,  xMax2=1.51,  yMin2=0,  yMax2=0.07)

MyHistogram_1(vector2=reduceColumn3_Average_week4,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-E-WT-week4-histogram",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=4.25,   xMin2=0,  xMax2=1.51,  yMin2=0,  yMax2=0.07)

MyHistogram_1(vector2=reduceColumn3_Average_week6,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-F-WT-week6-histogram",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=4.25,   xMin2=0,  xMax2=1.51,  yMin2=0,  yMax2=0.07)

MyHistogram_1(vector2=reduceColumn3_Average_week8,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-G-WT-week8-histogram",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=4.25,   xMin2=0,  xMax2=1.51,  yMin2=0,  yMax2=0.07)






MyHistogram_1A(vector2=reduceColumn3_Average_week0_EEDheto,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-H-CKO-week0-EEDheto-histogram",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=4.25,   xMin2=0,  xMax2=1.01,  yMin2=0,  yMax2=0.07)

MyHistogram_1A(vector2=reduceColumn3_Average_week0_EEDko,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-I-CKO-week0-EEDko-histogram",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=4.25,   xMin2=0,  xMax2=1.01,  yMin2=0,  yMax2=0.07)

MyHistogram_1A(vector2=reduceColumn3_Average_week4_EEDheto,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-J-CKO-week4-EEDheto-histogram",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=4.25,   xMin2=0,  xMax2=1.01,  yMin2=0,  yMax2=0.07)

MyHistogram_1A(vector2=reduceColumn3_Average_week4_EEDko,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-K-CKO-week4-EEDko-histogram",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=4.25,   xMin2=0,  xMax2=1.01,  yMin2=0,  yMax2=0.07)







MyHistogram_1B(vector2=reduceColumn3_Average_banding,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-L-TAC-banding-histogram",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=4.25,   xMin2=0,  xMax2=0.81,  yMin2=0,  yMax2=0.08)

MyHistogram_1B(vector2=reduceColumn3_Average_sham,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-M-TAC-sham-histogram",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=4.25,   xMin2=0,  xMax2=0.81,  yMin2=0,  yMax2=0.08)












MyHistogram_3(vector2=reduceColumn3_Average_H3,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-Z-A-WT-H3-density",  
              title2=myTitle_g,  
              xLab2="H3 signal",   
              height2=3.25,  width2=5.5,   xMin2=0,  xMax2=1.21,  yMin2=0,  yMax2=4)  ## height = 4cm,  width = 6cm

MyHistogram_3(vector2=reduceColumn3_Average_week0,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-Z-B-WT-week0-density",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=5.5,   xMin2=0,  xMax2=1.21,  yMin2=0,  yMax2=4)

MyHistogram_3(vector2=reduceColumn3_Average_week1,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-Z-C-WT-week1-density",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=5.5,   xMin2=0,  xMax2=1.21,  yMin2=0,  yMax2=4)

MyHistogram_3(vector2=reduceColumn3_Average_week2,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-Z-D-WT-week2-density",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=5.5,   xMin2=0,  xMax2=1.21,  yMin2=0,  yMax2=4)

MyHistogram_3(vector2=reduceColumn3_Average_week4,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-Z-E-WT-week4-density",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=5.5,   xMin2=0,  xMax2=1.21,  yMin2=0,  yMax2=4)

MyHistogram_3(vector2=reduceColumn3_Average_week6,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-Z-F-WT-week6-density",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=5.5,   xMin2=0,  xMax2=1.21,  yMin2=0,  yMax2=4)

MyHistogram_3(vector2=reduceColumn3_Average_week8,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-Z-G-WT-week8-density",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=5.5,   xMin2=0,  xMax2=1.21,  yMin2=0,  yMax2=4)






MyHistogram_3(vector2=reduceColumn3_Average_week0_EEDheto,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-Z-H-CKO-week0-EEDheto-density",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=5.5,   xMin2=0,  xMax2=1.21,  yMin2=0,  yMax2=4)

MyHistogram_3(vector2=reduceColumn3_Average_week0_EEDko,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-Z-I-CKO-week0-EEDko-density",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=5.5,   xMin2=0,  xMax2=1.21,  yMin2=0,  yMax2=4)

MyHistogram_3(vector2=reduceColumn3_Average_week4_EEDheto,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-Z-J-CKO-week4-EEDheto-density",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=5.5,   xMin2=0,  xMax2=1.21,  yMin2=0,  yMax2=4)

MyHistogram_3(vector2=reduceColumn3_Average_week4_EEDko,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-Z-K-CKO-week4-EEDko-density",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=5.5,   xMin2=0,  xMax2=1.21,  yMin2=0,  yMax2=4)







MyHistogram_3(vector2=reduceColumn3_Average_banding,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-Z-L-TAC-banding-density",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=5.5,   xMin2=0,  xMax2=1.21,  yMin2=0,  yMax2=4)

MyHistogram_3(vector2=reduceColumn3_Average_sham,  
              path2=subdir_8_part3,    
              fileName2="Part3-8-Z-M-TAC-sham-density",  
              title2=myTitle_g,  
              xLab2="H2BGFP signal",   
              height2=3.25,  width2=5.5,   xMin2=0,  xMax2=1.21,  yMin2=0,  yMax2=4)










pdf( file=paste(subdir_8_part3, "/", "Part3-8-A-2kb-all-densityCurve.pdf", sep="") )
thres_1 <- 1.5
par(mfrow=c(2, 1))
plot( density(reduceColumn2_Average_H3[ reduceColumn2_Average_H3 < thres_1]), xlim=c(0, thres_1)  )
plot( density(reduceColumn2_Average_week0[ reduceColumn2_Average_week0 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn2_Average_week0[ reduceColumn2_Average_week0 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn2_Average_week1[ reduceColumn2_Average_week1 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn2_Average_week2[ reduceColumn2_Average_week2 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn2_Average_week4[ reduceColumn2_Average_week4 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn2_Average_week6[ reduceColumn2_Average_week6 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn2_Average_week8[ reduceColumn2_Average_week8 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn2_Average_week0_EEDheto[ reduceColumn2_Average_week0_EEDheto < thres_1])  , xlim=c(0, thres_1)  )  
plot( density(reduceColumn2_Average_week0_EEDko[ reduceColumn2_Average_week0_EEDko < thres_1]) , xlim=c(0, thres_1)   )
plot( density(reduceColumn2_Average_week4_EEDheto[ reduceColumn2_Average_week4_EEDheto < thres_1]) , xlim=c(0, thres_1)   )
plot( density(reduceColumn2_Average_week4_EEDko[ reduceColumn2_Average_week4_EEDko < thres_1]) , xlim=c(0, thres_1)   )
plot( density(reduceColumn2_Average_banding[ reduceColumn2_Average_banding < thres_1]) , xlim=c(0, thres_1)   )
plot( density(reduceColumn2_Average_sham[ reduceColumn2_Average_sham < thres_1])  , xlim=c(0, thres_1)  )
dev.off()






pdf( file=paste(subdir_8_part3, "/", "Part3-8-B-1kb-all-densityCurve.pdf", sep="") )
thres_1 <- 1.5
par(mfrow=c(2, 1))
plot( density(reduceColumn3_Average_H3[ reduceColumn3_Average_H3 < thres_1]), xlim=c(0, thres_1)  )
plot( density(reduceColumn3_Average_week0[ reduceColumn3_Average_week0 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn3_Average_week0[ reduceColumn3_Average_week0 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn3_Average_week1[ reduceColumn3_Average_week1 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn3_Average_week2[ reduceColumn3_Average_week2 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn3_Average_week4[ reduceColumn3_Average_week4 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn3_Average_week6[ reduceColumn3_Average_week6 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn3_Average_week8[ reduceColumn3_Average_week8 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn3_Average_week0_EEDheto[ reduceColumn3_Average_week0_EEDheto < thres_1])  , xlim=c(0, thres_1)  )  
plot( density(reduceColumn3_Average_week0_EEDko[ reduceColumn3_Average_week0_EEDko < thres_1]) , xlim=c(0, thres_1)   )
plot( density(reduceColumn3_Average_week4_EEDheto[ reduceColumn3_Average_week4_EEDheto < thres_1]) , xlim=c(0, thres_1)   )
plot( density(reduceColumn3_Average_week4_EEDko[ reduceColumn3_Average_week4_EEDko < thres_1]) , xlim=c(0, thres_1)   )
plot( density(reduceColumn3_Average_banding[ reduceColumn3_Average_banding < thres_1]) , xlim=c(0, thres_1)   )
plot( density(reduceColumn3_Average_sham[ reduceColumn3_Average_sham < thres_1])  , xlim=c(0, thres_1)  )
dev.off()








pdf( file=paste(subdir_8_part3, "/", "Part3-8-C-500bp-all-densityCurve.pdf", sep="") )
thres_1 <- 1.5
par(mfrow=c(2, 1))
plot( density(reduceColumn4_Average_H3[ reduceColumn4_Average_H3 < thres_1]), xlim=c(0, thres_1)  )
plot( density(reduceColumn4_Average_week0[ reduceColumn4_Average_week0 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn4_Average_week0[ reduceColumn4_Average_week0 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn4_Average_week1[ reduceColumn4_Average_week1 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn4_Average_week2[ reduceColumn4_Average_week2 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn4_Average_week4[ reduceColumn4_Average_week4 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn4_Average_week6[ reduceColumn4_Average_week6 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn4_Average_week8[ reduceColumn4_Average_week8 < thres_1])  , xlim=c(0, thres_1)  )
plot( density(reduceColumn4_Average_week0_EEDheto[ reduceColumn4_Average_week0_EEDheto < thres_1])  , xlim=c(0, thres_1)  )  
plot( density(reduceColumn4_Average_week0_EEDko[ reduceColumn4_Average_week0_EEDko < thres_1]) , xlim=c(0, thres_1)   )
plot( density(reduceColumn4_Average_week4_EEDheto[ reduceColumn4_Average_week4_EEDheto < thres_1]) , xlim=c(0, thres_1)   )
plot( density(reduceColumn4_Average_week4_EEDko[ reduceColumn4_Average_week4_EEDko < thres_1]) , xlim=c(0, thres_1)   )
plot( density(reduceColumn4_Average_banding[ reduceColumn4_Average_banding < thres_1]) , xlim=c(0, thres_1)   )
plot( density(reduceColumn4_Average_sham[ reduceColumn4_Average_sham < thres_1])  , xlim=c(0, thres_1)  )
dev.off()
























###############################################################################
subdir_9_part3 <- paste(Part3_g,  "/9-averageRows-5kb-Cmp", sep = "")
if( ! file.exists(subdir_9_part3) ) { dir.create(subdir_9_part3) }


MyHistogram_4(vector2=c(row_Average_week0, row_Average_week1,  row_Average_week2,  row_Average_week4,  row_Average_week6, row_Average_week8  ),  
              sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                             rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ),  
              path2=subdir_9_part3,    
              fileName2="Part3-9-A-WT-6Samples-stacking",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2  )


MyHistogram_5(vector2=c(row_Average_week0, row_Average_week1,  row_Average_week2,  row_Average_week4,  row_Average_week6, row_Average_week8  ),  
              sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                             rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ),  
              path2=subdir_9_part3,    
              fileName2="Part3-9-B-WT-6Samples-dodge",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2  )



MyHistogram_6(vector2=c(row_Average_week0, row_Average_week1,  row_Average_week2,  row_Average_week4,  row_Average_week6, row_Average_week8  ),  
              sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                             rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ),  
              colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
              path2=subdir_9_part3,    
              fileName2="Part3-9-C-WT-6Samples-density",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=10 )


MyHistogram_7(vector2=c(row_Average_week0, row_Average_week1,  row_Average_week2,  row_Average_week4,  row_Average_week6, row_Average_week8  ),  
              sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                             rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ),  
              colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
              path2=subdir_9_part3,    
              fileName2="Part3-9-D-WT-6Samples-density2",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=0.8,  yMin2=0,  yMax2=10 )






MyHistogram_6(vector2=c(row_Average_week0, row_Average_H3  ),  
              sampleType2=c( rep("week0", numOfRows1),    rep("H3", numOfRows1)   ),   
              colours2=c( "week0"="red2",   "H3"="yellow2"  ),  
              path2=subdir_9_part3,    
              fileName2="Part3-9-E-WT-2Samples-density",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=10 )


MyHistogram_7(vector2=c(row_Average_week0, row_Average_H3 ),  
              sampleType2=c( rep("week0", numOfRows1),    rep("H3", numOfRows1)   ),  
              colours2=c( "week0"="red2",   "H3"="yellow2"  ), 
              path2=subdir_9_part3,    
              fileName2="Part3-9-F-WT-2Samples-density2",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.5,  yMin2=0,  yMax2=10 )








MyHistogram_4(vector2=c(row_Average_week0_EEDheto,  row_Average_week0_EEDko, row_Average_week4_EEDheto,  row_Average_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1),  
                             rep("week4_EEDheto", numOfRows1),   rep("week4_EEDko", numOfRows1) ),  
              path2=subdir_9_part3,    
              fileName2="Part3-9-G-CKO-4Samples-stacking",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=100,  alpha2=1)

MyHistogram_5(vector2=c(row_Average_week0_EEDheto,  row_Average_week0_EEDko, row_Average_week4_EEDheto,  row_Average_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1),  
                             rep("week4_EEDheto", numOfRows1),   rep("week4_EEDko", numOfRows1) ),  
              path2=subdir_9_part3,    
              fileName2="Part3-9-H-CKO-4Samples-dodge",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=100 )


MyHistogram_6(vector2=c(row_Average_week0_EEDheto,  row_Average_week0_EEDko, row_Average_week4_EEDheto,  row_Average_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1),  
                             rep("week4_EEDheto", numOfRows1),   rep("week4_EEDko", numOfRows1) ),  
              colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
              path2=subdir_9_part3,    
              fileName2="Part3-9-I-CKO-4Samples-density",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=10 )

MyHistogram_7(vector2=c(row_Average_week0_EEDheto,  row_Average_week0_EEDko, row_Average_week4_EEDheto,  row_Average_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1),  
                             rep("week4_EEDheto", numOfRows1),   rep("week4_EEDko", numOfRows1) ),  
              colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
              path2=subdir_9_part3,    
              fileName2="Part3-9-J-CKO-4Samples-density2",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=0.8,  yMin2=0,  yMax2=10 )









MyHistogram_6(vector2=c(row_Average_banding,  row_Average_sham ),
              sampleType2=c( rep("banding", numOfRows1),    rep("sham", numOfRows1)  ), 
              colours2=c( "banding"="blue",  "sham"="green4"  ), 
              path2=subdir_9_part3,    
              fileName2="Part3-9-K-TAC-4Samples-density",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=10 )

MyHistogram_7(vector2=c(row_Average_banding,  row_Average_sham ),
              sampleType2=c( rep("banding", numOfRows1),    rep("sham", numOfRows1)  ), 
              colours2=c( "banding"="blue",  "sham"="green4"  ),  
              path2=subdir_9_part3,    
              fileName2="Part3-9-L-TAC-4Samples-density2",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=0.8,  yMin2=0,  yMax2=10 )














###############################################################################
subdir_9A_part3 <- paste(Part3_g,  "/9A-averageRows-1kb-Cmp", sep = "")
if( ! file.exists(subdir_9A_part3) ) { dir.create(subdir_9A_part3) }



MyHistogram_4A(vector2=c(reduceColumn3_Average_week0, reduceColumn3_Average_week1,  reduceColumn3_Average_week2,  reduceColumn3_Average_week4,  reduceColumn3_Average_week6, reduceColumn3_Average_week8  ),  
              sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                             rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ),  
              path2=subdir_9A_part3,    
              fileName2="3-9A-1A-WT-6Samples-stacking",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.5,  yMin2=0,  yMax2=5,  alpha2=1 )


MyHistogram_5A(vector2=c(reduceColumn3_Average_week0, reduceColumn3_Average_week1,  reduceColumn3_Average_week2,  reduceColumn3_Average_week4,  reduceColumn3_Average_week6, reduceColumn3_Average_week8  ),  
              sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                             rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ),  
              path2=subdir_9A_part3,    
              fileName2="3-9A-1B-WT-6Samples-dodge",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.5,  yMin2=0,  yMax2=5 )



MyHistogram_6A(vector2=c(reduceColumn3_Average_week0, reduceColumn3_Average_week1,  reduceColumn3_Average_week2,  reduceColumn3_Average_week4,  reduceColumn3_Average_week6, reduceColumn3_Average_week8  ),  
              sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                             rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ),  
              colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
              path2=subdir_9A_part3,    
              fileName2="3-9A-1C-WT-6Samples-density",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.5,  yMin2=0,  yMax2=5 )


MyHistogram_7A(vector2=c(reduceColumn3_Average_week0, reduceColumn3_Average_week1,  reduceColumn3_Average_week2,  reduceColumn3_Average_week4,  reduceColumn3_Average_week6, reduceColumn3_Average_week8  ),  
              sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                             rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ),  
              colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
              path2=subdir_9A_part3,    
              fileName2="3-9A-1D-WT-6Samples-density2",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=2.1,  yMin2=0,  yMax2=5 )






MyHistogram_6A(vector2=c(reduceColumn3_Average_week0, reduceColumn3_Average_H3  ),  
              sampleType2=c( rep("week0", numOfRows1),    rep("H3", numOfRows1)   ),   
              colours2=c( "week0"="red2",   "H3"="yellow2"  ),  
              path2=subdir_9A_part3,    
              fileName2="3-9A-2A-WT-2Samples-density",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.5,  yMin2=0,  yMax2=5 )


MyHistogram_7A(vector2=c(reduceColumn3_Average_week0, reduceColumn3_Average_H3 ),  
              sampleType2=c( rep("week0", numOfRows1),    rep("H3", numOfRows1)   ),  
              colours2=c( "week0"="red2",   "H3"="yellow2"  ), 
              path2=subdir_9A_part3,    
              fileName2="3-9A-2B-WT-2Samples-density2",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=2.1,  yMin2=0,  yMax2=5 )








MyHistogram_4(vector2=c(reduceColumn3_Average_week0_EEDheto,  reduceColumn3_Average_week0_EEDko, reduceColumn3_Average_week4_EEDheto,  reduceColumn3_Average_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1),  
                             rep("week4_EEDheto", numOfRows1),   rep("week4_EEDko", numOfRows1) ),  
              path2=subdir_9A_part3,    
              fileName2="3-9A-3A-CKO-4Samples-stacking",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=100,  alpha2=1)

MyHistogram_5(vector2=c(reduceColumn3_Average_week0_EEDheto,  reduceColumn3_Average_week0_EEDko, reduceColumn3_Average_week4_EEDheto,  reduceColumn3_Average_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1),  
                             rep("week4_EEDheto", numOfRows1),   rep("week4_EEDko", numOfRows1) ),  
              path2=subdir_9A_part3,    
              fileName2="3-9A-3B-CKO-4Samples-dodge",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=100 )


MyHistogram_6(vector2=c(reduceColumn3_Average_week0_EEDheto,  reduceColumn3_Average_week0_EEDko, reduceColumn3_Average_week4_EEDheto,  reduceColumn3_Average_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1),  
                             rep("week4_EEDheto", numOfRows1),   rep("week4_EEDko", numOfRows1) ),  
              colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
              path2=subdir_9A_part3,    
              fileName2="3-9A-3C-CKO-4Samples-density",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=5 )

MyHistogram_7A(vector2=c(reduceColumn3_Average_week0_EEDheto,  reduceColumn3_Average_week0_EEDko, reduceColumn3_Average_week4_EEDheto,  reduceColumn3_Average_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1),  
                             rep("week4_EEDheto", numOfRows1),   rep("week4_EEDko", numOfRows1) ),  
              colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
              path2=subdir_9A_part3,    
              fileName2="3-9A-3D-CKO-4Samples-density2",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=2.1,  yMin2=0,  yMax2=5 )









MyHistogram_6(vector2=c(reduceColumn3_Average_banding,  reduceColumn3_Average_sham ),
              sampleType2=c( rep("banding", numOfRows1),    rep("sham", numOfRows1)  ), 
              colours2=c( "banding"="blue",  "sham"="green4"  ), 
              path2=subdir_9A_part3,    
              fileName2="3-9A-4A-TAC-4Samples-density",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=5 )

MyHistogram_7(vector2=c(reduceColumn3_Average_banding,  reduceColumn3_Average_sham ),
              sampleType2=c( rep("banding", numOfRows1),    rep("sham", numOfRows1)  ), 
              colours2=c( "banding"="blue",  "sham"="green4"  ),  
              path2=subdir_9A_part3,    
              fileName2="3-9A-4B-TAC-4Samples-density2",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=0.8,  yMin2=0,  yMax2=5 )























###############################################################################
subdir_9B_part3 <- paste(Part3_g,  "/9B-averageRows-500bp-Cmp", sep = "")
if( ! file.exists(subdir_9B_part3) ) { dir.create(subdir_9B_part3) }



MyHistogram_4A(vector2=c(reduceColumn4_Average_week0, reduceColumn4_Average_week1,  reduceColumn4_Average_week2,  reduceColumn4_Average_week4,  reduceColumn4_Average_week6, reduceColumn4_Average_week8  ),  
               sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                              rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ),  
               path2=subdir_9B_part3,    
               fileName2="3-9B-1A-WT-6Samples-stacking",   
               title2=myTitle_g,    xLab2="H2BGFP signal",  
               height2=3.5,  width2=7,   xMin2=0,  xMax2=1.5,  yMin2=0,  yMax2=5,  alpha2=1 )


MyHistogram_5A(vector2=c(reduceColumn4_Average_week0, reduceColumn4_Average_week1,  reduceColumn4_Average_week2,  reduceColumn4_Average_week4,  reduceColumn4_Average_week6, reduceColumn4_Average_week8  ),  
               sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                              rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ),  
               path2=subdir_9B_part3,    
               fileName2="3-9B-1B-WT-6Samples-dodge",   
               title2=myTitle_g,    xLab2="H2BGFP signal",  
               height2=3.5,  width2=7,   xMin2=0,  xMax2=1.5,  yMin2=0,  yMax2=5 )



MyHistogram_6A(vector2=c(reduceColumn4_Average_week0, reduceColumn4_Average_week1,  reduceColumn4_Average_week2,  reduceColumn4_Average_week4,  reduceColumn4_Average_week6, reduceColumn4_Average_week8  ),  
               sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                              rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ),  
               colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
               path2=subdir_9B_part3,    
               fileName2="3-9B-1C-WT-6Samples-density",   
               title2=myTitle_g,    xLab2="H2BGFP signal",  
               height2=3.5,  width2=7,   xMin2=0,  xMax2=1.5,  yMin2=0,  yMax2=5 )


MyHistogram_7A(vector2=c(reduceColumn4_Average_week0, reduceColumn4_Average_week1,  reduceColumn4_Average_week2,  reduceColumn4_Average_week4,  reduceColumn4_Average_week6, reduceColumn4_Average_week8  ),  
               sampleType2=c( rep("week0", numOfRows1),    rep("week1", numOfRows1),  rep("week2", numOfRows1),  
                              rep("week4", numOfRows1),    rep("week6", numOfRows1),  rep("week8",   numOfRows1)  ),  
               colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
               path2=subdir_9B_part3,    
               fileName2="3-9B-1D-WT-6Samples-density2",   
               title2=myTitle_g,    xLab2="H2BGFP signal",  
               height2=3.5,  width2=7,   xMin2=0,  xMax2=2.1,  yMin2=0,  yMax2=5 )






MyHistogram_6A(vector2=c(reduceColumn4_Average_week0, reduceColumn4_Average_H3  ),  
               sampleType2=c( rep("week0", numOfRows1),    rep("H3", numOfRows1)   ),   
               colours2=c( "week0"="red2",   "H3"="yellow2"  ),  
               path2=subdir_9B_part3,    
               fileName2="3-9B-2A-WT-2Samples-density",   
               title2=myTitle_g,    xLab2="H2BGFP signal",  
               height2=3.5,  width2=7,   xMin2=0,  xMax2=1.5,  yMin2=0,  yMax2=5 )


MyHistogram_7A(vector2=c(reduceColumn4_Average_week0, reduceColumn4_Average_H3 ),  
               sampleType2=c( rep("week0", numOfRows1),    rep("H3", numOfRows1)   ),  
               colours2=c( "week0"="red2",   "H3"="yellow2"  ), 
               path2=subdir_9B_part3,    
               fileName2="3-9B-2B-WT-2Samples-density2",   
               title2=myTitle_g,    xLab2="H2BGFP signal",  
               height2=3.5,  width2=7,   xMin2=0,  xMax2=2.1,  yMin2=0,  yMax2=5 )








MyHistogram_4(vector2=c(reduceColumn4_Average_week0_EEDheto,  reduceColumn4_Average_week0_EEDko, reduceColumn4_Average_week4_EEDheto,  reduceColumn4_Average_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1),  
                             rep("week4_EEDheto", numOfRows1),   rep("week4_EEDko", numOfRows1) ),  
              path2=subdir_9B_part3,    
              fileName2="3-9B-3A-CKO-4Samples-stacking",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=100,  alpha2=1)

MyHistogram_5(vector2=c(reduceColumn4_Average_week0_EEDheto,  reduceColumn4_Average_week0_EEDko, reduceColumn4_Average_week4_EEDheto,  reduceColumn4_Average_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1),  
                             rep("week4_EEDheto", numOfRows1),   rep("week4_EEDko", numOfRows1) ),  
              path2=subdir_9B_part3,    
              fileName2="3-9B-3B-CKO-4Samples-dodge",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=100 )


MyHistogram_6(vector2=c(reduceColumn4_Average_week0_EEDheto,  reduceColumn4_Average_week0_EEDko, reduceColumn4_Average_week4_EEDheto,  reduceColumn4_Average_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1),  
                             rep("week4_EEDheto", numOfRows1),   rep("week4_EEDko", numOfRows1) ),  
              colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
              path2=subdir_9B_part3,    
              fileName2="3-9B-3C-CKO-4Samples-density",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=5 )

MyHistogram_7A(vector2=c(reduceColumn4_Average_week0_EEDheto,  reduceColumn4_Average_week0_EEDko, reduceColumn4_Average_week4_EEDheto,  reduceColumn4_Average_week4_EEDko ),  
               sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1),  
                              rep("week4_EEDheto", numOfRows1),   rep("week4_EEDko", numOfRows1) ),  
               colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
               path2=subdir_9B_part3,    
               fileName2="3-9B-3D-CKO-4Samples-density2",   
               title2=myTitle_g,    xLab2="H2BGFP signal",  
               height2=3.5,  width2=7,   xMin2=0,  xMax2=2.1,  yMin2=0,  yMax2=5 )









MyHistogram_6(vector2=c(reduceColumn4_Average_banding,  reduceColumn4_Average_sham ),
              sampleType2=c( rep("banding", numOfRows1),    rep("sham", numOfRows1)  ), 
              colours2=c( "banding"="blue",  "sham"="green4"  ), 
              path2=subdir_9B_part3,    
              fileName2="3-9B-4A-TAC-4Samples-density",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=5 )

MyHistogram_7(vector2=c(reduceColumn4_Average_banding,  reduceColumn4_Average_sham ),
              sampleType2=c( rep("banding", numOfRows1),    rep("sham", numOfRows1)  ), 
              colours2=c( "banding"="blue",  "sham"="green4"  ),  
              path2=subdir_9B_part3,    
              fileName2="3-9B-4B-TAC-4Samples-density2",   
              title2=myTitle_g,    xLab2="H2BGFP signal",  
              height2=3.5,  width2=7,   xMin2=0,  xMax2=0.8,  yMin2=0,  yMax2=5 )


























##Lowest: 0%~20%
##Low: 20%~40% 
##Medium: 40%~60% 
##High: 60%~80% 
##Highest: 80%~100%
###############################################################################
subdir_10_part3 <- paste(Part3_g,  "/10-rows5Classes-curve", sep = "")
if( ! file.exists(subdir_10_part3) ) { dir.create(subdir_10_part3) }

dim(reduceRow1_Average_H3) 
dim(reduceRow1_Average_week0) 
dim(reduceRow1_Average_week1) 
dim(reduceRow1_Average_week2) 
dim(reduceRow1_Average_week4) 
dim(reduceRow1_Average_week6) 
dim(reduceRow1_Average_week8) 

MyAverageLines_1(vector2=c(reduceRow1_Average_H3[1, ],  reduceRow1_Average_H3[2, ], reduceRow1_Average_H3[3, ],  reduceRow1_Average_H3[4, ],  reduceRow1_Average_H3[5, ]  ),   
                 numSample2=5,   
                 sampleType2=c( rep("H3_1", numOfColumns1), rep("H3_2", numOfColumns1),    rep("H3_3", numOfColumns1),  rep("H3_4", numOfColumns1),  rep("H3_5", numOfColumns1)   ), 
                 sampleRank2=c( "H3_1", "H3_2",  "H3_3",   "H3_4",  "H3_5"  ),     
                 colours2=c( "H3_1"="yellow2",  "H3_2"="red2",   "H3_3"="cyan2",    "H3_4"="green2",  "H3_5"="blue2"   ), 
                 path2=subdir_10_part3,     fileName2="3-10-1A-WT-H3-averageCurve",  
                 title2="H3 Level",     xLab2="Relative distance (kb)",    yLab2="H3 signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(reduceRow1_Average_week0[1,],  reduceRow1_Average_week1[1,],  reduceRow1_Average_week2[1,],   
                           reduceRow1_Average_week4[1,],  reduceRow1_Average_week6[1,],  reduceRow1_Average_week8[1,]  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_10_part3,     fileName2="3-10-1B-WT-6Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(reduceRow1_Average_week0[2,],  reduceRow1_Average_week1[2,],  reduceRow1_Average_week2[2,],   
                           reduceRow1_Average_week4[2,],  reduceRow1_Average_week6[2,],  reduceRow1_Average_week8[2,]  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_10_part3,     fileName2="3-10-1C-WT-6Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(reduceRow1_Average_week0[3,],  reduceRow1_Average_week1[3,],  reduceRow1_Average_week2[3,],   
                           reduceRow1_Average_week4[3,],  reduceRow1_Average_week6[3,],  reduceRow1_Average_week8[3,]  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_10_part3,     fileName2="3-10-1D-WT-6Samples-averageCurve",  
                 title2="Genes (Medium TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(reduceRow1_Average_week0[4,],  reduceRow1_Average_week1[4,],  reduceRow1_Average_week2[4,],   
                           reduceRow1_Average_week4[4,],  reduceRow1_Average_week6[4,],  reduceRow1_Average_week8[4,]  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_10_part3,     fileName2="3-10-1E-WT-6Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(reduceRow1_Average_week0[5,],  reduceRow1_Average_week1[5,],  reduceRow1_Average_week2[5,],   
                           reduceRow1_Average_week4[5,],  reduceRow1_Average_week6[5,],  reduceRow1_Average_week8[5,]  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_10_part3,     fileName2="3-10-1F-WT-6Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g  )




dim(reduceRow1_Average_week0_EEDheto)   
dim(reduceRow1_Average_week0_EEDko)  
dim(reduceRow1_Average_week4_EEDheto) 
dim(reduceRow1_Average_week4_EEDko)

MyAverageLines_1(vector2=c(reduceRow1_Average_week0_EEDheto[1,],  reduceRow1_Average_week0_EEDko[1,],      
                           reduceRow1_Average_week4_EEDheto[1,],  reduceRow1_Average_week4_EEDko[1,]   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_10_part3,     fileName2="3-10-2A-CKO-4Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(reduceRow1_Average_week0_EEDheto[2,],  reduceRow1_Average_week0_EEDko[2,],      
                           reduceRow1_Average_week4_EEDheto[2,],  reduceRow1_Average_week4_EEDko[2,]   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_10_part3,     fileName2="3-10-2B-CKO-4Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(reduceRow1_Average_week0_EEDheto[3,],  reduceRow1_Average_week0_EEDko[3,],      
                           reduceRow1_Average_week4_EEDheto[3,],  reduceRow1_Average_week4_EEDko[3,]   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_10_part3,     fileName2="3-10-2C-CKO-4Samples-averageCurve",  
                 title2="Genes (Medium TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(reduceRow1_Average_week0_EEDheto[4,],  reduceRow1_Average_week0_EEDko[4,],      
                           reduceRow1_Average_week4_EEDheto[4,],  reduceRow1_Average_week4_EEDko[4,]   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_10_part3,     fileName2="3-10-2D-CKO-4Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(reduceRow1_Average_week0_EEDheto[5,],  reduceRow1_Average_week0_EEDko[5,],      
                           reduceRow1_Average_week4_EEDheto[5,],  reduceRow1_Average_week4_EEDko[5,]   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_10_part3,     fileName2="3-10-2E-CKO-4Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g  )





dim(reduceRow1_Average_banding)  
dim(reduceRow1_Average_sham)

MyAverageLines_1(vector2=c(reduceRow1_Average_banding[1,],  reduceRow1_Average_sham[1,]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_10_part3,     fileName2="3-10-3A-TAC-2Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(reduceRow1_Average_banding[2,],  reduceRow1_Average_sham[2,]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_10_part3,     fileName2="3-10-3B-TAC-2Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(reduceRow1_Average_banding[3,],  reduceRow1_Average_sham[3,]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_10_part3,     fileName2="3-10-3C-TAC-2Samples-averageCurve",  
                 title2="Genes (Medium TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(reduceRow1_Average_banding[4,],  reduceRow1_Average_sham[4,]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_10_part3,     fileName2="3-10-3D-TAC-2Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(reduceRow1_Average_banding[5,],  reduceRow1_Average_sham[5,]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_10_part3,     fileName2="3-10-3E-TAC-2Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g  )







SOME_COLUMNS   <-  c(200:300)
dim(Average_H3)
dim(Average_week0)
numRows_oneClass <- floor( nrow(Average_H3)/5 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
INDEX_4 <- seq(from = numRows_oneClass*3+1, to = numRows_oneClass*4, by =1 )
INDEX_5 <- seq(from = numRows_oneClass*4+1, to = numRows_oneClass*5, by =1 )
nrow(Average_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) + length(INDEX_4) + length(INDEX_5))

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_week0[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week0[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week0[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week0[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week0[INDEX_5, SOME_COLUMNS]),   
                            rowMeans(Average_week1[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week1[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week1[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week1[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week1[INDEX_5, SOME_COLUMNS]), 
                            rowMeans(Average_week2[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week2[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week2[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week2[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week2[INDEX_5, SOME_COLUMNS]), 
                            rowMeans(Average_week4[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week4[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week4[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week4[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week4[INDEX_5, SOME_COLUMNS]), 
                            rowMeans(Average_week6[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week6[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week6[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week6[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week6[INDEX_5, SOME_COLUMNS]), 
                            rowMeans(Average_week8[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week8[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week8[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week8[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week8[INDEX_5, SOME_COLUMNS])
                           ),   
                  sampleType2=c( rep("week0_1", length(INDEX_1)),    rep("week0_2", length(INDEX_2)),  rep("week0_3", length(INDEX_3)),     rep("week0_4", length(INDEX_4)),  rep("week0_5", length(INDEX_5)) ,    
                                 rep("week1_1", length(INDEX_1)),    rep("week1_2", length(INDEX_2)),  rep("week1_3", length(INDEX_3)),     rep("week1_4", length(INDEX_4)),  rep("week1_5", length(INDEX_5)) , 
                                 rep("week2_1", length(INDEX_1)),    rep("week2_2", length(INDEX_2)),  rep("week2_3", length(INDEX_3)),     rep("week2_4", length(INDEX_4)),  rep("week2_5", length(INDEX_5)) , 
                                 rep("week4_1", length(INDEX_1)),    rep("week4_2", length(INDEX_2)),  rep("week4_3", length(INDEX_3)),     rep("week4_4", length(INDEX_4)),  rep("week4_5", length(INDEX_5)) , 
                                 rep("week6_1", length(INDEX_1)),    rep("week6_2", length(INDEX_2)),  rep("week6_3", length(INDEX_3)),     rep("week6_4", length(INDEX_4)),  rep("week6_5", length(INDEX_5)) ,
                                 rep("week8_1", length(INDEX_1)),    rep("week8_2", length(INDEX_2)),  rep("week8_3", length(INDEX_3)),     rep("week8_4", length(INDEX_4)),  rep("week8_5", length(INDEX_5))
                                 ), 
                  sampleRank2=c( "week0_1",  "week0_2",   "week0_3",  "week0_4",  "week0_5",  
                                 "week1_1",  "week1_2",   "week1_3",  "week1_4",  "week1_5",  
                                 "week2_1",  "week2_2",   "week2_3",  "week2_4",  "week2_5",  
                                 "week4_1",  "week4_2",   "week4_3",  "week4_4",  "week4_5",  
                                 "week6_1",  "week6_2",   "week6_3",  "week6_4",  "week6_5",  
                                 "week8_1",  "week8_2",   "week8_3",  "week8_4",  "week8_5"   ),     
                  colours2=c( "week0_1"="red2",     "week0_2"="red2",      "week0_3"="red2",     "week0_4"="red2",     "week0_5"="red2",  
                              "week1_1"="cyan2",    "week1_2"="cyan2",     "week1_3"="cyan2",    "week1_4"="cyan2",    "week1_5"="cyan2",  
                              "week2_1"="green2",   "week2_2"="green2",    "week2_3"="green2",   "week2_4"="green2",   "week2_5"="green2",  
                              "week4_1"="blue2",    "week4_2"="blue2",     "week4_3"="blue2",    "week4_4"="blue2",    "week4_5"="blue2",  
                              "week6_1"="purple2",  "week6_2"="purple2",   "week6_3"="purple2",  "week6_4"="purple2",  "week6_5"="purple2",  
                              "week8_1"="orange2",  "week8_2"="orange2",   "week8_3"="orange2",  "week8_4"="orange2",  "week8_5"="orange2"  
                             ), 
                  path2=subdir_10_part3,   fileName2= paste("Part3-10-4A-WT-6Samples-BoxViolin",  "1kb",   mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4,   width2=16,   Ymin2=0, Ymax2=1.5)   ## width = 1 + 30*0.5, height=5cm 






MyBoxViolinPlot_1(vector2=c(rowMeans(Average_week0[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week0[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week0[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week0[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week0[INDEX_5, SOME_COLUMNS]),   
                            rowMeans(Average_week1[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week1[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week1[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week1[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week1[INDEX_5, SOME_COLUMNS]), 
                            rowMeans(Average_week2[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week2[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week2[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week2[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week2[INDEX_5, SOME_COLUMNS]), 
                            rowMeans(Average_week4[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week4[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week4[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week4[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week4[INDEX_5, SOME_COLUMNS]), 
                            rowMeans(Average_week6[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week6[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week6[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week6[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week6[INDEX_5, SOME_COLUMNS]), 
                            rowMeans(Average_week8[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week8[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week8[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week8[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week8[INDEX_5, SOME_COLUMNS])
                            ),   
            sampleType2=c( rep("week0_1", length(INDEX_1)),    rep("week0_2", length(INDEX_2)),  rep("week0_3", length(INDEX_3)),     rep("week0_4", length(INDEX_4)),  rep("week0_5", length(INDEX_5)) ,    
               rep("week1_1", length(INDEX_1)),    rep("week1_2", length(INDEX_2)),  rep("week1_3", length(INDEX_3)),     rep("week1_4", length(INDEX_4)),  rep("week1_5", length(INDEX_5)) , 
               rep("week2_1", length(INDEX_1)),    rep("week2_2", length(INDEX_2)),  rep("week2_3", length(INDEX_3)),     rep("week2_4", length(INDEX_4)),  rep("week2_5", length(INDEX_5)) , 
               rep("week4_1", length(INDEX_1)),    rep("week4_2", length(INDEX_2)),  rep("week4_3", length(INDEX_3)),     rep("week4_4", length(INDEX_4)),  rep("week4_5", length(INDEX_5)) , 
               rep("week6_1", length(INDEX_1)),    rep("week6_2", length(INDEX_2)),  rep("week6_3", length(INDEX_3)),     rep("week6_4", length(INDEX_4)),  rep("week6_5", length(INDEX_5)) ,
               rep("week8_1", length(INDEX_1)),    rep("week8_2", length(INDEX_2)),  rep("week8_3", length(INDEX_3)),     rep("week8_4", length(INDEX_4)),  rep("week8_5", length(INDEX_5))
            ), 
          sampleRank2=c( "week0_1",  "week1_1",   "week2_1",   "week4_1",  "week6_1", "week8_1",  
                         "week0_2",  "week1_2",   "week2_2",   "week4_2",  "week6_2", "week8_2",   
                         "week0_3",  "week1_3",   "week2_3",   "week4_3",  "week6_3", "week8_3",   
                         "week0_4",  "week1_4",   "week2_4",   "week4_4",  "week6_4", "week8_4",   
                         "week0_5",  "week1_5",   "week2_5",   "week4_5",  "week6_5", "week8_5"  
                         ),     
        colours2=c( "week0_1"="red2",  "week1_1"="cyan2",   "week2_1"="green2",   "week4_1"="blue2",  "week6_1"="purple2", "week8_1"="orange2",  
                    "week0_2"="red2",  "week1_2"="cyan2",   "week2_2"="green2",   "week4_2"="blue2",  "week6_2"="purple2", "week8_2"="orange2",   
                    "week0_3"="red2",  "week1_3"="cyan2",   "week2_3"="green2",   "week4_3"="blue2",  "week6_3"="purple2", "week8_3"="orange2",   
                    "week0_4"="red2",  "week1_4"="cyan2",   "week2_4"="green2",   "week4_4"="blue2",  "week6_4"="purple2", "week8_4"="orange2",   
                    "week0_5"="red2",  "week1_5"="cyan2",   "week2_5"="green2",   "week4_5"="blue2",  "week6_5"="purple2", "week8_5"="orange2"  
                  ),  
      path2=subdir_10_part3,   fileName2= paste("Part3-10-4B-WT-6Samples-BoxViolin",  "1kb",   mySample_g, sep = "_") ,  
      title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
      height2=4,   width2=16,   Ymin2=0, Ymax2=1.5)   ## width = 1 + 30*0.5, height=5cm 





MyBoxViolinPlot_1(vector2=c(rowMeans(Average_week0_EEDheto[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week0_EEDheto[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week0_EEDheto[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week0_EEDheto[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week0_EEDheto[INDEX_5, SOME_COLUMNS]),   
                            rowMeans(Average_week0_EEDko[INDEX_1, SOME_COLUMNS]),    rowMeans(Average_week0_EEDko[INDEX_2, SOME_COLUMNS]),    rowMeans(Average_week0_EEDko[INDEX_3, SOME_COLUMNS]),   rowMeans(Average_week0_EEDko[INDEX_4, SOME_COLUMNS]),    rowMeans(Average_week0_EEDko[INDEX_5, SOME_COLUMNS]), 
                            rowMeans(Average_week4_EEDheto[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week4_EEDheto[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week4_EEDheto[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week4_EEDheto[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week4_EEDheto[INDEX_5, SOME_COLUMNS]),   
                            rowMeans(Average_week4_EEDko[INDEX_1, SOME_COLUMNS]),    rowMeans(Average_week4_EEDko[INDEX_2, SOME_COLUMNS]),    rowMeans(Average_week4_EEDko[INDEX_3, SOME_COLUMNS]),   rowMeans(Average_week4_EEDko[INDEX_4, SOME_COLUMNS]),    rowMeans(Average_week4_EEDko[INDEX_5, SOME_COLUMNS]) 
                            ),   
                  sampleType2=c( rep("week0_EEDheto_1", length(INDEX_1)),    rep("week0_EEDheto_2", length(INDEX_2)),  rep("week0_EEDheto_3", length(INDEX_3)),     rep("week0_EEDheto_4", length(INDEX_4)),  rep("week0_EEDheto_5", length(INDEX_5)) , 
                                 rep("week0_EEDko_1", length(INDEX_1)),      rep("week0_EEDko_2", length(INDEX_2)),    rep("week0_EEDko_3", length(INDEX_3)),       rep("week0_EEDko_4", length(INDEX_4)),    rep("week0_EEDko_5", length(INDEX_5)) ,
                                 rep("week4_EEDheto_1", length(INDEX_1)),    rep("week4_EEDheto_2", length(INDEX_2)),  rep("week4_EEDheto_3", length(INDEX_3)),     rep("week4_EEDheto_4", length(INDEX_4)),  rep("week4_EEDheto_5", length(INDEX_5)) , 
                                 rep("week4_EEDko_1", length(INDEX_1)),      rep("week4_EEDko_2", length(INDEX_2)),    rep("week4_EEDko_3", length(INDEX_3)),       rep("week4_EEDko_4", length(INDEX_4)),    rep("week4_EEDko_5", length(INDEX_5)) 
                                 ), 
                  sampleRank2=c( "week0_EEDheto_1",  "week0_EEDheto_2",   "week0_EEDheto_3",  "week0_EEDheto_4",   "week0_EEDheto_5",
                                 "week0_EEDko_1",    "week0_EEDko_2",     "week0_EEDko_3",    "week0_EEDko_4",     "week0_EEDko_5",
                                 "week4_EEDheto_1",  "week4_EEDheto_2",   "week4_EEDheto_3",  "week4_EEDheto_4",   "week4_EEDheto_5",
                                 "week4_EEDko_1",    "week4_EEDko_2",     "week4_EEDko_3",    "week4_EEDko_4",     "week4_EEDko_5"
                                 ),  
                  colours2=c(  "week0_EEDheto_1"="red",     "week0_EEDheto_2"="red",        "week0_EEDheto_3"="red",      "week0_EEDheto_4"="red",       "week0_EEDheto_5"="red",
                               "week0_EEDko_1"="red4",       "week0_EEDko_2"="red4",        "week0_EEDko_3"="red4",       "week0_EEDko_4"="red4",        "week0_EEDko_5"="red4",
                               "week4_EEDheto_1"="skyblue",  "week4_EEDheto_2"="skyblue",   "week4_EEDheto_3"="skyblue",  "week4_EEDheto_4"="skyblue",   "week4_EEDheto_5"="skyblue",
                               "week4_EEDko_1"="blue",       "week4_EEDko_2"="blue",        "week4_EEDko_3"="blue",       "week4_EEDko_4"="blue",        "week4_EEDko_5"="blue"
                            ),  
                  path2=subdir_10_part3,  fileName2= paste("Part3-10-4C-EEDko-4Samples-BoxViolin", "1kb",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=4.4,   width2=11, Ymin2=0, Ymax2=1.0 )   ## width = 1 + 20*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_week0_EEDheto[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week0_EEDheto[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week0_EEDheto[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week0_EEDheto[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week0_EEDheto[INDEX_5, SOME_COLUMNS]),   
                            rowMeans(Average_week0_EEDko[INDEX_1, SOME_COLUMNS]),    rowMeans(Average_week0_EEDko[INDEX_2, SOME_COLUMNS]),    rowMeans(Average_week0_EEDko[INDEX_3, SOME_COLUMNS]),   rowMeans(Average_week0_EEDko[INDEX_4, SOME_COLUMNS]),    rowMeans(Average_week0_EEDko[INDEX_5, SOME_COLUMNS]), 
                            rowMeans(Average_week4_EEDheto[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week4_EEDheto[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week4_EEDheto[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week4_EEDheto[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_week4_EEDheto[INDEX_5, SOME_COLUMNS]),   
                            rowMeans(Average_week4_EEDko[INDEX_1, SOME_COLUMNS]),    rowMeans(Average_week4_EEDko[INDEX_2, SOME_COLUMNS]),    rowMeans(Average_week4_EEDko[INDEX_3, SOME_COLUMNS]),   rowMeans(Average_week4_EEDko[INDEX_4, SOME_COLUMNS]),    rowMeans(Average_week4_EEDko[INDEX_5, SOME_COLUMNS]) 
                           ),   
                  sampleType2=c( rep("week0_EEDheto_1", length(INDEX_1)),    rep("week0_EEDheto_2", length(INDEX_2)),  rep("week0_EEDheto_3", length(INDEX_3)),     rep("week0_EEDheto_4", length(INDEX_4)),  rep("week0_EEDheto_5", length(INDEX_5)) , 
                                 rep("week0_EEDko_1", length(INDEX_1)),      rep("week0_EEDko_2", length(INDEX_2)),    rep("week0_EEDko_3", length(INDEX_3)),       rep("week0_EEDko_4", length(INDEX_4)),    rep("week0_EEDko_5", length(INDEX_5)) ,
                                 rep("week4_EEDheto_1", length(INDEX_1)),    rep("week4_EEDheto_2", length(INDEX_2)),  rep("week4_EEDheto_3", length(INDEX_3)),     rep("week4_EEDheto_4", length(INDEX_4)),  rep("week4_EEDheto_5", length(INDEX_5)) , 
                                 rep("week4_EEDko_1", length(INDEX_1)),      rep("week4_EEDko_2", length(INDEX_2)),    rep("week4_EEDko_3", length(INDEX_3)),       rep("week4_EEDko_4", length(INDEX_4)),    rep("week4_EEDko_5", length(INDEX_5)) 
                  ), 
                  sampleRank2=c( "week0_EEDheto_1",  "week0_EEDko_1",    "week4_EEDheto_1",  "week4_EEDko_1", 
                                 "week0_EEDheto_2",  "week0_EEDko_2",    "week4_EEDheto_2",  "week4_EEDko_2", 
                                 "week0_EEDheto_3",  "week0_EEDko_3",    "week4_EEDheto_3",  "week4_EEDko_3", 
                                 "week0_EEDheto_4",  "week0_EEDko_4",    "week4_EEDheto_4",  "week4_EEDko_4", 
                                 "week0_EEDheto_5",  "week0_EEDko_5",    "week4_EEDheto_5",  "week4_EEDko_5"
                  ),  
                  colours2=c(  "week0_EEDheto_1"="red",  "week0_EEDko_1"="red4",    "week4_EEDheto_1"="skyblue",  "week4_EEDko_1"="blue", 
                               "week0_EEDheto_2"="red",  "week0_EEDko_2"="red4",    "week4_EEDheto_2"="skyblue",  "week4_EEDko_2"="blue", 
                               "week0_EEDheto_3"="red",  "week0_EEDko_3"="red4",    "week4_EEDheto_3"="skyblue",  "week4_EEDko_3"="blue", 
                               "week0_EEDheto_4"="red",  "week0_EEDko_4"="red4",    "week4_EEDheto_4"="skyblue",  "week4_EEDko_4"="blue", 
                               "week0_EEDheto_5"="red",  "week0_EEDko_5"="red4",    "week4_EEDheto_5"="skyblue",  "week4_EEDko_5"="blue"
                  ),  
                  path2=subdir_10_part3,  fileName2= paste("Part3-10-4D-EEDko-4Samples-BoxViolin", "1kb",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=4.4,   width2=11, Ymin2=0, Ymax2=1.0 )   ## width = 1 + 20*0.5, height=5cm 




MyBoxViolinPlot_1(vector2=c(rowMeans(Average_banding[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_banding[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_banding[INDEX_3, SOME_COLUMNS]), rowMeans(Average_banding[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_banding[INDEX_5, SOME_COLUMNS]),
                            rowMeans(Average_sham[INDEX_1, SOME_COLUMNS]),     rowMeans(Average_sham[INDEX_2, SOME_COLUMNS]),     rowMeans(Average_sham[INDEX_3, SOME_COLUMNS]),    rowMeans(Average_sham[INDEX_4, SOME_COLUMNS]),     rowMeans(Average_sham[INDEX_5, SOME_COLUMNS])
                            ),   
                  sampleType2=c( rep("banding_1", length(INDEX_1)),    rep("banding_2", length(INDEX_2)),  rep("banding_3", length(INDEX_3)),     rep("banding_4", length(INDEX_4)),  rep("banding_5", length(INDEX_5)) ,
                                 rep("sham_1", length(INDEX_1)),       rep("sham_2", length(INDEX_2)),     rep("sham_3", length(INDEX_3)),        rep("sham_4", length(INDEX_4)),     rep("sham_5", length(INDEX_5)) 
                                 ), 
                  sampleRank2=c( "banding_1",  "banding_2",   "banding_3",  "banding_4",   "banding_5",  
                                 "sham_1",     "sham_2",      "sham_3",     "sham_4",      "sham_5"
                                 ),     
                  colours2=c( "banding_1"="blue",    "banding_2"="blue",     "banding_3"="blue",    "banding_4"="blue",     "banding_5"="blue",  
                              "sham_1"="green4",     "sham_2"="green4",      "sham_3"="green4",     "sham_4"="green4",      "sham_5"="green4"
                             ), 
                  path2=subdir_10_part3,  fileName2= paste("Part3-10-4E-TAC-2Samples-BoxViolin", "1kb",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=4.2,   width2=6, Ymin2=0, Ymax2=0.8 )   ## width = 1 + 10*0.5, height=5cm 





MyBoxViolinPlot_1(vector2=c(rowMeans(Average_banding[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_banding[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_banding[INDEX_3, SOME_COLUMNS]), rowMeans(Average_banding[INDEX_4, SOME_COLUMNS]),  rowMeans(Average_banding[INDEX_5, SOME_COLUMNS]),
                            rowMeans(Average_sham[INDEX_1, SOME_COLUMNS]),     rowMeans(Average_sham[INDEX_2, SOME_COLUMNS]),     rowMeans(Average_sham[INDEX_3, SOME_COLUMNS]),    rowMeans(Average_sham[INDEX_4, SOME_COLUMNS]),     rowMeans(Average_sham[INDEX_5, SOME_COLUMNS])
                            ),   
                  sampleType2=c( rep("banding_1", length(INDEX_1)),    rep("banding_2", length(INDEX_2)),  rep("banding_3", length(INDEX_3)),     rep("banding_4", length(INDEX_4)),  rep("banding_5", length(INDEX_5)) ,
                                 rep("sham_1", length(INDEX_1)),       rep("sham_2", length(INDEX_2)),     rep("sham_3", length(INDEX_3)),        rep("sham_4", length(INDEX_4)),     rep("sham_5", length(INDEX_5)) 
                  ), 
                  sampleRank2=c( "banding_1",   "sham_1", 
                                 "banding_2",   "sham_2", 
                                 "banding_3",   "sham_3",
                                 "banding_4",   "sham_4", 
                                 "banding_5",   "sham_5"
                  ),     
                  colours2=c( "banding_1"="blue",   "sham_1"="green4", 
                              "banding_2"="blue",   "sham_2"="green4", 
                              "banding_3"="blue",   "sham_3"="green4",
                              "banding_4"="blue",   "sham_4"="green4", 
                              "banding_5"="blue",   "sham_5"="green4"
                  ), 
                  path2=subdir_10_part3,  fileName2= paste("Part3-10-4F-TAC-2Samples-BoxViolin", "1kb",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=4.2,   width2=6, Ymin2=0, Ymax2=0.8 )   ## width = 1 + 10*0.5, height=5cm 















##Lowest: 0%~25%
##Low: 25%~50% 
##High: 50%~75% 
##Highest: 75%~100%
###############################################################################
subdir_11_part3 <- paste(Part3_g,  "/11-rows4Classes-curve", sep = "")
if( ! file.exists(subdir_11_part3) ) { dir.create(subdir_11_part3) }

dim(reduceRow2_Average_H3) 
dim(reduceRow2_Average_week0) 
dim(reduceRow2_Average_week1) 
dim(reduceRow2_Average_week2) 
dim(reduceRow2_Average_week4) 
dim(reduceRow2_Average_week6) 
dim(reduceRow2_Average_week8) 

MyAverageLines_1(vector2=c(reduceRow2_Average_H3[1, ],  reduceRow2_Average_H3[2, ], reduceRow2_Average_H3[3, ],  reduceRow2_Average_H3[4, ] ),   
                 numSample2=4,   
                 sampleType2=c( rep("H3_1", numOfColumns1), rep("H3_2", numOfColumns1),    rep("H3_3", numOfColumns1),  rep("H3_4", numOfColumns1)  ), 
                 sampleRank2=c( "H3_1", "H3_2",  "H3_3",   "H3_4"   ),     
                 colours2=c( "H3_1"="yellow2",  "H3_2"="red2",   "H3_3"="cyan2",    "H3_4"="green2"    ), 
                 path2=subdir_11_part3,     fileName2="3-11-1A-WT-H3-averageCurve",  
                 title2="H3 Level",     xLab2="Relative distance (kb)",    yLab2="H3 signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(reduceRow2_Average_week0[1,],  reduceRow2_Average_week1[1,],  reduceRow2_Average_week2[1,],   
                           reduceRow2_Average_week4[1,],  reduceRow2_Average_week6[1,],  reduceRow2_Average_week8[1,]  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_11_part3,     fileName2="3-11-1B-WT-6Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(reduceRow2_Average_week0[2,],  reduceRow2_Average_week1[2,],  reduceRow2_Average_week2[2,],   
                           reduceRow2_Average_week4[2,],  reduceRow2_Average_week6[2,],  reduceRow2_Average_week8[2,]  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_11_part3,     fileName2="3-11-1C-WT-6Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(reduceRow2_Average_week0[3,],  reduceRow2_Average_week1[3,],  reduceRow2_Average_week2[3,],   
                           reduceRow2_Average_week4[3,],  reduceRow2_Average_week6[3,],  reduceRow2_Average_week8[3,]  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_11_part3,     fileName2="3-11-1D-WT-6Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(reduceRow2_Average_week0[4,],  reduceRow2_Average_week1[4,],  reduceRow2_Average_week2[4,],   
                           reduceRow2_Average_week4[4,],  reduceRow2_Average_week6[4,],  reduceRow2_Average_week8[4,]  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_11_part3,     fileName2="3-11-1E-WT-6Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )




dim(reduceRow2_Average_week0_EEDheto)   
dim(reduceRow2_Average_week0_EEDko)  
dim(reduceRow2_Average_week4_EEDheto) 
dim(reduceRow2_Average_week4_EEDko)

MyAverageLines_1(vector2=c(reduceRow2_Average_week0_EEDheto[1,],  reduceRow2_Average_week0_EEDko[1,],      
                           reduceRow2_Average_week4_EEDheto[1,],  reduceRow2_Average_week4_EEDko[1,]   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_11_part3,     fileName2="3-11-2A-CKO-4Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.6,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(reduceRow2_Average_week0_EEDheto[2,],  reduceRow2_Average_week0_EEDko[2,],      
                           reduceRow2_Average_week4_EEDheto[2,],  reduceRow2_Average_week4_EEDko[2,]   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_11_part3,     fileName2="3-11-2B-CKO-4Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.6,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(reduceRow2_Average_week0_EEDheto[3,],  reduceRow2_Average_week0_EEDko[3,],      
                           reduceRow2_Average_week4_EEDheto[3,],  reduceRow2_Average_week4_EEDko[3,]   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_11_part3,     fileName2="3-11-2C-CKO-4Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.6,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(reduceRow2_Average_week0_EEDheto[4,],  reduceRow2_Average_week0_EEDko[4,],      
                           reduceRow2_Average_week4_EEDheto[4,],  reduceRow2_Average_week4_EEDko[4,]   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_11_part3,     fileName2="3-11-2D-CKO-4Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.6,    height2=3.3,   width2=6.05 , center2=myCenter_g )





dim(reduceRow2_Average_banding)  
dim(reduceRow2_Average_sham)

MyAverageLines_1(vector2=c(reduceRow2_Average_banding[1,],  reduceRow2_Average_sham[1,]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_11_part3,     fileName2="3-11-3A-TAC-2Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(reduceRow2_Average_banding[2,],  reduceRow2_Average_sham[2,]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_11_part3,     fileName2="3-11-3B-TAC-2Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(reduceRow2_Average_banding[3,],  reduceRow2_Average_sham[3,]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_11_part3,     fileName2="3-11-3C-TAC-2Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(reduceRow2_Average_banding[4,],  reduceRow2_Average_sham[4,]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_11_part3,     fileName2="3-11-3D-TAC-2Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )











SOME_COLUMNS   <-  c(200:300)
dim(Average_H3)
dim(Average_week0)
numRows_oneClass <- floor( nrow(Average_H3)/4 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
INDEX_4 <- seq(from = numRows_oneClass*3+1, to = numRows_oneClass*4, by =1 )
nrow(Average_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) + length(INDEX_4)  )

MyBoxViolinPlot_1(
                 vector2=c( rowMeans(Average_week0[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week0[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week0[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week0[INDEX_4, SOME_COLUMNS]),      
                            rowMeans(Average_week1[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week1[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week1[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week1[INDEX_4, SOME_COLUMNS]),    
                            rowMeans(Average_week2[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week2[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week2[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week2[INDEX_4, SOME_COLUMNS]),  
                            rowMeans(Average_week4[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week4[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week4[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week4[INDEX_4, SOME_COLUMNS]),   
                            rowMeans(Average_week6[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week6[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week6[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week6[INDEX_4, SOME_COLUMNS]),   
                            rowMeans(Average_week8[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week8[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week8[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week8[INDEX_4, SOME_COLUMNS])  
                 ),   
                 sampleType2=c( rep("week0_1", length(INDEX_1)),    rep("week0_2", length(INDEX_2)),  rep("week0_3", length(INDEX_3)),     rep("week0_4", length(INDEX_4))   ,    
                                rep("week1_1", length(INDEX_1)),    rep("week1_2", length(INDEX_2)),  rep("week1_3", length(INDEX_3)),     rep("week1_4", length(INDEX_4))   , 
                                rep("week2_1", length(INDEX_1)),    rep("week2_2", length(INDEX_2)),  rep("week2_3", length(INDEX_3)),     rep("week2_4", length(INDEX_4))   , 
                                rep("week4_1", length(INDEX_1)),    rep("week4_2", length(INDEX_2)),  rep("week4_3", length(INDEX_3)),     rep("week4_4", length(INDEX_4))   , 
                                rep("week6_1", length(INDEX_1)),    rep("week6_2", length(INDEX_2)),  rep("week6_3", length(INDEX_3)),     rep("week6_4", length(INDEX_4))   ,
                                rep("week8_1", length(INDEX_1)),    rep("week8_2", length(INDEX_2)),  rep("week8_3", length(INDEX_3)),     rep("week8_4", length(INDEX_4))  
                 ), 
                 sampleRank2=c( "week0_1",  "week0_2",   "week0_3",  "week0_4",     
                                "week1_1",  "week1_2",   "week1_3",  "week1_4",     
                                "week2_1",  "week2_2",   "week2_3",  "week2_4",    
                                "week4_1",  "week4_2",   "week4_3",  "week4_4",     
                                "week6_1",  "week6_2",   "week6_3",  "week6_4",   
                                "week8_1",  "week8_2",   "week8_3",  "week8_4"    
                                ),      
                 colours2=c( "week0_1"="red2",     "week0_2"="red2",      "week0_3"="red2",     "week0_4"="red2",     
                             "week1_1"="cyan2",    "week1_2"="cyan2",     "week1_3"="cyan2",    "week1_4"="cyan2",      
                             "week2_1"="green2",   "week2_2"="green2",    "week2_3"="green2",   "week2_4"="green2",    
                             "week4_1"="blue2",    "week4_2"="blue2",     "week4_3"="blue2",    "week4_4"="blue2",     
                             "week6_1"="purple2",  "week6_2"="purple2",   "week6_3"="purple2",  "week6_4"="purple2",   
                             "week8_1"="orange2",  "week8_2"="orange2",   "week8_3"="orange2",  "week8_4"="orange2"  
                 ), 
                 path2=subdir_11_part3,   
                 fileName2= paste("Part3-11-4A-WT-6Samples-BoxViolin",  "1kb",   mySample_g, sep = "_") ,  
                 title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                 height2=4,   width2=13,   Ymin2=0, Ymax2=1.5)   ## width = 1 + 24*0.5, height=5cm 






MyBoxViolinPlot_1(
                  vector2=c(rowMeans(Average_week0[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week0[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week0[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week0[INDEX_4, SOME_COLUMNS]),   
                            rowMeans(Average_week1[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week1[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week1[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week1[INDEX_4, SOME_COLUMNS]),  
                            rowMeans(Average_week2[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week2[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week2[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week2[INDEX_4, SOME_COLUMNS]),   
                            rowMeans(Average_week4[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week4[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week4[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week4[INDEX_4, SOME_COLUMNS]),  
                            rowMeans(Average_week6[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week6[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week6[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week6[INDEX_4, SOME_COLUMNS]),   
                            rowMeans(Average_week8[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week8[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week8[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week8[INDEX_4, SOME_COLUMNS]) 
                  ),   
                  sampleType2=c( rep("week0_1", length(INDEX_1)),    rep("week0_2", length(INDEX_2)),  rep("week0_3", length(INDEX_3)),     rep("week0_4", length(INDEX_4)),      
                                 rep("week1_1", length(INDEX_1)),    rep("week1_2", length(INDEX_2)),  rep("week1_3", length(INDEX_3)),     rep("week1_4", length(INDEX_4)),   
                                 rep("week2_1", length(INDEX_1)),    rep("week2_2", length(INDEX_2)),  rep("week2_3", length(INDEX_3)),     rep("week2_4", length(INDEX_4)),   
                                 rep("week4_1", length(INDEX_1)),    rep("week4_2", length(INDEX_2)),  rep("week4_3", length(INDEX_3)),     rep("week4_4", length(INDEX_4)),   
                                 rep("week6_1", length(INDEX_1)),    rep("week6_2", length(INDEX_2)),  rep("week6_3", length(INDEX_3)),     rep("week6_4", length(INDEX_4)),   
                                 rep("week8_1", length(INDEX_1)),    rep("week8_2", length(INDEX_2)),  rep("week8_3", length(INDEX_3)),     rep("week8_4", length(INDEX_4)) 
                  ), 
                  sampleRank2=c( "week0_1",  "week1_1",   "week2_1",   "week4_1",  "week6_1", "week8_1",  
                                 "week0_2",  "week1_2",   "week2_2",   "week4_2",  "week6_2", "week8_2",   
                                 "week0_3",  "week1_3",   "week2_3",   "week4_3",  "week6_3", "week8_3",   
                                 "week0_4",  "week1_4",   "week2_4",   "week4_4",  "week6_4", "week8_4"                    
                                 ),     
                  colours2=c( "week0_1"="red2",  "week1_1"="cyan2",   "week2_1"="green2",   "week4_1"="blue2",  "week6_1"="purple2", "week8_1"="orange2",  
                              "week0_2"="red2",  "week1_2"="cyan2",   "week2_2"="green2",   "week4_2"="blue2",  "week6_2"="purple2", "week8_2"="orange2",   
                              "week0_3"="red2",  "week1_3"="cyan2",   "week2_3"="green2",   "week4_3"="blue2",  "week6_3"="purple2", "week8_3"="orange2",   
                              "week0_4"="red2",  "week1_4"="cyan2",   "week2_4"="green2",   "week4_4"="blue2",  "week6_4"="purple2", "week8_4"="orange2"    
                  ),  
                  path2=subdir_11_part3,   
                  fileName2= paste("Part3-11-4B-WT-6Samples-BoxViolin",  "1kb",   mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4,   width2=13,   Ymin2=0, Ymax2=1.5)   ## width = 1 + 30*0.5, height=5cm 





MyBoxViolinPlot_1(vector2=c(rowMeans(Average_week0_EEDheto[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week0_EEDheto[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week0_EEDheto[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week0_EEDheto[INDEX_4, SOME_COLUMNS]),     
                            rowMeans(Average_week0_EEDko[INDEX_1, SOME_COLUMNS]),    rowMeans(Average_week0_EEDko[INDEX_2, SOME_COLUMNS]),    rowMeans(Average_week0_EEDko[INDEX_3, SOME_COLUMNS]),   rowMeans(Average_week0_EEDko[INDEX_4, SOME_COLUMNS]),    
                            rowMeans(Average_week4_EEDheto[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week4_EEDheto[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week4_EEDheto[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week4_EEDheto[INDEX_4, SOME_COLUMNS]),     
                            rowMeans(Average_week4_EEDko[INDEX_1, SOME_COLUMNS]),    rowMeans(Average_week4_EEDko[INDEX_2, SOME_COLUMNS]),    rowMeans(Average_week4_EEDko[INDEX_3, SOME_COLUMNS]),   rowMeans(Average_week4_EEDko[INDEX_4, SOME_COLUMNS])      
                           ),   
                  sampleType2=c( rep("week0_EEDheto_1", length(INDEX_1)),    rep("week0_EEDheto_2", length(INDEX_2)),  rep("week0_EEDheto_3", length(INDEX_3)),     rep("week0_EEDheto_4", length(INDEX_4)),   
                                 rep("week0_EEDko_1", length(INDEX_1)),      rep("week0_EEDko_2", length(INDEX_2)),    rep("week0_EEDko_3", length(INDEX_3)),       rep("week0_EEDko_4", length(INDEX_4)),    
                                 rep("week4_EEDheto_1", length(INDEX_1)),    rep("week4_EEDheto_2", length(INDEX_2)),  rep("week4_EEDheto_3", length(INDEX_3)),     rep("week4_EEDheto_4", length(INDEX_4)),  
                                 rep("week4_EEDko_1", length(INDEX_1)),      rep("week4_EEDko_2", length(INDEX_2)),    rep("week4_EEDko_3", length(INDEX_3)),       rep("week4_EEDko_4", length(INDEX_4)) 
                  ), 
                  sampleRank2=c( "week0_EEDheto_1",  "week0_EEDheto_2",   "week0_EEDheto_3",  "week0_EEDheto_4" ,
                                 "week0_EEDko_1",    "week0_EEDko_2",     "week0_EEDko_3",    "week0_EEDko_4", 
                                 "week4_EEDheto_1",  "week4_EEDheto_2",   "week4_EEDheto_3",  "week4_EEDheto_4", 
                                 "week4_EEDko_1",    "week4_EEDko_2",     "week4_EEDko_3",    "week4_EEDko_4" 
                  ),  
                  colours2=c(  "week0_EEDheto_1"="red",     "week0_EEDheto_2"="red",        "week0_EEDheto_3"="red",      "week0_EEDheto_4"="red",       
                               "week0_EEDko_1"="red4",       "week0_EEDko_2"="red4",        "week0_EEDko_3"="red4",       "week0_EEDko_4"="red4",         
                               "week4_EEDheto_1"="skyblue",  "week4_EEDheto_2"="skyblue",   "week4_EEDheto_3"="skyblue",  "week4_EEDheto_4"="skyblue",   
                               "week4_EEDko_1"="blue",       "week4_EEDko_2"="blue",        "week4_EEDko_3"="blue",       "week4_EEDko_4"="blue" 
                  ),  
                  path2=subdir_11_part3,  fileName2= paste("Part3-11-4C-EEDko-4Samples-BoxViolin", "1kb",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=4.4,   width2=9, Ymin2=0, Ymax2=1.0 )   ## width = 1 + 20*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_week0_EEDheto[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week0_EEDheto[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week0_EEDheto[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week0_EEDheto[INDEX_4, SOME_COLUMNS]),     
                            rowMeans(Average_week0_EEDko[INDEX_1, SOME_COLUMNS]),    rowMeans(Average_week0_EEDko[INDEX_2, SOME_COLUMNS]),    rowMeans(Average_week0_EEDko[INDEX_3, SOME_COLUMNS]),   rowMeans(Average_week0_EEDko[INDEX_4, SOME_COLUMNS]),     
                            rowMeans(Average_week4_EEDheto[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week4_EEDheto[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week4_EEDheto[INDEX_3, SOME_COLUMNS]), rowMeans(Average_week4_EEDheto[INDEX_4, SOME_COLUMNS]),     
                            rowMeans(Average_week4_EEDko[INDEX_1, SOME_COLUMNS]),    rowMeans(Average_week4_EEDko[INDEX_2, SOME_COLUMNS]),    rowMeans(Average_week4_EEDko[INDEX_3, SOME_COLUMNS]),   rowMeans(Average_week4_EEDko[INDEX_4, SOME_COLUMNS])     
                 ),   
sampleType2=c( rep("week0_EEDheto_1", length(INDEX_1)),    rep("week0_EEDheto_2", length(INDEX_2)),  rep("week0_EEDheto_3", length(INDEX_3)),     rep("week0_EEDheto_4", length(INDEX_4)),   
               rep("week0_EEDko_1", length(INDEX_1)),      rep("week0_EEDko_2", length(INDEX_2)),    rep("week0_EEDko_3", length(INDEX_3)),       rep("week0_EEDko_4", length(INDEX_4)),     
               rep("week4_EEDheto_1", length(INDEX_1)),    rep("week4_EEDheto_2", length(INDEX_2)),  rep("week4_EEDheto_3", length(INDEX_3)),     rep("week4_EEDheto_4", length(INDEX_4)),   
               rep("week4_EEDko_1", length(INDEX_1)),      rep("week4_EEDko_2", length(INDEX_2)),    rep("week4_EEDko_3", length(INDEX_3)),       rep("week4_EEDko_4", length(INDEX_4)) 
), 
sampleRank2=c( "week0_EEDheto_1",  "week0_EEDko_1",    "week4_EEDheto_1",  "week4_EEDko_1", 
               "week0_EEDheto_2",  "week0_EEDko_2",    "week4_EEDheto_2",  "week4_EEDko_2", 
               "week0_EEDheto_3",  "week0_EEDko_3",    "week4_EEDheto_3",  "week4_EEDko_3", 
               "week0_EEDheto_4",  "week0_EEDko_4",    "week4_EEDheto_4",  "week4_EEDko_4"
),  
colours2=c(  "week0_EEDheto_1"="red",  "week0_EEDko_1"="red4",    "week4_EEDheto_1"="skyblue",  "week4_EEDko_1"="blue", 
             "week0_EEDheto_2"="red",  "week0_EEDko_2"="red4",    "week4_EEDheto_2"="skyblue",  "week4_EEDko_2"="blue", 
             "week0_EEDheto_3"="red",  "week0_EEDko_3"="red4",    "week4_EEDheto_3"="skyblue",  "week4_EEDko_3"="blue", 
             "week0_EEDheto_4"="red",  "week0_EEDko_4"="red4",    "week4_EEDheto_4"="skyblue",  "week4_EEDko_4"="blue" 
),  
path2=subdir_11_part3,  fileName2= paste("Part3-11-4D-EEDko-4Samples-BoxViolin", "1kb",  mySample_g, sep = "_") , 
title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
height2=4.4,   width2=9, Ymin2=0, Ymax2=1.0 )   ## width = 1 + 16*0.5, height=5cm 




MyBoxViolinPlot_1(vector2=c(rowMeans(Average_banding[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_banding[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_banding[INDEX_3, SOME_COLUMNS]), rowMeans(Average_banding[INDEX_4, SOME_COLUMNS]),  
                            rowMeans(Average_sham[INDEX_1, SOME_COLUMNS]),     rowMeans(Average_sham[INDEX_2, SOME_COLUMNS]),     rowMeans(Average_sham[INDEX_3, SOME_COLUMNS]),    rowMeans(Average_sham[INDEX_4, SOME_COLUMNS]) 
),   
sampleType2=c( rep("banding_1", length(INDEX_1)),    rep("banding_2", length(INDEX_2)),  rep("banding_3", length(INDEX_3)),     rep("banding_4", length(INDEX_4)) ,
               rep("sham_1", length(INDEX_1)),       rep("sham_2", length(INDEX_2)),     rep("sham_3", length(INDEX_3)),        rep("sham_4", length(INDEX_4)) 
), 
sampleRank2=c( "banding_1",  "banding_2",   "banding_3",  "banding_4",    
               "sham_1",     "sham_2",      "sham_3",     "sham_4" 
),     
colours2=c( "banding_1"="blue",    "banding_2"="blue",     "banding_3"="blue",    "banding_4"="blue",   
            "sham_1"="green4",     "sham_2"="green4",      "sham_3"="green4",     "sham_4"="green4" 
), 
path2=subdir_11_part3,  fileName2= paste("Part3-11-4E-TAC-2Samples-BoxViolin", "1kb",  mySample_g, sep = "_") , 
title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
height2=4.2,   width2=6, Ymin2=0, Ymax2=0.8 )   ## width = 1 + 10*0.5, height=5cm 





MyBoxViolinPlot_1(vector2=c(rowMeans(Average_banding[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_banding[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_banding[INDEX_3, SOME_COLUMNS]), rowMeans(Average_banding[INDEX_4, SOME_COLUMNS]),   
                            rowMeans(Average_sham[INDEX_1, SOME_COLUMNS]),     rowMeans(Average_sham[INDEX_2, SOME_COLUMNS]),     rowMeans(Average_sham[INDEX_3, SOME_COLUMNS]),    rowMeans(Average_sham[INDEX_4, SOME_COLUMNS]) 
),   
sampleType2=c( rep("banding_1", length(INDEX_1)),    rep("banding_2", length(INDEX_2)),  rep("banding_3", length(INDEX_3)),     rep("banding_4", length(INDEX_4)),   
               rep("sham_1", length(INDEX_1)),       rep("sham_2", length(INDEX_2)),     rep("sham_3", length(INDEX_3)),        rep("sham_4", length(INDEX_4)) 
), 
sampleRank2=c( "banding_1",   "sham_1", 
               "banding_2",   "sham_2", 
               "banding_3",   "sham_3",
               "banding_4",   "sham_4"  
),     
colours2=c( "banding_1"="blue",   "sham_1"="green4", 
            "banding_2"="blue",   "sham_2"="green4", 
            "banding_3"="blue",   "sham_3"="green4",
            "banding_4"="blue",   "sham_4"="green4" 
), 
path2=subdir_11_part3,  fileName2= paste("Part3-11-4F-TAC-2Samples-BoxViolin", "1kb",  mySample_g, sep = "_") , 
title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
height2=4.2,   width2=6, Ymin2=0, Ymax2=0.8 )   ## width = 1 + 10*0.5, height=5cm 
























##Low: 0%~33% 
##Medium: 33%~66%
##High: 66%~100% 

###############################################################################
subdir_12_part3 <- paste(Part3_g,  "/12-rows3Classes-curve", sep = "")
if( ! file.exists(subdir_12_part3) ) { dir.create(subdir_12_part3) }

dim(reduceRow3_Average_H3) 
dim(reduceRow3_Average_week0) 
dim(reduceRow3_Average_week1) 
dim(reduceRow3_Average_week2) 
dim(reduceRow3_Average_week4) 
dim(reduceRow3_Average_week6) 
dim(reduceRow3_Average_week8) 

MyAverageLines_1(vector2=c(reduceRow3_Average_H3[1, ],  reduceRow3_Average_H3[2, ], reduceRow3_Average_H3[3, ]  ),   
                 numSample2=3,   
                 sampleType2=c( rep("H3_1", numOfColumns1), rep("H3_2", numOfColumns1),    rep("H3_3", numOfColumns1)   ), 
                 sampleRank2=c( "H3_1", "H3_2",  "H3_3"   ),     
                 colours2=c( "H3_1"="yellow2",  "H3_2"="red2",   "H3_3"="cyan2"     ), 
                 path2=subdir_12_part3,     fileName2="3-12-1A-WT-H3-averageCurve",  
                 title2="H3 Level",     xLab2="Relative distance (kb)",    yLab2="H3 signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0   , center2=myCenter_g)

MyAverageLines_1(vector2=c(reduceRow3_Average_week0[1,],  reduceRow3_Average_week1[1,],  reduceRow3_Average_week2[1,],   
                           reduceRow3_Average_week4[1,],  reduceRow3_Average_week6[1,],  reduceRow3_Average_week8[1,]  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_12_part3,     fileName2="3-12-1B-WT-6Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0  , center2=myCenter_g )

MyAverageLines_1(vector2=c(reduceRow3_Average_week0[2,],  reduceRow3_Average_week1[2,],  reduceRow3_Average_week2[2,],   
                           reduceRow3_Average_week4[2,],  reduceRow3_Average_week6[2,],  reduceRow3_Average_week8[2,]  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_12_part3,     fileName2="3-12-1C-WT-6Samples-averageCurve",  
                 title2="Genes (Medium TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0  , center2=myCenter_g )

MyAverageLines_1(vector2=c(reduceRow3_Average_week0[3,],  reduceRow3_Average_week1[3,],  reduceRow3_Average_week2[3,],   
                           reduceRow3_Average_week4[3,],  reduceRow3_Average_week6[3,],  reduceRow3_Average_week8[3,]  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_12_part3,     fileName2="3-12-1D-WT-6Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0  , center2=myCenter_g )



dim(reduceRow3_Average_week0_EEDheto)   
dim(reduceRow3_Average_week0_EEDko)  
dim(reduceRow3_Average_week4_EEDheto) 
dim(reduceRow3_Average_week4_EEDko)

MyAverageLines_1(vector2=c(reduceRow3_Average_week0_EEDheto[1,],  reduceRow3_Average_week0_EEDko[1,],      
                           reduceRow3_Average_week4_EEDheto[1,],  reduceRow3_Average_week4_EEDko[1,]   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_12_part3,     fileName2="3-12-2A-CKO-4Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.6,    height2=3.3,   width2=6.05  , center2=myCenter_g )

MyAverageLines_1(vector2=c(reduceRow3_Average_week0_EEDheto[2,],  reduceRow3_Average_week0_EEDko[2,],      
                           reduceRow3_Average_week4_EEDheto[2,],  reduceRow3_Average_week4_EEDko[2,]   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_12_part3,     fileName2="3-12-2B-CKO-4Samples-averageCurve",  
                 title2="Genes (Medium TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.6,    height2=3.3,   width2=6.05  , center2=myCenter_g )

MyAverageLines_1(vector2=c(reduceRow3_Average_week0_EEDheto[3,],  reduceRow3_Average_week0_EEDko[3,],      
                           reduceRow3_Average_week4_EEDheto[3,],  reduceRow3_Average_week4_EEDko[3,]   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_12_part3,     fileName2="3-12-2C-CKO-4Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.6,    height2=3.3,   width2=6.05  , center2=myCenter_g )





dim(reduceRow3_Average_banding)  
dim(reduceRow3_Average_sham)

MyAverageLines_1(vector2=c(reduceRow3_Average_banding[1,],  reduceRow3_Average_sham[1,]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_12_part3,     fileName2="3-12-3A-TAC-2Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(reduceRow3_Average_banding[2,],  reduceRow3_Average_sham[2,]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_12_part3,     fileName2="3-12-3B-TAC-2Samples-averageCurve",  
                 title2="Genes (Medium TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15  , center2=myCenter_g )

MyAverageLines_1(vector2=c(reduceRow3_Average_banding[3,],  reduceRow3_Average_sham[3,]  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_12_part3,     fileName2="3-12-3C-TAC-2Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15  , center2=myCenter_g )
















SOME_COLUMNS   <-  c(200:300)
dim(Average_H3)
dim(Average_week0)
numRows_oneClass <- floor( nrow(Average_H3)/3 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
nrow(Average_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3)   )

MyBoxViolinPlot_1(
  vector2=c( rowMeans(Average_week0[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week0[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week0[INDEX_3, SOME_COLUMNS]),       
             rowMeans(Average_week1[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week1[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week1[INDEX_3, SOME_COLUMNS]),    
             rowMeans(Average_week2[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week2[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week2[INDEX_3, SOME_COLUMNS]),   
             rowMeans(Average_week4[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week4[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week4[INDEX_3, SOME_COLUMNS]),     
             rowMeans(Average_week6[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week6[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week6[INDEX_3, SOME_COLUMNS]),    
             rowMeans(Average_week8[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week8[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week8[INDEX_3, SOME_COLUMNS])  
  ),   
  sampleType2=c( rep("week0_1", length(INDEX_1)),    rep("week0_2", length(INDEX_2)),  rep("week0_3", length(INDEX_3)),         
                 rep("week1_1", length(INDEX_1)),    rep("week1_2", length(INDEX_2)),  rep("week1_3", length(INDEX_3)),       
                 rep("week2_1", length(INDEX_1)),    rep("week2_2", length(INDEX_2)),  rep("week2_3", length(INDEX_3)),       
                 rep("week4_1", length(INDEX_1)),    rep("week4_2", length(INDEX_2)),  rep("week4_3", length(INDEX_3)),      
                 rep("week6_1", length(INDEX_1)),    rep("week6_2", length(INDEX_2)),  rep("week6_3", length(INDEX_3)),     
                 rep("week8_1", length(INDEX_1)),    rep("week8_2", length(INDEX_2)),  rep("week8_3", length(INDEX_3))   
  ), 
  sampleRank2=c( "week0_1",  "week0_2",   "week0_3",      
                 "week1_1",  "week1_2",   "week1_3",       
                 "week2_1",  "week2_2",   "week2_3",     
                 "week4_1",  "week4_2",   "week4_3",      
                 "week6_1",  "week6_2",   "week6_3",     
                 "week8_1",  "week8_2",   "week8_3"     
  ),      
  colours2=c( "week0_1"="red2",     "week0_2"="red2",      "week0_3"="red2",          
              "week1_1"="cyan2",    "week1_2"="cyan2",     "week1_3"="cyan2",          
              "week2_1"="green2",   "week2_2"="green2",    "week2_3"="green2",      
              "week4_1"="blue2",    "week4_2"="blue2",     "week4_3"="blue2",          
              "week6_1"="purple2",  "week6_2"="purple2",   "week6_3"="purple2",     
              "week8_1"="orange2",  "week8_2"="orange2",   "week8_3"="orange2"   
  ), 
  path2=subdir_12_part3,   
  fileName2= paste("Part3-12-4A-WT-6Samples-BoxViolin",  "1kb",   mySample_g, sep = "_") ,  
  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
  height2=4,   width2=10,   Ymin2=0, Ymax2=1.5)   ## width = 1 + 18*0.5, height=5cm 






MyBoxViolinPlot_1(
  vector2=c(rowMeans(Average_week0[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week0[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week0[INDEX_3, SOME_COLUMNS]),    
            rowMeans(Average_week1[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week1[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week1[INDEX_3, SOME_COLUMNS]),  
            rowMeans(Average_week2[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week2[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week2[INDEX_3, SOME_COLUMNS]),    
            rowMeans(Average_week4[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week4[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week4[INDEX_3, SOME_COLUMNS]),  
            rowMeans(Average_week6[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week6[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week6[INDEX_3, SOME_COLUMNS]),     
            rowMeans(Average_week8[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week8[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week8[INDEX_3, SOME_COLUMNS])  
  ),   
  sampleType2=c( rep("week0_1", length(INDEX_1)),    rep("week0_2", length(INDEX_2)),  rep("week0_3", length(INDEX_3)),          
                 rep("week1_1", length(INDEX_1)),    rep("week1_2", length(INDEX_2)),  rep("week1_3", length(INDEX_3)),        
                 rep("week2_1", length(INDEX_1)),    rep("week2_2", length(INDEX_2)),  rep("week2_3", length(INDEX_3)),       
                 rep("week4_1", length(INDEX_1)),    rep("week4_2", length(INDEX_2)),  rep("week4_3", length(INDEX_3)),     
                 rep("week6_1", length(INDEX_1)),    rep("week6_2", length(INDEX_2)),  rep("week6_3", length(INDEX_3)),      
                 rep("week8_1", length(INDEX_1)),    rep("week8_2", length(INDEX_2)),  rep("week8_3", length(INDEX_3)) 
  ), 
  sampleRank2=c( "week0_1",  "week1_1",   "week2_1",   "week4_1",  "week6_1", "week8_1",  
                 "week0_2",  "week1_2",   "week2_2",   "week4_2",  "week6_2", "week8_2",   
                 "week0_3",  "week1_3",   "week2_3",   "week4_3",  "week6_3", "week8_3"    
  ),     
  colours2=c( "week0_1"="red2",  "week1_1"="cyan2",   "week2_1"="green2",   "week4_1"="blue2",  "week6_1"="purple2", "week8_1"="orange2",  
              "week0_2"="red2",  "week1_2"="cyan2",   "week2_2"="green2",   "week4_2"="blue2",  "week6_2"="purple2", "week8_2"="orange2",   
              "week0_3"="red2",  "week1_3"="cyan2",   "week2_3"="green2",   "week4_3"="blue2",  "week6_3"="purple2", "week8_3"="orange2"    
  ),  
  path2=subdir_12_part3,   
  fileName2= paste("Part3-12-4B-WT-6Samples-BoxViolin",  "1kb",   mySample_g, sep = "_") ,  
  title2=myTitle_g,   xLab2="Samples",    yLab2="H2BGFP signal",   
  height2=4,   width2=10,   Ymin2=0, Ymax2=1.5)   ## width = 1 + 18*0.5, height=5cm 





MyBoxViolinPlot_1(vector2=c(rowMeans(Average_week0_EEDheto[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week0_EEDheto[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week0_EEDheto[INDEX_3, SOME_COLUMNS]),      
                            rowMeans(Average_week0_EEDko[INDEX_1, SOME_COLUMNS]),    rowMeans(Average_week0_EEDko[INDEX_2, SOME_COLUMNS]),    rowMeans(Average_week0_EEDko[INDEX_3, SOME_COLUMNS]),       
                            rowMeans(Average_week4_EEDheto[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week4_EEDheto[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week4_EEDheto[INDEX_3, SOME_COLUMNS]),      
                            rowMeans(Average_week4_EEDko[INDEX_1, SOME_COLUMNS]),    rowMeans(Average_week4_EEDko[INDEX_2, SOME_COLUMNS]),    rowMeans(Average_week4_EEDko[INDEX_3, SOME_COLUMNS])      
),   
sampleType2=c( rep("week0_EEDheto_1", length(INDEX_1)),    rep("week0_EEDheto_2", length(INDEX_2)),  rep("week0_EEDheto_3", length(INDEX_3)),        
               rep("week0_EEDko_1", length(INDEX_1)),      rep("week0_EEDko_2", length(INDEX_2)),    rep("week0_EEDko_3", length(INDEX_3)),           
               rep("week4_EEDheto_1", length(INDEX_1)),    rep("week4_EEDheto_2", length(INDEX_2)),  rep("week4_EEDheto_3", length(INDEX_3)),      
               rep("week4_EEDko_1", length(INDEX_1)),      rep("week4_EEDko_2", length(INDEX_2)),    rep("week4_EEDko_3", length(INDEX_3)) 
), 
sampleRank2=c( "week0_EEDheto_1",  "week0_EEDheto_2",   "week0_EEDheto_3",  
               "week0_EEDko_1",    "week0_EEDko_2",     "week0_EEDko_3",     
               "week4_EEDheto_1",  "week4_EEDheto_2",   "week4_EEDheto_3",  
               "week4_EEDko_1",    "week4_EEDko_2",     "week4_EEDko_3" 
),  
colours2=c(  "week0_EEDheto_1"="red",     "week0_EEDheto_2"="red",        "week0_EEDheto_3"="red",         
             "week0_EEDko_1"="red4",       "week0_EEDko_2"="red4",        "week0_EEDko_3"="red4",           
             "week4_EEDheto_1"="skyblue",  "week4_EEDheto_2"="skyblue",   "week4_EEDheto_3"="skyblue",   
             "week4_EEDko_1"="blue",       "week4_EEDko_2"="blue",        "week4_EEDko_3"="blue" 
),  
path2=subdir_12_part3,  fileName2= paste("Part3-12-4C-EEDko-4Samples-BoxViolin", "1kb",  mySample_g, sep = "_") , 
title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
height2=4.4,   width2=7, Ymin2=0, Ymax2=1.0 )   ## width = 1 + 12*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_week0_EEDheto[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week0_EEDheto[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week0_EEDheto[INDEX_3, SOME_COLUMNS]),     
                            rowMeans(Average_week0_EEDko[INDEX_1, SOME_COLUMNS]),    rowMeans(Average_week0_EEDko[INDEX_2, SOME_COLUMNS]),    rowMeans(Average_week0_EEDko[INDEX_3, SOME_COLUMNS]),         
                            rowMeans(Average_week4_EEDheto[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_week4_EEDheto[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_week4_EEDheto[INDEX_3, SOME_COLUMNS]),       
                            rowMeans(Average_week4_EEDko[INDEX_1, SOME_COLUMNS]),    rowMeans(Average_week4_EEDko[INDEX_2, SOME_COLUMNS]),    rowMeans(Average_week4_EEDko[INDEX_3, SOME_COLUMNS])     
),   
sampleType2=c( rep("week0_EEDheto_1", length(INDEX_1)),    rep("week0_EEDheto_2", length(INDEX_2)),  rep("week0_EEDheto_3", length(INDEX_3)),       
               rep("week0_EEDko_1", length(INDEX_1)),      rep("week0_EEDko_2", length(INDEX_2)),    rep("week0_EEDko_3", length(INDEX_3)),            
               rep("week4_EEDheto_1", length(INDEX_1)),    rep("week4_EEDheto_2", length(INDEX_2)),  rep("week4_EEDheto_3", length(INDEX_3)),        
               rep("week4_EEDko_1", length(INDEX_1)),      rep("week4_EEDko_2", length(INDEX_2)),    rep("week4_EEDko_3", length(INDEX_3)) 
), 
sampleRank2=c( "week0_EEDheto_1",  "week0_EEDko_1",    "week4_EEDheto_1",  "week4_EEDko_1", 
               "week0_EEDheto_2",  "week0_EEDko_2",    "week4_EEDheto_2",  "week4_EEDko_2", 
               "week0_EEDheto_3",  "week0_EEDko_3",    "week4_EEDheto_3",  "week4_EEDko_3" 
),  
colours2=c(  "week0_EEDheto_1"="red",  "week0_EEDko_1"="red4",    "week4_EEDheto_1"="skyblue",  "week4_EEDko_1"="blue", 
             "week0_EEDheto_2"="red",  "week0_EEDko_2"="red4",    "week4_EEDheto_2"="skyblue",  "week4_EEDko_2"="blue", 
             "week0_EEDheto_3"="red",  "week0_EEDko_3"="red4",    "week4_EEDheto_3"="skyblue",  "week4_EEDko_3"="blue" 
),  
path2=subdir_12_part3,  fileName2= paste("Part3-12-4D-EEDko-4Samples-BoxViolin", "1kb",  mySample_g, sep = "_") , 
title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
height2=4.4,   width2=7, Ymin2=0, Ymax2=1.0 )   ## width = 1 + 20*0.5, height=5cm 




MyBoxViolinPlot_1(vector2=c(rowMeans(Average_banding[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_banding[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_banding[INDEX_3, SOME_COLUMNS]),  
                            rowMeans(Average_sham[INDEX_1, SOME_COLUMNS]),     rowMeans(Average_sham[INDEX_2, SOME_COLUMNS]),     rowMeans(Average_sham[INDEX_3, SOME_COLUMNS])  
),   
sampleType2=c( rep("banding_1", length(INDEX_1)),    rep("banding_2", length(INDEX_2)),  rep("banding_3", length(INDEX_3)),   
               rep("sham_1", length(INDEX_1)),       rep("sham_2", length(INDEX_2)),     rep("sham_3", length(INDEX_3)) 
), 
sampleRank2=c( "banding_1",  "banding_2",   "banding_3",      
               "sham_1",     "sham_2",      "sham_3"  
),     
colours2=c( "banding_1"="blue",    "banding_2"="blue",     "banding_3"="blue",   
            "sham_1"="green4",     "sham_2"="green4",      "sham_3"="green4" 
), 
path2=subdir_12_part3,  fileName2= paste("Part3-12-4E-TAC-2Samples-BoxViolin", "1kb",  mySample_g, sep = "_") , 
title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
height2=4.2,   width2=4, Ymin2=0, Ymax2=0.8 )   ## width = 1 + 6*0.5, height=5cm 





MyBoxViolinPlot_1(vector2=c(rowMeans(Average_banding[INDEX_1, SOME_COLUMNS]),  rowMeans(Average_banding[INDEX_2, SOME_COLUMNS]),  rowMeans(Average_banding[INDEX_3, SOME_COLUMNS]),    
                            rowMeans(Average_sham[INDEX_1, SOME_COLUMNS]),     rowMeans(Average_sham[INDEX_2, SOME_COLUMNS]),     rowMeans(Average_sham[INDEX_3, SOME_COLUMNS]) 
),   
sampleType2=c( rep("banding_1", length(INDEX_1)),    rep("banding_2", length(INDEX_2)),  rep("banding_3", length(INDEX_3)),      
               rep("sham_1", length(INDEX_1)),       rep("sham_2", length(INDEX_2)),     rep("sham_3", length(INDEX_3)) 
), 
sampleRank2=c( "banding_1",   "sham_1", 
               "banding_2",   "sham_2", 
               "banding_3",   "sham_3"
),     
colours2=c( "banding_1"="blue",   "sham_1"="green4", 
            "banding_2"="blue",   "sham_2"="green4", 
            "banding_3"="blue",   "sham_3"="green4"
), 
path2=subdir_12_part3,  fileName2= paste("Part3-12-4F-TAC-2Samples-BoxViolin", "1kb",  mySample_g, sep = "_") , 
title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
height2=4.2,   width2=4, Ymin2=0, Ymax2=0.8 )   ## width = 1 + 10*0.5, height=5cm 



























###############################################################################
subdir_13_part3 <- paste(Part3_g,  "/13-averageRows-1kb-Scatter", sep = "")
if( ! file.exists(subdir_13_part3) ) { dir.create(subdir_13_part3) }



MyScatterDiagram_1(vector2=reduceColumn3_Average_H3, 
                   path2=subdir_13_part3,   
                   fileName2="3-13-1A-WT-H3", 
                   xScale2=0.001,  xLab2="DNA regions (x1000)",   
                   yLab2="H2BGFP signal",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.0, alpha2=0.5)

MyScatterDiagram_1(vector2=reduceColumn3_Average_week0, 
                   path2=subdir_13_part3,   
                   fileName2="3-13-1B-WT-week0", 
                   xScale2=0.001,  xLab2="DNA regions (x1000)",   
                   yLab2="H2BGFP signal",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.0, alpha2=0.5)


MyScatterDiagram_1(vector2=reduceColumn3_Average_week1, 
                   path2=subdir_13_part3,   
                   fileName2="3-13-1C-WT-week1", 
                   xScale2=0.001,  xLab2="DNA regions (x1000)",   
                   yLab2="H2BGFP signal",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.0, alpha2=0.5)

MyScatterDiagram_1(vector2=reduceColumn3_Average_week2, 
                   path2=subdir_13_part3,   
                   fileName2="3-13-1D-WT-week2", 
                   xScale2=0.001,  xLab2="DNA regions (x1000)",   
                   yLab2="H2BGFP signal",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.0, alpha2=0.5)

MyScatterDiagram_1(vector2=reduceColumn3_Average_week4, 
                   path2=subdir_13_part3,   
                   fileName2="3-13-1E-WT-week4", 
                   xScale2=0.001,  xLab2="DNA regions (x1000)",   
                   yLab2="H2BGFP signal",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.0, alpha2=0.5)

MyScatterDiagram_1(vector2=reduceColumn3_Average_week6, 
                   path2=subdir_13_part3,   
                   fileName2="3-13-1F-WT-week6", 
                   xScale2=0.001,  xLab2="DNA regions (x1000)",   
                   yLab2="H2BGFP signal",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.0, alpha2=0.5)

MyScatterDiagram_1(vector2=reduceColumn3_Average_week8, 
                   path2=subdir_13_part3,   
                   fileName2="3-13-1G-WT-week8", 
                   xScale2=0.001,  xLab2="DNA regions (x1000)",   
                   yLab2="H2BGFP signal",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.0, alpha2=0.5)






MyScatterDiagram_1(vector2=reduceColumn3_Average_week0_EEDheto, 
                   path2=subdir_13_part3,   
                   fileName2="3-13-2A-CKO-week0-EEDheto-scatter", 
                   xScale2=0.001,  xLab2="DNA regions (x1000)",   
                   yLab2="H2BGFP signal",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.0, alpha2=0.5)

MyScatterDiagram_1(vector2=reduceColumn3_Average_week0_EEDko, 
                   path2=subdir_13_part3,   
                   fileName2="3-13-2B-CKO-week0-EEDko-scatter", 
                   xScale2=0.001,  xLab2="DNA regions (x1000)",   
                   yLab2="H2BGFP signal",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.0, alpha2=0.5)


MyScatterDiagram_1(vector2=reduceColumn3_Average_week4_EEDheto, 
                   path2=subdir_13_part3,   
                   fileName2="3-13-2C-CKO-week4-EEDheto-scatter", 
                   xScale2=0.001,  xLab2="DNA regions (x1000)",   
                   yLab2="H2BGFP signal",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.0, alpha2=0.5)

MyScatterDiagram_1(vector2=reduceColumn3_Average_week4_EEDko, 
                   path2=subdir_13_part3,   
                   fileName2="3-13-2D-CKO-week4-EEDko-scatter", 
                   xScale2=0.001,  xLab2="DNA regions (x1000)",   
                   yLab2="H2BGFP signal",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.0, alpha2=0.5)






MyScatterDiagram_1(vector2=reduceColumn3_Average_banding, 
                   path2=subdir_13_part3,   
                   fileName2="3-13-3A-TAC-banding-scatter", 
                   xScale2=0.001,  xLab2="DNA regions (x1000)",   
                   yLab2="H2BGFP signal",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.0, alpha2=0.5)

MyScatterDiagram_1(vector2=reduceColumn3_Average_sham, 
                   path2=subdir_13_part3,   
                   fileName2="3-13-3B-TAC-sham-scatter", 
                   xScale2=0.001,  xLab2="DNA regions (x1000)",   
                   yLab2="H2BGFP signal",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.0, alpha2=0.5)










###############################################################################
subdir_14_part3 <- paste(Part3_g,  "/14-averageRows-1kb-Scatter-Cmp", sep = "")
if( ! file.exists(subdir_14_part3) ) { dir.create(subdir_14_part3) }

MyScatterDiagram_2(vector2X=reduceColumn3_Average_H3, vector2Y=reduceColumn3_Average_week0, 
                   path2=subdir_14_part3,   fileName2="3-14-1A-WT_H3-vs-week0", 
                   xLab2="H3",   yLab2="week0",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.6,   xMin2=0, xMax2=1.6,  alpha2=0.5)

MyScatterDiagram_2(vector2X=reduceColumn3_Average_week0, vector2Y=reduceColumn3_Average_week1, 
                   path2=subdir_14_part3,   fileName2="3-14-1B-WT_week0-vs-week1", 
                   xLab2="week0",   yLab2="week1",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.6,   xMin2=0, xMax2=1.6,  alpha2=0.5)

MyScatterDiagram_2(vector2X=reduceColumn3_Average_week1, vector2Y=reduceColumn3_Average_week2, 
                   path2=subdir_14_part3,   fileName2="3-14-1C-WT_week1-vs-week2", 
                   xLab2="week1",   yLab2="week2",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.6,   xMin2=0, xMax2=1.6,  alpha2=0.5)

MyScatterDiagram_2(vector2X=reduceColumn3_Average_week2, vector2Y=reduceColumn3_Average_week4, 
                   path2=subdir_14_part3,   fileName2="3-14-1D-WT_week2-vs-week4", 
                   xLab2="week2",   yLab2="week4",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.6,   xMin2=0, xMax2=1.6,  alpha2=0.5)

MyScatterDiagram_2(vector2X=reduceColumn3_Average_week4, vector2Y=reduceColumn3_Average_week6, 
                   path2=subdir_14_part3,   fileName2="3-14-1E-WT_week4-vs-week6", 
                   xLab2="week4",   yLab2="week6",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.6,   xMin2=0, xMax2=1.6,  alpha2=0.5)

MyScatterDiagram_2(vector2X=reduceColumn3_Average_week6, vector2Y=reduceColumn3_Average_week8, 
                   path2=subdir_14_part3,   fileName2="3-14-1F-WT_week6-vs-week8", 
                   xLab2="week6",   yLab2="week8",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=1.6,   xMin2=0, xMax2=1.6,  alpha2=0.5)






MyScatterDiagram_2(vector2X=reduceColumn3_Average_week0_EEDheto, vector2Y=reduceColumn3_Average_week0_EEDko, 
                   path2=subdir_14_part3,   fileName2="3-14-2A-CKO_week0EEDheto-vs-week0EEDko", 
                   xLab2="week0_EEDheto",   yLab2="week0_EEDko",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=0.8,   xMin2=0, xMax2=0.8,  alpha2=0.5)

MyScatterDiagram_2(vector2X=reduceColumn3_Average_week0_EEDheto, vector2Y=reduceColumn3_Average_week4_EEDheto, 
                   path2=subdir_14_part3,   fileName2="3-14-2B-CKO_week0EEDheto-vs-week4EEDheto", 
                   xLab2="week0_EEDheto",   yLab2="week4_EEDheto",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=0.8,   xMin2=0, xMax2=0.8,  alpha2=0.5)

MyScatterDiagram_2(vector2X=reduceColumn3_Average_week0_EEDko, vector2Y=reduceColumn3_Average_week4_EEDko, 
                   path2=subdir_14_part3,   fileName2="3-14-2C-CKO_week0EEDko-vs-week4EEDko", 
                   xLab2="week0_EEDko",   yLab2="week4_EEDko",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=0.8,   xMin2=0, xMax2=0.8,  alpha2=0.5)

MyScatterDiagram_2(vector2X=reduceColumn3_Average_week4_EEDheto, vector2Y=reduceColumn3_Average_week4_EEDko, 
                   path2=subdir_14_part3,   fileName2="3-14-2D-CKO_week4EEDheto-vs-week4EEDko", 
                   xLab2="week4_EEDheto",   yLab2="week4_EEDko",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=0.8,   xMin2=0, xMax2=0.8,  alpha2=0.5)





MyScatterDiagram_2(vector2X=reduceColumn3_Average_banding, vector2Y=reduceColumn3_Average_sham, 
                   path2=subdir_14_part3,   fileName2="3-14-3A-TAC_banding-vs-sham", 
                   xLab2="banding",   yLab2="sham",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=0.7,   xMin2=0, xMax2=0.7,  alpha2=0.5)















###############################################################################
subdir_15_part3 <- paste(Part3_g,  "/15-averageRows-1kb-Scatter-Colour-Cmp", sep = "")
if( ! file.exists(subdir_15_part3) ) { dir.create(subdir_15_part3) }

MyScatterDiagram_3(vector2X=reduceColumn3_Average_H3, vector2Y=reduceColumn3_Average_week0, 
                   path2=subdir_15_part3,   fileName2="3-15-1A-WT_H3-vs-week0", 
                   xLab2="H3",   yLab2="week0",  title2=myTitle_g,  
                   height2=3.5,   width2=4.6, yMin2=0, yMax2=1.6,   xMin2=0, xMax2=1.6,  alpha2=0.4, 
                   ratioThres2=1.5, colours2=c("red", "blue", "purple"))

MyScatterDiagram_3(vector2X=reduceColumn3_Average_week0, vector2Y=reduceColumn3_Average_week1, 
                   path2=subdir_15_part3,   fileName2="3-15-1B-WT_week0-vs-week1", 
                   xLab2="week0",   yLab2="week1",  title2=myTitle_g,  
                   height2=3.5,   width2=4.6, yMin2=0, yMax2=1.6,   xMin2=0, xMax2=1.6,  alpha2=0.4,
                   ratioThres2=1.5, colours2=c("red", "blue", "purple"))

MyScatterDiagram_3(vector2X=reduceColumn3_Average_week1, vector2Y=reduceColumn3_Average_week2, 
                   path2=subdir_15_part3,   fileName2="3-15-1C-WT_week1-vs-week2", 
                   xLab2="week1",   yLab2="week2",  title2=myTitle_g,  
                   height2=3.5,   width2=4.6, yMin2=0, yMax2=1.6,   xMin2=0, xMax2=1.6,  alpha2=0.4,
                   ratioThres2=1.5, colours2=c("red", "blue", "purple"))

MyScatterDiagram_3(vector2X=reduceColumn3_Average_week2, vector2Y=reduceColumn3_Average_week4, 
                   path2=subdir_15_part3,   fileName2="3-15-1D-WT_week2-vs-week4", 
                   xLab2="week2",   yLab2="week4",  title2=myTitle_g,  
                   height2=3.5,   width2=4.6, yMin2=0, yMax2=1.6,   xMin2=0, xMax2=1.6,  alpha2=0.4,
                   ratioThres2=1.5, colours2=c("red", "blue", "purple"))

MyScatterDiagram_3(vector2X=reduceColumn3_Average_week4, vector2Y=reduceColumn3_Average_week6, 
                   path2=subdir_15_part3,   fileName2="3-15-1E-WT_week4-vs-week6", 
                   xLab2="week4",   yLab2="week6",  title2=myTitle_g,  
                   height2=3.5,   width2=4.6, yMin2=0, yMax2=1.6,   xMin2=0, xMax2=1.6,  alpha2=0.4,
                   ratioThres2=1.5, colours2=c("red", "blue", "purple"))

MyScatterDiagram_3(vector2X=reduceColumn3_Average_week6, vector2Y=reduceColumn3_Average_week8, 
                   path2=subdir_15_part3,   fileName2="3-15-1F-WT_week6-vs-week8", 
                   xLab2="week6",   yLab2="week8",  title2=myTitle_g,  
                   height2=3.5,   width2=4.6, yMin2=0, yMax2=1.6,   xMin2=0, xMax2=1.6,  alpha2=0.4,
                   ratioThres2=1.5, colours2=c("red", "blue", "purple"))






MyScatterDiagram_3(vector2X=reduceColumn3_Average_week0_EEDheto, vector2Y=reduceColumn3_Average_week0_EEDko, 
                   path2=subdir_15_part3,   fileName2="3-15-2A-CKO_week0EEDheto-vs-week0EEDko", 
                   xLab2="week0_EEDheto",   yLab2="week0_EEDko",  title2=myTitle_g,  
                   height2=3.5,   width2=4.6, yMin2=0, yMax2=0.8,   xMin2=0, xMax2=0.8,  alpha2=0.4,
                   ratioThres2=1.5, colours2=c("red", "blue", "purple"))

MyScatterDiagram_3(vector2X=reduceColumn3_Average_week0_EEDheto, vector2Y=reduceColumn3_Average_week4_EEDheto, 
                   path2=subdir_15_part3,   fileName2="3-15-2B-CKO_week0EEDheto-vs-week4EEDheto", 
                   xLab2="week0_EEDheto",   yLab2="week4_EEDheto",  title2=myTitle_g,  
                   height2=3.5,   width2=4.6, yMin2=0, yMax2=0.8,   xMin2=0, xMax2=0.8,  alpha2=0.4,
                   ratioThres2=1.5, colours2=c("red", "blue", "purple"))

MyScatterDiagram_3(vector2X=reduceColumn3_Average_week0_EEDko, vector2Y=reduceColumn3_Average_week4_EEDko, 
                   path2=subdir_15_part3,   fileName2="3-15-2C-CKO_week0EEDko-vs-week4EEDko", 
                   xLab2="week0_EEDko",   yLab2="week4_EEDko",  title2=myTitle_g,  
                   height2=3.5,   width2=4.6, yMin2=0, yMax2=0.8,   xMin2=0, xMax2=0.8,  alpha2=0.4,
                   ratioThres2=1.5, colours2=c("red", "blue", "purple"))

MyScatterDiagram_3(vector2X=reduceColumn3_Average_week4_EEDheto, vector2Y=reduceColumn3_Average_week4_EEDko, 
                   path2=subdir_15_part3,   fileName2="3-15-2D-CKO_week4EEDheto-vs-week4EEDko", 
                   xLab2="week4_EEDheto",   yLab2="week4_EEDko",  title2=myTitle_g,  
                   height2=3.5,   width2=4.6, yMin2=0, yMax2=0.8,   xMin2=0, xMax2=0.8,  alpha2=0.4,
                   ratioThres2=1.5, colours2=c("red", "blue", "purple"))





MyScatterDiagram_3(vector2X=reduceColumn3_Average_banding, vector2Y=reduceColumn3_Average_sham, 
                   path2=subdir_15_part3,   fileName2="3-15-3A-TAC_banding-vs-sham", 
                   xLab2="banding",   yLab2="sham",  title2=myTitle_g,  
                   height2=3.5,   width2=4.6, yMin2=0, yMax2=0.7,   xMin2=0, xMax2=0.7,  alpha2=0.4, 
                   ratioThres2=1.5, colours2=c("red", "blue", "purple"))


















####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################





















