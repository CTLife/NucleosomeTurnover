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


############################################################################
## Part  6:  Classify DNA regions based on NTR.
############################################################################





###########################################################################
subdir_1_part6 <- paste(Part6_g,  "/1-sort-NTR-WT", sep = "")
if( ! file.exists(subdir_1_part6) ) { dir.create(subdir_1_part6) }

sink( file=paste(subdir_1_part6,   "Part6-1-A.runLog",  sep = "/") )


howToSort_6 <- function(vector1) {
  valueToSort = as.numeric( vector1 ) 
  print(valueToSort)
  return(valueToSort)
}


length(myNTR_7_WT_1)
toSort_6A   <- howToSort_6(myNTR_7_WT_1)   
length(toSort_6A)
index_6A     <-   order(toSort_6A)   ##  rev( order(toSort1) ) 
length(index_6A)
toSort_6A[ index_6A[1:100] ]

length_eachClass_6A <- length(index_6A)/5
calss1 <- index_6A[(0*length_eachClass_6A+1):(1*length_eachClass_6A)]
calss2 <- index_6A[(1*length_eachClass_6A+1):(2*length_eachClass_6A)]
calss3 <- index_6A[(2*length_eachClass_6A+1):(3*length_eachClass_6A)]
calss4 <- index_6A[(3*length_eachClass_6A+1):(4*length_eachClass_6A)]
calss5 <- index_6A[(4*length_eachClass_6A+1):(5*length_eachClass_6A)]
length(index_6A) - length(calss1)  - length(calss2)  - length(calss3)  - length(calss4)  - length(calss5) 


## NOL
dim(Average_H3) 
dim(Average_week0) 
dim(Average_week1) 
dim(Average_week2) 
dim(Average_week4) 
dim(Average_week6) 
dim(Average_week8) 

MyAverageLines_1(vector2=c( colMeans(Average_H3[calss1, ]),  colMeans(Average_H3[calss2, ]), colMeans(Average_H3[calss3, ]),  colMeans(Average_H3[calss4, ])  ,  colMeans(Average_H3[calss5, ]) ),   
                 numSample2=5,   
                 sampleType2=c( rep("H3_1", 500 ), rep("H3_2",  500),    rep("H3_3",  500 ),  rep("H3_4",  500 )  ,  rep("H3_5",  500 )  ), 
                 sampleRank2=c( "H3_1", "H3_2",  "H3_3",   "H3_4"  ,   "H3_5"   ),     
                 colours2=c( "H3_1"="yellow2",  "H3_2"="red2",   "H3_3"="cyan2",    "H3_4"="green2"   ,    "H3_5"="blue2"  ), 
                 path2=subdir_1_part6,     fileName2="6-1-1A-WT-H3-averageCurve",  
                 title2="H3 Level",     xLab2="Relative distance (kb)",    yLab2="H3 signal",   
                 Ymin2=0,   Ymax2=2.0,    height2=3.3,   width2=5.0 , center2=myCenter_g )



MyBoxViolinPlot_1(vector2=c( rowMeans(Average_H3[calss1, c(200:300)]),  rowMeans(Average_H3[calss2,  c(200:300)]), rowMeans(Average_H3[calss3,  c(200:300)]),  rowMeans(Average_H3[calss4,  c(200:300)])  ,  rowMeans(Average_H3[calss5,  c(200:300)]) ),       
                  sampleType2=c(  rep("H3_1", length(calss1) ), rep("H3_2",  length(calss2)),    rep("H3_3",  length(calss3) ),  rep("H3_4",  length(calss4) )  ,  rep("H3_5",  length(calss5) ) ), 
                  sampleRank2=c( "H3_1", "H3_2",  "H3_3",   "H3_4"  ,   "H3_5"  ),     
                  colours2=c("H3_1"="yellow2",  "H3_2"="red2",   "H3_3"="cyan2",    "H3_4"="green2"   ,    "H3_5"="blue2"   ), 
                  path2=subdir_1_part6,     fileName2="6-1-1A-WT-H3-violinPlot",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="H3 signal",   
                  Ymin2=0,   Ymax2=2,    height2=5,   width2=3.5  )


MyAverageLines_1(vector2=c( colMeans(Average_week0[calss1, ]),  colMeans(Average_week0[calss2, ]), colMeans(Average_week0[calss3, ]),  colMeans(Average_week0[calss4, ])  ,  colMeans(Average_week0[calss5, ]) ),   
                 numSample2=5,   
                 sampleType2=c( rep("week0_1", 500 ), rep("week0_2",  500),    rep("week0_3",  500 ),  rep("week0_4",  500 )  ,  rep("week0_5",  500 )  ), 
                 sampleRank2=c( "week0_1", "week0_2",  "week0_3",   "week0_4"  ,   "week0_5"   ),     
                 colours2=c( "week0_1"="yellow2",  "week0_2"="red2",   "week0_3"="cyan2",    "week0_4"="green2"   ,    "week0_5"="blue2"  ), 
                 path2=subdir_1_part6,     fileName2="6-1-1A-WT-week0-averageCurve",  
                 title2="week0 Level",     xLab2="Relative distance (kb)",    yLab2="week0 signal",   
                 Ymin2=0,   Ymax2=2.0,    height2=3.3,   width2=5.0 , center2=myCenter_g )



MyBoxViolinPlot_1(vector2=c( rowMeans(Average_week0[calss1, c(200:300)]),  rowMeans(Average_week0[calss2,  c(200:300)]), rowMeans(Average_week0[calss3,  c(200:300)]),  rowMeans(Average_week0[calss4,  c(200:300)])  ,  rowMeans(Average_week0[calss5,  c(200:300)]) ),       
                  sampleType2=c(  rep("week0_1", length(calss1) ), rep("week0_2",  length(calss2)),    rep("week0_3",  length(calss3) ),  rep("week0_4",  length(calss4) )  ,  rep("week0_5",  length(calss5) ) ), 
                  sampleRank2=c( "week0_1", "week0_2",  "week0_3",   "week0_4"  ,   "week0_5"  ),     
                  colours2=c("week0_1"="yellow2",  "week0_2"="red2",   "week0_3"="cyan2",    "week0_4"="green2"   ,    "week0_5"="blue2"   ), 
                  path2=subdir_1_part6,     fileName2="6-1-1A-WT-week0-violinPlot",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="week0 signal",   
                  Ymin2=0,   Ymax2=2,    height2=5,   width2=3.5  )



MyAverageLines_1(vector2=c( colMeans(Average_week8[calss1, ]),  colMeans(Average_week8[calss2, ]), colMeans(Average_week8[calss3, ]),  colMeans(Average_week8[calss4, ])  ,  colMeans(Average_week8[calss5, ]) ),   
                 numSample2=5,   
                 sampleType2=c( rep("week8_1", 500 ), rep("week8_2",  500),    rep("week8_3",  500 ),  rep("week8_4",  500 )  ,  rep("week8_5",  500 )  ), 
                 sampleRank2=c( "week8_1", "week8_2",  "week8_3",   "week8_4"  ,   "week8_5"   ),     
                 colours2=c( "week8_1"="yellow2",  "week8_2"="red2",   "week8_3"="cyan2",    "week8_4"="green2"   ,    "week8_5"="blue2"  ), 
                 path2=subdir_1_part6,     fileName2="6-1-1A-WT-week8-averageCurve",  
                 title2="week8 Level",     xLab2="Relative distance (kb)",    yLab2="week8 signal",   
                 Ymin2=0,   Ymax2=2.0,    height2=3.3,   width2=5.0 , center2=myCenter_g )



MyBoxViolinPlot_1(vector2=c( rowMeans(Average_week8[calss1, c(200:300)]),  rowMeans(Average_week8[calss2,  c(200:300)]), rowMeans(Average_week8[calss3,  c(200:300)]),  rowMeans(Average_week8[calss4,  c(200:300)])  ,  rowMeans(Average_week8[calss5,  c(200:300)]) ),       
                  sampleType2=c(  rep("week8_1", length(calss1) ), rep("week8_2",  length(calss2)),    rep("week8_3",  length(calss3) ),  rep("week8_4",  length(calss4) )  ,  rep("week8_5",  length(calss5) ) ), 
                  sampleRank2=c( "week8_1", "week8_2",  "week8_3",   "week8_4"  ,   "week8_5"  ),     
                  colours2=c("week8_1"="yellow2",  "week8_2"="red2",   "week8_3"="cyan2",    "week8_4"="green2"   ,    "week8_5"="blue2"   ), 
                  path2=subdir_1_part6,     fileName2="6-1-1A-WT-week8-violinPlot",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="week8 signal",   
                  Ymin2=0,   Ymax2=1,    height2=5,   width2=3.5  )











MyAverageLines_1(vector2=c(colMeans(Average_week0[calss1,]),  colMeans(Average_week1[calss1,]),   colMeans( Average_week2[calss1,]),   
                           colMeans(Average_week4[calss1,]),  colMeans(Average_week6[calss1,]),   colMeans(Average_week8[calss1,])  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_1_part6,     fileName2="6-1-1B-WT-6Samples-averageCurve",  
                 title2="NOL (Lowest)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0[calss2,]),  colMeans(Average_week1[calss2,]),    colMeans(Average_week2[calss2,]),   
                           colMeans(Average_week4[calss2,]),  colMeans(Average_week6[calss2,]),    colMeans(Average_week8[calss2,])  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_1_part6,     fileName2="6-1-1C-WT-6Samples-averageCurve",  
                 title2="NOL (Low)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0[calss3,]),  colMeans(Average_week1[calss3,]),  colMeans(Average_week2[calss3,]),   
                           colMeans(Average_week4[calss3,]),  colMeans(Average_week6[calss3,]),  colMeans(Average_week8[calss3,] ) ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_1_part6,     fileName2="6-1-1D-WT-6Samples-averageCurve",  
                 title2="NOL (Medium)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0[calss4,]),  colMeans(Average_week1[calss4,]),  colMeans(Average_week2[calss4,]),   
                           colMeans(Average_week4[calss4,]),  colMeans(Average_week6[calss4,]),  colMeans(Average_week8[calss4,] )   ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_1_part6,     fileName2="6-1-1E-WT-6Samples-averageCurve",  
                 title2="NOL (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )


MyAverageLines_1(vector2=c(colMeans(Average_week0[calss5,]),  colMeans(Average_week1[calss5,]),  colMeans(Average_week2[calss5,]),   
                           colMeans(Average_week4[calss5,]),  colMeans(Average_week6[calss5,]),  colMeans(Average_week8[calss5,] )   ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_1_part6,     fileName2="6-1-1F-WT-6Samples-averageCurve",  
                 title2="NOL (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )




MyBoxViolinPlot_1(vector2=c( myNTR_7_WT_1[calss1], myNTR_7_WT_1[calss2],  myNTR_7_WT_1[calss3],  myNTR_7_WT_1[calss4],  myNTR_7_WT_1[calss5] ),   
                  sampleType2=c( rep("Lowest", length(calss1)),    rep("Low", length(calss2)),  rep("Medium", length(calss3)),  
                                 rep("High", length(calss4)),      rep("Highest", length(calss5))  ),    
                  sampleRank2=c( "Lowest",  "Low",   "Medium",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",    "Medium"="green2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_1_part6,   fileName2= paste("6-1-1G-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=3.5,   Ymin2=0, Ymax2=0.4)    ## width = 1 + 5*0.5, height=5cm 










dim(Average_week0_EEDheto)   
dim(Average_week0_EEDko)  
dim(Average_week4_EEDheto) 
dim(Average_week4_EEDko)

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss1,]), colMeans( Average_week0_EEDko[calss1,]),      
                           colMeans(Average_week4_EEDheto[calss1,]),  colMeans(Average_week4_EEDko[calss1,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_1_part6,     fileName2="6-1-2A-CKO-4Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss2,]),  colMeans(Average_week0_EEDko[calss2,]),      
                           colMeans(Average_week4_EEDheto[calss2,]),  colMeans(Average_week4_EEDko[calss2,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_1_part6,     fileName2="6-1-2B-CKO-4Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss3,]),  colMeans(Average_week0_EEDko[calss3,]),      
                           colMeans(Average_week4_EEDheto[calss3,]),  colMeans(Average_week4_EEDko[calss3,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_1_part6,     fileName2="6-1-2C-CKO-4Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss4,]),  colMeans(Average_week0_EEDko[calss4,]),      
                           colMeans(Average_week4_EEDheto[calss4,]),  colMeans(Average_week4_EEDko[calss4,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_1_part6,     fileName2="6-1-2D-CKO-4Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss5,]),  colMeans(Average_week0_EEDko[calss5,]),      
                           colMeans(Average_week4_EEDheto[calss5,]),  colMeans(Average_week4_EEDko[calss5,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_1_part6,     fileName2="6-1-2E-CKO-4Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )





MyBoxViolinPlot_1(vector2=c( myNTR_7_EEDheto_2A[calss1], myNTR_7_EEDheto_2A[calss2],  myNTR_7_EEDheto_2A[calss3],  myNTR_7_EEDheto_2A[calss4],  myNTR_7_EEDheto_2A[calss5],
                             myNTR_7_EEDko_2B[calss1],   myNTR_7_EEDko_2B[calss2],    myNTR_7_EEDko_2B[calss3],    myNTR_7_EEDko_2B[calss4],    myNTR_7_EEDko_2B[calss5] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(calss1)),  rep("EEDheto_Low", length(calss2)),  rep("EEDheto_Medium", length(calss3)),  rep("EEDheto_High", length(calss4)),  rep("EEDheto_Highest", length(calss5)), 
                                 rep("EEDko_Lowest", length(calss1)),    rep("EEDko_Low", length(calss2)),    rep("EEDko_Medium", length(calss3)),    rep("EEDko_High", length(calss4)),    rep("EEDko_Highest", length(calss5)) ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",     "EEDko_Medium",    "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2",    "EEDko_Highest"="purple2" ), 
                  path2=subdir_1_part6,   fileName2= paste("6-1-2F-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=0.3)    ## width = 1 + 10*0.5, height=5cm 








dim(Average_banding)  
dim(Average_sham)

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss1,]),  colMeans(Average_sham[calss1,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_1_part6,     fileName2="6-1-3A-TAC-2Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss2,]),  colMeans(Average_sham[calss2,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_1_part6,     fileName2="6-1-3B-TAC-2Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss3,]),  colMeans(Average_sham[calss3,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_1_part6,     fileName2="6-1-3C-TAC-2Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss4,]),  colMeans(Average_sham[calss4,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_1_part6,     fileName2="6-1-3D-TAC-2Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )


MyAverageLines_1(vector2=c(colMeans(Average_banding[calss5,]),  colMeans(Average_sham[calss5,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_1_part6,     fileName2="6-1-3E-TAC-2Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )






MyBoxViolinPlot_1(vector2=c( myNTR_7_banding_3A[calss1], myNTR_7_banding_3A[calss2],  myNTR_7_banding_3A[calss3],  myNTR_7_banding_3A[calss4],  myNTR_7_banding_3A[calss5],
                             myNTR_7_sham_3B[calss1],    myNTR_7_sham_3B[calss2],     myNTR_7_sham_3B[calss3],     myNTR_7_sham_3B[calss4],     myNTR_7_sham_3B[calss5]  ),   
                  sampleType2=c( rep("banding_Lowest", length(calss1)),  rep("banding_Low", length(calss2)),  rep("banding_Medium", length(calss3)),  rep("banding_High", length(calss4)),  rep("banding_Highest", length(calss5)), 
                                 rep("sham_Lowest", length(calss1)),     rep("sham_Low", length(calss2)),     rep("sham_Medium", length(calss3)),     rep("sham_High", length(calss4)),     rep("sham_Highest", length(calss5))  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",   "banding_Medium",  "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",     "sham_Low",      "sham_Medium",     "sham_High",     "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",      "sham_Low"="cyan2",       "sham_Medium"="green2",     "sham_High"="blue2",     "sham_Highest"="purple2"  ), 
                  path2=subdir_1_part6,   fileName2= paste("6-1-3F-TAC-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=1.5)    ## width = 1 + 10*0.5, height=5cm 


sink() 



















###########################################################################
subdir_2_part6 <- paste(Part6_g,  "/2-sort-NTR-EEDheto", sep = "")
if( ! file.exists(subdir_2_part6) ) { dir.create(subdir_2_part6) }

sink( file=paste(subdir_2_part6,   "Part6-2-A.runLog",  sep = "/") )


howToSort_6 <- function(vector1) {
  valueToSort = as.numeric( vector1 ) 
  print(valueToSort)
  return(valueToSort)
}


length(myNTR_7_EEDheto_2A)
toSort_6A   <- howToSort_6(myNTR_7_EEDheto_2A)   
length(toSort_6A)
index_6A     <-   order(toSort_6A)      ##  rev( order(toSort1) ) 
length(index_6A)
toSort_6A[ index_6A[1:100] ]

length_eachClass_6A <- length(index_6A)/5
calss1 <- index_6A[(0*length_eachClass_6A+1):(1*length_eachClass_6A)]
calss2 <- index_6A[(1*length_eachClass_6A+1):(2*length_eachClass_6A)]
calss3 <- index_6A[(2*length_eachClass_6A+1):(3*length_eachClass_6A)]
calss4 <- index_6A[(3*length_eachClass_6A+1):(4*length_eachClass_6A)]
calss5 <- index_6A[(4*length_eachClass_6A+1):(5*length_eachClass_6A)]
length(index_6A) - length(calss1)  - length(calss2)  - length(calss3)  - length(calss4)  - length(calss5) 


## NOL
dim(Average_H3) 
dim(Average_week0) 
dim(Average_week1) 
dim(Average_week2) 
dim(Average_week4) 
dim(Average_week6) 
dim(Average_week8) 

MyAverageLines_1(vector2=c( colMeans(Average_H3[calss1, ]),  colMeans(Average_H3[calss2, ]), colMeans(Average_H3[calss3, ]),  colMeans(Average_H3[calss4, ])  ,  colMeans(Average_H3[calss5, ]) ),   
                 numSample2=5,   
                 sampleType2=c( rep("H3_1", 500 ), rep("H3_2",  500),    rep("H3_3",  500 ),  rep("H3_4",  500 )  ,  rep("H3_5",  500 )  ), 
                 sampleRank2=c( "H3_1", "H3_2",  "H3_3",   "H3_4"  ,   "H3_5"   ),     
                 colours2=c( "H3_1"="yellow2",  "H3_2"="red2",   "H3_3"="cyan2",    "H3_4"="green2"   ,    "H3_5"="blue2"  ), 
                 path2=subdir_2_part6,     fileName2="6-2-1A-WT-H3-averageCurve",  
                 title2="H3 Level",     xLab2="Relative distance (kb)",    yLab2="H3 signal",   
                 Ymin2=0,   Ymax2=2.0,    height2=3.3,   width2=5.0 , center2=myCenter_g )



MyBoxViolinPlot_1(vector2=c( rowMeans(Average_H3[calss1, c(200:300)]),  rowMeans(Average_H3[calss2,  c(200:300)]), rowMeans(Average_H3[calss3,  c(200:300)]),  rowMeans(Average_H3[calss4,  c(200:300)])  ,  rowMeans(Average_H3[calss5,  c(200:300)]) ),       
                  sampleType2=c(  rep("H3_1", length(calss1) ), rep("H3_2",  length(calss2)),    rep("H3_3",  length(calss3) ),  rep("H3_4",  length(calss4) )  ,  rep("H3_5",  length(calss5) ) ), 
                  sampleRank2=c( "H3_1", "H3_2",  "H3_3",   "H3_4"  ,   "H3_5"  ),     
                  colours2=c("H3_1"="yellow2",  "H3_2"="red2",   "H3_3"="cyan2",    "H3_4"="green2"   ,    "H3_5"="blue2"   ), 
                  path2=subdir_2_part6,     fileName2="6-2-1A-WT-H3-violinPlot",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="H3 signal",   
                  Ymin2=0,   Ymax2=2,    height2=5,   width2=3.5  )


MyAverageLines_1(vector2=c( colMeans(Average_week0[calss1, ]),  colMeans(Average_week0[calss2, ]), colMeans(Average_week0[calss3, ]),  colMeans(Average_week0[calss4, ])  ,  colMeans(Average_week0[calss5, ]) ),   
                 numSample2=5,   
                 sampleType2=c( rep("week0_1", 500 ), rep("week0_2",  500),    rep("week0_3",  500 ),  rep("week0_4",  500 )  ,  rep("week0_5",  500 )  ), 
                 sampleRank2=c( "week0_1", "week0_2",  "week0_3",   "week0_4"  ,   "week0_5"   ),     
                 colours2=c( "week0_1"="yellow2",  "week0_2"="red2",   "week0_3"="cyan2",    "week0_4"="green2"   ,    "week0_5"="blue2"  ), 
                 path2=subdir_2_part6,     fileName2="6-2-1A-WT-week0-averageCurve",  
                 title2="week0 Level",     xLab2="Relative distance (kb)",    yLab2="week0 signal",   
                 Ymin2=0,   Ymax2=2.0,    height2=3.3,   width2=5.0 , center2=myCenter_g )



MyBoxViolinPlot_1(vector2=c( rowMeans(Average_week0[calss1, c(200:300)]),  rowMeans(Average_week0[calss2,  c(200:300)]), rowMeans(Average_week0[calss3,  c(200:300)]),  rowMeans(Average_week0[calss4,  c(200:300)])  ,  rowMeans(Average_week0[calss5,  c(200:300)]) ),       
                  sampleType2=c(  rep("week0_1", length(calss1) ), rep("week0_2",  length(calss2)),    rep("week0_3",  length(calss3) ),  rep("week0_4",  length(calss4) )  ,  rep("week0_5",  length(calss5) ) ), 
                  sampleRank2=c( "week0_1", "week0_2",  "week0_3",   "week0_4"  ,   "week0_5"  ),     
                  colours2=c("week0_1"="yellow2",  "week0_2"="red2",   "week0_3"="cyan2",    "week0_4"="green2"   ,    "week0_5"="blue2"   ), 
                  path2=subdir_2_part6,     fileName2="6-2-1A-WT-week0-violinPlot",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="week0 signal",   
                  Ymin2=0,   Ymax2=2,    height2=5,   width2=3.5  )



MyAverageLines_1(vector2=c( colMeans(Average_week8[calss1, ]),  colMeans(Average_week8[calss2, ]), colMeans(Average_week8[calss3, ]),  colMeans(Average_week8[calss4, ])  ,  colMeans(Average_week8[calss5, ]) ),   
                 numSample2=5,   
                 sampleType2=c( rep("week8_1", 500 ), rep("week8_2",  500),    rep("week8_3",  500 ),  rep("week8_4",  500 )  ,  rep("week8_5",  500 )  ), 
                 sampleRank2=c( "week8_1", "week8_2",  "week8_3",   "week8_4"  ,   "week8_5"   ),     
                 colours2=c( "week8_1"="yellow2",  "week8_2"="red2",   "week8_3"="cyan2",    "week8_4"="green2"   ,    "week8_5"="blue2"  ), 
                 path2=subdir_2_part6,     fileName2="6-2-1A-WT-week8-averageCurve",  
                 title2="week8 Level",     xLab2="Relative distance (kb)",    yLab2="week8 signal",   
                 Ymin2=0,   Ymax2=2.0,    height2=3.3,   width2=5.0 , center2=myCenter_g )



MyBoxViolinPlot_1(vector2=c( rowMeans(Average_week8[calss1, c(200:300)]),  rowMeans(Average_week8[calss2,  c(200:300)]), rowMeans(Average_week8[calss3,  c(200:300)]),  rowMeans(Average_week8[calss4,  c(200:300)])  ,  rowMeans(Average_week8[calss5,  c(200:300)]) ),       
                  sampleType2=c(  rep("week8_1", length(calss1) ), rep("week8_2",  length(calss2)),    rep("week8_3",  length(calss3) ),  rep("week8_4",  length(calss4) )  ,  rep("week8_5",  length(calss5) ) ), 
                  sampleRank2=c( "week8_1", "week8_2",  "week8_3",   "week8_4"  ,   "week8_5"  ),     
                  colours2=c("week8_1"="yellow2",  "week8_2"="red2",   "week8_3"="cyan2",    "week8_4"="green2"   ,    "week8_5"="blue2"   ), 
                  path2=subdir_2_part6,     fileName2="6-2-1A-WT-week8-violinPlot",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="week8 signal",   
                  Ymin2=0,   Ymax2=1,    height2=5,   width2=3.5  )











MyAverageLines_1(vector2=c(colMeans(Average_week0[calss1,]),  colMeans(Average_week1[calss1,]),   colMeans( Average_week2[calss1,]),   
                           colMeans(Average_week4[calss1,]),  colMeans(Average_week6[calss1,]),   colMeans(Average_week8[calss1,])  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_2_part6,     fileName2="6-2-1B-WT-6Samples-averageCurve",  
                 title2="NOL (Lowest)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0[calss2,]),  colMeans(Average_week1[calss2,]),    colMeans(Average_week2[calss2,]),   
                           colMeans(Average_week4[calss2,]),  colMeans(Average_week6[calss2,]),    colMeans(Average_week8[calss2,])  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_2_part6,     fileName2="6-2-1C-WT-6Samples-averageCurve",  
                 title2="NOL (Low)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0[calss3,]),  colMeans(Average_week1[calss3,]),  colMeans(Average_week2[calss3,]),   
                           colMeans(Average_week4[calss3,]),  colMeans(Average_week6[calss3,]),  colMeans(Average_week8[calss3,] ) ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_2_part6,     fileName2="6-2-1D-WT-6Samples-averageCurve",  
                 title2="NOL (Medium)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0[calss4,]),  colMeans(Average_week1[calss4,]),  colMeans(Average_week2[calss4,]),   
                           colMeans(Average_week4[calss4,]),  colMeans(Average_week6[calss4,]),  colMeans(Average_week8[calss4,] )   ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_2_part6,     fileName2="6-2-1E-WT-6Samples-averageCurve",  
                 title2="NOL (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )


MyAverageLines_1(vector2=c(colMeans(Average_week0[calss5,]),  colMeans(Average_week1[calss5,]),  colMeans(Average_week2[calss5,]),   
                           colMeans(Average_week4[calss5,]),  colMeans(Average_week6[calss5,]),  colMeans(Average_week8[calss5,] )   ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_2_part6,     fileName2="6-2-1F-WT-6Samples-averageCurve",  
                 title2="NOL (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )




MyBoxViolinPlot_1(vector2=c( myNTR_7_WT_1[calss1], myNTR_7_WT_1[calss2],  myNTR_7_WT_1[calss3],  myNTR_7_WT_1[calss4],  myNTR_7_WT_1[calss5] ),   
                  sampleType2=c( rep("Lowest", length(calss1)),    rep("Low", length(calss2)),  rep("Medium", length(calss3)),  
                                 rep("High", length(calss4)),      rep("Highest", length(calss5))  ),    
                  sampleRank2=c( "Lowest",  "Low",   "Medium",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",    "Medium"="green2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_2_part6,   fileName2= paste("6-2-1G-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=3.5,   Ymin2=0, Ymax2=0.4)    ## width = 1 + 5*0.5, height=5cm 










dim(Average_week0_EEDheto)   
dim(Average_week0_EEDko)  
dim(Average_week4_EEDheto) 
dim(Average_week4_EEDko)

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss1,]), colMeans( Average_week0_EEDko[calss1,]),      
                           colMeans(Average_week4_EEDheto[calss1,]),  colMeans(Average_week4_EEDko[calss1,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_2_part6,     fileName2="6-2-2A-CKO-4Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss2,]),  colMeans(Average_week0_EEDko[calss2,]),      
                           colMeans(Average_week4_EEDheto[calss2,]),  colMeans(Average_week4_EEDko[calss2,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_2_part6,     fileName2="6-2-2B-CKO-4Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss3,]),  colMeans(Average_week0_EEDko[calss3,]),      
                           colMeans(Average_week4_EEDheto[calss3,]),  colMeans(Average_week4_EEDko[calss3,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_2_part6,     fileName2="6-2-2C-CKO-4Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss4,]),  colMeans(Average_week0_EEDko[calss4,]),      
                           colMeans(Average_week4_EEDheto[calss4,]),  colMeans(Average_week4_EEDko[calss4,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_2_part6,     fileName2="6-2-2D-CKO-4Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss5,]),  colMeans(Average_week0_EEDko[calss5,]),      
                           colMeans(Average_week4_EEDheto[calss5,]),  colMeans(Average_week4_EEDko[calss5,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_2_part6,     fileName2="6-2-2E-CKO-4Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )





MyBoxViolinPlot_1(vector2=c( myNTR_7_EEDheto_2A[calss1], myNTR_7_EEDheto_2A[calss2],  myNTR_7_EEDheto_2A[calss3],  myNTR_7_EEDheto_2A[calss4],  myNTR_7_EEDheto_2A[calss5],
                             myNTR_7_EEDko_2B[calss1],   myNTR_7_EEDko_2B[calss2],    myNTR_7_EEDko_2B[calss3],    myNTR_7_EEDko_2B[calss4],    myNTR_7_EEDko_2B[calss5] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(calss1)),  rep("EEDheto_Low", length(calss2)),  rep("EEDheto_Medium", length(calss3)),  rep("EEDheto_High", length(calss4)),  rep("EEDheto_Highest", length(calss5)), 
                                 rep("EEDko_Lowest", length(calss1)),    rep("EEDko_Low", length(calss2)),    rep("EEDko_Medium", length(calss3)),    rep("EEDko_High", length(calss4)),    rep("EEDko_Highest", length(calss5)) ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",     "EEDko_Medium",    "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2",    "EEDko_Highest"="purple2" ), 
                  path2=subdir_2_part6,   fileName2= paste("6-2-2F-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=0.3)    ## width = 1 + 10*0.5, height=5cm 








dim(Average_banding)  
dim(Average_sham)

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss1,]),  colMeans(Average_sham[calss1,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_2_part6,     fileName2="6-2-3A-TAC-2Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss2,]),  colMeans(Average_sham[calss2,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_2_part6,     fileName2="6-2-3B-TAC-2Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss3,]),  colMeans(Average_sham[calss3,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_2_part6,     fileName2="6-2-3C-TAC-2Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss4,]),  colMeans(Average_sham[calss4,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_2_part6,     fileName2="6-2-3D-TAC-2Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )


MyAverageLines_1(vector2=c(colMeans(Average_banding[calss5,]),  colMeans(Average_sham[calss5,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_2_part6,     fileName2="6-2-3E-TAC-2Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )






MyBoxViolinPlot_1(vector2=c( myNTR_7_banding_3A[calss1], myNTR_7_banding_3A[calss2],  myNTR_7_banding_3A[calss3],  myNTR_7_banding_3A[calss4],  myNTR_7_banding_3A[calss5],
                             myNTR_7_sham_3B[calss1],    myNTR_7_sham_3B[calss2],     myNTR_7_sham_3B[calss3],     myNTR_7_sham_3B[calss4],     myNTR_7_sham_3B[calss5]  ),   
                  sampleType2=c( rep("banding_Lowest", length(calss1)),  rep("banding_Low", length(calss2)),  rep("banding_Medium", length(calss3)),  rep("banding_High", length(calss4)),  rep("banding_Highest", length(calss5)), 
                                 rep("sham_Lowest", length(calss1)),     rep("sham_Low", length(calss2)),     rep("sham_Medium", length(calss3)),     rep("sham_High", length(calss4)),     rep("sham_Highest", length(calss5))  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",   "banding_Medium",  "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",     "sham_Low",      "sham_Medium",     "sham_High",     "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",      "sham_Low"="cyan2",       "sham_Medium"="green2",     "sham_High"="blue2",     "sham_Highest"="purple2"  ), 
                  path2=subdir_2_part6,   fileName2= paste("6-2-3F-TAC-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=1.5)    ## width = 1 + 10*0.5, height=5cm 



sink() 















###########################################################################
subdir_3_part6 <- paste(Part6_g,  "/3-sort-NTR-EEDko", sep = "")
if( ! file.exists(subdir_3_part6) ) { dir.create(subdir_3_part6) }

sink( file=paste(subdir_3_part6,   "Part6-3-A.runLog",  sep = "/") )


howToSort_6 <- function(vector1) {
  valueToSort = as.numeric( vector1 ) 
  print(valueToSort)
  return(valueToSort)
}


length(myNTR_7_EEDko_2B)
toSort_6A   <- howToSort_6(myNTR_7_EEDko_2B)   
length(toSort_6A)
index_6A     <-   order(toSort_6A)   ##  rev( order(toSort1) ) 
length(index_6A)
toSort_6A[ index_6A[1:100] ]

length_eachClass_6A <- length(index_6A)/5
calss1 <- index_6A[(0*length_eachClass_6A+1):(1*length_eachClass_6A)]
calss2 <- index_6A[(1*length_eachClass_6A+1):(2*length_eachClass_6A)]
calss3 <- index_6A[(2*length_eachClass_6A+1):(3*length_eachClass_6A)]
calss4 <- index_6A[(3*length_eachClass_6A+1):(4*length_eachClass_6A)]
calss5 <- index_6A[(4*length_eachClass_6A+1):(5*length_eachClass_6A)]
length(index_6A) - length(calss1)  - length(calss2)  - length(calss3)  - length(calss4)  - length(calss5) 


## NOL
dim(Average_H3) 
dim(Average_week0) 
dim(Average_week1) 
dim(Average_week2) 
dim(Average_week4) 
dim(Average_week6) 
dim(Average_week8) 

MyAverageLines_1(vector2=c( colMeans(Average_H3[calss1, ]),  colMeans(Average_H3[calss2, ]), colMeans(Average_H3[calss3, ]),  colMeans(Average_H3[calss4, ])  ,  colMeans(Average_H3[calss5, ]) ),   
                 numSample2=5,   
                 sampleType2=c( rep("H3_1", 500 ), rep("H3_2",  500),    rep("H3_3",  500 ),  rep("H3_4",  500 )  ,  rep("H3_5",  500 )  ), 
                 sampleRank2=c( "H3_1", "H3_2",  "H3_3",   "H3_4"  ,   "H3_5"   ),     
                 colours2=c( "H3_1"="yellow2",  "H3_2"="red2",   "H3_3"="cyan2",    "H3_4"="green2"   ,    "H3_5"="blue2"  ), 
                 path2=subdir_3_part6,     fileName2="6-3-1A-WT-H3-averageCurve",  
                 title2="H3 Level",     xLab2="Relative distance (kb)",    yLab2="H3 signal",   
                 Ymin2=0,   Ymax2=2.0,    height2=3.3,   width2=5.0 , center2=myCenter_g )



MyBoxViolinPlot_1(vector2=c( rowMeans(Average_H3[calss1, c(200:300)]),  rowMeans(Average_H3[calss2,  c(200:300)]), rowMeans(Average_H3[calss3,  c(200:300)]),  rowMeans(Average_H3[calss4,  c(200:300)])  ,  rowMeans(Average_H3[calss5,  c(200:300)]) ),       
                  sampleType2=c(  rep("H3_1", length(calss1) ), rep("H3_2",  length(calss2)),    rep("H3_3",  length(calss3) ),  rep("H3_4",  length(calss4) )  ,  rep("H3_5",  length(calss5) ) ), 
                  sampleRank2=c( "H3_1", "H3_2",  "H3_3",   "H3_4"  ,   "H3_5"  ),     
                  colours2=c("H3_1"="yellow2",  "H3_2"="red2",   "H3_3"="cyan2",    "H3_4"="green2"   ,    "H3_5"="blue2"   ), 
                  path2=subdir_3_part6,     fileName2="6-3-1A-WT-H3-violinPlot",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="H3 signal",   
                  Ymin2=0,   Ymax2=2,    height2=5,   width2=3.5  )


MyAverageLines_1(vector2=c( colMeans(Average_week0[calss1, ]),  colMeans(Average_week0[calss2, ]), colMeans(Average_week0[calss3, ]),  colMeans(Average_week0[calss4, ])  ,  colMeans(Average_week0[calss5, ]) ),   
                 numSample2=5,   
                 sampleType2=c( rep("week0_1", 500 ), rep("week0_2",  500),    rep("week0_3",  500 ),  rep("week0_4",  500 )  ,  rep("week0_5",  500 )  ), 
                 sampleRank2=c( "week0_1", "week0_2",  "week0_3",   "week0_4"  ,   "week0_5"   ),     
                 colours2=c( "week0_1"="yellow2",  "week0_2"="red2",   "week0_3"="cyan2",    "week0_4"="green2"   ,    "week0_5"="blue2"  ), 
                 path2=subdir_3_part6,     fileName2="6-3-1A-WT-week0-averageCurve",  
                 title2="week0 Level",     xLab2="Relative distance (kb)",    yLab2="week0 signal",   
                 Ymin2=0,   Ymax2=2.0,    height2=3.3,   width2=5.0 , center2=myCenter_g )



MyBoxViolinPlot_1(vector2=c( rowMeans(Average_week0[calss1, c(200:300)]),  rowMeans(Average_week0[calss2,  c(200:300)]), rowMeans(Average_week0[calss3,  c(200:300)]),  rowMeans(Average_week0[calss4,  c(200:300)])  ,  rowMeans(Average_week0[calss5,  c(200:300)]) ),       
                  sampleType2=c(  rep("week0_1", length(calss1) ), rep("week0_2",  length(calss2)),    rep("week0_3",  length(calss3) ),  rep("week0_4",  length(calss4) )  ,  rep("week0_5",  length(calss5) ) ), 
                  sampleRank2=c( "week0_1", "week0_2",  "week0_3",   "week0_4"  ,   "week0_5"  ),     
                  colours2=c("week0_1"="yellow2",  "week0_2"="red2",   "week0_3"="cyan2",    "week0_4"="green2"   ,    "week0_5"="blue2"   ), 
                  path2=subdir_3_part6,     fileName2="6-3-1A-WT-week0-violinPlot",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="week0 signal",   
                  Ymin2=0,   Ymax2=2,    height2=5,   width2=3.5  )



MyAverageLines_1(vector2=c( colMeans(Average_week8[calss1, ]),  colMeans(Average_week8[calss2, ]), colMeans(Average_week8[calss3, ]),  colMeans(Average_week8[calss4, ])  ,  colMeans(Average_week8[calss5, ]) ),   
                 numSample2=5,   
                 sampleType2=c( rep("week8_1", 500 ), rep("week8_2",  500),    rep("week8_3",  500 ),  rep("week8_4",  500 )  ,  rep("week8_5",  500 )  ), 
                 sampleRank2=c( "week8_1", "week8_2",  "week8_3",   "week8_4"  ,   "week8_5"   ),     
                 colours2=c( "week8_1"="yellow2",  "week8_2"="red2",   "week8_3"="cyan2",    "week8_4"="green2"   ,    "week8_5"="blue2"  ), 
                 path2=subdir_3_part6,     fileName2="6-3-1A-WT-week8-averageCurve",  
                 title2="week8 Level",     xLab2="Relative distance (kb)",    yLab2="week8 signal",   
                 Ymin2=0,   Ymax2=2.0,    height2=3.3,   width2=5.0 , center2=myCenter_g )



MyBoxViolinPlot_1(vector2=c( rowMeans(Average_week8[calss1, c(200:300)]),  rowMeans(Average_week8[calss2,  c(200:300)]), rowMeans(Average_week8[calss3,  c(200:300)]),  rowMeans(Average_week8[calss4,  c(200:300)])  ,  rowMeans(Average_week8[calss5,  c(200:300)]) ),       
                  sampleType2=c(  rep("week8_1", length(calss1) ), rep("week8_2",  length(calss2)),    rep("week8_3",  length(calss3) ),  rep("week8_4",  length(calss4) )  ,  rep("week8_5",  length(calss5) ) ), 
                  sampleRank2=c( "week8_1", "week8_2",  "week8_3",   "week8_4"  ,   "week8_5"  ),     
                  colours2=c("week8_1"="yellow2",  "week8_2"="red2",   "week8_3"="cyan2",    "week8_4"="green2"   ,    "week8_5"="blue2"   ), 
                  path2=subdir_3_part6,     fileName2="6-3-1A-WT-week8-violinPlot",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="week8 signal",   
                  Ymin2=0,   Ymax2=1,    height2=5,   width2=3.5  )











MyAverageLines_1(vector2=c(colMeans(Average_week0[calss1,]),  colMeans(Average_week1[calss1,]),   colMeans( Average_week2[calss1,]),   
                           colMeans(Average_week4[calss1,]),  colMeans(Average_week6[calss1,]),   colMeans(Average_week8[calss1,])  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_3_part6,     fileName2="6-3-1B-WT-6Samples-averageCurve",  
                 title2="NOL (Lowest)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0[calss2,]),  colMeans(Average_week1[calss2,]),    colMeans(Average_week2[calss2,]),   
                           colMeans(Average_week4[calss2,]),  colMeans(Average_week6[calss2,]),    colMeans(Average_week8[calss2,])  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_3_part6,     fileName2="6-3-1C-WT-6Samples-averageCurve",  
                 title2="NOL (Low)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0[calss3,]),  colMeans(Average_week1[calss3,]),  colMeans(Average_week2[calss3,]),   
                           colMeans(Average_week4[calss3,]),  colMeans(Average_week6[calss3,]),  colMeans(Average_week8[calss3,] ) ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_3_part6,     fileName2="6-3-1D-WT-6Samples-averageCurve",  
                 title2="NOL (Medium)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0[calss4,]),  colMeans(Average_week1[calss4,]),  colMeans(Average_week2[calss4,]),   
                           colMeans(Average_week4[calss4,]),  colMeans(Average_week6[calss4,]),  colMeans(Average_week8[calss4,] )   ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_3_part6,     fileName2="6-3-1E-WT-6Samples-averageCurve",  
                 title2="NOL (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )


MyAverageLines_1(vector2=c(colMeans(Average_week0[calss5,]),  colMeans(Average_week1[calss5,]),  colMeans(Average_week2[calss5,]),   
                           colMeans(Average_week4[calss5,]),  colMeans(Average_week6[calss5,]),  colMeans(Average_week8[calss5,] )   ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_3_part6,     fileName2="6-3-1F-WT-6Samples-averageCurve",  
                 title2="NOL (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )




MyBoxViolinPlot_1(vector2=c( myNTR_7_WT_1[calss1], myNTR_7_WT_1[calss2],  myNTR_7_WT_1[calss3],  myNTR_7_WT_1[calss4],  myNTR_7_WT_1[calss5] ),   
                  sampleType2=c( rep("Lowest", length(calss1)),    rep("Low", length(calss2)),  rep("Medium", length(calss3)),  
                                 rep("High", length(calss4)),      rep("Highest", length(calss5))  ),    
                  sampleRank2=c( "Lowest",  "Low",   "Medium",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",    "Medium"="green2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_3_part6,   fileName2= paste("6-3-1G-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=3.5,   Ymin2=0, Ymax2=0.4)    ## width = 1 + 5*0.5, height=5cm 










dim(Average_week0_EEDheto)   
dim(Average_week0_EEDko)  
dim(Average_week4_EEDheto) 
dim(Average_week4_EEDko)

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss1,]), colMeans( Average_week0_EEDko[calss1,]),      
                           colMeans(Average_week4_EEDheto[calss1,]),  colMeans(Average_week4_EEDko[calss1,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_3_part6,     fileName2="6-3-2A-CKO-4Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss2,]),  colMeans(Average_week0_EEDko[calss2,]),      
                           colMeans(Average_week4_EEDheto[calss2,]),  colMeans(Average_week4_EEDko[calss2,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_3_part6,     fileName2="6-3-2B-CKO-4Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss3,]),  colMeans(Average_week0_EEDko[calss3,]),      
                           colMeans(Average_week4_EEDheto[calss3,]),  colMeans(Average_week4_EEDko[calss3,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_3_part6,     fileName2="6-3-2C-CKO-4Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss4,]),  colMeans(Average_week0_EEDko[calss4,]),      
                           colMeans(Average_week4_EEDheto[calss4,]),  colMeans(Average_week4_EEDko[calss4,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_3_part6,     fileName2="6-3-2D-CKO-4Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss5,]),  colMeans(Average_week0_EEDko[calss5,]),      
                           colMeans(Average_week4_EEDheto[calss5,]),  colMeans(Average_week4_EEDko[calss5,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_3_part6,     fileName2="6-3-2E-CKO-4Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )





MyBoxViolinPlot_1(vector2=c( myNTR_7_EEDheto_2A[calss1], myNTR_7_EEDheto_2A[calss2],  myNTR_7_EEDheto_2A[calss3],  myNTR_7_EEDheto_2A[calss4],  myNTR_7_EEDheto_2A[calss5],
                             myNTR_7_EEDko_2B[calss1],   myNTR_7_EEDko_2B[calss2],    myNTR_7_EEDko_2B[calss3],    myNTR_7_EEDko_2B[calss4],    myNTR_7_EEDko_2B[calss5] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(calss1)),  rep("EEDheto_Low", length(calss2)),  rep("EEDheto_Medium", length(calss3)),  rep("EEDheto_High", length(calss4)),  rep("EEDheto_Highest", length(calss5)), 
                                 rep("EEDko_Lowest", length(calss1)),    rep("EEDko_Low", length(calss2)),    rep("EEDko_Medium", length(calss3)),    rep("EEDko_High", length(calss4)),    rep("EEDko_Highest", length(calss5)) ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",     "EEDko_Medium",    "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2",    "EEDko_Highest"="purple2" ), 
                  path2=subdir_3_part6,   fileName2= paste("6-3-2F-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=0.3)    ## width = 1 + 10*0.5, height=5cm 








dim(Average_banding)  
dim(Average_sham)

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss1,]),  colMeans(Average_sham[calss1,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_3_part6,     fileName2="6-3-3A-TAC-2Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss2,]),  colMeans(Average_sham[calss2,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_3_part6,     fileName2="6-3-3B-TAC-2Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss3,]),  colMeans(Average_sham[calss3,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_3_part6,     fileName2="6-3-3C-TAC-2Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss4,]),  colMeans(Average_sham[calss4,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_3_part6,     fileName2="6-3-3D-TAC-2Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )


MyAverageLines_1(vector2=c(colMeans(Average_banding[calss5,]),  colMeans(Average_sham[calss5,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_3_part6,     fileName2="6-3-3E-TAC-2Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )






MyBoxViolinPlot_1(vector2=c( myNTR_7_banding_3A[calss1], myNTR_7_banding_3A[calss2],  myNTR_7_banding_3A[calss3],  myNTR_7_banding_3A[calss4],  myNTR_7_banding_3A[calss5],
                             myNTR_7_sham_3B[calss1],    myNTR_7_sham_3B[calss2],     myNTR_7_sham_3B[calss3],     myNTR_7_sham_3B[calss4],     myNTR_7_sham_3B[calss5]  ),   
                  sampleType2=c( rep("banding_Lowest", length(calss1)),  rep("banding_Low", length(calss2)),  rep("banding_Medium", length(calss3)),  rep("banding_High", length(calss4)),  rep("banding_Highest", length(calss5)), 
                                 rep("sham_Lowest", length(calss1)),     rep("sham_Low", length(calss2)),     rep("sham_Medium", length(calss3)),     rep("sham_High", length(calss4)),     rep("sham_Highest", length(calss5))  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",   "banding_Medium",  "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",     "sham_Low",      "sham_Medium",     "sham_High",     "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",      "sham_Low"="cyan2",       "sham_Medium"="green2",     "sham_High"="blue2",     "sham_Highest"="purple2"  ), 
                  path2=subdir_3_part6,   fileName2= paste("6-3-3F-TAC-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=1.5)    ## width = 1 + 10*0.5, height=5cm 


sink()



















###########################################################################
subdir_4_part6 <- paste(Part6_g,  "/4-sort-NTR-banding", sep = "")
if( ! file.exists(subdir_4_part6) ) { dir.create(subdir_4_part6) }

sink( file=paste(subdir_4_part6,   "Part6-4-A.runLog",  sep = "/") )


howToSort_6 <- function(vector1) {
  valueToSort = as.numeric( vector1 ) 
  print(valueToSort)
  return(valueToSort)
}


length(myNTR_7_banding_3A)
toSort_6A   <- howToSort_6(myNTR_7_banding_3A)   
length(toSort_6A)
index_6A     <-   order(toSort_6A)   ##  rev( order(toSort1) ) 
length(index_6A)
toSort_6A[ index_6A[1:100] ]

length_eachClass_6A <- length(index_6A)/5
calss1 <- index_6A[(0*length_eachClass_6A+1):(1*length_eachClass_6A)]
calss2 <- index_6A[(1*length_eachClass_6A+1):(2*length_eachClass_6A)]
calss3 <- index_6A[(2*length_eachClass_6A+1):(3*length_eachClass_6A)]
calss4 <- index_6A[(3*length_eachClass_6A+1):(4*length_eachClass_6A)]
calss5 <- index_6A[(4*length_eachClass_6A+1):(5*length_eachClass_6A)]
length(index_6A) - length(calss1)  - length(calss2)  - length(calss3)  - length(calss4)  - length(calss5) 


## NOL
dim(Average_H3) 
dim(Average_week0) 
dim(Average_week1) 
dim(Average_week2) 
dim(Average_week4) 
dim(Average_week6) 
dim(Average_week8) 

MyAverageLines_1(vector2=c( colMeans(Average_H3[calss1, ]),  colMeans(Average_H3[calss2, ]), colMeans(Average_H3[calss3, ]),  colMeans(Average_H3[calss4, ])  ,  colMeans(Average_H3[calss5, ]) ),   
                 numSample2=5,   
                 sampleType2=c( rep("H3_1", 500 ), rep("H3_2",  500),    rep("H3_3",  500 ),  rep("H3_4",  500 )  ,  rep("H3_5",  500 )  ), 
                 sampleRank2=c( "H3_1", "H3_2",  "H3_3",   "H3_4"  ,   "H3_5"   ),     
                 colours2=c( "H3_1"="yellow2",  "H3_2"="red2",   "H3_3"="cyan2",    "H3_4"="green2"   ,    "H3_5"="blue2"  ), 
                 path2=subdir_4_part6,     fileName2="6-4-1A-WT-H3-averageCurve",  
                 title2="H3 Level",     xLab2="Relative distance (kb)",    yLab2="H3 signal",   
                 Ymin2=0,   Ymax2=2.0,    height2=3.3,   width2=5.0 , center2=myCenter_g )



MyBoxViolinPlot_1(vector2=c( rowMeans(Average_H3[calss1, c(200:300)]),  rowMeans(Average_H3[calss2,  c(200:300)]), rowMeans(Average_H3[calss3,  c(200:300)]),  rowMeans(Average_H3[calss4,  c(200:300)])  ,  rowMeans(Average_H3[calss5,  c(200:300)]) ),       
                  sampleType2=c(  rep("H3_1", length(calss1) ), rep("H3_2",  length(calss2)),    rep("H3_3",  length(calss3) ),  rep("H3_4",  length(calss4) )  ,  rep("H3_5",  length(calss5) ) ), 
                  sampleRank2=c( "H3_1", "H3_2",  "H3_3",   "H3_4"  ,   "H3_5"  ),     
                  colours2=c("H3_1"="yellow2",  "H3_2"="red2",   "H3_3"="cyan2",    "H3_4"="green2"   ,    "H3_5"="blue2"   ), 
                  path2=subdir_4_part6,     fileName2="6-4-1A-WT-H3-violinPlot",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="H3 signal",   
                  Ymin2=0,   Ymax2=2,    height2=5,   width2=3.5  )


MyAverageLines_1(vector2=c( colMeans(Average_week0[calss1, ]),  colMeans(Average_week0[calss2, ]), colMeans(Average_week0[calss3, ]),  colMeans(Average_week0[calss4, ])  ,  colMeans(Average_week0[calss5, ]) ),   
                 numSample2=5,   
                 sampleType2=c( rep("week0_1", 500 ), rep("week0_2",  500),    rep("week0_3",  500 ),  rep("week0_4",  500 )  ,  rep("week0_5",  500 )  ), 
                 sampleRank2=c( "week0_1", "week0_2",  "week0_3",   "week0_4"  ,   "week0_5"   ),     
                 colours2=c( "week0_1"="yellow2",  "week0_2"="red2",   "week0_3"="cyan2",    "week0_4"="green2"   ,    "week0_5"="blue2"  ), 
                 path2=subdir_4_part6,     fileName2="6-4-1A-WT-week0-averageCurve",  
                 title2="week0 Level",     xLab2="Relative distance (kb)",    yLab2="week0 signal",   
                 Ymin2=0,   Ymax2=2.0,    height2=3.3,   width2=5.0 , center2=myCenter_g )



MyBoxViolinPlot_1(vector2=c( rowMeans(Average_week0[calss1, c(200:300)]),  rowMeans(Average_week0[calss2,  c(200:300)]), rowMeans(Average_week0[calss3,  c(200:300)]),  rowMeans(Average_week0[calss4,  c(200:300)])  ,  rowMeans(Average_week0[calss5,  c(200:300)]) ),       
                  sampleType2=c(  rep("week0_1", length(calss1) ), rep("week0_2",  length(calss2)),    rep("week0_3",  length(calss3) ),  rep("week0_4",  length(calss4) )  ,  rep("week0_5",  length(calss5) ) ), 
                  sampleRank2=c( "week0_1", "week0_2",  "week0_3",   "week0_4"  ,   "week0_5"  ),     
                  colours2=c("week0_1"="yellow2",  "week0_2"="red2",   "week0_3"="cyan2",    "week0_4"="green2"   ,    "week0_5"="blue2"   ), 
                  path2=subdir_4_part6,     fileName2="6-4-1A-WT-week0-violinPlot",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="week0 signal",   
                  Ymin2=0,   Ymax2=2,    height2=5,   width2=3.5  )



MyAverageLines_1(vector2=c( colMeans(Average_week8[calss1, ]),  colMeans(Average_week8[calss2, ]), colMeans(Average_week8[calss3, ]),  colMeans(Average_week8[calss4, ])  ,  colMeans(Average_week8[calss5, ]) ),   
                 numSample2=5,   
                 sampleType2=c( rep("week8_1", 500 ), rep("week8_2",  500),    rep("week8_3",  500 ),  rep("week8_4",  500 )  ,  rep("week8_5",  500 )  ), 
                 sampleRank2=c( "week8_1", "week8_2",  "week8_3",   "week8_4"  ,   "week8_5"   ),     
                 colours2=c( "week8_1"="yellow2",  "week8_2"="red2",   "week8_3"="cyan2",    "week8_4"="green2"   ,    "week8_5"="blue2"  ), 
                 path2=subdir_4_part6,     fileName2="6-4-1A-WT-week8-averageCurve",  
                 title2="week8 Level",     xLab2="Relative distance (kb)",    yLab2="week8 signal",   
                 Ymin2=0,   Ymax2=2.0,    height2=3.3,   width2=5.0 , center2=myCenter_g )



MyBoxViolinPlot_1(vector2=c( rowMeans(Average_week8[calss1, c(200:300)]),  rowMeans(Average_week8[calss2,  c(200:300)]), rowMeans(Average_week8[calss3,  c(200:300)]),  rowMeans(Average_week8[calss4,  c(200:300)])  ,  rowMeans(Average_week8[calss5,  c(200:300)]) ),       
                  sampleType2=c(  rep("week8_1", length(calss1) ), rep("week8_2",  length(calss2)),    rep("week8_3",  length(calss3) ),  rep("week8_4",  length(calss4) )  ,  rep("week8_5",  length(calss5) ) ), 
                  sampleRank2=c( "week8_1", "week8_2",  "week8_3",   "week8_4"  ,   "week8_5"  ),     
                  colours2=c("week8_1"="yellow2",  "week8_2"="red2",   "week8_3"="cyan2",    "week8_4"="green2"   ,    "week8_5"="blue2"   ), 
                  path2=subdir_4_part6,     fileName2="6-4-1A-WT-week8-violinPlot",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="week8 signal",   
                  Ymin2=0,   Ymax2=1,    height2=5,   width2=3.5  )











MyAverageLines_1(vector2=c(colMeans(Average_week0[calss1,]),  colMeans(Average_week1[calss1,]),   colMeans( Average_week2[calss1,]),   
                           colMeans(Average_week4[calss1,]),  colMeans(Average_week6[calss1,]),   colMeans(Average_week8[calss1,])  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_4_part6,     fileName2="6-4-1B-WT-6Samples-averageCurve",  
                 title2="NOL (Lowest)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0[calss2,]),  colMeans(Average_week1[calss2,]),    colMeans(Average_week2[calss2,]),   
                           colMeans(Average_week4[calss2,]),  colMeans(Average_week6[calss2,]),    colMeans(Average_week8[calss2,])  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_4_part6,     fileName2="6-4-1C-WT-6Samples-averageCurve",  
                 title2="NOL (Low)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0[calss3,]),  colMeans(Average_week1[calss3,]),  colMeans(Average_week2[calss3,]),   
                           colMeans(Average_week4[calss3,]),  colMeans(Average_week6[calss3,]),  colMeans(Average_week8[calss3,] ) ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_4_part6,     fileName2="6-4-1D-WT-6Samples-averageCurve",  
                 title2="NOL (Medium)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0[calss4,]),  colMeans(Average_week1[calss4,]),  colMeans(Average_week2[calss4,]),   
                           colMeans(Average_week4[calss4,]),  colMeans(Average_week6[calss4,]),  colMeans(Average_week8[calss4,] )   ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_4_part6,     fileName2="6-4-1E-WT-6Samples-averageCurve",  
                 title2="NOL (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )


MyAverageLines_1(vector2=c(colMeans(Average_week0[calss5,]),  colMeans(Average_week1[calss5,]),  colMeans(Average_week2[calss5,]),   
                           colMeans(Average_week4[calss5,]),  colMeans(Average_week6[calss5,]),  colMeans(Average_week8[calss5,] )   ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_4_part6,     fileName2="6-4-1F-WT-6Samples-averageCurve",  
                 title2="NOL (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )




MyBoxViolinPlot_1(vector2=c( myNTR_7_WT_1[calss1], myNTR_7_WT_1[calss2],  myNTR_7_WT_1[calss3],  myNTR_7_WT_1[calss4],  myNTR_7_WT_1[calss5] ),   
                  sampleType2=c( rep("Lowest", length(calss1)),    rep("Low", length(calss2)),  rep("Medium", length(calss3)),  
                                 rep("High", length(calss4)),      rep("Highest", length(calss5))  ),    
                  sampleRank2=c( "Lowest",  "Low",   "Medium",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",    "Medium"="green2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_4_part6,   fileName2= paste("6-4-1G-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=3.5,   Ymin2=0, Ymax2=0.4)    ## width = 1 + 5*0.5, height=5cm 










dim(Average_week0_EEDheto)   
dim(Average_week0_EEDko)  
dim(Average_week4_EEDheto) 
dim(Average_week4_EEDko)

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss1,]), colMeans( Average_week0_EEDko[calss1,]),      
                           colMeans(Average_week4_EEDheto[calss1,]),  colMeans(Average_week4_EEDko[calss1,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_4_part6,     fileName2="6-4-2A-CKO-4Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss2,]),  colMeans(Average_week0_EEDko[calss2,]),      
                           colMeans(Average_week4_EEDheto[calss2,]),  colMeans(Average_week4_EEDko[calss2,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_4_part6,     fileName2="6-4-2B-CKO-4Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss3,]),  colMeans(Average_week0_EEDko[calss3,]),      
                           colMeans(Average_week4_EEDheto[calss3,]),  colMeans(Average_week4_EEDko[calss3,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_4_part6,     fileName2="6-4-2C-CKO-4Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss4,]),  colMeans(Average_week0_EEDko[calss4,]),      
                           colMeans(Average_week4_EEDheto[calss4,]),  colMeans(Average_week4_EEDko[calss4,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_4_part6,     fileName2="6-4-2D-CKO-4Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss5,]),  colMeans(Average_week0_EEDko[calss5,]),      
                           colMeans(Average_week4_EEDheto[calss5,]),  colMeans(Average_week4_EEDko[calss5,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_4_part6,     fileName2="6-4-2E-CKO-4Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )





MyBoxViolinPlot_1(vector2=c( myNTR_7_EEDheto_2A[calss1], myNTR_7_EEDheto_2A[calss2],  myNTR_7_EEDheto_2A[calss3],  myNTR_7_EEDheto_2A[calss4],  myNTR_7_EEDheto_2A[calss5],
                             myNTR_7_EEDko_2B[calss1],   myNTR_7_EEDko_2B[calss2],    myNTR_7_EEDko_2B[calss3],    myNTR_7_EEDko_2B[calss4],    myNTR_7_EEDko_2B[calss5] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(calss1)),  rep("EEDheto_Low", length(calss2)),  rep("EEDheto_Medium", length(calss3)),  rep("EEDheto_High", length(calss4)),  rep("EEDheto_Highest", length(calss5)), 
                                 rep("EEDko_Lowest", length(calss1)),    rep("EEDko_Low", length(calss2)),    rep("EEDko_Medium", length(calss3)),    rep("EEDko_High", length(calss4)),    rep("EEDko_Highest", length(calss5)) ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",     "EEDko_Medium",    "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2",    "EEDko_Highest"="purple2" ), 
                  path2=subdir_4_part6,   fileName2= paste("6-4-2F-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=0.3)    ## width = 1 + 10*0.5, height=5cm 








dim(Average_banding)  
dim(Average_sham)

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss1,]),  colMeans(Average_sham[calss1,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_4_part6,     fileName2="6-4-3A-TAC-2Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss2,]),  colMeans(Average_sham[calss2,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_4_part6,     fileName2="6-4-3B-TAC-2Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss3,]),  colMeans(Average_sham[calss3,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_4_part6,     fileName2="6-4-3C-TAC-2Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss4,]),  colMeans(Average_sham[calss4,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_4_part6,     fileName2="6-4-3D-TAC-2Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )


MyAverageLines_1(vector2=c(colMeans(Average_banding[calss5,]),  colMeans(Average_sham[calss5,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_4_part6,     fileName2="6-4-3E-TAC-2Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )






MyBoxViolinPlot_1(vector2=c( myNTR_7_banding_3A[calss1], myNTR_7_banding_3A[calss2],  myNTR_7_banding_3A[calss3],  myNTR_7_banding_3A[calss4],  myNTR_7_banding_3A[calss5],
                             myNTR_7_sham_3B[calss1],    myNTR_7_sham_3B[calss2],     myNTR_7_sham_3B[calss3],     myNTR_7_sham_3B[calss4],     myNTR_7_sham_3B[calss5]  ),   
                  sampleType2=c( rep("banding_Lowest", length(calss1)),  rep("banding_Low", length(calss2)),  rep("banding_Medium", length(calss3)),  rep("banding_High", length(calss4)),  rep("banding_Highest", length(calss5)), 
                                 rep("sham_Lowest", length(calss1)),     rep("sham_Low", length(calss2)),     rep("sham_Medium", length(calss3)),     rep("sham_High", length(calss4)),     rep("sham_Highest", length(calss5))  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",   "banding_Medium",  "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",     "sham_Low",      "sham_Medium",     "sham_High",     "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",      "sham_Low"="cyan2",       "sham_Medium"="green2",     "sham_High"="blue2",     "sham_Highest"="purple2"  ), 
                  path2=subdir_4_part6,   fileName2= paste("6-4-3F-TAC-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=1.5)    ## width = 1 + 10*0.5, height=5cm 


sink() 


















###########################################################################
subdir_5_part6 <- paste(Part6_g,  "/5-sort-NTR-sham", sep = "")
if( ! file.exists(subdir_5_part6) ) { dir.create(subdir_5_part6) }

sink( file=paste(subdir_5_part6,   "Part6-5-A.runLog",  sep = "/") )


howToSort_6 <- function(vector1) {
  valueToSort = as.numeric( vector1 ) 
  print(valueToSort)
  return(valueToSort)
}


length( myNTR_7_sham_3B )
toSort_6A   <- howToSort_6( myNTR_7_sham_3B )   
length(toSort_6A)
index_6A     <-   order(toSort_6A)   ##  rev( order(toSort1) ) 
length(index_6A)
toSort_6A[ index_6A[1:100] ]

length_eachClass_6A <- length(index_6A)/5
calss1 <- index_6A[(0*length_eachClass_6A+1):(1*length_eachClass_6A)]
calss2 <- index_6A[(1*length_eachClass_6A+1):(2*length_eachClass_6A)]
calss3 <- index_6A[(2*length_eachClass_6A+1):(3*length_eachClass_6A)]
calss4 <- index_6A[(3*length_eachClass_6A+1):(4*length_eachClass_6A)]
calss5 <- index_6A[(4*length_eachClass_6A+1):(5*length_eachClass_6A)]
length(index_6A) - length(calss1)  - length(calss2)  - length(calss3)  - length(calss4)  - length(calss5) 


## NOL
dim(Average_H3) 
dim(Average_week0) 
dim(Average_week1) 
dim(Average_week2) 
dim(Average_week4) 
dim(Average_week6) 
dim(Average_week8) 

MyAverageLines_1(vector2=c( colMeans(Average_H3[calss1, ]),  colMeans(Average_H3[calss2, ]), colMeans(Average_H3[calss3, ]),  colMeans(Average_H3[calss4, ])  ,  colMeans(Average_H3[calss5, ]) ),   
                 numSample2=5,   
                 sampleType2=c( rep("H3_1", 500 ), rep("H3_2",  500),    rep("H3_3",  500 ),  rep("H3_4",  500 )  ,  rep("H3_5",  500 )  ), 
                 sampleRank2=c( "H3_1", "H3_2",  "H3_3",   "H3_4"  ,   "H3_5"   ),     
                 colours2=c( "H3_1"="yellow2",  "H3_2"="red2",   "H3_3"="cyan2",    "H3_4"="green2"   ,    "H3_5"="blue2"  ), 
                 path2=subdir_5_part6,     fileName2="6-5-1A-WT-H3-averageCurve",  
                 title2="H3 Level",     xLab2="Relative distance (kb)",    yLab2="H3 signal",   
                 Ymin2=0,   Ymax2=2.0,    height2=3.3,   width2=5.0 , center2=myCenter_g )



MyBoxViolinPlot_1(vector2=c( rowMeans(Average_H3[calss1, c(200:300)]),  rowMeans(Average_H3[calss2,  c(200:300)]), rowMeans(Average_H3[calss3,  c(200:300)]),  rowMeans(Average_H3[calss4,  c(200:300)])  ,  rowMeans(Average_H3[calss5,  c(200:300)]) ),       
                  sampleType2=c(  rep("H3_1", length(calss1) ), rep("H3_2",  length(calss2)),    rep("H3_3",  length(calss3) ),  rep("H3_4",  length(calss4) )  ,  rep("H3_5",  length(calss5) ) ), 
                  sampleRank2=c( "H3_1", "H3_2",  "H3_3",   "H3_4"  ,   "H3_5"  ),     
                  colours2=c("H3_1"="yellow2",  "H3_2"="red2",   "H3_3"="cyan2",    "H3_4"="green2"   ,    "H3_5"="blue2"   ), 
                  path2=subdir_5_part6,     fileName2="6-5-1A-WT-H3-violinPlot",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="H3 signal",   
                  Ymin2=0,   Ymax2=2,    height2=5,   width2=3.5  )


MyAverageLines_1(vector2=c( colMeans(Average_week0[calss1, ]),  colMeans(Average_week0[calss2, ]), colMeans(Average_week0[calss3, ]),  colMeans(Average_week0[calss4, ])  ,  colMeans(Average_week0[calss5, ]) ),   
                 numSample2=5,   
                 sampleType2=c( rep("week0_1", 500 ), rep("week0_2",  500),    rep("week0_3",  500 ),  rep("week0_4",  500 )  ,  rep("week0_5",  500 )  ), 
                 sampleRank2=c( "week0_1", "week0_2",  "week0_3",   "week0_4"  ,   "week0_5"   ),     
                 colours2=c( "week0_1"="yellow2",  "week0_2"="red2",   "week0_3"="cyan2",    "week0_4"="green2"   ,    "week0_5"="blue2"  ), 
                 path2=subdir_5_part6,     fileName2="6-5-1A-WT-week0-averageCurve",  
                 title2="week0 Level",     xLab2="Relative distance (kb)",    yLab2="week0 signal",   
                 Ymin2=0,   Ymax2=2.0,    height2=3.3,   width2=5.0 , center2=myCenter_g )



MyBoxViolinPlot_1(vector2=c( rowMeans(Average_week0[calss1, c(200:300)]),  rowMeans(Average_week0[calss2,  c(200:300)]), rowMeans(Average_week0[calss3,  c(200:300)]),  rowMeans(Average_week0[calss4,  c(200:300)])  ,  rowMeans(Average_week0[calss5,  c(200:300)]) ),       
                  sampleType2=c(  rep("week0_1", length(calss1) ), rep("week0_2",  length(calss2)),    rep("week0_3",  length(calss3) ),  rep("week0_4",  length(calss4) )  ,  rep("week0_5",  length(calss5) ) ), 
                  sampleRank2=c( "week0_1", "week0_2",  "week0_3",   "week0_4"  ,   "week0_5"  ),     
                  colours2=c("week0_1"="yellow2",  "week0_2"="red2",   "week0_3"="cyan2",    "week0_4"="green2"   ,    "week0_5"="blue2"   ), 
                  path2=subdir_5_part6,     fileName2="6-5-1A-WT-week0-violinPlot",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="week0 signal",   
                  Ymin2=0,   Ymax2=2,    height2=5,   width2=3.5  )



MyAverageLines_1(vector2=c( colMeans(Average_week8[calss1, ]),  colMeans(Average_week8[calss2, ]), colMeans(Average_week8[calss3, ]),  colMeans(Average_week8[calss4, ])  ,  colMeans(Average_week8[calss5, ]) ),   
                 numSample2=5,   
                 sampleType2=c( rep("week8_1", 500 ), rep("week8_2",  500),    rep("week8_3",  500 ),  rep("week8_4",  500 )  ,  rep("week8_5",  500 )  ), 
                 sampleRank2=c( "week8_1", "week8_2",  "week8_3",   "week8_4"  ,   "week8_5"   ),     
                 colours2=c( "week8_1"="yellow2",  "week8_2"="red2",   "week8_3"="cyan2",    "week8_4"="green2"   ,    "week8_5"="blue2"  ), 
                 path2=subdir_5_part6,     fileName2="6-5-1A-WT-week8-averageCurve",  
                 title2="week8 Level",     xLab2="Relative distance (kb)",    yLab2="week8 signal",   
                 Ymin2=0,   Ymax2=2.0,    height2=3.3,   width2=5.0 , center2=myCenter_g )



MyBoxViolinPlot_1(vector2=c( rowMeans(Average_week8[calss1, c(200:300)]),  rowMeans(Average_week8[calss2,  c(200:300)]), rowMeans(Average_week8[calss3,  c(200:300)]),  rowMeans(Average_week8[calss4,  c(200:300)])  ,  rowMeans(Average_week8[calss5,  c(200:300)]) ),       
                  sampleType2=c(  rep("week8_1", length(calss1) ), rep("week8_2",  length(calss2)),    rep("week8_3",  length(calss3) ),  rep("week8_4",  length(calss4) )  ,  rep("week8_5",  length(calss5) ) ), 
                  sampleRank2=c( "week8_1", "week8_2",  "week8_3",   "week8_4"  ,   "week8_5"  ),     
                  colours2=c("week8_1"="yellow2",  "week8_2"="red2",   "week8_3"="cyan2",    "week8_4"="green2"   ,    "week8_5"="blue2"   ), 
                  path2=subdir_5_part6,     fileName2="6-5-1A-WT-week8-violinPlot",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="week8 signal",   
                  Ymin2=0,   Ymax2=1,    height2=5,   width2=3.5  )











MyAverageLines_1(vector2=c(colMeans(Average_week0[calss1,]),  colMeans(Average_week1[calss1,]),   colMeans( Average_week2[calss1,]),   
                           colMeans(Average_week4[calss1,]),  colMeans(Average_week6[calss1,]),   colMeans(Average_week8[calss1,])  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_5_part6,     fileName2="6-5-1B-WT-6Samples-averageCurve",  
                 title2="NOL (Lowest)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0[calss2,]),  colMeans(Average_week1[calss2,]),    colMeans(Average_week2[calss2,]),   
                           colMeans(Average_week4[calss2,]),  colMeans(Average_week6[calss2,]),    colMeans(Average_week8[calss2,])  ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_5_part6,     fileName2="6-5-1C-WT-6Samples-averageCurve",  
                 title2="NOL (Low)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0[calss3,]),  colMeans(Average_week1[calss3,]),  colMeans(Average_week2[calss3,]),   
                           colMeans(Average_week4[calss3,]),  colMeans(Average_week6[calss3,]),  colMeans(Average_week8[calss3,] ) ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_5_part6,     fileName2="6-5-1D-WT-6Samples-averageCurve",  
                 title2="NOL (Medium)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0[calss4,]),  colMeans(Average_week1[calss4,]),  colMeans(Average_week2[calss4,]),   
                           colMeans(Average_week4[calss4,]),  colMeans(Average_week6[calss4,]),  colMeans(Average_week8[calss4,] )   ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_5_part6,     fileName2="6-5-1E-WT-6Samples-averageCurve",  
                 title2="NOL (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )


MyAverageLines_1(vector2=c(colMeans(Average_week0[calss5,]),  colMeans(Average_week1[calss5,]),  colMeans(Average_week2[calss5,]),   
                           colMeans(Average_week4[calss5,]),  colMeans(Average_week6[calss5,]),  colMeans(Average_week8[calss5,] )   ),   
                 numSample2=6,   
                 sampleType2=c( rep("week0", numOfColumns1),    rep("week1", numOfColumns1),  rep("week2", numOfColumns1),  
                                rep("week4", numOfColumns1),    rep("week6", numOfColumns1),  rep("week8", numOfColumns1)  ), 
                 sampleRank2=c("week0",  "week1",   "week2",  "week4",  "week6",   "week8" ),     
                 colours2=c( "week0"="red2",   "week1"="cyan2",    "week2"="green2",  "week4"="blue2",  "week6"="purple2",   "week8"="orange2" ), 
                 path2=subdir_5_part6,     fileName2="6-5-1F-WT-6Samples-averageCurve",  
                 title2="NOL (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=1.25,    height2=3.3,   width2=5.0 , center2=myCenter_g )




MyBoxViolinPlot_1(vector2=c( myNTR_7_WT_1[calss1], myNTR_7_WT_1[calss2],  myNTR_7_WT_1[calss3],  myNTR_7_WT_1[calss4],  myNTR_7_WT_1[calss5] ),   
                  sampleType2=c( rep("Lowest", length(calss1)),    rep("Low", length(calss2)),  rep("Medium", length(calss3)),  
                                 rep("High", length(calss4)),      rep("Highest", length(calss5))  ),    
                  sampleRank2=c( "Lowest",  "Low",   "Medium",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",    "Medium"="green2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_5_part6,   fileName2= paste("6-5-1G-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=3.5,   Ymin2=0, Ymax2=0.4)    ## width = 1 + 5*0.5, height=5cm 










dim(Average_week0_EEDheto)   
dim(Average_week0_EEDko)  
dim(Average_week4_EEDheto) 
dim(Average_week4_EEDko)

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss1,]), colMeans( Average_week0_EEDko[calss1,]),      
                           colMeans(Average_week4_EEDheto[calss1,]),  colMeans(Average_week4_EEDko[calss1,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_5_part6,     fileName2="6-5-2A-CKO-4Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss2,]),  colMeans(Average_week0_EEDko[calss2,]),      
                           colMeans(Average_week4_EEDheto[calss2,]),  colMeans(Average_week4_EEDko[calss2,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_5_part6,     fileName2="6-5-2B-CKO-4Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss3,]),  colMeans(Average_week0_EEDko[calss3,]),      
                           colMeans(Average_week4_EEDheto[calss3,]),  colMeans(Average_week4_EEDko[calss3,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_5_part6,     fileName2="6-5-2C-CKO-4Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss4,]),  colMeans(Average_week0_EEDko[calss4,]),      
                           colMeans(Average_week4_EEDheto[calss4,]),  colMeans(Average_week4_EEDko[calss4,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_5_part6,     fileName2="6-5-2D-CKO-4Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_week0_EEDheto[calss5,]),  colMeans(Average_week0_EEDko[calss5,]),      
                           colMeans(Average_week4_EEDheto[calss5,]),  colMeans(Average_week4_EEDko[calss5,])   ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"="red",  "week0_EEDko"="red4",   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_5_part6,     fileName2="6-5-2E-CKO-4Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.65,    height2=3.3,   width2=6.05 , center2=myCenter_g )





MyBoxViolinPlot_1(vector2=c( myNTR_7_EEDheto_2A[calss1], myNTR_7_EEDheto_2A[calss2],  myNTR_7_EEDheto_2A[calss3],  myNTR_7_EEDheto_2A[calss4],  myNTR_7_EEDheto_2A[calss5],
                             myNTR_7_EEDko_2B[calss1],   myNTR_7_EEDko_2B[calss2],    myNTR_7_EEDko_2B[calss3],    myNTR_7_EEDko_2B[calss4],    myNTR_7_EEDko_2B[calss5] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(calss1)),  rep("EEDheto_Low", length(calss2)),  rep("EEDheto_Medium", length(calss3)),  rep("EEDheto_High", length(calss4)),  rep("EEDheto_Highest", length(calss5)), 
                                 rep("EEDko_Lowest", length(calss1)),    rep("EEDko_Low", length(calss2)),    rep("EEDko_Medium", length(calss3)),    rep("EEDko_High", length(calss4)),    rep("EEDko_Highest", length(calss5)) ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",     "EEDko_Medium",    "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2",    "EEDko_Highest"="purple2" ), 
                  path2=subdir_5_part6,   fileName2= paste("6-5-2F-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=0.3)    ## width = 1 + 10*0.5, height=5cm 








dim(Average_banding)  
dim(Average_sham)

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss1,]),  colMeans(Average_sham[calss1,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_5_part6,     fileName2="6-5-3A-TAC-2Samples-averageCurve",  
                 title2="Genes (Lowest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss2,]),  colMeans(Average_sham[calss2,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_5_part6,     fileName2="6-5-3B-TAC-2Samples-averageCurve",  
                 title2="Genes (Low TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss3,]),  colMeans(Average_sham[calss3,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_5_part6,     fileName2="6-5-3C-TAC-2Samples-averageCurve",  
                 title2="Genes (High TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )

MyAverageLines_1(vector2=c(colMeans(Average_banding[calss4,]),  colMeans(Average_sham[calss4,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_5_part6,     fileName2="6-5-3D-TAC-2Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )


MyAverageLines_1(vector2=c(colMeans(Average_banding[calss5,]),  colMeans(Average_sham[calss5,])  ),   
                 numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1),    rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",   "sham"  ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_5_part6,     fileName2="6-5-3E-TAC-2Samples-averageCurve",  
                 title2="Genes (Highest TPM)",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0,   Ymax2=0.5,    height2=3.3,   width2=5.15 , center2=myCenter_g )






MyBoxViolinPlot_1(vector2=c( myNTR_7_banding_3A[calss1], myNTR_7_banding_3A[calss2],  myNTR_7_banding_3A[calss3],  myNTR_7_banding_3A[calss4],  myNTR_7_banding_3A[calss5],
                             myNTR_7_sham_3B[calss1],    myNTR_7_sham_3B[calss2],     myNTR_7_sham_3B[calss3],     myNTR_7_sham_3B[calss4],     myNTR_7_sham_3B[calss5]  ),   
                  sampleType2=c( rep("banding_Lowest", length(calss1)),  rep("banding_Low", length(calss2)),  rep("banding_Medium", length(calss3)),  rep("banding_High", length(calss4)),  rep("banding_Highest", length(calss5)), 
                                 rep("sham_Lowest", length(calss1)),     rep("sham_Low", length(calss2)),     rep("sham_Medium", length(calss3)),     rep("sham_High", length(calss4)),     rep("sham_Highest", length(calss5))  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",   "banding_Medium",  "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",     "sham_Low",      "sham_Medium",     "sham_High",     "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",      "sham_Low"="cyan2",       "sham_Medium"="green2",     "sham_High"="blue2",     "sham_Highest"="purple2"  ), 
                  path2=subdir_5_part6,   fileName2= paste("6-5-3F-TAC-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=1.5)    ## width = 1 + 10*0.5, height=5cm 


sink()  



####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################




