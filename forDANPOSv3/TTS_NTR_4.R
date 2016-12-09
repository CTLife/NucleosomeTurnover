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
## Part  4:  Figures about nucleosome turnover rate (NTR) and compute the correlation between NOL and NTR.
############################################################################







#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################



###############  average all rows for each columns
subdir_1_part4 <- paste(Part4_g,  "/1-NTR-averageAllRows", sep = "")
if( ! file.exists(subdir_1_part4) ) { dir.create(subdir_1_part4) }

# WT 6 samples
sink( file=paste(subdir_1_part4,  "/4-1-1A-LogLinearModel.runLog",      sep = "") )
myNTR_WT_1A  <-  c()
for (i in c(1:numOfColumns1)) {
    vec1 <- c(column_Average_week0[i], column_Average_week1[i], column_Average_week2[i], column_Average_week4[i], column_Average_week6[i], column_Average_week8[i])
    NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log(vec1), file2=paste(subdir_1_part4,  "/4-1-1B-LogLinearModel", sep="")  )   ## log(2.718281828459045) = 1
    myNTR_WT_1A <- c(myNTR_WT_1A, NTR_bin)
}
length(myNTR_WT_1A)
summary(myNTR_WT_1A)
sink()  

##Log-linear model is better than linear model
sink( file=paste(subdir_1_part4,  "/4-1-1C-LinearModel.runLog",      sep = "") )
myNTR_WT_1A2  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(column_Average_week0[i], column_Average_week1[i], column_Average_week2[i], column_Average_week4[i], column_Average_week6[i], column_Average_week8[i])
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=vec1 , file2=paste(subdir_1_part4,  "/4-1-1D-LinearModel", sep="")  )   ##
  myNTR_WT_1A2 <- c(myNTR_WT_1A2, NTR_bin)
}
length(myNTR_WT_1A2)
summary(myNTR_WT_1A2)
sink()  



HalfLife_WT_1A <- log(2)/(myNTR_WT_1A+0.0001)
length(HalfLife_WT_1A)
summary(HalfLife_WT_1A)

myNTR_WT_1A[myNTR_WT_1A<0]  <- 0
myNTR_WT_1A[myNTR_WT_1A>10]  <- 10
HalfLife_WT_1A[HalfLife_WT_1A<0]  <- 0
HalfLife_WT_1A[HalfLife_WT_1A>50]  <- 50


MyAverageLines_1(vector2=myNTR_WT_1A,    numSample2=1,   
                 sampleType2=c( rep("WT", numOfColumns1)  ), 
                 sampleRank2=c( "WT" ),     
                 colours2=c( "WT"="red" ), 
                 path2=subdir_1_part4,     fileName2="4-1-1A-WT-NTR-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.3,    height2=3.3,   width2=5 , center2=myCenter_g  )

MyAverageLines_1(vector2=HalfLife_WT_1A,    numSample2=1,   
                 sampleType2=c( rep("WT", numOfColumns1)  ), 
                 sampleRank2=c( "WT" ),     
                 colours2=c( "WT"="red" ), 
                 path2=subdir_1_part4,     fileName2="4-1-1B-WT-HalfLife-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half life (weeks)",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5 , center2=myCenter_g  )


myNTR_withEquation_1( xAxis2=c(0, 1, 2, 4, 6, 8),  
                      yAxis2=log( c( column_Average_week0[260], column_Average_week1[260], column_Average_week2[260], column_Average_week4[260], column_Average_week6[260], column_Average_week8[260] ) ),                               
                      path2=subdir_1_part4,  fileName2="4-1-1C-WT-NTR-260", 
                      titel2=myTitle_g )
  
myNTR_withEquation_1( xAxis2=c(0, 1, 2, 4, 6, 8),  
                      yAxis2=log( c( mean(column_Average_week0[255:265]), mean(column_Average_week1[255:265]), mean(column_Average_week2[255:265]), 
                                     mean(column_Average_week4[255:265]), mean(column_Average_week6[255:265]), mean(column_Average_week8[255:265]) ) ),                               
                      path2=subdir_1_part4,  fileName2="4-1-1D-WT-NTR-plusOneNucleosome", 
                      titel2=myTitle_g )

myNTR_withEquation_1( xAxis2=c(0, 1, 2, 4, 6, 8),  
                      yAxis2=log( c( mean(column_Average_week0[235:245]), mean(column_Average_week1[235:245]), mean(column_Average_week2[235:245]), 
                                     mean(column_Average_week4[235:245]), mean(column_Average_week6[235:245]), mean(column_Average_week8[235:245]) ) ),                               
                      path2=subdir_1_part4,  fileName2="4-1-1E-WT-NTR-minusOneNucleosome", 
                      titel2=myTitle_g )

  





# CKO 4 samples
sink( file=paste(subdir_1_part4,  "/4-1-2A-runLog.txt",      sep = "") )
myNTR_EEDheto_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(column_Average_week0_EEDheto[i],  column_Average_week4_EEDheto[i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ), file2=paste(subdir_1_part4,  "/4-1-2A-LogLinearModel", sep="")  )
  myNTR_EEDheto_2A <- c(myNTR_EEDheto_2A, NTR_bin)
}
sink()  
length(myNTR_EEDheto_2A)
summary(myNTR_EEDheto_2A)

HalfLife_EEDheto_2A <- log(2)/(myNTR_EEDheto_2A+0.0001)
length(HalfLife_EEDheto_2A)
summary(HalfLife_EEDheto_2A)

myNTR_EEDheto_2A[myNTR_EEDheto_2A<0]  <- 0
myNTR_EEDheto_2A[myNTR_EEDheto_2A>10]  <- 10
HalfLife_EEDheto_2A[HalfLife_EEDheto_2A<0]  <- 0
HalfLife_EEDheto_2A[HalfLife_EEDheto_2A>50]  <- 50

sink( file=paste(subdir_1_part4,  "/4-1-2B-runLog.txt",      sep = "") )
myNTR_EEDko_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(column_Average_week0_EEDko[i],  column_Average_week4_EEDko[i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_1_part4,  "/4-1-2B-LogLinearModel", sep="")  )
  myNTR_EEDko_2B <- c(myNTR_EEDko_2B, NTR_bin)
}
sink()  
length(myNTR_EEDko_2B)
summary(myNTR_EEDko_2B)

HalfLife_EEDko_2B <- log(2)/(myNTR_EEDko_2B+0.0001)
length(HalfLife_EEDko_2B)
summary(HalfLife_EEDko_2B)

myNTR_EEDko_2B[myNTR_EEDko_2B<0]  <- 0
myNTR_EEDko_2B[myNTR_EEDko_2B>10]  <- 10
HalfLife_EEDko_2B[HalfLife_EEDko_2B<0]  <- 0
HalfLife_EEDko_2B[HalfLife_EEDko_2B>50]  <- 50

MyAverageLines_1(vector2=c(myNTR_EEDheto_2A,  myNTR_EEDko_2B),    numSample2=2,   
                 sampleType2=c( rep("EEDheto", numOfColumns1)  ,   rep("EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "EEDheto",  "EEDko" ),     
                 colours2=c( "EEDheto"="red",  "EEDko"="blue"  ), 
                 path2=subdir_1_part4,     fileName2="4-1-2A-CKO-NTR-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.3,    height2=3.3,   width2=6.05 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(HalfLife_EEDheto_2A,  HalfLife_EEDko_2B),    numSample2=2,   
                 sampleType2=c( rep("EEDheto", numOfColumns1)  ,   rep("EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "EEDheto",  "EEDko" ),     
                 colours2=c( "EEDheto"="red",  "EEDko"="blue"  ), 
                 path2=subdir_1_part4,     fileName2="4-1-2B-CKO-HalfLife-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half life (weeks)",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=6.05 , center2=myCenter_g  )

# two different methods for computing NTR based on 2 time points.
#row_Average_one1_EEDheto_NTR2       <-  log( column_Average_week0_EEDheto/column_Average_week4_EEDheto)/4
#summary(row_Average_one1_EEDheto_NTR2)
#summary(myNTR_EEDheto_2A)






# TAC 2 samples
sink( file=paste(subdir_1_part4,  "/4-1-3A-runLog.txt",      sep = "") )
myNTR_banding_3A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(column_Average_week0[i],  column_Average_banding[i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  2.5),  yAxis2=log( vec1 ) , file2=paste(subdir_1_part4,  "/4-1-3A-LogLinearModel", sep="") )
  myNTR_banding_3A <- c(myNTR_banding_3A, NTR_bin)
}
sink()  
length(myNTR_banding_3A)
summary(myNTR_banding_3A)

HalfLife_banding_3A <- log(2)/(myNTR_banding_3A+0.0001)
length(HalfLife_banding_3A)
summary(HalfLife_banding_3A)

myNTR_banding_3A[myNTR_banding_3A<0]  <- 0
myNTR_banding_3A[myNTR_banding_3A>10]  <- 10
HalfLife_banding_3A[HalfLife_banding_3A<0]  <- 0
HalfLife_banding_3A[HalfLife_banding_3A>50]  <- 50

sink( file=paste(subdir_1_part4,  "/4-1-3B-runLog.txt",      sep = "") )
myNTR_sham_3B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(column_Average_week0[i],  column_Average_sham[i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  2.5),  yAxis2=log( vec1 ) , file2=paste(subdir_1_part4,  "/4-1-3B-LogLinearModel", sep="") )
  myNTR_sham_3B <- c(myNTR_sham_3B, NTR_bin)
}
sink()  
length(myNTR_sham_3B)
summary(myNTR_sham_3B)

HalfLife_sham_3B <- log(2)/(myNTR_sham_3B+0.0001)
length(HalfLife_sham_3B)
summary(HalfLife_sham_3B)

myNTR_sham_3B[myNTR_sham_3B<0]  <- 0
myNTR_sham_3B[myNTR_sham_3B>10]  <- 10
HalfLife_sham_3B[HalfLife_sham_3B<0]  <- 0
HalfLife_sham_3B[HalfLife_sham_3B>50]  <- 50

MyAverageLines_1(vector2=c(myNTR_banding_3A,  myNTR_sham_3B),    numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1)  ,   rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",  "sham" ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_1_part4,     fileName2="4-1-3A-TAC-NTR-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.8,    height2=3.3,   width2=5.15 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(HalfLife_banding_3A,  HalfLife_sham_3B),    numSample2=2,   
                 sampleType2=c( rep("banding", numOfColumns1)  ,   rep("sham", numOfColumns1)  ), 
                 sampleRank2=c( "banding",  "sham" ),     
                 colours2=c( "banding"="blue",  "sham"="green4"  ), 
                 path2=subdir_1_part4,     fileName2="4-1-3B-TAC-HalfLife-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half life (weeks)",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.15 , center2=myCenter_g  )

















###############  5 categaries based on rows
subdir_2_part4 <- paste(Part4_g,  "/2-NTR-rows5Classes", sep = "")
if( ! file.exists(subdir_2_part4) ) { dir.create(subdir_2_part4) }

# WT 6 samples
sink( file=paste(subdir_2_part4,  "/4-2-1A-runLog.txt",      sep = "") )
NTR_2_1_WT_1A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0[1, i], reduceRow1_Average_week1[1, i], reduceRow1_Average_week2[1, i], reduceRow1_Average_week4[1, i], reduceRow1_Average_week6[1, i], reduceRow1_Average_week8[1, i])
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-1A-LogLinearModel", sep="") )
  NTR_2_1_WT_1A <- c(NTR_2_1_WT_1A, NTR_bin)
}
sink()  
length(NTR_2_1_WT_1A)
summary(NTR_2_1_WT_1A)
HalfLife_2_1_WT_1A <- log(2)/(NTR_2_1_WT_1A+0.0001)
length(HalfLife_2_1_WT_1A)
summary(HalfLife_2_1_WT_1A)
NTR_2_1_WT_1A[NTR_2_1_WT_1A<0]  <- 0
NTR_2_1_WT_1A[NTR_2_1_WT_1A>10]  <- 10
HalfLife_2_1_WT_1A[HalfLife_2_1_WT_1A<0]  <- 0
HalfLife_2_1_WT_1A[HalfLife_2_1_WT_1A>50]  <- 50

sink( file=paste(subdir_2_part4,  "/4-2-1B-runLog.txt",      sep = "") )
NTR_2_2_WT_1A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0[2, i], reduceRow1_Average_week1[2, i], reduceRow1_Average_week2[2, i], reduceRow1_Average_week4[2, i], reduceRow1_Average_week6[2, i], reduceRow1_Average_week8[2, i])
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-1B-LogLinearModel", sep="") )
  NTR_2_2_WT_1A <- c(NTR_2_2_WT_1A, NTR_bin)
}
sink()  
length(NTR_2_2_WT_1A)
summary(NTR_2_2_WT_1A)
HalfLife_2_2_WT_1A <- log(2)/(NTR_2_2_WT_1A+0.0001)
length(HalfLife_2_2_WT_1A)
summary(HalfLife_2_2_WT_1A)
NTR_2_2_WT_1A[NTR_2_2_WT_1A<0]  <- 0
NTR_2_2_WT_1A[NTR_2_2_WT_1A>10]  <- 10
HalfLife_2_2_WT_1A[HalfLife_2_2_WT_1A<0]  <- 0
HalfLife_2_2_WT_1A[HalfLife_2_2_WT_1A>50]  <- 50

sink( file=paste(subdir_2_part4,  "/4-2-1C-runLog.txt",      sep = "") )
NTR_2_3_WT_1A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0[3, i], reduceRow1_Average_week1[3, i], reduceRow1_Average_week2[3, i], reduceRow1_Average_week4[3, i], reduceRow1_Average_week6[3, i], reduceRow1_Average_week8[3, i])
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-1C-LogLinearModel", sep="") )
  NTR_2_3_WT_1A <- c(NTR_2_3_WT_1A, NTR_bin)
}
sink()  
length(NTR_2_3_WT_1A)
summary(NTR_2_3_WT_1A)
HalfLife_2_3_WT_1A <- log(2)/(NTR_2_3_WT_1A+0.0001)
length(HalfLife_2_3_WT_1A)
summary(HalfLife_2_3_WT_1A)
NTR_2_3_WT_1A[NTR_2_3_WT_1A<0]  <- 0
NTR_2_3_WT_1A[NTR_2_3_WT_1A>10]  <- 10
HalfLife_2_3_WT_1A[HalfLife_2_3_WT_1A<0]  <- 0
HalfLife_2_3_WT_1A[HalfLife_2_3_WT_1A>50]  <- 50

sink( file=paste(subdir_2_part4,  "/4-2-1D-runLog.txt",      sep = "") )
NTR_2_4_WT_1A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0[4, i], reduceRow1_Average_week1[4, i], reduceRow1_Average_week2[4, i], reduceRow1_Average_week4[4, i], reduceRow1_Average_week6[4, i], reduceRow1_Average_week8[4, i])
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-1D-LogLinearModel", sep="") )
  NTR_2_4_WT_1A <- c(NTR_2_4_WT_1A, NTR_bin)
}
sink()  
length(NTR_2_4_WT_1A)
summary(NTR_2_4_WT_1A)
HalfLife_2_4_WT_1A <- log(2)/(NTR_2_4_WT_1A+0.0001)
length(HalfLife_2_4_WT_1A)
summary(HalfLife_2_4_WT_1A)
NTR_2_4_WT_1A[NTR_2_4_WT_1A<0]  <- 0
NTR_2_4_WT_1A[NTR_2_4_WT_1A>10]  <- 10
HalfLife_2_4_WT_1A[HalfLife_2_4_WT_1A<0]   <- 0
HalfLife_2_4_WT_1A[HalfLife_2_4_WT_1A>50]  <- 50

sink( file=paste(subdir_2_part4,  "/4-2-1E-runLog.txt",      sep = "") )
NTR_2_5_WT_1A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0[5, i], reduceRow1_Average_week1[5, i], reduceRow1_Average_week2[5, i], reduceRow1_Average_week4[5, i], reduceRow1_Average_week6[5, i], reduceRow1_Average_week8[5, i])
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-1E-LogLinearModel", sep="") )
  NTR_2_5_WT_1A <- c(NTR_2_5_WT_1A, NTR_bin)
}
sink()  
length(NTR_2_5_WT_1A)
summary(NTR_2_5_WT_1A)
HalfLife_2_5_WT_1A <- log(2)/(NTR_2_5_WT_1A+0.0001)
length(HalfLife_2_5_WT_1A)
summary(HalfLife_2_5_WT_1A)
NTR_2_5_WT_1A[NTR_2_5_WT_1A<0]  <- 0
NTR_2_5_WT_1A[NTR_2_5_WT_1A>10]  <- 10
HalfLife_2_5_WT_1A[HalfLife_2_5_WT_1A<0]  <- 0
HalfLife_2_5_WT_1A[HalfLife_2_5_WT_1A>50]  <- 50

MyAverageLines_1(vector2=c(NTR_2_1_WT_1A,  NTR_2_2_WT_1A,  NTR_2_3_WT_1A,  NTR_2_4_WT_1A,  NTR_2_5_WT_1A),    
                 numSample2=5,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("Medium", numOfColumns1),  rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",  "Medium",  "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",  "Medium"="red4",  "High"="olivedrab1",  "Highest"="olivedrab4" ), 
                 path2=subdir_2_part4,     fileName2="4-2-1A-WT-NTR-average-curve",  
                 title2="All Genes (WT)",     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.3,    height2=3.3,   width2=5.5  , center2=myCenter_g  )

MyAverageLines_1(vector2=c(HalfLife_2_1_WT_1A,  HalfLife_2_2_WT_1A,  HalfLife_2_3_WT_1A,  HalfLife_2_4_WT_1A,  HalfLife_2_5_WT_1A),    
                 numSample2=5,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("Medium", numOfColumns1),  rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",  "Medium",  "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",  "Medium"="red4",  "High"="olivedrab1",  "Highest"="olivedrab4" ), 
                 path2=subdir_2_part4,     fileName2="4-2-1B-WT-HalfLife-average-curve",  
                 title2="All Genes (WT)",     xLab2="Relative distance (kb)",    yLab2="Half life",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.5  , center2=myCenter_g  )











# CKO 4 samples
sink( file=paste(subdir_2_part4,  "/4-2-2A-runLog.txt",      sep = "") )
NTR_2_1_EEDheto_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0_EEDheto[1, i],  reduceRow1_Average_week4_EEDheto[1, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-2A-LogLinearModel", sep="") )
  NTR_2_1_EEDheto_2A <- c(NTR_2_1_EEDheto_2A, NTR_bin)
}
sink()  
length(NTR_2_1_EEDheto_2A)
summary(NTR_2_1_EEDheto_2A)
HalfLife_2_1_EEDheto_2A <- log(2)/(NTR_2_1_EEDheto_2A+0.0001)
length(HalfLife_2_1_EEDheto_2A)
summary(HalfLife_2_1_EEDheto_2A)
NTR_2_1_EEDheto_2A[NTR_2_1_EEDheto_2A<0]  <- 0
NTR_2_1_EEDheto_2A[NTR_2_1_EEDheto_2A>10]  <- 10
HalfLife_2_1_EEDheto_2A[HalfLife_2_1_EEDheto_2A<0]  <- 0
HalfLife_2_1_EEDheto_2A[HalfLife_2_1_EEDheto_2A>50]  <- 50

sink( file=paste(subdir_2_part4,  "/4-2-2B-runLog.txt",      sep = "") )
NTR_2_1_EEDko_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0_EEDko[1, i],  reduceRow1_Average_week4_EEDko[1, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-2B-LogLinearModel", sep=""))
  NTR_2_1_EEDko_2B <- c(NTR_2_1_EEDko_2B, NTR_bin)
}
sink()  
length(NTR_2_1_EEDko_2B)
summary(NTR_2_1_EEDko_2B)
HalfLife_2_1_EEDko_2B <- log(2)/(NTR_2_1_EEDko_2B+0.0001)
length(HalfLife_2_1_EEDko_2B)
summary(HalfLife_2_1_EEDko_2B)
NTR_2_1_EEDko_2B[NTR_2_1_EEDko_2B<0]   <- 0
NTR_2_1_EEDko_2B[NTR_2_1_EEDko_2B>10]  <- 10
HalfLife_2_1_EEDko_2B[HalfLife_2_1_EEDko_2B<0]  <- 0
HalfLife_2_1_EEDko_2B[HalfLife_2_1_EEDko_2B>50]  <- 50




sink( file=paste(subdir_2_part4,  "/4-2-2C-runLog.txt",      sep = "") )
NTR_2_2_EEDheto_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0_EEDheto[2, i],  reduceRow1_Average_week4_EEDheto[2, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-2C-LogLinearModel", sep=""))
  NTR_2_2_EEDheto_2A <- c(NTR_2_2_EEDheto_2A, NTR_bin)
}
sink()  
length(NTR_2_2_EEDheto_2A)
summary(NTR_2_2_EEDheto_2A)
HalfLife_2_2_EEDheto_2A <- log(2)/(NTR_2_2_EEDheto_2A+0.0001)
length(HalfLife_2_2_EEDheto_2A)
summary(HalfLife_2_2_EEDheto_2A)
NTR_2_2_EEDheto_2A[NTR_2_2_EEDheto_2A<0]  <- 0
NTR_2_2_EEDheto_2A[NTR_2_2_EEDheto_2A>10]  <- 10
HalfLife_2_2_EEDheto_2A[HalfLife_2_2_EEDheto_2A<0]  <- 0
HalfLife_2_2_EEDheto_2A[HalfLife_2_2_EEDheto_2A>50]  <- 50

sink( file=paste(subdir_2_part4,  "/4-2-2D-runLog.txt",      sep = "") )
NTR_2_2_EEDko_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0_EEDko[2, i],  reduceRow1_Average_week4_EEDko[2, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-2D-LogLinearModel", sep=""))
  NTR_2_2_EEDko_2B <- c(NTR_2_2_EEDko_2B, NTR_bin)
}
sink()  
length(NTR_2_2_EEDko_2B)
summary(NTR_2_2_EEDko_2B)
HalfLife_2_2_EEDko_2B <- log(2)/(NTR_2_2_EEDko_2B+0.0001)
length(HalfLife_2_2_EEDko_2B)
summary(HalfLife_2_2_EEDko_2B)
NTR_2_2_EEDko_2B[NTR_2_2_EEDko_2B<0]  <- 0
NTR_2_2_EEDko_2B[NTR_2_2_EEDko_2B>10]  <- 10
HalfLife_2_2_EEDko_2B[HalfLife_2_2_EEDko_2B<0]  <- 0
HalfLife_2_2_EEDko_2B[HalfLife_2_2_EEDko_2B>50]  <- 50




sink( file=paste(subdir_2_part4,  "/4-2-2E-runLog.txt",      sep = "") )
NTR_2_3_EEDheto_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0_EEDheto[3, i],  reduceRow1_Average_week4_EEDheto[3, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-2E-LogLinearModel", sep=""))
  NTR_2_3_EEDheto_2A <- c(NTR_2_3_EEDheto_2A, NTR_bin)
}
sink()  
length(NTR_2_3_EEDheto_2A)
summary(NTR_2_3_EEDheto_2A)
HalfLife_2_3_EEDheto_2A <- log(2)/(NTR_2_3_EEDheto_2A+0.0001)
length(HalfLife_2_3_EEDheto_2A)
summary(HalfLife_2_3_EEDheto_2A)
NTR_2_3_EEDheto_2A[NTR_2_3_EEDheto_2A<0]  <- 0
NTR_2_3_EEDheto_2A[NTR_2_3_EEDheto_2A>10]  <- 10
HalfLife_2_3_EEDheto_2A[HalfLife_2_3_EEDheto_2A<0]  <- 0
HalfLife_2_3_EEDheto_2A[HalfLife_2_3_EEDheto_2A>50]  <- 50

sink( file=paste(subdir_2_part4,  "/4-2-2F-runLog.txt",      sep = "") )
NTR_2_3_EEDko_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0_EEDko[3, i],  reduceRow1_Average_week4_EEDko[3, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-2F-LogLinearModel", sep=""))
  NTR_2_3_EEDko_2B <- c(NTR_2_3_EEDko_2B, NTR_bin)
}
sink()  
length(NTR_2_3_EEDko_2B)
summary(NTR_2_3_EEDko_2B)
HalfLife_2_3_EEDko_2B <- log(2)/(NTR_2_3_EEDko_2B+0.0001)
length(HalfLife_2_3_EEDko_2B)
summary(HalfLife_2_3_EEDko_2B)
NTR_2_3_EEDko_2B[NTR_2_3_EEDko_2B<0]  <- 0
NTR_2_3_EEDko_2B[NTR_2_3_EEDko_2B>10]  <- 10
HalfLife_2_3_EEDko_2B[HalfLife_2_3_EEDko_2B<0]  <- 0
HalfLife_2_3_EEDko_2B[HalfLife_2_3_EEDko_2B>50]  <- 50




sink( file=paste(subdir_2_part4,  "/4-2-2G-runLog.txt",      sep = "") )
NTR_2_4_EEDheto_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0_EEDheto[4, i],  reduceRow1_Average_week4_EEDheto[4, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-2G-LogLinearModel", sep=""))
  NTR_2_4_EEDheto_2A <- c(NTR_2_4_EEDheto_2A, NTR_bin)
}
sink()  
length(NTR_2_4_EEDheto_2A)
summary(NTR_2_4_EEDheto_2A)
HalfLife_2_4_EEDheto_2A <- log(2)/(NTR_2_4_EEDheto_2A+0.0001)
length(HalfLife_2_4_EEDheto_2A)
summary(HalfLife_2_4_EEDheto_2A)
NTR_2_4_EEDheto_2A[NTR_2_4_EEDheto_2A<0]  <- 0
NTR_2_4_EEDheto_2A[NTR_2_4_EEDheto_2A>10]  <- 10
HalfLife_2_4_EEDheto_2A[HalfLife_2_4_EEDheto_2A<0]  <- 0
HalfLife_2_4_EEDheto_2A[HalfLife_2_4_EEDheto_2A>50]  <- 50

sink( file=paste(subdir_2_part4,  "/4-2-2H-runLog.txt",      sep = "") )
NTR_2_4_EEDko_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0_EEDko[4, i],  reduceRow1_Average_week4_EEDko[4, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-2H-LogLinearModel", sep=""))
  NTR_2_4_EEDko_2B <- c(NTR_2_4_EEDko_2B, NTR_bin)
}
sink()  
length(NTR_2_4_EEDko_2B)
summary(NTR_2_4_EEDko_2B)
HalfLife_2_4_EEDko_2B <- log(2)/(NTR_2_4_EEDko_2B+0.0001)
length(HalfLife_2_4_EEDko_2B)
summary(HalfLife_2_4_EEDko_2B)
NTR_2_4_EEDko_2B[NTR_2_4_EEDko_2B<0]  <- 0
NTR_2_4_EEDko_2B[NTR_2_4_EEDko_2B<10]  <- 10
HalfLife_2_4_EEDko_2B[HalfLife_2_4_EEDko_2B<0]  <- 0
HalfLife_2_4_EEDko_2B[HalfLife_2_4_EEDko_2B<50]  <- 50






sink( file=paste(subdir_2_part4,  "/4-2-2I-runLog.txt",      sep = "") )
NTR_2_5_EEDheto_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0_EEDheto[5, i],  reduceRow1_Average_week4_EEDheto[5, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ), file2=paste(subdir_2_part4,  "/4-2-2I-LogLinearModel", sep="") )
  NTR_2_5_EEDheto_2A <- c(NTR_2_5_EEDheto_2A, NTR_bin)
}
sink()  
length(NTR_2_5_EEDheto_2A)
summary(NTR_2_5_EEDheto_2A)
HalfLife_2_5_EEDheto_2A <- log(2)/(NTR_2_5_EEDheto_2A+0.0001)
length(HalfLife_2_5_EEDheto_2A)
summary(HalfLife_2_5_EEDheto_2A)
NTR_2_5_EEDheto_2A[NTR_2_5_EEDheto_2A<0]  <- 0
NTR_2_5_EEDheto_2A[NTR_2_5_EEDheto_2A>10]  <- 10
HalfLife_2_5_EEDheto_2A[HalfLife_2_5_EEDheto_2A<0]  <- 0
HalfLife_2_5_EEDheto_2A[HalfLife_2_5_EEDheto_2A>50]  <- 50

sink( file=paste(subdir_2_part4,  "/4-2-2J-runLog.txt",      sep = "") )
NTR_2_5_EEDko_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0_EEDko[5, i],  reduceRow1_Average_week4_EEDko[5, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-2J-LogLinearModel", sep=""))
  NTR_2_5_EEDko_2B <- c(NTR_2_5_EEDko_2B, NTR_bin)
}
sink()  
length(NTR_2_5_EEDko_2B)
summary(NTR_2_5_EEDko_2B)
HalfLife_2_5_EEDko_2B <- log(2)/(NTR_2_5_EEDko_2B+0.0001)
length(HalfLife_2_5_EEDko_2B)
summary(HalfLife_2_5_EEDko_2B)
NTR_2_5_EEDko_2B[NTR_2_5_EEDko_2B<0]  <- 0
NTR_2_5_EEDko_2B[NTR_2_5_EEDko_2B>10]  <- 10
HalfLife_2_5_EEDko_2B[HalfLife_2_5_EEDko_2B<0]  <- 0
HalfLife_2_5_EEDko_2B[HalfLife_2_5_EEDko_2B>50]  <- 50


MyAverageLines_1(vector2=c(NTR_2_1_EEDheto_2A,  NTR_2_2_EEDheto_2A,  NTR_2_3_EEDheto_2A,  NTR_2_4_EEDheto_2A,  NTR_2_5_EEDheto_2A),    
                 numSample2=5,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("Medium", numOfColumns1),  rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",  "Medium",  "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",  "Medium"="red4",  "High"="olivedrab1",  "Highest"="olivedrab4" ), 
                 path2=subdir_2_part4,     fileName2="4-2-2A-EEDheto-NTR-average-curve",  
                 title2="All Genes (EEDheto)",     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.3,    height2=3.3,   width2=5.5  , center2=myCenter_g  )

MyAverageLines_1(vector2=c(NTR_2_1_EEDko_2B,  NTR_2_2_EEDko_2B,  NTR_2_3_EEDko_2B,  NTR_2_4_EEDko_2B,  NTR_2_5_EEDko_2B),    
                 numSample2=5,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("Medium", numOfColumns1),  rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",  "Medium",  "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",  "Medium"="red4",  "High"="olivedrab1",  "Highest"="olivedrab4" ), 
                 path2=subdir_2_part4,     fileName2="4-2-2B-EEDko-NTR-average-curve",  
                 title2="All Genes (EEDko)",     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.3,    height2=3.3,   width2=5.5  , center2=myCenter_g  )




MyAverageLines_1(vector2=c(HalfLife_2_1_EEDheto_2A,  HalfLife_2_2_EEDheto_2A,  HalfLife_2_3_EEDheto_2A,  HalfLife_2_4_EEDheto_2A,  HalfLife_2_5_EEDheto_2A),    
                 numSample2=5,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("Medium", numOfColumns1),  rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",  "Medium",  "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",  "Medium"="red4",  "High"="olivedrab1",  "Highest"="olivedrab4" ), 
                 path2=subdir_2_part4,     fileName2="4-2-2C-EEDheto-HalfLife-average-curve",  
                 title2="All Genes (EEDheto)",     xLab2="Relative distance (kb)",    yLab2="Half Life",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.5  , center2=myCenter_g  )

MyAverageLines_1(vector2=c(HalfLife_2_1_EEDko_2B,  HalfLife_2_2_EEDko_2B,  HalfLife_2_3_EEDko_2B,  HalfLife_2_4_EEDko_2B,  HalfLife_2_5_EEDko_2B),    
                 numSample2=5,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("Medium", numOfColumns1),  rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",  "Medium",  "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",  "Medium"="red4",  "High"="olivedrab1",  "Highest"="olivedrab4" ), 
                 path2=subdir_2_part4,     fileName2="4-2-2D-EEDko-HalfLife-average-curve",  
                 title2="All Genes (EEDko)",     xLab2="Relative distance (kb)",    yLab2="Half Life",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.5   , center2=myCenter_g )









# TAC 2 samples
sink( file=paste(subdir_2_part4,  "/4-2-3A-runLog.txt",      sep = "") )
NTR_2_1_banding_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0[1, i],  reduceRow1_Average_banding[1, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ), file2=paste(subdir_2_part4,  "/4-2-3A-LogLinearModel", sep="") )
  NTR_2_1_banding_2A <- c(NTR_2_1_banding_2A, NTR_bin)
}
sink()  
length(NTR_2_1_banding_2A)
summary(NTR_2_1_banding_2A)
HalfLife_2_1_banding_2A <- log(2)/(NTR_2_1_banding_2A+0.0001)
length(HalfLife_2_1_banding_2A)
summary(HalfLife_2_1_banding_2A)
NTR_2_1_banding_2A[NTR_2_1_banding_2A<0]  <- 0
NTR_2_1_banding_2A[NTR_2_1_banding_2A>10]  <- 10
HalfLife_2_1_banding_2A[HalfLife_2_1_banding_2A<0]  <- 0
HalfLife_2_1_banding_2A[HalfLife_2_1_banding_2A>50]  <- 50

sink( file=paste(subdir_2_part4,  "/4-2-3B-runLog.txt",      sep = "") )
NTR_2_1_sham_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0[1, i],  reduceRow1_Average_sham[1, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-3B-LogLinearModel", sep=""))
  NTR_2_1_sham_2B <- c(NTR_2_1_sham_2B, NTR_bin)
}
sink()  
length(NTR_2_1_sham_2B)
summary(NTR_2_1_sham_2B)
HalfLife_2_1_sham_2B <- log(2)/(NTR_2_1_sham_2B+0.0001)
length(HalfLife_2_1_sham_2B)
summary(HalfLife_2_1_sham_2B)
NTR_2_1_sham_2B[NTR_2_1_sham_2B<0]  <- 0
NTR_2_1_sham_2B[NTR_2_1_sham_2B>10]  <- 10
HalfLife_2_1_sham_2B[HalfLife_2_1_sham_2B<0]  <- 0
HalfLife_2_1_sham_2B[HalfLife_2_1_sham_2B>50]  <- 50




sink( file=paste(subdir_2_part4,  "/4-2-3C-runLog.txt",      sep = "") )
NTR_2_2_banding_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0[2, i],  reduceRow1_Average_banding[2, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-3C-LogLinearModel", sep=""))
  NTR_2_2_banding_2A <- c(NTR_2_2_banding_2A, NTR_bin)
}
sink()  
length(NTR_2_2_banding_2A)
summary(NTR_2_2_banding_2A)
HalfLife_2_2_banding_2A <- log(2)/(NTR_2_2_banding_2A+0.0001)
length(HalfLife_2_2_banding_2A)
summary(HalfLife_2_2_banding_2A)
NTR_2_2_banding_2A[NTR_2_2_banding_2A<0]  <- 0
NTR_2_2_banding_2A[NTR_2_2_banding_2A>10]  <- 10
HalfLife_2_2_banding_2A[HalfLife_2_2_banding_2A<0]  <- 0
HalfLife_2_2_banding_2A[HalfLife_2_2_banding_2A>50]  <- 50

sink( file=paste(subdir_2_part4,  "/4-2-3D-runLog.txt",      sep = "") )
NTR_2_2_sham_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0[2, i],  reduceRow1_Average_sham[2, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-3D-LogLinearModel", sep=""))
  NTR_2_2_sham_2B <- c(NTR_2_2_sham_2B, NTR_bin)
}
sink()  
length(NTR_2_2_sham_2B)
summary(NTR_2_2_sham_2B)
HalfLife_2_2_sham_2B <- log(2)/(NTR_2_2_sham_2B+0.0001)
length(HalfLife_2_2_sham_2B)
summary(HalfLife_2_2_sham_2B)
NTR_2_2_sham_2B[NTR_2_2_sham_2B<0]  <- 0
NTR_2_2_sham_2B[NTR_2_2_sham_2B>10]  <- 10
HalfLife_2_2_sham_2B[HalfLife_2_2_sham_2B<0]  <- 0
HalfLife_2_2_sham_2B[HalfLife_2_2_sham_2B>50]  <- 50




sink( file=paste(subdir_2_part4,  "/4-2-3E-runLog.txt",      sep = "") )
NTR_2_3_banding_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0[3, i],  reduceRow1_Average_banding[3, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-3E-LogLinearModel", sep=""))
  NTR_2_3_banding_2A <- c(NTR_2_3_banding_2A, NTR_bin)
}
sink()  
length(NTR_2_3_banding_2A)
summary(NTR_2_3_banding_2A)
HalfLife_2_3_banding_2A <- log(2)/(NTR_2_3_banding_2A+0.0001)
length(HalfLife_2_3_banding_2A)
summary(HalfLife_2_3_banding_2A)
NTR_2_3_banding_2A[NTR_2_3_banding_2A<0]  <- 0
NTR_2_3_banding_2A[NTR_2_3_banding_2A>10]  <- 10
HalfLife_2_3_banding_2A[HalfLife_2_3_banding_2A<0]  <- 0
HalfLife_2_3_banding_2A[HalfLife_2_3_banding_2A>50]  <- 50

sink( file=paste(subdir_2_part4,  "/4-2-3F-runLog.txt",      sep = "") )
NTR_2_3_sham_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0[3, i],  reduceRow1_Average_sham[3, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-3F-LogLinearModel", sep="") )
  NTR_2_3_sham_2B <- c(NTR_2_3_sham_2B, NTR_bin)
}
sink()  
length(NTR_2_3_sham_2B)
summary(NTR_2_3_sham_2B)
HalfLife_2_3_sham_2B <- log(2)/(NTR_2_3_sham_2B+0.0001)
length(HalfLife_2_3_sham_2B)
summary(HalfLife_2_3_sham_2B)
NTR_2_3_sham_2B[NTR_2_3_sham_2B<0]  <- 0
NTR_2_3_sham_2B[NTR_2_3_sham_2B>10]  <- 10
HalfLife_2_3_sham_2B[HalfLife_2_3_sham_2B<0]  <- 0
HalfLife_2_3_sham_2B[HalfLife_2_3_sham_2B>50]  <- 50




sink( file=paste(subdir_2_part4,  "/4-2-3G-runLog.txt",      sep = "") )
NTR_2_4_banding_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0[4, i],  reduceRow1_Average_banding[4, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-3G-LogLinearModel", sep=""))
  NTR_2_4_banding_2A <- c(NTR_2_4_banding_2A, NTR_bin)
}
sink()  
length(NTR_2_4_banding_2A)
summary(NTR_2_4_banding_2A)
HalfLife_2_4_banding_2A <- log(2)/(NTR_2_4_banding_2A+0.0001)
length(HalfLife_2_4_banding_2A)
summary(HalfLife_2_4_banding_2A)
NTR_2_4_banding_2A[NTR_2_4_banding_2A<0]  <- 0
NTR_2_4_banding_2A[NTR_2_4_banding_2A>10]  <- 10
HalfLife_2_4_banding_2A[HalfLife_2_4_banding_2A<0]  <- 0
HalfLife_2_4_banding_2A[HalfLife_2_4_banding_2A>50]  <- 50

sink( file=paste(subdir_2_part4,  "/4-2-3H-runLog.txt",      sep = "") )
NTR_2_4_sham_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0[4, i],  reduceRow1_Average_sham[4, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-3H-LogLinearModel", sep=""))
  NTR_2_4_sham_2B <- c(NTR_2_4_sham_2B, NTR_bin)
}
sink()  
length(NTR_2_4_sham_2B)
summary(NTR_2_4_sham_2B)
HalfLife_2_4_sham_2B <- log(2)/(NTR_2_4_sham_2B+0.0001)
length(HalfLife_2_4_sham_2B)
summary(HalfLife_2_4_sham_2B)
NTR_2_4_sham_2B[NTR_2_4_sham_2B<0]  <- 0
NTR_2_4_sham_2B[NTR_2_4_sham_2B>10]  <- 10
HalfLife_2_4_sham_2B[HalfLife_2_4_sham_2B<0]  <- 0
HalfLife_2_4_sham_2B[HalfLife_2_4_sham_2B>50]  <- 50






sink( file=paste(subdir_2_part4,  "/4-2-3I-runLog.txt",      sep = "") )
NTR_2_5_banding_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0[5, i],  reduceRow1_Average_banding[5, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_2_part4,  "/4-2-3I-LogLinearModel", sep=""))
  NTR_2_5_banding_2A <- c(NTR_2_5_banding_2A, NTR_bin)
}
sink()  
length(NTR_2_5_banding_2A)
summary(NTR_2_5_banding_2A)
HalfLife_2_5_banding_2A <- log(2)/(NTR_2_5_banding_2A+0.0001)
length(HalfLife_2_5_banding_2A)
summary(HalfLife_2_5_banding_2A)
NTR_2_5_banding_2A[NTR_2_5_banding_2A<0]  <- 0
NTR_2_5_banding_2A[NTR_2_5_banding_2A>10]  <- 10
HalfLife_2_5_banding_2A[HalfLife_2_5_banding_2A<0]  <- 0
HalfLife_2_5_banding_2A[HalfLife_2_5_banding_2A>50]  <- 50

sink( file=paste(subdir_2_part4,  "/4-2-3J-runLog.txt",      sep = "") )
NTR_2_5_sham_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow1_Average_week0[5, i],  reduceRow1_Average_sham[5, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ), file2=paste(subdir_2_part4,  "/4-2-3J-LogLinearModel", sep="") )
  NTR_2_5_sham_2B <- c(NTR_2_5_sham_2B, NTR_bin)
}
sink()  
length(NTR_2_5_sham_2B)
summary(NTR_2_5_sham_2B)
HalfLife_2_5_sham_2B <- log(2)/(NTR_2_5_sham_2B+0.0001)
length(HalfLife_2_5_sham_2B)
summary(HalfLife_2_5_sham_2B)
NTR_2_5_sham_2B[NTR_2_5_sham_2B<0]  <- 0
NTR_2_5_sham_2B[NTR_2_5_sham_2B>10]  <- 10
HalfLife_2_5_sham_2B[HalfLife_2_5_sham_2B<0]  <- 0
HalfLife_2_5_sham_2B[HalfLife_2_5_sham_2B>50]  <- 50


MyAverageLines_1(vector2=c(NTR_2_1_banding_2A,  NTR_2_2_banding_2A,  NTR_2_3_banding_2A,  NTR_2_4_banding_2A,  NTR_2_5_banding_2A),    
                 numSample2=5,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("Medium", numOfColumns1),  rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",  "Medium",  "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",  "Medium"="red4",  "High"="olivedrab1",  "Highest"="olivedrab4" ), 
                 path2=subdir_2_part4,     fileName2="4-2-3A-banding-NTR-average-curve",  
                 title2="All Genes (banding)",     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.8,    height2=3.3,   width2=5.5  , center2=myCenter_g  )

MyAverageLines_1(vector2=c(NTR_2_1_sham_2B,  NTR_2_2_sham_2B,  NTR_2_3_sham_2B,  NTR_2_4_sham_2B,  NTR_2_5_sham_2B),    
                 numSample2=5,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("Medium", numOfColumns1),  rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",  "Medium",  "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",  "Medium"="red4",  "High"="olivedrab1",  "Highest"="olivedrab4" ), 
                 path2=subdir_2_part4,     fileName2="4-2-3B-sham-NTR-average-curve",  
                 title2="All Genes (sham)",     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.8,    height2=3.3,   width2=5.5  , center2=myCenter_g  )


MyAverageLines_1(vector2=c(HalfLife_2_1_banding_2A,  HalfLife_2_2_banding_2A,  HalfLife_2_3_banding_2A,  HalfLife_2_4_banding_2A,  HalfLife_2_5_banding_2A),    
                 numSample2=5,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("Medium", numOfColumns1),  rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",  "Medium",  "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",  "Medium"="red4",  "High"="olivedrab1",  "Highest"="olivedrab4" ), 
                 path2=subdir_2_part4,     fileName2="4-2-3C-banding-HalfLife-average-curve",  
                 title2="All Genes (banding)",     xLab2="Relative distance (kb)",    yLab2="Half Life",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.5  , center2=myCenter_g  )

MyAverageLines_1(vector2=c(HalfLife_2_1_sham_2B,  HalfLife_2_2_sham_2B,  HalfLife_2_3_sham_2B,  HalfLife_2_4_sham_2B,  HalfLife_2_5_sham_2B),    
                 numSample2=5,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("Medium", numOfColumns1),  rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",  "Medium",  "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",  "Medium"="red4",  "High"="olivedrab1",  "Highest"="olivedrab4" ), 
                 path2=subdir_2_part4,     fileName2="4-2-3D-sham-HalfLife-average-curve",  
                 title2="All Genes (sham)",     xLab2="Relative distance (kb)",    yLab2="Half Life",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.5  , center2=myCenter_g  )
























###############  4 categaries based on rows
subdir_3_part4 <- paste(Part4_g,  "/3-NTR-rows4Classes", sep = "")
if( ! file.exists(subdir_3_part4) ) { dir.create(subdir_3_part4) }


# WT 6 samples
sink( file=paste(subdir_3_part4,  "/4-3-1A-runLog.txt",      sep = "") )
NTR_3_1_WT_1A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0[1, i], reduceRow2_Average_week1[1, i], reduceRow2_Average_week2[1, i], reduceRow2_Average_week4[1, i], reduceRow2_Average_week6[1, i], reduceRow2_Average_week8[1, i])
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log( vec1 ) , file2=paste(subdir_3_part4,  "/4-3-1A-LogLinearModel", sep=""))
  NTR_3_1_WT_1A <- c(NTR_3_1_WT_1A, NTR_bin)
}
sink()  
length(NTR_3_1_WT_1A)
summary(NTR_3_1_WT_1A)
HalfLife_3_1_WT_1A <- log(2)/(NTR_3_1_WT_1A+0.0001)
length(HalfLife_3_1_WT_1A)
summary(HalfLife_3_1_WT_1A)
NTR_3_1_WT_1A[NTR_3_1_WT_1A<0]  <- 0
NTR_3_1_WT_1A[NTR_3_1_WT_1A>10]  <- 10
HalfLife_3_1_WT_1A[HalfLife_3_1_WT_1A<0]  <- 0
HalfLife_3_1_WT_1A[HalfLife_3_1_WT_1A>50]  <- 50


sink( file=paste(subdir_3_part4,  "/4-3-1B-runLog.txt",      sep = "") )
NTR_3_2_WT_1A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0[2, i], reduceRow2_Average_week1[2, i], reduceRow2_Average_week2[2, i], reduceRow2_Average_week4[2, i], reduceRow2_Average_week6[2, i], reduceRow2_Average_week8[2, i])
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log( vec1 ) , file2=paste(subdir_3_part4,  "/4-3-1B-LogLinearModel", sep=""))
  NTR_3_2_WT_1A <- c(NTR_3_2_WT_1A, NTR_bin)
}
sink()  
length(NTR_3_2_WT_1A)
summary(NTR_3_2_WT_1A)
HalfLife_3_2_WT_1A <- log(2)/(NTR_3_2_WT_1A+0.0001)
length(HalfLife_3_2_WT_1A)
summary(HalfLife_3_2_WT_1A)
NTR_3_2_WT_1A[NTR_3_2_WT_1A<0]  <- 0
NTR_3_2_WT_1A[NTR_3_2_WT_1A>10]  <- 10
HalfLife_3_2_WT_1A[HalfLife_3_2_WT_1A<0]  <- 0
HalfLife_3_2_WT_1A[HalfLife_3_2_WT_1A>50]  <- 50


sink( file=paste(subdir_3_part4,  "/4-3-1C-runLog.txt",      sep = "") )
NTR_3_3_WT_1A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0[3, i], reduceRow2_Average_week1[3, i], reduceRow2_Average_week2[3, i], reduceRow2_Average_week4[3, i], reduceRow2_Average_week6[3, i], reduceRow2_Average_week8[3, i])
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log( vec1 ) , file2=paste(subdir_3_part4,  "/4-3-1C-LogLinearModel", sep="") )
  NTR_3_3_WT_1A <- c(NTR_3_3_WT_1A, NTR_bin)
}
sink()  
length(NTR_3_3_WT_1A)
summary(NTR_3_3_WT_1A)
HalfLife_3_3_WT_1A <- log(2)/(NTR_3_3_WT_1A+0.0001)
length(HalfLife_3_3_WT_1A)
summary(HalfLife_3_3_WT_1A)
NTR_3_3_WT_1A[NTR_3_3_WT_1A<0]  <- 0
NTR_3_3_WT_1A[NTR_3_3_WT_1A>10]  <- 10
HalfLife_3_3_WT_1A[HalfLife_3_3_WT_1A<0]  <- 0
HalfLife_3_3_WT_1A[HalfLife_3_3_WT_1A>50]  <- 50


sink( file=paste(subdir_3_part4,  "/4-3-1D-runLog.txt",      sep = "") )
NTR_3_4_WT_1A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0[4, i], reduceRow2_Average_week1[4, i], reduceRow2_Average_week2[4, i], reduceRow2_Average_week4[4, i], reduceRow2_Average_week6[4, i], reduceRow2_Average_week8[4, i])
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log( vec1 ) , file2=paste(subdir_3_part4,  "/4-3-1D-LogLinearModel", sep=""))
  NTR_3_4_WT_1A <- c(NTR_3_4_WT_1A, NTR_bin)
}
sink()  
length(NTR_3_4_WT_1A)
summary(NTR_3_4_WT_1A)
HalfLife_3_4_WT_1A <- log(2)/(NTR_3_4_WT_1A+0.0001)
length(HalfLife_3_4_WT_1A)
summary(HalfLife_3_4_WT_1A)
NTR_3_4_WT_1A[NTR_3_4_WT_1A<0]  <- 0
NTR_3_4_WT_1A[NTR_3_4_WT_1A>10]  <- 10
HalfLife_3_4_WT_1A[HalfLife_3_4_WT_1A<0]   <- 0
HalfLife_3_4_WT_1A[HalfLife_3_4_WT_1A>50]  <- 50


MyAverageLines_1(vector2=c(NTR_3_1_WT_1A,  NTR_3_2_WT_1A,  NTR_3_3_WT_1A,  NTR_3_4_WT_1A),    
                 numSample2=4,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                               rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",    "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",   "High"="red4",  "Highest"="olivedrab4" ), 
                 path2=subdir_3_part4,     fileName2="4-3-1A-WT-NTR-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.3,    height2=3.3,   width2=5.5  , center2=myCenter_g )


MyAverageLines_1(vector2=c(HalfLife_3_1_WT_1A,  HalfLife_3_2_WT_1A,  HalfLife_3_3_WT_1A,  HalfLife_3_4_WT_1A),    
                 numSample2=4,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",   "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",   "High"="red4",  "Highest"="olivedrab4" ), 
                 path2=subdir_3_part4,     fileName2="4-3-1B-WT-HalfLife-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half life",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.5   , center2=myCenter_g)





# CKO 4 samples
sink( file=paste(subdir_3_part4,  "/4-3-2A-runLog.txt",      sep = "") )
NTR_3_1_EEDheto_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0_EEDheto[1, i],  reduceRow2_Average_week4_EEDheto[1, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_3_part4,  "/4-3-2A-LogLinearModel", sep="") )
  NTR_3_1_EEDheto_2A <- c(NTR_3_1_EEDheto_2A, NTR_bin)
}
sink()  
length(NTR_3_1_EEDheto_2A)
summary(NTR_3_1_EEDheto_2A)
HalfLife_3_1_EEDheto_2A <- log(2)/(NTR_3_1_EEDheto_2A+0.0001)
length(HalfLife_3_1_EEDheto_2A)
summary(HalfLife_3_1_EEDheto_2A)
NTR_3_1_EEDheto_2A[NTR_3_1_EEDheto_2A<0]  <- 0
NTR_3_1_EEDheto_2A[NTR_3_1_EEDheto_2A>10]  <- 10
HalfLife_3_1_EEDheto_2A[HalfLife_3_1_EEDheto_2A<0]  <- 0
HalfLife_3_1_EEDheto_2A[HalfLife_3_1_EEDheto_2A>50]  <- 50

sink( file=paste(subdir_3_part4,  "/4-3-2B-runLog.txt",      sep = "") )
NTR_3_1_EEDko_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0_EEDko[1, i],  reduceRow2_Average_week4_EEDko[1, i]   )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_3_part4,  "/4-3-2B-LogLinearModel", sep=""))
  NTR_3_1_EEDko_2B <- c(NTR_3_1_EEDko_2B, NTR_bin)
}
sink()  
length(NTR_3_1_EEDko_2B)
summary(NTR_3_1_EEDko_2B)
HalfLife_3_1_EEDko_2B <- log(2)/(NTR_3_1_EEDko_2B+0.0001)
length(HalfLife_3_1_EEDko_2B)
summary(HalfLife_3_1_EEDko_2B)
NTR_3_1_EEDko_2B[NTR_3_1_EEDko_2B<0]   <- 0
NTR_3_1_EEDko_2B[NTR_3_1_EEDko_2B>10]  <- 10
HalfLife_3_1_EEDko_2B[HalfLife_3_1_EEDko_2B<0]  <- 0
HalfLife_3_1_EEDko_2B[HalfLife_3_1_EEDko_2B>50]  <- 50




sink( file=paste(subdir_3_part4,  "/4-3-2C-runLog.txt",      sep = "") )
NTR_3_2_EEDheto_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0_EEDheto[2, i],  reduceRow2_Average_week4_EEDheto[2, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_3_part4,  "/4-3-2C-LogLinearModel", sep="") )
  NTR_3_2_EEDheto_2A <- c(NTR_3_2_EEDheto_2A, NTR_bin)
}
sink()  
length(NTR_3_2_EEDheto_2A)
summary(NTR_3_2_EEDheto_2A)
HalfLife_3_2_EEDheto_2A <- log(2)/(NTR_3_2_EEDheto_2A+0.0001)
length(HalfLife_3_2_EEDheto_2A)
summary(HalfLife_3_2_EEDheto_2A)
NTR_3_2_EEDheto_2A[NTR_3_2_EEDheto_2A<0]  <- 0
NTR_3_2_EEDheto_2A[NTR_3_2_EEDheto_2A>10]  <- 10
HalfLife_3_2_EEDheto_2A[HalfLife_3_2_EEDheto_2A<0]  <- 0
HalfLife_3_2_EEDheto_2A[HalfLife_3_2_EEDheto_2A>50]  <- 50

sink( file=paste(subdir_3_part4,  "/4-3-2D-runLog.txt",      sep = "") )
NTR_3_2_EEDko_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0_EEDko[2, i],  reduceRow2_Average_week4_EEDko[2, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_3_part4,  "/4-3-2D-LogLinearModel", sep="") )
  NTR_3_2_EEDko_2B <- c(NTR_3_2_EEDko_2B, NTR_bin)
}
sink()  
length(NTR_3_2_EEDko_2B)
summary(NTR_3_2_EEDko_2B)
HalfLife_3_2_EEDko_2B <- log(2)/(NTR_3_2_EEDko_2B+0.0001)
length(HalfLife_3_2_EEDko_2B)
summary(HalfLife_3_2_EEDko_2B)
NTR_3_2_EEDko_2B[NTR_3_2_EEDko_2B<0]  <- 0
NTR_3_2_EEDko_2B[NTR_3_2_EEDko_2B>10]  <- 10
HalfLife_3_2_EEDko_2B[HalfLife_3_2_EEDko_2B<0]  <- 0
HalfLife_3_2_EEDko_2B[HalfLife_3_2_EEDko_2B>50]  <- 50




sink( file=paste(subdir_3_part4,  "/4-3-2E-runLog.txt",      sep = "") )
NTR_3_3_EEDheto_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0_EEDheto[3, i],  reduceRow2_Average_week4_EEDheto[3, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_3_part4,  "/4-3-2E-LogLinearModel", sep="")  )
  NTR_3_3_EEDheto_2A <- c(NTR_3_3_EEDheto_2A, NTR_bin)
}
sink()  
length(NTR_3_3_EEDheto_2A)
summary(NTR_3_3_EEDheto_2A)
HalfLife_3_3_EEDheto_2A <- log(2)/(NTR_3_3_EEDheto_2A+0.0001)
length(HalfLife_3_3_EEDheto_2A)
summary(HalfLife_3_3_EEDheto_2A)
NTR_3_3_EEDheto_2A[NTR_3_3_EEDheto_2A<0]  <- 0
NTR_3_3_EEDheto_2A[NTR_3_3_EEDheto_2A>10]  <- 10
HalfLife_3_3_EEDheto_2A[HalfLife_3_3_EEDheto_2A<0]  <- 0
HalfLife_3_3_EEDheto_2A[HalfLife_3_3_EEDheto_2A>50]  <- 50

sink( file=paste(subdir_3_part4,  "/4-3-2F-runLog.txt",      sep = "") )
NTR_3_3_EEDko_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0_EEDko[3, i],  reduceRow2_Average_week4_EEDko[3, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_3_part4,  "/4-3-2F-LogLinearModel", sep="")  )
  NTR_3_3_EEDko_2B <- c(NTR_3_3_EEDko_2B, NTR_bin)
}
sink()  
length(NTR_3_3_EEDko_2B)
summary(NTR_3_3_EEDko_2B)
HalfLife_3_3_EEDko_2B <- log(2)/(NTR_3_3_EEDko_2B+0.0001)
length(HalfLife_3_3_EEDko_2B)
summary(HalfLife_3_3_EEDko_2B)
NTR_3_3_EEDko_2B[NTR_3_3_EEDko_2B<0]  <- 0
NTR_3_3_EEDko_2B[NTR_3_3_EEDko_2B>10]  <- 10
HalfLife_3_3_EEDko_2B[HalfLife_3_3_EEDko_2B<0]  <- 0
HalfLife_3_3_EEDko_2B[HalfLife_3_3_EEDko_2B>50]  <- 50




sink( file=paste(subdir_3_part4,  "/4-3-2G-runLog.txt",      sep = "") )
NTR_3_4_EEDheto_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0_EEDheto[4, i],  reduceRow2_Average_week4_EEDheto[4, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_3_part4,  "/4-3-2G-LogLinearModel", sep="") )
  NTR_3_4_EEDheto_2A <- c(NTR_3_4_EEDheto_2A, NTR_bin)
}
sink()  
length(NTR_3_4_EEDheto_2A)
summary(NTR_3_4_EEDheto_2A)
HalfLife_3_4_EEDheto_2A <- log(2)/(NTR_3_4_EEDheto_2A+0.0001)
length(HalfLife_3_4_EEDheto_2A)
summary(HalfLife_3_4_EEDheto_2A)
NTR_3_4_EEDheto_2A[NTR_3_4_EEDheto_2A<0]  <- 0
NTR_3_4_EEDheto_2A[NTR_3_4_EEDheto_2A>10]  <- 10
HalfLife_3_4_EEDheto_2A[HalfLife_3_4_EEDheto_2A<0]  <- 0
HalfLife_3_4_EEDheto_2A[HalfLife_3_4_EEDheto_2A>50]  <- 50

sink( file=paste(subdir_3_part4,  "/4-3-2H-runLog.txt",      sep = "") )
NTR_3_4_EEDko_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0_EEDko[4, i],  reduceRow2_Average_week4_EEDko[4, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_3_part4,  "/4-3-2H-LogLinearModel", sep="") )
  NTR_3_4_EEDko_2B <- c(NTR_3_4_EEDko_2B, NTR_bin)
}
sink()  
length(NTR_3_4_EEDko_2B)
summary(NTR_3_4_EEDko_2B)
HalfLife_3_4_EEDko_2B <- log(2)/(NTR_3_4_EEDko_2B+0.0001)
length(HalfLife_3_4_EEDko_2B)
summary(HalfLife_3_4_EEDko_2B)
NTR_3_4_EEDko_2B[NTR_3_4_EEDko_2B<0]  <- 0
NTR_3_4_EEDko_2B[NTR_3_4_EEDko_2B<10]  <- 10
HalfLife_3_4_EEDko_2B[HalfLife_3_4_EEDko_2B<0]  <- 0
HalfLife_3_4_EEDko_2B[HalfLife_3_4_EEDko_2B<50]  <- 50


MyAverageLines_1(vector2=c(NTR_3_1_EEDheto_2A,  NTR_3_2_EEDheto_2A,  NTR_3_3_EEDheto_2A,  NTR_3_4_EEDheto_2A),    
                 numSample2=4,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",   "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",   "High"="red4",  "Highest"="olivedrab4" ), 
                 path2=subdir_3_part4,     fileName2="4-3-2A-EEDheto-NTR-average-curve",  
                 title2="All Genes (EEDheto)",     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.3,    height2=3.3,   width2=5.5  , center2=myCenter_g )

MyAverageLines_1(vector2=c(NTR_3_1_EEDko_2B,  NTR_3_2_EEDko_2B,  NTR_3_3_EEDko_2B,  NTR_3_4_EEDko_2B ),    
                 numSample2=4,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",  "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",   "High"="red4",  "Highest"="olivedrab4" ), 
                 path2=subdir_3_part4,     fileName2="4-3-2B-EEDko-NTR-average-curve",  
                 title2="All Genes (EEDko)",     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.3,    height2=3.3,   width2=5.5  , center2=myCenter_g )




MyAverageLines_1(vector2=c(HalfLife_3_1_EEDheto_2A,  HalfLife_3_2_EEDheto_2A,  HalfLife_3_3_EEDheto_2A,  HalfLife_3_4_EEDheto_2A ),    
                 numSample2=4,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",   "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",    "High"="red4",  "Highest"="olivedrab4" ), 
                 path2=subdir_3_part4,     fileName2="4-3-2C-EEDheto-HalfLife-average-curve",  
                 title2="All Genes (EEDheto)",     xLab2="Relative distance (kb)",    yLab2="Half Life",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.5  , center2=myCenter_g )

MyAverageLines_1(vector2=c(HalfLife_3_1_EEDko_2B,  HalfLife_3_2_EEDko_2B,  HalfLife_3_3_EEDko_2B,  HalfLife_3_4_EEDko_2B ),    
                 numSample2=4,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",  "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",   "High"="red4",  "Highest"="olivedrab4" ), 
                 path2=subdir_3_part4,     fileName2="4-3-2D-EEDko-HalfLife-average-curve",  
                 title2="All Genes (EEDko)",     xLab2="Relative distance (kb)",    yLab2="Half Life",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.5  , center2=myCenter_g )









# TAC 2 samples
sink( file=paste(subdir_3_part4,  "/4-3-3A-runLog.txt",      sep = "") )
NTR_3_1_banding_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0[1, i],  reduceRow2_Average_banding[1, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_3_part4,  "/4-3-3A-LogLinearModel", sep="") )
  NTR_3_1_banding_2A <- c(NTR_3_1_banding_2A, NTR_bin)
}
sink()  
length(NTR_3_1_banding_2A)
summary(NTR_3_1_banding_2A)
HalfLife_3_1_banding_2A <- log(2)/(NTR_3_1_banding_2A+0.0001)
length(HalfLife_3_1_banding_2A)
summary(HalfLife_3_1_banding_2A)
NTR_3_1_banding_2A[NTR_3_1_banding_2A<0]  <- 0
NTR_3_1_banding_2A[NTR_3_1_banding_2A>10]  <- 10
HalfLife_3_1_banding_2A[HalfLife_3_1_banding_2A<0]  <- 0
HalfLife_3_1_banding_2A[HalfLife_3_1_banding_2A>50]  <- 50

sink( file=paste(subdir_3_part4,  "/4-3-3B-runLog.txt",      sep = "") )
NTR_3_1_sham_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0[1, i],  reduceRow2_Average_sham[1, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_3_part4,  "/4-3-3B-LogLinearModel", sep=""))
  NTR_3_1_sham_2B <- c(NTR_3_1_sham_2B, NTR_bin)
}
sink()  
length(NTR_3_1_sham_2B)
summary(NTR_3_1_sham_2B)
HalfLife_3_1_sham_2B <- log(2)/(NTR_3_1_sham_2B+0.0001)
length(HalfLife_3_1_sham_2B)
summary(HalfLife_3_1_sham_2B)
NTR_3_1_sham_2B[NTR_3_1_sham_2B<0]  <- 0
NTR_3_1_sham_2B[NTR_3_1_sham_2B>10]  <- 10
HalfLife_3_1_sham_2B[HalfLife_3_1_sham_2B<0]  <- 0
HalfLife_3_1_sham_2B[HalfLife_3_1_sham_2B>50]  <- 50




sink( file=paste(subdir_3_part4,  "/4-3-3C-runLog.txt",      sep = "") )
NTR_3_2_banding_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0[2, i],  reduceRow2_Average_banding[2, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_3_part4,  "/4-3-3C-LogLinearModel", sep=""))
  NTR_3_2_banding_2A <- c(NTR_3_2_banding_2A, NTR_bin)
}
sink()  
length(NTR_3_2_banding_2A)
summary(NTR_3_2_banding_2A)
HalfLife_3_2_banding_2A <- log(2)/(NTR_3_2_banding_2A+0.0001)
length(HalfLife_3_2_banding_2A)
summary(HalfLife_3_2_banding_2A)
NTR_3_2_banding_2A[NTR_3_2_banding_2A<0]  <- 0
NTR_3_2_banding_2A[NTR_3_2_banding_2A>10]  <- 10
HalfLife_3_2_banding_2A[HalfLife_3_2_banding_2A<0]  <- 0
HalfLife_3_2_banding_2A[HalfLife_3_2_banding_2A>50]  <- 50

sink( file=paste(subdir_3_part4,  "/4-3-3D-runLog.txt",      sep = "") )
NTR_3_2_sham_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0[2, i],  reduceRow2_Average_sham[2, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_3_part4,  "/4-3-3D-LogLinearModel", sep="") )
  NTR_3_2_sham_2B <- c(NTR_3_2_sham_2B, NTR_bin)
}
sink()  
length(NTR_3_2_sham_2B)
summary(NTR_3_2_sham_2B)
HalfLife_3_2_sham_2B <- log(2)/(NTR_3_2_sham_2B+0.0001)
length(HalfLife_3_2_sham_2B)
summary(HalfLife_3_2_sham_2B)
NTR_3_2_sham_2B[NTR_3_2_sham_2B<0]  <- 0
NTR_3_2_sham_2B[NTR_3_2_sham_2B>10]  <- 10
HalfLife_3_2_sham_2B[HalfLife_3_2_sham_2B<0]  <- 0
HalfLife_3_2_sham_2B[HalfLife_3_2_sham_2B>50]  <- 50




sink( file=paste(subdir_3_part4,  "/4-3-3E-runLog.txt",      sep = "") )
NTR_3_3_banding_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0[3, i],  reduceRow2_Average_banding[3, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_3_part4,  "/4-3-3E-LogLinearModel", sep="") )
  NTR_3_3_banding_2A <- c(NTR_3_3_banding_2A, NTR_bin)
}
sink()  
length(NTR_3_3_banding_2A)
summary(NTR_3_3_banding_2A)
HalfLife_3_3_banding_2A <- log(2)/(NTR_3_3_banding_2A+0.0001)
length(HalfLife_3_3_banding_2A)
summary(HalfLife_3_3_banding_2A)
NTR_3_3_banding_2A[NTR_3_3_banding_2A<0]  <- 0
NTR_3_3_banding_2A[NTR_3_3_banding_2A>10]  <- 10
HalfLife_3_3_banding_2A[HalfLife_3_3_banding_2A<0]  <- 0
HalfLife_3_3_banding_2A[HalfLife_3_3_banding_2A>50]  <- 50

sink( file=paste(subdir_3_part4,  "/4-3-3F-runLog.txt",      sep = "") )
NTR_3_3_sham_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0[3, i],  reduceRow2_Average_sham[3, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_3_part4,  "/4-3-3F-LogLinearModel", sep="") )
  NTR_3_3_sham_2B <- c(NTR_3_3_sham_2B, NTR_bin)
}
sink()  
length(NTR_3_3_sham_2B)
summary(NTR_3_3_sham_2B)
HalfLife_3_3_sham_2B <- log(2)/(NTR_3_3_sham_2B+0.0001)
length(HalfLife_3_3_sham_2B)
summary(HalfLife_3_3_sham_2B)
NTR_3_3_sham_2B[NTR_3_3_sham_2B<0]  <- 0
NTR_3_3_sham_2B[NTR_3_3_sham_2B>10]  <- 10
HalfLife_3_3_sham_2B[HalfLife_3_3_sham_2B<0]  <- 0
HalfLife_3_3_sham_2B[HalfLife_3_3_sham_2B>50]  <- 50




sink( file=paste(subdir_3_part4,  "/4-3-3G-runLog.txt",      sep = "") )
NTR_3_4_banding_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0[4, i],  reduceRow2_Average_banding[4, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_3_part4,  "/4-3-3G-LogLinearModel", sep="") )
  NTR_3_4_banding_2A <- c(NTR_3_4_banding_2A, NTR_bin)
}
sink()  
length(NTR_3_4_banding_2A)
summary(NTR_3_4_banding_2A)
HalfLife_3_4_banding_2A <- log(2)/(NTR_3_4_banding_2A+0.0001)
length(HalfLife_3_4_banding_2A)
summary(HalfLife_3_4_banding_2A)
NTR_3_4_banding_2A[NTR_3_4_banding_2A<0]  <- 0
NTR_3_4_banding_2A[NTR_3_4_banding_2A>10]  <- 10
HalfLife_3_4_banding_2A[HalfLife_3_4_banding_2A<0]  <- 0
HalfLife_3_4_banding_2A[HalfLife_3_4_banding_2A>50]  <- 50

sink( file=paste(subdir_3_part4,  "/4-3-3H-runLog.txt",      sep = "") )
NTR_3_4_sham_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow2_Average_week0[4, i],  reduceRow2_Average_sham[4, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_3_part4,  "/4-3-3H-LogLinearModel", sep=""))
  NTR_3_4_sham_2B <- c(NTR_3_4_sham_2B, NTR_bin)
}
sink()  
length(NTR_3_4_sham_2B)
summary(NTR_3_4_sham_2B)
HalfLife_3_4_sham_2B <- log(2)/(NTR_3_4_sham_2B+0.0001)
length(HalfLife_3_4_sham_2B)
summary(HalfLife_3_4_sham_2B)
NTR_3_4_sham_2B[NTR_3_4_sham_2B<0]  <- 0
NTR_3_4_sham_2B[NTR_3_4_sham_2B>10]  <- 10
HalfLife_3_4_sham_2B[HalfLife_3_4_sham_2B<0]  <- 0
HalfLife_3_4_sham_2B[HalfLife_3_4_sham_2B>50]  <- 50




MyAverageLines_1(vector2=c(NTR_3_1_banding_2A,  NTR_3_2_banding_2A,  NTR_3_3_banding_2A,  NTR_3_4_banding_2A ),    
                 numSample2=4,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",    "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",     "High"="red4",  "Highest"="olivedrab4" ), 
                 path2=subdir_3_part4,     fileName2="4-3-3A-banding-NTR-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.8,    height2=3.3,   width2=5.5  , center2=myCenter_g )

MyAverageLines_1(vector2=c(NTR_3_1_sham_2B,  NTR_3_2_sham_2B,  NTR_3_3_sham_2B,  NTR_3_4_sham_2B ),    
                 numSample2=4,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                 rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",     "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",   "High"="red4",  "Highest"="olivedrab4" ), 
                 path2=subdir_3_part4,     fileName2="4-3-3B-sham-NTR-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.8,    height2=3.3,   width2=5.5  , center2=myCenter_g )


MyAverageLines_1(vector2=c(HalfLife_3_1_banding_2A,  HalfLife_3_2_banding_2A,  HalfLife_3_3_banding_2A,  HalfLife_3_4_banding_2A ),    
                 numSample2=4,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",   "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",    "High"="olivedrab1",  "Highest"="olivedrab4" ), 
                 path2=subdir_3_part4,     fileName2="4-3-3C-banding-HalfLife-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half Life",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.5 , center2=myCenter_g  )

MyAverageLines_1(vector2=c(HalfLife_3_1_sham_2B,  HalfLife_3_2_sham_2B,  HalfLife_3_3_sham_2B,  HalfLife_3_4_sham_2B ),    
                 numSample2=4,   
                 sampleType2=c( rep("Lowest", numOfColumns1),  rep("Low", numOfColumns1)  ,  
                                rep("High", numOfColumns1)  ,   rep("Highest", numOfColumns1)    ), 
                 sampleRank2=c( "Lowest",   "Low",     "High",  "Highest"),     
                 colours2=c( "Lowest"="pink2",   "Low"="red",    "High"="red4",  "Highest"="olivedrab4" ), 
                 path2=subdir_3_part4,     fileName2="4-3-3D-sham-HalfLife-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half Life",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.5  , center2=myCenter_g )
























###############  3 categaries based on rows
subdir_4_part4 <- paste(Part4_g,  "/4-NTR-rows3Classes", sep = "")
if( ! file.exists(subdir_4_part4) ) { dir.create(subdir_4_part4) }


# WT 6 samples
sink( file=paste(subdir_4_part4,  "/4-4-1A-runLog.txt",      sep = "") )
NTR_4_1_WT_1A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow3_Average_week0[1, i], reduceRow3_Average_week1[1, i], reduceRow3_Average_week2[1, i], reduceRow3_Average_week4[1, i], reduceRow3_Average_week6[1, i], reduceRow3_Average_week8[1, i])
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log( vec1 ) , file2=paste(subdir_4_part4,  "/4-4-1A-LogLinearModel", sep="") )
  NTR_4_1_WT_1A <- c(NTR_4_1_WT_1A, NTR_bin)
}
sink()  
length(NTR_4_1_WT_1A)
summary(NTR_4_1_WT_1A)
HalfLife_4_1_WT_1A <- log(2)/(NTR_4_1_WT_1A+0.0001)
length(HalfLife_4_1_WT_1A)
summary(HalfLife_4_1_WT_1A)
NTR_4_1_WT_1A[NTR_4_1_WT_1A<0]  <- 0
NTR_4_1_WT_1A[NTR_4_1_WT_1A>10]  <- 10
HalfLife_4_1_WT_1A[HalfLife_4_1_WT_1A<0]  <- 0
HalfLife_4_1_WT_1A[HalfLife_4_1_WT_1A>50]  <- 50


sink( file=paste(subdir_4_part4,  "/4-4-1B-runLog.txt",      sep = "") )
NTR_4_2_WT_1A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow3_Average_week0[2, i], reduceRow3_Average_week1[2, i], reduceRow3_Average_week2[2, i], reduceRow3_Average_week4[2, i], reduceRow3_Average_week6[2, i], reduceRow3_Average_week8[2, i])
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log( vec1 )  , file2=paste(subdir_4_part4,  "/4-4-1B-LogLinearModel", sep="") )
  NTR_4_2_WT_1A <- c(NTR_4_2_WT_1A, NTR_bin)
}
sink()  
length(NTR_4_2_WT_1A)
summary(NTR_4_2_WT_1A)
HalfLife_4_2_WT_1A <- log(2)/(NTR_4_2_WT_1A+0.0001)
length(HalfLife_4_2_WT_1A)
summary(HalfLife_4_2_WT_1A)
NTR_4_2_WT_1A[NTR_4_2_WT_1A<0]  <- 0
NTR_4_2_WT_1A[NTR_4_2_WT_1A>10]  <- 10
HalfLife_4_2_WT_1A[HalfLife_4_2_WT_1A<0]  <- 0
HalfLife_4_2_WT_1A[HalfLife_4_2_WT_1A>50]  <- 50


sink( file=paste(subdir_4_part4,  "/4-4-1C-runLog.txt",      sep = "") )
NTR_4_3_WT_1A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow3_Average_week0[3, i], reduceRow3_Average_week1[3, i], reduceRow3_Average_week2[3, i], reduceRow3_Average_week4[3, i], reduceRow3_Average_week6[3, i], reduceRow3_Average_week8[3, i])
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log( vec1 )  , file2=paste(subdir_4_part4,  "/4-4-1C-LogLinearModel", sep="") )
  NTR_4_3_WT_1A <- c(NTR_4_3_WT_1A, NTR_bin)
}
sink()  
length(NTR_4_3_WT_1A)
summary(NTR_4_3_WT_1A)
HalfLife_4_3_WT_1A <- log(2)/(NTR_4_3_WT_1A+0.0001)
length(HalfLife_4_3_WT_1A)
summary(HalfLife_4_3_WT_1A)
NTR_4_3_WT_1A[NTR_4_3_WT_1A<0]  <- 0
NTR_4_3_WT_1A[NTR_4_3_WT_1A>10]  <- 10
HalfLife_4_3_WT_1A[HalfLife_4_3_WT_1A<0]  <- 0
HalfLife_4_3_WT_1A[HalfLife_4_3_WT_1A>50]  <- 50


MyAverageLines_1(vector2=c(NTR_4_1_WT_1A,  NTR_4_2_WT_1A,  NTR_4_3_WT_1A ),    
                 numSample2=3,   
                 sampleType2=c( rep("Low", numOfColumns1),  rep("Medium", numOfColumns1)  ,   rep("High", numOfColumns1)      ), 
                 sampleRank2=c( "Low",    "Medium",  "High"),     
                 colours2=c( "Low"="pink2",   "Medium"="red",   "High"="red4" ), 
                 path2=subdir_4_part4,     fileName2="4-4-1A-WT-NTR-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.3,    height2=3.3,   width2=5.5  , center2=myCenter_g )


MyAverageLines_1(vector2=c(HalfLife_4_1_WT_1A,  HalfLife_4_2_WT_1A,  HalfLife_4_3_WT_1A ),    
                 numSample2=3,   
                 sampleType2=c( rep("Low", numOfColumns1),  rep("Medium", numOfColumns1)  ,   rep("High", numOfColumns1)      ), 
                 sampleRank2=c( "Low",   "Medium",   "High" ),     
                 colours2=c( "Low"="pink2",   "Medium"="red",   "High"="red4"  ), 
                 path2=subdir_4_part4,     fileName2="4-4-1B-WT-HalfLife-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half life",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.5  , center2=myCenter_g )





# CKO 4 samples
sink( file=paste(subdir_4_part4,  "/4-4-2A-runLog.txt",      sep = "") )
NTR_4_1_EEDheto_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow3_Average_week0_EEDheto[1, i],  reduceRow3_Average_week4_EEDheto[1, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_4_part4,  "/4-4-2A-LogLinearModel", sep="")  )
  NTR_4_1_EEDheto_2A <- c(NTR_4_1_EEDheto_2A, NTR_bin)
}
sink()  
length(NTR_4_1_EEDheto_2A)
summary(NTR_4_1_EEDheto_2A)
HalfLife_4_1_EEDheto_2A <- log(2)/(NTR_4_1_EEDheto_2A+0.0001)
length(HalfLife_4_1_EEDheto_2A)
summary(HalfLife_4_1_EEDheto_2A)
NTR_4_1_EEDheto_2A[NTR_4_1_EEDheto_2A<0]  <- 0
NTR_4_1_EEDheto_2A[NTR_4_1_EEDheto_2A>10]  <- 10
HalfLife_4_1_EEDheto_2A[HalfLife_4_1_EEDheto_2A<0]  <- 0
HalfLife_4_1_EEDheto_2A[HalfLife_4_1_EEDheto_2A>50]  <- 50

sink( file=paste(subdir_4_part4,  "/4-4-2B-runLog.txt",      sep = "") )
NTR_4_1_EEDko_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow3_Average_week0_EEDko[1, i],  reduceRow3_Average_week4_EEDko[1, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_4_part4,  "/4-4-2B-LogLinearModel", sep="") )
  NTR_4_1_EEDko_2B <- c(NTR_4_1_EEDko_2B, NTR_bin)
}
sink()  
length(NTR_4_1_EEDko_2B)
summary(NTR_4_1_EEDko_2B)
HalfLife_4_1_EEDko_2B <- log(2)/(NTR_4_1_EEDko_2B+0.0001)
length(HalfLife_4_1_EEDko_2B)
summary(HalfLife_4_1_EEDko_2B)
NTR_4_1_EEDko_2B[NTR_4_1_EEDko_2B<0]   <- 0
NTR_4_1_EEDko_2B[NTR_4_1_EEDko_2B>10]  <- 10
HalfLife_4_1_EEDko_2B[HalfLife_4_1_EEDko_2B<0]  <- 0
HalfLife_4_1_EEDko_2B[HalfLife_4_1_EEDko_2B>50]  <- 50




sink( file=paste(subdir_4_part4,  "/4-4-2C-runLog.txt",      sep = "") )
NTR_4_2_EEDheto_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow3_Average_week0_EEDheto[2, i],  reduceRow3_Average_week4_EEDheto[2, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_4_part4,  "/4-4-2C-LogLinearModel", sep="") )
  NTR_4_2_EEDheto_2A <- c(NTR_4_2_EEDheto_2A, NTR_bin)
}
sink()  
length(NTR_4_2_EEDheto_2A)
summary(NTR_4_2_EEDheto_2A)
HalfLife_4_2_EEDheto_2A <- log(2)/(NTR_4_2_EEDheto_2A+0.0001)
length(HalfLife_4_2_EEDheto_2A)
summary(HalfLife_4_2_EEDheto_2A)
NTR_4_2_EEDheto_2A[NTR_4_2_EEDheto_2A<0]  <- 0
NTR_4_2_EEDheto_2A[NTR_4_2_EEDheto_2A>10]  <- 10
HalfLife_4_2_EEDheto_2A[HalfLife_4_2_EEDheto_2A<0]  <- 0
HalfLife_4_2_EEDheto_2A[HalfLife_4_2_EEDheto_2A>50]  <- 50

sink( file=paste(subdir_4_part4,  "/4-4-2D-runLog.txt",      sep = "") )
NTR_4_2_EEDko_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow3_Average_week0_EEDko[2, i],  reduceRow3_Average_week4_EEDko[2, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_4_part4,  "/4-4-2D-LogLinearModel", sep="") )
  NTR_4_2_EEDko_2B <- c(NTR_4_2_EEDko_2B, NTR_bin)
}
sink()  
length(NTR_4_2_EEDko_2B)
summary(NTR_4_2_EEDko_2B)
HalfLife_4_2_EEDko_2B <- log(2)/(NTR_4_2_EEDko_2B+0.0001)
length(HalfLife_4_2_EEDko_2B)
summary(HalfLife_4_2_EEDko_2B)
NTR_4_2_EEDko_2B[NTR_4_2_EEDko_2B<0]  <- 0
NTR_4_2_EEDko_2B[NTR_4_2_EEDko_2B>10]  <- 10
HalfLife_4_2_EEDko_2B[HalfLife_4_2_EEDko_2B<0]  <- 0
HalfLife_4_2_EEDko_2B[HalfLife_4_2_EEDko_2B>50]  <- 50




sink( file=paste(subdir_4_part4,  "/4-4-2E-runLog.txt",      sep = "") )
NTR_4_3_EEDheto_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow3_Average_week0_EEDheto[3, i],  reduceRow3_Average_week4_EEDheto[3, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_4_part4,  "/4-4-2E-LogLinearModel", sep="") )
  NTR_4_3_EEDheto_2A <- c(NTR_4_3_EEDheto_2A, NTR_bin)
}
sink()  
length(NTR_4_3_EEDheto_2A)
summary(NTR_4_3_EEDheto_2A)
HalfLife_4_3_EEDheto_2A <- log(2)/(NTR_4_3_EEDheto_2A+0.0001)
length(HalfLife_4_3_EEDheto_2A)
summary(HalfLife_4_3_EEDheto_2A)
NTR_4_3_EEDheto_2A[NTR_4_3_EEDheto_2A<0]  <- 0
NTR_4_3_EEDheto_2A[NTR_4_3_EEDheto_2A>10]  <- 10
HalfLife_4_3_EEDheto_2A[HalfLife_4_3_EEDheto_2A<0]  <- 0
HalfLife_4_3_EEDheto_2A[HalfLife_4_3_EEDheto_2A>50]  <- 50

sink( file=paste(subdir_4_part4,  "/4-4-2F-runLog.txt",      sep = "") )
NTR_4_3_EEDko_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow3_Average_week0_EEDko[3, i],  reduceRow3_Average_week4_EEDko[3, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_4_part4,  "/4-4-2F-LogLinearModel", sep="") )
  NTR_4_3_EEDko_2B <- c(NTR_4_3_EEDko_2B, NTR_bin)
}
sink()  
length(NTR_4_3_EEDko_2B)
summary(NTR_4_3_EEDko_2B)
HalfLife_4_3_EEDko_2B <- log(2)/(NTR_4_3_EEDko_2B+0.0001)
length(HalfLife_4_3_EEDko_2B)
summary(HalfLife_4_3_EEDko_2B)
NTR_4_3_EEDko_2B[NTR_4_3_EEDko_2B<0]  <- 0
NTR_4_3_EEDko_2B[NTR_4_3_EEDko_2B>10]  <- 10
HalfLife_4_3_EEDko_2B[HalfLife_4_3_EEDko_2B<0]  <- 0
HalfLife_4_3_EEDko_2B[HalfLife_4_3_EEDko_2B>50]  <- 50




MyAverageLines_1(vector2=c(NTR_4_1_EEDheto_2A,  NTR_4_2_EEDheto_2A,  NTR_4_3_EEDheto_2A ),    
                 numSample2=3,   
                 sampleType2=c( rep("Low", numOfColumns1),  rep("Medium", numOfColumns1)  ,   rep("High", numOfColumns1)      ), 
                 sampleRank2=c( "Low",    "Medium",  "High"),     
                 colours2=c( "Low"="pink2",   "Medium"="red",   "High"="red4" ), 
                 path2=subdir_4_part4,     fileName2="4-4-2A-EEDheto-NTR-average-curve",  
                 title2="All Genes (EEDheto)",     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.3,    height2=3.3,   width2=5.5  , center2=myCenter_g )

MyAverageLines_1(vector2=c(NTR_4_1_EEDko_2B,  NTR_4_2_EEDko_2B,  NTR_4_3_EEDko_2B  ),    
                 numSample2=3,   
                 sampleType2=c( rep("Low", numOfColumns1),  rep("Medium", numOfColumns1)  ,   rep("High", numOfColumns1)      ), 
                 sampleRank2=c( "Low",    "Medium",  "High"),     
                 colours2=c( "Low"="pink2",   "Medium"="red",   "High"="red4" ), 
                 path2=subdir_4_part4,     fileName2="4-4-2B-EEDko-NTR-average-curve",  
                 title2="All Genes (EEDko)",     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.3,    height2=3.3,   width2=5.5  , center2=myCenter_g )




MyAverageLines_1(vector2=c(HalfLife_4_1_EEDheto_2A,  HalfLife_4_2_EEDheto_2A,  HalfLife_4_3_EEDheto_2A  ),    
                 numSample2=3,   
                 sampleType2=c( rep("Low", numOfColumns1),  rep("Medium", numOfColumns1)  ,   rep("High", numOfColumns1)      ), 
                 sampleRank2=c( "Low",    "Medium",  "High"),     
                 colours2=c( "Low"="pink2",   "Medium"="red",   "High"="red4" ), 
                 path2=subdir_4_part4,     fileName2="4-4-2C-EEDheto-HalfLife-average-curve",  
                 title2="All Genes (EEDheto)",     xLab2="Relative distance (kb)",    yLab2="Half Life",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.5  , center2=myCenter_g )

MyAverageLines_1(vector2=c(HalfLife_4_1_EEDko_2B,  HalfLife_4_2_EEDko_2B,  HalfLife_4_3_EEDko_2B  ),    
                 numSample2=3,   
                 sampleType2=c( rep("Low", numOfColumns1),  rep("Medium", numOfColumns1)  ,   rep("High", numOfColumns1)      ), 
                 sampleRank2=c( "Low",    "Medium",  "High"),     
                 colours2=c( "Low"="pink2",   "Medium"="red",   "High"="red4" ), 
                 path2=subdir_4_part4,     fileName2="4-4-2D-EEDko-HalfLife-average-curve",  
                 title2="All Genes (EEDko)",     xLab2="Relative distance (kb)",    yLab2="Half Life",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.5  , center2=myCenter_g )









# TAC 2 samples
sink( file=paste(subdir_4_part4,  "/4-4-3A-runLog.txt",      sep = "") )
NTR_4_1_banding_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow3_Average_week0[1, i],  reduceRow3_Average_banding[1, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_4_part4,  "/4-4-3A-LogLinearModel", sep="") )
  NTR_4_1_banding_2A <- c(NTR_4_1_banding_2A, NTR_bin)
}
sink()  
length(NTR_4_1_banding_2A)
summary(NTR_4_1_banding_2A)
HalfLife_4_1_banding_2A <- log(2)/(NTR_4_1_banding_2A+0.0001)
length(HalfLife_4_1_banding_2A)
summary(HalfLife_4_1_banding_2A)
NTR_4_1_banding_2A[NTR_4_1_banding_2A<0]  <- 0
NTR_4_1_banding_2A[NTR_4_1_banding_2A>10]  <- 10
HalfLife_4_1_banding_2A[HalfLife_4_1_banding_2A<0]  <- 0
HalfLife_4_1_banding_2A[HalfLife_4_1_banding_2A>50]  <- 50

sink( file=paste(subdir_4_part4,  "/4-4-3B-runLog.txt",      sep = "") )
NTR_4_1_sham_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow3_Average_week0[1, i],  reduceRow3_Average_sham[1, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 ) , file2=paste(subdir_4_part4,  "/4-4-3B-LogLinearModel", sep="")  )
  NTR_4_1_sham_2B <- c(NTR_4_1_sham_2B, NTR_bin)
}
sink()  
length(NTR_4_1_sham_2B)
summary(NTR_4_1_sham_2B)
HalfLife_4_1_sham_2B <- log(2)/(NTR_4_1_sham_2B+0.0001)
length(HalfLife_4_1_sham_2B)
summary(HalfLife_4_1_sham_2B)
NTR_4_1_sham_2B[NTR_4_1_sham_2B<0]  <- 0
NTR_4_1_sham_2B[NTR_4_1_sham_2B>10]  <- 10
HalfLife_4_1_sham_2B[HalfLife_4_1_sham_2B<0]  <- 0
HalfLife_4_1_sham_2B[HalfLife_4_1_sham_2B>50]  <- 50




sink( file=paste(subdir_4_part4,  "/4-4-3C-runLog.txt",      sep = "") )
NTR_4_2_banding_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow3_Average_week0[2, i],  reduceRow3_Average_banding[2, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_4_part4,  "/4-4-3C-LogLinearModel", sep="") )
  NTR_4_2_banding_2A <- c(NTR_4_2_banding_2A, NTR_bin)
}
sink()  
length(NTR_4_2_banding_2A)
summary(NTR_4_2_banding_2A)
HalfLife_4_2_banding_2A <- log(2)/(NTR_4_2_banding_2A+0.0001)
length(HalfLife_4_2_banding_2A)
summary(HalfLife_4_2_banding_2A)
NTR_4_2_banding_2A[NTR_4_2_banding_2A<0]  <- 0
NTR_4_2_banding_2A[NTR_4_2_banding_2A>10]  <- 10
HalfLife_4_2_banding_2A[HalfLife_4_2_banding_2A<0]  <- 0
HalfLife_4_2_banding_2A[HalfLife_4_2_banding_2A>50]  <- 50

sink( file=paste(subdir_4_part4,  "/4-4-3D-runLog.txt",      sep = "") )
NTR_4_2_sham_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow3_Average_week0[2, i],  reduceRow3_Average_sham[2, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_4_part4,  "/4-4-3D-LogLinearModel", sep="") )
  NTR_4_2_sham_2B <- c(NTR_4_2_sham_2B, NTR_bin)
}
sink()  
length(NTR_4_2_sham_2B)
summary(NTR_4_2_sham_2B)
HalfLife_4_2_sham_2B <- log(2)/(NTR_4_2_sham_2B+0.0001)
length(HalfLife_4_2_sham_2B)
summary(HalfLife_4_2_sham_2B)
NTR_4_2_sham_2B[NTR_4_2_sham_2B<0]  <- 0
NTR_4_2_sham_2B[NTR_4_2_sham_2B>10]  <- 10
HalfLife_4_2_sham_2B[HalfLife_4_2_sham_2B<0]  <- 0
HalfLife_4_2_sham_2B[HalfLife_4_2_sham_2B>50]  <- 50




sink( file=paste(subdir_4_part4,  "/4-4-3E-runLog.txt",      sep = "") )
NTR_4_3_banding_2A  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow3_Average_week0[3, i],  reduceRow3_Average_banding[3, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_4_part4,  "/4-4-3E-LogLinearModel", sep="") )
  NTR_4_3_banding_2A <- c(NTR_4_3_banding_2A, NTR_bin)
}
sink()  
length(NTR_4_3_banding_2A)
summary(NTR_4_3_banding_2A)
HalfLife_4_3_banding_2A <- log(2)/(NTR_4_3_banding_2A+0.0001)
length(HalfLife_4_3_banding_2A)
summary(HalfLife_4_3_banding_2A)
NTR_4_3_banding_2A[NTR_4_3_banding_2A<0]  <- 0
NTR_4_3_banding_2A[NTR_4_3_banding_2A>10]  <- 10
HalfLife_4_3_banding_2A[HalfLife_4_3_banding_2A<0]  <- 0
HalfLife_4_3_banding_2A[HalfLife_4_3_banding_2A>50]  <- 50

sink( file=paste(subdir_4_part4,  "/4-4-3F-runLog.txt",      sep = "") )
NTR_4_3_sham_2B  <-  c()
for (i in c(1:numOfColumns1)) {
  vec1 <- c(reduceRow3_Average_week0[3, i],  reduceRow3_Average_sham[3, i] )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_4_part4,  "/4-4-3F-LogLinearModel", sep="") )
  NTR_4_3_sham_2B <- c(NTR_4_3_sham_2B, NTR_bin)
}
sink()  
length(NTR_4_3_sham_2B)
summary(NTR_4_3_sham_2B)
HalfLife_4_3_sham_2B <- log(2)/(NTR_4_3_sham_2B+0.0001)
length(HalfLife_4_3_sham_2B)
summary(HalfLife_4_3_sham_2B)
NTR_4_3_sham_2B[NTR_4_3_sham_2B<0]  <- 0
NTR_4_3_sham_2B[NTR_4_3_sham_2B>10]  <- 10
HalfLife_4_3_sham_2B[HalfLife_4_3_sham_2B<0]  <- 0
HalfLife_4_3_sham_2B[HalfLife_4_3_sham_2B>50]  <- 50



MyAverageLines_1(vector2=c(NTR_4_1_banding_2A,  NTR_4_2_banding_2A,  NTR_4_3_banding_2A  ),    
                 numSample2=3,   
                 sampleType2=c( rep("Low", numOfColumns1),  rep("Medium", numOfColumns1)  ,   rep("High", numOfColumns1)      ), 
                 sampleRank2=c( "Low",    "Medium",  "High"),     
                 colours2=c( "Low"="pink2",   "Medium"="red",   "High"="red4" ), 
                 path2=subdir_4_part4,     fileName2="4-4-3A-banding-NTR-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.8,    height2=3.3,   width2=5.5   , center2=myCenter_g)

MyAverageLines_1(vector2=c(NTR_4_1_sham_2B,  NTR_4_2_sham_2B,  NTR_4_3_sham_2B  ),    
                 numSample2=3,   
                 sampleType2=c( rep("Low", numOfColumns1),  rep("Medium", numOfColumns1)  ,   rep("High", numOfColumns1)      ), 
                 sampleRank2=c( "Low",    "Medium",  "High"),     
                 colours2=c( "Low"="pink2",   "Medium"="red",   "High"="red4" ), 
                 path2=subdir_4_part4,     fileName2="4-4-3B-sham-NTR-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.8,    height2=3.3,   width2=5.5   , center2=myCenter_g)


MyAverageLines_1(vector2=c(HalfLife_4_1_banding_2A,  HalfLife_4_2_banding_2A,  HalfLife_4_3_banding_2A  ),    
                 numSample2=3,   
                 sampleType2=c( rep("Low", numOfColumns1),  rep("Medium", numOfColumns1)  ,   rep("High", numOfColumns1)      ), 
                 sampleRank2=c( "Low",    "Medium",  "High"),     
                 colours2=c( "Low"="pink2",   "Medium"="red",   "High"="red4" ), 
                 path2=subdir_4_part4,     fileName2="4-4-3C-banding-HalfLife-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half Life",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.5  , center2=myCenter_g )

MyAverageLines_1(vector2=c(HalfLife_4_1_sham_2B,  HalfLife_4_2_sham_2B,  HalfLife_4_3_sham_2B  ),    
                 numSample2=3,   
                 sampleType2=c( rep("Low", numOfColumns1),  rep("Medium", numOfColumns1)  ,   rep("High", numOfColumns1)      ), 
                 sampleRank2=c( "Low",    "Medium",  "High"),     
                 colours2=c( "Low"="pink2",   "Medium"="red",   "High"="red4" ), 
                 path2=subdir_4_part4,     fileName2="4-4-3D-sham-HalfLife-average-curve",  
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half Life",   
                 Ymin2=0,   Ymax2=50,    height2=3.3,   width2=5.5  , center2=myCenter_g )





















#########################################################################################################  average some columns for each row
subdir_5_part4 <- paste(Part4_g,  "/5-NTR-eachRow-5kb", sep = "")
if( ! file.exists(subdir_5_part4) ) { dir.create(subdir_5_part4) }


############################# WT 6 samples
sink( file=paste(subdir_5_part4,  "/4-5-1A-WT-runLog.txt",      sep = "") )
myNTR_5_WT_1  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(row_Average_week0[i]+0.0006, row_Average_week1[i]+0.0005, row_Average_week2[i]+0.0004, row_Average_week4[i]+0.0003, row_Average_week6[i]+0.0002, row_Average_week8[i]+0.0001)                        
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log(vec1)  , file2=paste(subdir_5_part4,  "/4-5-1B-WT-LogLinearModel", sep="") )
  myNTR_5_WT_1 <- c(myNTR_5_WT_1, NTR_bin)
}
sink()  

length(myNTR_5_WT_1)
summary(myNTR_5_WT_1)
HalfLife_5_WT_1 <- log(2)/(myNTR_5_WT_1+0.0001)
length(HalfLife_5_WT_1)
summary(HalfLife_5_WT_1)

write.table(x=myNTR_5_WT_1,      file = paste(subdir_5_part4,  "/4-5-1C-raw-NTR.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_5_WT_1,   file = paste(subdir_5_part4,  "/4-5-1D-raw-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_5_WT_1 <- data.frame( xAxis = c(1:length(myNTR_5_WT_1)),      yAxis = sort(myNTR_5_WT_1) )
FigureTemp_5_WT_1 <- ggplot(myframe_5_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_WT_1,  path1=subdir_5_part4, fileName1="4-5-1A-figure-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_5_WT_1 <- ggplot(myframe_5_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) +  
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_WT_1,  path1=subdir_5_part4, fileName1="4-5-1A-figure-Line-raw",  height1=4, width1=7)

myNTR_5_WT_1[myNTR_5_WT_1<0]  <- 0
myNTR_5_WT_1[myNTR_5_WT_1>1]  <- 1
HalfLife_5_WT_1[HalfLife_5_WT_1<0]  <- 0
HalfLife_5_WT_1[HalfLife_5_WT_1>50] <- 50

write.table(x=thisRowNames1,     file = paste(subdir_5_part4,   "/4-5-1E-RowNames.txt",  sep = ""),         append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=myNTR_5_WT_1,      file = paste(subdir_5_part4,  "/4-5-1F-noLess0-NTR.txt",    sep = ""),     append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_5_WT_1,   file = paste(subdir_5_part4,  "/4-5-1G-noLess0-HalfLife.txt", sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_5_WT_1 <- data.frame( xAxis = c(1:length(myNTR_5_WT_1)),      yAxis = sort(myNTR_5_WT_1) )
FigureTemp_5_WT_1 <- ggplot(myframe_5_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_WT_1,  path1=subdir_5_part4, fileName1="4-5-1B-figure-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_5_WT_1 <- ggplot(myframe_5_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) +  
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_WT_1,  path1=subdir_5_part4, fileName1="4-5-1B-figure-Line-noLess0",  height1=4, width1=7)

MyBoxViolinPlot_1(vector2=myNTR_5_WT_1,      sampleType2=c( rep("NTR", numOfRows1)  ), 
                  sampleRank2=c( "NTR" ),     colours2=c( "NTR"="red" ), 
                  path2=subdir_5_part4,     fileName2="4-5-1C-WT-NTR-average-curve",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="NTR",   
                  Ymin2=0,   Ymax2=0.4,    height2=5,   width2=1.5  )

MyBoxViolinPlot_1(vector2=HalfLife_5_WT_1,      sampleType2=c( rep("HalfLife", numOfRows1)  ), 
                  sampleRank2=c( "HalfLife" ),     colours2=c( "HalfLife"="red" ), 
                  path2=subdir_5_part4,     fileName2="4-5-1D-WT-HalfLife-average-curve",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="Half life (weeks)",   
                  Ymin2=0,   Ymax2=20,    height2=5,   width2=1.5  )


dim(Average_H3)
dim(Average_week0)
numRows_oneClass <- floor( nrow(Average_H3)/5 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
INDEX_4 <- seq(from = numRows_oneClass*3+1, to = numRows_oneClass*4, by =1 )
INDEX_5 <- seq(from = numRows_oneClass*4+1, to = numRows_oneClass*5, by =1 )
nrow(Average_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) + length(INDEX_4) + length(INDEX_5))

MyBoxViolinPlot_1(vector2=c( myNTR_5_WT_1[INDEX_1], myNTR_5_WT_1[INDEX_2],  myNTR_5_WT_1[INDEX_3],  myNTR_5_WT_1[INDEX_4],  myNTR_5_WT_1[INDEX_5] ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),  rep("Medium", length(INDEX_3)),  
                                 rep("High", length(INDEX_4)),      rep("Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "Lowest",  "Low",   "Medium",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",    "Medium"="green2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_5_part4,   fileName2= paste("4-5-1E-WT-NTR-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=3.5,   Ymin2=0, Ymax2=0.4)    ## width = 1 + 5*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c( HalfLife_5_WT_1[INDEX_1], HalfLife_5_WT_1[INDEX_2],  HalfLife_5_WT_1[INDEX_3],  HalfLife_5_WT_1[INDEX_4],  HalfLife_5_WT_1[INDEX_5] ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),  rep("Medium", length(INDEX_3)),  
                                 rep("High", length(INDEX_4)),      rep("Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "Lowest",  "Low",   "Medium",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",    "Medium"="green2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_5_part4,   fileName2= paste("4-5-1F-WT-HalfLife-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=3.5,   Ymin2=0, Ymax2=40)    ## width = 1 + 5*0.5, height=5cm 












############################# CKO 4 samples
sink( file=paste(subdir_5_part4,  "/4-5-2A-1-EEDheto-runLog.txt",      sep = "") )
myNTR_5_EEDheto_2A  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(row_Average_week0_EEDheto[i]+0.0001,  row_Average_week4_EEDheto[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_5_part4,  "/4-5-2A-2-EEDheto-LogLinearModel", sep="") )
  myNTR_5_EEDheto_2A <- c(myNTR_5_EEDheto_2A, NTR_bin)
}
sink()  

length(myNTR_5_EEDheto_2A)
summary(myNTR_5_EEDheto_2A)
HalfLife_5_EEDheto_2A <- log(2)/(myNTR_5_EEDheto_2A+0.0001)
length(HalfLife_5_EEDheto_2A)
summary(HalfLife_5_EEDheto_2A)

write.table(x=myNTR_5_EEDheto_2A,      file = paste(subdir_5_part4,  "/4-5-2A-3-raw-NTR_EEDheto.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_5_EEDheto_2A,   file = paste(subdir_5_part4,  "/4-5-2A-4-raw-HalfLife_EEDheto.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_5_EEDheto_2A <- data.frame( xAxis = c(1:length(myNTR_5_EEDheto_2A)),      yAxis = sort(myNTR_5_EEDheto_2A) )
FigureTemp_5_EEDheto_2A <- ggplot(myframe_5_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_EEDheto_2A,  path1=subdir_5_part4, fileName1="4-5-2A-1-figure-EEDheto-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_5_EEDheto_2A <- ggplot(myframe_5_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) +
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_EEDheto_2A,  path1=subdir_5_part4, fileName1="4-5-2A-1-figure-EEDheto-Line-raw",  height1=4, width1=7)

myNTR_5_EEDheto_2A[myNTR_5_EEDheto_2A<0]  <- 0
myNTR_5_EEDheto_2A[myNTR_5_EEDheto_2A>1]  <- 1
HalfLife_5_EEDheto_2A[HalfLife_5_EEDheto_2A<0]   <- 0
HalfLife_5_EEDheto_2A[HalfLife_5_EEDheto_2A>50]  <- 50

write.table(x=myNTR_5_EEDheto_2A,      file = paste(subdir_5_part4,  "/4-5-2A-5-noLess0-myNTR_EEDheto.txt",  sep = ""),      append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_5_EEDheto_2A,   file = paste(subdir_5_part4,  "/4-5-2A-6-noLess0-HalfLife_EEDheto.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_5_EEDheto_2A <- data.frame( xAxis = c(1:length(myNTR_5_EEDheto_2A)),      yAxis = sort(myNTR_5_EEDheto_2A) )
FigureTemp_5_EEDheto_2A <- ggplot(myframe_5_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_EEDheto_2A,  path1=subdir_5_part4, fileName1="4-5-2A-2-figure-EEDheto-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_5_EEDheto_2A <- ggplot(myframe_5_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_EEDheto_2A,  path1=subdir_5_part4, fileName1="4-5-2A-2-figure-EEDheto-Line-noLess0",  height1=4, width1=7)





sink( file=paste(subdir_5_part4,  "/4-5-2B-1-EEDko-runLog.txt",      sep = "") )
myNTR_5_EEDko_2B  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(row_Average_week0_EEDko[i]+0.0001,  row_Average_week4_EEDko[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_5_part4,  "/4-5-2B-2-EEDko-LogLinearModel", sep="") )
  myNTR_5_EEDko_2B <- c(myNTR_5_EEDko_2B, NTR_bin)
}
sink()  

length(myNTR_5_EEDko_2B)
summary(myNTR_5_EEDko_2B)
HalfLife_5_EEDko_2B <- log(2)/(myNTR_5_EEDko_2B+0.0001)
length(HalfLife_5_EEDko_2B)
summary(HalfLife_5_EEDko_2B)

write.table(x=myNTR_5_EEDko_2B,      file = paste(subdir_5_part4,  "/4-5-2B-3-raw-NTR_EEDko.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_5_EEDko_2B,   file = paste(subdir_5_part4,  "/4-5-2B-4-raw-HalfLife_EEDko.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_5_EEDko_2B <- data.frame( xAxis = c(1:length(myNTR_5_EEDko_2B)),      yAxis = sort(myNTR_5_EEDko_2B) )
FigureTemp_5_EEDko_2B <- ggplot(myframe_5_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_EEDko_2B,  path1=subdir_5_part4, fileName1="4-5-2B-1-figure-EEDko-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_5_EEDko_2B <- ggplot(myframe_5_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_EEDko_2B,  path1=subdir_5_part4, fileName1="4-5-2B-1-figure-EEDko-Line-raw",  height1=4, width1=7)

myNTR_5_EEDko_2B[myNTR_5_EEDko_2B<0]  <- 0
myNTR_5_EEDko_2B[myNTR_5_EEDko_2B>1]  <- 1
HalfLife_5_EEDko_2B[HalfLife_5_EEDko_2B<0]   <- 0
HalfLife_5_EEDko_2B[HalfLife_5_EEDko_2B>50]  <- 50

write.table(x=myNTR_5_EEDko_2B,      file = paste(subdir_5_part4,  "/4-5-2B-5-noLess0-myNTR_EEDko.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_5_EEDko_2B,   file = paste(subdir_5_part4,  "/4-5-2B-6-noLess0-HalfLife_EEDko.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_5_EEDko_2B <- data.frame( xAxis = c(1:length(myNTR_5_EEDko_2B)),      yAxis = sort(myNTR_5_EEDko_2B) )
FigureTemp_5_EEDko_2B <- ggplot(myframe_5_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_EEDko_2B,  path1=subdir_5_part4, fileName1="4-5-2B-2-figure-EEDko-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_5_EEDko_2B <- ggplot(myframe_5_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_EEDko_2B,  path1=subdir_5_part4, fileName1="4-5-2B-2-figure-EEDko-Line-noLess0",  height1=4, width1=7)



MyBoxViolinPlot_1(vector2=c(myNTR_5_EEDheto_2A,  myNTR_5_EEDko_2B),       
                  sampleType2=c( rep("EEDheto", numOfRows1)  ,   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "EEDheto",  "EEDko" ),     
                  colours2=c( "EEDheto"="red",  "EEDko"="blue"  ), 
                  path2=subdir_5_part4,     fileName2="4-5-2A-merge-CKO-NTR-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                  Ymin2=0,   Ymax2=0.2,    height2=5,   width2=2  )

MyBoxViolinPlot_1(vector2=c(HalfLife_5_EEDheto_2A,  HalfLife_5_EEDko_2B),       
                  sampleType2=c( rep("EEDheto", numOfRows1)  ,   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "EEDheto",  "EEDko" ),     
                  colours2=c( "EEDheto"="red",  "EEDko"="blue"  ), 
                  path2=subdir_5_part4,     fileName2="4-5-2B-merge-CKO-HalfLife-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half life (weeks)",   
                  Ymin2=0,   Ymax2=20,    height2=5,   width2=2  )



diff1_5_CKO = myNTR_5_EEDko_2B - myNTR_5_EEDheto_2A
diff2_5_CKO = log2( ( (myNTR_5_EEDko_2B+0.001) / (myNTR_5_EEDheto_2A+0.001))  )
diff3_5_CKO = log2( ( (myNTR_5_EEDheto_2A+0.001) / (myNTR_5_EEDko_2B+0.001))  )   ## 2^(-10) = 0.001
length(diff1_5_CKO)
length(diff2_5_CKO)
length(diff3_5_CKO)
summary(diff1_5_CKO)
summary(diff2_5_CKO)
summary(diff3_5_CKO)


myframe_5_CKO <- data.frame( xAxis = c(1:length(diff1_5_CKO)),      yAxis = sort(diff1_5_CKO) )

FigureTemp_5_CKO <- ggplot(myframe_5_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR (EEDko - EEDheto)") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.1, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -0.1, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_CKO,  path1=subdir_5_part4, fileName1="4-5-2C-merge-figure-CKO-Line-minus-limitY",  height1=4, width1=7)

FigureTemp_5_CKO <- ggplot(myframe_5_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR (EEDko - EEDheto)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.1, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -0.1, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_CKO,  path1=subdir_5_part4, fileName1="4-5-2C-merge-figure-CKO-Line-minus",  height1=4, width1=7)



myframe_5_CKO <- data.frame( xAxis = c(1:length(diff2_5_CKO)),      yAxis = sort(diff2_5_CKO) )

FigureTemp_5_CKO <- ggplot(myframe_5_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDko/EEDheto)") +  ggtitle(myTitle_g) + ylim(-10, 10 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_CKO,  path1=subdir_5_part4, fileName1="4-5-2D-merge-figure-CKO-Line-ratio-limitY",  height1=4, width1=7)

FigureTemp_5_CKO <- ggplot(myframe_5_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDko/EEDheto)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_CKO,  path1=subdir_5_part4, fileName1="4-5-2D-merge-figure-CKO-Line-ratio",  height1=4, width1=7)





myframe_5_CKO <- data.frame( xAxis = c(1:length(diff3_5_CKO)),      yAxis = sort(diff3_5_CKO) )

FigureTemp_5_CKO <- ggplot(myframe_5_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDheto/EEDko)") +  ggtitle(myTitle_g) + ylim(-10, 10 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_CKO,  path1=subdir_5_part4, fileName1="4-5-2E-merge-figure-CKO-Line-ratio-limitY",  height1=4, width1=7)

FigureTemp_5_CKO <- ggplot(myframe_5_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDheto/EEDko)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_CKO,  path1=subdir_5_part4, fileName1="4-5-2E-merge-figure-CKO-Line-ratio",  height1=4, width1=7)


MyBoxViolinPlot_1(vector2=c( myNTR_5_EEDheto_2A[INDEX_1], myNTR_5_EEDheto_2A[INDEX_2],  myNTR_5_EEDheto_2A[INDEX_3],  myNTR_5_EEDheto_2A[INDEX_4],  myNTR_5_EEDheto_2A[INDEX_5],
                             myNTR_5_EEDko_2B[INDEX_1],   myNTR_5_EEDko_2B[INDEX_2],    myNTR_5_EEDko_2B[INDEX_3],    myNTR_5_EEDko_2B[INDEX_4],    myNTR_5_EEDko_2B[INDEX_5] 
                             ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),  rep("EEDheto_Medium", length(INDEX_3)),  rep("EEDheto_High", length(INDEX_4)),  rep("EEDheto_Highest", length(INDEX_5)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),    rep("EEDko_Medium", length(INDEX_3)),    rep("EEDko_High", length(INDEX_4)),    rep("EEDko_Highest", length(INDEX_5))
                                 ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",     "EEDko_Medium",    "EEDko_High",    "EEDko_Highest"  
                                 ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2",    "EEDko_Highest"="purple2" 
                              ), 
                  path2=subdir_5_part4,   fileName2= paste("4-5-2F-WT-NTR-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=0.3)    ## width = 1 + 10*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_5_EEDheto_2A[INDEX_1], HalfLife_5_EEDheto_2A[INDEX_2],  HalfLife_5_EEDheto_2A[INDEX_3],  HalfLife_5_EEDheto_2A[INDEX_4],  HalfLife_5_EEDheto_2A[INDEX_5],
                             HalfLife_5_EEDko_2B[INDEX_1],   HalfLife_5_EEDko_2B[INDEX_2],    HalfLife_5_EEDko_2B[INDEX_3],    HalfLife_5_EEDko_2B[INDEX_4],    HalfLife_5_EEDko_2B[INDEX_5] 
                           ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),  rep("EEDheto_Medium", length(INDEX_3)),  rep("EEDheto_High", length(INDEX_4)),  rep("EEDheto_Highest", length(INDEX_5)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),    rep("EEDko_Medium", length(INDEX_3)),    rep("EEDko_High", length(INDEX_4)),    rep("EEDko_Highest", length(INDEX_5))
                  ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",     "EEDko_Medium",    "EEDko_High",    "EEDko_Highest"  
                  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2",    "EEDko_Highest"="purple2" 
                  ), 
                  path2=subdir_5_part4,   fileName2= paste("4-5-2G-WT-HalfLife-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=40)    ## width = 1 + 10*0.5, height=5cm 




MyBoxViolinPlot_1(vector2=c( myNTR_5_EEDheto_2A[INDEX_1], myNTR_5_EEDheto_2A[INDEX_2],  myNTR_5_EEDheto_2A[INDEX_3],  myNTR_5_EEDheto_2A[INDEX_4],  myNTR_5_EEDheto_2A[INDEX_5],
                             myNTR_5_EEDko_2B[INDEX_1],   myNTR_5_EEDko_2B[INDEX_2],    myNTR_5_EEDko_2B[INDEX_3],    myNTR_5_EEDko_2B[INDEX_4],    myNTR_5_EEDko_2B[INDEX_5] 
                           ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),  rep("EEDheto_Medium", length(INDEX_3)),  rep("EEDheto_High", length(INDEX_4)),  rep("EEDheto_Highest", length(INDEX_5)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),    rep("EEDko_Medium", length(INDEX_3)),    rep("EEDko_High", length(INDEX_4)),    rep("EEDko_Highest", length(INDEX_5))
                  ), 
                  sampleRank2=c( "EEDheto_Lowest",   "EEDko_Lowest",   
                                 "EEDheto_Low",      "EEDko_Low", 
                                 "EEDheto_Medium",   "EEDko_Medium", 
                                 "EEDheto_High",     "EEDko_High",  
                                 "EEDheto_Highest",  "EEDko_Highest" 
                  ),     
                  colours2=c( "EEDheto_Lowest"="red2",        "EEDko_Lowest"="red2",   
                              "EEDheto_Low"="cyan2",          "EEDko_Low"="cyan2",  
                              "EEDheto_Medium"="green2",      "EEDko_Medium"="green2", 
                              "EEDheto_High"="blue2",         "EEDko_High"="blue2",  
                              "EEDheto_Highest"="purple2" ,   "EEDko_Highest"="purple2" 
                  ), 
                  path2=subdir_5_part4,   fileName2= paste("4-5-2H-WT-NTR-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0,   Ymax2=0.3)    ## width = 1 + 10*0.5, height=5cm 















############################# TAC 2 samples
sink( file=paste(subdir_5_part4,  "/4-5-3A-1-banding-runLog.txt",      sep = "") )
myNTR_5_banding_3A  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(row_Average_week0[i]+0.0001,  row_Average_banding[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  2.5),  yAxis2=log( vec1 )  , file2=paste(subdir_5_part4,  "/4-5-3A-2-banding-LogLinearModel", sep="") )
  myNTR_5_banding_3A <- c(myNTR_5_banding_3A, NTR_bin)
}
sink()  
length(myNTR_5_banding_3A)
summary(myNTR_5_banding_3A)
HalfLife_5_banding_3A <- log(2)/(myNTR_5_banding_3A+0.0001)
length(HalfLife_5_banding_3A)
summary(HalfLife_5_banding_3A)

write.table(x=myNTR_5_banding_3A,      file = paste(subdir_5_part4,  "/4-5-3A-3-banding-raw-myNTR.txt",  sep = ""),      append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_5_banding_3A,   file = paste(subdir_5_part4,  "/4-5-3A-4-banding-raw-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_5_banding_3A <- data.frame( xAxis = c(1:length(myNTR_5_banding_3A)),      yAxis = sort(myNTR_5_banding_3A) )
FigureTemp_5_banding_3A <- ggplot(myframe_5_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_banding_3A,  path1=subdir_5_part4, fileName1="4-5-3A-1-figure-banding-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_5_banding_3A <- ggplot(myframe_5_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_banding_3A,  path1=subdir_5_part4, fileName1="4-5-3A-1-figure-banding-Line-raw",  height1=4, width1=7)

myNTR_5_banding_3A[myNTR_5_banding_3A<0]  <- 0
myNTR_5_banding_3A[myNTR_5_banding_3A>1]  <- 1
HalfLife_5_banding_3A[HalfLife_5_banding_3A<0]  <- 0
HalfLife_5_banding_3A[HalfLife_5_banding_3A>50]  <- 50
min(myNTR_5_banding_3A)

write.table(x=myNTR_5_banding_3A,      file = paste(subdir_5_part4,  "/4-5-3A-5-noLess0-myNTR_5_banding_3A.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_5_banding_3A,   file = paste(subdir_5_part4,  "/4-5-3A-6-noLess0-HalfLife_5_banding_3A.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_5_banding_3A <- data.frame( xAxis = c(1:length(myNTR_5_banding_3A)),      yAxis = sort(myNTR_5_banding_3A) )
FigureTemp_5_banding_3A <- ggplot(myframe_5_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_banding_3A,  path1=subdir_5_part4, fileName1="4-5-3A-2-figure-banding-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_5_banding_3A <- ggplot(myframe_5_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_banding_3A,  path1=subdir_5_part4, fileName1="4-5-3A-2-figure-banding-Line-noLess0",  height1=4, width1=7)



sink( file=paste(subdir_5_part4,  "/4-5-3B-1-sham-runLog.txt",      sep = "") )
myNTR_5_sham_3B  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(row_Average_week0[i]+0.0001,  row_Average_sham[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  2.5),  yAxis2=log( vec1 )  , file2=paste(subdir_5_part4,  "/4-5-3B-2-sham-LogLinearModel", sep="") )
  myNTR_5_sham_3B <- c(myNTR_5_sham_3B, NTR_bin)
}
sink()  
length(myNTR_5_sham_3B)
summary(myNTR_5_sham_3B)
HalfLife_5_sham_3B <- log(2)/(myNTR_5_sham_3B+0.0001)
length(HalfLife_5_sham_3B)
summary(HalfLife_5_sham_3B)

write.table(x=myNTR_5_sham_3B,      file = paste(subdir_5_part4,  "/4-5-3B-3-sham-raw-NTR.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_5_sham_3B,   file = paste(subdir_5_part4,  "/4-5-3B-4-sham-raw-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_5_sham_3B <- data.frame( xAxis = c(1:length(myNTR_5_sham_3B)),      yAxis = sort(myNTR_5_sham_3B) )
FigureTemp_5_sham_3B <- ggplot(myframe_5_sham_3B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_sham_3B,  path1=subdir_5_part4, fileName1="4-5-3B-1-figure-sham-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_5_sham_3B <- ggplot(myframe_5_sham_3B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_5_sham_3B,  path1=subdir_5_part4, fileName1="4-5-3B-1-figure-sham-Line-raw",  height1=4, width1=7)

myNTR_5_sham_3B[myNTR_5_sham_3B<0]  <- 0
myNTR_5_sham_3B[myNTR_5_sham_3B>1]  <- 1
HalfLife_5_sham_3B[HalfLife_5_sham_3B<0]   <- 0
HalfLife_5_sham_3B[HalfLife_5_sham_3B>50]  <- 50

write.table(x=myNTR_5_sham_3B,      file = paste(subdir_5_part4,  "/4-5-3B-5-sham-noLess0-NTR.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_5_sham_3B,   file = paste(subdir_5_part4,  "/4-5-3B-6-sham-noLess0-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
MyDistributionLine_1(x2=myNTR_5_sham_3B,  path2=subdir_5_part4,  fileName2="4-5-3B-2-figure-sham-Line-noLess0",  
                     title2=myTitle_g,  xLab2=myTitle_g,  yLab2="NTR (sham)",   
                     Ymin2=0,   Ymax2=1,   height2=4,  width2=7, yintercept2=c(0.2, 0.2, 0.2)  )


MyBoxViolinPlot_1(vector2=c(myNTR_5_banding_3A,  myNTR_5_sham_3B),       
                  sampleType2=c( rep("banding", numOfRows1)  ,   rep("sham", numOfRows1)  ), 
                  sampleRank2=c( "banding",  "sham" ),     
                  colours2=c( "banding"="blue",  "sham"="green4"  ), 
                  path2=subdir_5_part4,     fileName2="4-5-3A-merge-TAC-NTR-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                  Ymin2=0,   Ymax2=1.0,    height2=5,   width2=2  )

MyBoxViolinPlot_1(vector2=c(HalfLife_5_banding_3A,  HalfLife_5_sham_3B),       
                  sampleType2=c( rep("banding", numOfRows1)  ,   rep("sham", numOfRows1)  ), 
                  sampleRank2=c( "banding",  "sham" ),     
                  colours2=c( "banding"="blue",  "sham"="green4"  ), 
                  path2=subdir_5_part4,     fileName2="4-5-3B-merge-TAC-HalfLife-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half life (weeks)",   
                  Ymin2=0,   Ymax2=6,    height2=5,   width2=2  )




diff1_5_TAC = myNTR_5_banding_3A - myNTR_5_sham_3B 
diff2_5_TAC = log2( ( (myNTR_5_banding_3A+ 0.001) / (myNTR_5_sham_3B+0.001))  )
diff3_5_TAC = log2( ( (myNTR_5_sham_3B+ 0.001) / (myNTR_5_banding_3A+0.001))  )
length(diff1_5_TAC)
length(diff2_5_TAC)
length(diff3_5_TAC)

MyDistributionLine_1(x2=diff1_5_TAC,  path2=subdir_5_part4,  fileName2="4-5-3C-merge1-figure-TAC-Line-minus",  
                     title2=myTitle_g,  xLab2=myTitle_g,  yLab2="banding-sham",   
                     Ymin2=-1,   Ymax2=1,   height2=4,  width2=7, yintercept2=c(-0.2, 0, 0.2)  )

MyDistributionLine_1(x2=diff2_5_TAC,  path2=subdir_5_part4,  fileName2="4-5-3C-merge2-figure-TAC-Line-ratio",  
                     title2=myTitle_g,  xLab2=myTitle_g,  yLab2="log2(banding/sham)",   
                     Ymin2=-5,   Ymax2=5,   height2=4,  width2=7, yintercept2=c(-1, 0, 1)  )

MyDistributionLine_1(x2=diff3_5_TAC,  path2=subdir_5_part4,  fileName2="4-5-3C-merge3-figure-TAC-Line-ratio",  
                     title2=myTitle_g,  xLab2=myTitle_g,  yLab2="log2(sham/banding)",   
                     Ymin2=-5,   Ymax2=5,   height2=4,  width2=7, yintercept2=c(-1, 0, 1)  )





MyBoxViolinPlot_1(vector2=c( myNTR_5_banding_3A[INDEX_1], myNTR_5_banding_3A[INDEX_2],  myNTR_5_banding_3A[INDEX_3],  myNTR_5_banding_3A[INDEX_4],  myNTR_5_banding_3A[INDEX_5],
                             myNTR_5_sham_3B[INDEX_1],    myNTR_5_sham_3B[INDEX_2],     myNTR_5_sham_3B[INDEX_3],     myNTR_5_sham_3B[INDEX_4],     myNTR_5_sham_3B[INDEX_5] 
),   
sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),  rep("banding_Medium", length(INDEX_3)),  rep("banding_High", length(INDEX_4)),  rep("banding_Highest", length(INDEX_5)), 
               rep("sham_Lowest", length(INDEX_1)),     rep("sham_Low", length(INDEX_2)),     rep("sham_Medium", length(INDEX_3)),     rep("sham_High", length(INDEX_4)),     rep("sham_Highest", length(INDEX_5))
), 
sampleRank2=c( "banding_Lowest",  "banding_Low",   "banding_Medium",  "banding_High",  "banding_Highest"  ,
               "sham_Lowest",     "sham_Low",      "sham_Medium",     "sham_High",     "sham_Highest"  
),     
colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  "banding_Highest"="purple2" ,
            "sham_Lowest"="red2",      "sham_Low"="cyan2",       "sham_Medium"="green2",     "sham_High"="blue2",     "sham_Highest"="purple2" 
), 
path2=subdir_5_part4,   fileName2= paste("4-5-3F-WT-NTR-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
height2=4.0,   width2=6,   Ymin2=0, Ymax2=1)    ## width = 1 + 10*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_5_banding_3A[INDEX_1], HalfLife_5_banding_3A[INDEX_2],  HalfLife_5_banding_3A[INDEX_3],  HalfLife_5_banding_3A[INDEX_4],  HalfLife_5_banding_3A[INDEX_5],
                             HalfLife_5_sham_3B[INDEX_1],    HalfLife_5_sham_3B[INDEX_2],     HalfLife_5_sham_3B[INDEX_3],     HalfLife_5_sham_3B[INDEX_4],     HalfLife_5_sham_3B[INDEX_5] 
),   
sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),  rep("banding_Medium", length(INDEX_3)),  rep("banding_High", length(INDEX_4)),  rep("banding_Highest", length(INDEX_5)), 
               rep("sham_Lowest", length(INDEX_1)),    rep("sham_Low", length(INDEX_2)),    rep("sham_Medium", length(INDEX_3)),    rep("sham_High", length(INDEX_4)),    rep("sham_Highest", length(INDEX_5))
), 
sampleRank2=c( "banding_Lowest",  "banding_Low",   "banding_Medium",  "banding_High",  "banding_Highest"  ,
               "sham_Lowest",    "sham_Low",     "sham_Medium",    "sham_High",    "sham_Highest"  
),     
colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  "banding_Highest"="purple2" ,
            "sham_Lowest"="red2",     "sham_Low"="cyan2",      "sham_Medium"="green2",    "sham_High"="blue2",    "sham_Highest"="purple2" 
), 
path2=subdir_5_part4,   fileName2= paste("4-5-3G-WT-HalfLife-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
height2=4.0,   width2=6,   Ymin2=0, Ymax2=10)    ## width = 1 + 10*0.5, height=5cm 




MyBoxViolinPlot_1(vector2=c( myNTR_5_banding_3A[INDEX_1], myNTR_5_banding_3A[INDEX_2],  myNTR_5_banding_3A[INDEX_3],  myNTR_5_banding_3A[INDEX_4],  myNTR_5_banding_3A[INDEX_5],
                             myNTR_5_sham_3B[INDEX_1],   myNTR_5_sham_3B[INDEX_2],    myNTR_5_sham_3B[INDEX_3],    myNTR_5_sham_3B[INDEX_4],    myNTR_5_sham_3B[INDEX_5] 
),   
sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),  rep("banding_Medium", length(INDEX_3)),  rep("banding_High", length(INDEX_4)),  rep("banding_Highest", length(INDEX_5)), 
               rep("sham_Lowest", length(INDEX_1)),    rep("sham_Low", length(INDEX_2)),    rep("sham_Medium", length(INDEX_3)),    rep("sham_High", length(INDEX_4)),    rep("sham_Highest", length(INDEX_5))
), 
sampleRank2=c( "banding_Lowest",   "sham_Lowest",   
               "banding_Low",      "sham_Low", 
               "banding_Medium",   "sham_Medium", 
               "banding_High",     "sham_High",  
               "banding_Highest",  "sham_Highest" 
),     
colours2=c( "banding_Lowest"="red2",        "sham_Lowest"="red2",   
            "banding_Low"="cyan2",          "sham_Low"="cyan2",  
            "banding_Medium"="green2",      "sham_Medium"="green2", 
            "banding_High"="blue2",         "sham_High"="blue2",  
            "banding_Highest"="purple2" ,   "sham_Highest"="purple2" 
), 
path2=subdir_5_part4,   fileName2= paste("4-5-3H-WT-NTR-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
height2=4.0,   width2=6,   Ymin2=0,   Ymax2=1)    ## width = 1 + 10*0.5, height=5cm 

























#########################################################################
subdir_5A_part4 <- paste(Part4_g,  "/5A-NTR-eachRow-5kb", sep = "")
if( ! file.exists(subdir_5A_part4) ) { dir.create(subdir_5A_part4) }


dim(Average_H3)
dim(Average_week0)
numRows_oneClass <- floor( nrow(Average_H3)/4 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
INDEX_4 <- seq(from = numRows_oneClass*3+1, to = numRows_oneClass*4, by =1 )
nrow(Average_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) + length(INDEX_4) )

MyBoxViolinPlot_1(vector2=c( myNTR_5_WT_1[INDEX_1], myNTR_5_WT_1[INDEX_2],  myNTR_5_WT_1[INDEX_3],  myNTR_5_WT_1[INDEX_4]  ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),  rep("High", length(INDEX_3)),  
                                 rep("Highest", length(INDEX_4))   ), 
                  sampleRank2=c( "Lowest",  "Low",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",    "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_5A_part4,   fileName2= paste("4-5A-1E-WT-NTR-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=3.0,   Ymin2=0, Ymax2=0.4)    ## width = 1 + 4*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c( HalfLife_5_WT_1[INDEX_1], HalfLife_5_WT_1[INDEX_2],  HalfLife_5_WT_1[INDEX_3],  HalfLife_5_WT_1[INDEX_4]  ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),      
                                 rep("High", length(INDEX_3)),      rep("Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "Lowest",  "Low",   "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",    "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_5A_part4,   fileName2= paste("4-5A-1F-WT-HalfLife-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=3.0,   Ymin2=0, Ymax2=40)    ## width = 1 + 4*0.5, height=5cm 











MyBoxViolinPlot_1(vector2=c( myNTR_5_EEDheto_2A[INDEX_1], myNTR_5_EEDheto_2A[INDEX_2],  myNTR_5_EEDheto_2A[INDEX_3],  myNTR_5_EEDheto_2A[INDEX_4],  
                             myNTR_5_EEDko_2B[INDEX_1],   myNTR_5_EEDko_2B[INDEX_2],    myNTR_5_EEDko_2B[INDEX_3],    myNTR_5_EEDko_2B[INDEX_4] 
                           ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),  rep("EEDheto_High", length(INDEX_3)),  rep("EEDheto_Highest", length(INDEX_4)),  
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),    rep("EEDko_High", length(INDEX_3)),    rep("EEDko_Highest", length(INDEX_4)) 
                  ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",      "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",        "EEDko_High",    "EEDko_Highest"  
                  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",     "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",       "EEDko_High"="blue2",    "EEDko_Highest"="purple2" 
                  ), 
                  path2=subdir_5A_part4,   fileName2= paste("4-5A-2F-WT-NTR-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=0.3)    ## width = 1 + 8*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_5_EEDheto_2A[INDEX_1], HalfLife_5_EEDheto_2A[INDEX_2],  HalfLife_5_EEDheto_2A[INDEX_3],  HalfLife_5_EEDheto_2A[INDEX_4],  
                             HalfLife_5_EEDko_2B[INDEX_1],   HalfLife_5_EEDko_2B[INDEX_2],    HalfLife_5_EEDko_2B[INDEX_3],    HalfLife_5_EEDko_2B[INDEX_4] 
                           ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),       rep("EEDheto_High", length(INDEX_3)),  rep("EEDheto_Highest", length(INDEX_4)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),         rep("EEDko_High", length(INDEX_3)),    rep("EEDko_Highest", length(INDEX_4))
                  ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",      "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",        "EEDko_High",    "EEDko_Highest"  
                  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",        "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",          "EEDko_High"="blue2",    "EEDko_Highest"="purple2" 
                  ), 
                  path2=subdir_5A_part4,   fileName2= paste("4-5A-2G-WT-HalfLife-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=40)    ## width = 1 + 8*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( myNTR_5_EEDheto_2A[INDEX_1], myNTR_5_EEDheto_2A[INDEX_2],  myNTR_5_EEDheto_2A[INDEX_3],  myNTR_5_EEDheto_2A[INDEX_4],  
                             myNTR_5_EEDko_2B[INDEX_1],   myNTR_5_EEDko_2B[INDEX_2],    myNTR_5_EEDko_2B[INDEX_3],    myNTR_5_EEDko_2B[INDEX_4] 
                           ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),       rep("EEDheto_High", length(INDEX_3)),  rep("EEDheto_Highest", length(INDEX_4)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),         rep("EEDko_High", length(INDEX_3)),    rep("EEDko_Highest", length(INDEX_4))
                  ), 
                  sampleRank2=c( "EEDheto_Lowest",   "EEDko_Lowest",   
                                 "EEDheto_Low",      "EEDko_Low", 
                                 "EEDheto_High",     "EEDko_High",  
                                 "EEDheto_Highest",  "EEDko_Highest" 
                  ),     
                  colours2=c( "EEDheto_Lowest"="red2",        "EEDko_Lowest"="red2",   
                              "EEDheto_Low"="cyan2",          "EEDko_Low"="cyan2",  
                              "EEDheto_High"="blue2",         "EEDko_High"="blue2",  
                              "EEDheto_Highest"="purple2" ,   "EEDko_Highest"="purple2" 
                  ), 
                  path2=subdir_5A_part4,   fileName2= paste("4-5A-2H-WT-NTR-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0,   Ymax2=0.3)    ## width = 1 + 10*0.5, height=5cm 










MyBoxViolinPlot_1(vector2=c( myNTR_5_banding_3A[INDEX_1], myNTR_5_banding_3A[INDEX_2],  myNTR_5_banding_3A[INDEX_3],  myNTR_5_banding_3A[INDEX_4] ,
                             myNTR_5_sham_3B[INDEX_1],    myNTR_5_sham_3B[INDEX_2],     myNTR_5_sham_3B[INDEX_3],     myNTR_5_sham_3B[INDEX_4] 
                            ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),      rep("banding_High", length(INDEX_3)),  rep("banding_Highest", length(INDEX_4)), 
                                 rep("sham_Lowest", length(INDEX_1)),     rep("sham_Low", length(INDEX_2)),         rep("sham_High", length(INDEX_3)),     rep("sham_Highest", length(INDEX_4))
                  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",       "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",     "sham_Low",          "sham_High",     "sham_Highest"  
                  ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",       "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",      "sham_Low"="cyan2",          "sham_High"="blue2",     "sham_Highest"="purple2" 
                  ), 
                  path2=subdir_5A_part4,   fileName2= paste("4-5A-3F-WT-NTR-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=1)    ## width = 1 + 10*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_5_banding_3A[INDEX_1], HalfLife_5_banding_3A[INDEX_2],  HalfLife_5_banding_3A[INDEX_3],  HalfLife_5_banding_3A[INDEX_4], 
                             HalfLife_5_sham_3B[INDEX_1],    HalfLife_5_sham_3B[INDEX_2],     HalfLife_5_sham_3B[INDEX_3],     HalfLife_5_sham_3B[INDEX_4]
                          ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),     rep("banding_High", length(INDEX_3)),   rep("banding_Highest", length(INDEX_4)), 
                                 rep("sham_Lowest", length(INDEX_1)),     rep("sham_Low", length(INDEX_2)),        rep("sham_High", length(INDEX_3)),      rep("sham_Highest", length(INDEX_4))
                  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",      "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",    "sham_Low",          "sham_High",    "sham_Highest"  
                  ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",      "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",     "sham_Low"="cyan2",           "sham_High"="blue2",    "sham_Highest"="purple2" 
                  ), 
                  path2=subdir_5A_part4,   fileName2= paste("4-5A-3G-WT-HalfLife-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=10)    ## width = 1 + 10*0.5, height=5cm 




MyBoxViolinPlot_1(vector2=c( myNTR_5_banding_3A[INDEX_1], myNTR_5_banding_3A[INDEX_2],  myNTR_5_banding_3A[INDEX_3],  myNTR_5_banding_3A[INDEX_4],  
                             myNTR_5_sham_3B[INDEX_1],   myNTR_5_sham_3B[INDEX_2],    myNTR_5_sham_3B[INDEX_3],    myNTR_5_sham_3B[INDEX_4]  
                  ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),     rep("banding_High", length(INDEX_3)),  rep("banding_Highest", length(INDEX_4)), 
                                 rep("sham_Lowest", length(INDEX_1)),    rep("sham_Low", length(INDEX_2)),         rep("sham_High", length(INDEX_3)),    rep("sham_Highest", length(INDEX_4))
                  ), 
                  sampleRank2=c( "banding_Lowest",   "sham_Lowest",   
                                 "banding_Low",      "sham_Low", 
                                 "banding_High",     "sham_High",  
                                 "banding_Highest",  "sham_Highest" 
                  ),     
                  colours2=c( "banding_Lowest"="red2",        "sham_Lowest"="red2",   
                              "banding_Low"="cyan2",          "sham_Low"="cyan2",  
                              "banding_High"="blue2",         "sham_High"="blue2",  
                              "banding_Highest"="purple2" ,   "sham_Highest"="purple2" 
                  ), 
                  path2=subdir_5A_part4,   fileName2= paste("4-5A-3H-WT-NTR-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0,   Ymax2=1)    ## width = 1 + 10*0.5, height=5cm 





























#########################################################################
subdir_5B_part4 <- paste(Part4_g,  "/5B-NTR-eachRow-5kb", sep = "")
if( ! file.exists(subdir_5B_part4) ) { dir.create(subdir_5B_part4) }

dim(Average_H3)
dim(Average_week0)
numRows_oneClass <- floor( nrow(Average_H3)/3 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
nrow(Average_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3)  )

MyBoxViolinPlot_1(vector2=c( myNTR_5_WT_1[INDEX_1], myNTR_5_WT_1[INDEX_2],  myNTR_5_WT_1[INDEX_3]   ),   
                  sampleType2=c( rep("Low", length(INDEX_1)),    rep("Medium", length(INDEX_2)),  rep("High", length(INDEX_3)) ), 
                  sampleRank2=c(   "Low", "Medium",  "High"   ),     
                  colours2=c(  "Low"="cyan2",   "Medium"="red2", "High"="blue2"  ), 
                  path2=subdir_5B_part4,   fileName2= paste("4-5B-1E-WT-NTR-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=2.5,   Ymin2=0, Ymax2=0.4)    ## width = 1 + 3*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c( HalfLife_5_WT_1[INDEX_1], HalfLife_5_WT_1[INDEX_2],  HalfLife_5_WT_1[INDEX_3]   ),   
                  sampleType2=c( rep("Low", length(INDEX_1)),    rep("Medium", length(INDEX_2)),      rep("High", length(INDEX_3))   ), 
                  sampleRank2=c( "Low",   "Medium",  "High"  ),     
                  colours2=c( "Low"="cyan2",    "Medium"="red2",   "High"="blue2"  ), 
                  path2=subdir_5B_part4,   fileName2= paste("4-5B-1F-WT-HalfLife-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=2.5,   Ymin2=0, Ymax2=40)    ## width = 1 + 4*0.5, height=5cm 






MyBoxViolinPlot_1(vector2=c( myNTR_5_EEDheto_2A[INDEX_1], myNTR_5_EEDheto_2A[INDEX_2],  myNTR_5_EEDheto_2A[INDEX_3],  
                             myNTR_5_EEDko_2B[INDEX_1],   myNTR_5_EEDko_2B[INDEX_2],    myNTR_5_EEDko_2B[INDEX_3]     
                  ),   
                  sampleType2=c( rep("EEDheto_Low", length(INDEX_1)),  rep("EEDheto_Medium", length(INDEX_2)),  rep("EEDheto_High", length(INDEX_3)),    
                                 rep("EEDko_Low", length(INDEX_1)),    rep("EEDko_Medium", length(INDEX_2)),    rep("EEDko_High", length(INDEX_3)) 
                  ), 
                  sampleRank2=c( "EEDheto_Low",  "EEDheto_Medium",      "EEDheto_High",  
                                 "EEDko_Low",    "EEDko_Medium",        "EEDko_High"  
                  ),     
                  colours2=c( "EEDheto_Low"="red2",   "EEDheto_Medium"="cyan2",     "EEDheto_High"="blue2",  
                              "EEDko_Low"="red2",     "EEDko_Medium"="cyan2",       "EEDko_High"="blue2" 
                  ), 
                  path2=subdir_5B_part4,   fileName2= paste("4-5B-2F-WT-NTR-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=0.3)    ## width = 1 + 8*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_5_EEDheto_2A[INDEX_1], HalfLife_5_EEDheto_2A[INDEX_2],  HalfLife_5_EEDheto_2A[INDEX_3],    
                             HalfLife_5_EEDko_2B[INDEX_1],   HalfLife_5_EEDko_2B[INDEX_2],    HalfLife_5_EEDko_2B[INDEX_3] 
                  ),   
                  sampleType2=c( rep("EEDheto_Low", length(INDEX_1)),  rep("EEDheto_Medium", length(INDEX_2)),       rep("EEDheto_High", length(INDEX_3)),   
                                 rep("EEDko_Low", length(INDEX_1)),    rep("EEDko_Medium", length(INDEX_2)),         rep("EEDko_High", length(INDEX_3)) 
                  ), 
                  sampleRank2=c( "EEDheto_Low",  "EEDheto_Medium",      "EEDheto_High",  
                                 "EEDko_Low",    "EEDko_Medium",        "EEDko_High"  
                  ),     
                  colours2=c( "EEDheto_Low"="red2",   "EEDheto_Medium"="cyan2",        "EEDheto_High"="blue2",   
                              "EEDko_Low"="red2",     "EEDko_Medium"="cyan2",          "EEDko_High"="blue2" 
                  ), 
                  path2=subdir_5B_part4,   fileName2= paste("4-5B-2G-WT-HalfLife-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=40)    ## width = 1 + 8*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( myNTR_5_EEDheto_2A[INDEX_1], myNTR_5_EEDheto_2A[INDEX_2],  myNTR_5_EEDheto_2A[INDEX_3],   
                             myNTR_5_EEDko_2B[INDEX_1],   myNTR_5_EEDko_2B[INDEX_2],    myNTR_5_EEDko_2B[INDEX_3] 
                  ),   
                  sampleType2=c( rep("EEDheto_Low", length(INDEX_1)),  rep("EEDheto_Medium", length(INDEX_2)),       rep("EEDheto_High", length(INDEX_3)),   
                                 rep("EEDko_Low", length(INDEX_1)),    rep("EEDko_Medium", length(INDEX_2)),         rep("EEDko_High", length(INDEX_3)) 
                  ), 
                  sampleRank2=c( "EEDheto_Low",      "EEDko_Low",   
                                 "EEDheto_Medium",   "EEDko_Medium", 
                                 "EEDheto_High",     "EEDko_High"  
                  ),     
                  colours2=c( "EEDheto_Low"="red2",           "EEDko_Low"="red2",   
                              "EEDheto_Medium"="cyan2",       "EEDko_Medium"="cyan2",  
                              "EEDheto_High"="blue2",         "EEDko_High"="blue2" 
                  ), 
                  path2=subdir_5B_part4,   fileName2= paste("4-5B-2H-WT-NTR-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0,   Ymax2=0.3)    ## width = 1 + 10*0.5, height=5cm 






MyBoxViolinPlot_1(vector2=c( myNTR_5_banding_3A[INDEX_1], myNTR_5_banding_3A[INDEX_2],  myNTR_5_banding_3A[INDEX_3],   
                             myNTR_5_sham_3B[INDEX_1],    myNTR_5_sham_3B[INDEX_2],     myNTR_5_sham_3B[INDEX_3]  
                  ),   
                  sampleType2=c( rep("banding_Low", length(INDEX_1)),  rep("banding_Medium", length(INDEX_2)),      rep("banding_High", length(INDEX_3)),   
                                 rep("sham_Low", length(INDEX_1)),     rep("sham_Medium", length(INDEX_2)),         rep("sham_High", length(INDEX_3))
                  ), 
                  sampleRank2=c( "banding_Low",  "banding_Medium",       "banding_High",  
                                 "sham_Low",     "sham_Medium",          "sham_High"  
                  ),     
                  colours2=c(  "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2", 
                               "sham_Low"="cyan2",       "sham_Medium"="green2",     "sham_High"="blue2" 
                  ), 
                  path2=subdir_5B_part4,   fileName2= paste("4-5B-3F-WT-NTR-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=1)    ## width = 1 + 10*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_5_banding_3A[INDEX_1], HalfLife_5_banding_3A[INDEX_2],  HalfLife_5_banding_3A[INDEX_3],   
                             HalfLife_5_sham_3B[INDEX_1],    HalfLife_5_sham_3B[INDEX_2],     HalfLife_5_sham_3B[INDEX_3]
                   ),   
                  sampleType2=c( rep("banding_Low", length(INDEX_1)),  rep("banding_Medium", length(INDEX_2)),  rep("banding_High", length(INDEX_3)),  
                                 rep("sham_Low", length(INDEX_1)),    rep("sham_Medium", length(INDEX_2)),    rep("sham_High", length(INDEX_3))
                  ), 
                  sampleRank2=c( "banding_Low",   "banding_Medium",  "banding_High",  
                                 "sham_Low",     "sham_Medium",    "sham_High"  
                  ),     
                  colours2=c( "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  
                              "sham_Low"="cyan2",       "sham_Medium"="green2",     "sham_High"="blue2" 
                  ), 
                  path2=subdir_5B_part4,   fileName2= paste("4-5B-3G-WT-HalfLife-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=10)    ## width = 1 + 10*0.5, height=5cm 




MyBoxViolinPlot_1(vector2=c( myNTR_5_banding_3A[INDEX_1], myNTR_5_banding_3A[INDEX_2],  myNTR_5_banding_3A[INDEX_3],  
                             myNTR_5_sham_3B[INDEX_1],    myNTR_5_sham_3B[INDEX_2],     myNTR_5_sham_3B[INDEX_3] 
                  ),   
                  sampleType2=c( rep("banding_Low", length(INDEX_1)),    rep("banding_Medium", length(INDEX_2)),  rep("banding_High", length(INDEX_3)),  
                                 rep("sham_Low",    length(INDEX_1)),    rep("sham_Medium", length(INDEX_2)),     rep("sham_High",    length(INDEX_3))
                  ), 
                  sampleRank2=c( "banding_Low",      "sham_Low", 
                                 "banding_Medium",   "sham_Medium", 
                                 "banding_High",     "sham_High" 
                  ),     
                  colours2=c( "banding_Low"="cyan2",          "sham_Low"="cyan2",  
                              "banding_Medium"="green2",      "sham_Medium"="green2", 
                              "banding_High"="blue2",         "sham_High"="blue2"  
                  ), 
                  path2=subdir_5B_part4,   fileName2= paste("4-5B-3H-WT-NTR-BoxViolin",   "5kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0,   Ymax2=1)    ## width = 1 + 10*0.5, height=5cm 




















#########################################################################################################  average some columns for each row
subdir_6_part4 <- paste(Part4_g,  "/6-NTR-eachRow-2kb", sep = "")
if( ! file.exists(subdir_6_part4) ) { dir.create(subdir_6_part4) }


############################# WT 6 samples
sink( file=paste(subdir_6_part4,  "/4-6-1A-WT-runLog.txt",      sep = "") )
myNTR_6_WT_1  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(reduceColumn2_Average_week0[i]+0.0006, reduceColumn2_Average_week1[i]+0.0005, reduceColumn2_Average_week2[i]+0.0004, reduceColumn2_Average_week4[i]+0.0003, reduceColumn2_Average_week6[i]+0.0002, reduceColumn2_Average_week8[i]+0.0001)                        
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log(vec1)  , file2=paste(subdir_6_part4,  "/4-6-1B-WT-LogLinearModel", sep="") )
  myNTR_6_WT_1 <- c(myNTR_6_WT_1, NTR_bin)
}
sink()  

length(myNTR_6_WT_1)
summary(myNTR_6_WT_1)
HalfLife_6_WT_1 <- log(2)/(myNTR_6_WT_1+0.0001)
length(HalfLife_6_WT_1)
summary(HalfLife_6_WT_1)

write.table(x=myNTR_6_WT_1,      file = paste(subdir_6_part4,  "/4-6-1C-raw-NTR.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_6_WT_1,   file = paste(subdir_6_part4,  "/4-6-1D-raw-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_6_WT_1 <- data.frame( xAxis = c(1:length(myNTR_6_WT_1)),      yAxis = sort(myNTR_6_WT_1) )
FigureTemp_6_WT_1 <- ggplot(myframe_6_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_WT_1,  path1=subdir_6_part4, fileName1="4-6-1A-figure-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_6_WT_1 <- ggplot(myframe_6_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) +  
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_WT_1,  path1=subdir_6_part4, fileName1="4-6-1A-figure-Line-raw",  height1=4, width1=7)

myNTR_6_WT_1[myNTR_6_WT_1<0]  <- 0
myNTR_6_WT_1[myNTR_6_WT_1>1]  <- 1
HalfLife_6_WT_1[HalfLife_6_WT_1<0]  <- 0
HalfLife_6_WT_1[HalfLife_6_WT_1>50] <- 50

write.table(x=thisRowNames1,     file = paste(subdir_6_part4,   "/4-6-1E-RowNames.txt",  sep = ""),         append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=myNTR_6_WT_1,      file = paste(subdir_6_part4,  "/4-6-1F-noLess0-NTR.txt",    sep = ""),     append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_6_WT_1,   file = paste(subdir_6_part4,  "/4-6-1G-noLess0-HalfLife.txt", sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_6_WT_1 <- data.frame( xAxis = c(1:length(myNTR_6_WT_1)),      yAxis = sort(myNTR_6_WT_1) )
FigureTemp_6_WT_1 <- ggplot(myframe_6_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_WT_1,  path1=subdir_6_part4, fileName1="4-6-1B-figure-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_6_WT_1 <- ggplot(myframe_6_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) +  
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_WT_1,  path1=subdir_6_part4, fileName1="4-6-1B-figure-Line-noLess0",  height1=4, width1=7)

MyBoxViolinPlot_1(vector2=myNTR_6_WT_1,      sampleType2=c( rep("NTR", numOfRows1)  ), 
                  sampleRank2=c( "NTR" ),     colours2=c( "NTR"="red" ), 
                  path2=subdir_6_part4,     fileName2="4-6-1C-WT-NTR-average-curve",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="NTR",   
                  Ymin2=0,   Ymax2=0.4,    height2=5,   width2=1.5  )

MyBoxViolinPlot_1(vector2=HalfLife_6_WT_1,      sampleType2=c( rep("HalfLife", numOfRows1)  ), 
                  sampleRank2=c( "HalfLife" ),     colours2=c( "HalfLife"="red" ), 
                  path2=subdir_6_part4,     fileName2="4-6-1D-WT-HalfLife-average-curve",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="Half life (weeks)",   
                  Ymin2=0,   Ymax2=20,    height2=5,   width2=1.5  )










############################# CKO 4 samples
sink( file=paste(subdir_6_part4,  "/4-6-2A-1-EEDheto-runLog.txt",      sep = "") )
myNTR_6_EEDheto_2A  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(reduceColumn2_Average_week0_EEDheto[i]+0.0001,  reduceColumn2_Average_week4_EEDheto[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_6_part4,  "/4-6-2A-2-EEDheto-LogLinearModel", sep="") )
  myNTR_6_EEDheto_2A <- c(myNTR_6_EEDheto_2A, NTR_bin)
}
sink()  

length(myNTR_6_EEDheto_2A)
summary(myNTR_6_EEDheto_2A)
HalfLife_6_EEDheto_2A <- log(2)/(myNTR_6_EEDheto_2A+0.0001)
length(HalfLife_6_EEDheto_2A)
summary(HalfLife_6_EEDheto_2A)

write.table(x=myNTR_6_EEDheto_2A,      file = paste(subdir_6_part4,  "/4-6-2A-3-raw-NTR_EEDheto.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_6_EEDheto_2A,   file = paste(subdir_6_part4,  "/4-6-2A-4-raw-HalfLife_EEDheto.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_6_EEDheto_2A <- data.frame( xAxis = c(1:length(myNTR_6_EEDheto_2A)),      yAxis = sort(myNTR_6_EEDheto_2A) )
FigureTemp_6_EEDheto_2A <- ggplot(myframe_6_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_EEDheto_2A,  path1=subdir_6_part4, fileName1="4-6-2A-1-figure-EEDheto-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_6_EEDheto_2A <- ggplot(myframe_6_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) +
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_EEDheto_2A,  path1=subdir_6_part4, fileName1="4-6-2A-1-figure-EEDheto-Line-raw",  height1=4, width1=7)

myNTR_6_EEDheto_2A[myNTR_6_EEDheto_2A<0]  <- 0
myNTR_6_EEDheto_2A[myNTR_6_EEDheto_2A>1]  <- 1
HalfLife_6_EEDheto_2A[HalfLife_6_EEDheto_2A<0]   <- 0
HalfLife_6_EEDheto_2A[HalfLife_6_EEDheto_2A>50]  <- 50

write.table(x=myNTR_6_EEDheto_2A,      file = paste(subdir_6_part4,  "/4-6-2A-6-noLess0-myNTR_EEDheto.txt",  sep = ""),      append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_6_EEDheto_2A,   file = paste(subdir_6_part4,  "/4-6-2A-6-noLess0-HalfLife_EEDheto.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_6_EEDheto_2A <- data.frame( xAxis = c(1:length(myNTR_6_EEDheto_2A)),      yAxis = sort(myNTR_6_EEDheto_2A) )
FigureTemp_6_EEDheto_2A <- ggplot(myframe_6_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_EEDheto_2A,  path1=subdir_6_part4, fileName1="4-6-2A-2-figure-EEDheto-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_6_EEDheto_2A <- ggplot(myframe_6_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_EEDheto_2A,  path1=subdir_6_part4, fileName1="4-6-2A-2-figure-EEDheto-Line-noLess0",  height1=4, width1=7)





sink( file=paste(subdir_6_part4,  "/4-6-2B-1-EEDko-runLog.txt",      sep = "") )
myNTR_6_EEDko_2B  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(reduceColumn2_Average_week0_EEDko[i]+0.0001,  reduceColumn2_Average_week4_EEDko[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_6_part4,  "/4-6-2B-2-EEDko-LogLinearModel", sep="") )
  myNTR_6_EEDko_2B <- c(myNTR_6_EEDko_2B, NTR_bin)
}
sink()  

length(myNTR_6_EEDko_2B)
summary(myNTR_6_EEDko_2B)
HalfLife_6_EEDko_2B <- log(2)/(myNTR_6_EEDko_2B+0.0001)
length(HalfLife_6_EEDko_2B)
summary(HalfLife_6_EEDko_2B)

write.table(x=myNTR_6_EEDko_2B,      file = paste(subdir_6_part4,  "/4-6-2B-3-raw-NTR_EEDko.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_6_EEDko_2B,   file = paste(subdir_6_part4,  "/4-6-2B-4-raw-HalfLife_EEDko.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_6_EEDko_2B <- data.frame( xAxis = c(1:length(myNTR_6_EEDko_2B)),      yAxis = sort(myNTR_6_EEDko_2B) )
FigureTemp_6_EEDko_2B <- ggplot(myframe_6_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_EEDko_2B,  path1=subdir_6_part4, fileName1="4-6-2B-1-figure-EEDko-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_6_EEDko_2B <- ggplot(myframe_6_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_EEDko_2B,  path1=subdir_6_part4, fileName1="4-6-2B-1-figure-EEDko-Line-raw",  height1=4, width1=7)

myNTR_6_EEDko_2B[myNTR_6_EEDko_2B<0]  <- 0
myNTR_6_EEDko_2B[myNTR_6_EEDko_2B>1]  <- 1
HalfLife_6_EEDko_2B[HalfLife_6_EEDko_2B<0]   <- 0
HalfLife_6_EEDko_2B[HalfLife_6_EEDko_2B>50]  <- 50

write.table(x=myNTR_6_EEDko_2B,      file = paste(subdir_6_part4,  "/4-6-2B-6-noLess0-myNTR_EEDko.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_6_EEDko_2B,   file = paste(subdir_6_part4,  "/4-6-2B-6-noLess0-HalfLife_EEDko.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_6_EEDko_2B <- data.frame( xAxis = c(1:length(myNTR_6_EEDko_2B)),      yAxis = sort(myNTR_6_EEDko_2B) )
FigureTemp_6_EEDko_2B <- ggplot(myframe_6_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_EEDko_2B,  path1=subdir_6_part4, fileName1="4-6-2B-2-figure-EEDko-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_6_EEDko_2B <- ggplot(myframe_6_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_EEDko_2B,  path1=subdir_6_part4, fileName1="4-6-2B-2-figure-EEDko-Line-noLess0",  height1=4, width1=7)



MyBoxViolinPlot_1(vector2=c(myNTR_6_EEDheto_2A,  myNTR_6_EEDko_2B),       
                  sampleType2=c( rep("EEDheto", numOfRows1)  ,   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "EEDheto",  "EEDko" ),     
                  colours2=c( "EEDheto"="red",  "EEDko"="blue"  ), 
                  path2=subdir_6_part4,     fileName2="4-6-2A-merge-CKO-NTR-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                  Ymin2=0,   Ymax2=0.3,    height2=5,   width2=2  )

MyBoxViolinPlot_1(vector2=c(HalfLife_6_EEDheto_2A,  HalfLife_6_EEDko_2B),       
                  sampleType2=c( rep("EEDheto", numOfRows1)  ,   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "EEDheto",  "EEDko" ),     
                  colours2=c( "EEDheto"="red",  "EEDko"="blue"  ), 
                  path2=subdir_6_part4,     fileName2="4-6-2B-merge-CKO-HalfLife-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half life (weeks)",   
                  Ymin2=0,   Ymax2=20,    height2=5,   width2=2  )



diff1_6_CKO = myNTR_6_EEDko_2B - myNTR_6_EEDheto_2A
diff2_6_CKO = log2( ( (myNTR_6_EEDko_2B+0.001) / (myNTR_6_EEDheto_2A+0.001))  )
diff3_6_CKO = log2( ( (myNTR_6_EEDheto_2A+0.001) / (myNTR_6_EEDko_2B+0.001))  )   ## 2^(-10) = 0.001
length(diff1_6_CKO)
length(diff2_6_CKO)
length(diff3_6_CKO)
summary(diff1_6_CKO)
summary(diff2_6_CKO)
summary(diff3_6_CKO)


myframe_6_CKO <- data.frame( xAxis = c(1:length(diff1_6_CKO)),      yAxis = sort(diff1_6_CKO) )

FigureTemp_6_CKO <- ggplot(myframe_6_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR (EEDko - EEDheto)") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.1, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -0.1, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_CKO,  path1=subdir_6_part4, fileName1="4-6-2C-merge-figure-CKO-Line-minus-limitY",  height1=4, width1=7)

FigureTemp_6_CKO <- ggplot(myframe_6_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR (EEDko - EEDheto)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.1, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -0.1, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_CKO,  path1=subdir_6_part4, fileName1="4-6-2C-merge-figure-CKO-Line-minus",  height1=4, width1=7)



myframe_6_CKO <- data.frame( xAxis = c(1:length(diff2_6_CKO)),      yAxis = sort(diff2_6_CKO) )

FigureTemp_6_CKO <- ggplot(myframe_6_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDko/EEDheto)") +  ggtitle(myTitle_g) + ylim(-10, 10 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_CKO,  path1=subdir_6_part4, fileName1="4-6-2D-merge-figure-CKO-Line-ratio-limitY",  height1=4, width1=7)

FigureTemp_6_CKO <- ggplot(myframe_6_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDko/EEDheto)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_CKO,  path1=subdir_6_part4, fileName1="4-6-2D-merge-figure-CKO-Line-ratio",  height1=4, width1=7)





myframe_6_CKO <- data.frame( xAxis = c(1:length(diff3_6_CKO)),      yAxis = sort(diff3_6_CKO) )

FigureTemp_6_CKO <- ggplot(myframe_6_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDheto/EEDko)") +  ggtitle(myTitle_g) + ylim(-10, 10 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_CKO,  path1=subdir_6_part4, fileName1="4-6-2E-merge-figure-CKO-Line-ratio-limitY",  height1=4, width1=7)

FigureTemp_6_CKO <- ggplot(myframe_6_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDheto/EEDko)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_CKO,  path1=subdir_6_part4, fileName1="4-6-2E-merge-figure-CKO-Line-ratio",  height1=4, width1=7)









############################# TAC 2 samples
sink( file=paste(subdir_6_part4,  "/4-6-3A-1-banding-runLog.txt",      sep = "") )
myNTR_6_banding_3A  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(reduceColumn2_Average_week0[i]+0.0001,  reduceColumn2_Average_banding[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  2.5),  yAxis2=log( vec1 )  , file2=paste(subdir_6_part4,  "/4-6-3A-2-banding-LogLinearModel", sep="") )
  myNTR_6_banding_3A <- c(myNTR_6_banding_3A, NTR_bin)
}
sink()  
length(myNTR_6_banding_3A)
summary(myNTR_6_banding_3A)
HalfLife_6_banding_3A <- log(2)/(myNTR_6_banding_3A+0.0001)
length(HalfLife_6_banding_3A)
summary(HalfLife_6_banding_3A)

write.table(x=myNTR_6_banding_3A,      file = paste(subdir_6_part4,  "/4-6-3A-3-banding-raw-myNTR.txt",  sep = ""),      append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_6_banding_3A,   file = paste(subdir_6_part4,  "/4-6-3A-4-banding-raw-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_6_banding_3A <- data.frame( xAxis = c(1:length(myNTR_6_banding_3A)),      yAxis = sort(myNTR_6_banding_3A) )
FigureTemp_6_banding_3A <- ggplot(myframe_6_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_banding_3A,  path1=subdir_6_part4, fileName1="4-6-3A-1-figure-banding-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_6_banding_3A <- ggplot(myframe_6_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_banding_3A,  path1=subdir_6_part4, fileName1="4-6-3A-1-figure-banding-Line-raw",  height1=4, width1=7)

myNTR_6_banding_3A[myNTR_6_banding_3A<0]  <- 0
myNTR_6_banding_3A[myNTR_6_banding_3A>1]  <- 1
HalfLife_6_banding_3A[HalfLife_6_banding_3A<0]  <- 0
HalfLife_6_banding_3A[HalfLife_6_banding_3A>50]  <- 50
min(myNTR_6_banding_3A)

write.table(x=myNTR_6_banding_3A,      file = paste(subdir_6_part4,  "/4-6-3A-6-noLess0-myNTR_6_banding_3A.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_6_banding_3A,   file = paste(subdir_6_part4,  "/4-6-3A-6-noLess0-HalfLife_6_banding_3A.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_6_banding_3A <- data.frame( xAxis = c(1:length(myNTR_6_banding_3A)),      yAxis = sort(myNTR_6_banding_3A) )
FigureTemp_6_banding_3A <- ggplot(myframe_6_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_banding_3A,  path1=subdir_6_part4, fileName1="4-6-3A-2-figure-banding-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_6_banding_3A <- ggplot(myframe_6_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_banding_3A,  path1=subdir_6_part4, fileName1="4-6-3A-2-figure-banding-Line-noLess0",  height1=4, width1=7)



sink( file=paste(subdir_6_part4,  "/4-6-3B-1-sham-runLog.txt",      sep = "") )
myNTR_6_sham_3B  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(reduceColumn2_Average_week0[i]+0.0001,  reduceColumn2_Average_sham[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  2.5),  yAxis2=log( vec1 )  , file2=paste(subdir_6_part4,  "/4-6-3B-2-sham-LogLinearModel", sep="") )
  myNTR_6_sham_3B <- c(myNTR_6_sham_3B, NTR_bin)
}
sink()  
length(myNTR_6_sham_3B)
summary(myNTR_6_sham_3B)
HalfLife_6_sham_3B <- log(2)/(myNTR_6_sham_3B+0.0001)
length(HalfLife_6_sham_3B)
summary(HalfLife_6_sham_3B)

write.table(x=myNTR_6_sham_3B,      file = paste(subdir_6_part4,  "/4-6-3B-3-sham-raw-NTR.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_6_sham_3B,   file = paste(subdir_6_part4,  "/4-6-3B-4-sham-raw-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_6_sham_3B <- data.frame( xAxis = c(1:length(myNTR_6_sham_3B)),      yAxis = sort(myNTR_6_sham_3B) )
FigureTemp_6_sham_3B <- ggplot(myframe_6_sham_3B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_sham_3B,  path1=subdir_6_part4, fileName1="4-6-3B-1-figure-sham-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_6_sham_3B <- ggplot(myframe_6_sham_3B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_sham_3B,  path1=subdir_6_part4, fileName1="4-6-3B-1-figure-sham-Line-raw",  height1=4, width1=7)

myNTR_6_sham_3B[myNTR_6_sham_3B<0]  <- 0
myNTR_6_sham_3B[myNTR_6_sham_3B>1]  <- 1
HalfLife_6_sham_3B[HalfLife_6_sham_3B<0]   <- 0
HalfLife_6_sham_3B[HalfLife_6_sham_3B>50]  <- 50

write.table(x=myNTR_6_sham_3B,      file = paste(subdir_6_part4,  "/4-6-3B-6-sham-noLess0-NTR.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_6_sham_3B,   file = paste(subdir_6_part4,  "/4-6-3B-6-sham-noLess0-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_6_sham_3B <- data.frame( xAxis = c(1:length(myNTR_6_sham_3B)),      yAxis = sort(myNTR_6_sham_3B) )
FigureTemp_6_sham_3B <- ggplot(myframe_6_sham_3B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_sham_3B,  path1=subdir_6_part4, fileName1="4-6-3B-2-figure-sham-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_6_sham_3B <- ggplot(myframe_6_sham_3B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_sham_3B,  path1=subdir_6_part4, fileName1="4-6-3B-2-figure-sham-Line-noLess0",  height1=4, width1=7)



MyBoxViolinPlot_1(vector2=c(myNTR_6_banding_3A,  myNTR_6_sham_3B),       
                  sampleType2=c( rep("banding", numOfRows1)  ,   rep("sham", numOfRows1)  ), 
                  sampleRank2=c( "banding",  "sham" ),     
                  colours2=c( "banding"="blue",  "sham"="green4"  ), 
                  path2=subdir_6_part4,     fileName2="4-6-3A-merge-TAC-NTR-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                  Ymin2=0,   Ymax2=1.5,    height2=5,   width2=2  )

MyBoxViolinPlot_1(vector2=c(HalfLife_6_banding_3A,  HalfLife_6_sham_3B),       
                  sampleType2=c( rep("banding", numOfRows1)  ,   rep("sham", numOfRows1)  ), 
                  sampleRank2=c( "banding",  "sham" ),     
                  colours2=c( "banding"="blue",  "sham"="green4"  ), 
                  path2=subdir_6_part4,     fileName2="4-6-3B-merge-TAC-HalfLife-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half life (weeks)",   
                  Ymin2=0,   Ymax2=12,    height2=5,   width2=2  )




diff1_6_TAC = myNTR_6_banding_3A - myNTR_6_sham_3B 
diff2_6_TAC = log2( ( (myNTR_6_banding_3A+ 0.001) / (myNTR_6_sham_3B+0.001))  )
diff3_6_TAC = log2( ( (myNTR_6_sham_3B+ 0.001) / (myNTR_6_banding_3A+0.001))  )
length(diff1_6_TAC)
length(diff2_6_TAC)
length(diff3_6_TAC)

myframe_6_TAC <- data.frame( xAxis = c(1:length(diff1_6_TAC)),      yAxis = sort(diff1_6_TAC) )

FigureTemp_6_TAC <- ggplot(myframe_6_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR (banding - sham)") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.1, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -0.1, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_TAC,  path1=subdir_6_part4, fileName1="4-6-3C-merge1-figure-TAC-Line-minus-limitY",  height1=4, width1=7)

FigureTemp_6_TAC <- ggplot(myframe_6_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR (banding - sham)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.1, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -0.1, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_TAC,  path1=subdir_6_part4, fileName1="4-6-3C-merge1-figure-TAC-Line-minus",  height1=4, width1=7)





myframe_6_TAC <- data.frame( xAxis = c(1:length(diff2_6_TAC)),      yAxis = sort(diff2_6_TAC) )

FigureTemp_6_TAC <- ggplot(myframe_6_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(banding/sham)") +  ggtitle(myTitle_g) + ylim(-5, 5 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept =1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_TAC,  path1=subdir_6_part4, fileName1="4-6-3C-merge2-figure-TAC-Line-ratio-limitY",  height1=4, width1=7)

FigureTemp_6_TAC <- ggplot(myframe_6_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(banding/sham)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_TAC,  path1=subdir_6_part4, fileName1="4-6-3C-merge2-figure-TAC-Line-ratio",  height1=4, width1=7)




myframe_6_TAC <- data.frame( xAxis = c(1:length(diff3_6_TAC)),      yAxis = sort(diff3_6_TAC) )

FigureTemp_6_TAC <- ggplot(myframe_6_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(sham/banding)") +  ggtitle(myTitle_g) + ylim(-5, 5 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept =1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_TAC,  path1=subdir_6_part4, fileName1="4-6-3C-merge3-figure-TAC-Line-ratio-limitY",  height1=4, width1=7)

FigureTemp_6_TAC <- ggplot(myframe_6_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(sham/banding)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_6_TAC,  path1=subdir_6_part4, fileName1="4-6-3C-merge3-figure-TAC-Line-ratio",  height1=4, width1=7)






















##########################################################################   
subdir_6A_part4 <- paste(Part4_g,  "/6A-NTR-eachRow-2kb", sep = "")
if( ! file.exists(subdir_6A_part4) ) { dir.create(subdir_6A_part4) }


dim(Average_H3)
dim(Average_week0)
numRows_oneClass <- floor( nrow(Average_H3)/5 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
INDEX_4 <- seq(from = numRows_oneClass*3+1, to = numRows_oneClass*4, by =1 )
INDEX_5 <- seq(from = numRows_oneClass*4+1, to = numRows_oneClass*5, by =1 )
nrow(Average_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) + length(INDEX_4) + length(INDEX_5))

MyBoxViolinPlot_1(vector2=c( myNTR_6_WT_1[INDEX_1], myNTR_6_WT_1[INDEX_2],  myNTR_6_WT_1[INDEX_3],  myNTR_6_WT_1[INDEX_4],  myNTR_6_WT_1[INDEX_5] ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),  rep("Medium", length(INDEX_3)),  
                                 rep("High", length(INDEX_4)),      rep("Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "Lowest",  "Low",   "Medium",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",    "Medium"="green2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_6A_part4,   fileName2= paste("4-6A-1-WT-NTR-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=3.5,   Ymin2=0, Ymax2=0.4)    ## width = 1 + 5*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c( HalfLife_6_WT_1[INDEX_1], HalfLife_6_WT_1[INDEX_2],  HalfLife_6_WT_1[INDEX_3],  HalfLife_6_WT_1[INDEX_4],  HalfLife_6_WT_1[INDEX_5] ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),  rep("Medium", length(INDEX_3)),  
                                 rep("High", length(INDEX_4)),      rep("Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "Lowest",  "Low",   "Medium",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",    "Medium"="green2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_6A_part4,   fileName2= paste("4-6A-2-WT-HalfLife-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=3.5,   Ymin2=0, Ymax2=50)    ## width = 1 + 5*0.5, height=5cm 











MyBoxViolinPlot_1(vector2=c( myNTR_6_EEDheto_2A[INDEX_1], myNTR_6_EEDheto_2A[INDEX_2],  myNTR_6_EEDheto_2A[INDEX_3],  myNTR_6_EEDheto_2A[INDEX_4],  myNTR_6_EEDheto_2A[INDEX_5],
                             myNTR_6_EEDko_2B[INDEX_1],   myNTR_6_EEDko_2B[INDEX_2],    myNTR_6_EEDko_2B[INDEX_3],    myNTR_6_EEDko_2B[INDEX_4],    myNTR_6_EEDko_2B[INDEX_5] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),  rep("EEDheto_Medium", length(INDEX_3)),  rep("EEDheto_High", length(INDEX_4)),  rep("EEDheto_Highest", length(INDEX_5)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),    rep("EEDko_Medium", length(INDEX_3)),    rep("EEDko_High", length(INDEX_4)),    rep("EEDko_Highest", length(INDEX_5)) ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",     "EEDko_Medium",    "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2",    "EEDko_Highest"="purple2" ), 
                  path2=subdir_6A_part4,   fileName2= paste("4-6A-3-WT-NTR-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=0.3)    ## width = 1 + 10*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_6_EEDheto_2A[INDEX_1], HalfLife_6_EEDheto_2A[INDEX_2],  HalfLife_6_EEDheto_2A[INDEX_3],  HalfLife_6_EEDheto_2A[INDEX_4],  HalfLife_6_EEDheto_2A[INDEX_5],
                             HalfLife_6_EEDko_2B[INDEX_1],   HalfLife_6_EEDko_2B[INDEX_2],    HalfLife_6_EEDko_2B[INDEX_3],    HalfLife_6_EEDko_2B[INDEX_4],    HalfLife_6_EEDko_2B[INDEX_5] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),  rep("EEDheto_Medium", length(INDEX_3)),  rep("EEDheto_High", length(INDEX_4)),  rep("EEDheto_Highest", length(INDEX_5)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),    rep("EEDko_Medium", length(INDEX_3)),    rep("EEDko_High", length(INDEX_4)),    rep("EEDko_Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",     "EEDko_Medium",    "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2",    "EEDko_Highest"="purple2"  ), 
                  path2=subdir_6A_part4,   fileName2= paste("4-6A-4-WT-HalfLife-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=40)    ## width = 1 + 10*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( myNTR_6_EEDheto_2A[INDEX_1], myNTR_6_EEDheto_2A[INDEX_2],  myNTR_6_EEDheto_2A[INDEX_3],  myNTR_6_EEDheto_2A[INDEX_4],  myNTR_6_EEDheto_2A[INDEX_5],
                             myNTR_6_EEDko_2B[INDEX_1],   myNTR_6_EEDko_2B[INDEX_2],    myNTR_6_EEDko_2B[INDEX_3],    myNTR_6_EEDko_2B[INDEX_4],    myNTR_6_EEDko_2B[INDEX_5] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),  rep("EEDheto_Medium", length(INDEX_3)),  rep("EEDheto_High", length(INDEX_4)),  rep("EEDheto_Highest", length(INDEX_5)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),    rep("EEDko_Medium", length(INDEX_3)),    rep("EEDko_High", length(INDEX_4)),    rep("EEDko_Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "EEDheto_Lowest",   "EEDko_Lowest",   
                                 "EEDheto_Low",      "EEDko_Low", 
                                 "EEDheto_Medium",   "EEDko_Medium", 
                                 "EEDheto_High",     "EEDko_High",  
                                 "EEDheto_Highest",  "EEDko_Highest" ),     
                  colours2=c( "EEDheto_Lowest"="red2",        "EEDko_Lowest"="red2",   
                              "EEDheto_Low"="cyan2",          "EEDko_Low"="cyan2",  
                              "EEDheto_Medium"="green2",      "EEDko_Medium"="green2", 
                              "EEDheto_High"="blue2",         "EEDko_High"="blue2",  
                              "EEDheto_Highest"="purple2" ,   "EEDko_Highest"="purple2" ), 
                  path2=subdir_6A_part4,   fileName2= paste("4-6A-5-WT-NTR-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0,   Ymax2=0.3)    ## width = 1 + 10*0.5, height=5cm 









MyBoxViolinPlot_1(vector2=c( myNTR_6_banding_3A[INDEX_1], myNTR_6_banding_3A[INDEX_2],  myNTR_6_banding_3A[INDEX_3],  myNTR_6_banding_3A[INDEX_4],  myNTR_6_banding_3A[INDEX_5],
                             myNTR_6_sham_3B[INDEX_1],    myNTR_6_sham_3B[INDEX_2],     myNTR_6_sham_3B[INDEX_3],     myNTR_6_sham_3B[INDEX_4],     myNTR_6_sham_3B[INDEX_5]  ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),  rep("banding_Medium", length(INDEX_3)),  rep("banding_High", length(INDEX_4)),  rep("banding_Highest", length(INDEX_5)), 
                                 rep("sham_Lowest", length(INDEX_1)),     rep("sham_Low", length(INDEX_2)),     rep("sham_Medium", length(INDEX_3)),     rep("sham_High", length(INDEX_4)),     rep("sham_Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",   "banding_Medium",  "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",     "sham_Low",      "sham_Medium",     "sham_High",     "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",      "sham_Low"="cyan2",       "sham_Medium"="green2",     "sham_High"="blue2",     "sham_Highest"="purple2"  ), 
                  path2=subdir_6A_part4,   fileName2= paste("4-6A-6-WT-NTR-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=1)    ## width = 1 + 10*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_6_banding_3A[INDEX_1], HalfLife_6_banding_3A[INDEX_2],  HalfLife_6_banding_3A[INDEX_3],  HalfLife_6_banding_3A[INDEX_4],  HalfLife_6_banding_3A[INDEX_5],
                             HalfLife_6_sham_3B[INDEX_1],    HalfLife_6_sham_3B[INDEX_2],     HalfLife_6_sham_3B[INDEX_3],     HalfLife_6_sham_3B[INDEX_4],     HalfLife_6_sham_3B[INDEX_5] ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),  rep("banding_Medium", length(INDEX_3)),  rep("banding_High", length(INDEX_4)),  rep("banding_Highest", length(INDEX_5)), 
                                 rep("sham_Lowest", length(INDEX_1)),    rep("sham_Low", length(INDEX_2)),    rep("sham_Medium", length(INDEX_3)),    rep("sham_High", length(INDEX_4)),    rep("sham_Highest", length(INDEX_5))
                  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",   "banding_Medium",  "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",    "sham_Low",     "sham_Medium",    "sham_High",    "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",     "sham_Low"="cyan2",      "sham_Medium"="green2",    "sham_High"="blue2",    "sham_Highest"="purple2" ), 
                  path2=subdir_6A_part4,   fileName2= paste("4-6A-7-WT-HalfLife-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=10)    ## width = 1 + 10*0.5, height=5cm 



MyBoxViolinPlot_1(vector2=c( myNTR_6_banding_3A[INDEX_1], myNTR_6_banding_3A[INDEX_2],  myNTR_6_banding_3A[INDEX_3],  myNTR_6_banding_3A[INDEX_4],  myNTR_6_banding_3A[INDEX_5],
                             myNTR_6_sham_3B[INDEX_1],   myNTR_6_sham_3B[INDEX_2],    myNTR_6_sham_3B[INDEX_3],    myNTR_6_sham_3B[INDEX_4],    myNTR_6_sham_3B[INDEX_5] ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),  rep("banding_Medium", length(INDEX_3)),  rep("banding_High", length(INDEX_4)),  rep("banding_Highest", length(INDEX_5)), 
                                 rep("sham_Lowest", length(INDEX_1)),    rep("sham_Low", length(INDEX_2)),    rep("sham_Medium", length(INDEX_3)),    rep("sham_High", length(INDEX_4)),    rep("sham_Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "banding_Lowest",   "sham_Lowest",   
                                 "banding_Low",      "sham_Low", 
                                 "banding_Medium",   "sham_Medium", 
                                 "banding_High",     "sham_High",  
                                 "banding_Highest",  "sham_Highest" ),     
                  colours2=c( "banding_Lowest"="red2",        "sham_Lowest"="red2",   
                              "banding_Low"="cyan2",          "sham_Low"="cyan2",  
                              "banding_Medium"="green2",      "sham_Medium"="green2", 
                              "banding_High"="blue2",         "sham_High"="blue2",  
                              "banding_Highest"="purple2" ,   "sham_Highest"="purple2"  ), 
                  path2=subdir_6A_part4,   fileName2= paste("4-6A-8-WT-NTR-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0,   Ymax2=1)    ## width = 1 + 10*0.5, height=5cm 































##########################################################################   
subdir_6B_part4 <- paste(Part4_g,  "/6B-NTR-eachRow-2kb", sep = "")
if( ! file.exists(subdir_6B_part4) ) { dir.create(subdir_6B_part4) }


dim(Average_H3)
dim(Average_week0)
numRows_oneClass <- floor( nrow(Average_H3)/4 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
INDEX_4 <- seq(from = numRows_oneClass*3+1, to = numRows_oneClass*4, by =1 )
nrow(Average_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) + length(INDEX_4) )

MyBoxViolinPlot_1(vector2=c( myNTR_6_WT_1[INDEX_1], myNTR_6_WT_1[INDEX_2],  myNTR_6_WT_1[INDEX_3],  myNTR_6_WT_1[INDEX_4]  ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),  rep("High", length(INDEX_3)),      rep("Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "Lowest",  "Low",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_6B_part4,   fileName2= paste("4-6B-1-WT-NTR-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=3.0,   Ymin2=0, Ymax2=0.4)    ## width = 1 + 4*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c( HalfLife_6_WT_1[INDEX_1], HalfLife_6_WT_1[INDEX_2],  HalfLife_6_WT_1[INDEX_3],  HalfLife_6_WT_1[INDEX_4]  ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),  rep("High", length(INDEX_3)),      rep("Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "Lowest",  "Low",    "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",   "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_6B_part4,   fileName2= paste("4-6B-2-WT-HalfLife-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=3.0,   Ymin2=0, Ymax2=50)    ## width = 1 + 4*0.5, height=5cm 











MyBoxViolinPlot_1(vector2=c( myNTR_6_EEDheto_2A[INDEX_1], myNTR_6_EEDheto_2A[INDEX_2],  myNTR_6_EEDheto_2A[INDEX_3],  myNTR_6_EEDheto_2A[INDEX_4],  
                             myNTR_6_EEDko_2B[INDEX_1],   myNTR_6_EEDko_2B[INDEX_2],    myNTR_6_EEDko_2B[INDEX_3],    myNTR_6_EEDko_2B[INDEX_4] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),       rep("EEDheto_High", length(INDEX_3)),  rep("EEDheto_Highest", length(INDEX_4)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),         rep("EEDko_High", length(INDEX_3)),    rep("EEDko_Highest", length(INDEX_4)) ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",        "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",          "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",         "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",           "EEDko_High"="blue2",    "EEDko_Highest"="purple2" ), 
                  path2=subdir_6B_part4,   fileName2= paste("4-6B-3-WT-NTR-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=0.3)    ## width = 1 + 8*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_6_EEDheto_2A[INDEX_1], HalfLife_6_EEDheto_2A[INDEX_2],  HalfLife_6_EEDheto_2A[INDEX_3],  HalfLife_6_EEDheto_2A[INDEX_4],   
                             HalfLife_6_EEDko_2B[INDEX_1],   HalfLife_6_EEDko_2B[INDEX_2],    HalfLife_6_EEDko_2B[INDEX_3],    HalfLife_6_EEDko_2B[INDEX_4]  ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),       rep("EEDheto_High", length(INDEX_3)),  rep("EEDheto_Highest", length(INDEX_4)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),         rep("EEDko_High", length(INDEX_3)),    rep("EEDko_Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",      "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",        "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",     "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",       "EEDko_High"="blue2",    "EEDko_Highest"="purple2"  ), 
                  path2=subdir_6B_part4,   fileName2= paste("4-6B-4-WT-HalfLife-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=40)    ## width = 1 + 8*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( myNTR_6_EEDheto_2A[INDEX_1], myNTR_6_EEDheto_2A[INDEX_2],  myNTR_6_EEDheto_2A[INDEX_3],  myNTR_6_EEDheto_2A[INDEX_4],  
                             myNTR_6_EEDko_2B[INDEX_1],   myNTR_6_EEDko_2B[INDEX_2],    myNTR_6_EEDko_2B[INDEX_3],    myNTR_6_EEDko_2B[INDEX_4]  ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),    rep("EEDheto_High", length(INDEX_3)),  rep("EEDheto_Highest", length(INDEX_4)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),      rep("EEDko_High", length(INDEX_3)),    rep("EEDko_Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "EEDheto_Lowest",   "EEDko_Lowest",   
                                 "EEDheto_Low",      "EEDko_Low", 
                                 "EEDheto_High",     "EEDko_High",  
                                 "EEDheto_Highest",  "EEDko_Highest" ),     
                  colours2=c( "EEDheto_Lowest"="red2",        "EEDko_Lowest"="red2",   
                              "EEDheto_Low"="cyan2",          "EEDko_Low"="cyan2",  
                              "EEDheto_High"="blue2",         "EEDko_High"="blue2",  
                              "EEDheto_Highest"="purple2" ,   "EEDko_Highest"="purple2" ), 
                  path2=subdir_6B_part4,   fileName2= paste("4-6B-5-WT-NTR-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=5,   Ymin2=0,   Ymax2=0.3)    ## width = 1 + 8*0.5, height=5cm 









MyBoxViolinPlot_1(vector2=c( myNTR_6_banding_3A[INDEX_1], myNTR_6_banding_3A[INDEX_2],  myNTR_6_banding_3A[INDEX_3],  myNTR_6_banding_3A[INDEX_4],  
                             myNTR_6_sham_3B[INDEX_1],    myNTR_6_sham_3B[INDEX_2],     myNTR_6_sham_3B[INDEX_3],     myNTR_6_sham_3B[INDEX_4]   ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),   rep("banding_High", length(INDEX_3)),  rep("banding_Highest", length(INDEX_4)), 
                                 rep("sham_Lowest", length(INDEX_1)),     rep("sham_Low", length(INDEX_2)),      rep("sham_High", length(INDEX_3)),     rep("sham_Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",       "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",     "sham_Low",          "sham_High",     "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",        "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",      "sham_Low"="cyan2",            "sham_High"="blue2",     "sham_Highest"="purple2"  ), 
                  path2=subdir_6B_part4,   fileName2= paste("4-6B-6-WT-NTR-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=1)    ## width = 1 + 8*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_6_banding_3A[INDEX_1], HalfLife_6_banding_3A[INDEX_2],  HalfLife_6_banding_3A[INDEX_3],  HalfLife_6_banding_3A[INDEX_4],  
                             HalfLife_6_sham_3B[INDEX_1],    HalfLife_6_sham_3B[INDEX_2],     HalfLife_6_sham_3B[INDEX_3],     HalfLife_6_sham_3B[INDEX_4]  ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),     rep("banding_High", length(INDEX_3)),  rep("banding_Highest", length(INDEX_4)), 
                                 rep("sham_Lowest", length(INDEX_1)),    rep("sham_Low", length(INDEX_2)),         rep("sham_High", length(INDEX_3)),    rep("sham_Highest", length(INDEX_4))
                  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",     "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",    "sham_Low",          "sham_High",    "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",      "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",     "sham_Low"="cyan2",           "sham_High"="blue2",    "sham_Highest"="purple2" ), 
                  path2=subdir_6B_part4,   fileName2= paste("4-6B-7-WT-HalfLife-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=10)    ## width = 1 + 8*0.5, height=5cm 



MyBoxViolinPlot_1(vector2=c( myNTR_6_banding_3A[INDEX_1], myNTR_6_banding_3A[INDEX_2],  myNTR_6_banding_3A[INDEX_3],  myNTR_6_banding_3A[INDEX_4],   
                             myNTR_6_sham_3B[INDEX_1],   myNTR_6_sham_3B[INDEX_2],    myNTR_6_sham_3B[INDEX_3],    myNTR_6_sham_3B[INDEX_4]  ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),    rep("banding_High", length(INDEX_3)),  rep("banding_Highest", length(INDEX_4)), 
                                 rep("sham_Lowest", length(INDEX_1)),    rep("sham_Low", length(INDEX_2)),        rep("sham_High", length(INDEX_3)),    rep("sham_Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "banding_Lowest",   "sham_Lowest",   
                                 "banding_Low",      "sham_Low", 
                                 "banding_High",     "sham_High",  
                                 "banding_Highest",  "sham_Highest" ),     
                  colours2=c( "banding_Lowest"="red2",        "sham_Lowest"="red2",   
                              "banding_Low"="cyan2",          "sham_Low"="cyan2",  
                              "banding_High"="blue2",         "sham_High"="blue2",  
                              "banding_Highest"="purple2" ,   "sham_Highest"="purple2"  ), 
                  path2=subdir_6B_part4,   fileName2= paste("4-6B-8-WT-NTR-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=5,   Ymin2=0,   Ymax2=1)    ## width = 1 + 8*0.5, height=5cm 









































##########################################################################   
subdir_6C_part4 <- paste(Part4_g,  "/6C-NTR-eachRow-2kb", sep = "")
if( ! file.exists(subdir_6C_part4) ) { dir.create(subdir_6C_part4) }


dim(Average_H3)
dim(Average_week0)
numRows_oneClass <- floor( nrow(Average_H3)/3 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
nrow(Average_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) )

MyBoxViolinPlot_1(vector2=c( myNTR_6_WT_1[INDEX_1], myNTR_6_WT_1[INDEX_2],  myNTR_6_WT_1[INDEX_3]  ),   
                  sampleType2=c(  rep("Low", length(INDEX_1)),    rep("Medium", length(INDEX_2)),    rep("High", length(INDEX_3))   ), 
                  sampleRank2=c(  "Low",   "Medium",  "High"  ),     
                  colours2=c(  "Low"="cyan2",    "Medium"="green2",  "High"="blue2"   ), 
                  path2=subdir_6C_part4,   fileName2= paste("4-6C-1-WT-NTR-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=2.5,   Ymin2=0, Ymax2=0.4)    ## width = 1 + 3*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c( HalfLife_6_WT_1[INDEX_1], HalfLife_6_WT_1[INDEX_2],  HalfLife_6_WT_1[INDEX_3]  ),   
                  sampleType2=c(  rep("Low", length(INDEX_1)),  rep("Medium", length(INDEX_2)),   rep("High", length(INDEX_3))   ), 
                  sampleRank2=c(  "Low",   "Medium",  "High"  ),     
                  colours2=c( "Low"="cyan2",    "Medium"="green2",  "High"="blue2"   ), 
                  path2=subdir_6C_part4,   fileName2= paste("4-6C-2-WT-HalfLife-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=2.5,   Ymin2=0, Ymax2=50)    ## width = 1 + 3*0.5, height=5cm 











MyBoxViolinPlot_1(vector2=c( myNTR_6_EEDheto_2A[INDEX_1], myNTR_6_EEDheto_2A[INDEX_2],  myNTR_6_EEDheto_2A[INDEX_3],  
                             myNTR_6_EEDko_2B[INDEX_1],   myNTR_6_EEDko_2B[INDEX_2],    myNTR_6_EEDko_2B[INDEX_3]  ),   
                  sampleType2=c(   rep("EEDheto_Low", length(INDEX_1)),  rep("EEDheto_Medium", length(INDEX_2)),  rep("EEDheto_High", length(INDEX_3)),    
                                   rep("EEDko_Low", length(INDEX_1)),    rep("EEDko_Medium", length(INDEX_2)),    rep("EEDko_High", length(INDEX_3))  ), 
                  sampleRank2=c(    "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",   
                                    "EEDko_Low",     "EEDko_Medium",    "EEDko_High"   ),     
                  colours2=c(    "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",   
                                 "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2"  ), 
                  path2=subdir_6C_part4,   fileName2= paste("4-6C-3-WT-NTR-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=4,   Ymin2=0, Ymax2=0.3)    ## width = 1 + 6*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_6_EEDheto_2A[INDEX_1], HalfLife_6_EEDheto_2A[INDEX_2],  HalfLife_6_EEDheto_2A[INDEX_3],   
                             HalfLife_6_EEDko_2B[INDEX_1],   HalfLife_6_EEDko_2B[INDEX_2],    HalfLife_6_EEDko_2B[INDEX_3]  ),   
                  sampleType2=c(    rep("EEDheto_Low", length(INDEX_1)),  rep("EEDheto_Medium", length(INDEX_2)),  rep("EEDheto_High", length(INDEX_3)),  
                                    rep("EEDko_Low", length(INDEX_1)),    rep("EEDko_Medium", length(INDEX_2)),    rep("EEDko_High", length(INDEX_3))  ), 
                  sampleRank2=c(  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  
                                  "EEDko_Low",     "EEDko_Medium",    "EEDko_High"   ),     
                  colours2=c(  "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",   
                               "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2"   ), 
                  path2=subdir_6C_part4,   fileName2= paste("4-6C-4-WT-HalfLife-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=4,   Ymin2=0, Ymax2=40)    ## width = 1 + 6*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( myNTR_6_EEDheto_2A[INDEX_1], myNTR_6_EEDheto_2A[INDEX_2],  myNTR_6_EEDheto_2A[INDEX_3],   
                             myNTR_6_EEDko_2B[INDEX_1],   myNTR_6_EEDko_2B[INDEX_2],    myNTR_6_EEDko_2B[INDEX_3]  ),   
                  sampleType2=c(    rep("EEDheto_Low", length(INDEX_1)),  rep("EEDheto_Medium", length(INDEX_2)),  rep("EEDheto_High", length(INDEX_3)),   
                                    rep("EEDko_Low", length(INDEX_1)),    rep("EEDko_Medium", length(INDEX_2)),    rep("EEDko_High", length(INDEX_3))  ), 
                  sampleRank2=c( "EEDheto_Low",      "EEDko_Low", 
                                 "EEDheto_Medium",   "EEDko_Medium", 
                                 "EEDheto_High",     "EEDko_High" ),     
                  colours2=c( "EEDheto_Low"="cyan2",          "EEDko_Low"="cyan2",  
                              "EEDheto_Medium"="green2",      "EEDko_Medium"="green2", 
                              "EEDheto_High"="blue2",         "EEDko_High"="blue2" ), 
                  path2=subdir_6C_part4,   fileName2= paste("4-6C-5-WT-NTR-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=4,   Ymin2=0,   Ymax2=0.3)    ## width = 1 + 6*0.5, height=5cm 









MyBoxViolinPlot_1(vector2=c( myNTR_6_banding_3A[INDEX_1], myNTR_6_banding_3A[INDEX_2],  myNTR_6_banding_3A[INDEX_3],   
                             myNTR_6_sham_3B[INDEX_1],    myNTR_6_sham_3B[INDEX_2],     myNTR_6_sham_3B[INDEX_3]    ),   
                  sampleType2=c(  rep("banding_Low", length(INDEX_1)),  rep("banding_Medium", length(INDEX_2)),  rep("banding_High", length(INDEX_3)),   
                                  rep("sham_Low", length(INDEX_1)),     rep("sham_Medium", length(INDEX_2)),     rep("sham_High", length(INDEX_3))   ), 
                  sampleRank2=c(  "banding_Low",   "banding_Medium",  "banding_High",  
                                  "sham_Low",      "sham_Medium",     "sham_High"   ),     
                  colours2=c(  "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  
                               "sham_Low"="cyan2",       "sham_Medium"="green2",     "sham_High"="blue2"   ), 
                  path2=subdir_6C_part4,   fileName2= paste("4-6C-6-WT-NTR-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=4,   Ymin2=0, Ymax2=1)    ## width = 1 + 6*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_6_banding_3A[INDEX_1], HalfLife_6_banding_3A[INDEX_2],  HalfLife_6_banding_3A[INDEX_3],  
                             HalfLife_6_sham_3B[INDEX_1],    HalfLife_6_sham_3B[INDEX_2],     HalfLife_6_sham_3B[INDEX_3]  ),   
                  sampleType2=c(   rep("banding_Low", length(INDEX_1)),  rep("banding_Medium", length(INDEX_2)),  rep("banding_High", length(INDEX_3)),   
                                   rep("sham_Low", length(INDEX_1)),     rep("sham_Medium", length(INDEX_2)),     rep("sham_High", length(INDEX_3)) ), 
                  sampleRank2=c(   "banding_Low",   "banding_Medium",  "banding_High",   
                                   "sham_Low",     "sham_Medium",    "sham_High"    ),     
                  colours2=c(   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",   
                                "sham_Low"="cyan2",      "sham_Medium"="green2",    "sham_High"="blue2"  ), 
                  path2=subdir_6C_part4,   fileName2= paste("4-6C-7-WT-HalfLife-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=4,   Ymin2=0, Ymax2=10)    ## width = 1 + 6*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( myNTR_6_banding_3A[INDEX_1], myNTR_6_banding_3A[INDEX_2],  myNTR_6_banding_3A[INDEX_3],   
                             myNTR_6_sham_3B[INDEX_1],   myNTR_6_sham_3B[INDEX_2],    myNTR_6_sham_3B[INDEX_3]  ),   
                  sampleType2=c(   rep("banding_Low", length(INDEX_1)),  rep("banding_Medium", length(INDEX_2)),  rep("banding_High", length(INDEX_3)),   
                                   rep("sham_Low", length(INDEX_1)),     rep("sham_Medium", length(INDEX_2)),    rep("sham_High", length(INDEX_3))  ), 
                  sampleRank2=c( "banding_Low",      "sham_Low", 
                                 "banding_Medium",   "sham_Medium", 
                                 "banding_High",     "sham_High"  ),     
                  colours2=c( "banding_Low"="cyan2",          "sham_Low"="cyan2",  
                              "banding_Medium"="green2",      "sham_Medium"="green2", 
                              "banding_High"="blue2",         "sham_High"="blue2"   ), 
                  path2=subdir_6C_part4,   fileName2= paste("4-6C-8-WT-NTR-BoxViolin",   "2kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=4,   Ymin2=0,   Ymax2=1)    ## width = 1 + 6*0.5, height=5cm 































#########################################################################################################  average some columns for each row
subdir_7_part4 <- paste(Part4_g,  "/7-NTR-eachRow-1000bp", sep = "")
if( ! file.exists(subdir_7_part4) ) { dir.create(subdir_7_part4) }


############################# WT 6 samples
sink( file=paste(subdir_7_part4,  "/4-7-1A-WT-runLog.txt",      sep = "") )
myNTR_7_WT_1  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(reduceColumn3_Average_week0[i]+0.0006, reduceColumn3_Average_week1[i]+0.0005, reduceColumn3_Average_week2[i]+0.0004, reduceColumn3_Average_week4[i]+0.0003, reduceColumn3_Average_week6[i]+0.0002, reduceColumn3_Average_week8[i]+0.0001)                        
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log(vec1)  , file2=paste(subdir_7_part4,  "/4-7-1B-WT-LogLinearModel", sep="") )
  myNTR_7_WT_1 <- c(myNTR_7_WT_1, NTR_bin)
}
sink()  

length(myNTR_7_WT_1)
summary(myNTR_7_WT_1)
HalfLife_7_WT_1 <- log(2)/(myNTR_7_WT_1+0.0001)
length(HalfLife_7_WT_1)
summary(HalfLife_7_WT_1)

write.table(x=myNTR_7_WT_1,      file = paste(subdir_7_part4,  "/4-7-1C-raw-NTR.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_7_WT_1,   file = paste(subdir_7_part4,  "/4-7-1D-raw-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_7_WT_1 <- data.frame( xAxis = c(1:length(myNTR_7_WT_1)),      yAxis = sort(myNTR_7_WT_1) )
FigureTemp_7_WT_1 <- ggplot(myframe_7_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_WT_1,  path1=subdir_7_part4, fileName1="4-7-1A-figure-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_7_WT_1 <- ggplot(myframe_7_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) +  
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_WT_1,  path1=subdir_7_part4, fileName1="4-7-1A-figure-Line-raw",  height1=4, width1=7)

myNTR_7_WT_1[myNTR_7_WT_1<0]  <- 0
myNTR_7_WT_1[myNTR_7_WT_1>1]  <- 1
HalfLife_7_WT_1[HalfLife_7_WT_1<0]  <- 0
HalfLife_7_WT_1[HalfLife_7_WT_1>50] <- 50

write.table(x=thisRowNames1,     file = paste(subdir_7_part4,   "/4-7-1E-RowNames.txt",  sep = ""),         append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=myNTR_7_WT_1,      file = paste(subdir_7_part4,  "/4-7-1F-noLess0-NTR.txt",    sep = ""),     append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_7_WT_1,   file = paste(subdir_7_part4,  "/4-7-1G-noLess0-HalfLife.txt", sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_7_WT_1 <- data.frame( xAxis = c(1:length(myNTR_7_WT_1)),      yAxis = sort(myNTR_7_WT_1) )
FigureTemp_7_WT_1 <- ggplot(myframe_7_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_WT_1,  path1=subdir_7_part4, fileName1="4-7-1B-figure-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_7_WT_1 <- ggplot(myframe_7_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) +  
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_WT_1,  path1=subdir_7_part4, fileName1="4-7-1B-figure-Line-noLess0",  height1=4, width1=7)

MyBoxViolinPlot_1(vector2=myNTR_7_WT_1,      sampleType2=c( rep("NTR", numOfRows1)  ), 
                  sampleRank2=c( "NTR" ),     colours2=c( "NTR"="red" ), 
                  path2=subdir_7_part4,     fileName2="4-7-1C-WT-NTR-average-curve",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="NTR",   
                  Ymin2=0,   Ymax2=0.4,    height2=5,   width2=1.5  )

MyBoxViolinPlot_1(vector2=HalfLife_7_WT_1,      sampleType2=c( rep("HalfLife", numOfRows1)  ), 
                  sampleRank2=c( "HalfLife" ),     colours2=c( "HalfLife"="red" ), 
                  path2=subdir_7_part4,     fileName2="4-7-1D-WT-HalfLife-average-curve",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="Half life (weeks)",   
                  Ymin2=0,   Ymax2=20,    height2=5,   width2=1.5  )










############################# CKO 4 samples
sink( file=paste(subdir_7_part4,  "/4-7-2A-1-EEDheto-runLog.txt",      sep = "") )
myNTR_7_EEDheto_2A  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(reduceColumn3_Average_week0_EEDheto[i]+0.0001,  reduceColumn3_Average_week4_EEDheto[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_7_part4,  "/4-7-2A-2-EEDheto-LogLinearModel", sep="") )
  myNTR_7_EEDheto_2A <- c(myNTR_7_EEDheto_2A, NTR_bin)
}
sink()  

length(myNTR_7_EEDheto_2A)
summary(myNTR_7_EEDheto_2A)
HalfLife_7_EEDheto_2A <- log(2)/(myNTR_7_EEDheto_2A+0.0001)
length(HalfLife_7_EEDheto_2A)
summary(HalfLife_7_EEDheto_2A)

write.table(x=myNTR_7_EEDheto_2A,      file = paste(subdir_7_part4,  "/4-7-2A-3-raw-NTR_EEDheto.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_7_EEDheto_2A,   file = paste(subdir_7_part4,  "/4-7-2A-4-raw-HalfLife_EEDheto.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_7_EEDheto_2A <- data.frame( xAxis = c(1:length(myNTR_7_EEDheto_2A)),      yAxis = sort(myNTR_7_EEDheto_2A) )
FigureTemp_7_EEDheto_2A <- ggplot(myframe_7_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_EEDheto_2A,  path1=subdir_7_part4, fileName1="4-7-2A-1-figure-EEDheto-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_7_EEDheto_2A <- ggplot(myframe_7_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) +
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_EEDheto_2A,  path1=subdir_7_part4, fileName1="4-7-2A-1-figure-EEDheto-Line-raw",  height1=4, width1=7)

myNTR_7_EEDheto_2A[myNTR_7_EEDheto_2A<0]  <- 0
myNTR_7_EEDheto_2A[myNTR_7_EEDheto_2A>1]  <- 1
HalfLife_7_EEDheto_2A[HalfLife_7_EEDheto_2A<0]   <- 0
HalfLife_7_EEDheto_2A[HalfLife_7_EEDheto_2A>50]  <- 50

write.table(x=myNTR_7_EEDheto_2A,      file = paste(subdir_7_part4,  "/4-7-2A-7-noLess0-myNTR_EEDheto.txt",  sep = ""),      append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_7_EEDheto_2A,   file = paste(subdir_7_part4,  "/4-7-2A-7-noLess0-HalfLife_EEDheto.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_7_EEDheto_2A <- data.frame( xAxis = c(1:length(myNTR_7_EEDheto_2A)),      yAxis = sort(myNTR_7_EEDheto_2A) )
FigureTemp_7_EEDheto_2A <- ggplot(myframe_7_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_EEDheto_2A,  path1=subdir_7_part4, fileName1="4-7-2A-2-figure-EEDheto-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_7_EEDheto_2A <- ggplot(myframe_7_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_EEDheto_2A,  path1=subdir_7_part4, fileName1="4-7-2A-2-figure-EEDheto-Line-noLess0",  height1=4, width1=7)





sink( file=paste(subdir_7_part4,  "/4-7-2B-1-EEDko-runLog.txt",      sep = "") )
myNTR_7_EEDko_2B  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(reduceColumn3_Average_week0_EEDko[i]+0.0001,  reduceColumn3_Average_week4_EEDko[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_7_part4,  "/4-7-2B-2-EEDko-LogLinearModel", sep="") )
  myNTR_7_EEDko_2B <- c(myNTR_7_EEDko_2B, NTR_bin)
}
sink()  

length(myNTR_7_EEDko_2B)
summary(myNTR_7_EEDko_2B)
HalfLife_7_EEDko_2B <- log(2)/(myNTR_7_EEDko_2B+0.0001)
length(HalfLife_7_EEDko_2B)
summary(HalfLife_7_EEDko_2B)

write.table(x=myNTR_7_EEDko_2B,      file = paste(subdir_7_part4,  "/4-7-2B-3-raw-NTR_EEDko.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_7_EEDko_2B,   file = paste(subdir_7_part4,  "/4-7-2B-4-raw-HalfLife_EEDko.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_7_EEDko_2B <- data.frame( xAxis = c(1:length(myNTR_7_EEDko_2B)),      yAxis = sort(myNTR_7_EEDko_2B) )
FigureTemp_7_EEDko_2B <- ggplot(myframe_7_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_EEDko_2B,  path1=subdir_7_part4, fileName1="4-7-2B-1-figure-EEDko-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_7_EEDko_2B <- ggplot(myframe_7_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_EEDko_2B,  path1=subdir_7_part4, fileName1="4-7-2B-1-figure-EEDko-Line-raw",  height1=4, width1=7)

myNTR_7_EEDko_2B[myNTR_7_EEDko_2B<0]  <- 0
myNTR_7_EEDko_2B[myNTR_7_EEDko_2B>1]  <- 1
HalfLife_7_EEDko_2B[HalfLife_7_EEDko_2B<0]   <- 0
HalfLife_7_EEDko_2B[HalfLife_7_EEDko_2B>50]  <- 50

write.table(x=myNTR_7_EEDko_2B,      file = paste(subdir_7_part4,  "/4-7-2B-7-noLess0-myNTR_EEDko.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_7_EEDko_2B,   file = paste(subdir_7_part4,  "/4-7-2B-7-noLess0-HalfLife_EEDko.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_7_EEDko_2B <- data.frame( xAxis = c(1:length(myNTR_7_EEDko_2B)),      yAxis = sort(myNTR_7_EEDko_2B) )
FigureTemp_7_EEDko_2B <- ggplot(myframe_7_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_EEDko_2B,  path1=subdir_7_part4, fileName1="4-7-2B-2-figure-EEDko-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_7_EEDko_2B <- ggplot(myframe_7_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_EEDko_2B,  path1=subdir_7_part4, fileName1="4-7-2B-2-figure-EEDko-Line-noLess0",  height1=4, width1=7)



MyBoxViolinPlot_1(vector2=c(myNTR_7_EEDheto_2A,  myNTR_7_EEDko_2B),       
                  sampleType2=c( rep("EEDheto", numOfRows1)  ,   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "EEDheto",  "EEDko" ),     
                  colours2=c( "EEDheto"="red",  "EEDko"="blue"  ), 
                  path2=subdir_7_part4,     fileName2="4-7-2A-merge-CKO-NTR-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                  Ymin2=0,   Ymax2=0.3,    height2=5,   width2=2  )

MyBoxViolinPlot_1(vector2=c(HalfLife_7_EEDheto_2A,  HalfLife_7_EEDko_2B),       
                  sampleType2=c( rep("EEDheto", numOfRows1)  ,   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "EEDheto",  "EEDko" ),     
                  colours2=c( "EEDheto"="red",  "EEDko"="blue"  ), 
                  path2=subdir_7_part4,     fileName2="4-7-2B-merge-CKO-HalfLife-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half life (weeks)",   
                  Ymin2=0,   Ymax2=20,    height2=5,   width2=2  )



diff1_7_CKO = myNTR_7_EEDko_2B - myNTR_7_EEDheto_2A
diff2_7_CKO = log2( ( (myNTR_7_EEDko_2B+0.001) / (myNTR_7_EEDheto_2A+0.001))  )
diff3_7_CKO = log2( ( (myNTR_7_EEDheto_2A+0.001) / (myNTR_7_EEDko_2B+0.001))  )   ## 2^(-10) = 0.001
length(diff1_7_CKO)
length(diff2_7_CKO)
length(diff3_7_CKO)
summary(diff1_7_CKO)
summary(diff2_7_CKO)
summary(diff3_7_CKO)


myframe_7_CKO <- data.frame( xAxis = c(1:length(diff1_7_CKO)),      yAxis = sort(diff1_7_CKO) )

FigureTemp_7_CKO <- ggplot(myframe_7_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR (EEDko - EEDheto)") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.1, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -0.1, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_CKO,  path1=subdir_7_part4, fileName1="4-7-2C-merge-figure-CKO-Line-minus-limitY",  height1=4, width1=7)

FigureTemp_7_CKO <- ggplot(myframe_7_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR (EEDko - EEDheto)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.1, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -0.1, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_CKO,  path1=subdir_7_part4, fileName1="4-7-2C-merge-figure-CKO-Line-minus",  height1=4, width1=7)



myframe_7_CKO <- data.frame( xAxis = c(1:length(diff2_7_CKO)),      yAxis = sort(diff2_7_CKO) )

FigureTemp_7_CKO <- ggplot(myframe_7_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDko/EEDheto)") +  ggtitle(myTitle_g) + ylim(-10, 10 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_CKO,  path1=subdir_7_part4, fileName1="4-7-2D-merge-figure-CKO-Line-ratio-limitY",  height1=4, width1=7)

FigureTemp_7_CKO <- ggplot(myframe_7_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDko/EEDheto)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_CKO,  path1=subdir_7_part4, fileName1="4-7-2D-merge-figure-CKO-Line-ratio",  height1=4, width1=7)





myframe_7_CKO <- data.frame( xAxis = c(1:length(diff3_7_CKO)),      yAxis = sort(diff3_7_CKO) )

FigureTemp_7_CKO <- ggplot(myframe_7_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDheto/EEDko)") +  ggtitle(myTitle_g) + ylim(-10, 10 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_CKO,  path1=subdir_7_part4, fileName1="4-7-2E-merge-figure-CKO-Line-ratio-limitY",  height1=4, width1=7)

FigureTemp_7_CKO <- ggplot(myframe_7_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDheto/EEDko)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_CKO,  path1=subdir_7_part4, fileName1="4-7-2E-merge-figure-CKO-Line-ratio",  height1=4, width1=7)









############################# TAC 2 samples
sink( file=paste(subdir_7_part4,  "/4-7-3A-1-banding-runLog.txt",      sep = "") )
myNTR_7_banding_3A  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(reduceColumn3_Average_week0[i]+0.0001,  reduceColumn3_Average_banding[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  2.5),  yAxis2=log( vec1 )  , file2=paste(subdir_7_part4,  "/4-7-3A-2-banding-LogLinearModel", sep="") )
  myNTR_7_banding_3A <- c(myNTR_7_banding_3A, NTR_bin)
}
sink()  
length(myNTR_7_banding_3A)
summary(myNTR_7_banding_3A)
HalfLife_7_banding_3A <- log(2)/(myNTR_7_banding_3A+0.0001)
length(HalfLife_7_banding_3A)
summary(HalfLife_7_banding_3A)

write.table(x=myNTR_7_banding_3A,      file = paste(subdir_7_part4,  "/4-7-3A-3-banding-raw-myNTR.txt",  sep = ""),      append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_7_banding_3A,   file = paste(subdir_7_part4,  "/4-7-3A-4-banding-raw-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_7_banding_3A <- data.frame( xAxis = c(1:length(myNTR_7_banding_3A)),      yAxis = sort(myNTR_7_banding_3A) )
FigureTemp_7_banding_3A <- ggplot(myframe_7_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_banding_3A,  path1=subdir_7_part4, fileName1="4-7-3A-1-figure-banding-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_7_banding_3A <- ggplot(myframe_7_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_banding_3A,  path1=subdir_7_part4, fileName1="4-7-3A-1-figure-banding-Line-raw",  height1=4, width1=7)

myNTR_7_banding_3A[myNTR_7_banding_3A<0]  <- 0
myNTR_7_banding_3A[myNTR_7_banding_3A>10]  <- 10
HalfLife_7_banding_3A[HalfLife_7_banding_3A<0]  <- 0
HalfLife_7_banding_3A[HalfLife_7_banding_3A>50]  <- 50
min(myNTR_7_banding_3A)

write.table(x=myNTR_7_banding_3A,      file = paste(subdir_7_part4,  "/4-7-3A-7-noLess0-myNTR_7_banding_3A.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_7_banding_3A,   file = paste(subdir_7_part4,  "/4-7-3A-7-noLess0-HalfLife_7_banding_3A.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_7_banding_3A <- data.frame( xAxis = c(1:length(myNTR_7_banding_3A)),      yAxis = sort(myNTR_7_banding_3A) )
FigureTemp_7_banding_3A <- ggplot(myframe_7_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_banding_3A,  path1=subdir_7_part4, fileName1="4-7-3A-2-figure-banding-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_7_banding_3A <- ggplot(myframe_7_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_banding_3A,  path1=subdir_7_part4, fileName1="4-7-3A-2-figure-banding-Line-noLess0",  height1=4, width1=7)



sink( file=paste(subdir_7_part4,  "/4-7-3B-1-sham-runLog.txt",      sep = "") )
myNTR_7_sham_3B  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(reduceColumn3_Average_week0[i]+0.0001,  reduceColumn3_Average_sham[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  2.5),  yAxis2=log( vec1 )  , file2=paste(subdir_7_part4,  "/4-7-3B-2-sham-LogLinearModel", sep="") )
  myNTR_7_sham_3B <- c(myNTR_7_sham_3B, NTR_bin)
}
sink()  
length(myNTR_7_sham_3B)
summary(myNTR_7_sham_3B)
HalfLife_7_sham_3B <- log(2)/(myNTR_7_sham_3B+0.0001)
length(HalfLife_7_sham_3B)
summary(HalfLife_7_sham_3B)

write.table(x=myNTR_7_sham_3B,      file = paste(subdir_7_part4,  "/4-7-3B-3-sham-raw-NTR.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_7_sham_3B,   file = paste(subdir_7_part4,  "/4-7-3B-4-sham-raw-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_7_sham_3B <- data.frame( xAxis = c(1:length(myNTR_7_sham_3B)),      yAxis = sort(myNTR_7_sham_3B) )
FigureTemp_7_sham_3B <- ggplot(myframe_7_sham_3B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_sham_3B,  path1=subdir_7_part4, fileName1="4-7-3B-1-figure-sham-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_7_sham_3B <- ggplot(myframe_7_sham_3B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_sham_3B,  path1=subdir_7_part4, fileName1="4-7-3B-1-figure-sham-Line-raw",  height1=4, width1=7)

myNTR_7_sham_3B[myNTR_7_sham_3B<0]  <- 0
myNTR_7_sham_3B[myNTR_7_sham_3B>10]  <- 10
HalfLife_7_sham_3B[HalfLife_7_sham_3B<0]   <- 0
HalfLife_7_sham_3B[HalfLife_7_sham_3B>50]  <- 50

write.table(x=myNTR_7_sham_3B,      file = paste(subdir_7_part4,  "/4-7-3B-7-sham-noLess0-NTR.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_7_sham_3B,   file = paste(subdir_7_part4,  "/4-7-3B-7-sham-noLess0-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_7_sham_3B <- data.frame( xAxis = c(1:length(myNTR_7_sham_3B)),      yAxis = sort(myNTR_7_sham_3B) )
FigureTemp_7_sham_3B <- ggplot(myframe_7_sham_3B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_sham_3B,  path1=subdir_7_part4, fileName1="4-7-3B-2-figure-sham-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_7_sham_3B <- ggplot(myframe_7_sham_3B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_sham_3B,  path1=subdir_7_part4, fileName1="4-7-3B-2-figure-sham-Line-noLess0",  height1=4, width1=7)



MyBoxViolinPlot_1(vector2=c(myNTR_7_banding_3A,  myNTR_7_sham_3B),       
                  sampleType2=c( rep("banding", numOfRows1)  ,   rep("sham", numOfRows1)  ), 
                  sampleRank2=c( "banding",  "sham" ),     
                  colours2=c( "banding"="blue",  "sham"="green4"  ), 
                  path2=subdir_7_part4,     fileName2="4-7-3A-merge-TAC-NTR-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                  Ymin2=0,   Ymax2=1.5,    height2=5,   width2=2  )

MyBoxViolinPlot_1(vector2=c(HalfLife_7_banding_3A,  HalfLife_7_sham_3B),       
                  sampleType2=c( rep("banding", numOfRows1)  ,   rep("sham", numOfRows1)  ), 
                  sampleRank2=c( "banding",  "sham" ),     
                  colours2=c( "banding"="blue",  "sham"="green4"  ), 
                  path2=subdir_7_part4,     fileName2="4-7-3B-merge-TAC-HalfLife-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half life (weeks)",   
                  Ymin2=0,   Ymax2=12,    height2=5,   width2=2  )




diff1_7_TAC = myNTR_7_banding_3A - myNTR_7_sham_3B 
diff2_7_TAC = log2( ( (myNTR_7_banding_3A+ 0.001) / (myNTR_7_sham_3B+0.001))  )
diff3_7_TAC = log2( ( (myNTR_7_sham_3B+ 0.001) / (myNTR_7_banding_3A+0.001))  )
length(diff1_7_TAC)
length(diff2_7_TAC)
length(diff3_7_TAC)

myframe_7_TAC <- data.frame( xAxis = c(1:length(diff1_7_TAC)),      yAxis = sort(diff1_7_TAC) )

FigureTemp_7_TAC <- ggplot(myframe_7_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR (banding - sham)") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.1, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -0.1, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_TAC,  path1=subdir_7_part4, fileName1="4-7-3C-merge1-figure-TAC-Line-minus-limitY",  height1=4, width1=7)

FigureTemp_7_TAC <- ggplot(myframe_7_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR (banding - sham)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.1, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -0.1, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_TAC,  path1=subdir_7_part4, fileName1="4-7-3C-merge1-figure-TAC-Line-minus",  height1=4, width1=7)





myframe_7_TAC <- data.frame( xAxis = c(1:length(diff2_7_TAC)),      yAxis = sort(diff2_7_TAC) )

FigureTemp_7_TAC <- ggplot(myframe_7_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(banding/sham)") +  ggtitle(myTitle_g) + ylim(-5, 5 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept =1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_TAC,  path1=subdir_7_part4, fileName1="4-7-3C-merge2-figure-TAC-Line-ratio-limitY",  height1=4, width1=7)

FigureTemp_7_TAC <- ggplot(myframe_7_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(banding/sham)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_TAC,  path1=subdir_7_part4, fileName1="4-7-3C-merge2-figure-TAC-Line-ratio",  height1=4, width1=7)




myframe_7_TAC <- data.frame( xAxis = c(1:length(diff3_7_TAC)),      yAxis = sort(diff3_7_TAC) )

FigureTemp_7_TAC <- ggplot(myframe_7_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(sham/banding)") +  ggtitle(myTitle_g) + ylim(-5, 5 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept =1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_TAC,  path1=subdir_7_part4, fileName1="4-7-3C-merge3-figure-TAC-Line-ratio-limitY",  height1=4, width1=7)

FigureTemp_7_TAC <- ggplot(myframe_7_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(sham/banding)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_7_TAC,  path1=subdir_7_part4, fileName1="4-7-3C-merge3-figure-TAC-Line-ratio",  height1=4, width1=7)





























##########################################################################   
subdir_7A_part4 <- paste(Part4_g,  "/7A-NTR-eachRow-1kb", sep = "")
if( ! file.exists(subdir_7A_part4) ) { dir.create(subdir_7A_part4) }


dim(Average_H3)
dim(Average_week0)
numRows_oneClass <- floor( nrow(Average_H3)/5 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
INDEX_4 <- seq(from = numRows_oneClass*3+1, to = numRows_oneClass*4, by =1 )
INDEX_5 <- seq(from = numRows_oneClass*4+1, to = numRows_oneClass*5, by =1 )
nrow(Average_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) + length(INDEX_4) + length(INDEX_5))

MyBoxViolinPlot_1(vector2=c( myNTR_7_WT_1[INDEX_1], myNTR_7_WT_1[INDEX_2],  myNTR_7_WT_1[INDEX_3],  myNTR_7_WT_1[INDEX_4],  myNTR_7_WT_1[INDEX_5] ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),  rep("Medium", length(INDEX_3)),  
                                 rep("High", length(INDEX_4)),      rep("Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "Lowest",  "Low",   "Medium",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",    "Medium"="green2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_7A_part4,   fileName2= paste("4-7A-1-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=3.5,   Ymin2=0, Ymax2=0.4)    ## width = 1 + 5*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c( HalfLife_7_WT_1[INDEX_1], HalfLife_7_WT_1[INDEX_2],  HalfLife_7_WT_1[INDEX_3],  HalfLife_7_WT_1[INDEX_4],  HalfLife_7_WT_1[INDEX_5] ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),  rep("Medium", length(INDEX_3)),  
                                 rep("High", length(INDEX_4)),      rep("Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "Lowest",  "Low",   "Medium",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",    "Medium"="green2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_7A_part4,   fileName2= paste("4-7A-2-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=3.5,   Ymin2=0, Ymax2=50)    ## width = 1 + 5*0.5, height=5cm 











MyBoxViolinPlot_1(vector2=c( myNTR_7_EEDheto_2A[INDEX_1], myNTR_7_EEDheto_2A[INDEX_2],  myNTR_7_EEDheto_2A[INDEX_3],  myNTR_7_EEDheto_2A[INDEX_4],  myNTR_7_EEDheto_2A[INDEX_5],
                             myNTR_7_EEDko_2B[INDEX_1],   myNTR_7_EEDko_2B[INDEX_2],    myNTR_7_EEDko_2B[INDEX_3],    myNTR_7_EEDko_2B[INDEX_4],    myNTR_7_EEDko_2B[INDEX_5] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),  rep("EEDheto_Medium", length(INDEX_3)),  rep("EEDheto_High", length(INDEX_4)),  rep("EEDheto_Highest", length(INDEX_5)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),    rep("EEDko_Medium", length(INDEX_3)),    rep("EEDko_High", length(INDEX_4)),    rep("EEDko_Highest", length(INDEX_5)) ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",     "EEDko_Medium",    "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2",    "EEDko_Highest"="purple2" ), 
                  path2=subdir_7A_part4,   fileName2= paste("4-7A-3-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=0.3)    ## width = 1 + 10*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_7_EEDheto_2A[INDEX_1], HalfLife_7_EEDheto_2A[INDEX_2],  HalfLife_7_EEDheto_2A[INDEX_3],  HalfLife_7_EEDheto_2A[INDEX_4],  HalfLife_7_EEDheto_2A[INDEX_5],
                             HalfLife_7_EEDko_2B[INDEX_1],   HalfLife_7_EEDko_2B[INDEX_2],    HalfLife_7_EEDko_2B[INDEX_3],    HalfLife_7_EEDko_2B[INDEX_4],    HalfLife_7_EEDko_2B[INDEX_5] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),  rep("EEDheto_Medium", length(INDEX_3)),  rep("EEDheto_High", length(INDEX_4)),  rep("EEDheto_Highest", length(INDEX_5)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),    rep("EEDko_Medium", length(INDEX_3)),    rep("EEDko_High", length(INDEX_4)),    rep("EEDko_Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",     "EEDko_Medium",    "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2",    "EEDko_Highest"="purple2"  ), 
                  path2=subdir_7A_part4,   fileName2= paste("4-7A-4-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=40)    ## width = 1 + 10*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( myNTR_7_EEDheto_2A[INDEX_1], myNTR_7_EEDheto_2A[INDEX_2],  myNTR_7_EEDheto_2A[INDEX_3],  myNTR_7_EEDheto_2A[INDEX_4],  myNTR_7_EEDheto_2A[INDEX_5],
                             myNTR_7_EEDko_2B[INDEX_1],   myNTR_7_EEDko_2B[INDEX_2],    myNTR_7_EEDko_2B[INDEX_3],    myNTR_7_EEDko_2B[INDEX_4],    myNTR_7_EEDko_2B[INDEX_5] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),  rep("EEDheto_Medium", length(INDEX_3)),  rep("EEDheto_High", length(INDEX_4)),  rep("EEDheto_Highest", length(INDEX_5)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),    rep("EEDko_Medium", length(INDEX_3)),    rep("EEDko_High", length(INDEX_4)),    rep("EEDko_Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "EEDheto_Lowest",   "EEDko_Lowest",   
                                 "EEDheto_Low",      "EEDko_Low", 
                                 "EEDheto_Medium",   "EEDko_Medium", 
                                 "EEDheto_High",     "EEDko_High",  
                                 "EEDheto_Highest",  "EEDko_Highest" ),     
                  colours2=c( "EEDheto_Lowest"="red2",        "EEDko_Lowest"="red2",   
                              "EEDheto_Low"="cyan2",          "EEDko_Low"="cyan2",  
                              "EEDheto_Medium"="green2",      "EEDko_Medium"="green2", 
                              "EEDheto_High"="blue2",         "EEDko_High"="blue2",  
                              "EEDheto_Highest"="purple2" ,   "EEDko_Highest"="purple2" ), 
                  path2=subdir_7A_part4,   fileName2= paste("4-7A-5-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0,   Ymax2=0.3)    ## width = 1 + 10*0.5, height=5cm 









MyBoxViolinPlot_1(vector2=c( myNTR_7_banding_3A[INDEX_1], myNTR_7_banding_3A[INDEX_2],  myNTR_7_banding_3A[INDEX_3],  myNTR_7_banding_3A[INDEX_4],  myNTR_7_banding_3A[INDEX_5],
                             myNTR_7_sham_3B[INDEX_1],    myNTR_7_sham_3B[INDEX_2],     myNTR_7_sham_3B[INDEX_3],     myNTR_7_sham_3B[INDEX_4],     myNTR_7_sham_3B[INDEX_5]  ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),  rep("banding_Medium", length(INDEX_3)),  rep("banding_High", length(INDEX_4)),  rep("banding_Highest", length(INDEX_5)), 
                                 rep("sham_Lowest", length(INDEX_1)),     rep("sham_Low", length(INDEX_2)),     rep("sham_Medium", length(INDEX_3)),     rep("sham_High", length(INDEX_4)),     rep("sham_Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",   "banding_Medium",  "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",     "sham_Low",      "sham_Medium",     "sham_High",     "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",      "sham_Low"="cyan2",       "sham_Medium"="green2",     "sham_High"="blue2",     "sham_Highest"="purple2"  ), 
                  path2=subdir_7A_part4,   fileName2= paste("4-7A-6-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=1.5)    ## width = 1 + 10*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_7_banding_3A[INDEX_1], HalfLife_7_banding_3A[INDEX_2],  HalfLife_7_banding_3A[INDEX_3],  HalfLife_7_banding_3A[INDEX_4],  HalfLife_7_banding_3A[INDEX_5],
                             HalfLife_7_sham_3B[INDEX_1],    HalfLife_7_sham_3B[INDEX_2],     HalfLife_7_sham_3B[INDEX_3],     HalfLife_7_sham_3B[INDEX_4],     HalfLife_7_sham_3B[INDEX_5] ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),  rep("banding_Medium", length(INDEX_3)),  rep("banding_High", length(INDEX_4)),  rep("banding_Highest", length(INDEX_5)), 
                                 rep("sham_Lowest", length(INDEX_1)),    rep("sham_Low", length(INDEX_2)),    rep("sham_Medium", length(INDEX_3)),    rep("sham_High", length(INDEX_4)),    rep("sham_Highest", length(INDEX_5))
                  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",   "banding_Medium",  "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",    "sham_Low",     "sham_Medium",    "sham_High",    "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",     "sham_Low"="cyan2",      "sham_Medium"="green2",    "sham_High"="blue2",    "sham_Highest"="purple2" ), 
                  path2=subdir_7A_part4,   fileName2= paste("4-7A-7-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=10)    ## width = 1 + 10*0.5, height=5cm 



MyBoxViolinPlot_1(vector2=c( myNTR_7_banding_3A[INDEX_1], myNTR_7_banding_3A[INDEX_2],  myNTR_7_banding_3A[INDEX_3],  myNTR_7_banding_3A[INDEX_4],  myNTR_7_banding_3A[INDEX_5],
                             myNTR_7_sham_3B[INDEX_1],   myNTR_7_sham_3B[INDEX_2],    myNTR_7_sham_3B[INDEX_3],    myNTR_7_sham_3B[INDEX_4],    myNTR_7_sham_3B[INDEX_5] ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),  rep("banding_Medium", length(INDEX_3)),  rep("banding_High", length(INDEX_4)),  rep("banding_Highest", length(INDEX_5)), 
                                 rep("sham_Lowest", length(INDEX_1)),    rep("sham_Low", length(INDEX_2)),    rep("sham_Medium", length(INDEX_3)),    rep("sham_High", length(INDEX_4)),    rep("sham_Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "banding_Lowest",   "sham_Lowest",   
                                 "banding_Low",      "sham_Low", 
                                 "banding_Medium",   "sham_Medium", 
                                 "banding_High",     "sham_High",  
                                 "banding_Highest",  "sham_Highest" ),     
                  colours2=c( "banding_Lowest"="red2",        "sham_Lowest"="red2",   
                              "banding_Low"="cyan2",          "sham_Low"="cyan2",  
                              "banding_Medium"="green2",      "sham_Medium"="green2", 
                              "banding_High"="blue2",         "sham_High"="blue2",  
                              "banding_Highest"="purple2" ,   "sham_Highest"="purple2"  ), 
                  path2=subdir_7A_part4,   fileName2= paste("4-7A-8-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0,   Ymax2=1.5)    ## width = 1 + 10*0.5, height=5cm 































##########################################################################   
subdir_7B_part4 <- paste(Part4_g,  "/7B-NTR-eachRow-1kb", sep = "")
if( ! file.exists(subdir_7B_part4) ) { dir.create(subdir_7B_part4) }


dim(Average_H3)
dim(Average_week0)
numRows_oneClass <- floor( nrow(Average_H3)/4 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
INDEX_4 <- seq(from = numRows_oneClass*3+1, to = numRows_oneClass*4, by =1 )
nrow(Average_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) + length(INDEX_4) )

MyBoxViolinPlot_1(vector2=c( myNTR_7_WT_1[INDEX_1], myNTR_7_WT_1[INDEX_2],  myNTR_7_WT_1[INDEX_3],  myNTR_7_WT_1[INDEX_4]  ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),  rep("High", length(INDEX_3)),      rep("Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "Lowest",  "Low",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-7B-1-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=3.0,   Ymin2=0, Ymax2=0.4)    ## width = 1 + 4*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c( HalfLife_7_WT_1[INDEX_1], HalfLife_7_WT_1[INDEX_2],  HalfLife_7_WT_1[INDEX_3],  HalfLife_7_WT_1[INDEX_4]  ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),  rep("High", length(INDEX_3)),      rep("Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "Lowest",  "Low",    "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",   "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-7B-2-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=3.0,   Ymin2=0, Ymax2=50)    ## width = 1 + 4*0.5, height=5cm 











MyBoxViolinPlot_1(vector2=c( myNTR_7_EEDheto_2A[INDEX_1], myNTR_7_EEDheto_2A[INDEX_2],  myNTR_7_EEDheto_2A[INDEX_3],  myNTR_7_EEDheto_2A[INDEX_4],  
                             myNTR_7_EEDko_2B[INDEX_1],   myNTR_7_EEDko_2B[INDEX_2],    myNTR_7_EEDko_2B[INDEX_3],    myNTR_7_EEDko_2B[INDEX_4] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),       rep("EEDheto_High", length(INDEX_3)),  rep("EEDheto_Highest", length(INDEX_4)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),         rep("EEDko_High", length(INDEX_3)),    rep("EEDko_Highest", length(INDEX_4)) ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",        "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",          "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",         "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",           "EEDko_High"="blue2",    "EEDko_Highest"="purple2" ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-7B-3-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=0.3)    ## width = 1 + 8*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_7_EEDheto_2A[INDEX_1], HalfLife_7_EEDheto_2A[INDEX_2],  HalfLife_7_EEDheto_2A[INDEX_3],  HalfLife_7_EEDheto_2A[INDEX_4],   
                             HalfLife_7_EEDko_2B[INDEX_1],   HalfLife_7_EEDko_2B[INDEX_2],    HalfLife_7_EEDko_2B[INDEX_3],    HalfLife_7_EEDko_2B[INDEX_4]  ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),       rep("EEDheto_High", length(INDEX_3)),  rep("EEDheto_Highest", length(INDEX_4)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),         rep("EEDko_High", length(INDEX_3)),    rep("EEDko_Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",      "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",        "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",     "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",       "EEDko_High"="blue2",    "EEDko_Highest"="purple2"  ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-7B-4-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=40)    ## width = 1 + 8*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( myNTR_7_EEDheto_2A[INDEX_1], myNTR_7_EEDheto_2A[INDEX_2],  myNTR_7_EEDheto_2A[INDEX_3],  myNTR_7_EEDheto_2A[INDEX_4],  
                             myNTR_7_EEDko_2B[INDEX_1],   myNTR_7_EEDko_2B[INDEX_2],    myNTR_7_EEDko_2B[INDEX_3],    myNTR_7_EEDko_2B[INDEX_4]  ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),    rep("EEDheto_High", length(INDEX_3)),  rep("EEDheto_Highest", length(INDEX_4)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),      rep("EEDko_High", length(INDEX_3)),    rep("EEDko_Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "EEDheto_Lowest",   "EEDko_Lowest",   
                                 "EEDheto_Low",      "EEDko_Low", 
                                 "EEDheto_High",     "EEDko_High",  
                                 "EEDheto_Highest",  "EEDko_Highest" ),     
                  colours2=c( "EEDheto_Lowest"="red2",        "EEDko_Lowest"="red2",   
                              "EEDheto_Low"="cyan2",          "EEDko_Low"="cyan2",  
                              "EEDheto_High"="blue2",         "EEDko_High"="blue2",  
                              "EEDheto_Highest"="purple2" ,   "EEDko_Highest"="purple2" ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-7B-5-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=5,   Ymin2=0,   Ymax2=0.3)    ## width = 1 + 8*0.5, height=5cm 









MyBoxViolinPlot_1(vector2=c( myNTR_7_banding_3A[INDEX_1], myNTR_7_banding_3A[INDEX_2],  myNTR_7_banding_3A[INDEX_3],  myNTR_7_banding_3A[INDEX_4],  
                             myNTR_7_sham_3B[INDEX_1],    myNTR_7_sham_3B[INDEX_2],     myNTR_7_sham_3B[INDEX_3],     myNTR_7_sham_3B[INDEX_4]   ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),   rep("banding_High", length(INDEX_3)),  rep("banding_Highest", length(INDEX_4)), 
                                 rep("sham_Lowest", length(INDEX_1)),     rep("sham_Low", length(INDEX_2)),      rep("sham_High", length(INDEX_3)),     rep("sham_Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",       "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",     "sham_Low",          "sham_High",     "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",        "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",      "sham_Low"="cyan2",            "sham_High"="blue2",     "sham_Highest"="purple2"  ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-7B-6-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=1.5)    ## width = 1 + 8*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_7_banding_3A[INDEX_1], HalfLife_7_banding_3A[INDEX_2],  HalfLife_7_banding_3A[INDEX_3],  HalfLife_7_banding_3A[INDEX_4],  
                             HalfLife_7_sham_3B[INDEX_1],    HalfLife_7_sham_3B[INDEX_2],     HalfLife_7_sham_3B[INDEX_3],     HalfLife_7_sham_3B[INDEX_4]  ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),     rep("banding_High", length(INDEX_3)),  rep("banding_Highest", length(INDEX_4)), 
                                 rep("sham_Lowest", length(INDEX_1)),    rep("sham_Low", length(INDEX_2)),         rep("sham_High", length(INDEX_3)),    rep("sham_Highest", length(INDEX_4))
                  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",     "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",    "sham_Low",          "sham_High",    "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",      "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",     "sham_Low"="cyan2",           "sham_High"="blue2",    "sham_Highest"="purple2" ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-7B-7-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=10)    ## width = 1 + 8*0.5, height=5cm 



MyBoxViolinPlot_1(vector2=c( myNTR_7_banding_3A[INDEX_1], myNTR_7_banding_3A[INDEX_2],  myNTR_7_banding_3A[INDEX_3],  myNTR_7_banding_3A[INDEX_4],   
                             myNTR_7_sham_3B[INDEX_1],   myNTR_7_sham_3B[INDEX_2],    myNTR_7_sham_3B[INDEX_3],    myNTR_7_sham_3B[INDEX_4]  ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),    rep("banding_High", length(INDEX_3)),  rep("banding_Highest", length(INDEX_4)), 
                                 rep("sham_Lowest", length(INDEX_1)),    rep("sham_Low", length(INDEX_2)),        rep("sham_High", length(INDEX_3)),    rep("sham_Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "banding_Lowest",   "sham_Lowest",   
                                 "banding_Low",      "sham_Low", 
                                 "banding_High",     "sham_High",  
                                 "banding_Highest",  "sham_Highest" ),     
                  colours2=c( "banding_Lowest"="red2",        "sham_Lowest"="red2",   
                              "banding_Low"="cyan2",          "sham_Low"="cyan2",  
                              "banding_High"="blue2",         "sham_High"="blue2",  
                              "banding_Highest"="purple2" ,   "sham_Highest"="purple2"  ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-7B-8-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=5,   Ymin2=0,   Ymax2=1.5)    ## width = 1 + 8*0.5, height=5cm 









































##########################################################################   
subdir_7C_part4 <- paste(Part4_g,  "/7C-NTR-eachRow-1kb", sep = "")
if( ! file.exists(subdir_7C_part4) ) { dir.create(subdir_7C_part4) }


dim(Average_H3)
dim(Average_week0)
numRows_oneClass <- floor( nrow(Average_H3)/3 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
nrow(Average_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) )

MyBoxViolinPlot_1(vector2=c( myNTR_7_WT_1[INDEX_1], myNTR_7_WT_1[INDEX_2],  myNTR_7_WT_1[INDEX_3]  ),   
                  sampleType2=c(  rep("Low", length(INDEX_1)),    rep("Medium", length(INDEX_2)),    rep("High", length(INDEX_3))   ), 
                  sampleRank2=c(  "Low",   "Medium",  "High"  ),     
                  colours2=c(  "Low"="cyan2",    "Medium"="green2",  "High"="blue2"   ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-7C-1-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=2.5,   Ymin2=0, Ymax2=0.4)    ## width = 1 + 3*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c( HalfLife_7_WT_1[INDEX_1], HalfLife_7_WT_1[INDEX_2],  HalfLife_7_WT_1[INDEX_3]  ),   
                  sampleType2=c(  rep("Low", length(INDEX_1)),  rep("Medium", length(INDEX_2)),   rep("High", length(INDEX_3))   ), 
                  sampleRank2=c(  "Low",   "Medium",  "High"  ),     
                  colours2=c( "Low"="cyan2",    "Medium"="green2",  "High"="blue2"   ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-7C-2-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=2.5,   Ymin2=0, Ymax2=50)    ## width = 1 + 3*0.5, height=5cm 











MyBoxViolinPlot_1(vector2=c( myNTR_7_EEDheto_2A[INDEX_1], myNTR_7_EEDheto_2A[INDEX_2],  myNTR_7_EEDheto_2A[INDEX_3],  
                             myNTR_7_EEDko_2B[INDEX_1],   myNTR_7_EEDko_2B[INDEX_2],    myNTR_7_EEDko_2B[INDEX_3]  ),   
                  sampleType2=c(   rep("EEDheto_Low", length(INDEX_1)),  rep("EEDheto_Medium", length(INDEX_2)),  rep("EEDheto_High", length(INDEX_3)),    
                                   rep("EEDko_Low", length(INDEX_1)),    rep("EEDko_Medium", length(INDEX_2)),    rep("EEDko_High", length(INDEX_3))  ), 
                  sampleRank2=c(    "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",   
                                    "EEDko_Low",     "EEDko_Medium",    "EEDko_High"   ),     
                  colours2=c(    "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",   
                                 "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2"  ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-7C-3-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=4,   Ymin2=0, Ymax2=0.3)    ## width = 1 + 6*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_7_EEDheto_2A[INDEX_1], HalfLife_7_EEDheto_2A[INDEX_2],  HalfLife_7_EEDheto_2A[INDEX_3],   
                             HalfLife_7_EEDko_2B[INDEX_1],   HalfLife_7_EEDko_2B[INDEX_2],    HalfLife_7_EEDko_2B[INDEX_3]  ),   
                  sampleType2=c(    rep("EEDheto_Low", length(INDEX_1)),  rep("EEDheto_Medium", length(INDEX_2)),  rep("EEDheto_High", length(INDEX_3)),  
                                    rep("EEDko_Low", length(INDEX_1)),    rep("EEDko_Medium", length(INDEX_2)),    rep("EEDko_High", length(INDEX_3))  ), 
                  sampleRank2=c(  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  
                                  "EEDko_Low",     "EEDko_Medium",    "EEDko_High"   ),     
                  colours2=c(  "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",   
                               "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2"   ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-7C-4-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=4,   Ymin2=0, Ymax2=40)    ## width = 1 + 6*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( myNTR_7_EEDheto_2A[INDEX_1], myNTR_7_EEDheto_2A[INDEX_2],  myNTR_7_EEDheto_2A[INDEX_3],   
                             myNTR_7_EEDko_2B[INDEX_1],   myNTR_7_EEDko_2B[INDEX_2],    myNTR_7_EEDko_2B[INDEX_3]  ),   
                  sampleType2=c(    rep("EEDheto_Low", length(INDEX_1)),  rep("EEDheto_Medium", length(INDEX_2)),  rep("EEDheto_High", length(INDEX_3)),   
                                    rep("EEDko_Low", length(INDEX_1)),    rep("EEDko_Medium", length(INDEX_2)),    rep("EEDko_High", length(INDEX_3))  ), 
                  sampleRank2=c( "EEDheto_Low",      "EEDko_Low", 
                                 "EEDheto_Medium",   "EEDko_Medium", 
                                 "EEDheto_High",     "EEDko_High" ),     
                  colours2=c( "EEDheto_Low"="cyan2",          "EEDko_Low"="cyan2",  
                              "EEDheto_Medium"="green2",      "EEDko_Medium"="green2", 
                              "EEDheto_High"="blue2",         "EEDko_High"="blue2" ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-7C-5-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=4,   Ymin2=0,   Ymax2=0.3)    ## width = 1 + 6*0.5, height=5cm 









MyBoxViolinPlot_1(vector2=c( myNTR_7_banding_3A[INDEX_1], myNTR_7_banding_3A[INDEX_2],  myNTR_7_banding_3A[INDEX_3],   
                             myNTR_7_sham_3B[INDEX_1],    myNTR_7_sham_3B[INDEX_2],     myNTR_7_sham_3B[INDEX_3]    ),   
                  sampleType2=c(  rep("banding_Low", length(INDEX_1)),  rep("banding_Medium", length(INDEX_2)),  rep("banding_High", length(INDEX_3)),   
                                  rep("sham_Low", length(INDEX_1)),     rep("sham_Medium", length(INDEX_2)),     rep("sham_High", length(INDEX_3))   ), 
                  sampleRank2=c(  "banding_Low",   "banding_Medium",  "banding_High",  
                                  "sham_Low",      "sham_Medium",     "sham_High"   ),     
                  colours2=c(  "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  
                               "sham_Low"="cyan2",       "sham_Medium"="green2",     "sham_High"="blue2"   ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-7C-6-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=4,   Ymin2=0, Ymax2=1)    ## width = 1 + 6*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_7_banding_3A[INDEX_1], HalfLife_7_banding_3A[INDEX_2],  HalfLife_7_banding_3A[INDEX_3],  
                             HalfLife_7_sham_3B[INDEX_1],    HalfLife_7_sham_3B[INDEX_2],     HalfLife_7_sham_3B[INDEX_3]  ),   
                  sampleType2=c(   rep("banding_Low", length(INDEX_1)),  rep("banding_Medium", length(INDEX_2)),  rep("banding_High", length(INDEX_3)),   
                                   rep("sham_Low", length(INDEX_1)),     rep("sham_Medium", length(INDEX_2)),     rep("sham_High", length(INDEX_3)) ), 
                  sampleRank2=c(   "banding_Low",   "banding_Medium",  "banding_High",   
                                   "sham_Low",     "sham_Medium",    "sham_High"    ),     
                  colours2=c(   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",   
                                "sham_Low"="cyan2",      "sham_Medium"="green2",    "sham_High"="blue2"  ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-7C-7-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=4,   Ymin2=0, Ymax2=10)    ## width = 1 + 6*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( myNTR_7_banding_3A[INDEX_1], myNTR_7_banding_3A[INDEX_2],  myNTR_7_banding_3A[INDEX_3],   
                             myNTR_7_sham_3B[INDEX_1],   myNTR_7_sham_3B[INDEX_2],    myNTR_7_sham_3B[INDEX_3]  ),   
                  sampleType2=c(   rep("banding_Low", length(INDEX_1)),  rep("banding_Medium", length(INDEX_2)),  rep("banding_High", length(INDEX_3)),   
                                   rep("sham_Low", length(INDEX_1)),     rep("sham_Medium", length(INDEX_2)),    rep("sham_High", length(INDEX_3))  ), 
                  sampleRank2=c( "banding_Low",      "sham_Low", 
                                 "banding_Medium",   "sham_Medium", 
                                 "banding_High",     "sham_High"  ),     
                  colours2=c( "banding_Low"="cyan2",          "sham_Low"="cyan2",  
                              "banding_Medium"="green2",      "sham_Medium"="green2", 
                              "banding_High"="blue2",         "sham_High"="blue2"   ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-7C-8-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=4,   Ymin2=0,   Ymax2=1)    ## width = 1 + 6*0.5, height=5cm 































#########################################################################################################  average some columns for each row
subdir_8_part4 <- paste(Part4_g,  "/8-NTR-eachRow-500bp", sep = "")
if( ! file.exists(subdir_8_part4) ) { dir.create(subdir_8_part4) }


############################# WT 6 samples
sink( file=paste(subdir_8_part4,  "/4-8-1A-WT-runLog.txt",      sep = "") )
myNTR_8_WT_1  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(reduceColumn4_Average_week0[i]+0.0006, reduceColumn4_Average_week1[i]+0.0005, reduceColumn4_Average_week2[i]+0.0004, reduceColumn4_Average_week4[i]+0.0003, reduceColumn4_Average_week6[i]+0.0002, reduceColumn4_Average_week8[i]+0.0001)                        
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log(vec1)  , file2=paste(subdir_8_part4,  "/4-8-1B-WT-LogLinearModel", sep="") )
  myNTR_8_WT_1 <- c(myNTR_8_WT_1, NTR_bin)
}
sink()  

length(myNTR_8_WT_1)
summary(myNTR_8_WT_1)
HalfLife_8_WT_1 <- log(2)/(myNTR_8_WT_1+0.0001)
length(HalfLife_8_WT_1)
summary(HalfLife_8_WT_1)

write.table(x=myNTR_8_WT_1,      file = paste(subdir_8_part4,  "/4-8-1C-raw-NTR.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_8_WT_1,   file = paste(subdir_8_part4,  "/4-8-1D-raw-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_8_WT_1 <- data.frame( xAxis = c(1:length(myNTR_8_WT_1)),      yAxis = sort(myNTR_8_WT_1) )
FigureTemp_8_WT_1 <- ggplot(myframe_8_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_WT_1,  path1=subdir_8_part4, fileName1="4-8-1A-figure-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_8_WT_1 <- ggplot(myframe_8_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) +  
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_WT_1,  path1=subdir_8_part4, fileName1="4-8-1A-figure-Line-raw",  height1=4, width1=7)

myNTR_8_WT_1[myNTR_8_WT_1<0]  <- 0
myNTR_8_WT_1[myNTR_8_WT_1>1]  <- 1
HalfLife_8_WT_1[HalfLife_8_WT_1<0]  <- 0
HalfLife_8_WT_1[HalfLife_8_WT_1>50] <- 50

write.table(x=thisRowNames1,     file = paste(subdir_8_part4,   "/4-8-1E-RowNames.txt",  sep = ""),         append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=myNTR_8_WT_1,      file = paste(subdir_8_part4,  "/4-8-1F-noLess0-NTR.txt",    sep = ""),     append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_8_WT_1,   file = paste(subdir_8_part4,  "/4-8-1G-noLess0-HalfLife.txt", sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_8_WT_1 <- data.frame( xAxis = c(1:length(myNTR_8_WT_1)),      yAxis = sort(myNTR_8_WT_1) )
FigureTemp_8_WT_1 <- ggplot(myframe_8_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_WT_1,  path1=subdir_8_part4, fileName1="4-8-1B-figure-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_8_WT_1 <- ggplot(myframe_8_WT_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) +  
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_WT_1,  path1=subdir_8_part4, fileName1="4-8-1B-figure-Line-noLess0",  height1=4, width1=7)

MyBoxViolinPlot_1(vector2=myNTR_8_WT_1,      sampleType2=c( rep("NTR", numOfRows1)  ), 
                  sampleRank2=c( "NTR" ),     colours2=c( "NTR"="red" ), 
                  path2=subdir_8_part4,     fileName2="4-8-1C-WT-NTR-average-curve",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="NTR",   
                  Ymin2=0,   Ymax2=0.7,    height2=5,   width2=1.5  )

MyBoxViolinPlot_1(vector2=HalfLife_8_WT_1,      sampleType2=c( rep("HalfLife", numOfRows1)  ), 
                  sampleRank2=c( "HalfLife" ),     colours2=c( "HalfLife"="red" ), 
                  path2=subdir_8_part4,     fileName2="4-8-1D-WT-HalfLife-average-curve",  
                  title2=myTitle_g,     xLab2="Samples",    yLab2="Half life (weeks)",   
                  Ymin2=0,   Ymax2=15,    height2=5,   width2=1.5  )










############################# CKO 4 samples
sink( file=paste(subdir_8_part4,  "/4-8-2A-1-EEDheto-runLog.txt",      sep = "") )
myNTR_8_EEDheto_2A  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(reduceColumn4_Average_week0_EEDheto[i]+0.0001,  reduceColumn4_Average_week4_EEDheto[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_8_part4,  "/4-8-2A-2-EEDheto-LogLinearModel", sep="") )
  myNTR_8_EEDheto_2A <- c(myNTR_8_EEDheto_2A, NTR_bin)
}
sink()  

length(myNTR_8_EEDheto_2A)
summary(myNTR_8_EEDheto_2A)
HalfLife_8_EEDheto_2A <- log(2)/(myNTR_8_EEDheto_2A+0.0001)
length(HalfLife_8_EEDheto_2A)
summary(HalfLife_8_EEDheto_2A)

write.table(x=myNTR_8_EEDheto_2A,      file = paste(subdir_8_part4,  "/4-8-2A-3-raw-NTR_EEDheto.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_8_EEDheto_2A,   file = paste(subdir_8_part4,  "/4-8-2A-4-raw-HalfLife_EEDheto.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_8_EEDheto_2A <- data.frame( xAxis = c(1:length(myNTR_8_EEDheto_2A)),      yAxis = sort(myNTR_8_EEDheto_2A) )
FigureTemp_8_EEDheto_2A <- ggplot(myframe_8_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_EEDheto_2A,  path1=subdir_8_part4, fileName1="4-8-2A-1-figure-EEDheto-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_8_EEDheto_2A <- ggplot(myframe_8_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) +
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_EEDheto_2A,  path1=subdir_8_part4, fileName1="4-8-2A-1-figure-EEDheto-Line-raw",  height1=4, width1=7)

myNTR_8_EEDheto_2A[myNTR_8_EEDheto_2A<0]  <- 0
myNTR_8_EEDheto_2A[myNTR_8_EEDheto_2A>1]  <- 1
HalfLife_8_EEDheto_2A[HalfLife_8_EEDheto_2A<0]   <- 0
HalfLife_8_EEDheto_2A[HalfLife_8_EEDheto_2A>50]  <- 50

write.table(x=myNTR_8_EEDheto_2A,      file = paste(subdir_8_part4,  "/4-8-2A-8-noLess0-myNTR_EEDheto.txt",  sep = ""),      append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_8_EEDheto_2A,   file = paste(subdir_8_part4,  "/4-8-2A-8-noLess0-HalfLife_EEDheto.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_8_EEDheto_2A <- data.frame( xAxis = c(1:length(myNTR_8_EEDheto_2A)),      yAxis = sort(myNTR_8_EEDheto_2A) )
FigureTemp_8_EEDheto_2A <- ggplot(myframe_8_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_EEDheto_2A,  path1=subdir_8_part4, fileName1="4-8-2A-2-figure-EEDheto-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_8_EEDheto_2A <- ggplot(myframe_8_EEDheto_2A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_EEDheto_2A,  path1=subdir_8_part4, fileName1="4-8-2A-2-figure-EEDheto-Line-noLess0",  height1=4, width1=7)





sink( file=paste(subdir_8_part4,  "/4-8-2B-1-EEDko-runLog.txt",      sep = "") )
myNTR_8_EEDko_2B  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(reduceColumn4_Average_week0_EEDko[i]+0.0001,  reduceColumn4_Average_week4_EEDko[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  4),  yAxis2=log( vec1 )  , file2=paste(subdir_8_part4,  "/4-8-2B-2-EEDko-LogLinearModel", sep="") )
  myNTR_8_EEDko_2B <- c(myNTR_8_EEDko_2B, NTR_bin)
}
sink()  

length(myNTR_8_EEDko_2B)
summary(myNTR_8_EEDko_2B)
HalfLife_8_EEDko_2B <- log(2)/(myNTR_8_EEDko_2B+0.0001)
length(HalfLife_8_EEDko_2B)
summary(HalfLife_8_EEDko_2B)

write.table(x=myNTR_8_EEDko_2B,      file = paste(subdir_8_part4,  "/4-8-2B-3-raw-NTR_EEDko.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_8_EEDko_2B,   file = paste(subdir_8_part4,  "/4-8-2B-4-raw-HalfLife_EEDko.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_8_EEDko_2B <- data.frame( xAxis = c(1:length(myNTR_8_EEDko_2B)),      yAxis = sort(myNTR_8_EEDko_2B) )
FigureTemp_8_EEDko_2B <- ggplot(myframe_8_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_EEDko_2B,  path1=subdir_8_part4, fileName1="4-8-2B-1-figure-EEDko-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_8_EEDko_2B <- ggplot(myframe_8_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +    geom_hline(yintercept = 0.0, colour="blue") +   geom_hline(yintercept = -0.2, colour="blue") +  
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_EEDko_2B,  path1=subdir_8_part4, fileName1="4-8-2B-1-figure-EEDko-Line-raw",  height1=4, width1=7)

myNTR_8_EEDko_2B[myNTR_8_EEDko_2B<0]  <- 0
myNTR_8_EEDko_2B[myNTR_8_EEDko_2B>1]  <- 1
HalfLife_8_EEDko_2B[HalfLife_8_EEDko_2B<0]   <- 0
HalfLife_8_EEDko_2B[HalfLife_8_EEDko_2B>50]  <- 50

write.table(x=myNTR_8_EEDko_2B,      file = paste(subdir_8_part4,  "/4-8-2B-8-noLess0-myNTR_EEDko.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_8_EEDko_2B,   file = paste(subdir_8_part4,  "/4-8-2B-8-noLess0-HalfLife_EEDko.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_8_EEDko_2B <- data.frame( xAxis = c(1:length(myNTR_8_EEDko_2B)),      yAxis = sort(myNTR_8_EEDko_2B) )
FigureTemp_8_EEDko_2B <- ggplot(myframe_8_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_EEDko_2B,  path1=subdir_8_part4, fileName1="4-8-2B-2-figure-EEDko-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_8_EEDko_2B <- ggplot(myframe_8_EEDko_2B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_EEDko_2B,  path1=subdir_8_part4, fileName1="4-8-2B-2-figure-EEDko-Line-noLess0",  height1=4, width1=7)



MyBoxViolinPlot_1(vector2=c(myNTR_8_EEDheto_2A,  myNTR_8_EEDko_2B),       
                  sampleType2=c( rep("EEDheto", numOfRows1)  ,   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "EEDheto",  "EEDko" ),     
                  colours2=c( "EEDheto"="red",  "EEDko"="blue"  ), 
                  path2=subdir_8_part4,     fileName2="4-8-2A-merge-CKO-NTR-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                  Ymin2=0,   Ymax2=0.5,    height2=5,   width2=2  )

MyBoxViolinPlot_1(vector2=c(HalfLife_8_EEDheto_2A,  HalfLife_8_EEDko_2B),       
                  sampleType2=c( rep("EEDheto", numOfRows1)  ,   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "EEDheto",  "EEDko" ),     
                  colours2=c( "EEDheto"="red",  "EEDko"="blue"  ), 
                  path2=subdir_8_part4,     fileName2="4-8-2B-merge-CKO-HalfLife-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half life (weeks)",   
                  Ymin2=0,   Ymax2=15,    height2=5,   width2=2  )



diff1_8_CKO = myNTR_8_EEDko_2B - myNTR_8_EEDheto_2A
diff2_8_CKO = log2( ( (myNTR_8_EEDko_2B+0.001) / (myNTR_8_EEDheto_2A+0.001))  )
diff3_8_CKO = log2( ( (myNTR_8_EEDheto_2A+0.001) / (myNTR_8_EEDko_2B+0.001))  )   ## 2^(-10) = 0.001
length(diff1_8_CKO)
length(diff2_8_CKO)
length(diff3_8_CKO)
summary(diff1_8_CKO)
summary(diff2_8_CKO)
summary(diff3_8_CKO)


myframe_8_CKO <- data.frame( xAxis = c(1:length(diff1_8_CKO)),      yAxis = sort(diff1_8_CKO) )

FigureTemp_8_CKO <- ggplot(myframe_8_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR (EEDko - EEDheto)") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.1, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -0.1, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_CKO,  path1=subdir_8_part4, fileName1="4-8-2C-merge-figure-CKO-Line-minus-limitY",  height1=4, width1=7)

FigureTemp_8_CKO <- ggplot(myframe_8_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR (EEDko - EEDheto)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.1, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -0.1, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_CKO,  path1=subdir_8_part4, fileName1="4-8-2C-merge-figure-CKO-Line-minus",  height1=4, width1=7)



myframe_8_CKO <- data.frame( xAxis = c(1:length(diff2_8_CKO)),      yAxis = sort(diff2_8_CKO) )

FigureTemp_8_CKO <- ggplot(myframe_8_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDko/EEDheto)") +  ggtitle(myTitle_g) + ylim(-10, 10 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_CKO,  path1=subdir_8_part4, fileName1="4-8-2D-merge-figure-CKO-Line-ratio-limitY",  height1=4, width1=7)

FigureTemp_8_CKO <- ggplot(myframe_8_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDko/EEDheto)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_CKO,  path1=subdir_8_part4, fileName1="4-8-2D-merge-figure-CKO-Line-ratio",  height1=4, width1=7)





myframe_8_CKO <- data.frame( xAxis = c(1:length(diff3_8_CKO)),      yAxis = sort(diff3_8_CKO) )

FigureTemp_8_CKO <- ggplot(myframe_8_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDheto/EEDko)") +  ggtitle(myTitle_g) + ylim(-10, 10 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_CKO,  path1=subdir_8_part4, fileName1="4-8-2E-merge-figure-CKO-Line-ratio-limitY",  height1=4, width1=7)

FigureTemp_8_CKO <- ggplot(myframe_8_CKO,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(EEDheto/EEDko)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_CKO,  path1=subdir_8_part4, fileName1="4-8-2E-merge-figure-CKO-Line-ratio",  height1=4, width1=7)









############################# TAC 2 samples
sink( file=paste(subdir_8_part4,  "/4-8-3A-1-banding-runLog.txt",      sep = "") )
myNTR_8_banding_3A  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(reduceColumn4_Average_week0[i]+0.0001,  reduceColumn4_Average_banding[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  2.5),  yAxis2=log( vec1 )  , file2=paste(subdir_8_part4,  "/4-8-3A-2-banding-LogLinearModel", sep="") )
  myNTR_8_banding_3A <- c(myNTR_8_banding_3A, NTR_bin)
}
sink()  
length(myNTR_8_banding_3A)
summary(myNTR_8_banding_3A)
HalfLife_8_banding_3A <- log(2)/(myNTR_8_banding_3A+0.0001)
length(HalfLife_8_banding_3A)
summary(HalfLife_8_banding_3A)

write.table(x=myNTR_8_banding_3A,      file = paste(subdir_8_part4,  "/4-8-3A-3-banding-raw-myNTR.txt",  sep = ""),      append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_8_banding_3A,   file = paste(subdir_8_part4,  "/4-8-3A-4-banding-raw-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_8_banding_3A <- data.frame( xAxis = c(1:length(myNTR_8_banding_3A)),      yAxis = sort(myNTR_8_banding_3A) )
FigureTemp_8_banding_3A <- ggplot(myframe_8_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_banding_3A,  path1=subdir_8_part4, fileName1="4-8-3A-1-figure-banding-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_8_banding_3A <- ggplot(myframe_8_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_banding_3A,  path1=subdir_8_part4, fileName1="4-8-3A-1-figure-banding-Line-raw",  height1=4, width1=7)

myNTR_8_banding_3A[myNTR_8_banding_3A<0]  <- 0
myNTR_8_banding_3A[myNTR_8_banding_3A>10]  <- 10
HalfLife_8_banding_3A[HalfLife_8_banding_3A<0]  <- 0
HalfLife_8_banding_3A[HalfLife_8_banding_3A>50]  <- 50
min(myNTR_8_banding_3A)

write.table(x=myNTR_8_banding_3A,      file = paste(subdir_8_part4,  "/4-8-3A-8-noLess0-myNTR_8_banding_3A.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_8_banding_3A,   file = paste(subdir_8_part4,  "/4-8-3A-8-noLess0-HalfLife_8_banding_3A.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_8_banding_3A <- data.frame( xAxis = c(1:length(myNTR_8_banding_3A)),      yAxis = sort(myNTR_8_banding_3A) )
FigureTemp_8_banding_3A <- ggplot(myframe_8_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_banding_3A,  path1=subdir_8_part4, fileName1="4-8-3A-2-figure-banding-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_8_banding_3A <- ggplot(myframe_8_banding_3A,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_banding_3A,  path1=subdir_8_part4, fileName1="4-8-3A-2-figure-banding-Line-noLess0",  height1=4, width1=7)



sink( file=paste(subdir_8_part4,  "/4-8-3B-1-sham-runLog.txt",      sep = "") )
myNTR_8_sham_3B  <-  c()
for (i in c(1:numOfRows1)) {
  vec1 <- c(reduceColumn4_Average_week0[i]+0.0001,  reduceColumn4_Average_sham[i]+0.0001 )
  NTR_bin <- MyTurnoverRate_1( xAxis2=c(0,  2.5),  yAxis2=log( vec1 )  , file2=paste(subdir_8_part4,  "/4-8-3B-2-sham-LogLinearModel", sep="") )
  myNTR_8_sham_3B <- c(myNTR_8_sham_3B, NTR_bin)
}
sink()  
length(myNTR_8_sham_3B)
summary(myNTR_8_sham_3B)
HalfLife_8_sham_3B <- log(2)/(myNTR_8_sham_3B+0.0001)
length(HalfLife_8_sham_3B)
summary(HalfLife_8_sham_3B)

write.table(x=myNTR_8_sham_3B,      file = paste(subdir_8_part4,  "/4-8-3B-3-sham-raw-NTR.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_8_sham_3B,   file = paste(subdir_8_part4,  "/4-8-3B-4-sham-raw-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_8_sham_3B <- data.frame( xAxis = c(1:length(myNTR_8_sham_3B)),      yAxis = sort(myNTR_8_sham_3B) )
FigureTemp_8_sham_3B <- ggplot(myframe_8_sham_3B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_sham_3B,  path1=subdir_8_part4, fileName1="4-8-3B-1-figure-sham-Line-raw-limitY",  height1=4, width1=7)
FigureTemp_8_sham_3B <- ggplot(myframe_8_sham_3B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + geom_hline(yintercept = 0.0, colour="blue") + geom_hline(yintercept = -0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_sham_3B,  path1=subdir_8_part4, fileName1="4-8-3B-1-figure-sham-Line-raw",  height1=4, width1=7)

myNTR_8_sham_3B[myNTR_8_sham_3B<0]  <- 0
myNTR_8_sham_3B[myNTR_8_sham_3B>10]  <- 10
HalfLife_8_sham_3B[HalfLife_8_sham_3B<0]   <- 0
HalfLife_8_sham_3B[HalfLife_8_sham_3B>50]  <- 50

write.table(x=myNTR_8_sham_3B,      file = paste(subdir_8_part4,  "/4-8-3B-8-sham-noLess0-NTR.txt",  sep = ""),        append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
write.table(x=HalfLife_8_sham_3B,   file = paste(subdir_8_part4,  "/4-8-3B-8-sham-noLess0-HalfLife.txt",  sep = ""),   append = FALSE, quote = FALSE, sep = "\t",  row.names = TRUE)
myframe_8_sham_3B <- data.frame( xAxis = c(1:length(myNTR_8_sham_3B)),      yAxis = sort(myNTR_8_sham_3B) )
FigureTemp_8_sham_3B <- ggplot(myframe_8_sham_3B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + ylim(0, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_sham_3B,  path1=subdir_8_part4, fileName1="4-8-3B-2-figure-sham-Line-noLess0-limitY",  height1=4, width1=7)
FigureTemp_8_sham_3B <- ggplot(myframe_8_sham_3B,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.2, colour="blue") +
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_sham_3B,  path1=subdir_8_part4, fileName1="4-8-3B-2-figure-sham-Line-noLess0",  height1=4, width1=7)



MyBoxViolinPlot_1(vector2=c(myNTR_8_banding_3A,  myNTR_8_sham_3B),       
                  sampleType2=c( rep("banding", numOfRows1)  ,   rep("sham", numOfRows1)  ), 
                  sampleRank2=c( "banding",  "sham" ),     
                  colours2=c( "banding"="blue",  "sham"="green4"  ), 
                  path2=subdir_8_part4,     fileName2="4-8-3A-merge-TAC-NTR-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="NTR",   
                  Ymin2=0,   Ymax2=2.0,    height2=5,   width2=2  )

MyBoxViolinPlot_1(vector2=c(HalfLife_8_banding_3A,  HalfLife_8_sham_3B),       
                  sampleType2=c( rep("banding", numOfRows1)  ,   rep("sham", numOfRows1)  ), 
                  sampleRank2=c( "banding",  "sham" ),     
                  colours2=c( "banding"="blue",  "sham"="green4"  ), 
                  path2=subdir_8_part4,     fileName2="4-8-3B-merge-TAC-HalfLife-average-curve",  
                  title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="Half life (weeks)",   
                  Ymin2=0,   Ymax2=5,    height2=5,   width2=2  )




diff1_8_TAC = myNTR_8_banding_3A - myNTR_8_sham_3B 
diff2_8_TAC = log2( ( (myNTR_8_banding_3A+ 0.001) / (myNTR_8_sham_3B+0.001))  )
diff3_8_TAC = log2( ( (myNTR_8_sham_3B+ 0.001) / (myNTR_8_banding_3A+0.001))  )
length(diff1_8_TAC)
length(diff2_8_TAC)
length(diff3_8_TAC)

myframe_8_TAC <- data.frame( xAxis = c(1:length(diff1_8_TAC)),      yAxis = sort(diff1_8_TAC) )

FigureTemp_8_TAC <- ggplot(myframe_8_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR (banding - sham)") +  ggtitle(myTitle_g) + ylim(-1, 1 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.1, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -0.1, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_TAC,  path1=subdir_8_part4, fileName1="4-8-3C-merge1-figure-TAC-Line-minus-limitY",  height1=4, width1=7)

FigureTemp_8_TAC <- ggplot(myframe_8_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("NTR (banding - sham)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 0.1, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -0.1, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_TAC,  path1=subdir_8_part4, fileName1="4-8-3C-merge1-figure-TAC-Line-minus",  height1=4, width1=7)





myframe_8_TAC <- data.frame( xAxis = c(1:length(diff2_8_TAC)),      yAxis = sort(diff2_8_TAC) )

FigureTemp_8_TAC <- ggplot(myframe_8_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(banding/sham)") +  ggtitle(myTitle_g) + ylim(-5, 5 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept =1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_TAC,  path1=subdir_8_part4, fileName1="4-8-3C-merge2-figure-TAC-Line-ratio-limitY",  height1=4, width1=7)

FigureTemp_8_TAC <- ggplot(myframe_8_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(banding/sham)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_TAC,  path1=subdir_8_part4, fileName1="4-8-3C-merge2-figure-TAC-Line-ratio",  height1=4, width1=7)




myframe_8_TAC <- data.frame( xAxis = c(1:length(diff3_8_TAC)),      yAxis = sort(diff3_8_TAC) )

FigureTemp_8_TAC <- ggplot(myframe_8_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(sham/banding)") +  ggtitle(myTitle_g) + ylim(-5, 5 ) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept =1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_TAC,  path1=subdir_8_part4, fileName1="4-8-3C-merge3-figure-TAC-Line-ratio-limitY",  height1=4, width1=7)

FigureTemp_8_TAC <- ggplot(myframe_8_TAC,   aes(x=xAxis, y=yAxis )  )  +  xlab(myTitle_g) +  ylab("log2(sham/banding)") +  ggtitle(myTitle_g) + 
  geom_line(size=0.5, colour="red") + geom_hline(yintercept = 1.5, colour="blue") + geom_hline(yintercept = 0, colour="blue") + geom_hline(yintercept = -1.5, colour="blue") + 
  MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
MySaveGgplot2_1(ggplot2Figure1=FigureTemp_8_TAC,  path1=subdir_8_part4, fileName1="4-8-3C-merge3-figure-TAC-Line-ratio",  height1=4, width1=7)
































##########################################################################   
subdir_8A_part4 <- paste(Part4_g,  "/8A-NTR-eachRow-500bp", sep = "")
if( ! file.exists(subdir_8A_part4) ) { dir.create(subdir_8A_part4) }


dim(Average_H3)
dim(Average_week0)
numRows_oneClass <- floor( nrow(Average_H3)/5 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
INDEX_4 <- seq(from = numRows_oneClass*3+1, to = numRows_oneClass*4, by =1 )
INDEX_5 <- seq(from = numRows_oneClass*4+1, to = numRows_oneClass*5, by =1 )
nrow(Average_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) + length(INDEX_4) + length(INDEX_5))

MyBoxViolinPlot_1(vector2=c( myNTR_8_WT_1[INDEX_1], myNTR_8_WT_1[INDEX_2],  myNTR_8_WT_1[INDEX_3],  myNTR_8_WT_1[INDEX_4],  myNTR_8_WT_1[INDEX_5] ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),  rep("Medium", length(INDEX_3)),  
                                 rep("High", length(INDEX_4)),      rep("Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "Lowest",  "Low",   "Medium",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",    "Medium"="green2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_8A_part4,   fileName2= paste("4-8A-1-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=3.5,   Ymin2=0, Ymax2=0.7)    ## width = 1 + 5*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c( HalfLife_8_WT_1[INDEX_1], HalfLife_8_WT_1[INDEX_2],  HalfLife_8_WT_1[INDEX_3],  HalfLife_8_WT_1[INDEX_4],  HalfLife_8_WT_1[INDEX_5] ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),  rep("Medium", length(INDEX_3)),  
                                 rep("High", length(INDEX_4)),      rep("Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "Lowest",  "Low",   "Medium",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",    "Medium"="green2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_8A_part4,   fileName2= paste("4-8A-2-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=3.5,   Ymin2=0, Ymax2=15)    ## width = 1 + 5*0.5, height=5cm 











MyBoxViolinPlot_1(vector2=c( myNTR_8_EEDheto_2A[INDEX_1], myNTR_8_EEDheto_2A[INDEX_2],  myNTR_8_EEDheto_2A[INDEX_3],  myNTR_8_EEDheto_2A[INDEX_4],  myNTR_8_EEDheto_2A[INDEX_5],
                             myNTR_8_EEDko_2B[INDEX_1],   myNTR_8_EEDko_2B[INDEX_2],    myNTR_8_EEDko_2B[INDEX_3],    myNTR_8_EEDko_2B[INDEX_4],    myNTR_8_EEDko_2B[INDEX_5] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),  rep("EEDheto_Medium", length(INDEX_3)),  rep("EEDheto_High", length(INDEX_4)),  rep("EEDheto_Highest", length(INDEX_5)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),    rep("EEDko_Medium", length(INDEX_3)),    rep("EEDko_High", length(INDEX_4)),    rep("EEDko_Highest", length(INDEX_5)) ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",     "EEDko_Medium",    "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2",    "EEDko_Highest"="purple2" ), 
                  path2=subdir_8A_part4,   fileName2= paste("4-8A-3-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=0.5)    ## width = 1 + 10*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_8_EEDheto_2A[INDEX_1], HalfLife_8_EEDheto_2A[INDEX_2],  HalfLife_8_EEDheto_2A[INDEX_3],  HalfLife_8_EEDheto_2A[INDEX_4],  HalfLife_8_EEDheto_2A[INDEX_5],
                             HalfLife_8_EEDko_2B[INDEX_1],   HalfLife_8_EEDko_2B[INDEX_2],    HalfLife_8_EEDko_2B[INDEX_3],    HalfLife_8_EEDko_2B[INDEX_4],    HalfLife_8_EEDko_2B[INDEX_5] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),  rep("EEDheto_Medium", length(INDEX_3)),  rep("EEDheto_High", length(INDEX_4)),  rep("EEDheto_Highest", length(INDEX_5)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),    rep("EEDko_Medium", length(INDEX_3)),    rep("EEDko_High", length(INDEX_4)),    rep("EEDko_Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",     "EEDko_Medium",    "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2",    "EEDko_Highest"="purple2"  ), 
                  path2=subdir_8A_part4,   fileName2= paste("4-8A-4-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=15)    ## width = 1 + 10*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( myNTR_8_EEDheto_2A[INDEX_1], myNTR_8_EEDheto_2A[INDEX_2],  myNTR_8_EEDheto_2A[INDEX_3],  myNTR_8_EEDheto_2A[INDEX_4],  myNTR_8_EEDheto_2A[INDEX_5],
                             myNTR_8_EEDko_2B[INDEX_1],   myNTR_8_EEDko_2B[INDEX_2],    myNTR_8_EEDko_2B[INDEX_3],    myNTR_8_EEDko_2B[INDEX_4],    myNTR_8_EEDko_2B[INDEX_5] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),  rep("EEDheto_Medium", length(INDEX_3)),  rep("EEDheto_High", length(INDEX_4)),  rep("EEDheto_Highest", length(INDEX_5)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),    rep("EEDko_Medium", length(INDEX_3)),    rep("EEDko_High", length(INDEX_4)),    rep("EEDko_Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "EEDheto_Lowest",   "EEDko_Lowest",   
                                 "EEDheto_Low",      "EEDko_Low", 
                                 "EEDheto_Medium",   "EEDko_Medium", 
                                 "EEDheto_High",     "EEDko_High",  
                                 "EEDheto_Highest",  "EEDko_Highest" ),     
                  colours2=c( "EEDheto_Lowest"="red2",        "EEDko_Lowest"="red2",   
                              "EEDheto_Low"="cyan2",          "EEDko_Low"="cyan2",  
                              "EEDheto_Medium"="green2",      "EEDko_Medium"="green2", 
                              "EEDheto_High"="blue2",         "EEDko_High"="blue2",  
                              "EEDheto_Highest"="purple2" ,   "EEDko_Highest"="purple2" ), 
                  path2=subdir_8A_part4,   fileName2= paste("4-8A-5-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0,   Ymax2=0.5)    ## width = 1 + 10*0.5, height=5cm 









MyBoxViolinPlot_1(vector2=c( myNTR_8_banding_3A[INDEX_1], myNTR_8_banding_3A[INDEX_2],  myNTR_8_banding_3A[INDEX_3],  myNTR_8_banding_3A[INDEX_4],  myNTR_8_banding_3A[INDEX_5],
                             myNTR_8_sham_3B[INDEX_1],    myNTR_8_sham_3B[INDEX_2],     myNTR_8_sham_3B[INDEX_3],     myNTR_8_sham_3B[INDEX_4],     myNTR_8_sham_3B[INDEX_5]  ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),  rep("banding_Medium", length(INDEX_3)),  rep("banding_High", length(INDEX_4)),  rep("banding_Highest", length(INDEX_5)), 
                                 rep("sham_Lowest", length(INDEX_1)),     rep("sham_Low", length(INDEX_2)),     rep("sham_Medium", length(INDEX_3)),     rep("sham_High", length(INDEX_4)),     rep("sham_Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",   "banding_Medium",  "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",     "sham_Low",      "sham_Medium",     "sham_High",     "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",      "sham_Low"="cyan2",       "sham_Medium"="green2",     "sham_High"="blue2",     "sham_Highest"="purple2"  ), 
                  path2=subdir_8A_part4,   fileName2= paste("4-8A-6-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=2)    ## width = 1 + 10*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_8_banding_3A[INDEX_1], HalfLife_8_banding_3A[INDEX_2],  HalfLife_8_banding_3A[INDEX_3],  HalfLife_8_banding_3A[INDEX_4],  HalfLife_8_banding_3A[INDEX_5],
                             HalfLife_8_sham_3B[INDEX_1],    HalfLife_8_sham_3B[INDEX_2],     HalfLife_8_sham_3B[INDEX_3],     HalfLife_8_sham_3B[INDEX_4],     HalfLife_8_sham_3B[INDEX_5] ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),  rep("banding_Medium", length(INDEX_3)),  rep("banding_High", length(INDEX_4)),  rep("banding_Highest", length(INDEX_5)), 
                                 rep("sham_Lowest", length(INDEX_1)),    rep("sham_Low", length(INDEX_2)),    rep("sham_Medium", length(INDEX_3)),    rep("sham_High", length(INDEX_4)),    rep("sham_Highest", length(INDEX_5))
                  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",   "banding_Medium",  "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",    "sham_Low",     "sham_Medium",    "sham_High",    "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",     "sham_Low"="cyan2",      "sham_Medium"="green2",    "sham_High"="blue2",    "sham_Highest"="purple2" ), 
                  path2=subdir_8A_part4,   fileName2= paste("4-8A-7-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=6,   Ymin2=0, Ymax2=5)    ## width = 1 + 10*0.5, height=5cm 



MyBoxViolinPlot_1(vector2=c( myNTR_8_banding_3A[INDEX_1], myNTR_8_banding_3A[INDEX_2],  myNTR_8_banding_3A[INDEX_3],  myNTR_8_banding_3A[INDEX_4],  myNTR_8_banding_3A[INDEX_5],
                             myNTR_8_sham_3B[INDEX_1],   myNTR_8_sham_3B[INDEX_2],    myNTR_8_sham_3B[INDEX_3],    myNTR_8_sham_3B[INDEX_4],    myNTR_8_sham_3B[INDEX_5] ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),  rep("banding_Medium", length(INDEX_3)),  rep("banding_High", length(INDEX_4)),  rep("banding_Highest", length(INDEX_5)), 
                                 rep("sham_Lowest", length(INDEX_1)),    rep("sham_Low", length(INDEX_2)),    rep("sham_Medium", length(INDEX_3)),    rep("sham_High", length(INDEX_4)),    rep("sham_Highest", length(INDEX_5))  ), 
                  sampleRank2=c( "banding_Lowest",   "sham_Lowest",   
                                 "banding_Low",      "sham_Low", 
                                 "banding_Medium",   "sham_Medium", 
                                 "banding_High",     "sham_High",  
                                 "banding_Highest",  "sham_Highest" ),     
                  colours2=c( "banding_Lowest"="red2",        "sham_Lowest"="red2",   
                              "banding_Low"="cyan2",          "sham_Low"="cyan2",  
                              "banding_Medium"="green2",      "sham_Medium"="green2", 
                              "banding_High"="blue2",         "sham_High"="blue2",  
                              "banding_Highest"="purple2" ,   "sham_Highest"="purple2"  ), 
                  path2=subdir_8A_part4,   fileName2= paste("4-8A-8-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=6,   Ymin2=0,   Ymax2=2)    ## width = 1 + 10*0.5, height=5cm 































##########################################################################   
subdir_7B_part4 <- paste(Part4_g,  "/8B-NTR-eachRow-500bp", sep = "")
if( ! file.exists(subdir_7B_part4) ) { dir.create(subdir_7B_part4) }


dim(Average_H3)
dim(Average_week0)
numRows_oneClass <- floor( nrow(Average_H3)/4 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
INDEX_4 <- seq(from = numRows_oneClass*3+1, to = numRows_oneClass*4, by =1 )
nrow(Average_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) + length(INDEX_4) )

MyBoxViolinPlot_1(vector2=c( myNTR_8_WT_1[INDEX_1], myNTR_8_WT_1[INDEX_2],  myNTR_8_WT_1[INDEX_3],  myNTR_8_WT_1[INDEX_4]  ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),  rep("High", length(INDEX_3)),      rep("Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "Lowest",  "Low",  "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",  "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-8B-1-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=3.0,   Ymin2=0, Ymax2=0.7)    ## width = 1 + 4*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c( HalfLife_8_WT_1[INDEX_1], HalfLife_8_WT_1[INDEX_2],  HalfLife_8_WT_1[INDEX_3],  HalfLife_8_WT_1[INDEX_4]  ),   
                  sampleType2=c( rep("Lowest", length(INDEX_1)),    rep("Low", length(INDEX_2)),  rep("High", length(INDEX_3)),      rep("Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "Lowest",  "Low",    "High",  "Highest"  ),     
                  colours2=c( "Lowest"="red2",   "Low"="cyan2",   "High"="blue2",  "Highest"="purple2"  ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-8B-2-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=3.0,   Ymin2=0, Ymax2=15)    ## width = 1 + 4*0.5, height=5cm 











MyBoxViolinPlot_1(vector2=c( myNTR_8_EEDheto_2A[INDEX_1], myNTR_8_EEDheto_2A[INDEX_2],  myNTR_8_EEDheto_2A[INDEX_3],  myNTR_8_EEDheto_2A[INDEX_4],  
                             myNTR_8_EEDko_2B[INDEX_1],   myNTR_8_EEDko_2B[INDEX_2],    myNTR_8_EEDko_2B[INDEX_3],    myNTR_8_EEDko_2B[INDEX_4] ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),       rep("EEDheto_High", length(INDEX_3)),  rep("EEDheto_Highest", length(INDEX_4)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),         rep("EEDko_High", length(INDEX_3)),    rep("EEDko_Highest", length(INDEX_4)) ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",        "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",          "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",         "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",           "EEDko_High"="blue2",    "EEDko_Highest"="purple2" ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-8B-3-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=0.5)    ## width = 1 + 8*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_8_EEDheto_2A[INDEX_1], HalfLife_8_EEDheto_2A[INDEX_2],  HalfLife_8_EEDheto_2A[INDEX_3],  HalfLife_8_EEDheto_2A[INDEX_4],   
                             HalfLife_8_EEDko_2B[INDEX_1],   HalfLife_8_EEDko_2B[INDEX_2],    HalfLife_8_EEDko_2B[INDEX_3],    HalfLife_8_EEDko_2B[INDEX_4]  ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),       rep("EEDheto_High", length(INDEX_3)),  rep("EEDheto_Highest", length(INDEX_4)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),         rep("EEDko_High", length(INDEX_3)),    rep("EEDko_Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "EEDheto_Lowest",  "EEDheto_Low",      "EEDheto_High",  "EEDheto_Highest"  ,
                                 "EEDko_Lowest",    "EEDko_Low",        "EEDko_High",    "EEDko_Highest"  ),     
                  colours2=c( "EEDheto_Lowest"="red2",   "EEDheto_Low"="cyan2",     "EEDheto_High"="blue2",  "EEDheto_Highest"="purple2" ,
                              "EEDko_Lowest"="red2",     "EEDko_Low"="cyan2",       "EEDko_High"="blue2",    "EEDko_Highest"="purple2"  ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-8B-4-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=15)    ## width = 1 + 8*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( myNTR_8_EEDheto_2A[INDEX_1], myNTR_8_EEDheto_2A[INDEX_2],  myNTR_8_EEDheto_2A[INDEX_3],  myNTR_8_EEDheto_2A[INDEX_4],  
                             myNTR_8_EEDko_2B[INDEX_1],   myNTR_8_EEDko_2B[INDEX_2],    myNTR_8_EEDko_2B[INDEX_3],    myNTR_8_EEDko_2B[INDEX_4]  ),   
                  sampleType2=c( rep("EEDheto_Lowest", length(INDEX_1)),  rep("EEDheto_Low", length(INDEX_2)),    rep("EEDheto_High", length(INDEX_3)),  rep("EEDheto_Highest", length(INDEX_4)), 
                                 rep("EEDko_Lowest", length(INDEX_1)),    rep("EEDko_Low", length(INDEX_2)),      rep("EEDko_High", length(INDEX_3)),    rep("EEDko_Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "EEDheto_Lowest",   "EEDko_Lowest",   
                                 "EEDheto_Low",      "EEDko_Low", 
                                 "EEDheto_High",     "EEDko_High",  
                                 "EEDheto_Highest",  "EEDko_Highest" ),     
                  colours2=c( "EEDheto_Lowest"="red2",        "EEDko_Lowest"="red2",   
                              "EEDheto_Low"="cyan2",          "EEDko_Low"="cyan2",  
                              "EEDheto_High"="blue2",         "EEDko_High"="blue2",  
                              "EEDheto_Highest"="purple2" ,   "EEDko_Highest"="purple2" ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-8B-5-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=5,   Ymin2=0,   Ymax2=0.5)    ## width = 1 + 8*0.5, height=5cm 









MyBoxViolinPlot_1(vector2=c( myNTR_8_banding_3A[INDEX_1], myNTR_8_banding_3A[INDEX_2],  myNTR_8_banding_3A[INDEX_3],  myNTR_8_banding_3A[INDEX_4],  
                             myNTR_8_sham_3B[INDEX_1],    myNTR_8_sham_3B[INDEX_2],     myNTR_8_sham_3B[INDEX_3],     myNTR_8_sham_3B[INDEX_4]   ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),   rep("banding_High", length(INDEX_3)),  rep("banding_Highest", length(INDEX_4)), 
                                 rep("sham_Lowest", length(INDEX_1)),     rep("sham_Low", length(INDEX_2)),      rep("sham_High", length(INDEX_3)),     rep("sham_Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",       "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",     "sham_Low",          "sham_High",     "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",        "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",      "sham_Low"="cyan2",            "sham_High"="blue2",     "sham_Highest"="purple2"  ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-8B-6-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=2)    ## width = 1 + 8*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_8_banding_3A[INDEX_1], HalfLife_8_banding_3A[INDEX_2],  HalfLife_8_banding_3A[INDEX_3],  HalfLife_8_banding_3A[INDEX_4],  
                             HalfLife_8_sham_3B[INDEX_1],    HalfLife_8_sham_3B[INDEX_2],     HalfLife_8_sham_3B[INDEX_3],     HalfLife_8_sham_3B[INDEX_4]  ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),     rep("banding_High", length(INDEX_3)),  rep("banding_Highest", length(INDEX_4)), 
                                 rep("sham_Lowest", length(INDEX_1)),    rep("sham_Low", length(INDEX_2)),         rep("sham_High", length(INDEX_3)),    rep("sham_Highest", length(INDEX_4))
                  ), 
                  sampleRank2=c( "banding_Lowest",  "banding_Low",     "banding_High",  "banding_Highest"  ,
                                 "sham_Lowest",    "sham_Low",          "sham_High",    "sham_Highest"   ),     
                  colours2=c( "banding_Lowest"="red2",   "banding_Low"="cyan2",      "banding_High"="blue2",  "banding_Highest"="purple2" ,
                              "sham_Lowest"="red2",     "sham_Low"="cyan2",           "sham_High"="blue2",    "sham_Highest"="purple2" ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-8B-7-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=5,   Ymin2=0, Ymax2=5)    ## width = 1 + 8*0.5, height=5cm 



MyBoxViolinPlot_1(vector2=c( myNTR_8_banding_3A[INDEX_1], myNTR_8_banding_3A[INDEX_2],  myNTR_8_banding_3A[INDEX_3],  myNTR_8_banding_3A[INDEX_4],   
                             myNTR_8_sham_3B[INDEX_1],   myNTR_8_sham_3B[INDEX_2],    myNTR_8_sham_3B[INDEX_3],    myNTR_8_sham_3B[INDEX_4]  ),   
                  sampleType2=c( rep("banding_Lowest", length(INDEX_1)),  rep("banding_Low", length(INDEX_2)),    rep("banding_High", length(INDEX_3)),  rep("banding_Highest", length(INDEX_4)), 
                                 rep("sham_Lowest", length(INDEX_1)),    rep("sham_Low", length(INDEX_2)),        rep("sham_High", length(INDEX_3)),    rep("sham_Highest", length(INDEX_4))  ), 
                  sampleRank2=c( "banding_Lowest",   "sham_Lowest",   
                                 "banding_Low",      "sham_Low", 
                                 "banding_High",     "sham_High",  
                                 "banding_Highest",  "sham_Highest" ),     
                  colours2=c( "banding_Lowest"="red2",        "sham_Lowest"="red2",   
                              "banding_Low"="cyan2",          "sham_Low"="cyan2",  
                              "banding_High"="blue2",         "sham_High"="blue2",  
                              "banding_Highest"="purple2" ,   "sham_Highest"="purple2"  ), 
                  path2=subdir_7B_part4,   fileName2= paste("4-8B-8-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=5,   Ymin2=0,   Ymax2=2)    ## width = 1 + 8*0.5, height=5cm 









































##########################################################################   
subdir_7C_part4 <- paste(Part4_g,  "/8C-NTR-eachRow-500bp", sep = "")
if( ! file.exists(subdir_7C_part4) ) { dir.create(subdir_7C_part4) }


dim(Average_H3)
dim(Average_week0)
numRows_oneClass <- floor( nrow(Average_H3)/3 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
nrow(Average_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) )

MyBoxViolinPlot_1(vector2=c( myNTR_8_WT_1[INDEX_1], myNTR_8_WT_1[INDEX_2],  myNTR_8_WT_1[INDEX_3]  ),   
                  sampleType2=c(  rep("Low", length(INDEX_1)),    rep("Medium", length(INDEX_2)),    rep("High", length(INDEX_3))   ), 
                  sampleRank2=c(  "Low",   "Medium",  "High"  ),     
                  colours2=c(  "Low"="cyan2",    "Medium"="green2",  "High"="blue2"   ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-8C-1-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=2.5,   Ymin2=0, Ymax2=0.7)    ## width = 1 + 3*0.5, height=5cm 

MyBoxViolinPlot_1(vector2=c( HalfLife_8_WT_1[INDEX_1], HalfLife_8_WT_1[INDEX_2],  HalfLife_8_WT_1[INDEX_3]  ),   
                  sampleType2=c(  rep("Low", length(INDEX_1)),  rep("Medium", length(INDEX_2)),   rep("High", length(INDEX_3))   ), 
                  sampleRank2=c(  "Low",   "Medium",  "High"  ),     
                  colours2=c( "Low"="cyan2",    "Medium"="green2",  "High"="blue2"   ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-8C-2-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=2.5,   Ymin2=0, Ymax2=15)    ## width = 1 + 3*0.5, height=5cm 











MyBoxViolinPlot_1(vector2=c( myNTR_8_EEDheto_2A[INDEX_1], myNTR_8_EEDheto_2A[INDEX_2],  myNTR_8_EEDheto_2A[INDEX_3],  
                             myNTR_8_EEDko_2B[INDEX_1],   myNTR_8_EEDko_2B[INDEX_2],    myNTR_8_EEDko_2B[INDEX_3]  ),   
                  sampleType2=c(   rep("EEDheto_Low", length(INDEX_1)),  rep("EEDheto_Medium", length(INDEX_2)),  rep("EEDheto_High", length(INDEX_3)),    
                                   rep("EEDko_Low", length(INDEX_1)),    rep("EEDko_Medium", length(INDEX_2)),    rep("EEDko_High", length(INDEX_3))  ), 
                  sampleRank2=c(    "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",   
                                    "EEDko_Low",     "EEDko_Medium",    "EEDko_High"   ),     
                  colours2=c(    "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",   
                                 "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2"  ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-8C-3-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=4,   Ymin2=0, Ymax2=0.5)    ## width = 1 + 6*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_8_EEDheto_2A[INDEX_1], HalfLife_8_EEDheto_2A[INDEX_2],  HalfLife_8_EEDheto_2A[INDEX_3],   
                             HalfLife_8_EEDko_2B[INDEX_1],   HalfLife_8_EEDko_2B[INDEX_2],    HalfLife_8_EEDko_2B[INDEX_3]  ),   
                  sampleType2=c(    rep("EEDheto_Low", length(INDEX_1)),  rep("EEDheto_Medium", length(INDEX_2)),  rep("EEDheto_High", length(INDEX_3)),  
                                    rep("EEDko_Low", length(INDEX_1)),    rep("EEDko_Medium", length(INDEX_2)),    rep("EEDko_High", length(INDEX_3))  ), 
                  sampleRank2=c(  "EEDheto_Low",   "EEDheto_Medium",  "EEDheto_High",  
                                  "EEDko_Low",     "EEDko_Medium",    "EEDko_High"   ),     
                  colours2=c(  "EEDheto_Low"="cyan2",    "EEDheto_Medium"="green2",  "EEDheto_High"="blue2",   
                               "EEDko_Low"="cyan2",      "EEDko_Medium"="green2",    "EEDko_High"="blue2"   ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-8C-4-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=4,   Ymin2=0, Ymax2=15)    ## width = 1 + 6*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( myNTR_8_EEDheto_2A[INDEX_1], myNTR_8_EEDheto_2A[INDEX_2],  myNTR_8_EEDheto_2A[INDEX_3],   
                             myNTR_8_EEDko_2B[INDEX_1],   myNTR_8_EEDko_2B[INDEX_2],    myNTR_8_EEDko_2B[INDEX_3]  ),   
                  sampleType2=c(    rep("EEDheto_Low", length(INDEX_1)),  rep("EEDheto_Medium", length(INDEX_2)),  rep("EEDheto_High", length(INDEX_3)),   
                                    rep("EEDko_Low", length(INDEX_1)),    rep("EEDko_Medium", length(INDEX_2)),    rep("EEDko_High", length(INDEX_3))  ), 
                  sampleRank2=c( "EEDheto_Low",      "EEDko_Low", 
                                 "EEDheto_Medium",   "EEDko_Medium", 
                                 "EEDheto_High",     "EEDko_High" ),     
                  colours2=c( "EEDheto_Low"="cyan2",          "EEDko_Low"="cyan2",  
                              "EEDheto_Medium"="green2",      "EEDko_Medium"="green2", 
                              "EEDheto_High"="blue2",         "EEDko_High"="blue2" ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-8C-5-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=4,   Ymin2=0,   Ymax2=0.5)    ## width = 1 + 6*0.5, height=5cm 









MyBoxViolinPlot_1(vector2=c( myNTR_8_banding_3A[INDEX_1], myNTR_8_banding_3A[INDEX_2],  myNTR_8_banding_3A[INDEX_3],   
                             myNTR_8_sham_3B[INDEX_1],    myNTR_8_sham_3B[INDEX_2],     myNTR_8_sham_3B[INDEX_3]    ),   
                  sampleType2=c(  rep("banding_Low", length(INDEX_1)),  rep("banding_Medium", length(INDEX_2)),  rep("banding_High", length(INDEX_3)),   
                                  rep("sham_Low", length(INDEX_1)),     rep("sham_Medium", length(INDEX_2)),     rep("sham_High", length(INDEX_3))   ), 
                  sampleRank2=c(  "banding_Low",   "banding_Medium",  "banding_High",  
                                  "sham_Low",      "sham_Medium",     "sham_High"   ),     
                  colours2=c(  "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",  
                               "sham_Low"="cyan2",       "sham_Medium"="green2",     "sham_High"="blue2"   ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-8C-6-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=4,   Ymin2=0, Ymax2=2)    ## width = 1 + 6*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( HalfLife_8_banding_3A[INDEX_1], HalfLife_8_banding_3A[INDEX_2],  HalfLife_8_banding_3A[INDEX_3],  
                             HalfLife_8_sham_3B[INDEX_1],    HalfLife_8_sham_3B[INDEX_2],     HalfLife_8_sham_3B[INDEX_3]  ),   
                  sampleType2=c(   rep("banding_Low", length(INDEX_1)),  rep("banding_Medium", length(INDEX_2)),  rep("banding_High", length(INDEX_3)),   
                                   rep("sham_Low", length(INDEX_1)),     rep("sham_Medium", length(INDEX_2)),     rep("sham_High", length(INDEX_3)) ), 
                  sampleRank2=c(   "banding_Low",   "banding_Medium",  "banding_High",   
                                   "sham_Low",     "sham_Medium",    "sham_High"    ),     
                  colours2=c(   "banding_Low"="cyan2",    "banding_Medium"="green2",  "banding_High"="blue2",   
                                "sham_Low"="cyan2",      "sham_Medium"="green2",    "sham_High"="blue2"  ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-8C-7-WT-HalfLife-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="Half Life",   
                  height2=4.0,   width2=4,   Ymin2=0, Ymax2=5)    ## width = 1 + 6*0.5, height=5cm 


MyBoxViolinPlot_1(vector2=c( myNTR_8_banding_3A[INDEX_1], myNTR_8_banding_3A[INDEX_2],  myNTR_8_banding_3A[INDEX_3],   
                             myNTR_8_sham_3B[INDEX_1],   myNTR_8_sham_3B[INDEX_2],    myNTR_8_sham_3B[INDEX_3]  ),   
                  sampleType2=c(   rep("banding_Low", length(INDEX_1)),  rep("banding_Medium", length(INDEX_2)),  rep("banding_High", length(INDEX_3)),   
                                   rep("sham_Low", length(INDEX_1)),     rep("sham_Medium", length(INDEX_2)),    rep("sham_High", length(INDEX_3))  ), 
                  sampleRank2=c( "banding_Low",      "sham_Low", 
                                 "banding_Medium",   "sham_Medium", 
                                 "banding_High",     "sham_High"  ),     
                  colours2=c( "banding_Low"="cyan2",          "sham_Low"="cyan2",  
                              "banding_Medium"="green2",      "sham_Medium"="green2", 
                              "banding_High"="blue2",         "sham_High"="blue2"   ), 
                  path2=subdir_7C_part4,   fileName2= paste("4-8C-8-WT-NTR-BoxViolin",   "1kb",  mySample_g, sep = "_") ,  
                  title2=myTitle_g,   xLab2="Samples",    yLab2="NTR",   
                  height2=4.0,   width2=4,   Ymin2=0,   Ymax2=2)    ## width = 1 + 6*0.5, height=5cm 


























          
############################### Correlations between NOL and NTR
subdir_9_part4 <- paste(Part4_g,  "/9-NTR-NOL-Correlations-5kb", sep = "")
if( ! file.exists(subdir_9_part4) ) { dir.create(subdir_9_part4) }


############################# WT 6 samples
length(row_Average_H3) 
length(row_Average_week0)    
length(myNTR_5_WT_1)
MyCor2Vars_4(vector1=row_Average_H3,    vector2=myNTR_5_WT_1, file1=paste(subdir_9_part4,  "/4-9-1A-WT-NTR-H3.txt",      sep="") )
MyCor2Vars_4(vector1=row_Average_week0, vector2=myNTR_5_WT_1, file1=paste(subdir_9_part4,  "/4-9-1B-WT-NTR-week0.txt",   sep="") )

MyCor2Vars_4(vector1=log(row_Average_H3+0.001),    vector2=myNTR_5_WT_1, file1=paste(subdir_9_part4,  "/4-9-1C-WT-NTR-H3-log.txt",      sep="") )
MyCor2Vars_4(vector1=log(row_Average_week0+0.001), vector2=myNTR_5_WT_1, file1=paste(subdir_9_part4,  "/4-9-1D-WT-NTR-week0-log.txt",   sep="") )


myNTR_withEquation_3( xAxis2=row_Average_H3,  yAxis2=myNTR_5_WT_1,                               
                      path2=subdir_9_part4,  fileName2="4-9-1A-WT-NTR-H3_1", 
                      title2=myTitle_g , xLab2="H3", yLab2="NTR", yMin2=0, yMax2=0.3,   xMin2=0, xMax2=1.4 , height2=3,  width2=4)

myNTR_withEquation_3( xAxis2= row_Average_H3 ,  yAxis2=myNTR_5_WT_1,                               
                      path2=subdir_9_part4,  fileName2="4-9-1A-WT-NTR-H3_2", 
                      title2=myTitle_g , xLab2="H3", yLab2="NTR", yMin2=0.00, yMax2=0.3,   xMin2=0.25, xMax2=1.4 , height2=3,  width2=4)


myNTR_withEquation_3( xAxis2=row_Average_week0,  yAxis2=myNTR_5_WT_1,                               
                      path2=subdir_9_part4,  fileName2="4-9-1B-WT-NTR-week0_1", 
                      title2=myTitle_g , xLab2="H3", yLab2="NTR", yMin2=0, yMax2=0.3,   xMin2=0, xMax2=1.4 , height2=3,  width2=4)

myNTR_withEquation_3( xAxis2= row_Average_week0 ,  yAxis2=myNTR_5_WT_1,                               
                      path2=subdir_9_part4,  fileName2="4-9-1B-WT-NTR-week0_2", 
                      title2=myTitle_g , xLab2="H3", yLab2="NTR", yMin2=0.00, yMax2=0.3,   xMin2=0.25, xMax2=1.4 , height2=3,  width2=4)

















############################### Correlations between NOL and NTR
subdir_10_part4 <- paste(Part4_g,  "/10-NTR-NOL-Correlations-2kb", sep = "")
if( ! file.exists(subdir_10_part4) ) { dir.create(subdir_10_part4) }


############################# WT 6 samples
length(reduceColumn2_Average_H3) 
length(reduceColumn2_Average_week0)    
length(myNTR_6_WT_1)
MyCor2Vars_4(vector1=reduceColumn2_Average_H3,    vector2=myNTR_6_WT_1, file1=paste(subdir_10_part4,  "/4-10-1A-WT-NTR-H3.txt",      sep="") )
MyCor2Vars_4(vector1=reduceColumn2_Average_week0, vector2=myNTR_6_WT_1, file1=paste(subdir_10_part4,  "/4-10-1B-WT-NTR-week0.txt",   sep="") )

MyCor2Vars_4(vector1=log(reduceColumn2_Average_H3+0.001),    vector2=myNTR_6_WT_1, file1=paste(subdir_10_part4,  "/4-10-1C-WT-NTR-H3-log.txt",      sep="") )
MyCor2Vars_4(vector1=log(reduceColumn2_Average_week0+0.001), vector2=myNTR_6_WT_1, file1=paste(subdir_10_part4,  "/4-10-1D-WT-NTR-week0-log.txt",   sep="") )



myNTR_withEquation_3( xAxis2=reduceColumn2_Average_H3,  yAxis2=myNTR_6_WT_1,                               
                      path2=subdir_10_part4,  fileName2="4-10-1A-WT-NTR-H3_1", 
                      title2=myTitle_g , xLab2="H3", yLab2="NTR", yMin2=0, yMax2=0.3,   xMin2=0, xMax2=1.4 , height2=3,  width2=4)

myNTR_withEquation_3( xAxis2= reduceColumn2_Average_H3 ,  yAxis2=myNTR_6_WT_1,                               
                      path2=subdir_10_part4,  fileName2="4-10-1A-WT-NTR-H3_2", 
                      title2=myTitle_g , xLab2="H3", yLab2="NTR", yMin2=0.00, yMax2=0.3,   xMin2=0.25, xMax2=1.4 , height2=3,  width2=4)


myNTR_withEquation_3( xAxis2=reduceColumn2_Average_week0,  yAxis2=myNTR_6_WT_1,                               
                      path2=subdir_10_part4,  fileName2="4-10-1B-WT-NTR-week0_1", 
                      title2=myTitle_g , xLab2="H3", yLab2="NTR", yMin2=0, yMax2=0.3,   xMin2=0, xMax2=1.4 , height2=3,  width2=4)

myNTR_withEquation_3( xAxis2= reduceColumn2_Average_week0 ,  yAxis2=myNTR_6_WT_1,                               
                      path2=subdir_10_part4,  fileName2="4-10-1B-WT-NTR-week0_2", 
                      title2=myTitle_g , xLab2="H3", yLab2="NTR", yMin2=0.00, yMax2=0.3,   xMin2=0.25, xMax2=1.4 , height2=3,  width2=4)




####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################




