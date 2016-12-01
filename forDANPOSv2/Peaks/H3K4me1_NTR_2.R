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
##   We should get a n*m matrix that means it contains n rows (n DNA fragment) and m columns (m bin, 1 bin=20bp).
##   The element of matrix is reads density: ChIP-input (If ChIP-input <0, we let it = 0)
#############################################################################################################################


############################################################################
## Part  2:  Read in the input matrix data, and reduce the dimensions (rows or columns) of matrix, and save the information of input files.
############################################################################





#################################################################### Start ##########################################################################################################################################
H3_Rep1     <- read.table("H3K4me1-WT/1_H3_Rep1_center_heatmap/H3K4me1-WT.1_H3_Rep1.wig.heatmap.xls",                        header=TRUE,   sep="",   quote = "",   comment.char = "")  
week0_Rep1  <- read.table("H3K4me1-WT/2_H2bGFP-week0_Rep1_center_heatmap/H3K4me1-WT.2_H2bGFP-week0_Rep1.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week0_Rep2  <- read.table("H3K4me1-WT/2_H2bGFP-week0_Rep2_center_heatmap/H3K4me1-WT.2_H2bGFP-week0_Rep2.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week1_Rep1  <- read.table("H3K4me1-WT/1_H2bGFP-week1_Rep1_center_heatmap/H3K4me1-WT.1_H2bGFP-week1_Rep1.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week1_Rep2  <- read.table("H3K4me1-WT/1_H2bGFP-week1_Rep2_center_heatmap/H3K4me1-WT.1_H2bGFP-week1_Rep2.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week2_Rep1  <- read.table("H3K4me1-WT/1_H2bGFP-week2_Rep1_center_heatmap/H3K4me1-WT.1_H2bGFP-week2_Rep1.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week2_Rep2  <- read.table("H3K4me1-WT/1_H2bGFP-week2_Rep2_center_heatmap/H3K4me1-WT.1_H2bGFP-week2_Rep2.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week4_Rep1  <- read.table("H3K4me1-WT/1_H2bGFP-week4_Rep1_center_heatmap/H3K4me1-WT.1_H2bGFP-week4_Rep1.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week4_Rep2  <- read.table("H3K4me1-WT/2_H2BGFP-CM-week4-2015_Rep1_center_heatmap/H3K4me1-WT.2_H2BGFP-CM-week4-2015_Rep1.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week4_Rep3  <- read.table("H3K4me1-WT/2_H2BGFP-CM-week4-2015_Rep2_center_heatmap/H3K4me1-WT.2_H2BGFP-CM-week4-2015_Rep2.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week6_Rep1  <- read.table("H3K4me1-WT/1_H2bGFP-week6_Rep1_center_heatmap/H3K4me1-WT.1_H2bGFP-week6_Rep1.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week6_Rep2  <- read.table("H3K4me1-WT/1_H2bGFP-week6_Rep2_center_heatmap/H3K4me1-WT.1_H2bGFP-week6_Rep2.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week8_Rep1  <- read.table("H3K4me1-WT/1_H2bGFP-week8_Rep1_center_heatmap/H3K4me1-WT.1_H2bGFP-week8_Rep1.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week0_EEDheto_Rep1  <- read.table("H3K4me1-WT/3_H2BGFP-CM-week0-EEDheto_Rep1_center_heatmap/H3K4me1-WT.3_H2BGFP-CM-week0-EEDheto_Rep1.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week0_EEDheto_Rep2  <- read.table("H3K4me1-WT/3_H2BGFP-CM-week0-EEDheto_Rep2_center_heatmap/H3K4me1-WT.3_H2BGFP-CM-week0-EEDheto_Rep2.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week0_EEDko_Rep1    <- read.table("H3K4me1-WT/3_H2BGFP-CM-week0-EEDko_Rep1_center_heatmap/H3K4me1-WT.3_H2BGFP-CM-week0-EEDko_Rep1.wig.heatmap.xls",        header=TRUE,   sep="",   quote = "",   comment.char = "")  
week4_EEDheto_Rep1  <- read.table("H3K4me1-WT/4_H2BGFP-CM-week4-EEDheto_Rep1_center_heatmap/H3K4me1-WT.4_H2BGFP-CM-week4-EEDheto_Rep1.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week4_EEDheto_Rep2  <- read.table("H3K4me1-WT/4_H2BGFP-CM-week4-EEDheto_Rep2_center_heatmap/H3K4me1-WT.4_H2BGFP-CM-week4-EEDheto_Rep2.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
week4_EEDko_Rep1    <- read.table("H3K4me1-WT/4_H2BGFP-CM-week4-EEDko_Rep1_center_heatmap/H3K4me1-WT.4_H2BGFP-CM-week4-EEDko_Rep1.wig.heatmap.xls",        header=TRUE,   sep="",   quote = "",   comment.char = "")  
week4_EEDko_Rep2    <- read.table("H3K4me1-WT/4_H2BGFP-CM-week4-EEDko_Rep2_center_heatmap/H3K4me1-WT.4_H2BGFP-CM-week4-EEDko_Rep2.wig.heatmap.xls",        header=TRUE,   sep="",   quote = "",   comment.char = "")  
banding_Rep1  <- read.table("H3K4me1-WT/5_H2BGFP-CM-week2.5-banding_Rep1_center_heatmap/H3K4me1-WT.5_H2BGFP-CM-week2.5-banding_Rep1.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
banding_Rep2  <- read.table("H3K4me1-WT/5_H2BGFP-CM-week2.5-banding_Rep2_center_heatmap/H3K4me1-WT.5_H2BGFP-CM-week2.5-banding_Rep2.wig.heatmap.xls",    header=TRUE,   sep="",   quote = "",   comment.char = "")  
sham_Rep1     <- read.table("H3K4me1-WT/5_H2BGFP-CM-week2.5-sham_Rep1_center_heatmap/H3K4me1-WT.5_H2BGFP-CM-week2.5-sham_Rep1.wig.heatmap.xls",          header=TRUE,   sep="",   quote = "",   comment.char = "")  
sham_Rep2     <- read.table("H3K4me1-WT/5_H2BGFP-CM-week2.5-sham_Rep2_center_heatmap/H3K4me1-WT.5_H2BGFP-CM-week2.5-sham_Rep2.wig.heatmap.xls",          header=TRUE,   sep="",   quote = "",   comment.char = "")  

dim(H3_Rep1)
dim(week0_Rep1)
dim(week0_Rep2) 
dim(week1_Rep1) 
dim(week1_Rep2) 
dim(week2_Rep1) 
dim(week2_Rep2)
dim(week4_Rep1)
dim(week4_Rep2) 
dim(week4_Rep3) 
dim(week6_Rep1)
dim(week6_Rep2) 
dim(week8_Rep1) 
dim(week0_EEDheto_Rep1)  
dim(week0_EEDheto_Rep2)  
dim(week0_EEDko_Rep1)  
dim(week4_EEDheto_Rep1)  
dim(week4_EEDheto_Rep2)  
dim(week4_EEDko_Rep1)  
dim(week4_EEDko_Rep2)  
dim(banding_Rep1)  
dim(banding_Rep2)  
dim(sham_Rep1)  
dim(sham_Rep2) 


identical(H3_Rep1[, 1],  H3_Rep1[, 1] )
identical(H3_Rep1[, 1],  week0_Rep1[, 1] )
identical(H3_Rep1[, 1],  week0_Rep2[, 1] )
identical(H3_Rep1[, 1],  week1_Rep1[, 1] )
identical(H3_Rep1[, 1],  week1_Rep2[, 1] )
identical(H3_Rep1[, 1],  week2_Rep1[, 1] )
identical(H3_Rep1[, 1],  week2_Rep2[, 1] )
identical(H3_Rep1[, 1],  week4_Rep1[, 1] )
identical(H3_Rep1[, 1],  week4_Rep2[, 1] )
identical(H3_Rep1[, 1],  week4_Rep3[, 1] )
identical(H3_Rep1[, 1],  week6_Rep1[, 1] )
identical(H3_Rep1[, 1],  week6_Rep2[, 1] )
identical(H3_Rep1[, 1],  week8_Rep1[, 1] )
identical(H3_Rep1[, 1],  week0_EEDheto_Rep1[, 1] )
identical(H3_Rep1[, 1],  week0_EEDheto_Rep2[, 1] )
identical(H3_Rep1[, 1],  week0_EEDko_Rep1[, 1] )
identical(H3_Rep1[, 1],  week4_EEDheto_Rep1[, 1] )
identical(H3_Rep1[, 1],  week4_EEDheto_Rep2[, 1] )
identical(H3_Rep1[, 1],  week4_EEDko_Rep1[, 1] )
identical(H3_Rep1[, 1],  week4_EEDko_Rep2[, 1] )
identical(H3_Rep1[, 1],  banding_Rep1[, 1] )
identical(H3_Rep1[, 1],  banding_Rep2[, 1] )
identical(H3_Rep1[, 1],  sham_Rep1[, 1] )
identical(H3_Rep1[, 1],  sham_Rep2[, 1] )

thisRowNames1 <- H3_Rep1[, 1]
thisRowNames2 <- week8_Rep1[, 1]
thisRowNames3 <- banding_Rep2[, 1]
length(thisRowNames1)
length(thisRowNames2)
length(thisRowNames3)
thisRowNames1[100]
thisRowNames2[100]
thisRowNames3[100]

write.table(x=thisRowNames1,  file=paste(Part2_g, "/Part2-1-A-rowNames.txt", sep="") , quote = FALSE, sep="\t")                                                        
write.table(x=thisRowNames2,  file=paste(Part2_g, "/Part2-1-B-rowNames.txt", sep="") , quote = FALSE, sep="\t")                                                        
write.table(x=thisRowNames3,  file=paste(Part2_g, "/Part2-1-C-rowNames.txt", sep="") , quote = FALSE, sep="\t")                                                        

					
			

			
######################################################### remove 1st-4th columns of each matrix
H3_Rep1    <- as.matrix(H3_Rep1[, -c(1:4)]) 
week0_Rep1 <- as.matrix(week0_Rep1[, -c(1:4)]) * 1.08
week0_Rep2 <- as.matrix(week0_Rep2[, -c(1:4)]) * 1.08
week1_Rep1 <- as.matrix(week1_Rep1[, -c(1:4)]) 
week1_Rep2 <- as.matrix(week1_Rep2[, -c(1:4)]) 
week2_Rep1 <- as.matrix(week2_Rep1[, -c(1:4)]) 
week2_Rep2 <- as.matrix(week2_Rep2[, -c(1:4)]) 
week4_Rep1 <- as.matrix(week4_Rep1[, -c(1:4)]) 
week4_Rep2 <- as.matrix(week4_Rep2[, -c(1:4)])
week4_Rep3 <- as.matrix(week4_Rep3[, -c(1:4)]) 
week6_Rep1 <- as.matrix(week6_Rep1[, -c(1:4)]) 
week6_Rep2 <- as.matrix(week6_Rep2[, -c(1:4)]) 
week8_Rep1 <- as.matrix(week8_Rep1[, -c(1:4)]) 
week0_EEDheto_Rep1 <- as.matrix(week0_EEDheto_Rep1[, -c(1:4)]) 
week0_EEDheto_Rep2 <- as.matrix(week0_EEDheto_Rep2[, -c(1:4)]) 
week0_EEDko_Rep1   <- as.matrix(week0_EEDko_Rep1[, -c(1:4)]) 
week4_EEDheto_Rep1 <- as.matrix(week4_EEDheto_Rep1[, -c(1:4)]) 
week4_EEDheto_Rep2 <- as.matrix(week4_EEDheto_Rep2[, -c(1:4)]) 
week4_EEDko_Rep1   <- as.matrix(week4_EEDko_Rep1[, -c(1:4)]) 
week4_EEDko_Rep2   <- as.matrix(week4_EEDko_Rep2[, -c(1:4)]) 
banding_Rep1 <- as.matrix(banding_Rep1[, -c(1:4)]) 
banding_Rep2 <- as.matrix(banding_Rep2[, -c(1:4)]) 
sham_Rep1    <- as.matrix(sham_Rep1[, -c(1:4)]) 
sham_Rep2    <- as.matrix(sham_Rep2[, -c(1:4)]) 

sink( file=paste(Part2_g,  "/Part2-2-A-dimersion.txt",      sep = "") )
dim(H3_Rep1)
dim(week0_Rep1)
dim(week0_Rep2) 
dim(week1_Rep1) 
dim(week1_Rep2) 
dim(week2_Rep1) 
dim(week2_Rep2)
dim(week4_Rep1)
dim(week4_Rep2) 
dim(week4_Rep3) 
dim(week6_Rep1)
dim(week6_Rep2) 
dim(week8_Rep1) 
dim(week0_EEDheto_Rep1)  
dim(week0_EEDheto_Rep2)  
dim(week0_EEDko_Rep1)  
dim(week4_EEDheto_Rep1)  
dim(week4_EEDheto_Rep2)  
dim(week4_EEDko_Rep1)  
dim(week4_EEDko_Rep2)  
dim(banding_Rep1)  
dim(banding_Rep2)  
dim(sham_Rep1)  
dim(sham_Rep2)  
sink()

sink( file=paste(Part2_g,  "/Part2-2-B-summaryForRawMatrix.txt",      sep = "") )
summary( H3_Rep1)  
summary( week0_Rep1)
summary( week0_Rep2) 
summary( week1_Rep1) 
summary( week1_Rep2) 
summary( week2_Rep1) 
summary( week2_Rep2)
summary( week4_Rep1)
summary( week4_Rep2) 
summary( week4_Rep3) 
summary( week6_Rep1)
summary( week6_Rep2) 
summary( week8_Rep1) 
summary( week0_EEDheto_Rep1)  
summary( week0_EEDheto_Rep2)  
summary( week0_EEDko_Rep1)  
summary( week4_EEDheto_Rep1)  
summary( week4_EEDheto_Rep2)  
summary( week4_EEDko_Rep1)  
summary( week4_EEDko_Rep2)  
summary( banding_Rep1)  
summary( banding_Rep2)  
summary( sham_Rep1)  
summary( sham_Rep2)  
sink() 

sink( file=paste(Part2_g,  "/Part2-2-C-summaryForRawMatrix-asVector.txt",      sep = "") )
summary( as.vector(H3_Rep1))
summary( as.vector(week0_Rep1))
summary( as.vector(week0_Rep2)) 
summary( as.vector(week1_Rep1)) 
summary( as.vector(week1_Rep2)) 
summary( as.vector(week2_Rep1)) 
summary( as.vector(week2_Rep2))
summary( as.vector(week4_Rep1))
summary( as.vector(week4_Rep2)) 
summary( as.vector(week4_Rep3)) 
summary( as.vector(week6_Rep1))
summary( as.vector(week6_Rep2)) 
summary( as.vector(week8_Rep1)) 
summary( as.vector(week0_EEDheto_Rep1))  
summary( as.vector(week0_EEDheto_Rep2))  
summary( as.vector(week0_EEDko_Rep1))  
summary( as.vector(week4_EEDheto_Rep1))  
summary( as.vector(week4_EEDheto_Rep2))  
summary( as.vector(week4_EEDko_Rep1))  
summary( as.vector(week4_EEDko_Rep2))  
summary( as.vector(banding_Rep1))  
summary( as.vector(banding_Rep2))  
summary( as.vector(sham_Rep1))  
summary( as.vector(sham_Rep2))  
sink() 





################################################## average biological replicates
Average_H3     <- H3_Rep1
Average_week0  <- (week0_Rep1+week0_Rep2)/2
Average_week1  <- (week1_Rep1+week1_Rep2)/2 
Average_week2  <- (week2_Rep1+week2_Rep2)/2 
Average_week4  <- (week4_Rep2+week4_Rep3)/2
Average_week6  <- week6_Rep2
Average_week8  <- week8_Rep1 
Average_week0_EEDheto  <- (week0_EEDheto_Rep1+week0_EEDheto_Rep2)/2  
Average_week0_EEDko    <- week0_EEDko_Rep1  
Average_week4_EEDheto  <- week4_EEDheto_Rep2 
Average_week4_EEDko    <- (week4_EEDko_Rep1+week4_EEDko_Rep2)/2  
Average_banding  <- (banding_Rep1+banding_Rep2)/2  
Average_sham     <- (sham_Rep1+sham_Rep2)/2  

sink( file=paste(Part2_g,  "/Part2-3-A-dimersion-AverageBioRep.txt",      sep = "") )
dim(Average_H3)
dim(Average_week0)
dim(Average_week1)
dim(Average_week2)
dim(Average_week4)
dim(Average_week6)
dim(Average_week8)
dim(Average_week0_EEDheto)  
dim(Average_week0_EEDko) 
dim(Average_week4_EEDheto)
dim(Average_week4_EEDko) 
dim(Average_banding) 
dim(Average_sham)
sink()   

sink( file=paste(Part2_g,  "/Part2-3-B-summaryAverage-asMatrix.txt",      sep = "") )
summary( Average_H3)
summary( Average_week0)
summary( Average_week1)
summary( Average_week2)
summary( Average_week4)
summary( Average_week6)
summary( Average_week8)
summary( Average_week0_EEDheto)  
summary( Average_week0_EEDko) 
summary( Average_week4_EEDheto)
summary( Average_week4_EEDko) 
summary( Average_banding) 
summary( Average_sham)
sink() 

sink( file=paste(Part2_g,  "/Part2-3-C-summaryAverage-asVector.txt",      sep = "") )
summary( as.vector(Average_H3))
summary( as.vector(Average_week0))
summary( as.vector(Average_week1))
summary( as.vector(Average_week2))
summary( as.vector(Average_week4))
summary( as.vector(Average_week6))
summary( as.vector(Average_week8))
summary( as.vector(Average_week0_EEDheto))  
summary( as.vector(Average_week0_EEDko)) 
summary( as.vector(Average_week4_EEDheto))
summary( as.vector(Average_week4_EEDko)) 
summary( as.vector(Average_banding)) 
summary( as.vector(Average_sham))
sink() 





#####################################################################  average every row , 计算每一行的平均
row_H3_Rep1     <- rowMeans(H3_Rep1)
row_week0_Rep1  <- rowMeans(week0_Rep1)
row_week0_Rep2  <- rowMeans(week0_Rep2) 
row_week1_Rep1  <- rowMeans(week1_Rep1) 
row_week1_Rep2  <- rowMeans(week1_Rep2) 
row_week2_Rep1  <- rowMeans(week2_Rep1) 
row_week2_Rep2  <- rowMeans(week2_Rep2)
row_week4_Rep1  <- rowMeans(week4_Rep1)
row_week4_Rep2  <- rowMeans(week4_Rep2) 
row_week4_Rep3  <- rowMeans(week4_Rep3) 
row_week6_Rep1  <- rowMeans(week6_Rep1)
row_week6_Rep2  <- rowMeans(week6_Rep2) 
row_week8_Rep1  <- rowMeans(week8_Rep1) 
row_week0_EEDheto_Rep1  <- rowMeans(week0_EEDheto_Rep1)  
row_week0_EEDheto_Rep2  <- rowMeans(week0_EEDheto_Rep2)  
row_week0_EEDko_Rep1    <- rowMeans(week0_EEDko_Rep1)  
row_week4_EEDheto_Rep1  <- rowMeans(week4_EEDheto_Rep1)  
row_week4_EEDheto_Rep2  <- rowMeans(week4_EEDheto_Rep2)  
row_week4_EEDko_Rep1    <- rowMeans(week4_EEDko_Rep1)  
row_week4_EEDko_Rep2    <- rowMeans(week4_EEDko_Rep2)  
row_banding_Rep1  <- rowMeans(banding_Rep1)  
row_banding_Rep2  <- rowMeans(banding_Rep2)  
row_sham_Rep1     <- rowMeans(sham_Rep1)  
row_sham_Rep2     <- rowMeans(sham_Rep2)  

length(row_H3_Rep1)
length(row_week0_Rep1)
length(row_week0_Rep2)
length(row_week1_Rep1) 
length(row_week1_Rep2) 
length(row_week2_Rep1)
length(row_week2_Rep2)
length(row_week4_Rep1)
length(row_week4_Rep2) 
length(row_week4_Rep3)
length(row_week6_Rep1)
length(row_week6_Rep2)
length(row_week8_Rep1)
length(row_week0_EEDheto_Rep1)  
length(row_week0_EEDheto_Rep2)  
length(row_week0_EEDko_Rep1) 
length(row_week4_EEDheto_Rep1)
length(row_week4_EEDheto_Rep2) 
length(row_week4_EEDko_Rep1) 
length(row_week4_EEDko_Rep2) 
length(row_banding_Rep1) 
length(row_banding_Rep2)
length(row_sham_Rep1)  
length(row_sham_Rep2) 

rowAverage_matrix  <- cbind(row_H3_Rep1, row_week0_Rep1, row_week0_Rep2,  row_week1_Rep1,  row_week1_Rep2, row_week2_Rep1, row_week2_Rep2, row_week4_Rep1, row_week4_Rep2, row_week4_Rep3, 
                            row_week6_Rep1, row_week6_Rep2, row_week8_Rep1, row_week0_EEDheto_Rep1, row_week0_EEDheto_Rep2, row_week0_EEDko_Rep1, row_week4_EEDheto_Rep1, 
                            row_week4_EEDheto_Rep2, row_week4_EEDko_Rep1, row_week4_EEDko_Rep2, row_banding_Rep1, row_banding_Rep2, row_sham_Rep1, row_sham_Rep2) 
colnames(rowAverage_matrix) <- cbind("row_H3_Rep1", "row_week0_Rep1", "row_week0_Rep2",  "row_week1_Rep1",  "row_week1_Rep2", "row_week2_Rep1", "row_week2_Rep2", "row_week4_Rep1", "row_week4_Rep2", "row_week4_Rep3", 
                                     "row_week6_Rep1", "row_week6_Rep2", "row_week8_Rep1", "row_week0_EEDheto_Rep1", "row_week0_EEDheto_Rep2", "row_week0_EEDko_Rep1", "row_week4_EEDheto_Rep1", 
                                     "row_week4_EEDheto_Rep2", "row_week4_EEDko_Rep1", "row_week4_EEDko_Rep2", "row_banding_Rep1", "row_banding_Rep2", "row_sham_Rep1", "row_sham_Rep2") 
rownames(rowAverage_matrix) <- thisRowNames1
write.table(x=rowAverage_matrix,  file=paste(Part2_g, "/Part2-4-A-rowAverage-withReps.txt", sep="") , quote = FALSE, sep="\t")     

sink( file=paste(Part2_g,  "/Part2-4-B-summary-averageRow-withReps.txt",      sep = "") )

cat("\n\n\n")
print("row_H3_Rep1:")
sd(row_H3_Rep1)
mean(row_H3_Rep1)
summary( as.vector(row_H3_Rep1))

cat("\n\n\n")
print("row_week0_Rep1:")
sd(row_week0_Rep1)
mean(row_week0_Rep1)
summary( as.vector(row_week0_Rep1))

cat("\n\n\n")
print("row_week0_Rep2:")
sd(row_week0_Rep2)
mean(row_week0_Rep2)
summary( as.vector(row_week0_Rep2))

cat("\n\n\n")
print("row_week1_Rep1:")
sd(row_week1_Rep1)
mean(row_week1_Rep1)
summary( as.vector(row_week1_Rep1)) 

cat("\n\n\n")
print("row_week1_Rep2:")
sd(row_week1_Rep2)
mean(row_week1_Rep2)
summary( as.vector(row_week1_Rep2)) 

cat("\n\n\n")
print("row_week2_Rep1:")
sd(row_week2_Rep1)
mean(row_week2_Rep1)
summary( as.vector(row_week2_Rep1))

cat("\n\n\n")
print("row_week2_Rep2:")
sd(row_week2_Rep2)
mean(row_week2_Rep2)
summary( as.vector(row_week2_Rep2))

cat("\n\n\n")
print("row_week4_Rep1:")
sd(row_week4_Rep1)
mean(row_week4_Rep1)
summary( as.vector(row_week4_Rep1))

cat("\n\n\n")
print("row_week4_Rep2:")
sd(row_week4_Rep2)
mean(row_week4_Rep2)
summary( as.vector(row_week4_Rep2)) 

cat("\n\n\n")
print("row_week4_Rep3:")
sd(row_week4_Rep3)
mean(row_week4_Rep3)
summary( as.vector(row_week4_Rep3))

cat("\n\n\n")
print("row_week6_Rep1:")
sd(row_week6_Rep1)
mean(row_week6_Rep1)
summary( as.vector(row_week6_Rep1))

cat("\n\n\n")
print("row_week6_Rep2:")
sd(row_week6_Rep2)
mean(row_week6_Rep2)
summary( as.vector(row_week6_Rep2))

cat("\n\n\n")
print("row_week8_Rep1:")
sd(row_week8_Rep1)
mean(row_week8_Rep1)
summary( as.vector(row_week8_Rep1))

cat("\n\n\n")
print("row_week0_EEDheto_Rep1:")
sd(row_week0_EEDheto_Rep1)
mean(row_week0_EEDheto_Rep1)
summary( as.vector(row_week0_EEDheto_Rep1))  

cat("\n\n\n")
print("row_week0_EEDheto_Rep2:")
sd(row_week0_EEDheto_Rep2)
mean(row_week0_EEDheto_Rep2)
summary( as.vector(row_week0_EEDheto_Rep2))  

cat("\n\n\n")
print("row_week0_EEDko_Rep1:")
sd(row_week0_EEDko_Rep1)
mean(row_week0_EEDko_Rep1)
summary( as.vector(row_week0_EEDko_Rep1)) 

cat("\n\n\n")
print("row_week4_EEDheto_Rep1:")
sd(row_week4_EEDheto_Rep1)
mean(row_week4_EEDheto_Rep1)
summary( as.vector(row_week4_EEDheto_Rep1))

cat("\n\n\n")
print("row_week4_EEDheto_Rep2:")
sd(row_week4_EEDheto_Rep2)
mean(row_week4_EEDheto_Rep2)
summary( as.vector(row_week4_EEDheto_Rep2)) 

cat("\n\n\n")
print("row_week4_EEDko_Rep1:")
sd(row_week4_EEDko_Rep1)
mean(row_week4_EEDko_Rep1)
summary( as.vector(row_week4_EEDko_Rep1)) 

cat("\n\n\n")
print("row_week4_EEDko_Rep2:")
sd(row_week4_EEDko_Rep2)
mean(row_week4_EEDko_Rep2)
summary( as.vector(row_week4_EEDko_Rep2)) 

cat("\n\n\n")
print("row_banding_Rep1:")
sd(row_banding_Rep1)
mean(row_banding_Rep1)
summary( as.vector(row_banding_Rep1)) 

cat("\n\n\n")
print("row_banding_Rep2:")
sd(row_banding_Rep2)
mean(row_banding_Rep2)
summary( as.vector(row_banding_Rep2))

cat("\n\n\n")
print("row_sham_Rep1:")
sd(row_sham_Rep1)
mean(row_sham_Rep1)
summary( as.vector(row_sham_Rep1))  

cat("\n\n\n")
print("row_sham_Rep2:")
sd(row_sham_Rep2)
mean(row_sham_Rep2)
summary( as.vector(row_sham_Rep2)) 

sink()


row_Average_H3     <- rowMeans(Average_H3)
row_Average_week0  <- rowMeans(Average_week0)
row_Average_week1  <- rowMeans(Average_week1)
row_Average_week2  <- rowMeans(Average_week2)
row_Average_week4  <- rowMeans(Average_week4)
row_Average_week6  <- rowMeans(Average_week6)
row_Average_week8  <- rowMeans(Average_week8)
row_Average_week0_EEDheto  <- rowMeans(Average_week0_EEDheto)  
row_Average_week0_EEDko    <- rowMeans(Average_week0_EEDko) 
row_Average_week4_EEDheto  <- rowMeans(Average_week4_EEDheto)
row_Average_week4_EEDko    <- rowMeans(Average_week4_EEDko) 
row_Average_banding  <- rowMeans(Average_banding) 
row_Average_sham     <- rowMeans(Average_sham)

length(row_Average_H3)
length(row_Average_week0)
length(row_Average_week1)
length(row_Average_week2)
length(row_Average_week4)
length(row_Average_week6)
length(row_Average_week8)
length(row_Average_week0_EEDheto)  
length(row_Average_week0_EEDko)
length(row_Average_week4_EEDheto )
length(row_Average_week4_EEDko)
length(row_Average_banding)
length(row_Average_sham)

rowAverage_matrix2  <- cbind(row_Average_H3, row_Average_week0,  row_Average_week1, row_Average_week2, row_Average_week4, row_Average_week6, row_Average_week8, 
                            row_Average_week0_EEDheto, row_Average_week0_EEDko, row_Average_week4_EEDheto, row_Average_week4_EEDko, row_Average_banding, row_Average_sham ) 
colnames(rowAverage_matrix2) <- cbind("H3", "week0",  "week1", "week2", "week4", "week6", "week8", "week0_EEDheto", "week0_EEDko", "week4_EEDheto", "week4_EEDko", "banding", "sham" ) 
rownames(rowAverage_matrix2) <- thisRowNames1
write.table(x=rowAverage_matrix2,  file=paste(Part2_g, "/Part2-4-C-rowAverage.txt", sep="") , quote = FALSE, sep="\t")     

sink( file=paste(Part2_g,  "/Part2-4-D-summary-averageRow.txt",      sep = "") )

cat("\n\n\n")
print("row_Average_H3:")
sd(row_Average_H3)
mean(row_Average_H3)
summary( as.vector(row_Average_H3))

cat("\n\n\n")
print("row_Average_week0:")
sd(row_Average_week0)
mean(row_Average_week0)
summary( as.vector(row_Average_week0))


cat("\n\n\n")
print("row_Average_week1:")
sd(row_Average_week1)
mean(row_Average_week1)
summary( as.vector(row_Average_week1)) 


cat("\n\n\n")
print("row_Average_week2:")
sd(row_Average_week2)
mean(row_Average_week2)
summary( as.vector(row_Average_week2))


cat("\n\n\n")
print("row_Average_week4:")
sd(row_Average_week4)
mean(row_Average_week4)
summary( as.vector(row_Average_week4))

cat("\n\n\n")
print("row_Average_week6:")
sd(row_Average_week6)
mean(row_Average_week6)
summary( as.vector(row_Average_week6))

cat("\n\n\n")
print("row_Average_week8:")
sd(row_Average_week8)
mean(row_Average_week8)
summary( as.vector(row_Average_week8))

cat("\n\n\n")
print("row_Average_week0_EEDheto:")
sd(row_Average_week0_EEDheto)
mean(row_Average_week0_EEDheto)
summary( as.vector(row_Average_week0_EEDheto))  


cat("\n\n\n")
print("row_Average_week0_EEDko:")
sd(row_Average_week0_EEDko)
mean(row_Average_week0_EEDko)
summary( as.vector(row_Average_week0_EEDko)) 

cat("\n\n\n")
print("row_Average_week4_EEDheto:")
sd(row_Average_week4_EEDheto)
mean(row_Average_week4_EEDheto)
summary( as.vector(row_Average_week4_EEDheto))

cat("\n\n\n")
print("row_Average_week4_EEDko:")
sd(row_Average_week4_EEDko)
mean(row_Average_week4_EEDko)
summary( as.vector(row_Average_week4_EEDko)) 


cat("\n\n\n")
print("row_Average_banding:")
sd(row_Average_banding)
mean(row_Average_banding)
summary( as.vector(row_Average_banding)) 


cat("\n\n\n")
print("row_Average_sham:")
sd(row_Average_sham)
mean(row_Average_sham)
summary( as.vector(row_Average_sham))  

sink()





#####################################################   average every column,  计算每一列的平均 
column_H3_Rep1     <- colMeans(H3_Rep1)
column_week0_Rep1  <- colMeans(week0_Rep1)
column_week0_Rep2  <- colMeans(week0_Rep2) 
column_week1_Rep1  <- colMeans(week1_Rep1) 
column_week1_Rep2  <- colMeans(week1_Rep2) 
column_week2_Rep1  <- colMeans(week2_Rep1) 
column_week2_Rep2  <- colMeans(week2_Rep2)
column_week4_Rep1  <- colMeans(week4_Rep1)
column_week4_Rep2  <- colMeans(week4_Rep2) 
column_week4_Rep3  <- colMeans(week4_Rep3) 
column_week6_Rep1  <- colMeans(week6_Rep1)
column_week6_Rep2  <- colMeans(week6_Rep2) 
column_week8_Rep1  <- colMeans(week8_Rep1) 
column_week0_EEDheto_Rep1  <- colMeans(week0_EEDheto_Rep1)  
column_week0_EEDheto_Rep2  <- colMeans(week0_EEDheto_Rep2)  
column_week0_EEDko_Rep1    <- colMeans(week0_EEDko_Rep1)  
column_week4_EEDheto_Rep1  <- colMeans(week4_EEDheto_Rep1)  
column_week4_EEDheto_Rep2  <- colMeans(week4_EEDheto_Rep2)  
column_week4_EEDko_Rep1    <- colMeans(week4_EEDko_Rep1)  
column_week4_EEDko_Rep2    <- colMeans(week4_EEDko_Rep2)  
column_banding_Rep1  <- colMeans(banding_Rep1)  
column_banding_Rep2  <- colMeans(banding_Rep2)  
column_sham_Rep1     <- colMeans(sham_Rep1)  
column_sham_Rep2     <- colMeans(sham_Rep2)  

length(column_H3_Rep1)
length(column_week0_Rep1)
length(column_week0_Rep2)
length(column_week1_Rep1) 
length(column_week1_Rep2) 
length(column_week2_Rep1)
length(column_week2_Rep2)
length(column_week4_Rep1)
length(column_week4_Rep2) 
length(column_week4_Rep3)
length(column_week6_Rep1)
length(column_week6_Rep2)
length(column_week8_Rep1)
length(column_week0_EEDheto_Rep1)  
length(column_week0_EEDheto_Rep2)  
length(column_week0_EEDko_Rep1) 
length(column_week4_EEDheto_Rep1)
length(column_week4_EEDheto_Rep2) 
length(column_week4_EEDko_Rep1) 
length(column_week4_EEDko_Rep2) 
length(column_banding_Rep1) 
length(column_banding_Rep2)
length(column_sham_Rep1)  
length(column_sham_Rep2) 

columnAverage_matrix  <- rbind(column_H3_Rep1, column_week0_Rep1, column_week0_Rep2,  column_week1_Rep1,  column_week1_Rep2, column_week2_Rep1, column_week2_Rep2, column_week4_Rep1, column_week4_Rep2, column_week4_Rep3, 
                            column_week6_Rep1, column_week6_Rep2, column_week8_Rep1, column_week0_EEDheto_Rep1, column_week0_EEDheto_Rep2, column_week0_EEDko_Rep1, column_week4_EEDheto_Rep1, 
                            column_week4_EEDheto_Rep2, column_week4_EEDko_Rep1, column_week4_EEDko_Rep2, column_banding_Rep1, column_banding_Rep2, column_sham_Rep1, column_sham_Rep2) 
rownames(columnAverage_matrix) <- cbind("column_H3_Rep1", "column_week0_Rep1", "column_week0_Rep2",  "column_week1_Rep1",  "column_week1_Rep2", "column_week2_Rep1", "column_week2_Rep2", "column_week4_Rep1", "column_week4_Rep2", "column_week4_Rep3", 
                                     "column_week6_Rep1", "column_week6_Rep2", "column_week8_Rep1", "column_week0_EEDheto_Rep1", "column_week0_EEDheto_Rep2", "column_week0_EEDko_Rep1", "column_week4_EEDheto_Rep1", 
                                     "column_week4_EEDheto_Rep2", "column_week4_EEDko_Rep1", "column_week4_EEDko_Rep2", "column_banding_Rep1", "column_banding_Rep2", "column_sham_Rep1", "column_sham_Rep2") 
write.table(x=columnAverage_matrix,  file=paste(Part2_g, "/Part2-5-A-colAverage-withReps.txt", sep="") , quote = FALSE, sep="\t")     

sink( file=paste(Part2_g,  "/Part2-5-B-summary-averageCol-withReps.txt",      sep = "") )

cat("\n\n\n")
print("column_H3_Rep1:")
sd(column_H3_Rep1)
mean(column_H3_Rep1)
summary( as.vector(column_H3_Rep1))

cat("\n\n\n")
print("column_week0_Rep1:")
sd(column_week0_Rep1)
mean(column_week0_Rep1)
summary( as.vector(column_week0_Rep1))

cat("\n\n\n")
print("column_week0_Rep2:")
sd(column_week0_Rep2)
mean(column_week0_Rep2)
summary( as.vector(column_week0_Rep2))

cat("\n\n\n")
print("column_week1_Rep1:")
sd(column_week1_Rep1)
mean(column_week1_Rep1)
summary( as.vector(column_week1_Rep1)) 

cat("\n\n\n")
print("column_week1_Rep2:")
sd(column_week1_Rep2)
mean(column_week1_Rep2)
summary( as.vector(column_week1_Rep2)) 

cat("\n\n\n")
print("column_week2_Rep1:")
sd(column_week2_Rep1)
mean(column_week2_Rep1)
summary( as.vector(column_week2_Rep1))

cat("\n\n\n")
print("column_week2_Rep2:")
sd(column_week2_Rep2)
mean(column_week2_Rep2)
summary( as.vector(column_week2_Rep2))

cat("\n\n\n")
print("column_week4_Rep1:")
sd(column_week4_Rep1)
mean(column_week4_Rep1)
summary( as.vector(column_week4_Rep1))

cat("\n\n\n")
print("column_week4_Rep2:")
sd(column_week4_Rep2)
mean(column_week4_Rep2)
summary( as.vector(column_week4_Rep2)) 

cat("\n\n\n")
print("column_week4_Rep3:")
sd(column_week4_Rep3)
mean(column_week4_Rep3)
summary( as.vector(column_week4_Rep3))

cat("\n\n\n")
print("column_week6_Rep1:")
sd(column_week6_Rep1)
mean(column_week6_Rep1)
summary( as.vector(column_week6_Rep1))

cat("\n\n\n")
print("column_week6_Rep2:")
sd(column_week6_Rep2)
mean(column_week6_Rep2)
summary( as.vector(column_week6_Rep2))

cat("\n\n\n")
print("column_week8_Rep1:")
sd(column_week8_Rep1)
mean(column_week8_Rep1)
summary( as.vector(column_week8_Rep1))

cat("\n\n\n")
print("column_week0_EEDheto_Rep1:")
sd(column_week0_EEDheto_Rep1)
mean(column_week0_EEDheto_Rep1)
summary( as.vector(column_week0_EEDheto_Rep1))  

cat("\n\n\n")
print("column_week0_EEDheto_Rep2:")
sd(column_week0_EEDheto_Rep2)
mean(column_week0_EEDheto_Rep2)
summary( as.vector(column_week0_EEDheto_Rep2))  

cat("\n\n\n")
print("column_week0_EEDko_Rep1:")
sd(column_week0_EEDko_Rep1)
mean(column_week0_EEDko_Rep1)
summary( as.vector(column_week0_EEDko_Rep1)) 

cat("\n\n\n")
print("column_week4_EEDheto_Rep1:")
sd(column_week4_EEDheto_Rep1)
mean(column_week4_EEDheto_Rep1)
summary( as.vector(column_week4_EEDheto_Rep1))

cat("\n\n\n")
print("column_week4_EEDheto_Rep2:")
sd(column_week4_EEDheto_Rep2)
mean(column_week4_EEDheto_Rep2)
summary( as.vector(column_week4_EEDheto_Rep2)) 

cat("\n\n\n")
print("column_week4_EEDko_Rep1:")
sd(column_week4_EEDko_Rep1)
mean(column_week4_EEDko_Rep1)
summary( as.vector(column_week4_EEDko_Rep1)) 

cat("\n\n\n")
print("column_week4_EEDko_Rep2:")
sd(column_week4_EEDko_Rep2)
mean(column_week4_EEDko_Rep2)
summary( as.vector(column_week4_EEDko_Rep2)) 

cat("\n\n\n")
print("column_banding_Rep1:")
sd(column_banding_Rep1)
mean(column_banding_Rep1)
summary( as.vector(column_banding_Rep1)) 

cat("\n\n\n")
print("column_banding_Rep2:")
sd(column_banding_Rep2)
mean(column_banding_Rep2)
summary( as.vector(column_banding_Rep2))

cat("\n\n\n")
print("column_sham_Rep1:")
sd(column_sham_Rep1)
mean(column_sham_Rep1)
summary( as.vector(column_sham_Rep1))  

cat("\n\n\n")
print("column_sham_Rep2:")
sd(column_sham_Rep2)
mean(column_sham_Rep2)
summary( as.vector(column_sham_Rep2)) 

sink()


column_Average_H3     <- colMeans(Average_H3)
column_Average_week0  <- colMeans(Average_week0)
column_Average_week1  <- colMeans(Average_week1)
column_Average_week2  <- colMeans(Average_week2)
column_Average_week4  <- colMeans(Average_week4)
column_Average_week6  <- colMeans(Average_week6)
column_Average_week8  <- colMeans(Average_week8)
column_Average_week0_EEDheto  <- colMeans(Average_week0_EEDheto)  
column_Average_week0_EEDko    <- colMeans(Average_week0_EEDko) 
column_Average_week4_EEDheto  <- colMeans(Average_week4_EEDheto)
column_Average_week4_EEDko    <- colMeans(Average_week4_EEDko) 
column_Average_banding  <- colMeans(Average_banding) 
column_Average_sham     <- colMeans(Average_sham)

length(column_Average_H3)
length(column_Average_week0)
length(column_Average_week1)
length(column_Average_week2)
length(column_Average_week4)
length(column_Average_week6)
length(column_Average_week8)
length(column_Average_week0_EEDheto)  
length(column_Average_week0_EEDko)
length(column_Average_week4_EEDheto )
length(column_Average_week4_EEDko)
length(column_Average_banding)
length(column_Average_sham)

columnAverage_matrix2  <- rbind(column_Average_H3, column_Average_week0,  column_Average_week1, column_Average_week2, column_Average_week4, column_Average_week6, column_Average_week8, 
                             column_Average_week0_EEDheto, column_Average_week0_EEDko, column_Average_week4_EEDheto, column_Average_week4_EEDko, column_Average_banding, column_Average_sham ) 
rownames(columnAverage_matrix2) <- cbind("H3", "week0",  "week1", "week2", "week4", "week6", "week8",  "week0_EEDheto", "week0_EEDko", "week4_EEDheto", "week4_EEDko", "banding", "sham" ) 
write.table(x=columnAverage_matrix2,  file=paste(Part2_g, "/Part2-5-C-columnAverage.txt", sep="") , quote = FALSE, sep="\t")     

sink( file=paste(Part2_g,  "/Part2-5-D-summary-averageCol.txt",      sep = "") )

cat("\n\n\n")
print("column_Average_H3:")
sd(column_Average_H3)
mean(column_Average_H3)
summary( as.vector(column_Average_H3))

cat("\n\n\n")
print("column_Average_week0:")
sd(column_Average_week0)
mean(column_Average_week0)
summary( as.vector(column_Average_week0))

cat("\n\n\n")
print("column_Average_week1:")
sd(column_Average_week1)
mean(column_Average_week1)
summary( as.vector(column_Average_week1)) 

cat("\n\n\n")
print("column_Average_week2:")
sd(column_Average_week2)
mean(column_Average_week2)
summary( as.vector(column_Average_week2))

cat("\n\n\n")
print("column_Average_week4:")
sd(column_Average_week4)
mean(column_Average_week4)
summary( as.vector(column_Average_week4))

cat("\n\n\n")
print("column_Average_week6:")
sd(column_Average_week6)
mean(column_Average_week6)
summary( as.vector(column_Average_week6))

cat("\n\n\n")
print("column_Average_week8:")
sd(column_Average_week8)
mean(column_Average_week8)
summary( as.vector(column_Average_week8))

cat("\n\n\n")
print("column_Average_week0_EEDheto:")
sd(column_Average_week0_EEDheto)
mean(column_Average_week0_EEDheto)
summary( as.vector(column_Average_week0_EEDheto))  

cat("\n\n\n")
print("column_Average_week0_EEDko:")
sd(column_Average_week0_EEDko)
mean(column_Average_week0_EEDko)
summary( as.vector(column_Average_week0_EEDko)) 

cat("\n\n\n")
print("column_Average_week4_EEDheto:")
sd(column_Average_week4_EEDheto)
mean(column_Average_week4_EEDheto)
summary( as.vector(column_Average_week4_EEDheto))

cat("\n\n\n")
print("column_Average_week4_EEDko:")
sd(column_Average_week4_EEDko)
mean(column_Average_week4_EEDko)
summary( as.vector(column_Average_week4_EEDko)) 

cat("\n\n\n")
print("column_Average_banding:")
sd(column_Average_banding)
mean(column_Average_banding)
summary( as.vector(column_Average_banding)) 

cat("\n\n\n")
print("column_Average_sham:")
sd(column_Average_sham)
mean(column_Average_sham)
summary( as.vector(column_Average_sham))  

sink()










##########  number of rows and columns  
numOfColumns1 <- ncol(Average_H3)
numOfRows1    <- nrow(Average_H3)  
numOfColumns1
numOfRows1

sink(paste(Part2_g,  "/Part2-6-A-numberOf-columns-rows.txt",           sep = ""))
cat("number Of Columns:", numOfColumns1, "\n")
cat("number Of Rows:",    numOfRows1, "\n\n\n")
sink() 










##########    The standard error of the mean (SEM) 
SEM_Average_H3 <- apply(Average_H3,    2, sd)/sqrt( nrow(Average_H3)-1 )  
SEM_Average_week0 <- apply(Average_week0,    2, sd)/sqrt( nrow(Average_week0)-1 )  
SEM_Average_week1 <- apply(Average_week1,    2, sd)/sqrt( nrow(Average_week1)-1 )  
SEM_Average_week2 <- apply(Average_week2,    2, sd)/sqrt( nrow(Average_week2)-1 )  
SEM_Average_week4 <- apply(Average_week4,    2, sd)/sqrt( nrow(Average_week4)-1 )  
SEM_Average_week6 <- apply(Average_week6,    2, sd)/sqrt( nrow(Average_week6)-1 )  
SEM_Average_week8 <- apply(Average_week8,    2, sd)/sqrt( nrow(Average_week8)-1 )  
SEM_Average_week0_EEDheto <- apply(Average_week0_EEDheto,    2, sd)/sqrt( nrow(Average_week0_EEDheto)-1 )
SEM_Average_week0_EEDko   <- apply(Average_week0_EEDko,    2, sd)/sqrt( nrow(Average_week0_EEDko)-1 )
SEM_Average_week4_EEDheto <- apply(Average_week4_EEDheto,    2, sd)/sqrt( nrow(Average_week4_EEDheto)-1 )
SEM_Average_week4_EEDko   <- apply(Average_week4_EEDko,    2, sd)/sqrt( nrow(Average_week4_EEDko)-1 )
SEM_Average_banding <- apply(Average_banding,    2, sd)/sqrt( nrow(Average_banding)-1 )
SEM_Average_sham    <- apply(Average_sham,    2, sd)/sqrt( nrow(Average_sham)-1 )

#NonZero_one1_1    <- seq(from = 0,  by=20,  length.out=numOfColumns1/20)
#NonZero_one1_1[1] <- 1

#SEM_Average_H3[-NonZero_one1_1]      <- NA 
#SEM_Average_week0[-NonZero_one1_1]   <- NA   
#SEM_Average_week1[-NonZero_one1_1]   <- NA  
#SEM_Average_week2[-NonZero_one1_1]   <- NA   
#SEM_Average_week4[-NonZero_one1_1]   <- NA   
#SEM_Average_week6[-NonZero_one1_1]   <- NA   
#SEM_Average_week8[-NonZero_one1_1]   <- NA   
#SEM_Average_week0_EEDheto[-NonZero_one1_1]   <- NA   
#SEM_Average_week0_EEDko[-NonZero_one1_1]     <- NA   
#SEM_Average_week4_EEDheto[-NonZero_one1_1]   <- NA   
#SEM_Average_week4_EEDko[-NonZero_one1_1]     <- NA   
#SEM_Average_banding[-NonZero_one1_1]   <- NA   
#SEM_Average_sham[-NonZero_one1_1]      <- NA   

length(SEM_Average_H3)










##########   reduce rows 
reduceRow1_Average_H3     <- reduceMatrixRow(Average_H3,       rowNum_1=5 ) 
reduceRow1_Average_week0  <- reduceMatrixRow(Average_week0,    rowNum_1=5 ) 
reduceRow1_Average_week1  <- reduceMatrixRow(Average_week1,    rowNum_1=5 ) 
reduceRow1_Average_week2  <- reduceMatrixRow(Average_week2,    rowNum_1=5 ) 
reduceRow1_Average_week4  <- reduceMatrixRow(Average_week4,    rowNum_1=5 ) 
reduceRow1_Average_week6  <- reduceMatrixRow(Average_week6,    rowNum_1=5 ) 
reduceRow1_Average_week8  <- reduceMatrixRow(Average_week8,    rowNum_1=5 ) 
reduceRow1_Average_week0_EEDheto  <- reduceMatrixRow(Average_week0_EEDheto,    rowNum_1=5 )   
reduceRow1_Average_week0_EEDko    <- reduceMatrixRow(Average_week0_EEDko,    rowNum_1=5 )  
reduceRow1_Average_week4_EEDheto  <- reduceMatrixRow(Average_week4_EEDheto,    rowNum_1=5 ) 
reduceRow1_Average_week4_EEDko    <- reduceMatrixRow(Average_week4_EEDko,    rowNum_1=5 )  
reduceRow1_Average_banding  <- reduceMatrixRow(Average_banding,    rowNum_1=5 )  
reduceRow1_Average_sham     <- reduceMatrixRow(Average_sham,    rowNum_1=5 ) 
dim(reduceRow1_Average_H3) 
dim(reduceRow1_Average_week0) 
dim(reduceRow1_Average_week1) 
dim(reduceRow1_Average_week2) 
dim(reduceRow1_Average_week4) 
dim(reduceRow1_Average_week6) 
dim(reduceRow1_Average_week8) 
dim(reduceRow1_Average_week0_EEDheto)   
dim(reduceRow1_Average_week0_EEDko)  
dim(reduceRow1_Average_week4_EEDheto) 
dim(reduceRow1_Average_week4_EEDko)  
dim(reduceRow1_Average_banding)  
dim(reduceRow1_Average_sham) 



reduceRow2_Average_H3     <- reduceMatrixRow(Average_H3,       rowNum_1=4 ) 
reduceRow2_Average_week0  <- reduceMatrixRow(Average_week0,    rowNum_1=4 ) 
reduceRow2_Average_week1  <- reduceMatrixRow(Average_week1,    rowNum_1=4 ) 
reduceRow2_Average_week2  <- reduceMatrixRow(Average_week2,    rowNum_1=4 ) 
reduceRow2_Average_week4  <- reduceMatrixRow(Average_week4,    rowNum_1=4 ) 
reduceRow2_Average_week6  <- reduceMatrixRow(Average_week6,    rowNum_1=4 ) 
reduceRow2_Average_week8  <- reduceMatrixRow(Average_week8,    rowNum_1=4 ) 
reduceRow2_Average_week0_EEDheto  <- reduceMatrixRow(Average_week0_EEDheto,    rowNum_1=4 )   
reduceRow2_Average_week0_EEDko    <- reduceMatrixRow(Average_week0_EEDko,    rowNum_1=4 )  
reduceRow2_Average_week4_EEDheto  <- reduceMatrixRow(Average_week4_EEDheto,    rowNum_1=4 ) 
reduceRow2_Average_week4_EEDko    <- reduceMatrixRow(Average_week4_EEDko,    rowNum_1=4 )  
reduceRow2_Average_banding  <- reduceMatrixRow(Average_banding,    rowNum_1=4 )  
reduceRow2_Average_sham     <- reduceMatrixRow(Average_sham,    rowNum_1=4 ) 
dim(reduceRow2_Average_H3) 
dim(reduceRow2_Average_week0) 
dim(reduceRow2_Average_week1) 
dim(reduceRow2_Average_week2) 
dim(reduceRow2_Average_week4) 
dim(reduceRow2_Average_week6) 
dim(reduceRow2_Average_week8) 
dim(reduceRow2_Average_week0_EEDheto)   
dim(reduceRow2_Average_week0_EEDko)  
dim(reduceRow2_Average_week4_EEDheto) 
dim(reduceRow2_Average_week4_EEDko)  
dim(reduceRow2_Average_banding)  
dim(reduceRow2_Average_sham) 



reduceRow3_Average_H3     <- reduceMatrixRow(Average_H3,       rowNum_1=3 ) 
reduceRow3_Average_week0  <- reduceMatrixRow(Average_week0,    rowNum_1=3 ) 
reduceRow3_Average_week1  <- reduceMatrixRow(Average_week1,    rowNum_1=3 ) 
reduceRow3_Average_week2  <- reduceMatrixRow(Average_week2,    rowNum_1=3 ) 
reduceRow3_Average_week4  <- reduceMatrixRow(Average_week4,    rowNum_1=3 ) 
reduceRow3_Average_week6  <- reduceMatrixRow(Average_week6,    rowNum_1=3 ) 
reduceRow3_Average_week8  <- reduceMatrixRow(Average_week8,    rowNum_1=3 ) 
reduceRow3_Average_week0_EEDheto  <- reduceMatrixRow(Average_week0_EEDheto,    rowNum_1=3 )   
reduceRow3_Average_week0_EEDko    <- reduceMatrixRow(Average_week0_EEDko,    rowNum_1=3 )  
reduceRow3_Average_week4_EEDheto  <- reduceMatrixRow(Average_week4_EEDheto,    rowNum_1=3 ) 
reduceRow3_Average_week4_EEDko    <- reduceMatrixRow(Average_week4_EEDko,    rowNum_1=3 )  
reduceRow3_Average_banding  <- reduceMatrixRow(Average_banding,    rowNum_1=3 )  
reduceRow3_Average_sham     <- reduceMatrixRow(Average_sham,    rowNum_1=3 ) 
dim(reduceRow3_Average_H3) 
dim(reduceRow3_Average_week0) 
dim(reduceRow3_Average_week1) 
dim(reduceRow3_Average_week2) 
dim(reduceRow3_Average_week4) 
dim(reduceRow3_Average_week6) 
dim(reduceRow3_Average_week8) 
dim(reduceRow3_Average_week0_EEDheto)   
dim(reduceRow3_Average_week0_EEDko)  
dim(reduceRow3_Average_week4_EEDheto) 
dim(reduceRow3_Average_week4_EEDko)  
dim(reduceRow3_Average_banding)  
dim(reduceRow3_Average_sham) 










####################################################################################################
##########   reduce columns:
reduceColumn1_Average_H3     <- reduceMatrixCol(Average_H3,       colNum_1=5 ) 
reduceColumn1_Average_week0  <- reduceMatrixCol(Average_week0,    colNum_1=5 ) 
reduceColumn1_Average_week1  <- reduceMatrixCol(Average_week1,    colNum_1=5 ) 
reduceColumn1_Average_week2  <- reduceMatrixCol(Average_week2,    colNum_1=5 ) 
reduceColumn1_Average_week4  <- reduceMatrixCol(Average_week4,    colNum_1=5 ) 
reduceColumn1_Average_week6  <- reduceMatrixCol(Average_week6,    colNum_1=5 ) 
reduceColumn1_Average_week8  <- reduceMatrixCol(Average_week8,    colNum_1=5 ) 
reduceColumn1_Average_week0_EEDheto  <- reduceMatrixCol(Average_week0_EEDheto,    colNum_1=5 )   
reduceColumn1_Average_week0_EEDko    <- reduceMatrixCol(Average_week0_EEDko,    colNum_1=5 )  
reduceColumn1_Average_week4_EEDheto  <- reduceMatrixCol(Average_week4_EEDheto,    colNum_1=5 ) 
reduceColumn1_Average_week4_EEDko    <- reduceMatrixCol(Average_week4_EEDko,    colNum_1=5 )  
reduceColumn1_Average_banding  <- reduceMatrixCol(Average_banding,    colNum_1=5 )  
reduceColumn1_Average_sham     <- reduceMatrixCol(Average_sham,    colNum_1=5 ) 
dim(reduceColumn1_Average_H3) 
dim(reduceColumn1_Average_week0) 
dim(reduceColumn1_Average_week1) 
dim(reduceColumn1_Average_week2) 
dim(reduceColumn1_Average_week4) 
dim(reduceColumn1_Average_week6) 
dim(reduceColumn1_Average_week8) 
dim(reduceColumn1_Average_week0_EEDheto)   
dim(reduceColumn1_Average_week0_EEDko)  
dim(reduceColumn1_Average_week4_EEDheto) 
dim(reduceColumn1_Average_week4_EEDko)  
dim(reduceColumn1_Average_banding)  
dim(reduceColumn1_Average_sham) 





## only one column
speColumns1 <- c(150:350)   ## center ± 2k
reduceColumn2_Average_H3     <- rowMeans(Average_H3[, speColumns1])
reduceColumn2_Average_week0  <- rowMeans(Average_week0[, speColumns1])
reduceColumn2_Average_week1  <- rowMeans(Average_week1[, speColumns1])
reduceColumn2_Average_week2  <- rowMeans(Average_week2[, speColumns1])
reduceColumn2_Average_week4  <- rowMeans(Average_week4[, speColumns1])
reduceColumn2_Average_week6  <- rowMeans(Average_week6[, speColumns1])
reduceColumn2_Average_week8  <- rowMeans(Average_week8[, speColumns1])
reduceColumn2_Average_week0_EEDheto  <- rowMeans(Average_week0_EEDheto[, speColumns1])  
reduceColumn2_Average_week0_EEDko    <- rowMeans(Average_week0_EEDko[, speColumns1]) 
reduceColumn2_Average_week4_EEDheto  <- rowMeans(Average_week4_EEDheto[, speColumns1])
reduceColumn2_Average_week4_EEDko    <- rowMeans(Average_week4_EEDko[, speColumns1]) 
reduceColumn2_Average_banding  <- rowMeans(Average_banding[, speColumns1]) 
reduceColumn2_Average_sham     <- rowMeans(Average_sham[, speColumns1])
length(reduceColumn2_Average_H3)
length(reduceColumn2_Average_week0)
length(reduceColumn2_Average_week1)
length(reduceColumn2_Average_week2)
length(reduceColumn2_Average_week4)
length(reduceColumn2_Average_week6)
length(reduceColumn2_Average_week8)
length(reduceColumn2_Average_week0_EEDheto)  
length(reduceColumn2_Average_week0_EEDko)
length(reduceColumn2_Average_week4_EEDheto )
length(reduceColumn2_Average_week4_EEDko)
length(reduceColumn2_Average_banding)
length(reduceColumn2_Average_sham)

write.table(x=reduceColumn2_Average_H3,    file = paste(Part2_g,  "Part2-7-A-H3.txt", sep = "/"),    append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn2_Average_week0, file = paste(Part2_g,  "Part2-7-B-week0.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn2_Average_week1, file = paste(Part2_g,  "Part2-7-C-week1.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn2_Average_week2, file = paste(Part2_g,  "Part2-7-D-week2.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn2_Average_week4, file = paste(Part2_g,  "Part2-7-E-week4.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn2_Average_week6, file = paste(Part2_g,  "Part2-7-F-week6.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn2_Average_week8, file = paste(Part2_g,  "Part2-7-G-week8.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn2_Average_week0_EEDheto, file = paste(Part2_g,  "Part2-7-H-week0_EEDheto.txt", sep = "/"),   append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn2_Average_week0_EEDko,   file = paste(Part2_g,  "Part2-7-I-week0_EEDko.txt", sep = "/"),     append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn2_Average_week4_EEDheto, file = paste(Part2_g,  "Part2-7-J-week4_EEDheto.txt", sep = "/"),   append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn2_Average_week4_EEDko,   file = paste(Part2_g,  "Part2-7-K-week4_EEDko.txt", sep = "/"),     append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn2_Average_banding,   file = paste(Part2_g,  "Part2-7-L-banding.txt", sep = "/"),  append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn2_Average_sham,      file = paste(Part2_g,  "Part2-7-M-sham.txt", sep = "/"),     append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        

MyMergeMatrix_1 <- cbind(reduceColumn2_Average_H3, reduceColumn2_Average_week0, reduceColumn2_Average_week1, reduceColumn2_Average_week2, reduceColumn2_Average_week4, 
                         reduceColumn2_Average_week6, reduceColumn2_Average_week8, reduceColumn2_Average_week0_EEDheto, reduceColumn2_Average_week0_EEDko, 
                         reduceColumn2_Average_week4_EEDheto, reduceColumn2_Average_week4_EEDko, reduceColumn2_Average_banding, reduceColumn2_Average_sham)
rownames(MyMergeMatrix_1) <- thisRowNames1
colnames(MyMergeMatrix_1) <- c("H3", "week0", "week1", "week2", "week4",  "week6", "week8", "week0_EEDheto", "week0_EEDko",  "week4_EEDheto", "week4_EEDko", "banding", "sham")  

write.table(x=MyMergeMatrix_1,      file = paste(Part2_g,  "Part2-7-N-all-merge.txt", sep = "/"),     append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")                                                        





speColumns2 <- c(200:300)   ## center ± 1kb
reduceColumn3_Average_H3     <- rowMeans(Average_H3[, speColumns2])
reduceColumn3_Average_week0  <- rowMeans(Average_week0[, speColumns2])
reduceColumn3_Average_week1  <- rowMeans(Average_week1[, speColumns2])
reduceColumn3_Average_week2  <- rowMeans(Average_week2[, speColumns2])
reduceColumn3_Average_week4  <- rowMeans(Average_week4[, speColumns2])
reduceColumn3_Average_week6  <- rowMeans(Average_week6[, speColumns2])
reduceColumn3_Average_week8  <- rowMeans(Average_week8[, speColumns2])
reduceColumn3_Average_week0_EEDheto  <- rowMeans(Average_week0_EEDheto[, speColumns2])  
reduceColumn3_Average_week0_EEDko    <- rowMeans(Average_week0_EEDko[, speColumns2]) 
reduceColumn3_Average_week4_EEDheto  <- rowMeans(Average_week4_EEDheto[, speColumns2])
reduceColumn3_Average_week4_EEDko    <- rowMeans(Average_week4_EEDko[, speColumns2]) 
reduceColumn3_Average_banding  <- rowMeans(Average_banding[, speColumns2]) 
reduceColumn3_Average_sham     <- rowMeans(Average_sham[, speColumns2])
length(reduceColumn3_Average_H3)
length(reduceColumn3_Average_week0)
length(reduceColumn3_Average_week1)
length(reduceColumn3_Average_week2)
length(reduceColumn3_Average_week4)
length(reduceColumn3_Average_week6)
length(reduceColumn3_Average_week8)
length(reduceColumn3_Average_week0_EEDheto)  
length(reduceColumn3_Average_week0_EEDko)
length(reduceColumn3_Average_week4_EEDheto )
length(reduceColumn3_Average_week4_EEDko)
length(reduceColumn3_Average_banding)
length(reduceColumn3_Average_sham)

write.table(x=reduceColumn3_Average_H3,    file = paste(Part2_g,  "Part2-8-A-H3.txt", sep = "/"),    append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn3_Average_week0, file = paste(Part2_g,  "Part2-8-B-week0.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn3_Average_week1, file = paste(Part2_g,  "Part2-8-C-week1.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn3_Average_week2, file = paste(Part2_g,  "Part2-8-D-week2.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn3_Average_week4, file = paste(Part2_g,  "Part2-8-E-week4.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn3_Average_week6, file = paste(Part2_g,  "Part2-8-F-week6.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn3_Average_week8, file = paste(Part2_g,  "Part2-8-G-week8.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn3_Average_week0_EEDheto, file = paste(Part2_g,  "Part2-8-H-week0_EEDheto.txt", sep = "/"),   append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn3_Average_week0_EEDko,   file = paste(Part2_g,  "Part2-8-I-week0_EEDko.txt", sep = "/"),     append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn3_Average_week4_EEDheto, file = paste(Part2_g,  "Part2-8-J-week4_EEDheto.txt", sep = "/"),   append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn3_Average_week4_EEDko,   file = paste(Part2_g,  "Part2-8-K-week4_EEDko.txt", sep = "/"),     append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn3_Average_banding,   file = paste(Part2_g,  "Part2-8-L-banding.txt", sep = "/"),  append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn3_Average_sham,      file = paste(Part2_g,  "Part2-8-M-sham.txt", sep = "/"),     append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        

MyMergeMatrix_2 <- cbind(reduceColumn3_Average_H3, reduceColumn3_Average_week0, reduceColumn3_Average_week1, reduceColumn3_Average_week2, reduceColumn3_Average_week4, 
                         reduceColumn3_Average_week6, reduceColumn3_Average_week8, reduceColumn3_Average_week0_EEDheto, reduceColumn3_Average_week0_EEDko, 
                         reduceColumn3_Average_week4_EEDheto, reduceColumn3_Average_week4_EEDko, reduceColumn3_Average_banding, reduceColumn3_Average_sham)
rownames(MyMergeMatrix_2) <- thisRowNames1
colnames(MyMergeMatrix_2) <- c("H3", "week0", "week1", "week2", "week4",  "week6", "week8", "week0_EEDheto", "week0_EEDko",  "week4_EEDheto", "week4_EEDko", "banding", "sham")  

write.table(x=MyMergeMatrix_2,      file = paste(Part2_g,  "Part2-8-N-all-merge.txt", sep = "/"),     append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")                                                        





speColumns3 <- c(225:275)   ## center ± 500bp
reduceColumn4_Average_H3     <- rowMeans(Average_H3[, speColumns3])
reduceColumn4_Average_week0  <- rowMeans(Average_week0[, speColumns3])
reduceColumn4_Average_week1  <- rowMeans(Average_week1[, speColumns3])
reduceColumn4_Average_week2  <- rowMeans(Average_week2[, speColumns3])
reduceColumn4_Average_week4  <- rowMeans(Average_week4[, speColumns3])
reduceColumn4_Average_week6  <- rowMeans(Average_week6[, speColumns3])
reduceColumn4_Average_week8  <- rowMeans(Average_week8[, speColumns3])
reduceColumn4_Average_week0_EEDheto  <- rowMeans(Average_week0_EEDheto[, speColumns3])  
reduceColumn4_Average_week0_EEDko    <- rowMeans(Average_week0_EEDko[, speColumns3]) 
reduceColumn4_Average_week4_EEDheto  <- rowMeans(Average_week4_EEDheto[, speColumns3])
reduceColumn4_Average_week4_EEDko    <- rowMeans(Average_week4_EEDko[, speColumns3]) 
reduceColumn4_Average_banding  <- rowMeans(Average_banding[, speColumns3]) 
reduceColumn4_Average_sham     <- rowMeans(Average_sham[, speColumns3])
length(reduceColumn4_Average_H3)
length(reduceColumn4_Average_week0)
length(reduceColumn4_Average_week1)
length(reduceColumn4_Average_week2)
length(reduceColumn4_Average_week4)
length(reduceColumn4_Average_week6)
length(reduceColumn4_Average_week8)
length(reduceColumn4_Average_week0_EEDheto)  
length(reduceColumn4_Average_week0_EEDko)
length(reduceColumn4_Average_week4_EEDheto )
length(reduceColumn4_Average_week4_EEDko)
length(reduceColumn4_Average_banding)
length(reduceColumn4_Average_sham)

write.table(x=reduceColumn4_Average_H3,    file = paste(Part2_g,  "Part2-9-A-H3.txt", sep = "/"),    append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn4_Average_week0, file = paste(Part2_g,  "Part2-9-B-week0.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn4_Average_week1, file = paste(Part2_g,  "Part2-9-C-week1.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn4_Average_week2, file = paste(Part2_g,  "Part2-9-D-week2.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn4_Average_week4, file = paste(Part2_g,  "Part2-9-E-week4.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn4_Average_week6, file = paste(Part2_g,  "Part2-9-F-week6.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn4_Average_week8, file = paste(Part2_g,  "Part2-9-G-week8.txt", sep = "/"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn4_Average_week0_EEDheto, file = paste(Part2_g,  "Part2-9-H-week0_EEDheto.txt", sep = "/"),   append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn4_Average_week0_EEDko,   file = paste(Part2_g,  "Part2-9-I-week0_EEDko.txt", sep = "/"),     append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn4_Average_week4_EEDheto, file = paste(Part2_g,  "Part2-9-J-week4_EEDheto.txt", sep = "/"),   append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn4_Average_week4_EEDko,   file = paste(Part2_g,  "Part2-9-K-week4_EEDko.txt", sep = "/"),     append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn4_Average_banding,   file = paste(Part2_g,  "Part2-9-L-banding.txt", sep = "/"),  append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
write.table(x=reduceColumn4_Average_sham,      file = paste(Part2_g,  "Part2-9-M-sham.txt", sep = "/"),     append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        

MyMergeMatrix_3 <- cbind(reduceColumn4_Average_H3, reduceColumn4_Average_week0, reduceColumn4_Average_week1, reduceColumn4_Average_week2, reduceColumn4_Average_week4, 
                         reduceColumn4_Average_week6, reduceColumn4_Average_week8, reduceColumn4_Average_week0_EEDheto, reduceColumn4_Average_week0_EEDko, 
                         reduceColumn4_Average_week4_EEDheto, reduceColumn4_Average_week4_EEDko, reduceColumn4_Average_banding, reduceColumn4_Average_sham)
rownames(MyMergeMatrix_3) <- thisRowNames1
colnames(MyMergeMatrix_3) <- c("H3", "week0", "week1", "week2", "week4",  "week6", "week8", "week0_EEDheto", "week0_EEDko",  "week4_EEDheto", "week4_EEDko", "banding", "sham")  

write.table(x=MyMergeMatrix_3,      file = paste(Part2_g,  "Part2-9-N-all-merge.txt", sep = "/"),     append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")                                                        





save.image( file = "H3K4me1-WT_Center.RData"  )  
#################################################################### End ##########################################################################################################################################




