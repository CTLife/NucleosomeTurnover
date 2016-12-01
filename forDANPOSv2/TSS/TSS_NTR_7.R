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
## Part  7:  Figures about nucleosome contribution value (NCV).  
############################################################################


library(scatterplot3d)


subdir_1_part7 <- paste(Part7_g,  "/1-NTR-NOL-TPM", sep = "")
if( ! file.exists(subdir_1_part7) ) { dir.create(subdir_1_part7) }


length(myNTR_7_WT_1)
dim(targetRegions)  
length(reduceColumn3_Average_H3)      
length(reduceColumn3_Average_week0)

TPM_WT_7 <- targetRegions[, 22]
TPM_WT_7[24410:24413]

MyCor2Vars_2(vector1=row_Average_H3,     vector2=TPM_WT_7,  file1=paste(subdir_1_part7,  "/Part7-1-1A-TPM-H3.txt",      sep="") )
MyCor2Vars_2(vector1=row_Average_week0,  vector2=TPM_WT_7,  file1=paste(subdir_1_part7,  "/Part7-1-1B-TPM-week0.txt",   sep="") )
MyCor2Vars_2(vector1=myNTR_7_WT_1,       vector2=TPM_WT_7,  file1=paste(subdir_1_part7,  "/Part7-1-1C-TPM-NTR.txt",     sep="") )




MyScatterDiagram_2A(vector2X=reduceColumn3_Average_H3, vector2Y=log2(TPM_WT_7+1), 
                   path2=subdir_1_part7,   fileName2="Part7-1-A-H3-vs-TPM", 
                   xLab2="H3",   yLab2="TPM",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=10,   xMin2=0, xMax2=1.6,  alpha2=0.5)

MyScatterDiagram_2A(vector2X=reduceColumn3_Average_week0, vector2Y=log2(TPM_WT_7+1), 
                   path2=subdir_1_part7,   fileName2="Part7-1-B-week0H2BGFP-vs-TPM", 
                   xLab2="week0 H2BGFP",   yLab2="TPM",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=10,   xMin2=0, xMax2=1.2,  alpha2=0.5)


MyScatterDiagram_2A(vector2X=myNTR_7_WT_1, vector2Y=log2(TPM_WT_7+1), 
                   path2=subdir_1_part7,   fileName2="Part7-1-C-NTR-vs-TPM", 
                   xLab2="NTR",   yLab2="TPM",  title2=myTitle_g,  
                   height2=3.5,   width2=4, yMin2=0, yMax2=10,   xMin2=0, xMax2=0.4,  alpha2=0.5)




pdf( file=paste(subdir_1_part7, "Part7-1-D-H3-NTR-TPM.pdf", sep="/") )
scatterplot3d(x=reduceColumn3_Average_H3,  y=myNTR_7_WT_1,   z=log2(TPM_WT_7+1) , color="red", 
              main=NULL,   sub=NULL,    xlim=c(0, 1.6),   ylim=c(0, 0.4),   zlim=c(0, 10),
              xlab=NULL,   ylab=NULL,   zlab=NULL,   scale.y=1,   angle=40,
              axis=TRUE,   tick.marks=TRUE,   label.tick.marks=TRUE,
              x.ticklabs=NULL,   y.ticklabs=NULL,    z.ticklabs=NULL,
              y.margin.add=0,   grid=TRUE,   box=TRUE,
              type="p", pch=20 )
dev.off()









#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################



###############  average all rows for each columns
subdir_1_part4 <- paste(Part4_g,  "/1-NTR-averageAllRows", sep = "")
if( ! file.exists(subdir_1_part4) ) { dir.create(subdir_1_part4) }





####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################




