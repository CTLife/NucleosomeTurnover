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


subdir_1_part7 <- paste(Part7_g,  "/1-NTR-NOL-randIndex-5kb", sep = "")
if( ! file.exists(subdir_1_part7) ) { dir.create(subdir_1_part7) }


rankIndex <- seq(from = 1, to = length(myNTR_5_WT_1), by = 1)

length(myNTR_5_WT_1)
length(rankIndex)
length(reduceColumn3_Average_H3)      
length(reduceColumn3_Average_week0)



MyCor2Vars_4(vector1=row_Average_H3,     vector2=rankIndex,  file1=paste(subdir_1_part7,  "/Part7-1-1A-rankIndex-H3.txt",      sep="") )
MyCor2Vars_4(vector1=row_Average_week0,  vector2=rankIndex,  file1=paste(subdir_1_part7,  "/Part7-1-1B-rankIndex-week0.txt",   sep="") )
MyCor2Vars_4(vector1=myNTR_5_WT_1,       vector2=rankIndex,  file1=paste(subdir_1_part7,  "/Part7-1-1C-rankIndex-NTR.txt",     sep="") )



MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= row_Average_H3  , 
                   path2=subdir_1_part7,   fileName2="Part7-1-A-H3-vs-rankIndex", 
                   xLab2="rankIndex",   yLab2="H3",  title2=myTitle_g,  
                   height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.4 )
 
MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= row_Average_week0  , 
                    path2=subdir_1_part7,   fileName2="Part7-1-B-week0-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week0",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.1 )

MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= myNTR_5_WT_1  , 
                    path2=subdir_1_part7,   fileName2="Part7-1-C-NTR-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="NTR",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=0.3 )

MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= row_Average_week1  , 
                    path2=subdir_1_part7,   fileName2="Part7-1-D-week1-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week1",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.1 )

MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= row_Average_week2  , 
                    path2=subdir_1_part7,   fileName2="Part7-1-E-week2-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week2",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.1 )


MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= row_Average_week4  , 
                    path2=subdir_1_part7,   fileName2="Part7-1-F-week4-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week4",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.1 )


MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= row_Average_week6  , 
                    path2=subdir_1_part7,   fileName2="Part7-1-G-week6-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week6",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.1 )


MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= row_Average_week8  , 
                    path2=subdir_1_part7,   fileName2="Part7-1-H-week8-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week8",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.1 )








pdf( file=paste(subdir_1_part7, "Part7-1-D-H3-NTR-rankIndex.pdf", sep="/") )
scatterplot3d(x=row_Average_H3,  y=rankIndex ,  z=myNTR_5_WT_1,  color="red",   
              main=NULL,   sub=NULL,    xlim=c(0, 1.4),   ylim=c(1, length(rankIndex)),   zlim=c(0, 0.3) ,
              xlab=NULL,   ylab=NULL,   zlab=NULL,   scale.y=1,   angle=40,
              axis=TRUE,   tick.marks=TRUE,   label.tick.marks=TRUE,
              x.ticklabs=NULL,   y.ticklabs=NULL,    z.ticklabs=NULL,
              y.margin.add=0,   grid=TRUE,   box=TRUE,
              type="p", pch="." )

scatterplot3d(x=row_Average_week0,  y=rankIndex ,  z=myNTR_5_WT_1,  color="red",   
              main=NULL,   sub=NULL,    xlim=c(0, 1.4),   ylim=c(1, length(rankIndex)),   zlim=c(0, 0.3) ,
              xlab=NULL,   ylab=NULL,   zlab=NULL,   scale.y=1,   angle=40,
              axis=TRUE,   tick.marks=TRUE,   label.tick.marks=TRUE,
              x.ticklabs=NULL,   y.ticklabs=NULL,    z.ticklabs=NULL,
              y.margin.add=0,   grid=TRUE,   box=TRUE,
              type="p", pch="." )

scatterplot3d(x=row_Average_week8,  y=rankIndex ,  z=myNTR_5_WT_1,  color="red",   
              main=NULL,   sub=NULL,    xlim=c(0, 1.4),   ylim=c(1, length(rankIndex)),   zlim=c(0, 0.3) ,
              xlab=NULL,   ylab=NULL,   zlab=NULL,   scale.y=1,   angle=40,
              axis=TRUE,   tick.marks=TRUE,   label.tick.marks=TRUE,
              x.ticklabs=NULL,   y.ticklabs=NULL,    z.ticklabs=NULL,
              y.margin.add=0,   grid=TRUE,   box=TRUE,
              type="p", pch="." )
dev.off()




















########################################################################################
subdir_2_part7 <- paste(Part7_g,  "/2-NTR-NOL-randIndex-2kb", sep = "")
if( ! file.exists(subdir_2_part7) ) { dir.create(subdir_2_part7) }


rankIndex <- seq(from = 1, to = length(myNTR_6_WT_1), by = 1)

length(myNTR_6_WT_1)
length(rankIndex)
length(reduceColumn3_Average_H3)      
length(reduceColumn3_Average_week0)



MyCor2Vars_4(vector1=reduceColumn2_Average_H3,     vector2=rankIndex,  file1=paste(subdir_2_part7,  "/Part7-2-1A-rankIndex-H3.txt",      sep="") )
MyCor2Vars_4(vector1=reduceColumn2_Average_week0,  vector2=rankIndex,  file1=paste(subdir_2_part7,  "/Part7-2-1B-rankIndex-week0.txt",   sep="") )
MyCor2Vars_4(vector1=myNTR_6_WT_1,       vector2=rankIndex,  file1=paste(subdir_2_part7,  "/Part7-2-1C-rankIndex-NTR.txt",     sep="") )



MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= reduceColumn2_Average_H3  , 
                    path2=subdir_2_part7,   fileName2="Part7-2-A-H3-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="H3",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.4 )

MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= reduceColumn2_Average_week0  , 
                    path2=subdir_2_part7,   fileName2="Part7-2-B-week0-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week0",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.1 )

MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= myNTR_6_WT_1  , 
                    path2=subdir_2_part7,   fileName2="Part7-2-C-NTR-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="NTR",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=0.3 )

MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= reduceColumn2_Average_week1  , 
                    path2=subdir_2_part7,   fileName2="Part7-2-D-week1-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week1",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.1 )

MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= reduceColumn2_Average_week2  , 
                    path2=subdir_2_part7,   fileName2="Part7-2-E-week2-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week2",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.1 )


MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= reduceColumn2_Average_week4  , 
                    path2=subdir_2_part7,   fileName2="Part7-2-F-week4-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week4",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.1 )


MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= reduceColumn2_Average_week6  , 
                    path2=subdir_2_part7,   fileName2="Part7-2-G-week6-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week6",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.1 )


MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= reduceColumn2_Average_week8  , 
                    path2=subdir_2_part7,   fileName2="Part7-2-H-week8-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week8",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.1 )








pdf( file=paste(subdir_2_part7, "Part7-2-D-H3-NTR-rankIndex.pdf", sep="/") )
scatterplot3d(x=reduceColumn2_Average_H3,  y=rankIndex ,  z=myNTR_6_WT_1,  color="red",   
              main=NULL,   sub=NULL,    xlim=c(0, 1.4),   ylim=c(1, length(rankIndex)),   zlim=c(0, 0.3) ,
              xlab=NULL,   ylab=NULL,   zlab=NULL,   scale.y=1,   angle=40,
              axis=TRUE,   tick.marks=TRUE,   label.tick.marks=TRUE,
              x.ticklabs=NULL,   y.ticklabs=NULL,    z.ticklabs=NULL,
              y.margin.add=0,   grid=TRUE,   box=TRUE,
              type="p", pch="." )

scatterplot3d(x=reduceColumn2_Average_week0,  y=rankIndex ,  z=myNTR_6_WT_1,  color="red",   
              main=NULL,   sub=NULL,    xlim=c(0, 1.4),   ylim=c(1, length(rankIndex)),   zlim=c(0, 0.3) ,
              xlab=NULL,   ylab=NULL,   zlab=NULL,   scale.y=1,   angle=40,
              axis=TRUE,   tick.marks=TRUE,   label.tick.marks=TRUE,
              x.ticklabs=NULL,   y.ticklabs=NULL,    z.ticklabs=NULL,
              y.margin.add=0,   grid=TRUE,   box=TRUE,
              type="p", pch="." )

scatterplot3d(x=reduceColumn2_Average_week8,  y=rankIndex ,  z=myNTR_6_WT_1,  color="red",   
              main=NULL,   sub=NULL,    xlim=c(0, 1.4),   ylim=c(1, length(rankIndex)),   zlim=c(0, 0.3) ,
              xlab=NULL,   ylab=NULL,   zlab=NULL,   scale.y=1,   angle=40,
              axis=TRUE,   tick.marks=TRUE,   label.tick.marks=TRUE,
              x.ticklabs=NULL,   y.ticklabs=NULL,    z.ticklabs=NULL,
              y.margin.add=0,   grid=TRUE,   box=TRUE,
              type="p", pch="." )
dev.off()






















########################################################################################
subdir_3_part7 <- paste(Part7_g,  "/3-NTR-NOL-randIndex-1kb", sep = "")
if( ! file.exists(subdir_3_part7) ) { dir.create(subdir_3_part7) }


rankIndex <- seq(from = 1, to = length(myNTR_7_WT_1), by = 1)

length(myNTR_7_WT_1)
length(rankIndex)
length(reduceColumn3_Average_H3)      
length(reduceColumn3_Average_week0)



MyCor2Vars_4(vector1=reduceColumn3_Average_H3,     vector2=rankIndex,  file1=paste(subdir_3_part7,  "/Part7-3-1A-rankIndex-H3.txt",      sep="") )
MyCor2Vars_4(vector1=reduceColumn3_Average_week0,  vector2=rankIndex,  file1=paste(subdir_3_part7,  "/Part7-3-1B-rankIndex-week0.txt",   sep="") )
MyCor2Vars_4(vector1=myNTR_7_WT_1,       vector2=rankIndex,  file1=paste(subdir_3_part7,  "/Part7-3-1C-rankIndex-NTR.txt",     sep="") )



MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= reduceColumn3_Average_H3  , 
                    path2=subdir_3_part7,   fileName2="Part7-3-A-H3-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="H3",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.6 )

MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= reduceColumn3_Average_week0  , 
                    path2=subdir_3_part7,   fileName2="Part7-3-B-week0-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week0",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.2 )

MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= myNTR_7_WT_1  , 
                    path2=subdir_3_part7,   fileName2="Part7-3-C-NTR-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="NTR",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=0.4 )

MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= reduceColumn3_Average_week1  , 
                    path2=subdir_3_part7,   fileName2="Part7-3-D-week1-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week1",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.2 )

MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= reduceColumn3_Average_week2  , 
                    path2=subdir_3_part7,   fileName2="Part7-3-E-week2-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week2",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.2 )


MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= reduceColumn3_Average_week4  , 
                    path2=subdir_3_part7,   fileName2="Part7-3-F-week4-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week4",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.2 )


MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= reduceColumn3_Average_week6  , 
                    path2=subdir_3_part7,   fileName2="Part7-3-G-week6-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week6",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.2 )


MyScatterDiagram_5A(vector2X=rankIndex , vector2Y= reduceColumn3_Average_week8  , 
                    path2=subdir_3_part7,   fileName2="Part7-3-H-week8-vs-rankIndex", 
                    xLab2="rankIndex",   yLab2="week8",  title2=myTitle_g,  
                    height2=3.5,   width2=6, xMin2=1, xMax2=length(rankIndex),   yMin2=0, yMax2=1.2 )








pdf( file=paste(subdir_3_part7, "Part7-3-D-H3-NTR-rankIndex.pdf", sep="/") )
scatterplot3d(x=reduceColumn3_Average_H3,  y=rankIndex ,  z=myNTR_7_WT_1,  color="red",   
              main=NULL,   sub=NULL,    xlim=c(0, 1.6),   ylim=c(1, length(rankIndex)),   zlim=c(0, 0.4) ,
              xlab=NULL,   ylab=NULL,   zlab=NULL,   scale.y=1,   angle=40,
              axis=TRUE,   tick.marks=TRUE,   label.tick.marks=TRUE,
              x.ticklabs=NULL,   y.ticklabs=NULL,    z.ticklabs=NULL,
              y.margin.add=0,   grid=TRUE,   box=TRUE,
              type="p", pch="." )

scatterplot3d(x=reduceColumn3_Average_week0,  y=rankIndex ,  z=myNTR_7_WT_1,  color="red",   
              main=NULL,   sub=NULL,    xlim=c(0, 1.6),   ylim=c(1, length(rankIndex)),   zlim=c(0, 0.4) ,
              xlab=NULL,   ylab=NULL,   zlab=NULL,   scale.y=1,   angle=40,
              axis=TRUE,   tick.marks=TRUE,   label.tick.marks=TRUE,
              x.ticklabs=NULL,   y.ticklabs=NULL,    z.ticklabs=NULL,
              y.margin.add=0,   grid=TRUE,   box=TRUE,
              type="p", pch="." )

scatterplot3d(x=reduceColumn3_Average_week8,  y=rankIndex ,  z=myNTR_7_WT_1,  color="red",   
              main=NULL,   sub=NULL,    xlim=c(0, 1.6),   ylim=c(1, length(rankIndex)),   zlim=c(0, 0.4) ,
              xlab=NULL,   ylab=NULL,   zlab=NULL,   scale.y=1,   angle=40,
              axis=TRUE,   tick.marks=TRUE,   label.tick.marks=TRUE,
              x.ticklabs=NULL,   y.ticklabs=NULL,    z.ticklabs=NULL,
              y.margin.add=0,   grid=TRUE,   box=TRUE,
              type="p", pch="." )
dev.off()








####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################




