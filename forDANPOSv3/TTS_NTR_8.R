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
## Part  8:  Figures about a smaller matrix, so we can compute NTR at a single DNA region. Finally, a NTR matrix should be generated.  
############################################################################

subdir_1_part8 <- paste(Part8_g,  "/1-5columns", sep = "")
if( ! file.exists(subdir_1_part8) ) { dir.create(subdir_1_part8) }




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


sink( file=paste(subdir_1_part8,  "/8-1-1A-LogLinearModel.runLog",      sep = "") )
myNTR_WT_1A  <-  matrix(data = NA, nrow = nrow(reduceColumn1_Average_H3), ncol = 5 )
for (i in c(1:nrow(reduceColumn1_Average_H3))) {
    for (j in c(1:5)) { 
        vec1 <- c(reduceColumn1_Average_week0[i,j]+0.00001, reduceColumn1_Average_week1[i,j]+0.00001, reduceColumn1_Average_week2[i,j]+0.00001, 
                  reduceColumn1_Average_week4[i,j]+0.00001, reduceColumn1_Average_week6[i,j]+0.00001, reduceColumn1_Average_week8[i,j]+0.00001)
        NTR_bin <- MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2=log(vec1), file2=paste(subdir_1_part8,  "/8-1-1B-LogLinearModel", sep="")  )    
        myNTR_WT_1A[i,j] <- NTR_bin
    }
}
dim( myNTR_WT_1A )
sink()  


write.table(x=myNTR_WT_1A,    file = paste(Part8_g,  "Part8-1-A-NTR.txt", sep = "/"),    
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        



MyHeatmap_1(matrix2=myNTR_WT_1A,  
            path2=subdir_1_part8,  fileName2="Part8-1-A-heatmap",  
            title2=myTitle_g,  yLab2="2kb regions",  xLab2="genes or peaks",   
            height2=3, width2=5, midpoint2=0.15,  limits2=c(0, 0.3)  )   
  

####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################




