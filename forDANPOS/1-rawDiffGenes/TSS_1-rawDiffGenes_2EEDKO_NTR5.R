#############################################################################################################################
## Part 5:  Classify DNA regions based on NOL.
#############################################################################################################################





#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################


subdir_1_part5 <- paste(Part5_g,  "/1-one1-diff-WEEK0signals", sep = "")
if( ! file.exists(subdir_1_part5) ) { dir.create(subdir_1_part5) }


##两样本的wilcoxon检验，更多的时候，它被称之为Mann-Whitney U检验，不过在R中，都是wilcox.test函数来计算，同样基于把数据用rank替换(不区分样本)，然后计算其中一个样本的秩和，
##这样子就把问题简化成从1到n1+n2中无放回取出n1个值的样本。零假设是两个样本来自于同一个分布，那么在这种情况下两个样本的秩是随机的，p值计算的是实际数据中秩分布的可能性有多大。
##如果零假设是真的，那么两份数据拥有相同的位置; 而备择假设则认为两份数据的位置有显著的位移。
##和t检验其实是不同的，t检验所检验的是单个的数值（均值），而Mann-Whitney U检验的是数据的分布，而不是单一的值。  

##闵可夫斯基距离(Minkowski Distance): d(X, Y) = [sum|Xi-Yi|^p]^(1/p)
##其中p是一个变参数。
##当p=1时，就是曼哈顿距离 (Manhattan Distance)
##当p=2时，就是欧氏距离 (Euclidean Distance)
##当p→∞时，就是切比雪夫距离 ( Chebyshev Distance ) 
##根据变参数的不同，闵氏距离可以表示一类的距离。
## dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)  
##  method: the distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
##  P: The power of the Minkowski distance.


dim(Average_one1_week0_EEDheto)
dim(Average_one1_week0_EEDko)

one1_week0_numcols <- ncol(Average_one1_week0_EEDko)
one1_week0_numrows <- nrow(Average_one1_week0_EEDko)
one1_week0_pvalues <- vector(length = one1_week0_numrows)
one1_week0_dist    <- vector(length = one1_week0_numrows)

for ( i in 1:one1_week0_numrows  )  {
        wilcoxTemp <- wilcox.test(x=Average_one1_week0_EEDheto[i,], y=Average_one1_week0_EEDko[i, ], alternative="two.sided",  mu=0,   paired=TRUE,   exact=NULL, correct=TRUE,  conf.int=FALSE,  conf.level=0.95)
        one1_week0_pvalues[i] <- wilcoxTemp$p.value
        if( is.na(one1_week0_pvalues[i]) ) {one1_week0_pvalues[i] = 0}
        ##one1_week0_dist[i]    <- dist( rbind(vector1, vector2), method="minkowski", diag=FALSE, upper=FALSE, p=1)
        one1_week0_dist[i]    <- mean( x=( Average_one1_week0_EEDheto[i,201:400] - Average_one1_week0_EEDko[i, 201:400] ),  na.rm = TRUE )
}

summary(one1_week0_dist)    
one1_week0_dist[1]
summary(one1_week0_dist)
summary(one1_week0_pvalues)


length(one1_week0_dist)
length(one1_week0_pvalues)
length( one1_week0_pvalues[one1_week0_pvalues>=0.05] )
length( one1_week0_pvalues[one1_week0_pvalues < 0.05] )
length( one1_week0_dist[one1_week0_dist < -0.08] )  ## one1_week0_EEDheto is less than one1_week0_EEDko
length( one1_week0_dist[one1_week0_dist >  0.08] )  ## one1_week0_EEDheto is more than one1_week0_EEDko


pValue_threshold=0.05
distance_threshold=0.08


## 差异不显著
one1_NOL_diffIndex_1 <- ( (one1_week0_pvalues>pValue_threshold)  &  (abs(one1_week0_dist)<distance_threshold) )
length(one1_NOL_diffIndex_1[one1_NOL_diffIndex_1])
write.table(x=one1_NOL_diffIndex_1,   file = paste(subdir_1_part5, "/1-one1-NoDiff.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")


## 差异显著  one1_week0_EEDheto is less than one1_week0_EEDko
one1_NOL_diffIndex_2 <- ( (one1_week0_pvalues<pValue_threshold)  &  (one1_week0_dist < -distance_threshold) )
length(one1_NOL_diffIndex_2[one1_NOL_diffIndex_2])
write.table(x=one1_NOL_diffIndex_2,   file = paste(subdir_1_part5, "/1-one1-EEDheto-Less-Than-EEDko.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")


## 差异显著  one1_week0_EEDheto is more than one1_week0_EEDko
one1_NOL_diffIndex_3 <- ( (one1_week0_pvalues<pValue_threshold)  &  (one1_week0_dist > distance_threshold) )
length(one1_NOL_diffIndex_3[one1_NOL_diffIndex_3])
write.table(x=one1_NOL_diffIndex_3,   file = paste(subdir_1_part5, "/1-one1-EEDheto-More-Than-EEDko.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")






MyAverageLines_3(vector2=c(colMeans(Average_one1_week0_EEDheto[one1_NOL_diffIndex_1, ]),  colMeans(Average_one1_week0_EEDko[one1_NOL_diffIndex_1, ]), 
                           colMeans(Average_one1_week4_EEDheto[one1_NOL_diffIndex_1, ]),  colMeans(Average_one1_week4_EEDko[one1_NOL_diffIndex_1, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part5,     fileName2="A-one1-NoDiff",  
                 title2="Down-regulated Genes (No Diff)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min,   Ymax2=col_NOL_max,    height2=3.2,   width2=5.55 , center2="TSS" )


one1_len_A <- nrow(Average_one1_week0_EEDheto[one1_NOL_diffIndex_1, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_one1_week0_EEDheto[one1_NOL_diffIndex_1, ]),  rowMeans(Average_one1_week0_EEDko[one1_NOL_diffIndex_1, ]), 
                            rowMeans(Average_one1_week4_EEDheto[one1_NOL_diffIndex_1, ]),  rowMeans(Average_one1_week4_EEDko[one1_NOL_diffIndex_1, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", one1_len_A),   rep("week0_EEDko",  one1_len_A), 
                                 rep("week4_EEDheto", one1_len_A),   rep("week4_EEDko",  one1_len_A) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_1_part5,  fileName2="A-one1-all-NoDiff-BoxViolin",  
                  title2="Down-regulated Genes  (No Diff)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_one1_week0_EEDheto[one1_NOL_diffIndex_1, 201:400]),  rowMeans(Average_one1_week0_EEDko[one1_NOL_diffIndex_1, 201:400]), 
                            rowMeans(Average_one1_week4_EEDheto[one1_NOL_diffIndex_1, 201:400]),  rowMeans(Average_one1_week4_EEDko[one1_NOL_diffIndex_1, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", one1_len_A),   rep("week0_EEDko",  one1_len_A), 
                                 rep("week4_EEDheto", one1_len_A),   rep("week4_EEDko",  one1_len_A) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_1_part5,  fileName2="A-one1-2kbRegion-NoDiff-BoxViolin",  
                  title2="Down-regulated Genes  (No Diff)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )










MyAverageLines_3(vector2=c(colMeans(Average_one1_week0_EEDheto[one1_NOL_diffIndex_2, ]),  colMeans(Average_one1_week0_EEDko[one1_NOL_diffIndex_2, ]), 
                           colMeans(Average_one1_week4_EEDheto[one1_NOL_diffIndex_2, ]),  colMeans(Average_one1_week4_EEDko[one1_NOL_diffIndex_2, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part5,     fileName2="B-one1-EEDhetoLess",  
                 title2="Down-regulated Genes (EEDheto Less)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min,   Ymax2=col_NOL_max,    height2=3.2,   width2=5.55 , center2="TSS" )


one1_len_B <- nrow(Average_one1_week0_EEDheto[one1_NOL_diffIndex_2, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_one1_week0_EEDheto[one1_NOL_diffIndex_2, ]),  rowMeans(Average_one1_week0_EEDko[one1_NOL_diffIndex_2, ]), 
                            rowMeans(Average_one1_week4_EEDheto[one1_NOL_diffIndex_2, ]),  rowMeans(Average_one1_week4_EEDko[one1_NOL_diffIndex_2, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", one1_len_B),   rep("week0_EEDko",  one1_len_B), 
                                 rep("week4_EEDheto", one1_len_B),   rep("week4_EEDko",  one1_len_B) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_1_part5,  fileName2="B-one1-all-EEDhetoLess-BoxViolin",  
                  title2="Down-regulated Genes  (EEDheto Less)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_one1_week0_EEDheto[one1_NOL_diffIndex_2, 201:400]),  rowMeans(Average_one1_week0_EEDko[one1_NOL_diffIndex_2, 201:400]), 
                            rowMeans(Average_one1_week4_EEDheto[one1_NOL_diffIndex_2, 201:400]),  rowMeans(Average_one1_week4_EEDko[one1_NOL_diffIndex_2, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", one1_len_B),   rep("week0_EEDko",  one1_len_B), 
                                 rep("week4_EEDheto", one1_len_B),   rep("week4_EEDko",  one1_len_B) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_1_part5,  fileName2="B-one1-2kbRegion-EEDhetoLess-BoxViolin",  
                  title2="Down-regulated Genes  (EEDheto Less)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )












MyAverageLines_3(vector2=c(colMeans(Average_one1_week0_EEDheto[one1_NOL_diffIndex_3, ]),  colMeans(Average_one1_week0_EEDko[one1_NOL_diffIndex_3, ]), 
                           colMeans(Average_one1_week4_EEDheto[one1_NOL_diffIndex_3, ]),  colMeans(Average_one1_week4_EEDko[one1_NOL_diffIndex_3, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part5,     fileName2="C-one1-EEDhetoMore",  
                 title2="Down-regulated Genes (EEDheto More)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min,   Ymax2=col_NOL_max,    height2=3.2,   width2=5.55 , center2="TSS" )


one1_len_C <- nrow(Average_one1_week0_EEDheto[one1_NOL_diffIndex_3, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_one1_week0_EEDheto[one1_NOL_diffIndex_3, ]),  rowMeans(Average_one1_week0_EEDko[one1_NOL_diffIndex_3, ]), 
                            rowMeans(Average_one1_week4_EEDheto[one1_NOL_diffIndex_3, ]),  rowMeans(Average_one1_week4_EEDko[one1_NOL_diffIndex_3, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", one1_len_C),   rep("week0_EEDko",  one1_len_C), 
                                 rep("week4_EEDheto", one1_len_C),   rep("week4_EEDko",  one1_len_C) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_1_part5,  fileName2="C-one1-all-EEDhetoMore-BoxViolin",  
                  title2="Down-regulated Genes  (EEDheto More)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_one1_week0_EEDheto[one1_NOL_diffIndex_3, 201:400]),  rowMeans(Average_one1_week0_EEDko[one1_NOL_diffIndex_3, 201:400]), 
                            rowMeans(Average_one1_week4_EEDheto[one1_NOL_diffIndex_3, 201:400]),  rowMeans(Average_one1_week4_EEDko[one1_NOL_diffIndex_3, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", one1_len_C),   rep("week0_EEDko",  one1_len_C), 
                                 rep("week4_EEDheto", one1_len_C),   rep("week4_EEDko",  one1_len_C) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_1_part5,  fileName2="C-one1-2kbRegion-EEDhetoMore-BoxViolin",  
                  title2="Down-regulated Genes  (EEDheto More)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )











column_NTR_one1_EEDheto_same   <-  log( colMeans(Average_one1_week0_EEDheto[one1_NOL_diffIndex_1, ]) / colMeans(Average_one1_week4_EEDheto[one1_NOL_diffIndex_1, ]) )/4
column_NTR_one1_EEDko_same     <-  log( colMeans(Average_one1_week0_EEDko[one1_NOL_diffIndex_1, ])   / colMeans(Average_one1_week4_EEDko[one1_NOL_diffIndex_1, ])   )/4
column_NTR_one1_EEDheto_less   <-  log( colMeans(Average_one1_week0_EEDheto[one1_NOL_diffIndex_2, ]) / colMeans(Average_one1_week4_EEDheto[one1_NOL_diffIndex_2, ]) )/4
column_NTR_one1_EEDko_less     <-  log( colMeans(Average_one1_week0_EEDko[one1_NOL_diffIndex_2, ])   / colMeans(Average_one1_week4_EEDko[one1_NOL_diffIndex_2, ])   )/4
column_NTR_one1_EEDheto_more   <-  log( colMeans(Average_one1_week0_EEDheto[one1_NOL_diffIndex_3, ]) / colMeans(Average_one1_week4_EEDheto[one1_NOL_diffIndex_3, ]) )/4
column_NTR_one1_EEDko_more     <-  log( colMeans(Average_one1_week0_EEDko[one1_NOL_diffIndex_3, ])   / colMeans(Average_one1_week4_EEDko[one1_NOL_diffIndex_3, ])   )/4

row_NTR_one1_EEDheto_same   <-  log( row_Average_one1_week0_EEDheto[one1_NOL_diffIndex_1] / row_Average_one1_week4_EEDheto[one1_NOL_diffIndex_1] )/4
row_NTR_one1_EEDko_same     <-  log( row_Average_one1_week0_EEDko[one1_NOL_diffIndex_1]   / row_Average_one1_week4_EEDko[one1_NOL_diffIndex_1]   )/4
row_NTR_one1_EEDheto_less   <-  log( row_Average_one1_week0_EEDheto[one1_NOL_diffIndex_2] / row_Average_one1_week4_EEDheto[one1_NOL_diffIndex_2] )/4
row_NTR_one1_EEDko_less     <-  log( row_Average_one1_week0_EEDko[one1_NOL_diffIndex_2]   / row_Average_one1_week4_EEDko[one1_NOL_diffIndex_2]   )/4
row_NTR_one1_EEDheto_more   <-  log( row_Average_one1_week0_EEDheto[one1_NOL_diffIndex_3] / row_Average_one1_week4_EEDheto[one1_NOL_diffIndex_3] )/4
row_NTR_one1_EEDko_more     <-  log( row_Average_one1_week0_EEDko[one1_NOL_diffIndex_3]   / row_Average_one1_week4_EEDko[one1_NOL_diffIndex_3]   )/4


MyAverageLines_3(vector2=c(column_NTR_one1_EEDheto_same,  column_NTR_one1_EEDko_same,  column_NTR_one1_EEDheto_less, column_NTR_one1_EEDko_less,
                           column_NTR_one1_EEDheto_more, column_NTR_one1_EEDko_more ),   
                 numSample2=6,   
                 sampleType2=c( rep("EEDheto_same", numOfColumns1),   rep("EEDko_same", numOfColumns1),  rep("EEDheto_less", numOfColumns1),   rep("EEDko_less", numOfColumns1), 
                                rep("EEDheto_more", numOfColumns1),   rep("EEDko_more", numOfColumns1)  ), 
                 sampleRank2=c( "EEDheto_same",   "EEDko_same",  "EEDheto_less",  "EEDko_less",  "EEDheto_more",  "EEDko_more"  ),  
                 colours2=c("EEDheto_same"="red",   "EEDko_same"="red4",  "EEDheto_less"="blue",  "EEDko_less"="blue4",  "EEDheto_more"="green",  "EEDko_more"="green4"),   
                 path2=subdir_1_part5,     fileName2="1-one1-NTR",  
                 title2="Down-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR_one1_EEDheto_same,  column_NTR_one1_EEDko_same ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_same", numOfColumns1),   rep("EEDko_same", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_same",   "EEDko_same"  ),  
                 colours2=c("EEDheto_same"="red",   "EEDko_same"="red4" ),   
                 path2=subdir_1_part5,     fileName2="1-one1-NTR-same",  
                 title2="Down-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR_one1_EEDheto_less,  column_NTR_one1_EEDko_less ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_less", numOfColumns1),   rep("EEDko_less", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_less",   "EEDko_less"  ),  
                 colours2=c("EEDheto_less"="blue",   "EEDko_less"="blue4" ),   
                 path2=subdir_1_part5,     fileName2="1-one1-NTR-less",  
                 title2="Down-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR_one1_EEDheto_more,  column_NTR_one1_EEDko_more ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_more", numOfColumns1),   rep("EEDko_more", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_more",   "EEDko_more"  ),  
                 colours2=c("EEDheto_more"="blue",   "EEDko_more"="blue4" ),   
                 path2=subdir_1_part5,     fileName2="1-one1-NTR-more",  
                 title2="Down-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyBoxViolinPlot_1(vector2=c(row_NTR_one1_EEDheto_same,  row_NTR_one1_EEDko_same,  row_NTR_one1_EEDheto_less, row_NTR_one1_EEDko_less,
                            row_NTR_one1_EEDheto_more,  row_NTR_one1_EEDko_more  ),  
                  sampleType2=c( rep("EEDheto_same", one1_len_A),   rep("EEDko_same", one1_len_A),  rep("EEDheto_less", one1_len_B),   rep("EEDko_less", one1_len_B), 
                                 rep("EEDheto_more", one1_len_C),   rep("EEDko_more", one1_len_C)  ), 
                  sampleRank2=c( "EEDheto_same",   "EEDko_same",  "EEDheto_less",  "EEDko_less",  "EEDheto_more",  "EEDko_more"  ),  
                  colours2=c("red",  "red4",   "blue",    "blue4", "green", "green4" ),   
                  path2=subdir_1_part5,  fileName2="1-one1-NTR-boxViolin",  
                  title2="Down-regulated Genes  (EEDheto More)",  xLab2="Samples",  yLab2="NTR",   
                  height2=3.88,   width2=4, Ymin2=-0.3, Ymax2=1 )





















################################################################################
subdir_2_part5 <- paste(Part5_g,  "/2-two2-diff-WEEK0signals", sep = "")
if( ! file.exists(subdir_2_part5) ) { dir.create(subdir_2_part5) }


dim(Average_two2_week0_EEDheto)
dim(Average_two2_week0_EEDko)

two2_week0_numcols <- ncol(Average_two2_week0_EEDko)
two2_week0_numrows <- nrow(Average_two2_week0_EEDko)
two2_week0_pvalues <- vector(length = two2_week0_numrows)
two2_week0_dist    <- vector(length = two2_week0_numrows)

for ( i in 1:two2_week0_numrows  )  {
        wilcoxTemp <- wilcox.test(x=Average_two2_week0_EEDheto[i,], y=Average_two2_week0_EEDko[i, ], alternative="two.sided",  mu=0,   paired=TRUE,   exact=NULL, correct=TRUE,  conf.int=FALSE,  conf.level=0.95)
        two2_week0_pvalues[i] <- wilcoxTemp$p.value
        if( is.na(two2_week0_pvalues[i]) ) {two2_week0_pvalues[i] = 0}
        ##two2_week0_dist[i]    <- dist( rbind(vector1, vector2), method="minkowski", diag=FALSE, upper=FALSE, p=1)
        two2_week0_dist[i]    <- mean( x=( Average_two2_week0_EEDheto[i,201:400] - Average_two2_week0_EEDko[i, 201:400] ),  na.rm = TRUE )
}

summary(two2_week0_dist)    
two2_week0_dist[1]
summary(two2_week0_pvalues)
summary(two2_week0_dist)

length(two2_week0_dist)
length(two2_week0_pvalues)
length( two2_week0_pvalues[two2_week0_pvalues>=0.05] )
length( two2_week0_pvalues[two2_week0_pvalues < 0.05] )
length( two2_week0_dist[two2_week0_dist < -0.08] )  ## two2_week0_EEDheto is less than two2_week0_EEDko
length( two2_week0_dist[two2_week0_dist >  0.08] )  ## two2_week0_EEDheto is more than two2_week0_EEDko


pValue_threshold=0.05
distance_threshold=0.08


## 差异不显著
two2_NOL_diffIndex_1 <- ( (two2_week0_pvalues>pValue_threshold)  &  (abs(two2_week0_dist)<distance_threshold) )
length(two2_NOL_diffIndex_1[two2_NOL_diffIndex_1])
write.table(x=two2_NOL_diffIndex_1,   file = paste(subdir_2_part5, "/1-two2-NoDiff.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")


## 差异显著  two2_week0_EEDheto is less than two2_week0_EEDko
two2_NOL_diffIndex_2 <- ( (two2_week0_pvalues<pValue_threshold)  &  (two2_week0_dist < -distance_threshold) )
length(two2_NOL_diffIndex_2[two2_NOL_diffIndex_2])
write.table(x=two2_NOL_diffIndex_2,   file = paste(subdir_2_part5, "/1-two2-EEDheto-Less-Than-EEDko.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")


## 差异显著  two2_week0_EEDheto is more than two2_week0_EEDko
two2_NOL_diffIndex_3 <- ( (two2_week0_pvalues<pValue_threshold)  &  (two2_week0_dist > distance_threshold) )
length(two2_NOL_diffIndex_3[two2_NOL_diffIndex_3])
write.table(x=two2_NOL_diffIndex_3,   file = paste(subdir_2_part5, "/1-two2-EEDheto-More-Than-EEDko.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")






MyAverageLines_3(vector2=c(colMeans(Average_two2_week0_EEDheto[two2_NOL_diffIndex_1, ]),  colMeans(Average_two2_week0_EEDko[two2_NOL_diffIndex_1, ]), 
                           colMeans(Average_two2_week4_EEDheto[two2_NOL_diffIndex_1, ]),  colMeans(Average_two2_week4_EEDko[two2_NOL_diffIndex_1, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_2_part5,     fileName2="A-two2-NoDiff",  
                 title2="Unchanged Genes (No Diff)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min,   Ymax2=col_NOL_max,    height2=3.2,   width2=5.55 , center2="TSS" )


two2_len_A <- nrow(Average_two2_week0_EEDheto[two2_NOL_diffIndex_1, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_two2_week0_EEDheto[two2_NOL_diffIndex_1, ]),  rowMeans(Average_two2_week0_EEDko[two2_NOL_diffIndex_1, ]), 
                            rowMeans(Average_two2_week4_EEDheto[two2_NOL_diffIndex_1, ]),  rowMeans(Average_two2_week4_EEDko[two2_NOL_diffIndex_1, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", two2_len_A),   rep("week0_EEDko",  two2_len_A), 
                                 rep("week4_EEDheto", two2_len_A),   rep("week4_EEDko",  two2_len_A) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_2_part5,  fileName2="A-two2-all-NoDiff-BoxViolin",  
                  title2="Unchanged Genes  (No Diff)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_two2_week0_EEDheto[two2_NOL_diffIndex_1, 201:400]),  rowMeans(Average_two2_week0_EEDko[two2_NOL_diffIndex_1, 201:400]), 
                            rowMeans(Average_two2_week4_EEDheto[two2_NOL_diffIndex_1, 201:400]),  rowMeans(Average_two2_week4_EEDko[two2_NOL_diffIndex_1, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", two2_len_A),   rep("week0_EEDko",  two2_len_A), 
                                 rep("week4_EEDheto", two2_len_A),   rep("week4_EEDko",  two2_len_A) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_2_part5,  fileName2="A-two2-2kbRegion-NoDiff-BoxViolin",  
                  title2="Unchanged Genes  (No Diff)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )










MyAverageLines_3(vector2=c(colMeans(Average_two2_week0_EEDheto[two2_NOL_diffIndex_2, ]),  colMeans(Average_two2_week0_EEDko[two2_NOL_diffIndex_2, ]), 
                           colMeans(Average_two2_week4_EEDheto[two2_NOL_diffIndex_2, ]),  colMeans(Average_two2_week4_EEDko[two2_NOL_diffIndex_2, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_2_part5,     fileName2="B-two2-EEDhetoLess",  
                 title2="Unchanged Genes (EEDheto Less)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min,   Ymax2=col_NOL_max,    height2=3.2,   width2=5.55 , center2="TSS" )


two2_len_B <- nrow(Average_two2_week0_EEDheto[two2_NOL_diffIndex_2, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_two2_week0_EEDheto[two2_NOL_diffIndex_2, ]),  rowMeans(Average_two2_week0_EEDko[two2_NOL_diffIndex_2, ]), 
                            rowMeans(Average_two2_week4_EEDheto[two2_NOL_diffIndex_2, ]),  rowMeans(Average_two2_week4_EEDko[two2_NOL_diffIndex_2, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", two2_len_B),   rep("week0_EEDko",  two2_len_B), 
                                 rep("week4_EEDheto", two2_len_B),   rep("week4_EEDko",  two2_len_B) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_2_part5,  fileName2="B-two2-all-EEDhetoLess-BoxViolin",  
                  title2="Unchanged Genes  (EEDheto Less)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_two2_week0_EEDheto[two2_NOL_diffIndex_2, 201:400]),  rowMeans(Average_two2_week0_EEDko[two2_NOL_diffIndex_2, 201:400]), 
                            rowMeans(Average_two2_week4_EEDheto[two2_NOL_diffIndex_2, 201:400]),  rowMeans(Average_two2_week4_EEDko[two2_NOL_diffIndex_2, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", two2_len_B),   rep("week0_EEDko",  two2_len_B), 
                                 rep("week4_EEDheto", two2_len_B),   rep("week4_EEDko",  two2_len_B) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_2_part5,  fileName2="B-two2-2kbRegion-EEDhetoLess-BoxViolin",  
                  title2="Unchanged Genes  (EEDheto Less)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )












MyAverageLines_3(vector2=c(colMeans(Average_two2_week0_EEDheto[two2_NOL_diffIndex_3, ]),  colMeans(Average_two2_week0_EEDko[two2_NOL_diffIndex_3, ]), 
                           colMeans(Average_two2_week4_EEDheto[two2_NOL_diffIndex_3, ]),  colMeans(Average_two2_week4_EEDko[two2_NOL_diffIndex_3, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_2_part5,     fileName2="C-two2-EEDhetoMore",  
                 title2="Unchanged Genes Genes (EEDheto More)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min,   Ymax2=col_NOL_max,    height2=3.2,   width2=5.55 , center2="TSS" )


two2_len_C <- nrow(Average_two2_week0_EEDheto[two2_NOL_diffIndex_3, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_two2_week0_EEDheto[two2_NOL_diffIndex_3, ]),  rowMeans(Average_two2_week0_EEDko[two2_NOL_diffIndex_3, ]), 
                            rowMeans(Average_two2_week4_EEDheto[two2_NOL_diffIndex_3, ]),  rowMeans(Average_two2_week4_EEDko[two2_NOL_diffIndex_3, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", two2_len_C),   rep("week0_EEDko",  two2_len_C), 
                                 rep("week4_EEDheto", two2_len_C),   rep("week4_EEDko",  two2_len_C) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_2_part5,  fileName2="C-two2-all-EEDhetoMore-BoxViolin",  
                  title2="Unchanged Genes Genes  (EEDheto More)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_two2_week0_EEDheto[two2_NOL_diffIndex_3, 201:400]),  rowMeans(Average_two2_week0_EEDko[two2_NOL_diffIndex_3, 201:400]), 
                            rowMeans(Average_two2_week4_EEDheto[two2_NOL_diffIndex_3, 201:400]),  rowMeans(Average_two2_week4_EEDko[two2_NOL_diffIndex_3, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", two2_len_C),   rep("week0_EEDko",  two2_len_C), 
                                 rep("week4_EEDheto", two2_len_C),   rep("week4_EEDko",  two2_len_C) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_2_part5,  fileName2="C-two2-2kbRegion-EEDhetoMore-BoxViolin",  
                  title2="Unchanged Genes Genes  (EEDheto More)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )












column_NTR_two2_EEDheto_same   <-  log( colMeans(Average_two2_week0_EEDheto[two2_NOL_diffIndex_1, ]) / colMeans(Average_two2_week4_EEDheto[two2_NOL_diffIndex_1, ]) )/4
column_NTR_two2_EEDko_same     <-  log( colMeans(Average_two2_week0_EEDko[two2_NOL_diffIndex_1, ])   / colMeans(Average_two2_week4_EEDko[two2_NOL_diffIndex_1, ])   )/4
column_NTR_two2_EEDheto_less   <-  log( colMeans(Average_two2_week0_EEDheto[two2_NOL_diffIndex_2, ]) / colMeans(Average_two2_week4_EEDheto[two2_NOL_diffIndex_2, ]) )/4
column_NTR_two2_EEDko_less     <-  log( colMeans(Average_two2_week0_EEDko[two2_NOL_diffIndex_2, ])   / colMeans(Average_two2_week4_EEDko[two2_NOL_diffIndex_2, ])   )/4
column_NTR_two2_EEDheto_more   <-  log( colMeans(Average_two2_week0_EEDheto[two2_NOL_diffIndex_3, ]) / colMeans(Average_two2_week4_EEDheto[two2_NOL_diffIndex_3, ]) )/4
column_NTR_two2_EEDko_more     <-  log( colMeans(Average_two2_week0_EEDko[two2_NOL_diffIndex_3, ])   / colMeans(Average_two2_week4_EEDko[two2_NOL_diffIndex_3, ])   )/4

row_NTR_two2_EEDheto_same   <-  log( row_Average_two2_week0_EEDheto[two2_NOL_diffIndex_1] / row_Average_two2_week4_EEDheto[two2_NOL_diffIndex_1] )/4
row_NTR_two2_EEDko_same     <-  log( row_Average_two2_week0_EEDko[two2_NOL_diffIndex_1]   / row_Average_two2_week4_EEDko[two2_NOL_diffIndex_1]   )/4
row_NTR_two2_EEDheto_less   <-  log( row_Average_two2_week0_EEDheto[two2_NOL_diffIndex_2] / row_Average_two2_week4_EEDheto[two2_NOL_diffIndex_2] )/4
row_NTR_two2_EEDko_less     <-  log( row_Average_two2_week0_EEDko[two2_NOL_diffIndex_2]   / row_Average_two2_week4_EEDko[two2_NOL_diffIndex_2]   )/4
row_NTR_two2_EEDheto_more   <-  log( row_Average_two2_week0_EEDheto[two2_NOL_diffIndex_3] / row_Average_two2_week4_EEDheto[two2_NOL_diffIndex_3] )/4
row_NTR_two2_EEDko_more     <-  log( row_Average_two2_week0_EEDko[two2_NOL_diffIndex_3]   / row_Average_two2_week4_EEDko[two2_NOL_diffIndex_3]   )/4


MyAverageLines_3(vector2=c(column_NTR_two2_EEDheto_same,  column_NTR_two2_EEDko_same,  column_NTR_two2_EEDheto_less, column_NTR_two2_EEDko_less,
                           column_NTR_two2_EEDheto_more, column_NTR_two2_EEDko_more ),   
                 numSample2=6,   
                 sampleType2=c( rep("EEDheto_same", numOfColumns1),   rep("EEDko_same", numOfColumns1),  rep("EEDheto_less", numOfColumns1),   rep("EEDko_less", numOfColumns1), 
                                rep("EEDheto_more", numOfColumns1),   rep("EEDko_more", numOfColumns1)  ), 
                 sampleRank2=c( "EEDheto_same",   "EEDko_same",  "EEDheto_less",  "EEDko_less",  "EEDheto_more",  "EEDko_more"  ),  
                 colours2=c("EEDheto_same"="red",   "EEDko_same"="red4",  "EEDheto_less"="blue",  "EEDko_less"="blue4",  "EEDheto_more"="green",  "EEDko_more"="green4"),   
                 path2=subdir_2_part5,     fileName2="1-two2-NTR",  
                 title2="Unchanged Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR_two2_EEDheto_same,  column_NTR_two2_EEDko_same ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_same", numOfColumns1),   rep("EEDko_same", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_same",   "EEDko_same"  ),  
                 colours2=c("EEDheto_same"="red",   "EEDko_same"="red4" ),   
                 path2=subdir_2_part5,     fileName2="1-two2-NTR-same",  
                 title2="Unchanged Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR_two2_EEDheto_less,  column_NTR_two2_EEDko_less ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_less", numOfColumns1),   rep("EEDko_less", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_less",   "EEDko_less"  ),  
                 colours2=c("EEDheto_less"="blue",   "EEDko_less"="blue4" ),   
                 path2=subdir_2_part5,     fileName2="1-two2-NTR-less",  
                 title2="Unchanged Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR_two2_EEDheto_more,  column_NTR_two2_EEDko_more ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_more", numOfColumns1),   rep("EEDko_more", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_more",   "EEDko_more"  ),  
                 colours2=c("EEDheto_more"="blue",   "EEDko_more"="blue4" ),   
                 path2=subdir_2_part5,     fileName2="1-two2-NTR-more",  
                 title2="Unchanged Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyBoxViolinPlot_1(vector2=c(row_NTR_two2_EEDheto_same,  row_NTR_two2_EEDko_same,  row_NTR_two2_EEDheto_less, row_NTR_two2_EEDko_less,
                            row_NTR_two2_EEDheto_more,  row_NTR_two2_EEDko_more  ),  
                  sampleType2=c( rep("EEDheto_same", two2_len_A),   rep("EEDko_same", two2_len_A),  rep("EEDheto_less", two2_len_B),   rep("EEDko_less", two2_len_B), 
                                 rep("EEDheto_more", two2_len_C),   rep("EEDko_more", two2_len_C)  ), 
                  sampleRank2=c( "EEDheto_same",   "EEDko_same",  "EEDheto_less",  "EEDko_less",  "EEDheto_more",  "EEDko_more"  ),  
                  colours2=c("red",  "red4",   "blue",    "blue4", "green", "green4" ),   
                  path2=subdir_2_part5,  fileName2="1-two2-NTR-boxViolin",  
                  title2="Unchanged Genes",  xLab2="Samples",  yLab2="NTR",   
                  height2=3.88,   width2=4, Ymin2=-0.3, Ymax2=1 )






























################################################################################
subdir_3_part5 <- paste(Part5_g,  "/3-three3-diff-WEEK0signals", sep = "")
if( ! file.exists(subdir_3_part5) ) { dir.create(subdir_3_part5) }


dim(Average_three3_week0_EEDheto)
dim(Average_three3_week0_EEDko)

three3_week0_numcols <- ncol(Average_three3_week0_EEDko)
three3_week0_numrows <- nrow(Average_three3_week0_EEDko)
three3_week0_pvalues <- vector(length = three3_week0_numrows)
three3_week0_dist    <- vector(length = three3_week0_numrows)

for ( i in 1:three3_week0_numrows  )  {
        wilcoxTemp <- wilcox.test(x=Average_three3_week0_EEDheto[i,], y=Average_three3_week0_EEDko[i, ], alternative="two.sided",  mu=0,   paired=TRUE,   exact=NULL, correct=TRUE,  conf.int=FALSE,  conf.level=0.95)
        three3_week0_pvalues[i] <- wilcoxTemp$p.value
        if( is.na(three3_week0_pvalues[i]) ) {three3_week0_pvalues[i] = 0}
        ##three3_week0_dist[i]    <- dist( rbind(vector1, vector2), method="minkowski", diag=FALSE, upper=FALSE, p=1)
        three3_week0_dist[i]    <- mean( x=( Average_three3_week0_EEDheto[i,201:400] - Average_three3_week0_EEDko[i, 201:400] ),  na.rm = TRUE )
}

summary(three3_week0_dist)    
three3_week0_dist[1]
summary(three3_week0_pvalues)
summary(three3_week0_dist)

length(three3_week0_dist)
length(three3_week0_pvalues)
length( three3_week0_pvalues[three3_week0_pvalues>=0.05] )
length( three3_week0_pvalues[three3_week0_pvalues < 0.05] )
length( three3_week0_dist[three3_week0_dist < -0.08] )  ## three3_week0_EEDheto is less than three3_week0_EEDko
length( three3_week0_dist[three3_week0_dist >  0.08] )  ## three3_week0_EEDheto is more than three3_week0_EEDko


pValue_threshold=0.05
distance_threshold=0.08


## 差异不显著
three3_NOL_diffIndex_1 <- ( (three3_week0_pvalues>pValue_threshold)  &  (abs(three3_week0_dist)<distance_threshold) )
length(three3_NOL_diffIndex_1[three3_NOL_diffIndex_1])
write.table(x=three3_NOL_diffIndex_1,   file = paste(subdir_3_part5, "/1-three3-NoDiff.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")


## 差异显著  three3_week0_EEDheto is less than three3_week0_EEDko
three3_NOL_diffIndex_2 <- ( (three3_week0_pvalues<pValue_threshold)  &  (three3_week0_dist < -distance_threshold) )
length(three3_NOL_diffIndex_2[three3_NOL_diffIndex_2])
write.table(x=three3_NOL_diffIndex_2,   file = paste(subdir_3_part5, "/1-three3-EEDheto-Less-Than-EEDko.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")


## 差异显著  three3_week0_EEDheto is more than three3_week0_EEDko
three3_NOL_diffIndex_3 <- ( (three3_week0_pvalues<pValue_threshold)  &  (three3_week0_dist > distance_threshold) )
length(three3_NOL_diffIndex_3[three3_NOL_diffIndex_3])
write.table(x=three3_NOL_diffIndex_3,   file = paste(subdir_3_part5, "/1-three3-EEDheto-More-Than-EEDko.txt", sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE, qmethod = c("escape", "double"),  fileEncoding = "")






MyAverageLines_3(vector2=c(colMeans(Average_three3_week0_EEDheto[three3_NOL_diffIndex_1, ]),  colMeans(Average_three3_week0_EEDko[three3_NOL_diffIndex_1, ]), 
                           colMeans(Average_three3_week4_EEDheto[three3_NOL_diffIndex_1, ]),  colMeans(Average_three3_week4_EEDko[three3_NOL_diffIndex_1, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_3_part5,     fileName2="A-three3-NoDiff",  
                 title2="Up-regulated Genes (No Diff)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min,   Ymax2=col_NOL_max,    height2=3.2,   width2=5.55 , center2="TSS" )


three3_len_A <- nrow(Average_three3_week0_EEDheto[three3_NOL_diffIndex_1, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_three3_week0_EEDheto[three3_NOL_diffIndex_1, ]),  rowMeans(Average_three3_week0_EEDko[three3_NOL_diffIndex_1, ]), 
                            rowMeans(Average_three3_week4_EEDheto[three3_NOL_diffIndex_1, ]),  rowMeans(Average_three3_week4_EEDko[three3_NOL_diffIndex_1, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", three3_len_A),   rep("week0_EEDko",  three3_len_A), 
                                 rep("week4_EEDheto", three3_len_A),   rep("week4_EEDko",  three3_len_A) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part5,  fileName2="A-three3-all-NoDiff-BoxViolin",  
                  title2="Up-regulated Genes  (No Diff)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_three3_week0_EEDheto[three3_NOL_diffIndex_1, 201:400]),  rowMeans(Average_three3_week0_EEDko[three3_NOL_diffIndex_1, 201:400]), 
                            rowMeans(Average_three3_week4_EEDheto[three3_NOL_diffIndex_1, 201:400]),  rowMeans(Average_three3_week4_EEDko[three3_NOL_diffIndex_1, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", three3_len_A),   rep("week0_EEDko",  three3_len_A), 
                                 rep("week4_EEDheto", three3_len_A),   rep("week4_EEDko",  three3_len_A) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part5,  fileName2="A-three3-2kbRegion-NoDiff-BoxViolin",  
                  title2="Up-regulated Genes  (No Diff)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )










MyAverageLines_3(vector2=c(colMeans(Average_three3_week0_EEDheto[three3_NOL_diffIndex_2, ]),  colMeans(Average_three3_week0_EEDko[three3_NOL_diffIndex_2, ]), 
                           colMeans(Average_three3_week4_EEDheto[three3_NOL_diffIndex_2, ]),  colMeans(Average_three3_week4_EEDko[three3_NOL_diffIndex_2, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_3_part5,     fileName2="B-three3-EEDhetoLess",  
                 title2="Up-regulated Genes (EEDheto Less)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min,   Ymax2=col_NOL_max,    height2=3.2,   width2=5.55 , center2="TSS" )


three3_len_B <- nrow(Average_three3_week0_EEDheto[three3_NOL_diffIndex_2, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_three3_week0_EEDheto[three3_NOL_diffIndex_2, ]),  rowMeans(Average_three3_week0_EEDko[three3_NOL_diffIndex_2, ]), 
                            rowMeans(Average_three3_week4_EEDheto[three3_NOL_diffIndex_2, ]),  rowMeans(Average_three3_week4_EEDko[three3_NOL_diffIndex_2, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", three3_len_B),   rep("week0_EEDko",  three3_len_B), 
                                 rep("week4_EEDheto", three3_len_B),   rep("week4_EEDko",  three3_len_B) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part5,  fileName2="B-three3-all-EEDhetoLess-BoxViolin",  
                  title2="Up-regulated Genes  (EEDheto Less)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_three3_week0_EEDheto[three3_NOL_diffIndex_2, 201:400]),  rowMeans(Average_three3_week0_EEDko[three3_NOL_diffIndex_2, 201:400]), 
                            rowMeans(Average_three3_week4_EEDheto[three3_NOL_diffIndex_2, 201:400]),  rowMeans(Average_three3_week4_EEDko[three3_NOL_diffIndex_2, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", three3_len_B),   rep("week0_EEDko",  three3_len_B), 
                                 rep("week4_EEDheto", three3_len_B),   rep("week4_EEDko",  three3_len_B) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part5,  fileName2="B-three3-2kbRegion-EEDhetoLess-BoxViolin",  
                  title2="Up-regulated Genes  (EEDheto Less)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )












MyAverageLines_3(vector2=c(colMeans(Average_three3_week0_EEDheto[three3_NOL_diffIndex_3, ]),  colMeans(Average_three3_week0_EEDko[three3_NOL_diffIndex_3, ]), 
                           colMeans(Average_three3_week4_EEDheto[three3_NOL_diffIndex_3, ]),  colMeans(Average_three3_week4_EEDko[three3_NOL_diffIndex_3, ]) ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   
                                rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_3_part5,     fileName2="C-three3-EEDhetoMore",  
                 title2="Up-regulated Genes (EEDheto More)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min,   Ymax2=col_NOL_max,    height2=3.2,   width2=5.55 , center2="TSS" )


three3_len_C <- nrow(Average_three3_week0_EEDheto[three3_NOL_diffIndex_3, ])

MyBoxViolinPlot_1(vector2=c(rowMeans(Average_three3_week0_EEDheto[three3_NOL_diffIndex_3, ]),  rowMeans(Average_three3_week0_EEDko[three3_NOL_diffIndex_3, ]), 
                            rowMeans(Average_three3_week4_EEDheto[three3_NOL_diffIndex_3, ]),  rowMeans(Average_three3_week4_EEDko[three3_NOL_diffIndex_3, ]) ),  
                  sampleType2=c( rep("week0_EEDheto", three3_len_C),   rep("week0_EEDko",  three3_len_C), 
                                 rep("week4_EEDheto", three3_len_C),   rep("week4_EEDko",  three3_len_C) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part5,  fileName2="C-three3-all-EEDhetoMore-BoxViolin",  
                  title2="Up-regulated Genes  (EEDheto More)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


MyBoxViolinPlot_1(vector2=c(rowMeans(Average_three3_week0_EEDheto[three3_NOL_diffIndex_3, 201:400]),  rowMeans(Average_three3_week0_EEDko[three3_NOL_diffIndex_3, 201:400]), 
                            rowMeans(Average_three3_week4_EEDheto[three3_NOL_diffIndex_3, 201:400]),  rowMeans(Average_three3_week4_EEDko[three3_NOL_diffIndex_3, 201:400]) ),  
                  sampleType2=c( rep("week0_EEDheto", three3_len_C),   rep("week0_EEDko",  three3_len_C), 
                                 rep("week4_EEDheto", three3_len_C),   rep("week4_EEDko",  three3_len_C) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part5,  fileName2="C-three3-2kbRegion-EEDhetoMore-BoxViolin",  
                  title2="Up-regulated Genes  (EEDheto More)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1 )


















column_NTR_three3_EEDheto_same   <-  log( colMeans(Average_three3_week0_EEDheto[three3_NOL_diffIndex_1, ]) / colMeans(Average_three3_week4_EEDheto[three3_NOL_diffIndex_1, ]) )/4
column_NTR_three3_EEDko_same     <-  log( colMeans(Average_three3_week0_EEDko[three3_NOL_diffIndex_1, ])   / colMeans(Average_three3_week4_EEDko[three3_NOL_diffIndex_1, ])   )/4
column_NTR_three3_EEDheto_less   <-  log( colMeans(Average_three3_week0_EEDheto[three3_NOL_diffIndex_2, ]) / colMeans(Average_three3_week4_EEDheto[three3_NOL_diffIndex_2, ]) )/4
column_NTR_three3_EEDko_less     <-  log( colMeans(Average_three3_week0_EEDko[three3_NOL_diffIndex_2, ])   / colMeans(Average_three3_week4_EEDko[three3_NOL_diffIndex_2, ])   )/4
column_NTR_three3_EEDheto_more   <-  log( colMeans(Average_three3_week0_EEDheto[three3_NOL_diffIndex_3, ]) / colMeans(Average_three3_week4_EEDheto[three3_NOL_diffIndex_3, ]) )/4
column_NTR_three3_EEDko_more     <-  log( colMeans(Average_three3_week0_EEDko[three3_NOL_diffIndex_3, ])   / colMeans(Average_three3_week4_EEDko[three3_NOL_diffIndex_3, ])   )/4

row_NTR_three3_EEDheto_same   <-  log( row_Average_three3_week0_EEDheto[three3_NOL_diffIndex_1] / row_Average_three3_week4_EEDheto[three3_NOL_diffIndex_1] )/4
row_NTR_three3_EEDko_same     <-  log( row_Average_three3_week0_EEDko[three3_NOL_diffIndex_1]   / row_Average_three3_week4_EEDko[three3_NOL_diffIndex_1]   )/4
row_NTR_three3_EEDheto_less   <-  log( row_Average_three3_week0_EEDheto[three3_NOL_diffIndex_2] / row_Average_three3_week4_EEDheto[three3_NOL_diffIndex_2] )/4
row_NTR_three3_EEDko_less     <-  log( row_Average_three3_week0_EEDko[three3_NOL_diffIndex_2]   / row_Average_three3_week4_EEDko[three3_NOL_diffIndex_2]   )/4
row_NTR_three3_EEDheto_more   <-  log( row_Average_three3_week0_EEDheto[three3_NOL_diffIndex_3] / row_Average_three3_week4_EEDheto[three3_NOL_diffIndex_3] )/4
row_NTR_three3_EEDko_more     <-  log( row_Average_three3_week0_EEDko[three3_NOL_diffIndex_3]   / row_Average_three3_week4_EEDko[three3_NOL_diffIndex_3]   )/4


MyAverageLines_3(vector2=c(column_NTR_three3_EEDheto_same,  column_NTR_three3_EEDko_same,  column_NTR_three3_EEDheto_less, column_NTR_three3_EEDko_less,
                           column_NTR_three3_EEDheto_more, column_NTR_three3_EEDko_more ),   
                 numSample2=6,   
                 sampleType2=c( rep("EEDheto_same", numOfColumns1),   rep("EEDko_same", numOfColumns1),  rep("EEDheto_less", numOfColumns1),   rep("EEDko_less", numOfColumns1), 
                                rep("EEDheto_more", numOfColumns1),   rep("EEDko_more", numOfColumns1)  ), 
                 sampleRank2=c( "EEDheto_same",   "EEDko_same",  "EEDheto_less",  "EEDko_less",  "EEDheto_more",  "EEDko_more"  ),  
                 colours2=c("EEDheto_same"="red",   "EEDko_same"="red4",  "EEDheto_less"="blue",  "EEDko_less"="blue4",  "EEDheto_more"="green",  "EEDko_more"="green4"),   
                 path2=subdir_3_part5,     fileName2="1-three3-NTR",  
                 title2="Up-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR_three3_EEDheto_same,  column_NTR_three3_EEDko_same ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_same", numOfColumns1),   rep("EEDko_same", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_same",   "EEDko_same"  ),  
                 colours2=c("EEDheto_same"="red",   "EEDko_same"="red4" ),   
                 path2=subdir_3_part5,     fileName2="1-three3-NTR-same",  
                 title2="Up-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR_three3_EEDheto_less,  column_NTR_three3_EEDko_less ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_less", numOfColumns1),   rep("EEDko_less", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_less",   "EEDko_less"  ),  
                 colours2=c("EEDheto_less"="blue",   "EEDko_less"="blue4" ),   
                 path2=subdir_3_part5,     fileName2="1-three3-NTR-less",  
                 title2="Up-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyAverageLines_3(vector2=c(column_NTR_three3_EEDheto_more,  column_NTR_three3_EEDko_more ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto_more", numOfColumns1),   rep("EEDko_more", numOfColumns1) ), 
                 sampleRank2=c( "EEDheto_more",   "EEDko_more"  ),  
                 colours2=c("EEDheto_more"="blue",   "EEDko_more"="blue4" ),   
                 path2=subdir_3_part5,     fileName2="1-three3-NTR-more",  
                 title2="Up-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="NTR",   
                 Ymin2=-0.3,   Ymax2=0.8,    height2=3.2,   width2=5.55 , center2="TSS" )


MyBoxViolinPlot_1(vector2=c(row_NTR_three3_EEDheto_same,  row_NTR_three3_EEDko_same,  row_NTR_three3_EEDheto_less, row_NTR_three3_EEDko_less,
                            row_NTR_three3_EEDheto_more,  row_NTR_three3_EEDko_more  ),  
                  sampleType2=c( rep("EEDheto_same", three3_len_A),   rep("EEDko_same", three3_len_A),  rep("EEDheto_less", three3_len_B),   rep("EEDko_less", three3_len_B), 
                                 rep("EEDheto_more", three3_len_C),   rep("EEDko_more", three3_len_C)  ), 
                  sampleRank2=c( "EEDheto_same",   "EEDko_same",  "EEDheto_less",  "EEDko_less",  "EEDheto_more",  "EEDko_more"  ),  
                  colours2=c("red",  "red4",   "blue",    "blue4", "green", "green4" ),   
                  path2=subdir_3_part5,  fileName2="1-three3-NTR-boxViolin",  
                  title2="Up-regulated Genes",  xLab2="Samples",  yLab2="NTR",   
                  height2=3.88,   width2=4, Ymin2=-0.3, Ymax2=1 )











####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################








