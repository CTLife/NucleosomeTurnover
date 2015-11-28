#############################################################################################################################
## Part 3:  Figures about nucleosome occupancy level (NOL).
#############################################################################################################################


 



#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
subdir_1_part3 <- paste(Part3_g,  "/1-AveCol-AllRow", sep = "")
if( ! file.exists(subdir_1_part3) ) { dir.create(subdir_1_part3) }


MyAverageLines_3(vector2=c(column_one1_week0_EEDheto_rep1, column_one1_week0_EEDheto_rep2,  column_one1_week0_EEDko_rep1, 
                           column_one1_week4_EEDheto_rep1, column_one1_week4_EEDheto_rep2,  column_one1_week4_EEDko_rep1,  column_one1_week4_EEDko_rep2 ),   
             numSample2=7,   
             sampleType2=c( rep("week0_EEDheto_rep1", numOfColumns1),    rep("week0_EEDheto_rep2", numOfColumns1),  rep("week0_EEDko_rep1", numOfColumns1),  
                            rep("week4_EEDheto_rep1", numOfColumns1),    rep("week4_EEDheto_rep2", numOfColumns1),  rep("week4_EEDko_rep1",   numOfColumns1),    rep("week4_EEDko_rep2", numOfColumns1) ), 
             sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                            "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),     
             colours2=c( "week0_EEDheto_rep1"="red",   "week0_EEDheto_rep2"="red4",    "week0_EEDko_rep1"="purple", 
                         "week4_EEDheto_rep1"="blue",  "week4_EEDheto_rep2"="blue4",   "week4_EEDko_rep1"="green",   "week4_EEDko_rep2"="green4" ), 
             path2=subdir_1_part3,     fileName2="1A-one1-7Samples-curve",  
             title2="Down-regulated Genes",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
             Ymin2=col_NOL_min,   Ymax2=col_NOL_max,    height2=3.2,   width2=6.0,   center2="TSS" )

MyAverageLines_3(vector2=c(column_Average_one1_week0_EEDheto,  column_Average_one1_week0_EEDko, column_Average_one1_week4_EEDheto, column_Average_one1_week4_EEDko),   
             numSample2=4,   
             sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
             sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
             colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
             path2=subdir_1_part3,     fileName2="1A-one1-4Samples-curve",  
             title2="Down-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
             Ymin2=col_NOL_min,   Ymax2=col_NOL_max,    height2=3.2,   width2=5.55 , center2="TSS" )





MyAverageLines_3(vector2=c(column_two2_week0_EEDheto_rep1, column_two2_week0_EEDheto_rep2,  column_two2_week0_EEDko_rep1, 
                           column_two2_week4_EEDheto_rep1, column_two2_week4_EEDheto_rep2,  column_two2_week4_EEDko_rep1,  column_two2_week4_EEDko_rep2 ),   
                 numSample2=7,   
                 sampleType2=c( rep("week0_EEDheto_rep1", numOfColumns1),    rep("week0_EEDheto_rep2", numOfColumns1),  rep("week0_EEDko_rep1", numOfColumns1),  
                                rep("week4_EEDheto_rep1", numOfColumns1),    rep("week4_EEDheto_rep2", numOfColumns1), 
                                rep("week4_EEDko_rep1",   numOfColumns1),    rep("week4_EEDko_rep2", numOfColumns1) ), 
                 sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),     
                 colours2=c( "week0_EEDheto_rep1"="red",   "week0_EEDheto_rep2"="red4",    "week0_EEDko_rep1"="purple", 
                             "week4_EEDheto_rep1"="blue",  "week4_EEDheto_rep2"="blue4",   "week4_EEDko_rep1"="green",   "week4_EEDko_rep2"="green4" ), 
                 path2=subdir_1_part3,     fileName2="1B-two2-7Samples-curve",  
                 title2="Unchanged Genes",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min,   Ymax2=col_NOL_max,    height2=3.2,   width2=6.0,   center2="TSS" )

MyAverageLines_3(vector2=c(column_Average_two2_week0_EEDheto,  column_Average_two2_week0_EEDko, column_Average_two2_week4_EEDheto, column_Average_two2_week4_EEDko),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part3,     fileName2="1B-two2-4Samples-curve",  
                 title2="Unchanged Genes",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min,   Ymax2=col_NOL_max,    height2=3.2,   width2=5.55 , center2="TSS" )





MyAverageLines_3(vector2=c(column_three3_week0_EEDheto_rep1, column_three3_week0_EEDheto_rep2,  column_three3_week0_EEDko_rep1, 
                           column_three3_week4_EEDheto_rep1, column_three3_week4_EEDheto_rep2,  column_three3_week4_EEDko_rep1,  column_three3_week4_EEDko_rep2 ),   
                 numSample2=7,   
                 sampleType2=c( rep("week0_EEDheto_rep1", numOfColumns1),    rep("week0_EEDheto_rep2", numOfColumns1),  rep("week0_EEDko_rep1", numOfColumns1),  
                                rep("week4_EEDheto_rep1", numOfColumns1),    rep("week4_EEDheto_rep2", numOfColumns1), 
                                rep("week4_EEDko_rep1",   numOfColumns1),    rep("week4_EEDko_rep2", numOfColumns1) ), 
                 sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),     
                 colours2=c( "week0_EEDheto_rep1"="red",   "week0_EEDheto_rep2"="red4",    "week0_EEDko_rep1"="purple", 
                             "week4_EEDheto_rep1"="blue",  "week4_EEDheto_rep2"="blue4",   "week4_EEDko_rep1"="green",   "week4_EEDko_rep2"="green4" ), 
                 path2=subdir_1_part3,     fileName2="1C-three3-7Samples-curve",  
                 title2="Up-regulated Genes",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min,   Ymax2=col_NOL_max,    height2=3.2,   width2=6.0 , center2="TSS" )

MyAverageLines_3(vector2=c(column_Average_three3_week0_EEDheto,  column_Average_three3_week0_EEDko, column_Average_three3_week4_EEDheto, column_Average_three3_week4_EEDko),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part3,     fileName2="1C-three3-4Samples-curve",  
                 title2="Up-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min,   Ymax2=col_NOL_max,    height2=3.2,   width2=5.55 , center2="TSS" )







## divided by H3

MyAverageLines_3(vector2=c(column_one1_week0_EEDheto_rep1/column_one1_H3_normal_rep1, column_one1_week0_EEDheto_rep2/column_one1_H3_normal_rep1,  column_one1_week0_EEDko_rep1/column_one1_H3_normal_rep1, 
                           column_one1_week4_EEDheto_rep1/column_one1_H3_normal_rep1, column_one1_week4_EEDheto_rep2/column_one1_H3_normal_rep1,  column_one1_week4_EEDko_rep1/column_one1_H3_normal_rep1,  column_one1_week4_EEDko_rep2/column_one1_H3_normal_rep1 ),                                      
                 numSample2=7,   
                 sampleType2=c( rep("week0_EEDheto_rep1", numOfColumns1),    rep("week0_EEDheto_rep2", numOfColumns1),  rep("week0_EEDko_rep1", numOfColumns1),  
                                rep("week4_EEDheto_rep1", numOfColumns1),    rep("week4_EEDheto_rep2", numOfColumns1),  rep("week4_EEDko_rep1",   numOfColumns1),    rep("week4_EEDko_rep2", numOfColumns1) ), 
                 sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),     
                 colours2=c( "week0_EEDheto_rep1"="red",   "week0_EEDheto_rep2"="red4",    "week0_EEDko_rep1"="purple", 
                             "week4_EEDheto_rep1"="blue",  "week4_EEDheto_rep2"="blue4",   "week4_EEDko_rep1"="green",   "week4_EEDko_rep2"="green4" ), 
                 path2=subdir_1_part3,     fileName2="2A-one1-7Samples-curve",  
                 title2="Down-regulated Genes",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal/H3",   
                 Ymin2=0,   Ymax2=1.5,    height2=3.2,   width2=6.0,   center2="TSS" )

MyAverageLines_3(vector2=c(column_Average_one1_week0_EEDheto/column_Average_one1_H3_normal,  column_Average_one1_week0_EEDko/column_Average_one1_H3_normal, 
                           column_Average_one1_week4_EEDheto/column_Average_one1_H3_normal, column_Average_one1_week4_EEDko/column_Average_one1_H3_normal),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part3,     fileName2="2A-one1-4Samples-curve",  
                 title2="Down-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal/H3",   
                 Ymin2=0,   Ymax2=1.5,    height2=3.2,   width2=5.55 , center2="TSS" )





MyAverageLines_3(vector2=c(column_two2_week0_EEDheto_rep1/column_two2_H3_normal_rep1, column_two2_week0_EEDheto_rep2/column_two2_H3_normal_rep1,  column_two2_week0_EEDko_rep1/column_two2_H3_normal_rep1, 
                           column_two2_week4_EEDheto_rep1/column_two2_H3_normal_rep1, column_two2_week4_EEDheto_rep2/column_two2_H3_normal_rep1,  column_two2_week4_EEDko_rep1/column_two2_H3_normal_rep1,  column_two2_week4_EEDko_rep2/column_two2_H3_normal_rep1 ),   
                 numSample2=7,   
                 sampleType2=c( rep("week0_EEDheto_rep1", numOfColumns1),    rep("week0_EEDheto_rep2", numOfColumns1),  rep("week0_EEDko_rep1", numOfColumns1),  
                                rep("week4_EEDheto_rep1", numOfColumns1),    rep("week4_EEDheto_rep2", numOfColumns1), 
                                rep("week4_EEDko_rep1",   numOfColumns1),    rep("week4_EEDko_rep2", numOfColumns1) ), 
                 sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),     
                 colours2=c( "week0_EEDheto_rep1"="red",   "week0_EEDheto_rep2"="red4",    "week0_EEDko_rep1"="purple", 
                             "week4_EEDheto_rep1"="blue",  "week4_EEDheto_rep2"="blue4",   "week4_EEDko_rep1"="green",   "week4_EEDko_rep2"="green4" ), 
                 path2=subdir_1_part3,     fileName2="2B-two2-7Samples-curve",  
                 title2="Unchanged Genes",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal/H3",   
                 Ymin2=0,   Ymax2=1.5,     height2=3.2,   width2=6.0,   center2="TSS" )

MyAverageLines_3(vector2=c(column_Average_two2_week0_EEDheto/column_Average_two2_H3_normal,  column_Average_two2_week0_EEDko/column_Average_two2_H3_normal, 
                           column_Average_two2_week4_EEDheto/column_Average_two2_H3_normal, column_Average_two2_week4_EEDko/column_Average_two2_H3_normal),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part3,     fileName2="2B-two2-4Samples-curve",  
                 title2="Unchanged Genes",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal/H3",   
                 Ymin2=0,   Ymax2=1.5,     height2=3.2,   width2=5.55 , center2="TSS" )





MyAverageLines_3(vector2=c(column_three3_week0_EEDheto_rep1/column_three3_H3_normal_rep1, column_three3_week0_EEDheto_rep2/column_three3_H3_normal_rep1,  column_three3_week0_EEDko_rep1/column_three3_H3_normal_rep1, 
                           column_three3_week4_EEDheto_rep1/column_three3_H3_normal_rep1, column_three3_week4_EEDheto_rep2/column_three3_H3_normal_rep1,  column_three3_week4_EEDko_rep1/column_three3_H3_normal_rep1,  column_three3_week4_EEDko_rep2/column_three3_H3_normal_rep1 ),   
                 numSample2=7,   
                 sampleType2=c( rep("week0_EEDheto_rep1", numOfColumns1),    rep("week0_EEDheto_rep2", numOfColumns1),  rep("week0_EEDko_rep1", numOfColumns1),  
                                rep("week4_EEDheto_rep1", numOfColumns1),    rep("week4_EEDheto_rep2", numOfColumns1), 
                                rep("week4_EEDko_rep1",   numOfColumns1),    rep("week4_EEDko_rep2", numOfColumns1) ), 
                 sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),     
                 colours2=c( "week0_EEDheto_rep1"="red",   "week0_EEDheto_rep2"="red4",    "week0_EEDko_rep1"="purple", 
                             "week4_EEDheto_rep1"="blue",  "week4_EEDheto_rep2"="blue4",   "week4_EEDko_rep1"="green",   "week4_EEDko_rep2"="green4" ), 
                 path2=subdir_1_part3,     fileName2="2C-three3-7Samples-curve",  
                 title2="Up-regulated Genes",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal/H3",   
                 Ymin2=0,   Ymax2=1.5,    height2=3.2,   width2=6.0 , center2="TSS" )

MyAverageLines_3(vector2=c(column_Average_three3_week0_EEDheto/column_Average_three3_H3_normal,  column_Average_three3_week0_EEDko/column_Average_three3_H3_normal, 
                           column_Average_three3_week4_EEDheto/column_Average_three3_H3_normal, column_Average_three3_week4_EEDko/column_Average_three3_H3_normal ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part3,     fileName2="2C-three3-4Samples-curve",  
                 title2="Up-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal/H3",   
                 Ymin2=0,   Ymax2=1.5,     height2=3.2,   width2=5.55 , center2="TSS" )












############## smooth
MyAverageLines_3(vector2=c(column_one1_week0_EEDheto_rep1_smooth, column_one1_week0_EEDheto_rep2_smooth,  column_one1_week0_EEDko_rep1_smooth, 
                           column_one1_week4_EEDheto_rep1_smooth, column_one1_week4_EEDheto_rep2_smooth,  column_one1_week4_EEDko_rep1_smooth,  column_one1_week4_EEDko_rep2_smooth ),   
                 numSample2=7,   
                 sampleType2=c( rep("week0_EEDheto_rep1", numOfColumns1),    rep("week0_EEDheto_rep2", numOfColumns1),  rep("week0_EEDko_rep1", numOfColumns1),  
                                rep("week4_EEDheto_rep1", numOfColumns1),    rep("week4_EEDheto_rep2", numOfColumns1), 
                                rep("week4_EEDko_rep1",   numOfColumns1),    rep("week4_EEDko_rep2", numOfColumns1) ), 
                 sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),     
                 colours2=c( "week0_EEDheto_rep1"="red",   "week0_EEDheto_rep2"="red4",    "week0_EEDko_rep1"="purple", 
                             "week4_EEDheto_rep1"="blue",  "week4_EEDheto_rep2"="blue4",   "week4_EEDko_rep1"="green",   "week4_EEDko_rep2"="green4" ), 
                 path2=subdir_1_part3,     fileName2="2A-one1-7Samples-curve-smooth",  
                 title2="Down-regulated Genes",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min2,   Ymax2=col_NOL_max2,    height2=3.2,   width2=6.0 , center2="TSS" )

MyAverageLines_3(vector2=c(column_Average_one1_week0_EEDheto_smooth,  column_Average_one1_week0_EEDko_smooth, column_Average_one1_week4_EEDheto_smooth, column_Average_one1_week4_EEDko_smooth),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part3,     fileName2="2A-one1-4Samples-curve-smooth",  
                 title2="Down-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min2,   Ymax2=col_NOL_max2,    height2=3.2,   width2=5.55 , center2="TSS" )





MyAverageLines_3(vector2=c(column_two2_week0_EEDheto_rep1_smooth, column_two2_week0_EEDheto_rep2_smooth,  column_two2_week0_EEDko_rep1_smooth, 
                           column_two2_week4_EEDheto_rep1_smooth, column_two2_week4_EEDheto_rep2_smooth,  column_two2_week4_EEDko_rep1_smooth,  column_two2_week4_EEDko_rep2_smooth ),   
                 numSample2=7,   
                 sampleType2=c( rep("week0_EEDheto_rep1", numOfColumns1),    rep("week0_EEDheto_rep2", numOfColumns1),  rep("week0_EEDko_rep1", numOfColumns1),  
                                rep("week4_EEDheto_rep1", numOfColumns1),    rep("week4_EEDheto_rep2", numOfColumns1), 
                                rep("week4_EEDko_rep1",   numOfColumns1),    rep("week4_EEDko_rep2", numOfColumns1) ), 
                 sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2_smooth" ),     
                 colours2=c( "week0_EEDheto_rep1"="red",   "week0_EEDheto_rep2"="red4",    "week0_EEDko_rep1"="purple", 
                             "week4_EEDheto_rep1"="blue",  "week4_EEDheto_rep2"="blue4",   "week4_EEDko_rep1"="green",   "week4_EEDko_rep2"="green4" ), 
                 path2=subdir_1_part3,     fileName2="2B-two2-7Samples-curve-smooth",  
                 title2="Unchanged Genes",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min2,   Ymax2=col_NOL_max2,    height2=3.2,   width2=6.0 , center2="TSS" )

MyAverageLines_3(vector2=c(column_Average_two2_week0_EEDheto_smooth,  column_Average_two2_week0_EEDko_smooth, column_Average_two2_week4_EEDheto_smooth, column_Average_two2_week4_EEDko_smooth),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part3,     fileName2="2B-two2-4Samples-curve-smooth",  
                 title2="Unchanged Genes",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min2,   Ymax2=col_NOL_max2,    height2=3.2,   width2=5.55 , center2="TSS" )





MyAverageLines_3(vector2=c(column_three3_week0_EEDheto_rep1_smooth, column_three3_week0_EEDheto_rep2_smooth,  column_three3_week0_EEDko_rep1_smooth, 
                           column_three3_week4_EEDheto_rep1_smooth, column_three3_week4_EEDheto_rep2_smooth,  column_three3_week4_EEDko_rep1_smooth,  column_three3_week4_EEDko_rep2_smooth ),   
                 numSample2=7,   
                 sampleType2=c( rep("week0_EEDheto_rep1", numOfColumns1),    rep("week0_EEDheto_rep2", numOfColumns1),  rep("week0_EEDko_rep1", numOfColumns1),  
                                rep("week4_EEDheto_rep1", numOfColumns1),    rep("week4_EEDheto_rep2", numOfColumns1), 
                                rep("week4_EEDko_rep1",   numOfColumns1),    rep("week4_EEDko_rep2", numOfColumns1) ), 
                 sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2_smooth" ),     
                 colours2=c( "week0_EEDheto_rep1"="red",   "week0_EEDheto_rep2"="red4",    "week0_EEDko_rep1"="purple", 
                             "week4_EEDheto_rep1"="blue",  "week4_EEDheto_rep2"="blue4",   "week4_EEDko_rep1"="green",   "week4_EEDko_rep2"="green4" ), 
                 path2=subdir_1_part3,     fileName2="2C-three3-7Samples-curve-smooth",  
                 title2="Up-regulated Genes",     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min2,   Ymax2=col_NOL_max2,    height2=3.2,   width2=6.0 , center2="TSS" )

MyAverageLines_3(vector2=c(column_Average_three3_week0_EEDheto_smooth,  column_Average_three3_week0_EEDko_smooth, column_Average_three3_week4_EEDheto_smooth, column_Average_three3_week4_EEDko_smooth),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part3,     fileName2="2C-three3-4Samples-curve-smooth",  
                 title2="Up-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min2,   Ymax2=col_NOL_max2,    height2=3.2,   width2=5.55 , center2="TSS" )






 
 


## Error Bar
MyAverageLines_4(vector2=c(column_Average_one1_week0_EEDheto_smooth,  column_Average_one1_week0_EEDko_smooth, column_Average_one1_week4_EEDheto_smooth, column_Average_one1_week4_EEDko_smooth),   
                 SEM2 = c( SEM_one1_week0_EEDheto,  SEM_one1_week0_EEDko, SEM_one1_week4_EEDheto,  SEM_one1_week4_EEDko ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part3,     fileName2="3A-one1-4Samples-curve-smooth-SEM",  
                 title2="Down-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min2,   Ymax2=col_NOL_max2+0.05,    height2=3.2,   width2=5.55 , center2="TSS" )



MyAverageLines_4(vector2=c(column_Average_two2_week0_EEDheto_smooth,  column_Average_two2_week0_EEDko_smooth, column_Average_two2_week4_EEDheto_smooth, column_Average_two2_week4_EEDko_smooth),   
                 SEM2 = c( SEM_two2_week0_EEDheto,  SEM_two2_week0_EEDko, SEM_two2_week4_EEDheto,  SEM_two2_week4_EEDko ),
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part3,     fileName2="3B-two2-4Samples-curve-smooth-SEM",  
                 title2="Unchanged Genes",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min2,   Ymax2=col_NOL_max2+0.05,    height2=3.2,   width2=5.55 , center2="TSS" )



MyAverageLines_4(vector2=c(column_Average_three3_week0_EEDheto_smooth,  column_Average_three3_week0_EEDko_smooth, column_Average_three3_week4_EEDheto_smooth, column_Average_three3_week4_EEDko_smooth), 
                 SEM2 = c( SEM_three3_week0_EEDheto,  SEM_three3_week0_EEDko, SEM_three3_week4_EEDheto,  SEM_three3_week4_EEDko ),
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_1_part3,     fileName2="3C-three3-4Samples-curve-smooth-SEM",  
                 title2="Up-regulated Genes",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min2,   Ymax2=col_NOL_max2+0.05,    height2=3.2,   width2=5.55 , center2="TSS" )


















## hypothesis test
################################################################################################
subdir_2_part3 <- paste(Part3_g,  "/2-averageColumnsOfAllRows-HypTest", sep = "")
if( ! file.exists(subdir_2_part3) ) { dir.create(subdir_2_part3) }


MyHypothesisTest_2(vector1=column_Average_one1_week0_EEDheto, vector2=column_Average_one1_week0_EEDko, 
                   file1=paste(subdir_2_part3,  "/1A-column_one1_week0_EEDheto-VS-column_one1_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=column_Average_one1_week0_EEDheto, vector2=column_Average_one1_week4_EEDheto, 
                   file1=paste(subdir_2_part3,  "/1A-column_one1_week0_EEDheto-VS-column_one1_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=column_Average_one1_week0_EEDko, vector2=column_Average_one1_week4_EEDko, 
                   file1=paste(subdir_2_part3,  "/1A-column_one1_week0_EEDko-VS-column_one1_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=column_Average_one1_week4_EEDheto, vector2=column_Average_one1_week4_EEDko, 
                   file1=paste(subdir_2_part3,  "/1A-column_one1_week4_EEDheto-VS-column_one1_week4_EEDheto.txt", sep = "")  
)   

MyHypothesisTest_3(column_Average_one1_week0_EEDheto_smooth,  c(301:400),  file2=paste(subdir_2_part3,  "/1B-column_one1_week0_EEDheto", sep = "")  )   
MyHypothesisTest_3(column_Average_one1_week0_EEDko_smooth,    c(301:400),  file2=paste(subdir_2_part3,  "/1B-column_one1_week0_EEDko",   sep = "")  )  
MyHypothesisTest_3(column_Average_one1_week4_EEDheto_smooth,  c(301:400),  file2=paste(subdir_2_part3,  "/1B-column_one1_week4_EEDheto", sep = "")  )   
MyHypothesisTest_3(column_Average_one1_week4_EEDko_smooth,    c(301:400),  file2=paste(subdir_2_part3,  "/1B-column_one1_week4_EEDko",   sep = "")  )  







MyHypothesisTest_2(vector1=column_Average_two2_week0_EEDheto, vector2=column_Average_two2_week0_EEDko, 
                   file1=paste(subdir_2_part3,  "/2A-column_two2_week0_EEDheto-VS-column_two2_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=column_Average_two2_week0_EEDheto, vector2=column_Average_two2_week4_EEDheto, 
                   file1=paste(subdir_2_part3,  "/2A-column_two2_week0_EEDheto-VS-column_two2_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=column_Average_two2_week0_EEDko, vector2=column_Average_two2_week4_EEDko, 
                   file1=paste(subdir_2_part3,  "/2A-column_two2_week0_EEDko-VS-column_two2_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=column_Average_two2_week4_EEDheto, vector2=column_Average_two2_week4_EEDko, 
                   file1=paste(subdir_2_part3,  "/2A-column_two2_week4_EEDheto-VS-column_two2_week4_EEDheto.txt", sep = "")  
)   

MyHypothesisTest_3(column_Average_two2_week0_EEDheto_smooth,  c(301:400),  file2=paste(subdir_2_part3,  "/2B-column_two2_week0_EEDheto", sep = "")  )   
MyHypothesisTest_3(column_Average_two2_week0_EEDko_smooth,    c(301:400),  file2=paste(subdir_2_part3,  "/2B-column_two2_week0_EEDko",   sep = "")  )  
MyHypothesisTest_3(column_Average_two2_week4_EEDheto_smooth,  c(301:400),  file2=paste(subdir_2_part3,  "/2B-column_two2_week4_EEDheto", sep = "")  )   
MyHypothesisTest_3(column_Average_two2_week4_EEDko_smooth,    c(301:400),  file2=paste(subdir_2_part3,  "/2B-column_two2_week4_EEDko",   sep = "")  )  







MyHypothesisTest_2(vector1=column_Average_three3_week0_EEDheto, vector2=column_Average_three3_week0_EEDko, 
                   file1=paste(subdir_2_part3,  "/3A-column_three3_week0_EEDheto-VS-column_three3_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=column_Average_three3_week0_EEDheto, vector2=column_Average_three3_week4_EEDheto, 
                   file1=paste(subdir_2_part3,  "/3A-column_three3_week0_EEDheto-VS-column_three3_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=column_Average_three3_week0_EEDko, vector2=column_Average_three3_week4_EEDko, 
                   file1=paste(subdir_2_part3,  "/3A-column_three3_week0_EEDko-VS-column_three3_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=column_Average_three3_week4_EEDheto, vector2=column_Average_three3_week4_EEDko, 
                   file1=paste(subdir_2_part3,  "/3A-column_three3_week4_EEDheto-VS-column_three3_week4_EEDheto.txt", sep = "")  
)   

MyHypothesisTest_3(column_Average_three3_week0_EEDheto_smooth,  c(301:400),  file2=paste(subdir_2_part3,  "/3B-column_three3_week0_EEDheto", sep = "")  )   
MyHypothesisTest_3(column_Average_three3_week0_EEDko_smooth,    c(301:400),  file2=paste(subdir_2_part3,  "/3B-column_three3_week0_EEDko",   sep = "")  )  
MyHypothesisTest_3(column_Average_three3_week4_EEDheto_smooth,  c(301:400),  file2=paste(subdir_2_part3,  "/3B-column_three3_week4_EEDheto", sep = "")  )   
MyHypothesisTest_3(column_Average_three3_week4_EEDko_smooth,    c(301:400),  file2=paste(subdir_2_part3,  "/3B-column_three3_week4_EEDko",   sep = "")  )  















MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto, vector2=row_Average_one1_week0_EEDko, 
                   file1=paste(subdir_2_part3,  "/4-row_one1_week0_EEDheto-VS-row_one1_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto, vector2=row_Average_one1_week4_EEDheto, 
                   file1=paste(subdir_2_part3,  "/4-row_one1_week0_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDko, vector2=row_Average_one1_week4_EEDko, 
                   file1=paste(subdir_2_part3,  "/4-row_one1_week0_EEDko-VS-row_one1_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week4_EEDheto, vector2=row_Average_one1_week4_EEDko, 
                   file1=paste(subdir_2_part3,  "/4-row_one1_week4_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "")  
)   






MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto, vector2=row_Average_two2_week0_EEDko, 
                   file1=paste(subdir_2_part3,  "/5-row_two2_week0_EEDheto-VS-row_two2_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto, vector2=row_Average_two2_week4_EEDheto, 
                   file1=paste(subdir_2_part3,  "/5-row_two2_week0_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDko, vector2=row_Average_two2_week4_EEDko, 
                   file1=paste(subdir_2_part3,  "/5-row_two2_week0_EEDko-VS-row_two2_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week4_EEDheto, vector2=row_Average_two2_week4_EEDko, 
                   file1=paste(subdir_2_part3,  "/5-row_two2_week4_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "")  
)   





MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto, vector2=row_Average_three3_week0_EEDko, 
                   file1=paste(subdir_2_part3,  "/6-row_three3_week0_EEDheto-VS-row_three3_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto, vector2=row_Average_three3_week4_EEDheto, 
                   file1=paste(subdir_2_part3,  "/6-row_three3_week0_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDko, vector2=row_Average_three3_week4_EEDko, 
                   file1=paste(subdir_2_part3,  "/6-row_three3_week0_EEDko-VS-row_three3_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week4_EEDheto, vector2=row_Average_three3_week4_EEDko, 
                   file1=paste(subdir_2_part3,  "/6-row_three3_week4_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "")  
)   






























###############################################################################  Box-Violin-Plot
subdir_3_part3 <- paste(Part3_g,  "/3-AveRows-CenterCol", sep = "")
if( ! file.exists(subdir_3_part3) ) { dir.create(subdir_3_part3) }


MyBoxViolinPlot_1(vector2=c(row_one1_week0_EEDheto_rep1, row_one1_week0_EEDheto_rep2, row_one1_week0_EEDko_rep1,   
                            row_one1_week4_EEDheto_rep1, row_one1_week4_EEDheto_rep2, row_one1_week4_EEDko_rep1, row_one1_week4_EEDko_rep2 ),   
                  sampleType2=c( rep("week0_EEDheto_rep1", numOfRows1),  rep("week0_EEDheto_rep2", numOfRows1),    rep("week0_EEDko_rep1", numOfRows1),  
                                 rep("week4_EEDheto_rep1", numOfRows1),  rep("week4_EEDheto_rep2", numOfRows1),  
                                 rep("week4_EEDko_rep1", numOfRows1),  rep("week4_EEDko_rep2", numOfRows1)  ), 
                  sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                 "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),  
                  colours2=c( "week0_EEDheto_rep1"="red",  "week0_EEDheto_rep2"="red4",  "week0_EEDko_rep1"="purple", 
                              "week4_EEDheto_rep1"="blue",  "week4_EEDheto_rep2"="blue4",   "week4_EEDko_rep1"="green",   "week4_EEDko_rep2"="green4" ), 
                  path2=subdir_3_part3,   fileName2="1-one1-7Samples-BoxViolin",  
                  title2="Down-regulated Genes",   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4.1,   width2=4.2,   Ymin2=0, Ymax2=1.25)  

MyBoxViolinPlot_1(vector2=c(row_Average_one1_week0_EEDheto,  row_Average_one1_week0_EEDko, row_Average_one1_week4_EEDheto,  row_Average_one1_week4_EEDko ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1), rep("week4_EEDheto", numOfRows1),  rep("week4_EEDko", numOfRows1) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part3,  fileName2="1-one1-4Samples-BoxViolin",  
                  title2="Down-regulated Genes",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )





MyBoxViolinPlot_1(vector2=c(row_two2_week0_EEDheto_rep1, row_two2_week0_EEDheto_rep2, row_two2_week0_EEDko_rep1,   
                            row_two2_week4_EEDheto_rep1, row_two2_week4_EEDheto_rep2, row_two2_week4_EEDko_rep1, row_two2_week4_EEDko_rep2 ),   
                  sampleType2=c( rep("week0_EEDheto_rep1", numOfRows2),  rep("week0_EEDheto_rep2", numOfRows2),    rep("week0_EEDko_rep1", numOfRows2),  
                                 rep("week4_EEDheto_rep1", numOfRows2),  rep("week4_EEDheto_rep2", numOfRows2),    rep("week4_EEDko_rep1", numOfRows2),  rep("week4_EEDko_rep2", numOfRows2)  ), 
                  sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                 "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),  
                  colours2=c( "week0_EEDheto_rep1"="red",  "week0_EEDheto_rep2"="red4",  "week0_EEDko_rep1"="purple", 
                              "week4_EEDheto_rep1"="blue",  "week4_EEDheto_rep2"="blue4",   "week4_EEDko_rep1"="green",   "week4_EEDko_rep2"="green4" ), 
                  path2=subdir_3_part3,   fileName2="1-two2-7Samples-BoxViolin",  
                  title2="Unchanged Genes",   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4.1,   width2=4.2,   Ymin2=0, Ymax2=1.25) 

MyBoxViolinPlot_1(vector2=c(row_Average_two2_week0_EEDheto,  row_Average_two2_week0_EEDko, row_Average_two2_week4_EEDheto,  row_Average_two2_week4_EEDko ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows2),   rep("week0_EEDko", numOfRows2), rep("week4_EEDheto", numOfRows2),  rep("week4_EEDko", numOfRows2) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part3,  fileName2="1-two2-4Samples-BoxViolin",  
                  title2="Unchanged Genes",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )





MyBoxViolinPlot_1(vector2=c(row_three3_week0_EEDheto_rep1, row_three3_week0_EEDheto_rep2, row_three3_week0_EEDko_rep1,   
                            row_three3_week4_EEDheto_rep1, row_three3_week4_EEDheto_rep2, row_three3_week4_EEDko_rep1, row_three3_week4_EEDko_rep2 ),   
                  sampleType2=c( rep("week0_EEDheto_rep1", numOfRows3),  rep("week0_EEDheto_rep2", numOfRows3),    rep("week0_EEDko_rep1", numOfRows3),  
                                 rep("week4_EEDheto_rep1", numOfRows3),  rep("week4_EEDheto_rep2", numOfRows3),  
                                 rep("week4_EEDko_rep1", numOfRows3),  rep("week4_EEDko_rep2", numOfRows3)  ), 
                  sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                 "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),  
                  colours2=c( "week0_EEDheto_rep1"="red",  "week0_EEDheto_rep2"="red4",  "week0_EEDko_rep1"="purple", 
                              "week4_EEDheto_rep1"="blue",  "week4_EEDheto_rep2"="blue4",   "week4_EEDko_rep1"="green",   "week4_EEDko_rep2"="green4" ), 
                  path2=subdir_3_part3,   fileName2="1-three3-7Samples-BoxViolin",  
                  title2="Up-regulated Genes",   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4.1,   width2=4.2, Ymin2=0, Ymax2=1.25  )

MyBoxViolinPlot_1(vector2=c(row_Average_three3_week0_EEDheto,  row_Average_three3_week0_EEDko, row_Average_three3_week4_EEDheto,  row_Average_three3_week4_EEDko ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows3),   rep("week0_EEDko", numOfRows3), rep("week4_EEDheto", numOfRows3),  rep("week4_EEDko", numOfRows3) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part3,  fileName2="1-three3-4Samples-BoxViolin",  
                  title2="Up-regulated Genes",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8,  Ymin2=0, Ymax2=1.25 )





###############################################################################  Box-Violin-Scatter-Plot

MyBoxViolinPlot_2(vector2=c(row_one1_week0_EEDheto_rep1, row_one1_week0_EEDheto_rep2, row_one1_week0_EEDko_rep1,   
                            row_one1_week4_EEDheto_rep1, row_one1_week4_EEDheto_rep2, row_one1_week4_EEDko_rep1, row_one1_week4_EEDko_rep2 ),   
                  sampleType2=c( rep("week0_EEDheto_rep1", numOfRows1),  rep("week0_EEDheto_rep2", numOfRows1),    rep("week0_EEDko_rep1", numOfRows1),  
                                 rep("week4_EEDheto_rep1", numOfRows1),  rep("week4_EEDheto_rep2", numOfRows1),  
                                 rep("week4_EEDko_rep1", numOfRows1),  rep("week4_EEDko_rep2", numOfRows1)  ), 
                  sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                 "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),  
                  path2=subdir_3_part3,  fileName2="2-one1-7Samples-EEDKO-BoxViolinScatter",  
                  title2="Down-regulated genes",  
                  xLab2="Samples",  
                  yLab2="H2BGFP signal",   
                  height2=4.1,   width2=4.2,   Ymin2=0, Ymax2=1.25,  alpha2=0.8 )

MyBoxViolinPlot_2(vector2=c(row_Average_one1_week0_EEDheto,  row_Average_one1_week0_EEDko, row_Average_one1_week4_EEDheto,  row_Average_one1_week4_EEDko ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1),  
                                 rep("week4_EEDheto", numOfRows1),  rep("week4_EEDko", numOfRows1) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),     
                  path2=subdir_3_part3,  
                  fileName2="2-one1-4Samples-EEDKO-BoxViolinScatter",  
                  title2="Down-regulated genes",  
                  xLab2="Samples",  
                  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8,  Ymin2=0, Ymax2=1.25, alpha2=0.8 )





MyBoxViolinPlot_2(vector2=c(row_two2_week0_EEDheto_rep1, row_two2_week0_EEDheto_rep2, row_two2_week0_EEDko_rep1,   
                            row_two2_week4_EEDheto_rep1, row_two2_week4_EEDheto_rep2, row_two2_week4_EEDko_rep1, row_two2_week4_EEDko_rep2 ),   
                  sampleType2=c( rep("week0_EEDheto_rep1", numOfRows2),  rep("week0_EEDheto_rep2", numOfRows2),    rep("week0_EEDko_rep1", numOfRows2),  
                                 rep("week4_EEDheto_rep1", numOfRows2),  rep("week4_EEDheto_rep2", numOfRows2),  
                                 rep("week4_EEDko_rep1", numOfRows2),  rep("week4_EEDko_rep2", numOfRows2)  ), 
                  sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                 "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),  
                  path2=subdir_3_part3,  
                  fileName2="2-two2-7Samples-EEDKO-BoxViolin",  
                  title2="Unchanged genes",  
                  xLab2="Samples",  
                  yLab2="H2BGFP signal",   
                  height2=4.1,   width2=4.2,  Ymin2=0, Ymax2=1.25,  alpha2=0.2 )

MyBoxViolinPlot_2(vector2=c(row_Average_two2_week0_EEDheto,  row_Average_two2_week0_EEDko, row_Average_two2_week4_EEDheto,  row_Average_two2_week4_EEDko ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows2),   rep("week0_EEDko", numOfRows2),  
                                 rep("week4_EEDheto", numOfRows2),  rep("week4_EEDko", numOfRows2) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),     
                  path2=subdir_3_part3,  
                  fileName2="2-two2-4Samples-EEDKO-BoxViolin",  
                  title2="Unchanged genes",  
                  xLab2="Samples",  
                  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25,  alpha2=0.2 )






MyBoxViolinPlot_2(vector2=c(row_three3_week0_EEDheto_rep1, row_three3_week0_EEDheto_rep2, row_three3_week0_EEDko_rep1,   
                            row_three3_week4_EEDheto_rep1, row_three3_week4_EEDheto_rep2, row_three3_week4_EEDko_rep1, row_three3_week4_EEDko_rep2 ),   
                  sampleType2=c( rep("week0_EEDheto_rep1", numOfRows3),  rep("week0_EEDheto_rep2", numOfRows3),    rep("week0_EEDko_rep1", numOfRows3),  
                                 rep("week4_EEDheto_rep1", numOfRows3),  rep("week4_EEDheto_rep2", numOfRows3),  
                                 rep("week4_EEDko_rep1", numOfRows3),  rep("week4_EEDko_rep2", numOfRows3)  ), 
                  sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                 "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),  
                  path2=subdir_3_part3,  
                  fileName2="2-three3-7Samples-EEDKO-BoxViolin",  
                  title2="Up-regulated genes",  
                  xLab2="Samples",  
                  yLab2="H2BGFP signal",   
                  height2=4.1,   width2=4.2, Ymin2=0, Ymax2=1.5,  alpha2=0.8 )

MyBoxViolinPlot_2(vector2=c(row_Average_three3_week0_EEDheto,  row_Average_three3_week0_EEDko, row_Average_three3_week4_EEDheto,  row_Average_three3_week4_EEDko ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows3),   rep("week0_EEDko", numOfRows3),  
                                 rep("week4_EEDheto", numOfRows3),  rep("week4_EEDko", numOfRows3) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),     
                  path2=subdir_3_part3,  
                  fileName2="2-three3-4Samples-EEDKO-BoxViolin",  
                  title2="Up-regulated genes",  
                  xLab2="Samples",  
                  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.5,  alpha2=0.8 )







###############################################################################  with path

MyBoxViolinPlot_3_s7(vector2=c(row_one1_week0_EEDheto_rep1, row_one1_week0_EEDheto_rep2, row_one1_week0_EEDko_rep1,   
                            row_one1_week4_EEDheto_rep1, row_one1_week4_EEDheto_rep2, row_one1_week4_EEDko_rep1, row_one1_week4_EEDko_rep2 ), 
                  sampleType2=c( rep("week0_EEDheto_rep1", numOfRows1),  rep("week0_EEDheto_rep2", numOfRows1),    rep("week0_EEDko_rep1", numOfRows1),  
                                 rep("week4_EEDheto_rep1", numOfRows1),  rep("week4_EEDheto_rep2", numOfRows1),  
                                 rep("week4_EEDko_rep1", numOfRows1),  rep("week4_EEDko_rep2", numOfRows1)  ), 
                  sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                 "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),  
                  colours2=c( "week0_EEDheto_rep1"="red",  "week0_EEDheto_rep2"="red4",  "week0_EEDko_rep1"="purple", 
                              "week4_EEDheto_rep1"="blue",  "week4_EEDheto_rep2"="blue4",   "week4_EEDko_rep1"="green",   "week4_EEDko_rep2"="green4" ), 
                  path2=subdir_3_part3,   fileName2="3-one1-7Samples-BoxViolin",  
                  title2="Down-regulated Genes",   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4.1,   width2=5.5,   Ymin2=0, Ymax2=1.25)  

MyBoxViolinPlot_3_s4(vector2=c(row_Average_one1_week0_EEDheto,  row_Average_one1_week0_EEDko, row_Average_one1_week4_EEDheto,  row_Average_one1_week4_EEDko ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1), rep("week4_EEDheto", numOfRows1),  rep("week4_EEDko", numOfRows1) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part3,  fileName2="3-one1-4Samples-BoxViolin",  
                  title2="Down-regulated Genes",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=4, Ymin2=0, Ymax2=1.25 )








MyBoxViolinPlot_3_s7(vector2=c(row_two2_week0_EEDheto_rep1, row_two2_week0_EEDheto_rep2, row_two2_week0_EEDko_rep1,   
                            row_two2_week4_EEDheto_rep1, row_two2_week4_EEDheto_rep2, row_two2_week4_EEDko_rep1, row_two2_week4_EEDko_rep2 ),   
                  sampleType2=c( rep("week0_EEDheto_rep1", numOfRows2),  rep("week0_EEDheto_rep2", numOfRows2),    rep("week0_EEDko_rep1", numOfRows2),  
                                 rep("week4_EEDheto_rep1", numOfRows2),  rep("week4_EEDheto_rep2", numOfRows2),    rep("week4_EEDko_rep1", numOfRows2),  rep("week4_EEDko_rep2", numOfRows2)  ), 
                  sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                 "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),  
                  colours2=c( "week0_EEDheto_rep1"="red",  "week0_EEDheto_rep2"="red4",  "week0_EEDko_rep1"="purple", 
                              "week4_EEDheto_rep1"="blue",  "week4_EEDheto_rep2"="blue4",   "week4_EEDko_rep1"="green",   "week4_EEDko_rep2"="green4" ), 
                  path2=subdir_3_part3,   fileName2="3-two2-7Samples-BoxViolin",  
                  title2="Unchanged Genes",   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4.1,   width2=5.5,   Ymin2=0, Ymax2=1.25) 

MyBoxViolinPlot_3_s4(vector2=c(row_Average_two2_week0_EEDheto,  row_Average_two2_week0_EEDko, row_Average_two2_week4_EEDheto,  row_Average_two2_week4_EEDko ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows2),   rep("week0_EEDko", numOfRows2), rep("week4_EEDheto", numOfRows2),  rep("week4_EEDko", numOfRows2) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part3,  fileName2="3-two2-4Samples-BoxViolin",  
                  title2="Unchanged Genes",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=4, Ymin2=0, Ymax2=1.25 )








MyBoxViolinPlot_3_s7(vector2=c(row_three3_week0_EEDheto_rep1, row_three3_week0_EEDheto_rep2, row_three3_week0_EEDko_rep1,   
                            row_three3_week4_EEDheto_rep1, row_three3_week4_EEDheto_rep2, row_three3_week4_EEDko_rep1, row_three3_week4_EEDko_rep2 ),   
                  sampleType2=c( rep("week0_EEDheto_rep1", numOfRows3),  rep("week0_EEDheto_rep2", numOfRows3),    rep("week0_EEDko_rep1", numOfRows3),  
                                 rep("week4_EEDheto_rep1", numOfRows3),  rep("week4_EEDheto_rep2", numOfRows3),  
                                 rep("week4_EEDko_rep1", numOfRows3),  rep("week4_EEDko_rep2", numOfRows3)  ), 
                  sampleRank2=c( "week0_EEDheto_rep1",  "week0_EEDheto_rep2",  "week0_EEDko_rep1", 
                                 "week4_EEDheto_rep1",  "week4_EEDheto_rep2",   "week4_EEDko_rep1",   "week4_EEDko_rep2" ),  
                  colours2=c( "week0_EEDheto_rep1"="red",  "week0_EEDheto_rep2"="red4",  "week0_EEDko_rep1"="purple", 
                              "week4_EEDheto_rep1"="blue",  "week4_EEDheto_rep2"="blue4",   "week4_EEDko_rep1"="green",   "week4_EEDko_rep2"="green4" ), 
                  path2=subdir_3_part3,   fileName2="3-three3-7Samples-BoxViolin",  
                  title2="Up-regulated Genes",   xLab2="Samples",    yLab2="H2BGFP signal",   
                  height2=4.1,   width2=5.5, Ymin2=0, Ymax2=1.25  )

MyBoxViolinPlot_3_s4(vector2=c(row_Average_three3_week0_EEDheto,  row_Average_three3_week0_EEDko, row_Average_three3_week4_EEDheto,  row_Average_three3_week4_EEDko ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows3),   rep("week0_EEDko", numOfRows3), rep("week4_EEDheto", numOfRows3),  rep("week4_EEDko", numOfRows3) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_3_part3,  fileName2="3-three3-4Samples-BoxViolin",  
                  title2="Up-regulated Genes",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=4,  Ymin2=0, Ymax2=1.25 )





















###############################################################################
subdir_4_part3 <- paste(Part3_g,  "/4-averageRowsOfCenterColumns-Scatter", sep = "")
if( ! file.exists(subdir_4_part3) ) { dir.create(subdir_4_part3) }


MyScatterDiagram_1(vector2=row_Average_one1_week0_EEDheto, 
                   path2=subdir_4_part3,   
                   fileName2="1-one1-week0-EEDheto-scatter", 
                   xScale2=0.001,  xLab2="Genes (x1000)",   
                   yLab2="H2BGFP signal",  title2="Down-regulated Genes (week0_EEDheto)",  
                   height2=3.2,   width2=3.52, yMin2=0, yMax2=1.5, alpha2=1.0)


MyScatterDiagram_1(vector2=row_Average_one1_week0_EEDko, 
                   path2=subdir_4_part3,   
                   fileName2="1-one1-week0-EEDko-scatter", 
                   xScale2=0.001,  xLab2="Genes (x1000)",   
                   yLab2="H2BGFP signal",  title2="Down-regulated Genes (week0_EEDko)",  
                   height2=3.2,   width2=3.52, yMin2=0, yMax2=1.5, alpha2=1.0)


MyScatterDiagram_1(vector2=row_Average_one1_week4_EEDheto, 
                   path2=subdir_4_part3,   
                   fileName2="1-one1-week4-EEDheto-scatter", 
                   xScale2=0.001,  xLab2="Genes (x1000)",   
                   yLab2="H2BGFP signal",  title2="Down-regulated Genes (week4_EEDheto)",  
                   height2=3.2,   width2=3.52, yMin2=0, yMax2=1.5, alpha2=1.0)


MyScatterDiagram_1(vector2=row_Average_one1_week4_EEDko, 
                   path2=subdir_4_part3,   
                   fileName2="1-one1-week4-EEDko-scatter", 
                   xScale2=0.001,  xLab2="Genes (x1000)",   
                   yLab2="H2BGFP signal",  title2="Down-regulated Genes (week4_EEDko)",  
                   height2=3.2,   width2=3.52, yMin2=0, yMax2=1.5, alpha2=1.0)










MyScatterDiagram_1(vector2=row_Average_two2_week0_EEDheto, 
                   path2=subdir_4_part3,   
                   fileName2="2-two2-week0-EEDheto-scatter", 
                   xScale2=0.001,  xLab2="Genes (x1000)",   
                   yLab2="H2BGFP signal",  title2="Unchanged Genes (week0_EEDheto)",  
                   height2=3.2,   width2=3.52, yMin2=0, yMax2=1.5, alpha2=0.2)


MyScatterDiagram_1(vector2=row_Average_two2_week0_EEDko, 
                   path2=subdir_4_part3,   
                   fileName2="2-two2-week0-EEDko-scatter", 
                   xScale2=0.001,  xLab2="Genes (x1000)",   
                   yLab2="H2BGFP signal",  title2="Unchanged Genes (week0_EEDko)",  
                   height2=3.2,   width2=3.52, yMin2=0, yMax2=1.5, alpha2=0.2)


MyScatterDiagram_1(vector2=row_Average_two2_week4_EEDheto, 
                   path2=subdir_4_part3,   
                   fileName2="2-two2-week4-EEDheto-scatter", 
                   xScale2=0.001,  xLab2="Genes (x1000)",   
                   yLab2="H2BGFP signal",  title2="Unchanged Genes (week4_EEDheto)",  
                   height2=3.2,   width2=3.52, yMin2=0, yMax2=1.5, alpha2=0.2)


MyScatterDiagram_1(vector2=row_Average_two2_week4_EEDko, 
                   path2=subdir_4_part3,   
                   fileName2="2-two2-week4-EEDko-scatter", 
                   xScale2=0.001,  xLab2="Genes (x1000)",   
                   yLab2="H2BGFP signal",  title2="Unchanged Genes (week4_EEDko)",  
                   height2=3.2,   width2=3.52, yMin2=0, yMax2=1.5, alpha2=0.2)







MyScatterDiagram_1(vector2=row_Average_three3_week0_EEDheto, 
                   path2=subdir_4_part3,   
                   fileName2="3-three3-week0-EEDheto-scatter", 
                   xScale2=0.001,  xLab2="Genes (x1000)",   
                   yLab2="H2BGFP signal",  title2="Up-regulated Genes (week0_EEDheto)",  
                   height2=3.2,   width2=3.52, yMin2=0, yMax2=1.5, alpha2=1.0)


MyScatterDiagram_1(vector2=row_Average_three3_week0_EEDko, 
                   path2=subdir_4_part3,   
                   fileName2="3-three3-week0-EEDko-scatter", 
                   xScale2=0.001,  xLab2="Genes (x1000)",   
                   yLab2="H2BGFP signal",  title2="Up-regulated Genes (week0_EEDko)",  
                   height2=3.2,   width2=3.52, yMin2=0, yMax2=1.5, alpha2=1.0)


MyScatterDiagram_1(vector2=row_Average_three3_week4_EEDheto, 
                   path2=subdir_4_part3,   
                   fileName2="3-three3-week4-EEDheto-scatter", 
                   xScale2=0.001,  xLab2="Genes (x1000)",   
                   yLab2="H2BGFP signal",  title2="Up-regulated Genes (week4_EEDheto)",  
                   height2=3.2,   width2=3.52, yMin2=0, yMax2=1.5, alpha2=1.0)


MyScatterDiagram_1(vector2=row_Average_three3_week4_EEDko, 
                   path2=subdir_4_part3,   
                   fileName2="3-three3-week4-EEDko-scatter", 
                   xScale2=0.001,  xLab2="Genes (x1000)",   
                   yLab2="H2BGFP signal",  title2="Up-regulated Genes (week4_EEDko)",  
                   height2=3.2,   width2=3.52, yMin2=0, yMax2=1.5, alpha2=1.0)














###############################################################################
subdir_5_part3 <- paste(Part3_g,  "/5-averageRowsOfCenterColumns-histogram", sep = "")
if( ! file.exists(subdir_5_part3) ) { dir.create(subdir_5_part3) }



MyHistogram_1(vector2=row_Average_one1_week0_EEDheto,  
              path2=subdir_5_part3,    
              fileName2="1-one1-week0-EEDheto-histogram",  
              title2="Down-regulated Genes (week0 EEDheto)",  
              xLab2="H2BGFP signal",   
              height2=3.2,  width2=3.5,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=0.25)

MyHistogram_1(vector2=row_Average_one1_week0_EEDko,  
              path2=subdir_5_part3,    
              fileName2="1-one1-week0-EEDko-histogram",  
              title2="Down-regulated Genes (week0 EEDko)",  
              xLab2="H2BGFP signal",   
              height2=3.2,  width2=3.5,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=0.25)

MyHistogram_1(vector2=row_Average_one1_week4_EEDheto,  
              path2=subdir_5_part3,    
              fileName2="1-one1-week4-EEDheto-histogram",  
              title2="Down-regulated Genes (week4 EEDheto)",  
              xLab2="H2BGFP signal",   
              height2=3.2,  width2=3.5,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=0.25)

MyHistogram_1(vector2=row_Average_one1_week4_EEDko,  
              path2=subdir_5_part3,    
              fileName2="1-one1-week4-EEDko-histogram",  
              title2="Down-regulated Genes (week4 EEDko)",  
              xLab2="H2BGFP signal",   
              height2=3.2,  width2=3.5,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=0.25)







MyHistogram_1(vector2=row_Average_two2_week0_EEDheto,  
            path2=subdir_5_part3,    
            fileName2="2-two2-week0-EEDheto-histogram",  
            title2="Unchanged Genes (week0 EEDheto)",  
            xLab2="H2BGFP signal",   
            height2=3.2,  width2=3.5,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=0.25)

MyHistogram_1(vector2=row_Average_two2_week0_EEDko,  
              path2=subdir_5_part3,    
              fileName2="2-two2-week0-EEDko-histogram",  
              title2="Unchanged Genes (week0 EEDko)",  
              xLab2="H2BGFP signal",   
              height2=3.2,  width2=3.5,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=0.25)

MyHistogram_1(vector2=row_Average_two2_week4_EEDheto,  
              path2=subdir_5_part3,    
              fileName2="2-two2-week4-EEDheto-histogram",  
              title2="Unchanged Genes (week4 EEDheto)",  
              xLab2="H2BGFP signal",   
              height2=3.2,  width2=3.5,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=0.25)

MyHistogram_1(vector2=row_Average_two2_week4_EEDko,  
              path2=subdir_5_part3,    
              fileName2="2-two2-week4-EEDko-histogram",  
              title2="Unchanged Genes (week4 EEDko)",  
              xLab2="H2BGFP signal",   
              height2=3.2,  width2=3.5,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=0.25)


  






MyHistogram_1(vector2=row_Average_three3_week0_EEDheto,  
              path2=subdir_5_part3,    
              fileName2="3-three3-week0-EEDheto-histogram",  
              title2="Up-regulated genes (week0 EEDheto)",  
              xLab2="H2BGFP signal",   
              height2=3.2,  width2=3.5,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=0.25)

MyHistogram_1(vector2=row_Average_three3_week0_EEDko,  
              path2=subdir_5_part3,    
              fileName2="3-three3-week0-EEDko-histogram",  
              title2="Up-regulated genes (week0 EEDko)",  
              xLab2="H2BGFP signal",   
              height2=3.2,  width2=3.5,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=0.25)

MyHistogram_1(vector2=row_Average_three3_week4_EEDheto,  
              path2=subdir_5_part3,    
              fileName2="3-three3-week4-EEDheto-histogram",  
              title2="Up-regulated genes (week4 EEDheto)",  
              xLab2="H2BGFP signal",   
              height2=3.2,  width2=3.5,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=0.25)

MyHistogram_1(vector2=row_Average_three3_week4_EEDko,  
              path2=subdir_5_part3,    
              fileName2="3-three3-week4-EEDko-histogram",  
              title2="Up-regulated genes (week4 EEDko)",  
              xLab2="H2BGFP signal",   
              height2=3.2,  width2=3.5,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=0.25)










MyHistogram_2(vector2=row_Average_three3_week4_EEDko,  
              path2=subdir_5_part3,    
              fileName2="4-three3-week4-EEDko-histogram",  
              title2="Up-regulated genes (week4 EEDko)",  
              xLab2="H2BGFP signal",   
              height2=3.2,  width2=3.5,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=100)





MyHistogram_3(vector2=row_Average_three3_week0_EEDko,  
              path2=subdir_5_part3,    
              fileName2="7-three3-week0-EEDko-histogram",  
              title2="Up-regulated genes (week0 EEDko)",  
              xLab2="H2BGFP signal",   
              height2=3.2,  width2=4.3,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=10)


MyHistogram_3(vector2=row_Average_three3_week4_EEDko,  
              path2=subdir_5_part3,    
              fileName2="7-three3-week4-EEDko-histogram",  
              title2="Up-regulated genes (week4 EEDko)",  
              xLab2="H2BGFP signal",   
              height2=3.2,  width2=4.3,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=10)







###############################################################################
subdir_6_part3 <- paste(Part3_g,  "/6-averageRowsOfCenterColumns-Cmp", sep = "")
if( ! file.exists(subdir_6_part3) ) { dir.create(subdir_6_part3) }

MyHistogram_4(vector2=c(row_Average_three3_week0_EEDheto,  row_Average_three3_week0_EEDko, row_Average_three3_week4_EEDheto,  row_Average_three3_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows3),   rep("week0_EEDko", numOfRows3),  
                             rep("week4_EEDheto", numOfRows3),   rep("week4_EEDko", numOfRows3) ),  
              path2=subdir_6_part3,    
              fileName2="1-three3-4Samples-EEDKO-BoxViolin",   
              title2="Up-regulated genes",    xLab2="H2BGFP signal",  
              height2=3.2,  width2=5.6,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=100,  alpha2=1)





MyHistogram_5(vector2=c(row_Average_three3_week0_EEDheto,  row_Average_three3_week0_EEDko, row_Average_three3_week4_EEDheto,  row_Average_three3_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows3),   rep("week0_EEDko", numOfRows3),  
                             rep("week4_EEDheto", numOfRows3),   rep("week4_EEDko", numOfRows3) ),  
              path2=subdir_6_part3,    
              fileName2="2-three3-4Samples-EEDKO-BoxViolin",   
              title2="Up-regulated genes",    xLab2="H2BGFP signal",  
              height2=3.2,  width2=5.5,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=100)

















MyHistogram_6(vector2=c(row_one1_week0_EEDheto_rep1, row_one1_week0_EEDheto_rep2, row_one1_week0_EEDko_rep1,   
                        row_one1_week4_EEDheto_rep1, row_one1_week4_EEDheto_rep2, row_one1_week4_EEDko_rep1, row_one1_week4_EEDko_rep2 ),   
               sampleType2=c( rep("week0_EEDheto_rep1", numOfRows1),  rep("week0_EEDheto_rep2", numOfRows1),    rep("week0_EEDko_rep1", numOfRows1),  
                             rep("week4_EEDheto_rep1", numOfRows1),  rep("week4_EEDheto_rep2", numOfRows1),  
                             rep("week4_EEDko_rep1", numOfRows1),  rep("week4_EEDko_rep2", numOfRows1)  ), 
               colours2=c("red",  "red4",  "purple",   "blue",  "blue4", "green", "green4"),   
               path2=subdir_6_part3,    
               fileName2="7-one1-7Samples-BoxViolin",   
               title2="Down-regulated Genes",    xLab2="H2BGFP signal",  
               height2=3.2,  width2=6,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=5)

MyHistogram_6(vector2=c(row_Average_one1_week0_EEDheto,  row_Average_one1_week0_EEDko, row_Average_one1_week4_EEDheto,  row_Average_one1_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1),  
                             rep("week4_EEDheto", numOfRows1),   rep("week4_EEDko", numOfRows1) ),  
              colours2=c("week0_EEDheto"="red",  "week0_EEDko"="orange",  "week4_EEDheto"="blue",  "week4_EEDko"="green"),   
              path2=subdir_6_part3,    
              fileName2="7-one1-4Samples-BoxViolin",   
              title2="Down-regulated Genes",    xLab2="H2BGFP signal",  
              height2=3.2,  width2=5.6,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=5)





MyHistogram_6(vector2=c(row_two2_week0_EEDheto_rep1, row_two2_week0_EEDheto_rep2, row_two2_week0_EEDko_rep1,   
                        row_two2_week4_EEDheto_rep1, row_two2_week4_EEDheto_rep2, row_two2_week4_EEDko_rep1, row_two2_week4_EEDko_rep2 ),   
              sampleType2=c( rep("week0_EEDheto_rep1", numOfRows2),  rep("week0_EEDheto_rep2", numOfRows2),    rep("week0_EEDko_rep1", numOfRows2),  
                             rep("week4_EEDheto_rep1", numOfRows2),  rep("week4_EEDheto_rep2", numOfRows2),  
                             rep("week4_EEDko_rep1", numOfRows2),  rep("week4_EEDko_rep2", numOfRows2)  ), 
              colours2=c("red",  "red4",  "purple",   "blue",  "blue4", "green", "green4"),   
              path2=subdir_6_part3,    
              fileName2="8-two2-7Samples-BoxViolin",   
              title2="Unchanged Genes",    xLab2="H2BGFP signal",  
              height2=3.2,  width2=6,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=5)

MyHistogram_6(vector2=c(row_Average_two2_week0_EEDheto,  row_Average_two2_week0_EEDko, row_Average_two2_week4_EEDheto,  row_Average_two2_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows2),   rep("week0_EEDko", numOfRows2),  
                             rep("week4_EEDheto", numOfRows2),   rep("week4_EEDko", numOfRows2) ),  
              colours2=c("week0_EEDheto"="red",  "week0_EEDko"="orange",  "week4_EEDheto"="blue",  "week4_EEDko"="green"),   
              path2=subdir_6_part3,    
              fileName2="8-two2-4Samples-BoxViolin",   
              title2="Unchanged Genes",    xLab2="H2BGFP signal",  
              height2=3.2,  width2=5.6,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=5)





MyHistogram_6(vector2=c(row_three3_week0_EEDheto_rep1, row_three3_week0_EEDheto_rep2, row_three3_week0_EEDko_rep1,   
                        row_three3_week4_EEDheto_rep1, row_three3_week4_EEDheto_rep2, row_three3_week4_EEDko_rep1, row_three3_week4_EEDko_rep2 ),   
              sampleType2=c( rep("week0_EEDheto_rep1", numOfRows3),  rep("week0_EEDheto_rep2", numOfRows3),    rep("week0_EEDko_rep1", numOfRows3),  
                             rep("week4_EEDheto_rep1", numOfRows3),  rep("week4_EEDheto_rep2", numOfRows3),  
                             rep("week4_EEDko_rep1", numOfRows3),  rep("week4_EEDko_rep2", numOfRows3)  ), 
              colours2=c("red",  "red4",  "purple",   "blue",  "blue4", "green", "green4"),   
              path2=subdir_6_part3,    
              fileName2="9-three3-7Samples-BoxViolin",   
              title2="Up-regulated genes",    xLab2="H2BGFP signal",  
              height2=3.2,  width2=6,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=5)

MyHistogram_6(vector2=c(row_Average_three3_week0_EEDheto,  row_Average_three3_week0_EEDko, row_Average_three3_week4_EEDheto,  row_Average_three3_week4_EEDko ),  
               sampleType2=c( rep("week0_EEDheto", numOfRows3),   rep("week0_EEDko", numOfRows3),  
                             rep("week4_EEDheto", numOfRows3),   rep("week4_EEDko", numOfRows3) ),  
               colours2=c("week0_EEDheto"="red",  "week0_EEDko"="orange",  "week4_EEDheto"="blue",  "week4_EEDko"="green"),   
               path2=subdir_6_part3,    
               fileName2="9-three3-4Samples-BoxViolin",   
               title2="Up-regulated genes",    xLab2="H2BGFP signal",  
               height2=3.2,  width2=5.6,   xMin2=0,  xMax2=1.2,  yMin2=0,  yMax2=5)





  




MyHistogram_7(vector2=c(row_Average_three3_week0_EEDheto,  row_Average_three3_week0_EEDko, row_Average_three3_week4_EEDheto,  row_Average_three3_week4_EEDko ),  
              sampleType2=c( rep("week0_EEDheto", numOfRows3),   rep("week0_EEDko", numOfRows3),  
                             rep("week4_EEDheto", numOfRows3),   rep("week4_EEDko", numOfRows3) ),  
              colours2=c("week0_EEDheto"="red",  "week0_EEDko"="orange",  "week4_EEDheto"="blue",  "week4_EEDko"="green"),   
              path2=subdir_6_part3,    
              fileName2="10-three3-4Samples-BoxViolin",   
              title2="Up-regulated genes",    xLab2="H2BGFP signal",  
              height2=3.2,  width2=5.6,   xMin2=0,  xMax2=0.6,  yMin2=0,  yMax2=5)


















###############################################################################
subdir_7_part3 <- paste(Part3_g,  "/7-rows5Classes-curve", sep = "")
if( ! file.exists(subdir_7_part3) ) { dir.create(subdir_7_part3) }





dim(reduceRow_Average_one1_week0_EEDheto) 
dim(reduceRow_Average_one1_week0_EEDko)   
dim(reduceRow_Average_one1_week4_EEDheto) 
dim(reduceRow_Average_one1_week4_EEDko)    

MyAverageLines_3(vector2=  c( colMeans(reduceRow_Average_one1_week0_EEDheto[c(1,2), ]),   colMeans(reduceRow_Average_one1_week0_EEDko[c(1,2),]),  
                              colMeans(reduceRow_Average_one1_week4_EEDheto[c(1,2), ]),   colMeans(reduceRow_Average_one1_week4_EEDko[c(1,2),] ) ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_7_part3,     fileName2="1A-one1-4Samples-curve",  
                 title2="Down-regulated Genes (0%~20%)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min-0.1,   Ymax2=col_NOL_max+0.1,    height2=3.2,   width2=5.65 , center2="TSS" )

MyAverageLines_3(vector2=  c( colMeans(reduceRow_Average_one1_week0_EEDheto[c(3,4), ]),   colMeans(reduceRow_Average_one1_week0_EEDko[c(3,4),]),  
                              colMeans(reduceRow_Average_one1_week4_EEDheto[c(3,4), ]),   colMeans(reduceRow_Average_one1_week4_EEDko[c(3,4),] ) ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_7_part3,     fileName2="1B-one1-4Samples-curve",  
                 title2="Down-regulated Genes (20%~40%)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min-0.1,   Ymax2=col_NOL_max+0.1,    height2=3.2,   width2=5.65 , center2="TSS" )

MyAverageLines_3(vector2=  c( colMeans(reduceRow_Average_one1_week0_EEDheto[c(5,6), ]),   colMeans(reduceRow_Average_one1_week0_EEDko[c(5,6),]),  
                              colMeans(reduceRow_Average_one1_week4_EEDheto[c(5,6), ]),   colMeans(reduceRow_Average_one1_week4_EEDko[c(5,6),] ) ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_7_part3,     fileName2="1C-one1-4Samples-curve",  
                 title2="Down-regulated Genes (40%~60%)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min-0.1,   Ymax2=col_NOL_max+0.1,    height2=3.2,   width2=5.65 , center2="TSS" )

MyAverageLines_3(vector2=  c( colMeans(reduceRow_Average_one1_week0_EEDheto[c(7,8), ]),   colMeans(reduceRow_Average_one1_week0_EEDko[c(7,8),]),  
                              colMeans(reduceRow_Average_one1_week4_EEDheto[c(7,8), ]),   colMeans(reduceRow_Average_one1_week4_EEDko[c(7,8),] ) ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_7_part3,     fileName2="1D-one1-4Samples-curve",  
                 title2="Down-regulated Genes (60%~80%)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min-0.1,   Ymax2=col_NOL_max+0.1,    height2=3.2,   width2=5.65 , center2="TSS" )


MyAverageLines_3(vector2=  c( colMeans(reduceRow_Average_one1_week0_EEDheto[c(9,10), ]),   colMeans(reduceRow_Average_one1_week0_EEDko[c(9,10),]),  
                              colMeans(reduceRow_Average_one1_week4_EEDheto[c(9,10), ]),   colMeans(reduceRow_Average_one1_week4_EEDko[c(9,10),] ) ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_7_part3,     fileName2="1E-one1-4Samples-curve",  
                 title2="Down-regulated Genes (80%~100%)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min-0.1,   Ymax2=col_NOL_max+0.1,    height2=3.2,   width2=5.65 , center2="TSS" )













dim(reduceRow_Average_two2_week0_EEDheto) 
dim(reduceRow_Average_two2_week0_EEDko)   
dim(reduceRow_Average_two2_week4_EEDheto) 
dim(reduceRow_Average_two2_week4_EEDko)    

MyAverageLines_3(vector2=  c( colMeans(reduceRow_Average_two2_week0_EEDheto[c(1,2), ]),   colMeans(reduceRow_Average_two2_week0_EEDko[c(1,2),]),  
                              colMeans(reduceRow_Average_two2_week4_EEDheto[c(1,2), ]),   colMeans(reduceRow_Average_two2_week4_EEDko[c(1,2),] ) ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_7_part3,     fileName2="2A-two2-4Samples-curve",  
                 title2="Unchanged Genes (0%~20%)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min-0.1,   Ymax2=col_NOL_max+0.1,    height2=3.2,   width2=5.65 , center2="TSS" )

MyAverageLines_3(vector2=  c( colMeans(reduceRow_Average_two2_week0_EEDheto[c(3,4), ]),   colMeans(reduceRow_Average_two2_week0_EEDko[c(3,4),]),  
                              colMeans(reduceRow_Average_two2_week4_EEDheto[c(3,4), ]),   colMeans(reduceRow_Average_two2_week4_EEDko[c(3,4),] ) ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_7_part3,     fileName2="2B-two2-4Samples-curve",  
                 title2="Unchanged Genes (20%~40%)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min-0.1,   Ymax2=col_NOL_max+0.1,    height2=3.2,   width2=5.65 , center2="TSS" )

MyAverageLines_3(vector2=  c( colMeans(reduceRow_Average_two2_week0_EEDheto[c(5,6), ]),   colMeans(reduceRow_Average_two2_week0_EEDko[c(5,6),]),  
                              colMeans(reduceRow_Average_two2_week4_EEDheto[c(5,6), ]),   colMeans(reduceRow_Average_two2_week4_EEDko[c(5,6),] ) ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_7_part3,     fileName2="2C-two2-4Samples-curve",  
                 title2="Unchanged Genes (40%~60%)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min-0.1,   Ymax2=col_NOL_max+0.1,    height2=3.2,   width2=5.65 , center2="TSS" )

MyAverageLines_3(vector2=  c( colMeans(reduceRow_Average_two2_week0_EEDheto[c(7,8), ]),   colMeans(reduceRow_Average_two2_week0_EEDko[c(7,8),]),  
                              colMeans(reduceRow_Average_two2_week4_EEDheto[c(7,8), ]),   colMeans(reduceRow_Average_two2_week4_EEDko[c(7,8),] ) ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_7_part3,     fileName2="2D-two2-4Samples-curve",  
                 title2="Unchanged Genes (60%~80%)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min-0.1,   Ymax2=col_NOL_max+0.1,    height2=3.2,   width2=5.65 , center2="TSS" )


MyAverageLines_3(vector2=  c( colMeans(reduceRow_Average_two2_week0_EEDheto[c(9,10), ]),   colMeans(reduceRow_Average_two2_week0_EEDko[c(9,10),]),  
                              colMeans(reduceRow_Average_two2_week4_EEDheto[c(9,10), ]),   colMeans(reduceRow_Average_two2_week4_EEDko[c(9,10),] ) ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_7_part3,     fileName2="2E-two2-4Samples-curve",  
                 title2="Unchanged Genes (80%~100%)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min-0.1,   Ymax2=col_NOL_max+0.1,    height2=3.2,   width2=5.65 , center2="TSS" )















dim(reduceRow_Average_three3_week0_EEDheto) 
dim(reduceRow_Average_three3_week0_EEDko)   
dim(reduceRow_Average_three3_week4_EEDheto) 
dim(reduceRow_Average_three3_week4_EEDko)    

MyAverageLines_3(vector2=  c( colMeans(reduceRow_Average_three3_week0_EEDheto[c(1,2), ]),   colMeans(reduceRow_Average_three3_week0_EEDko[c(1,2),]),  
                              colMeans(reduceRow_Average_three3_week4_EEDheto[c(1,2), ]),   colMeans(reduceRow_Average_three3_week4_EEDko[c(1,2),] ) ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_7_part3,     fileName2="3A-three3-4Samples-curve",  
                 title2="Up-regulated Genes (0%~20%)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min-0.1,   Ymax2=col_NOL_max+0.1,    height2=3.2,   width2=5.65 , center2="TSS" )

MyAverageLines_3(vector2=  c( colMeans(reduceRow_Average_three3_week0_EEDheto[c(3,4), ]),   colMeans(reduceRow_Average_three3_week0_EEDko[c(3,4),]),  
                              colMeans(reduceRow_Average_three3_week4_EEDheto[c(3,4), ]),   colMeans(reduceRow_Average_three3_week4_EEDko[c(3,4),] ) ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_7_part3,     fileName2="3B-three3-4Samples-curve",  
                 title2="Up-regulated Genes (20%~40%)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min-0.1,   Ymax2=col_NOL_max+0.1,    height2=3.2,   width2=5.65 , center2="TSS" )

MyAverageLines_3(vector2=  c( colMeans(reduceRow_Average_three3_week0_EEDheto[c(5,6), ]),   colMeans(reduceRow_Average_three3_week0_EEDko[c(5,6),]),  
                              colMeans(reduceRow_Average_three3_week4_EEDheto[c(5,6), ]),   colMeans(reduceRow_Average_three3_week4_EEDko[c(5,6),] ) ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_7_part3,     fileName2="3C-three3-4Samples-curve",  
                 title2="Up-regulated Genes (40%~60%)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min-0.1,   Ymax2=col_NOL_max+0.1,    height2=3.2,   width2=5.65 , center2="TSS" )

MyAverageLines_3(vector2=  c( colMeans(reduceRow_Average_three3_week0_EEDheto[c(7,8), ]),   colMeans(reduceRow_Average_three3_week0_EEDko[c(7,8),]),  
                              colMeans(reduceRow_Average_three3_week4_EEDheto[c(7,8), ]),   colMeans(reduceRow_Average_three3_week4_EEDko[c(7,8),] ) ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_7_part3,     fileName2="3D-three3-4Samples-curve",  
                 title2="Up-regulated Genes (60%~80%)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min-0.1,   Ymax2=col_NOL_max+0.1,    height2=3.2,   width2=5.65 , center2="TSS" )


MyAverageLines_3(vector2=  c( colMeans(reduceRow_Average_three3_week0_EEDheto[c(9,10), ]),   colMeans(reduceRow_Average_three3_week0_EEDko[c(9,10),]),  
                              colMeans(reduceRow_Average_three3_week4_EEDheto[c(9,10), ]),   colMeans(reduceRow_Average_three3_week4_EEDko[c(9,10),] ) ), 
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),   rep("week0_EEDko", numOfColumns1),  rep("week4_EEDheto", numOfColumns1),   rep("week4_EEDko", numOfColumns1)  ), 
                 sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                 colours2=c("week0_EEDheto"="red",  "week0_EEDko"="red4",  "week4_EEDheto"="blue",  "week4_EEDko"="blue4"),   
                 path2=subdir_7_part3,     fileName2="3E-three3-4Samples-curve",  
                 title2="Up-regulated Genes (80%~100%)",      xLab2="Relative distance (kb)",   yLab2="H2BGFP signal",   
                 Ymin2=col_NOL_min-0.1,   Ymax2=col_NOL_max+0.1,    height2=3.2,   width2=5.65 , center2="TSS" )



















###############################################################################
subdir_8_part3 <- paste(Part3_g,  "/8-rows5Classes-ViolinPlot", sep = "")
if( ! file.exists(subdir_8_part3) ) { dir.create(subdir_8_part3) }


row_Average_one1 <- numClasses( c(1:numOfRows1), binNum1=5)
row_Average_one1_len1 <- length(row_Average_one1[[1]])
row_Average_one1_len2 <- length(row_Average_one1[[2]])
row_Average_one1_len3 <- length(row_Average_one1[[3]])
row_Average_one1_len4 <- length(row_Average_one1[[4]])
row_Average_one1_len5 <- length(row_Average_one1[[5]])

MyBoxViolinPlot_1(vector2=c(row_Average_one1_week0_EEDheto[ row_Average_one1[[1]] ],  row_Average_one1_week0_EEDko[ row_Average_one1[[1]] ], 
                            row_Average_one1_week4_EEDheto[ row_Average_one1[[1]] ],  row_Average_one1_week4_EEDko[ row_Average_one1[[1]] ] ),   
                  sampleType2=c( rep("week0_EEDheto", row_Average_one1_len1),   rep("week0_EEDko", row_Average_one1_len1), 
                                 rep("week4_EEDheto", row_Average_one1_len1),   rep("week4_EEDko", row_Average_one1_len1) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_8_part3,  fileName2="1A-one1-BoxViolin",  
                  title2="Down-regulated Genes (0%~20%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )

MyBoxViolinPlot_1(vector2=c(row_Average_one1_week0_EEDheto[ row_Average_one1[[2]] ],  row_Average_one1_week0_EEDko[ row_Average_one1[[2]] ], 
                            row_Average_one1_week4_EEDheto[ row_Average_one1[[2]] ],  row_Average_one1_week4_EEDko[ row_Average_one1[[2]] ] ),   
                  sampleType2=c( rep("week0_EEDheto", row_Average_one1_len2),   rep("week0_EEDko", row_Average_one1_len2), 
                                 rep("week4_EEDheto", row_Average_one1_len2),   rep("week4_EEDko", row_Average_one1_len2) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_8_part3,  fileName2="1B-one1-BoxViolin",  
                  title2="Down-regulated Genes (20%~40%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )

MyBoxViolinPlot_1(vector2=c(row_Average_one1_week0_EEDheto[ row_Average_one1[[3]] ],  row_Average_one1_week0_EEDko[ row_Average_one1[[3]] ], 
                            row_Average_one1_week4_EEDheto[ row_Average_one1[[3]] ],  row_Average_one1_week4_EEDko[ row_Average_one1[[3]] ] ),   
                  sampleType2=c( rep("week0_EEDheto", row_Average_one1_len3),   rep("week0_EEDko", row_Average_one1_len3), 
                                 rep("week4_EEDheto", row_Average_one1_len3),   rep("week4_EEDko", row_Average_one1_len3) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_8_part3,  fileName2="1C-one1-BoxViolin",  
                  title2="Down-regulated Genes (40%~60%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )


MyBoxViolinPlot_1(vector2=c(row_Average_one1_week0_EEDheto[ row_Average_one1[[4]] ],  row_Average_one1_week0_EEDko[ row_Average_one1[[4]] ], 
                            row_Average_one1_week4_EEDheto[ row_Average_one1[[4]] ],  row_Average_one1_week4_EEDko[ row_Average_one1[[4]] ] ),   
                  sampleType2=c( rep("week0_EEDheto", row_Average_one1_len4),   rep("week0_EEDko", row_Average_one1_len4), 
                                 rep("week4_EEDheto", row_Average_one1_len4),   rep("week4_EEDko", row_Average_one1_len4) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_8_part3,  fileName2="1D-one1-BoxViolin",  
                  title2="Down-regulated Genes (60%~80%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )


MyBoxViolinPlot_1(vector2=c(row_Average_one1_week0_EEDheto[ row_Average_one1[[5]] ],  row_Average_one1_week0_EEDko[ row_Average_one1[[5]] ], 
                            row_Average_one1_week4_EEDheto[ row_Average_one1[[5]] ],  row_Average_one1_week4_EEDko[ row_Average_one1[[5]] ] ),   
                  sampleType2=c( rep("week0_EEDheto", row_Average_one1_len5),   rep("week0_EEDko", row_Average_one1_len5), 
                                 rep("week4_EEDheto", row_Average_one1_len5),   rep("week4_EEDko", row_Average_one1_len5) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_8_part3,  fileName2="1E-one1-BoxViolin",  
                  title2="Down-regulated Genes (80%~100%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )



MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[1]] ], vector2=row_Average_one1_week0_EEDko[ row_Average_one1[[1]] ], 
                   file1=paste(subdir_8_part3,  "/1A-row_one1_week0_EEDheto-VS-row_one1_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[1]] ], vector2=row_Average_one1_week4_EEDheto[ row_Average_one1[[1]] ], 
                   file1=paste(subdir_8_part3,  "/1A-row_one1_week0_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDko[ row_Average_one1[[1]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[1]] ], 
                   file1=paste(subdir_8_part3,  "/1A-row_one1_week0_EEDko-VS-row_one1_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week4_EEDheto[ row_Average_one1[[1]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[1]] ], 
                   file1=paste(subdir_8_part3,  "/1A-row_one1_week4_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[2]] ], vector2=row_Average_one1_week0_EEDko[ row_Average_one1[[2]] ], 
                   file1=paste(subdir_8_part3,  "/1B-row_one1_week0_EEDheto-VS-row_one1_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[2]] ], vector2=row_Average_one1_week4_EEDheto[ row_Average_one1[[2]] ], 
                   file1=paste(subdir_8_part3,  "/1B-row_one1_week0_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDko[ row_Average_one1[[2]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[2]] ], 
                   file1=paste(subdir_8_part3,  "/1B-row_one1_week0_EEDko-VS-row_one1_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week4_EEDheto[ row_Average_one1[[2]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[2]] ], 
                   file1=paste(subdir_8_part3,  "/1B-row_one1_week4_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[3]] ], vector2=row_Average_one1_week0_EEDko[ row_Average_one1[[3]] ], 
                   file1=paste(subdir_8_part3,  "/1C-row_one1_week0_EEDheto-VS-row_one1_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[3]] ], vector2=row_Average_one1_week4_EEDheto[ row_Average_one1[[3]] ], 
                   file1=paste(subdir_8_part3,  "/1C-row_one1_week0_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDko[ row_Average_one1[[3]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[3]] ], 
                   file1=paste(subdir_8_part3,  "/1C-row_one1_week0_EEDko-VS-row_one1_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week4_EEDheto[ row_Average_one1[[3]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[3]] ], 
                   file1=paste(subdir_8_part3,  "/1C-row_one1_week4_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[4]] ], vector2=row_Average_one1_week0_EEDko[ row_Average_one1[[4]] ], 
                   file1=paste(subdir_8_part3,  "/1D-row_one1_week0_EEDheto-VS-row_one1_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[4]] ], vector2=row_Average_one1_week4_EEDheto[ row_Average_one1[[4]] ], 
                   file1=paste(subdir_8_part3,  "/1D-row_one1_week0_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDko[ row_Average_one1[[4]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[4]] ], 
                   file1=paste(subdir_8_part3,  "/1D-row_one1_week0_EEDko-VS-row_one1_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week4_EEDheto[ row_Average_one1[[4]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[4]] ], 
                   file1=paste(subdir_8_part3,  "/1D-row_one1_week4_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "")  
)   



MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[5]] ], vector2=row_Average_one1_week0_EEDko[ row_Average_one1[[5]] ], 
                   file1=paste(subdir_8_part3,  "/1E-row_one1_week0_EEDheto-VS-row_one1_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[5]] ], vector2=row_Average_one1_week4_EEDheto[ row_Average_one1[[5]] ], 
                   file1=paste(subdir_8_part3,  "/1E-row_one1_week0_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDko[ row_Average_one1[[5]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[5]] ], 
                   file1=paste(subdir_8_part3,  "/1E-row_one1_week0_EEDko-VS-row_one1_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week4_EEDheto[ row_Average_one1[[5]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[5]] ], 
                   file1=paste(subdir_8_part3,  "/1E-row_one1_week4_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "")  
)   
















row_Average_two2 <- numClasses( c(1:numOfRows1), binNum1=5)
row_Average_two2_len1 <- length(row_Average_two2[[1]])
row_Average_two2_len2 <- length(row_Average_two2[[2]])
row_Average_two2_len3 <- length(row_Average_two2[[3]])
row_Average_two2_len4 <- length(row_Average_two2[[4]])
row_Average_two2_len5 <- length(row_Average_two2[[5]])

MyBoxViolinPlot_1(vector2=c(row_Average_two2_week0_EEDheto[ row_Average_two2[[1]] ],  row_Average_two2_week0_EEDko[ row_Average_two2[[1]] ], 
                            row_Average_two2_week4_EEDheto[ row_Average_two2[[1]] ],  row_Average_two2_week4_EEDko[ row_Average_two2[[1]] ] ),   
                  sampleType2=c( rep("week0_EEDheto", row_Average_two2_len1),   rep("week0_EEDko", row_Average_two2_len1), 
                                 rep("week4_EEDheto", row_Average_two2_len1),   rep("week4_EEDko", row_Average_two2_len1) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_8_part3,  fileName2="2A-two2-BoxViolin",  
                  title2="Down-regulated Genes (0%~20%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )

MyBoxViolinPlot_1(vector2=c(row_Average_two2_week0_EEDheto[ row_Average_two2[[2]] ],  row_Average_two2_week0_EEDko[ row_Average_two2[[2]] ], 
                            row_Average_two2_week4_EEDheto[ row_Average_two2[[2]] ],  row_Average_two2_week4_EEDko[ row_Average_two2[[2]] ] ),   
                  sampleType2=c( rep("week0_EEDheto", row_Average_two2_len2),   rep("week0_EEDko", row_Average_two2_len2), 
                                 rep("week4_EEDheto", row_Average_two2_len2),   rep("week4_EEDko", row_Average_two2_len2) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_8_part3,  fileName2="2B-two2-BoxViolin",  
                  title2="Down-regulated Genes (20%~40%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )

MyBoxViolinPlot_1(vector2=c(row_Average_two2_week0_EEDheto[ row_Average_two2[[3]] ],  row_Average_two2_week0_EEDko[ row_Average_two2[[3]] ], 
                            row_Average_two2_week4_EEDheto[ row_Average_two2[[3]] ],  row_Average_two2_week4_EEDko[ row_Average_two2[[3]] ] ),   
                  sampleType2=c( rep("week0_EEDheto", row_Average_two2_len3),   rep("week0_EEDko", row_Average_two2_len3), 
                                 rep("week4_EEDheto", row_Average_two2_len3),   rep("week4_EEDko", row_Average_two2_len3) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_8_part3,  fileName2="2C-two2-BoxViolin",  
                  title2="Down-regulated Genes (40%~60%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )


MyBoxViolinPlot_1(vector2=c(row_Average_two2_week0_EEDheto[ row_Average_two2[[4]] ],  row_Average_two2_week0_EEDko[ row_Average_two2[[4]] ], 
                            row_Average_two2_week4_EEDheto[ row_Average_two2[[4]] ],  row_Average_two2_week4_EEDko[ row_Average_two2[[4]] ] ),   
                  sampleType2=c( rep("week0_EEDheto", row_Average_two2_len4),   rep("week0_EEDko", row_Average_two2_len4), 
                                 rep("week4_EEDheto", row_Average_two2_len4),   rep("week4_EEDko", row_Average_two2_len4) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_8_part3,  fileName2="2D-two2-BoxViolin",  
                  title2="Down-regulated Genes (60%~80%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )


MyBoxViolinPlot_1(vector2=c(row_Average_two2_week0_EEDheto[ row_Average_two2[[5]] ],  row_Average_two2_week0_EEDko[ row_Average_two2[[5]] ], 
                            row_Average_two2_week4_EEDheto[ row_Average_two2[[5]] ],  row_Average_two2_week4_EEDko[ row_Average_two2[[5]] ] ),   
                  sampleType2=c( rep("week0_EEDheto", row_Average_two2_len5),   rep("week0_EEDko", row_Average_two2_len5), 
                                 rep("week4_EEDheto", row_Average_two2_len5),   rep("week4_EEDko", row_Average_two2_len5) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_8_part3,  fileName2="2E-two2-BoxViolin",  
                  title2="Down-regulated Genes (80%~100%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )



MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[1]] ], vector2=row_Average_two2_week0_EEDko[ row_Average_two2[[1]] ], 
                   file1=paste(subdir_8_part3,  "/2A-row_two2_week0_EEDheto-VS-row_two2_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[1]] ], vector2=row_Average_two2_week4_EEDheto[ row_Average_two2[[1]] ], 
                   file1=paste(subdir_8_part3,  "/2A-row_two2_week0_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDko[ row_Average_two2[[1]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[1]] ], 
                   file1=paste(subdir_8_part3,  "/2A-row_two2_week0_EEDko-VS-row_two2_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week4_EEDheto[ row_Average_two2[[1]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[1]] ], 
                   file1=paste(subdir_8_part3,  "/2A-row_two2_week4_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[2]] ], vector2=row_Average_two2_week0_EEDko[ row_Average_two2[[2]] ], 
                   file1=paste(subdir_8_part3,  "/2B-row_two2_week0_EEDheto-VS-row_two2_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[2]] ], vector2=row_Average_two2_week4_EEDheto[ row_Average_two2[[2]] ], 
                   file1=paste(subdir_8_part3,  "/2B-row_two2_week0_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDko[ row_Average_two2[[2]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[2]] ], 
                   file1=paste(subdir_8_part3,  "/2B-row_two2_week0_EEDko-VS-row_two2_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week4_EEDheto[ row_Average_two2[[2]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[2]] ], 
                   file1=paste(subdir_8_part3,  "/2B-row_two2_week4_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[3]] ], vector2=row_Average_two2_week0_EEDko[ row_Average_two2[[3]] ], 
                   file1=paste(subdir_8_part3,  "/2C-row_two2_week0_EEDheto-VS-row_two2_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[3]] ], vector2=row_Average_two2_week4_EEDheto[ row_Average_two2[[3]] ], 
                   file1=paste(subdir_8_part3,  "/2C-row_two2_week0_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDko[ row_Average_two2[[3]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[3]] ], 
                   file1=paste(subdir_8_part3,  "/2C-row_two2_week0_EEDko-VS-row_two2_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week4_EEDheto[ row_Average_two2[[3]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[3]] ], 
                   file1=paste(subdir_8_part3,  "/2C-row_two2_week4_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[4]] ], vector2=row_Average_two2_week0_EEDko[ row_Average_two2[[4]] ], 
                   file1=paste(subdir_8_part3,  "/2D-row_two2_week0_EEDheto-VS-row_two2_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[4]] ], vector2=row_Average_two2_week4_EEDheto[ row_Average_two2[[4]] ], 
                   file1=paste(subdir_8_part3,  "/2D-row_two2_week0_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDko[ row_Average_two2[[4]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[4]] ], 
                   file1=paste(subdir_8_part3,  "/2D-row_two2_week0_EEDko-VS-row_two2_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week4_EEDheto[ row_Average_two2[[4]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[4]] ], 
                   file1=paste(subdir_8_part3,  "/2D-row_two2_week4_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[5]] ], vector2=row_Average_two2_week0_EEDko[ row_Average_two2[[5]] ], 
                   file1=paste(subdir_8_part3,  "/2E-row_two2_week0_EEDheto-VS-row_two2_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[5]] ], vector2=row_Average_two2_week4_EEDheto[ row_Average_two2[[5]] ], 
                   file1=paste(subdir_8_part3,  "/2E-row_two2_week0_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDko[ row_Average_two2[[5]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[5]] ], 
                   file1=paste(subdir_8_part3,  "/2E-row_two2_week0_EEDko-VS-row_two2_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week4_EEDheto[ row_Average_two2[[5]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[5]] ], 
                   file1=paste(subdir_8_part3,  "/2E-row_two2_week4_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "")  
)   















row_Average_three3 <- numClasses( c(1:numOfRows1), binNum1=5)
row_Average_three3_len1 <- length(row_Average_three3[[1]])
row_Average_three3_len2 <- length(row_Average_three3[[2]])
row_Average_three3_len3 <- length(row_Average_three3[[3]])
row_Average_three3_len4 <- length(row_Average_three3[[4]])
row_Average_three3_len5 <- length(row_Average_three3[[5]])

MyBoxViolinPlot_1(vector2=c(row_Average_three3_week0_EEDheto[ row_Average_three3[[1]] ],  row_Average_three3_week0_EEDko[ row_Average_three3[[1]] ], 
                            row_Average_three3_week4_EEDheto[ row_Average_three3[[1]] ],  row_Average_three3_week4_EEDko[ row_Average_three3[[1]] ] ),   
                  sampleType2=c( rep("week0_EEDheto", row_Average_three3_len1),   rep("week0_EEDko", row_Average_three3_len1), 
                                 rep("week4_EEDheto", row_Average_three3_len1),   rep("week4_EEDko", row_Average_three3_len1) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_8_part3,  fileName2="3A-three3-BoxViolin",  
                  title2="Down-regulated Genes (0%~20%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )

MyBoxViolinPlot_1(vector2=c(row_Average_three3_week0_EEDheto[ row_Average_three3[[2]] ],  row_Average_three3_week0_EEDko[ row_Average_three3[[2]] ], 
                            row_Average_three3_week4_EEDheto[ row_Average_three3[[2]] ],  row_Average_three3_week4_EEDko[ row_Average_three3[[2]] ] ),   
                  sampleType2=c( rep("week0_EEDheto", row_Average_three3_len2),   rep("week0_EEDko", row_Average_three3_len2), 
                                 rep("week4_EEDheto", row_Average_three3_len2),   rep("week4_EEDko", row_Average_three3_len2) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_8_part3,  fileName2="3B-three3-BoxViolin",  
                  title2="Down-regulated Genes (20%~40%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )

MyBoxViolinPlot_1(vector2=c(row_Average_three3_week0_EEDheto[ row_Average_three3[[3]] ],  row_Average_three3_week0_EEDko[ row_Average_three3[[3]] ], 
                            row_Average_three3_week4_EEDheto[ row_Average_three3[[3]] ],  row_Average_three3_week4_EEDko[ row_Average_three3[[3]] ] ),   
                  sampleType2=c( rep("week0_EEDheto", row_Average_three3_len3),   rep("week0_EEDko", row_Average_three3_len3), 
                                 rep("week4_EEDheto", row_Average_three3_len3),   rep("week4_EEDko", row_Average_three3_len3) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_8_part3,  fileName2="3C-three3-BoxViolin",  
                  title2="Down-regulated Genes (40%~60%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )


MyBoxViolinPlot_1(vector2=c(row_Average_three3_week0_EEDheto[ row_Average_three3[[4]] ],  row_Average_three3_week0_EEDko[ row_Average_three3[[4]] ], 
                            row_Average_three3_week4_EEDheto[ row_Average_three3[[4]] ],  row_Average_three3_week4_EEDko[ row_Average_three3[[4]] ] ),   
                  sampleType2=c( rep("week0_EEDheto", row_Average_three3_len4),   rep("week0_EEDko", row_Average_three3_len4), 
                                 rep("week4_EEDheto", row_Average_three3_len4),   rep("week4_EEDko", row_Average_three3_len4) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_8_part3,  fileName2="3D-three3-BoxViolin",  
                  title2="Down-regulated Genes (60%~80%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )


MyBoxViolinPlot_1(vector2=c(row_Average_three3_week0_EEDheto[ row_Average_three3[[5]] ],  row_Average_three3_week0_EEDko[ row_Average_three3[[5]] ], 
                            row_Average_three3_week4_EEDheto[ row_Average_three3[[5]] ],  row_Average_three3_week4_EEDko[ row_Average_three3[[5]] ] ),   
                  sampleType2=c( rep("week0_EEDheto", row_Average_three3_len5),   rep("week0_EEDko", row_Average_three3_len5), 
                                 rep("week4_EEDheto", row_Average_three3_len5),   rep("week4_EEDko", row_Average_three3_len5) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c("red",  "red4",   "blue",    "blue4" ),   
                  path2=subdir_8_part3,  fileName2="3E-three3-BoxViolin",  
                  title2="Down-regulated Genes (80%~100%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )



MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[1]] ], vector2=row_Average_three3_week0_EEDko[ row_Average_three3[[1]] ], 
                   file1=paste(subdir_8_part3,  "/3A-row_three3_week0_EEDheto-VS-row_three3_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[1]] ], vector2=row_Average_three3_week4_EEDheto[ row_Average_three3[[1]] ], 
                   file1=paste(subdir_8_part3,  "/3A-row_three3_week0_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDko[ row_Average_three3[[1]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[1]] ], 
                   file1=paste(subdir_8_part3,  "/3A-row_three3_week0_EEDko-VS-row_three3_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week4_EEDheto[ row_Average_three3[[1]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[1]] ], 
                   file1=paste(subdir_8_part3,  "/3A-row_three3_week4_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[2]] ], vector2=row_Average_three3_week0_EEDko[ row_Average_three3[[2]] ], 
                   file1=paste(subdir_8_part3,  "/3B-row_three3_week0_EEDheto-VS-row_three3_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[2]] ], vector2=row_Average_three3_week4_EEDheto[ row_Average_three3[[2]] ], 
                   file1=paste(subdir_8_part3,  "/3B-row_three3_week0_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDko[ row_Average_three3[[2]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[2]] ], 
                   file1=paste(subdir_8_part3,  "/3B-row_three3_week0_EEDko-VS-row_three3_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week4_EEDheto[ row_Average_three3[[2]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[2]] ], 
                   file1=paste(subdir_8_part3,  "/3B-row_three3_week4_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[3]] ], vector2=row_Average_three3_week0_EEDko[ row_Average_three3[[3]] ], 
                   file1=paste(subdir_8_part3,  "/3C-row_three3_week0_EEDheto-VS-row_three3_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[3]] ], vector2=row_Average_three3_week4_EEDheto[ row_Average_three3[[3]] ], 
                   file1=paste(subdir_8_part3,  "/3C-row_three3_week0_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDko[ row_Average_three3[[3]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[3]] ], 
                   file1=paste(subdir_8_part3,  "/3C-row_three3_week0_EEDko-VS-row_three3_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week4_EEDheto[ row_Average_three3[[3]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[3]] ], 
                   file1=paste(subdir_8_part3,  "/3C-row_three3_week4_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[4]] ], vector2=row_Average_three3_week0_EEDko[ row_Average_three3[[4]] ], 
                   file1=paste(subdir_8_part3,  "/3D-row_three3_week0_EEDheto-VS-row_three3_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[4]] ], vector2=row_Average_three3_week4_EEDheto[ row_Average_three3[[4]] ], 
                   file1=paste(subdir_8_part3,  "/3D-row_three3_week0_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDko[ row_Average_three3[[4]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[4]] ], 
                   file1=paste(subdir_8_part3,  "/3D-row_three3_week0_EEDko-VS-row_three3_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week4_EEDheto[ row_Average_three3[[4]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[4]] ], 
                   file1=paste(subdir_8_part3,  "/3D-row_three3_week4_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "")  
)   






MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[5]] ], vector2=row_Average_three3_week0_EEDko[ row_Average_three3[[5]] ], 
                   file1=paste(subdir_8_part3,  "/3E-row_three3_week0_EEDheto-VS-row_three3_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[5]] ], vector2=row_Average_three3_week4_EEDheto[ row_Average_three3[[5]] ], 
                   file1=paste(subdir_8_part3,  "/3E-row_three3_week0_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDko[ row_Average_three3[[5]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[5]] ], 
                   file1=paste(subdir_8_part3,  "/3E-row_three3_week0_EEDko-VS-row_three3_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week4_EEDheto[ row_Average_three3[[5]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[5]] ], 
                   file1=paste(subdir_8_part3,  "/3E-row_three3_week4_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "")  
)   






















###############################################################################
subdir_9_part3 <- paste(Part3_g,  "/9-rows5Classes-path", sep = "")
if( ! file.exists(subdir_9_part3) ) { dir.create(subdir_9_part3) }


row_Average_one1 <- numClasses( c(1:numOfRows1), binNum1=5)
row_Average_one1_len1 <- length(row_Average_one1[[1]])
row_Average_one1_len2 <- length(row_Average_one1[[2]])
row_Average_one1_len3 <- length(row_Average_one1[[3]])
row_Average_one1_len4 <- length(row_Average_one1[[4]])
row_Average_one1_len5 <- length(row_Average_one1[[5]])

MyBoxViolinPlot_3_s4(vector2=c(row_Average_one1_week0_EEDheto[ row_Average_one1[[1]] ],  row_Average_one1_week0_EEDko[ row_Average_one1[[1]] ], 
                               row_Average_one1_week4_EEDheto[ row_Average_one1[[1]] ],  row_Average_one1_week4_EEDko[ row_Average_one1[[1]] ] ),   
                     sampleType2=c( rep("week0_EEDheto", row_Average_one1_len1),   rep("week0_EEDko", row_Average_one1_len1), 
                                    rep("week4_EEDheto", row_Average_one1_len1),   rep("week4_EEDko", row_Average_one1_len1) ), 
                     sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                     colours2=c("red",  "red4",   "blue",    "blue4" ),   
                     path2=subdir_9_part3,  fileName2="1A-one1-BoxViolin",  
                     title2="Down-regulated Genes (0%~20%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                     height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )

MyBoxViolinPlot_3_s4(vector2=c(row_Average_one1_week0_EEDheto[ row_Average_one1[[2]] ],  row_Average_one1_week0_EEDko[ row_Average_one1[[2]] ], 
                               row_Average_one1_week4_EEDheto[ row_Average_one1[[2]] ],  row_Average_one1_week4_EEDko[ row_Average_one1[[2]] ] ),   
                     sampleType2=c( rep("week0_EEDheto", row_Average_one1_len2),   rep("week0_EEDko", row_Average_one1_len2), 
                                    rep("week4_EEDheto", row_Average_one1_len2),   rep("week4_EEDko", row_Average_one1_len2) ), 
                     sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                     colours2=c("red",  "red4",   "blue",    "blue4" ),   
                     path2=subdir_9_part3,  fileName2="1B-one1-BoxViolin",  
                     title2="Down-regulated Genes (20%~40%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                     height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )

MyBoxViolinPlot_3_s4(vector2=c(row_Average_one1_week0_EEDheto[ row_Average_one1[[3]] ],  row_Average_one1_week0_EEDko[ row_Average_one1[[3]] ], 
                               row_Average_one1_week4_EEDheto[ row_Average_one1[[3]] ],  row_Average_one1_week4_EEDko[ row_Average_one1[[3]] ] ),   
                     sampleType2=c( rep("week0_EEDheto", row_Average_one1_len3),   rep("week0_EEDko", row_Average_one1_len3), 
                                    rep("week4_EEDheto", row_Average_one1_len3),   rep("week4_EEDko", row_Average_one1_len3) ), 
                     sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                     colours2=c("red",  "red4",   "blue",    "blue4" ),   
                     path2=subdir_9_part3,  fileName2="1C-one1-BoxViolin",  
                     title2="Down-regulated Genes (40%~60%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                     height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )


MyBoxViolinPlot_3_s4(vector2=c(row_Average_one1_week0_EEDheto[ row_Average_one1[[4]] ],  row_Average_one1_week0_EEDko[ row_Average_one1[[4]] ], 
                               row_Average_one1_week4_EEDheto[ row_Average_one1[[4]] ],  row_Average_one1_week4_EEDko[ row_Average_one1[[4]] ] ),   
                     sampleType2=c( rep("week0_EEDheto", row_Average_one1_len4),   rep("week0_EEDko", row_Average_one1_len4), 
                                    rep("week4_EEDheto", row_Average_one1_len4),   rep("week4_EEDko", row_Average_one1_len4) ), 
                     sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                     colours2=c("red",  "red4",   "blue",    "blue4" ),   
                     path2=subdir_9_part3,  fileName2="1D-one1-BoxViolin",  
                     title2="Down-regulated Genes (60%~80%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                     height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )


MyBoxViolinPlot_3_s4(vector2=c(row_Average_one1_week0_EEDheto[ row_Average_one1[[5]] ],  row_Average_one1_week0_EEDko[ row_Average_one1[[5]] ], 
                               row_Average_one1_week4_EEDheto[ row_Average_one1[[5]] ],  row_Average_one1_week4_EEDko[ row_Average_one1[[5]] ] ),   
                     sampleType2=c( rep("week0_EEDheto", row_Average_one1_len5),   rep("week0_EEDko", row_Average_one1_len5), 
                                    rep("week4_EEDheto", row_Average_one1_len5),   rep("week4_EEDko", row_Average_one1_len5) ), 
                     sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                     colours2=c("red",  "red4",   "blue",    "blue4" ),   
                     path2=subdir_9_part3,  fileName2="1E-one1-BoxViolin",  
                     title2="Down-regulated Genes (80%~100%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                     height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )



MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[1]] ], vector2=row_Average_one1_week0_EEDko[ row_Average_one1[[1]] ], 
                   file1=paste(subdir_9_part3,  "/1A-row_one1_week0_EEDheto-VS-row_one1_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[1]] ], vector2=row_Average_one1_week4_EEDheto[ row_Average_one1[[1]] ], 
                   file1=paste(subdir_9_part3,  "/1A-row_one1_week0_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDko[ row_Average_one1[[1]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[1]] ], 
                   file1=paste(subdir_9_part3,  "/1A-row_one1_week0_EEDko-VS-row_one1_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week4_EEDheto[ row_Average_one1[[1]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[1]] ], 
                   file1=paste(subdir_9_part3,  "/1A-row_one1_week4_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[2]] ], vector2=row_Average_one1_week0_EEDko[ row_Average_one1[[2]] ], 
                   file1=paste(subdir_9_part3,  "/1B-row_one1_week0_EEDheto-VS-row_one1_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[2]] ], vector2=row_Average_one1_week4_EEDheto[ row_Average_one1[[2]] ], 
                   file1=paste(subdir_9_part3,  "/1B-row_one1_week0_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDko[ row_Average_one1[[2]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[2]] ], 
                   file1=paste(subdir_9_part3,  "/1B-row_one1_week0_EEDko-VS-row_one1_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week4_EEDheto[ row_Average_one1[[2]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[2]] ], 
                   file1=paste(subdir_9_part3,  "/1B-row_one1_week4_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[3]] ], vector2=row_Average_one1_week0_EEDko[ row_Average_one1[[3]] ], 
                   file1=paste(subdir_9_part3,  "/1C-row_one1_week0_EEDheto-VS-row_one1_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[3]] ], vector2=row_Average_one1_week4_EEDheto[ row_Average_one1[[3]] ], 
                   file1=paste(subdir_9_part3,  "/1C-row_one1_week0_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDko[ row_Average_one1[[3]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[3]] ], 
                   file1=paste(subdir_9_part3,  "/1C-row_one1_week0_EEDko-VS-row_one1_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week4_EEDheto[ row_Average_one1[[3]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[3]] ], 
                   file1=paste(subdir_9_part3,  "/1C-row_one1_week4_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[4]] ], vector2=row_Average_one1_week0_EEDko[ row_Average_one1[[4]] ], 
                   file1=paste(subdir_9_part3,  "/1D-row_one1_week0_EEDheto-VS-row_one1_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[4]] ], vector2=row_Average_one1_week4_EEDheto[ row_Average_one1[[4]] ], 
                   file1=paste(subdir_9_part3,  "/1D-row_one1_week0_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDko[ row_Average_one1[[4]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[4]] ], 
                   file1=paste(subdir_9_part3,  "/1D-row_one1_week0_EEDko-VS-row_one1_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week4_EEDheto[ row_Average_one1[[4]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[4]] ], 
                   file1=paste(subdir_9_part3,  "/1D-row_one1_week4_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "")  
)   



MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[5]] ], vector2=row_Average_one1_week0_EEDko[ row_Average_one1[[5]] ], 
                   file1=paste(subdir_9_part3,  "/1E-row_one1_week0_EEDheto-VS-row_one1_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDheto[ row_Average_one1[[5]] ], vector2=row_Average_one1_week4_EEDheto[ row_Average_one1[[5]] ], 
                   file1=paste(subdir_9_part3,  "/1E-row_one1_week0_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week0_EEDko[ row_Average_one1[[5]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[5]] ], 
                   file1=paste(subdir_9_part3,  "/1E-row_one1_week0_EEDko-VS-row_one1_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_one1_week4_EEDheto[ row_Average_one1[[5]] ], vector2=row_Average_one1_week4_EEDko[ row_Average_one1[[5]] ], 
                   file1=paste(subdir_9_part3,  "/1E-row_one1_week4_EEDheto-VS-row_one1_week4_EEDheto.txt", sep = "")  
)   
















row_Average_two2 <- numClasses( c(1:numOfRows1), binNum1=5)
row_Average_two2_len1 <- length(row_Average_two2[[1]])
row_Average_two2_len2 <- length(row_Average_two2[[2]])
row_Average_two2_len3 <- length(row_Average_two2[[3]])
row_Average_two2_len4 <- length(row_Average_two2[[4]])
row_Average_two2_len5 <- length(row_Average_two2[[5]])

MyBoxViolinPlot_3_s4(vector2=c(row_Average_two2_week0_EEDheto[ row_Average_two2[[1]] ],  row_Average_two2_week0_EEDko[ row_Average_two2[[1]] ], 
                               row_Average_two2_week4_EEDheto[ row_Average_two2[[1]] ],  row_Average_two2_week4_EEDko[ row_Average_two2[[1]] ] ),   
                     sampleType2=c( rep("week0_EEDheto", row_Average_two2_len1),   rep("week0_EEDko", row_Average_two2_len1), 
                                    rep("week4_EEDheto", row_Average_two2_len1),   rep("week4_EEDko", row_Average_two2_len1) ), 
                     sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                     colours2=c("red",  "red4",   "blue",    "blue4" ),   
                     path2=subdir_9_part3,  fileName2="2A-two2-BoxViolin",  
                     title2="Down-regulated Genes (0%~20%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                     height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )

MyBoxViolinPlot_3_s4(vector2=c(row_Average_two2_week0_EEDheto[ row_Average_two2[[2]] ],  row_Average_two2_week0_EEDko[ row_Average_two2[[2]] ], 
                               row_Average_two2_week4_EEDheto[ row_Average_two2[[2]] ],  row_Average_two2_week4_EEDko[ row_Average_two2[[2]] ] ),   
                     sampleType2=c( rep("week0_EEDheto", row_Average_two2_len2),   rep("week0_EEDko", row_Average_two2_len2), 
                                    rep("week4_EEDheto", row_Average_two2_len2),   rep("week4_EEDko", row_Average_two2_len2) ), 
                     sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                     colours2=c("red",  "red4",   "blue",    "blue4" ),   
                     path2=subdir_9_part3,  fileName2="2B-two2-BoxViolin",  
                     title2="Down-regulated Genes (20%~40%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                     height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )

MyBoxViolinPlot_3_s4(vector2=c(row_Average_two2_week0_EEDheto[ row_Average_two2[[3]] ],  row_Average_two2_week0_EEDko[ row_Average_two2[[3]] ], 
                               row_Average_two2_week4_EEDheto[ row_Average_two2[[3]] ],  row_Average_two2_week4_EEDko[ row_Average_two2[[3]] ] ),   
                     sampleType2=c( rep("week0_EEDheto", row_Average_two2_len3),   rep("week0_EEDko", row_Average_two2_len3), 
                                    rep("week4_EEDheto", row_Average_two2_len3),   rep("week4_EEDko", row_Average_two2_len3) ), 
                     sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                     colours2=c("red",  "red4",   "blue",    "blue4" ),   
                     path2=subdir_9_part3,  fileName2="2C-two2-BoxViolin",  
                     title2="Down-regulated Genes (40%~60%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                     height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )


MyBoxViolinPlot_3_s4(vector2=c(row_Average_two2_week0_EEDheto[ row_Average_two2[[4]] ],  row_Average_two2_week0_EEDko[ row_Average_two2[[4]] ], 
                               row_Average_two2_week4_EEDheto[ row_Average_two2[[4]] ],  row_Average_two2_week4_EEDko[ row_Average_two2[[4]] ] ),   
                     sampleType2=c( rep("week0_EEDheto", row_Average_two2_len4),   rep("week0_EEDko", row_Average_two2_len4), 
                                    rep("week4_EEDheto", row_Average_two2_len4),   rep("week4_EEDko", row_Average_two2_len4) ), 
                     sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                     colours2=c("red",  "red4",   "blue",    "blue4" ),   
                     path2=subdir_9_part3,  fileName2="2D-two2-BoxViolin",  
                     title2="Down-regulated Genes (60%~80%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                     height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )


MyBoxViolinPlot_3_s4(vector2=c(row_Average_two2_week0_EEDheto[ row_Average_two2[[5]] ],  row_Average_two2_week0_EEDko[ row_Average_two2[[5]] ], 
                               row_Average_two2_week4_EEDheto[ row_Average_two2[[5]] ],  row_Average_two2_week4_EEDko[ row_Average_two2[[5]] ] ),   
                     sampleType2=c( rep("week0_EEDheto", row_Average_two2_len5),   rep("week0_EEDko", row_Average_two2_len5), 
                                    rep("week4_EEDheto", row_Average_two2_len5),   rep("week4_EEDko", row_Average_two2_len5) ), 
                     sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                     colours2=c("red",  "red4",   "blue",    "blue4" ),   
                     path2=subdir_9_part3,  fileName2="2E-two2-BoxViolin",  
                     title2="Down-regulated Genes (80%~100%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                     height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )



MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[1]] ], vector2=row_Average_two2_week0_EEDko[ row_Average_two2[[1]] ], 
                   file1=paste(subdir_9_part3,  "/2A-row_two2_week0_EEDheto-VS-row_two2_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[1]] ], vector2=row_Average_two2_week4_EEDheto[ row_Average_two2[[1]] ], 
                   file1=paste(subdir_9_part3,  "/2A-row_two2_week0_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDko[ row_Average_two2[[1]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[1]] ], 
                   file1=paste(subdir_9_part3,  "/2A-row_two2_week0_EEDko-VS-row_two2_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week4_EEDheto[ row_Average_two2[[1]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[1]] ], 
                   file1=paste(subdir_9_part3,  "/2A-row_two2_week4_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[2]] ], vector2=row_Average_two2_week0_EEDko[ row_Average_two2[[2]] ], 
                   file1=paste(subdir_9_part3,  "/2B-row_two2_week0_EEDheto-VS-row_two2_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[2]] ], vector2=row_Average_two2_week4_EEDheto[ row_Average_two2[[2]] ], 
                   file1=paste(subdir_9_part3,  "/2B-row_two2_week0_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDko[ row_Average_two2[[2]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[2]] ], 
                   file1=paste(subdir_9_part3,  "/2B-row_two2_week0_EEDko-VS-row_two2_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week4_EEDheto[ row_Average_two2[[2]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[2]] ], 
                   file1=paste(subdir_9_part3,  "/2B-row_two2_week4_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[3]] ], vector2=row_Average_two2_week0_EEDko[ row_Average_two2[[3]] ], 
                   file1=paste(subdir_9_part3,  "/2C-row_two2_week0_EEDheto-VS-row_two2_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[3]] ], vector2=row_Average_two2_week4_EEDheto[ row_Average_two2[[3]] ], 
                   file1=paste(subdir_9_part3,  "/2C-row_two2_week0_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDko[ row_Average_two2[[3]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[3]] ], 
                   file1=paste(subdir_9_part3,  "/2C-row_two2_week0_EEDko-VS-row_two2_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week4_EEDheto[ row_Average_two2[[3]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[3]] ], 
                   file1=paste(subdir_9_part3,  "/2C-row_two2_week4_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[4]] ], vector2=row_Average_two2_week0_EEDko[ row_Average_two2[[4]] ], 
                   file1=paste(subdir_9_part3,  "/2D-row_two2_week0_EEDheto-VS-row_two2_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[4]] ], vector2=row_Average_two2_week4_EEDheto[ row_Average_two2[[4]] ], 
                   file1=paste(subdir_9_part3,  "/2D-row_two2_week0_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDko[ row_Average_two2[[4]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[4]] ], 
                   file1=paste(subdir_9_part3,  "/2D-row_two2_week0_EEDko-VS-row_two2_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week4_EEDheto[ row_Average_two2[[4]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[4]] ], 
                   file1=paste(subdir_9_part3,  "/2D-row_two2_week4_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[5]] ], vector2=row_Average_two2_week0_EEDko[ row_Average_two2[[5]] ], 
                   file1=paste(subdir_9_part3,  "/2E-row_two2_week0_EEDheto-VS-row_two2_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDheto[ row_Average_two2[[5]] ], vector2=row_Average_two2_week4_EEDheto[ row_Average_two2[[5]] ], 
                   file1=paste(subdir_9_part3,  "/2E-row_two2_week0_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week0_EEDko[ row_Average_two2[[5]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[5]] ], 
                   file1=paste(subdir_9_part3,  "/2E-row_two2_week0_EEDko-VS-row_two2_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_two2_week4_EEDheto[ row_Average_two2[[5]] ], vector2=row_Average_two2_week4_EEDko[ row_Average_two2[[5]] ], 
                   file1=paste(subdir_9_part3,  "/2E-row_two2_week4_EEDheto-VS-row_two2_week4_EEDheto.txt", sep = "")  
)   















row_Average_three3 <- numClasses( c(1:numOfRows1), binNum1=5)
row_Average_three3_len1 <- length(row_Average_three3[[1]])
row_Average_three3_len2 <- length(row_Average_three3[[2]])
row_Average_three3_len3 <- length(row_Average_three3[[3]])
row_Average_three3_len4 <- length(row_Average_three3[[4]])
row_Average_three3_len5 <- length(row_Average_three3[[5]])

MyBoxViolinPlot_3_s4(vector2=c(row_Average_three3_week0_EEDheto[ row_Average_three3[[1]] ],  row_Average_three3_week0_EEDko[ row_Average_three3[[1]] ], 
                               row_Average_three3_week4_EEDheto[ row_Average_three3[[1]] ],  row_Average_three3_week4_EEDko[ row_Average_three3[[1]] ] ),   
                     sampleType2=c( rep("week0_EEDheto", row_Average_three3_len1),   rep("week0_EEDko", row_Average_three3_len1), 
                                    rep("week4_EEDheto", row_Average_three3_len1),   rep("week4_EEDko", row_Average_three3_len1) ), 
                     sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                     colours2=c("red",  "red4",   "blue",    "blue4" ),   
                     path2=subdir_9_part3,  fileName2="3A-three3-BoxViolin",  
                     title2="Down-regulated Genes (0%~20%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                     height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )

MyBoxViolinPlot_3_s4(vector2=c(row_Average_three3_week0_EEDheto[ row_Average_three3[[2]] ],  row_Average_three3_week0_EEDko[ row_Average_three3[[2]] ], 
                               row_Average_three3_week4_EEDheto[ row_Average_three3[[2]] ],  row_Average_three3_week4_EEDko[ row_Average_three3[[2]] ] ),   
                     sampleType2=c( rep("week0_EEDheto", row_Average_three3_len2),   rep("week0_EEDko", row_Average_three3_len2), 
                                    rep("week4_EEDheto", row_Average_three3_len2),   rep("week4_EEDko", row_Average_three3_len2) ), 
                     sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                     colours2=c("red",  "red4",   "blue",    "blue4" ),   
                     path2=subdir_9_part3,  fileName2="3B-three3-BoxViolin",  
                     title2="Down-regulated Genes (20%~40%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                     height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )

MyBoxViolinPlot_3_s4(vector2=c(row_Average_three3_week0_EEDheto[ row_Average_three3[[3]] ],  row_Average_three3_week0_EEDko[ row_Average_three3[[3]] ], 
                               row_Average_three3_week4_EEDheto[ row_Average_three3[[3]] ],  row_Average_three3_week4_EEDko[ row_Average_three3[[3]] ] ),   
                     sampleType2=c( rep("week0_EEDheto", row_Average_three3_len3),   rep("week0_EEDko", row_Average_three3_len3), 
                                    rep("week4_EEDheto", row_Average_three3_len3),   rep("week4_EEDko", row_Average_three3_len3) ), 
                     sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                     colours2=c("red",  "red4",   "blue",    "blue4" ),   
                     path2=subdir_9_part3,  fileName2="3C-three3-BoxViolin",  
                     title2="Down-regulated Genes (40%~60%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                     height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )


MyBoxViolinPlot_3_s4(vector2=c(row_Average_three3_week0_EEDheto[ row_Average_three3[[4]] ],  row_Average_three3_week0_EEDko[ row_Average_three3[[4]] ], 
                               row_Average_three3_week4_EEDheto[ row_Average_three3[[4]] ],  row_Average_three3_week4_EEDko[ row_Average_three3[[4]] ] ),   
                     sampleType2=c( rep("week0_EEDheto", row_Average_three3_len4),   rep("week0_EEDko", row_Average_three3_len4), 
                                    rep("week4_EEDheto", row_Average_three3_len4),   rep("week4_EEDko", row_Average_three3_len4) ), 
                     sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                     colours2=c("red",  "red4",   "blue",    "blue4" ),   
                     path2=subdir_9_part3,  fileName2="3D-three3-BoxViolin",  
                     title2="Down-regulated Genes (60%~80%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                     height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )


MyBoxViolinPlot_3_s4(vector2=c(row_Average_three3_week0_EEDheto[ row_Average_three3[[5]] ],  row_Average_three3_week0_EEDko[ row_Average_three3[[5]] ], 
                               row_Average_three3_week4_EEDheto[ row_Average_three3[[5]] ],  row_Average_three3_week4_EEDko[ row_Average_three3[[5]] ] ),   
                     sampleType2=c( rep("week0_EEDheto", row_Average_three3_len5),   rep("week0_EEDko", row_Average_three3_len5), 
                                    rep("week4_EEDheto", row_Average_three3_len5),   rep("week4_EEDko", row_Average_three3_len5) ), 
                     sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                     colours2=c("red",  "red4",   "blue",    "blue4" ),   
                     path2=subdir_9_part3,  fileName2="3E-three3-BoxViolin",  
                     title2="Down-regulated Genes (80%~100%)",  xLab2="Samples",  yLab2="H2BGFP signal",   
                     height2=3.88,   width2=2.8, Ymin2=0, Ymax2=1.25 )



MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[1]] ], vector2=row_Average_three3_week0_EEDko[ row_Average_three3[[1]] ], 
                   file1=paste(subdir_9_part3,  "/3A-row_three3_week0_EEDheto-VS-row_three3_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[1]] ], vector2=row_Average_three3_week4_EEDheto[ row_Average_three3[[1]] ], 
                   file1=paste(subdir_9_part3,  "/3A-row_three3_week0_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDko[ row_Average_three3[[1]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[1]] ], 
                   file1=paste(subdir_9_part3,  "/3A-row_three3_week0_EEDko-VS-row_three3_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week4_EEDheto[ row_Average_three3[[1]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[1]] ], 
                   file1=paste(subdir_9_part3,  "/3A-row_three3_week4_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[2]] ], vector2=row_Average_three3_week0_EEDko[ row_Average_three3[[2]] ], 
                   file1=paste(subdir_9_part3,  "/3B-row_three3_week0_EEDheto-VS-row_three3_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[2]] ], vector2=row_Average_three3_week4_EEDheto[ row_Average_three3[[2]] ], 
                   file1=paste(subdir_9_part3,  "/3B-row_three3_week0_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDko[ row_Average_three3[[2]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[2]] ], 
                   file1=paste(subdir_9_part3,  "/3B-row_three3_week0_EEDko-VS-row_three3_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week4_EEDheto[ row_Average_three3[[2]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[2]] ], 
                   file1=paste(subdir_9_part3,  "/3B-row_three3_week4_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[3]] ], vector2=row_Average_three3_week0_EEDko[ row_Average_three3[[3]] ], 
                   file1=paste(subdir_9_part3,  "/3C-row_three3_week0_EEDheto-VS-row_three3_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[3]] ], vector2=row_Average_three3_week4_EEDheto[ row_Average_three3[[3]] ], 
                   file1=paste(subdir_9_part3,  "/3C-row_three3_week0_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDko[ row_Average_three3[[3]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[3]] ], 
                   file1=paste(subdir_9_part3,  "/3C-row_three3_week0_EEDko-VS-row_three3_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week4_EEDheto[ row_Average_three3[[3]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[3]] ], 
                   file1=paste(subdir_9_part3,  "/3C-row_three3_week4_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "")  
)   




MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[4]] ], vector2=row_Average_three3_week0_EEDko[ row_Average_three3[[4]] ], 
                   file1=paste(subdir_9_part3,  "/3D-row_three3_week0_EEDheto-VS-row_three3_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[4]] ], vector2=row_Average_three3_week4_EEDheto[ row_Average_three3[[4]] ], 
                   file1=paste(subdir_9_part3,  "/3D-row_three3_week0_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDko[ row_Average_three3[[4]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[4]] ], 
                   file1=paste(subdir_9_part3,  "/3D-row_three3_week0_EEDko-VS-row_three3_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week4_EEDheto[ row_Average_three3[[4]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[4]] ], 
                   file1=paste(subdir_9_part3,  "/3D-row_three3_week4_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "")  
)   






MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[5]] ], vector2=row_Average_three3_week0_EEDko[ row_Average_three3[[5]] ], 
                   file1=paste(subdir_9_part3,  "/3E-row_three3_week0_EEDheto-VS-row_three3_week0_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDheto[ row_Average_three3[[5]] ], vector2=row_Average_three3_week4_EEDheto[ row_Average_three3[[5]] ], 
                   file1=paste(subdir_9_part3,  "/3E-row_three3_week0_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week0_EEDko[ row_Average_three3[[5]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[5]] ], 
                   file1=paste(subdir_9_part3,  "/3E-row_three3_week0_EEDko-VS-row_three3_week4_EEDko.txt", sep = "") 
)   

MyHypothesisTest_2(vector1=row_Average_three3_week4_EEDheto[ row_Average_three3[[5]] ], vector2=row_Average_three3_week4_EEDko[ row_Average_three3[[5]] ], 
                   file1=paste(subdir_9_part3,  "/3E-row_three3_week4_EEDheto-VS-row_three3_week4_EEDheto.txt", sep = "")  
)   













###############################################################################
subdir_10_part3 <- paste(Part3_g,  "/10-HeatMap", sep = "")
if( ! file.exists(subdir_10_part3) ) { dir.create(subdir_10_part3) }





column_Average_one1 <-   rbind(column_Average_one1_week0_EEDheto,  column_Average_one1_week0_EEDko, column_Average_one1_week4_EEDheto, column_Average_one1_week4_EEDko)
dim(column_Average_one1)
rownames(column_Average_one1) <- c("week0_EEDheto", "week0_EEDko", "week4_EEDheto" , "week4_EEDko"  )
length( column_Average_one1[column_Average_one1> -1] )
length( column_Average_one1[column_Average_one1>0.5] )
length( column_Average_one1[column_Average_one1<0.3] )

column_Average_one1[column_Average_one1>0.6]  <- 0.6

column_heatmap_one1_min <- min(column_Average_one1)
column_heatmap_one1_max <- max(column_Average_one1)
column_heatmap_one1_ave <- (column_heatmap_one1_max + column_heatmap_one1_min)/2

zp1_A <- ggplot( melt(column_Average_one1), aes(x = as.numeric(Var2), y = as.factor(Var1), fill = value) )   ## as.numeric  or  as.factor
zp1_A <- zp1_A + geom_tile()  + xlab("Relative distance (kb)") + ylab("Samples") +  ggtitle("Down-regulated Genes") 
zp1_A <- zp1_A + scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint =column_heatmap_one1_ave,  limits=c(column_heatmap_one1_min, column_heatmap_one1_max) ) 
zp1_A <- zp1_A + scale_y_discrete(expand = c(0, 0))  + scale_x_continuous( breaks=c(1, 100, 200, 300, 400, 500, 600 ), labels=c("-3",  "-2",  "-1",   "TSS",    "1",   "2",  "3" ) ,  expand = c(0, 0)  )     ## discrete  or  continuous  
zp1_A <- zp1_A +   MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )
MySaveGgplot2_1(ggplot2Figure1=zp1_A,  path1=subdir_10_part3, fileName1="1A-one1-4Samples-heatmap-1",  height1=3, width1=6)

zp1_B <- ggplot( melt(column_Average_one1), aes(x = as.numeric(Var2), y = as.factor(Var1), fill = value) )   ## as.numeric  or  as.factor
zp1_B <- zp1_B + geom_tile()  + xlab("Relative distance (kb)") + ylab("Samples") +  ggtitle("Down-regulated Genes")
zp1_B <- zp1_B + scale_fill_gradient2( low = "yellow", mid = "yellowgreen", high = "green", midpoint =column_heatmap_one1_ave,  limits=c(column_heatmap_one1_min, column_heatmap_one1_max) ) 
zp1_B <- zp1_B + scale_y_discrete(expand = c(0, 0))  + scale_x_continuous( breaks=c(1, 100, 200, 300, 400, 500, 600 ), labels=c("-3",  "-2",  "-1",   "TSS",    "1",   "2",  "3" ) ,  expand = c(0, 0)  )     ## discrete  or  continuous  
zp1_B <- zp1_B +   MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) 
MySaveGgplot2_1(ggplot2Figure1=zp1_B,  path1=subdir_10_part3, fileName1="1B-one1-4Samples-heatmap-2",  height1=3, width1=6)














column_Average_two2 <-   rbind(column_Average_two2_week0_EEDheto,  column_Average_two2_week0_EEDko, column_Average_two2_week4_EEDheto, column_Average_two2_week4_EEDko)
dim(column_Average_two2)
rownames(column_Average_two2) <- c("week0_EEDheto", "week0_EEDko", "week4_EEDheto" , "week4_EEDko"  )
length( column_Average_two2[column_Average_two2> -1] )
length( column_Average_two2[column_Average_two2>0.5] )
length( column_Average_two2[column_Average_two2<0.3] )

column_Average_two2[column_Average_two2>0.6]  <- 0.6

column_heatmap_two2_min <- min(column_Average_two2)
column_heatmap_two2_max <- max(column_Average_two2)
column_heatmap_two2_ave <- (column_heatmap_two2_max + column_heatmap_two2_min)/2

zp1_A <- ggplot( melt(column_Average_two2), aes(x = as.numeric(Var2), y = as.factor(Var1), fill = value) )   ## as.numeric  or  as.factor
zp1_A <- zp1_A + geom_tile()  + xlab("Relative distance (kb)") + ylab("Samples") +  ggtitle("Unchanged Genes") 
zp1_A <- zp1_A + scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint =column_heatmap_two2_ave,  limits=c(column_heatmap_two2_min, column_heatmap_two2_max) ) 
zp1_A <- zp1_A + scale_y_discrete(expand = c(0, 0))  + scale_x_continuous( breaks=c(1, 100, 200, 300, 400, 500, 600 ), labels=c("-3",  "-2",  "-1",   "TSS",    "1",   "2",  "3" ) ,  expand = c(0, 0)  )     ## discrete  or  continuous  
zp1_A <- zp1_A +   MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )
MySaveGgplot2_1(ggplot2Figure1=zp1_A,  path1=subdir_10_part3, fileName1="2A-two2-4Samples-heatmap-1",  height1=3, width1=6)

zp1_B <- ggplot( melt(column_Average_two2), aes(x = as.numeric(Var2), y = as.factor(Var1), fill = value) )   ## as.numeric  or  as.factor
zp1_B <- zp1_B + geom_tile()  + xlab("Relative distance (kb)") + ylab("Samples") +  ggtitle("Unchanged Genes")
zp1_B <- zp1_B + scale_fill_gradient2( low = "yellow", mid = "yellowgreen", high = "green", midpoint =column_heatmap_two2_ave,  limits=c(column_heatmap_two2_min, column_heatmap_two2_max) ) 
zp1_B <- zp1_B + scale_y_discrete(expand = c(0, 0))  + scale_x_continuous( breaks=c(1, 100, 200, 300, 400, 500, 600 ), labels=c("-3",  "-2",  "-1",   "TSS",    "1",   "2",  "3" ) ,  expand = c(0, 0)  )     ## discrete  or  continuous  
zp1_B <- zp1_B +   MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) 
MySaveGgplot2_1(ggplot2Figure1=zp1_B,  path1=subdir_10_part3, fileName1="2B-two2-4Samples-heatmap-2",  height1=3, width1=6)











column_Average_three3 <-   rbind(column_Average_three3_week0_EEDheto,  column_Average_three3_week0_EEDko, column_Average_three3_week4_EEDheto, column_Average_three3_week4_EEDko)
dim(column_Average_three3)
rownames(column_Average_three3) <- c("week0_EEDheto", "week0_EEDko", "week4_EEDheto" , "week4_EEDko"  )
length( column_Average_three3[column_Average_three3> -1] )
length( column_Average_three3[column_Average_three3>0.5] )
length( column_Average_three3[column_Average_three3<0.3] )

column_Average_three3[column_Average_three3>0.6]  <- 0.6

column_heatmap_three3_min <- min(column_Average_three3)
column_heatmap_three3_max <- max(column_Average_three3)
column_heatmap_three3_ave <- (column_heatmap_three3_max + column_heatmap_three3_min)/2

zp1_A <- ggplot( melt(column_Average_three3), aes(x = as.numeric(Var2), y = as.factor(Var1), fill = value) )   ## as.numeric  or  as.factor
zp1_A <- zp1_A + geom_tile()  + xlab("Relative distance (kb)") + ylab("Samples") +  ggtitle("Up-regulated Genes") 
zp1_A <- zp1_A + scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint =column_heatmap_three3_ave,  limits=c(column_heatmap_three3_min, column_heatmap_three3_max) ) 
zp1_A <- zp1_A + scale_y_discrete(expand = c(0, 0))  + scale_x_continuous( breaks=c(1, 100, 200, 300, 400, 500, 600 ), labels=c("-3",  "-2",  "-1",   "TSS",    "1",   "2",  "3" ) ,  expand = c(0, 0)  )     ## discrete  or  continuous  
zp1_A <- zp1_A +   MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )
MySaveGgplot2_1(ggplot2Figure1=zp1_A,  path1=subdir_10_part3, fileName1="3A-three3-4Samples-heatmap-1",  height1=3, width1=6)

zp1_B <- ggplot( melt(column_Average_three3), aes(x = as.numeric(Var2), y = as.factor(Var1), fill = value) )   ## as.numeric  or  as.factor
zp1_B <- zp1_B + geom_tile()  + xlab("Relative distance (kb)") + ylab("Samples") +  ggtitle("Up-regulated Genes")
zp1_B <- zp1_B + scale_fill_gradient2( low = "yellow", mid = "yellowgreen", high = "green", midpoint =column_heatmap_three3_ave,  limits=c(column_heatmap_three3_min, column_heatmap_three3_max) ) 
zp1_B <- zp1_B + scale_y_discrete(expand = c(0, 0))  + scale_x_continuous( breaks=c(1, 100, 200, 300, 400, 500, 600 ), labels=c("-3",  "-2",  "-1",   "TSS",    "1",   "2",  "3" ) ,  expand = c(0, 0)  )     ## discrete  or  continuous  
zp1_B <- zp1_B +   MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) 
MySaveGgplot2_1(ggplot2Figure1=zp1_B,  path1=subdir_10_part3, fileName1="3B-three3-4Samples-heatmap-2",  height1=3, width1=6)







################## for one sample
dim(reduceColumn_Average_one1_week0_EEDheto)
length( reduceColumn_Average_one1_week0_EEDheto[reduceColumn_Average_one1_week0_EEDheto> -1] )
length( reduceColumn_Average_one1_week0_EEDheto[reduceColumn_Average_one1_week0_EEDheto>1] )
length( reduceColumn_Average_one1_week0_EEDheto[reduceColumn_Average_one1_week0_EEDheto<0.1] )

reduceColumn_Average_one1_week0_EEDheto[reduceColumn_Average_one1_week0_EEDheto>1]  <- 1
reduceColumn_Average_one1_week0_EEDheto[reduceColumn_Average_one1_week0_EEDheto<0.1]  <- 0.1

reduceColumn_Average_one1_week0_EEDheto_min <- min(reduceColumn_Average_one1_week0_EEDheto)
reduceColumn_Average_one1_week0_EEDheto_max <- max(reduceColumn_Average_one1_week0_EEDheto)
reduceColumn_Average_one1_week0_EEDheto_ave <- (reduceColumn_Average_one1_week0_EEDheto_max + reduceColumn_Average_one1_week0_EEDheto)/2

zp1_A <- ggplot( melt(reduceColumn_Average_one1_week0_EEDheto), aes(x = as.numeric(Var2), y = as.numeric(Var1), fill = value) )   ## as.numeric  or  as.factor
zp1_A <- zp1_A + geom_tile()  + xlab("Relative distance (kb)") + ylab("Samples") +  ggtitle("Down-regulated Genes") 
zp1_A <- zp1_A + scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint =column_heatmap_three3_ave,  limits=c(column_heatmap_three3_min, column_heatmap_three3_max) ) 
zp1_A <- zp1_A + scale_y_continuous(expand = c(0, 0))  + scale_x_continuous(  expand = c(0, 0)  )     ## discrete  or  continuous  
zp1_A <- zp1_A +   MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )
MySaveGgplot2_1(ggplot2Figure1=zp1_A,  path1=subdir_10_part3, fileName1="10A-one1_week0_EEDheto",  height1=3, width1=6)

zp1_B <- ggplot( melt(reduceColumn_Average_one1_week0_EEDheto), aes(x = as.numeric(Var2), y = as.numeric(Var1), fill = value) )   ## as.numeric  or  as.factor
zp1_B <- zp1_B + geom_tile()  + xlab("Relative distance (kb)") + ylab("Samples") +  ggtitle("Down-regulated Genes")
zp1_B <- zp1_B + scale_fill_gradient2( low = "yellow", mid = "yellowgreen", high = "green", midpoint =column_heatmap_three3_ave,  limits=c(column_heatmap_three3_min, column_heatmap_three3_max) ) 
zp1_B <- zp1_B + scale_y_continuous(expand = c(0, 0))  + scale_x_continuous( expand = c(0, 0)  )     ## discrete  or  continuous  
zp1_B <- zp1_B +   MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) 
MySaveGgplot2_1(ggplot2Figure1=zp1_B,  path1=subdir_10_part3, fileName1="10B-one1_week0_EEDheto",  height1=3, width1=6)







####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################





















