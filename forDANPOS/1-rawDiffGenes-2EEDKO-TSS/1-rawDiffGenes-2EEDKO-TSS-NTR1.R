#############################################################################################################################
## There are ten parts for analysis the DANPOS2.2.2 processsed results of nucleosome turnover ChIP-seq data:
## Part 1:  Load R libraries and define some functions.
## Part 2:  Reads in the data (Read in the input matrix data), and reduce the dimension (rows or columns) of matrix, 
##          and statistics the information of input files.
## Part 3:  Figures about nucleosome occupancy level (NOL).
## Part 4:  Figures about nucleosome turnover rate (NTR). Correlations between NOL and NTR.
## Part 5:  Same NOL DNA regions, but their NTRs are different.
## Part 6:  Same NTR DNA regions, but their NOLs are different.
## Part 7:  Figures about nucleosome contribution value (NCV).  3D-graph: x=NOL, y=NTR, z=biological function. (based on data)
##          NCV=f(NOL, NTR)=b1/[1+(NOL/k1)^n1] + b2/[1+(k2/NTR)^n2]   (based on Modelling)
## Part 8:  Figures about a smaller matrix, so we can compute NTR at a single DNA region.
## Part 9:  NOL, NTR, NCV and Biological Functions. (3 indexes and biological functions)
## Part 10: Other analysis.
#############################################################################################################################





AllResults_g <- "z-2EEDKO-TSS"
if( ! file.exists(AllResults_g) ) { dir.create(AllResults_g) }

Part1_g   <- paste(AllResults_g,  "/Part1-Functions",      sep = "")
Part2_g   <- paste(AllResults_g,  "/Part2-ReadData",       sep = "")
Part3_g   <- paste(AllResults_g,  "/Part3-NOL",            sep = "")
Part4_g   <- paste(AllResults_g,  "/Part4-NTR",            sep = "")
Part5_g   <- paste(AllResults_g,  "/Part5-sameNOL-NTR",    sep = "")
Part6_g   <- paste(AllResults_g,  "/Part6-sameNTR-NOL",    sep = "")
Part7_g   <- paste(AllResults_g,  "/Part7-NCV",            sep = "")
Part8_g   <- paste(AllResults_g,  "/Part8-SmallMatrix",    sep = "")
Part9_g   <- paste(AllResults_g,  "/Part9-3index-BioFun",  sep = "")
Part10_g  <- paste(AllResults_g,  "/Part10-Others",        sep = "")
if( ! file.exists(Part1_g) )  { dir.create(Part1_g) }
if( ! file.exists(Part2_g) )  { dir.create(Part2_g) }
if( ! file.exists(Part3_g) )  { dir.create(Part3_g) }
if( ! file.exists(Part4_g) )  { dir.create(Part4_g) }
if( ! file.exists(Part5_g) )  { dir.create(Part5_g) }
if( ! file.exists(Part6_g) )  { dir.create(Part6_g) }
if( ! file.exists(Part7_g) )  { dir.create(Part7_g) }
if( ! file.exists(Part8_g) )  { dir.create(Part8_g) }
if( ! file.exists(Part9_g) )  { dir.create(Part9_g) }
if( ! file.exists(Part10_g))  { dir.create(Part10_g)}


library("reshape2")
library("ggplot2") 
library("grid")
library("Cairo")
library("RColorBrewer")
library("gplots")  
library("stats")
library("KernSmooth")
library("psych")
library("minerva")
library("matrixStats")

library("extrafont")
font_import()
fonts()
fonttable()
loadfonts()
loadfonts(device="postscript")
names(postscriptFonts())





MyTheme_1 <- function(textSize1=14, hjust1=NULL, vjust1=NULL,  angle1=NULL) {    #hjust=1,    vjust=1,      angle=30
        theme(  
                line  = element_line(colour="black",  size=1.0,   linetype=1,      lineend=NULL),                                                              ## all line elements.          局部优先总体,下面3个也是,只对非局部设置有效.   所有线属性.
                rect  = element_rect(colour="black",  size=1.0,   linetype=1,      fill="transparent" ),                                                       ## all rectangluar elements.    hjust=1: 靠右对齐.   所有矩形区域属性.
                text  = element_text(family="serif",  face=NULL,  colour="black",  size=textSize1, hjust=NULL, vjust=NULL,   angle=NULL, lineheight=NULL),     ## all text elements.           "serif" for a serif font. 所有文本相关属性.
                title = element_text(family="serif",  face=NULL,  colour="black",  size=textSize1, hjust=NULL, vjust=NULL,   angle=NULL, lineheight=NULL),     ## all title elements: plot, axes, legends.    hjust:水平对齐的方向.  所有标题属性.
    
                axis.title        = element_text(family="serif", face=NULL, colour="black", size=textSize1,    hjust=NULL,   vjust=NULL,   angle=NULL,   lineheight=NULL),       ## label of axes (element_text; inherits from text).  horizontal: 水平的, 水平线 
                axis.title.x      = element_text(family="serif", face=NULL, colour="black", size=textSize1,    hjust=NULL,   vjust=-0.5,   angle=NULL,   lineheight=NULL),       ## x axis label (element_text; inherits from axis.title)
                axis.title.y      = element_text(family="serif", face=NULL, colour="black", size=textSize1,    hjust=NULL,   vjust=1.5,    angle=NULL,   lineheight=NULL),       ## y axis label (element_text; inherits from axis.title)
                axis.text         = element_text(family="serif", face=NULL, colour="black", size=textSize1,    hjust=NULL,   vjust=NULL,   angle=NULL,   lineheight=NULL),       ## tick labels along axes (element_text; inherits from text). 坐标轴刻度的标签的属性.                                                         
                axis.text.x       = element_text(family="serif", face=NULL, colour="black", size=textSize1,    hjust=hjust1, vjust=vjust1, angle=angle1, lineheight=angle1),     ## x axis tick labels (element_text; inherits from axis.text)
                axis.text.y       = element_text(family="serif", face=NULL, colour="black", size=textSize1,    hjust=NULL,   vjust=NULL,   angle=NULL,   lineheight=NULL),       ## y axis tick labels (element_text; inherits from axis.text)
                axis.ticks        = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),                                                                            ## tick marks along axes (element_line; inherits from line). 坐标轴刻度线.
                axis.ticks.x      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),                                                                            ## x axis tick marks (element_line; inherits from axis.ticks)
                axis.ticks.y      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),                                                                            ## y axis tick marks (element_line; inherits from axis.ticks)
                axis.ticks.length = grid::unit(2.0,   "mm",   data=NULL),                                                                                                        ## length of tick marks (unit), ‘"mm"’ Millimetres.  10 mm = 1 cm.  刻度线长度
                axis.ticks.margin = grid::unit(1.0,   "mm",   data=NULL),  	                                                                                                 ## space between tick mark and tick label (unit),  ‘"mm"’ Millimetres.  10 mm = 1 cm. 刻度线和刻度标签之间的间距.                                                                           
                axis.line         = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	                                                                 ## lines along axes (element_line; inherits from line). 坐标轴线
                axis.line.x       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	                                                                 ## line along x axis (element_line; inherits from axis.line)
                axis.line.y       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),	                                                                 ## line along y axis (element_line; inherits from axis.line)
    
                legend.background    = element_rect(colour="transparent", size=1, linetype=1, fill="transparent" ), 	      ## background of legend (element_rect; inherits from rect)
                legend.margin        = grid::unit(1, "mm", data=NULL), 	                                                      ## extra space added around legend (unit). linetype=1指的是矩形边框的类型.
                legend.key           = element_rect(colour="transparent", size=2, linetype=1, fill="transparent" ), 	      ## background underneath legend keys. 图例符号. size=1指的是矩形边框的大小.
                legend.key.size      = grid::unit(6,   "mm", data=NULL) , 	                                              ## size of legend keys   (unit; inherits from legend.key.size)
                legend.key.height    = grid::unit(6.5, "mm", data=NULL) , 	                                              ## key background height (unit; inherits from legend.key.size)
                legend.key.width     = grid::unit(8,   "mm", data=NULL) ,                                                     ## key background width  (unit; inherits from legend.key.size)
                legend.text          = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	##legend item labels. 图例文字标签.
                legend.text.align    = 0, 	                    ## alignment of legend labels (number from 0 (left) to 1 (right))
                legend.title         = element_blank(),   	    ## title of legend (element_text; inherits from title)
                legend.title.align   = 0, 	                    ## alignment of legend title (number from 0 (left) to 1 (right))
                legend.position      = "right", 	            ## the position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
                legend.direction     = "vertical",        	    ## layout of items in legends  ("horizontal" or "vertical")   图例排列方向
                legend.justification = "center",      	            ## anchor point for positioning legend inside plot ("center" or two-element numeric vector)  图例居中方式
                legend.box           = NULL, 	                    ## arrangement of multiple legends ("horizontal" or "vertical")  多图例的排列方式
                legend.box.just      = NULL, 	                    ## justification of each legend within the overall bounding box, when there are multiple legends ("top", "bottom", "left", or "right")  多图例的居中方式
    
                panel.background   = element_rect(colour="transparent", size=0.0, linetype=1, fill="transparent" ),   	## background of plotting area, drawn underneath plot (element_rect; inherits from rect)
                panel.border       = element_rect(colour="black", size=0.5, linetype=1, fill=NA ), 	                ## border around plotting area, drawn on top of plot so that it covers tick marks and grid lines. This should be used with fill=NA (element_rect; inherits from rect)
                panel.margin       = grid::unit(1, "mm", data=NULL) , 	                                                ## margin around facet panels (unit)  分面绘图区之间的边距
                panel.grid         = element_blank(), 	                                                                ## grid lines (element_line; inherits from line)  绘图区网格线
                panel.grid.major   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	## major grid lines (element_line; inherits from panel.grid)  主网格线
                panel.grid.minor   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## minor grid lines (element_line; inherits from panel.grid)  次网格线
                panel.grid.major.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	## vertical major grid lines (element_line; inherits from panel.grid.major)
                panel.grid.major.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal major grid lines (element_line; inherits from panel.grid.major)
                panel.grid.minor.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## vertical minor grid lines (element_line; inherits from panel.grid.minor)
                panel.grid.minor.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal minor grid lines (element_line; inherits from panel.grid.minor)
    
                plot.background	= element_rect(colour="transparent", size=NULL, linetype=NULL, fill="transparent" ),                                                ## background of the entire plot (element_rect; inherits from rect)  整个图形的背景
                plot.title      = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=NULL, vjust=1.5,   angle=NULL, lineheight=NULL),    ## plot title (text appearance) (element_text; inherits from title)  图形标题
                plot.margin     = grid::unit(c(5, 5, 5, 5), "mm", data=NULL), 	                                                                                    ## margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
    
                strip.background = element_rect(colour=NULL, size=NULL, linetype=NULL, fill=NULL ), 	                                                        ## background of facet labels (element_rect; inherits from rect)  分面标签背景
                strip.text       = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	## facet labels (element_text; inherits from text)
                strip.text.x     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	## facet labels along horizontal direction (element_text; inherits from strip.text)
                strip.text.y     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL)   	## facet labels along vertical direction (element_text; inherits from strip.text) 
        ) 
} 







MySaveGgplot2_1 <- function(ggplot2Figure1,  path1, fileName1,  height1, width1) {
        SVG1 <- paste(path1,  "/",  "SVG",  sep = "",  collapse = NULL)
        PNG1 <- paste(path1,  "/",  "PNG",  sep = "",  collapse = NULL)
        PDF1 <- paste(path1,  "/",  "PDF",  sep = "",  collapse = NULL)
        EPS1 <- paste(path1,  "/",  "EPS",  sep = "",  collapse = NULL)
        EPS2 <- paste(path1,  "/",  "EPS2", sep = "",  collapse = NULL)
        if( ! file.exists(SVG1) ) { dir.create(SVG1) }
        if( ! file.exists(PNG1) ) { dir.create(PNG1) }
        if( ! file.exists(PDF1) ) { dir.create(PDF1) }
        if( ! file.exists(EPS1) ) { dir.create(EPS1) }
        if( ! file.exists(EPS2) ) { dir.create(EPS2) }
        postscript( file=paste(EPS2,  "/",  fileName1,  ".eps",  sep="",  collapse=NULL),   height =height1, width = width1,  family = "serif",  paper = "special",  onefile = FALSE,  horizontal = FALSE)      
        print( ggplot2Figure1 )
        dev.off()  
        ggsave( filename = paste(SVG1,  "/",  fileName1,  ".svg",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1000 )
        ggsave( filename = paste(PNG1,  "/",  fileName1,  ".png",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1000 )
        ggsave( filename = paste(PDF1,  "/",  fileName1,  ".pdf",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1000 )
        ggsave( filename = paste(EPS1,  "/",  fileName1,  ".eps",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1000,   device=cairo_ps)         
}
                            
  

                       
                                                                                                                                                                                                        


## Postscript function is removed, because of it doesn't support transparency.
MySaveGgplot2_2 <- function(ggplot2Figure1,  path1, fileName1,  height1, width1) {
        SVG1 <- paste(path1,  "/",  "SVG",  sep = "",  collapse = NULL)
        PNG1 <- paste(path1,  "/",  "PNG",  sep = "",  collapse = NULL)
        PDF1 <- paste(path1,  "/",  "PDF",  sep = "",  collapse = NULL)
        EPS1 <- paste(path1,  "/",  "EPS",  sep = "",  collapse = NULL)
        if( ! file.exists(SVG1) ) { dir.create(SVG1) }
        if( ! file.exists(PNG1) ) { dir.create(PNG1) }
        if( ! file.exists(PDF1) ) { dir.create(PDF1) }
        if( ! file.exists(EPS1) ) { dir.create(EPS1) }
        ggsave( filename = paste(SVG1,  "/",  fileName1,  ".svg",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1000 )
        ggsave( filename = paste(PNG1,  "/",  fileName1,  ".png",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1000 )
        ggsave( filename = paste(PDF1,  "/",  fileName1,  ".pdf",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1000 )
        ggsave( filename = paste(EPS1,  "/",  fileName1,  ".eps",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1000,   device=cairo_ps)         
}








MyAverageLines_1 <- function(vector2,  numSample2,  sampleType2, sampleRank2,   colours2,  path2,  fileName2,  title2,  xLab2,  yLab2,   Ymin2=0,   Ymax2=10,   height2=3, width2=4, center2="0") {  
        binLen <- 10    ## Bin size is 10 bp
        binNum <- 1000  ## 1000 * 10 = 10000 bp
        Position_Local   <-  seq(from = -5,  by=0.01,  length.out=binNum)   ## unit is kb
        DataFrame_Local  <-  data.frame( xAxis = c( rep(Position_Local, numSample2) ),      yAxis = vector2,    SampleType = sampleType2 ) 
  
        FigureTemp1 <- ggplot(DataFrame_Local,   aes(x=xAxis, y=yAxis,  color=SampleType) )  +  xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + 
                scale_colour_manual( values=colours2, breaks=sampleRank2  ) + geom_line(size=0.5) +  geom_vline(xintercept=0, lty=2, col="gray45", size=0.5) +   
                scale_x_continuous( breaks=c(-5,  -2.5,  0,  2.5,   5  ), labels=c("-5",  "-2.5", center2,  "2.5",  "5") ) +  
                MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) +  guides( colour = guide_legend(override.aes = list(size=1.5)) )  

        FigureTemp2 <- ggplot(DataFrame_Local,   aes(x=xAxis, y=yAxis,  color=SampleType) )  + xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + ylim(Ymin2, Ymax2 ) +
                scale_colour_manual( values=colours2, breaks=sampleRank2  ) +  geom_line(size=0.5) +  geom_vline(xintercept=0, lty=2, col="gray45", size=0.5) +   
                scale_x_continuous( breaks=c(-5,  -2.5,  0,  2.5,   5  ), labels=c("-5",  "-2.5",  center2,  "2.5",  "5") ) +  
                MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) +  guides( colour = guide_legend(override.aes = list(size=1.5)) )  

        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=fileName2,                                             height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-limitY",  sep="",  collapse=NULL),  height1=height2, width1=width2)
}  





##Error bar一般是standard error，即标准误，反映的是抽样引起的误差。统计量（样本的函数）的标准差。The standard error of the mean (SEM) 
##有些时候也用standard deviation，即标准差，反映的是样本数据之间的离散程度（或整齐程度）。
##两者的联系: S.E.=S.D./sqrt(N-1),  S.E.受样本数量影响很大，而S.D.则受样本数据齐性影响更大。 
## with error-bar:
MyAverageLines_2 <- function(vector2,  SEM2, numSample2,  sampleType2, sampleRank2,   colours2,  path2,  fileName2,  title2,  xLab2,  yLab2,   Ymin2=0,   Ymax2=10,   height2=3, width2=4, center2="0") {  
        binLen <- 10    ## Bin size is 10 bp
        binNum <- 1000  ## 1000 * 10 = 10000 bp
        Position_Local   <-  seq(from = -5,  by=0.01,  length.out=binNum)   ## unit is kb
        DataFrame_Local  <-  data.frame( xAxis = c( rep(Position_Local, numSample2) ),      yAxis = vector2,    SampleType = sampleType2,  Error = SEM2 ) 
        Limits_Local     <- aes( ymax = yAxis+Error,  ymin = yAxis-Error )
  
        FigureTemp1 <- ggplot(DataFrame_Local,   aes(x=xAxis, y=yAxis,  color=SampleType) )  +  xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + 
                scale_colour_manual( values=colours2, breaks=sampleRank2  ) + geom_line(size=0.5) +  geom_vline(xintercept=0, lty=2, col="gray45", size=0.5) +   
                scale_x_continuous( breaks=c(-5,  -2.5,  0,  2.5,   5  ), labels=c("-5",  "-2.5", center2,  "2.5",  "5") ) +  geom_errorbar(Limits_Local,   width=0.03) + 
                MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) +  guides( colour = guide_legend(override.aes = list(size=1.5)) )  
  
        FigureTemp2 <- ggplot(DataFrame_Local,   aes(x=xAxis, y=yAxis,  color=SampleType) )  + xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + ylim(Ymin2, Ymax2 ) +
                scale_colour_manual( values=colours2, breaks=sampleRank2  ) +  geom_line(size=0.5) +  geom_vline(xintercept=0, lty=2, col="gray45", size=0.5) +   
                scale_x_continuous( breaks=c(-5,  -2.5,  0,  2.5,   5  ), labels=c("-5",  "-2.5",  center2,  "2.5",  "5") ) +  geom_errorbar(Limits_Local,   width=0.03) + 
                MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) +  guides( colour = guide_legend(override.aes = list(size=1.5)) )  
  
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=fileName2,                                             height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-limitY",  sep="",  collapse=NULL),  height1=height2, width1=width2)
}  







MyAverageLines_3 <- function(vector2,  numSample2,  sampleType2, sampleRank2,   colours2,  path2,  fileName2,  title2,  xLab2,  yLab2,   Ymin2=0,   Ymax2=10,   height2=3, width2=4, center2="0") {  
        binLen <- 10    ## Bin size is 10 bp
        binNum <- 600   ## 600 * 10 = 6000 bp
        Position_Local   <-  seq(from = -3,  by=0.01,  length.out=binNum)   ## unit is kb
        DataFrame_Local  <-  data.frame( xAxis = c( rep(Position_Local, numSample2) ),      yAxis = vector2,    SampleType = sampleType2 ) 
  
        FigureTemp1 <- ggplot(DataFrame_Local,   aes(x=xAxis, y=yAxis,  color=SampleType) )  +  xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + 
                scale_colour_manual( values=colours2, breaks=sampleRank2  ) + geom_line(size=0.5) +  geom_vline(xintercept=0, lty=2, col="gray45", size=0.5) +   
                scale_x_continuous( breaks=c(-3,  -2, -1, 0,  1, 2,   3  ), labels=c("-3",  "-2", "-1", center2, "1", "2",  "3") ) +  
                MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) +  guides( colour = guide_legend(override.aes = list(size=1.5)) )  
  
        FigureTemp2 <- ggplot(DataFrame_Local,   aes(x=xAxis, y=yAxis,  color=SampleType) )  + xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + ylim(Ymin2, Ymax2 ) +
                scale_colour_manual( values=colours2, breaks=sampleRank2  ) +  geom_line(size=0.5) +  geom_vline(xintercept=0, lty=2, col="gray45", size=0.5) +   
                scale_x_continuous( breaks=c(-3,  -2, -1, 0,  1, 2,   3  ), labels=c("-3",  "-2", "-1", center2, "1", "2",  "3") ) +  
                MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) +  guides( colour = guide_legend(override.aes = list(size=1.5)) )  
  
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=fileName2,                                             height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-limitY",  sep="",  collapse=NULL),  height1=height2, width1=width2)
}  






MyAverageLines_4 <- function(vector2,  SEM2, numSample2,  sampleType2, sampleRank2,   colours2,  path2,  fileName2,  title2,  xLab2,  yLab2,   Ymin2=0,   Ymax2=10,   height2=3, width2=4, center2="0") {  
        binLen <- 10    ## Bin size is 10 bp
        binNum <- 600   ## 600 * 10 = 6000 bp
        Position_Local   <-  seq(from = -3,  by=0.01,  length.out=binNum)   ## unit is kb
        DataFrame_Local  <-  data.frame( xAxis = c( rep(Position_Local, numSample2) ),      yAxis = vector2,    SampleType = sampleType2,  Error = SEM2 ) 
        Limits_Local     <- aes( ymax = yAxis+Error,  ymin = yAxis-Error )
  
        FigureTemp1 <- ggplot(DataFrame_Local,   aes(x=xAxis, y=yAxis,  color=SampleType) )  +  xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + 
                scale_colour_manual( values=colours2, breaks=sampleRank2  ) + geom_line(size=0.5) +  geom_vline(xintercept=0, lty=2, col="gray45", size=0.5) +   
                scale_x_continuous( breaks=c(-3,  -2, -1, 0,  1, 2,   3  ), labels=c("-3",  "-2", "-1", center2, "1", "2",  "3") ) +  geom_errorbar(Limits_Local,   width=0.03) + 
                MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) +  guides( colour = guide_legend(override.aes = list(size=1.5)) )  
  
        FigureTemp2 <- ggplot(DataFrame_Local,   aes(x=xAxis, y=yAxis,  color=SampleType) )  + xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + ylim(Ymin2, Ymax2 ) +
                scale_colour_manual( values=colours2, breaks=sampleRank2  ) +  geom_line(size=0.5) +  geom_vline(xintercept=0, lty=2, col="gray45", size=0.5) +   
                scale_x_continuous( breaks=c(-3,  -2, -1, 0,  1, 2,   3  ), labels=c("-3",  "-2", "-1", center2, "1", "2",  "3") ) +  geom_errorbar(Limits_Local,   width=0.03) + 
                MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) +  guides( colour = guide_legend(override.aes = list(size=1.5)) )  
  
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=fileName2,                                             height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-limitY",  sep="",  collapse=NULL),  height1=height2, width1=width2)
}  









## df contains two columns, the first column (cond_col=1) is sample type, the second column (val_col=2) is value. (must be).
whisk_1 <- function(df, cond_col=1, val_col=2) {  
        require(reshape2)
        condname <- names(df)[cond_col]  ## save the name of the first column.
        names(df)[cond_col] <- "cond" 
        names(df)[val_col]  <- "value"
        b   <- boxplot(value~cond, data=df, plot=FALSE)   
        df2 <- cbind(as.data.frame(b$stats), c("min","lq","m","uq","max"))
        names(df2) <- c(levels(df$cond), "pos")
        df2 <- melt(df2, id="pos", variable.name="cond")
        df2 <- dcast(df2,cond~pos)   
        names(df2)[1] <- condname 
        print(df2)
        df2
}





MyBoxViolinPlot_1 <- function(vector2,   sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,   Ymin2=0, Ymax2=3) { 
        vector2[vector2>Ymax2] <- Ymax2
        vector2[vector2<Ymin2] <- Ymin2
        DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2    ) 
  
        FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType) ) + scale_x_discrete(limits=sampleRank2)  +
                geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1(DataFrame_Local),   width=0.2, size=0.5 ) +
                geom_boxplot( width=0.6,   aes(y=yAxis), outlier.colour="gray45",  outlier.shape=1, outlier.size=0.1, size=0.5, fill=colours2 ) +    
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )

        FigureTemp2 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) +   scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )

        FigureTemp3 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray")   + scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),    width = 0.1, size=0.3, colour = "black") +
                geom_boxplot( aes(y=yAxis),   width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )

        FigureTemp4 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )

        FigureTemp5 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-boxPlot",                     sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-ViolinPlot-3adjust",          sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "-ViolinPlot-noAdjust",         sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "-ViolinPlot-noAdjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "-ViolinPlot-3Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
  
}  







## Add scatter plot
MyBoxViolinPlot_2 <- function(vector2,   sampleType2,  sampleRank2,    path2,   fileName2,    title2,  xLab2,  yLab2,   height2=4,  width2=4,   Ymin2=0, Ymax2=3, alpha2=0.4) {                                                     
        vector2[vector2>Ymax2] <- Ymax2
        vector2[vector2<Ymin2] <- Ymin2
        DataFrame1  <- data.frame(   sampleType=sampleType2,   yAxis=vector2 )  
  
        FigureTemp1 <- ggplot( DataFrame1, aes(x=sampleType) ) +  scale_x_discrete(limits=sampleRank2)  + 
                geom_jitter(aes(y=yAxis), size=1, colour="grey1", alpha=alpha2, position = position_jitter(width=0.25) ) +
                geom_boxplot( alpha=0, width=0.7,   aes(y=yAxis),  outlier.size=0,  size=0.7,  fill=NA, colour="red") +    
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 )  + ylim(Ymin2, Ymax2 )

        FigureTemp2 <- ggplot(DataFrame1, aes(x=sampleType) ) +  scale_x_discrete(limits=sampleRank2)  + 
                geom_jitter(aes(y=yAxis), size=1, colour="grey45", alpha=alpha2, position = position_jitter(width=0.25) ) +
                geom_violin(aes(y=yAxis), fill = NA, colour = "blue", adjust = 3,  alpha=0, size=0.6) +  
                geom_boxplot( alpha=0, width=0.3,   aes(y=yAxis),  outlier.size=0,  size=0.6,  fill=NA, colour="red" ) +    
                stat_summary(   aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=0.0, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 )  + ylim(Ymin2, Ymax2 )
  
        FigureTemp3 <- ggplot(DataFrame1, aes(x=sampleType) ) + scale_x_discrete(limits=sampleRank2)  +  
                geom_jitter(aes(y=yAxis), size=1, colour="grey45", alpha=alpha2, position = position_jitter(width=0.25) ) +
                geom_violin(aes(y=yAxis), fill = NA, colour = "blue",   alpha=0, size=0.6) +  
                geom_boxplot( alpha=0, width=0.3,   aes(y=yAxis),  outlier.size=0,  size=0.6,  fill=NA, colour="red" ) +    
                stat_summary(   aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=0.0, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 )  + ylim(Ymin2, Ymax2 )
  
        MySaveGgplot2_2(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-ScatterBoxPlot",        sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_2(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-ViolinPlot-3adjust",    sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_2(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "-ViolinPlot-noAdjust",   sep="",  collapse=NULL),  height1=height2, width1=width2)
 
}  






## Add path for 11 samples
MyBoxViolinPlot_3_s11 <- function(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3) { 
        vector2[vector2>Ymax2] <- Ymax2
        vector2[vector2<Ymin2] <- Ymin2
        DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2    ) 
        
        dataframePath1 <- data.frame( x=c(1.05, 1.05, 1.95, 1.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath2 <- data.frame( x=c(2.05, 2.05, 2.95, 2.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath3 <- data.frame( x=c(3.05, 3.05, 3.95, 3.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath4 <- data.frame( x=c(4.05, 4.05, 4.95, 4.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath5 <- data.frame( x=c(5.05, 5.05, 5.95, 5.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath6 <- data.frame( x=c(6.05, 6.05, 6.95, 6.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath7 <- data.frame( x=c(7.05, 7.05, 7.95, 7.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath8 <- data.frame( x=c(8.05, 8.05, 8.95, 8.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath9 <- data.frame( x=c(9.05, 9.05, 9.95, 9.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath10<- data.frame( x=c(10.05, 10.05, 10.95, 10.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        
        
        #Path1 = 'geom_path( data = dataframePath1, aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7)'
        #if(numSample2>2) {
        #        Path1 = paste(Path1, ' + geom_path( data = dataframePath2, aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) ',   sep="",  collapse=NULL) 
        #}
        #if(numSample2>3) {
        #        Path1 = paste(Path1, ' + geom_path( data = dataframePath3, aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) ',   sep="",  collapse=NULL) 
        #}
        #if(numSample2>4) {
        #        Path1 = paste(Path1, ' + geom_path( data = dataframePath4, aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) ',   sep="",  collapse=NULL) 
        # }
        #if(numSample2>5) {
        #        Path1 = paste(Path1, ' + geom_path( data = dataframePath5, aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) ',   sep="",  collapse=NULL) 
        #}
        #if(numSample2>6) {
        #       Path1 = paste(Path1, ' + geom_path( data = dataframePath6, aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) ',   sep="",  collapse=NULL)
        #}
        #if(numSample2>7) {
        #        Path1 = paste(Path1, ' + geom_path( data = dataframePath7, aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) ',   sep="",  collapse=NULL)
        #}
        #if(numSample2>8) {
        #        Path1 = paste(Path1, ' + geom_path( data = dataframePath8, aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) ',   sep="",  collapse=NULL)
        #}
        #if(numSample2>9) {
        #        Path1 = paste(Path1, ' + geom_path( data = dataframePath9, aes(x = x, y = y), size=0.3 ) + annotate("text", x=9.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) ',   sep="",  collapse=NULL) 
        #}
        
        FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType) ) + scale_x_discrete(limits=sampleRank2)  +
                geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1(DataFrame_Local),   width=0.2, size=0.5 ) +
                geom_boxplot( width=0.6,   aes(y=yAxis), outlier.colour="gray45",  outlier.shape=1, outlier.size=0.1, size=0.5, fill=colours2 ) +    
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath8,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath9,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=9.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath10, aes(x = x, y = y), size=0.3 ) + annotate("text", x=10.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        
        
        FigureTemp2 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) +   scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath8,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath9,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=9.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath10, aes(x = x, y = y), size=0.3 ) + annotate("text", x=10.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp3 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray")   + scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),    width = 0.1, size=0.3, colour = "black") +
                geom_boxplot( aes(y=yAxis),   width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath8,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath9,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=9.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath10, aes(x = x, y = y), size=0.3 ) + annotate("text", x=10.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp4 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath8,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath9,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=9.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath10, aes(x = x, y = y), size=0.3 ) + annotate("text", x=10.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp5 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath8,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath9,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=9.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath10, aes(x = x, y = y), size=0.3 ) + annotate("text", x=10.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-path-boxPlot",                     sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3adjust",          sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust",         sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
        
}  


## Add path for 10 samples
MyBoxViolinPlot_3_s10 <- function(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3) { 
        vector2[vector2>Ymax2] <- Ymax2
        vector2[vector2<Ymin2] <- Ymin2
        DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2    ) 
        
        dataframePath1 <- data.frame( x=c(1.05, 1.05, 1.95, 1.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath2 <- data.frame( x=c(2.05, 2.05, 2.95, 2.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath3 <- data.frame( x=c(3.05, 3.05, 3.95, 3.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath4 <- data.frame( x=c(4.05, 4.05, 4.95, 4.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath5 <- data.frame( x=c(5.05, 5.05, 5.95, 5.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath6 <- data.frame( x=c(6.05, 6.05, 6.95, 6.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath7 <- data.frame( x=c(7.05, 7.05, 7.95, 7.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath8 <- data.frame( x=c(8.05, 8.05, 8.95, 8.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath9 <- data.frame( x=c(9.05, 9.05, 9.95, 9.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        
        FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType) ) + scale_x_discrete(limits=sampleRank2)  +
                geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1(DataFrame_Local),   width=0.2, size=0.5 ) +
                geom_boxplot( width=0.6,   aes(y=yAxis), outlier.colour="gray45",  outlier.shape=1, outlier.size=0.1, size=0.5, fill=colours2 ) +    
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath8,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath9,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=9.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        
        
        FigureTemp2 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) +   scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath8,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath9,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=9.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp3 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray")   + scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),    width = 0.1, size=0.3, colour = "black") +
                geom_boxplot( aes(y=yAxis),   width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath8,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath9,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=9.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp4 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath8,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath9,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=9.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp5 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath8,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath9,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=9.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-path-boxPlot",                     sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3adjust",          sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust",         sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
}  


## Add path for 9 samples
MyBoxViolinPlot_3_s9 <- function(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3) { 
        vector2[vector2>Ymax2] <- Ymax2
        vector2[vector2<Ymin2] <- Ymin2
        DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2    ) 
        
        dataframePath1 <- data.frame( x=c(1.05, 1.05, 1.95, 1.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath2 <- data.frame( x=c(2.05, 2.05, 2.95, 2.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath3 <- data.frame( x=c(3.05, 3.05, 3.95, 3.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath4 <- data.frame( x=c(4.05, 4.05, 4.95, 4.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath5 <- data.frame( x=c(5.05, 5.05, 5.95, 5.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath6 <- data.frame( x=c(6.05, 6.05, 6.95, 6.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath7 <- data.frame( x=c(7.05, 7.05, 7.95, 7.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath8 <- data.frame( x=c(8.05, 8.05, 8.95, 8.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        
        FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType) ) + scale_x_discrete(limits=sampleRank2)  +
                geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1(DataFrame_Local),   width=0.2, size=0.5 ) +
                geom_boxplot( width=0.6,   aes(y=yAxis), outlier.colour="gray45",  outlier.shape=1, outlier.size=0.1, size=0.5, fill=colours2 ) +    
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath8,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        
        
        FigureTemp2 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) +   scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath8,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp3 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray")   + scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),    width = 0.1, size=0.3, colour = "black") +
                geom_boxplot( aes(y=yAxis),   width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath8,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp4 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath8,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp5 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath8,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=8.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-path-boxPlot",                     sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3adjust",          sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust",         sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
}  


## Add path for 8 samples
MyBoxViolinPlot_3_s8 <- function(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3) { 
        vector2[vector2>Ymax2] <- Ymax2
        vector2[vector2<Ymin2] <- Ymin2
        DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2    ) 
        
        dataframePath1 <- data.frame( x=c(1.05, 1.05, 1.95, 1.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath2 <- data.frame( x=c(2.05, 2.05, 2.95, 2.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath3 <- data.frame( x=c(3.05, 3.05, 3.95, 3.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath4 <- data.frame( x=c(4.05, 4.05, 4.95, 4.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath5 <- data.frame( x=c(5.05, 5.05, 5.95, 5.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath6 <- data.frame( x=c(6.05, 6.05, 6.95, 6.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath7 <- data.frame( x=c(7.05, 7.05, 7.95, 7.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        
        FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType) ) + scale_x_discrete(limits=sampleRank2)  +
                geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1(DataFrame_Local),   width=0.2, size=0.5 ) +
                geom_boxplot( width=0.6,   aes(y=yAxis), outlier.colour="gray45",  outlier.shape=1, outlier.size=0.1, size=0.5, fill=colours2 ) +    
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        
        FigureTemp2 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) +   scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp3 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray")   + scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),    width = 0.1, size=0.3, colour = "black") +
                geom_boxplot( aes(y=yAxis),   width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp4 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp5 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath7,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=7.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-path-boxPlot",                     sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3adjust",          sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust",         sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
}  


## Add path for 7 samples
MyBoxViolinPlot_3_s7 <- function(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3) { 
        vector2[vector2>Ymax2] <- Ymax2
        vector2[vector2<Ymin2] <- Ymin2
        DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2    ) 
        
        dataframePath1 <- data.frame( x=c(1.05, 1.05, 1.95, 1.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath2 <- data.frame( x=c(2.05, 2.05, 2.95, 2.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath3 <- data.frame( x=c(3.05, 3.05, 3.95, 3.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath4 <- data.frame( x=c(4.05, 4.05, 4.95, 4.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath5 <- data.frame( x=c(5.05, 5.05, 5.95, 5.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath6 <- data.frame( x=c(6.05, 6.05, 6.95, 6.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        
        FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType) ) + scale_x_discrete(limits=sampleRank2)  +
                geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1(DataFrame_Local),   width=0.2, size=0.5 ) +
                geom_boxplot( width=0.6,   aes(y=yAxis), outlier.colour="gray45",  outlier.shape=1, outlier.size=0.1, size=0.5, fill=colours2 ) +    
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp2 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) +   scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp3 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray")   + scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),    width = 0.1, size=0.3, colour = "black") +
                geom_boxplot( aes(y=yAxis),   width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp4 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp5 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath6,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=6.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-path-boxPlot",                     sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3adjust",          sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust",         sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
}  



## Add path for 6 samples
MyBoxViolinPlot_3_s6 <- function(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3) { 
        vector2[vector2>Ymax2] <- Ymax2
        vector2[vector2<Ymin2] <- Ymin2
        DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2    ) 
        
        dataframePath1 <- data.frame( x=c(1.05, 1.05, 1.95, 1.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath2 <- data.frame( x=c(2.05, 2.05, 2.95, 2.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath3 <- data.frame( x=c(3.05, 3.05, 3.95, 3.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath4 <- data.frame( x=c(4.05, 4.05, 4.95, 4.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath5 <- data.frame( x=c(5.05, 5.05, 5.95, 5.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        
        FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType) ) + scale_x_discrete(limits=sampleRank2)  +
                geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1(DataFrame_Local),   width=0.2, size=0.5 ) +
                geom_boxplot( width=0.6,   aes(y=yAxis), outlier.colour="gray45",  outlier.shape=1, outlier.size=0.1, size=0.5, fill=colours2 ) +    
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp2 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) +   scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp3 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray")   + scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),    width = 0.1, size=0.3, colour = "black") +
                geom_boxplot( aes(y=yAxis),   width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp4 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp5 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath5,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=5.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-path-boxPlot",                     sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3adjust",          sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust",         sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
}  



## Add path for 5 samples
MyBoxViolinPlot_3_s5 <- function(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3) { 
        vector2[vector2>Ymax2] <- Ymax2
        vector2[vector2<Ymin2] <- Ymin2
        DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2    ) 
        
        dataframePath1 <- data.frame( x=c(1.05, 1.05, 1.95, 1.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath2 <- data.frame( x=c(2.05, 2.05, 2.95, 2.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath3 <- data.frame( x=c(3.05, 3.05, 3.95, 3.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath4 <- data.frame( x=c(4.05, 4.05, 4.95, 4.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        
        FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType) ) + scale_x_discrete(limits=sampleRank2)  +
                geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1(DataFrame_Local),   width=0.2, size=0.5 ) +
                geom_boxplot( width=0.6,   aes(y=yAxis), outlier.colour="gray45",  outlier.shape=1, outlier.size=0.1, size=0.5, fill=colours2 ) +    
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7)
        
        FigureTemp2 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) +   scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7)
        
        FigureTemp3 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray")   + scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),    width = 0.1, size=0.3, colour = "black") +
                geom_boxplot( aes(y=yAxis),   width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp4 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp5 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath4,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=4.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-path-boxPlot",                     sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3adjust",          sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust",         sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
}  



## Add path for 4 samples
MyBoxViolinPlot_3_s4 <- function(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3) { 
        vector2[vector2>Ymax2] <- Ymax2
        vector2[vector2<Ymin2] <- Ymin2
        DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2    ) 
        
        dataframePath1 <- data.frame( x=c(1.05, 1.05, 1.95, 1.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath2 <- data.frame( x=c(2.05, 2.05, 2.95, 2.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath3 <- data.frame( x=c(3.05, 3.05, 3.95, 3.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        
        FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType) ) + scale_x_discrete(limits=sampleRank2)  +
                geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1(DataFrame_Local),   width=0.2, size=0.5 ) +
                geom_boxplot( width=0.6,   aes(y=yAxis), outlier.colour="gray45",  outlier.shape=1, outlier.size=0.1, size=0.5, fill=colours2 ) +    
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp2 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) +   scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp3 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray")   + scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),    width = 0.1, size=0.3, colour = "black") +
                geom_boxplot( aes(y=yAxis),   width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp4 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp5 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath3,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=3.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-path-boxPlot",                     sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3adjust",          sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust",         sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
}  



## Add path for 3 samples
MyBoxViolinPlot_3_s3 <- function(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3) { 
        vector2[vector2>Ymax2] <- Ymax2
        vector2[vector2<Ymin2] <- Ymin2
        DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2    ) 
        
        dataframePath1 <- data.frame( x=c(1.05, 1.05, 1.95, 1.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        dataframePath2 <- data.frame( x=c(2.05, 2.05, 2.95, 2.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        
        FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType) ) + scale_x_discrete(limits=sampleRank2)  +
                geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1(DataFrame_Local),   width=0.2, size=0.5 ) +
                geom_boxplot( width=0.6,   aes(y=yAxis), outlier.colour="gray45",  outlier.shape=1, outlier.size=0.1, size=0.5, fill=colours2 ) +    
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp2 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) +   scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp3 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray")   + scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),    width = 0.1, size=0.3, colour = "black") +
                geom_boxplot( aes(y=yAxis),   width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp4 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp5 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) +
                geom_path( data = dataframePath2,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=2.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-path-boxPlot",                     sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3adjust",          sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust",         sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
}  



## Add path for 2 samples
MyBoxViolinPlot_3_s2 <- function(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3) { 
        vector2[vector2>Ymax2] <- Ymax2
        vector2[vector2<Ymin2] <- Ymin2
        DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2    ) 
        
        dataframePath1 <- data.frame( x=c(1.05, 1.05, 1.95, 1.95),  y=c(Ymax2+0.05, Ymax2+0.1, Ymax2+0.1, Ymax2+0.05) )  
        
        FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType) ) + scale_x_discrete(limits=sampleRank2)  +
                geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1(DataFrame_Local),   width=0.2, size=0.5 ) +
                geom_boxplot( width=0.6,   aes(y=yAxis), outlier.colour="gray45",  outlier.shape=1, outlier.size=0.1, size=0.5, fill=colours2 ) +    
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp2 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) +   scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp3 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray")   + scale_x_discrete(limits=sampleRank2)  +
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),    width = 0.1, size=0.3, colour = "black") +
                geom_boxplot( aes(y=yAxis),   width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        FigureTemp4 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7)
        
        FigureTemp5 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) + scale_x_discrete(limits=sampleRank2)  +  
                #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show_guide = FALSE) + 
                xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1( hjust1=1, vjust1=1,  angle1=30,   textSize=14 ) + 
                geom_path( data = dataframePath1,  aes(x = x, y = y), size=0.3 ) + annotate("text", x=1.5, y=Ymax2+0.16, label="p<3e-16", size=2.7) 
        
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-path-boxPlot",                     sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3adjust",          sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust",         sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-noAdjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "-path-ViolinPlot-3Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
}  






MyNormalizeLinear_1 <- function(vector2, lower2=-1, upper2=1) {
        max1 = max(vector2)  
        min1 = min(vector2)  
        vector2 = lower2 + (upper2-lower2)*(vector2-min1)/(max1-min1)
        return(vector2)
}
        




MyLmEquation_1 <- function(m) { 
        l <- list( a = format(coef(m)[1], digits = 3, nsmall=3),   b = format(abs(coef(m)[2]), digits = 3, nsmall=4),   r2 = format(summary(m)$r.squared, digits = 3, nsmall=3), pv = format(as.numeric(unlist(summary(m)$coefficients[,4][2] )), digits = 3, nsmall=5) )                                                          
        if (coef(m)[2] >= 0)  {
                eq <- substitute(y == a~+~b* "x," ~~r^2~"="~r2* ","  ~~p~"="~pv, l)  ## substitute函数中需要替换的变量用列表参数方式给出。~ is space.
        } else {
                eq <- substitute(y == a~-~b* "x," ~~r^2~"="~r2* ","  ~~p~"="~pv, l)    
        }
        as.character(as.expression(eq));                 
}





MyTurnoverRate_1 <- function( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2 ) {  ## for any number of time points
        regression_A <- lm(yAxis2 ~ xAxis2, model=TRUE, x=TRUE, y=TRUE, qr=TRUE)    
        print("####################################################################")
        print(regression_A)
        print(summary(regression_A))
        print("####################################################################")
        b = -(coef(regression_A)[2])   ## b is turnover rate, b=0 for no turnover, b>0 for decreasing, b<0 for increasing.
        return(b)
}





MyComputeNCV_1 <- function(x1, x2,  b1, b2,  k1, k2,  n1, n2) {  ## x1 is NOL, x2 is NTR. To consider NOL or NTR only, we can let b2=0 or b1=0.
        temp1 <- (x1/k1)**n1
        temp2 <- (k2/x2)**n2
        NCV1 <- b1/(1+temp1) + b2/(1+temp2)
        return(NCV1)
}
## example: computeNCV_1(x1=1.7, x2=0.9,  b1=2.1, b2=3.0,  k1=0.77, k2=0.76,  n1=1.3, n2=1.1)





MyCor2Vars_1 <- function(vector1, vector2, file1) {
        dataFrame <- data.frame(var1 <- vector1,   var2 <- vector2)
        sink(file=file1)
  
        print("Pearson product-moment correlation coefficient (pairwise) :")
        print( corr.test(x=dataFrame, y = dataFrame,   use = "pairwise", method="pearson",  adjust="holm", alpha=.05, ci=TRUE) ) ## pearson, spearman, kendall
        cat("\n\n\n\n\n")  
  
        print("Pearson product-moment correlation coefficient (complete) :")
        print( corr.test(x=dataFrame, y = dataFrame,   use = "complete", method="pearson",  adjust="holm", alpha=.05, ci=TRUE) ) ## pearson, spearman, kendall
        cat("\n\n\n\n\n")  
  
        print("Spearman's rank correlation coefficient or Spearman's rho (pairwise) :")
        print( corr.test(x=dataFrame, y = dataFrame,   use = "pairwise", method="spearman", adjust="holm", alpha=.05, ci=TRUE) )  ## pearson, spearman, kendall
        cat("\n\n\n\n\n")
  
        print("Spearman's rank correlation coefficient or Spearman's rho (complete) ::")
        print( corr.test(x=dataFrame, y = dataFrame,   use = "complete", method="spearman", adjust="holm", alpha=.05, ci=TRUE) )  ## pearson, spearman, kendall
        cat("\n\n\n\n\n")
  
        print("Kendall rank correlation coefficient, commonly referred to as Kendall's tau (τ) coefficient (pairwise) :")
        print( corr.test(x=dataFrame, y = dataFrame,   use = "pairwise", method="kendall",  adjust="holm", alpha=.05, ci=TRUE) ) ## pearson, spearman, kendall
        cat("\n\n\n\n\n")
  
        print("Kendall rank correlation coefficient, commonly referred to as Kendall's tau (τ) coefficient (complete) :")
        print( corr.test(x=dataFrame, y = dataFrame,   use = "complete", method="kendall",  adjust="holm", alpha=.05, ci=TRUE) ) ## pearson, spearman, kendall
        cat("\n\n\n\n\n")

        print("maximal information coefficient (MIC) in maximal information-based nonparametric exploration (MINE) statistics (alpha=0.6, default):")
        mine(x=vector1, y = vector2,   master=NULL,  alpha=0.6,  C=15,  n.cores=16,  var.thr=1e-5,  eps=NULL)
        cat("\n\n\n\n\n")
  
        print("maximal information coefficient (MIC) in maximal information-based nonparametric exploration (MINE) statistics (alpha=1.0):")
        mine(x=vector1, y = vector2,   master=NULL,  alpha=1.0,  C=15,  n.cores=16,  var.thr=1e-5,  eps=NULL)
        cat("\n\n\n\n\n")
  
        sink()  
}








## T test and Wilcoxon test  (unpaired),   very slowly
MyHypothesisTest_1 <- function(vector1, vector2, file1) {
        sink(file=file1)
        print("######################## Apply continuity correction in the normal approximation for the p-value. ###############################################")

        print("##################################################################################################################################")
        wilcoxTest_2 <- wilcox.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=FALSE,  exact=TRUE, correct=TRUE,  conf.int=TRUE,  conf.level=0.95)
        print( wilcoxTest_2  )
        cat( "\n\nExact p-value:", wilcoxTest_2$p.value, "\n\n\n\n\n" ) 

        print("##################################################################################################################################")
        wilcoxTest_4 <- wilcox.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=FALSE,  exact=TRUE, correct=TRUE,  conf.int=TRUE,  conf.level=0.95)
        print( wilcoxTest_4  )
        cat( "\n\nExact p-value:", wilcoxTest_4$p.value, "\n\n\n\n\n" )
        print("##################################################################################################################################")

        print("##################################################################################################################################")
        wilcoxTest_6 <- wilcox.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=FALSE,  exact=TRUE, correct=TRUE,  conf.int=TRUE,  conf.level=0.95)
        print( wilcoxTest_6  )
        cat( "\n\nExact p-value:", wilcoxTest_6$p.value, "\n\n\n\n\n\n" )
  
  
        print("######################## Don't apply continuity correction in the normal approximation for the p-value. ###############################################")
        print("##################################################################################################################################")

        print("##################################################################################################################################")
        wilcoxTestB_2 <- wilcox.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=FALSE,  exact=TRUE, correct=FALSE,  conf.int=TRUE,  conf.level=0.95)
        print( wilcoxTestB_2  )
        cat( "\n\nExact p-value:", wilcoxTestB_2$p.value, "\n\n\n\n\n" ) 
        print("##################################################################################################################################")

        print("##################################################################################################################################")
        wilcoxTestB_4 <- wilcox.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=FALSE,  exact=TRUE, correct=FALSE,  conf.int=TRUE,  conf.level=0.95)
        print( wilcoxTestB_4  )
        cat( "\n\nExact p-value:", wilcoxTestB_4$p.value, "\n\n\n\n\n" )
        print("##################################################################################################################################")

        print("##################################################################################################################################")
        wilcoxTestB_6 <- wilcox.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=FALSE,  exact=TRUE, correct=FALSE,  conf.int=TRUE,  conf.level=0.95)
        print( wilcoxTestB_6  )
        cat( "\n\nExact p-value:", wilcoxTestB_6$p.value, "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" )
  
  
  
  
  
        print("######################## T-test, var.equal=FALSE. ###############################################")
        print("##################################################################################################################################")

        print("##################################################################################################################################")
        tTest_2 <- t.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=FALSE,   var.equal=FALSE,  conf.level=0.95)
        print( tTest_2  )
        cat( "\n\nExact p-value:", tTest_2$p.value, "\n\n\n\n\n" )
        print("##################################################################################################################################")

        print("##################################################################################################################################")
        tTest_4 <- t.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=FALSE,   var.equal=FALSE,  conf.level=0.95)
        print( tTest_4  )
        cat( "\n\nExact p-value:", tTest_4$p.value, "\n\n\n\n\n" )
        print("##################################################################################################################################")

        print("##################################################################################################################################")
        tTest_6 <- t.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=FALSE,   var.equal=FALSE,  conf.level=0.95)
        print( tTest_6  )
        cat( "\n\nExact p-value:", tTest_6$p.value, "\n\n\n\n\n\n" )
  
  
        print("######################## T-test, var.equal=TRUE. ##############################################################################################")

        print("##################################################################################################################################")
        tTestB_2 <- t.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=FALSE,   var.equal=TRUE,  conf.level=0.95)
        print( tTestB_2  )
        cat( "\n\nExact p-value:", tTestB_2$p.value, "\n\n\n\n\n" )
        print("##################################################################################################################################")

        print("##################################################################################################################################")
        tTestB_4 <- t.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=FALSE,   var.equal=TRUE,  conf.level=0.95)
        print( tTestB_4 )
        cat( "\n\nExact p-value:", tTestB_4$p.value, "\n\n\n\n\n" )

        print("##################################################################################################################################")
        tTestB_6 <- t.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=FALSE,   var.equal=TRUE,  conf.level=0.95)
        print( tTestB_6  )
        cat( "\n\nExact p-value:", tTestB_6$p.value, "\n\n\n\n\n" )
  
        sink()
}







## T test and Wilcoxon test  (paired)
MyHypothesisTest_2 <- function(vector1, vector2, file1) {
        sink(file=file1)
  
        print("######################## Apply continuity correction in the normal approximation for the p-value. ###############################################")
        print("##################################################################################################################################")
        wilcoxTest_1 <- wilcox.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=TRUE,   exact=TRUE, correct=TRUE,  conf.int=TRUE,  conf.level=0.95)
        print( wilcoxTest_1  )
        cat( "\n\nExact p-value:", wilcoxTest_1$p.value, "\n\n\n\n\n" )

        print("##################################################################################################################################")
        wilcoxTest_3 <- wilcox.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=TRUE,   exact=TRUE, correct=TRUE,  conf.int=TRUE,  conf.level=0.95)
        print( wilcoxTest_3  )
        cat( "\n\nExact p-value:", wilcoxTest_3$p.value, "\n\n\n\n\n" ) 

        print("##################################################################################################################################")
        wilcoxTest_5 <- wilcox.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=TRUE,   exact=TRUE, correct=TRUE,  conf.int=TRUE,  conf.level=0.95)
        print( wilcoxTest_5  )
        cat( "\n\nExact p-value:", wilcoxTest_5$p.value, "\n\n\n\n\n" )

  
        print("######################## Don't apply continuity correction in the normal approximation for the p-value. ###############################################")
        print("##################################################################################################################################")
        wilcoxTestB_1 <- wilcox.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=TRUE,   exact=TRUE, correct=FALSE,  conf.int=TRUE,  conf.level=0.95)
        print( wilcoxTestB_1  )
        cat( "\n\nExact p-value:", wilcoxTestB_1$p.value, "\n\n\n\n\n" )

        print("##################################################################################################################################")
        wilcoxTestB_3 <- wilcox.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=TRUE,   exact=TRUE, correct=FALSE,  conf.int=TRUE,  conf.level=0.95)
        print( wilcoxTestB_3 )
        cat( "\n\nExact p-value:", wilcoxTestB_3$p.value, "\n\n\n\n\n" ) 

        print("##################################################################################################################################")
        wilcoxTestB_5 <- wilcox.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=TRUE,   exact=TRUE, correct=FALSE,  conf.int=TRUE,  conf.level=0.95)
        print( wilcoxTestB_5  )
        cat( "\n\nExact p-value:", wilcoxTestB_5$p.value, "\n\n\n\n\n" )

  
  
  
        print("######################## T-test, var.equal=FALSE. ###############################################")
        print("##################################################################################################################################")
        tTest_1 <- t.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=TRUE,    var.equal=FALSE,  conf.level=0.95)
        print( tTest_1  )
        cat( "\n\nExact p-value:", tTest_1$p.value, "\n\n\n\n\n" )

        print("##################################################################################################################################")
        tTest_3 <- t.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=TRUE,    var.equal=FALSE,  conf.level=0.95)
        print( tTest_3  )
        cat( "\n\nExact p-value:", tTest_3$p.value, "\n\n\n\n\n" )

        print("##################################################################################################################################")
        tTest_5 <- t.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=TRUE,    var.equal=FALSE,  conf.level=0.95)
        print( tTest_5  )
        cat( "\n\nExact p-value:", tTest_5$p.value, "\n\n\n\n\n" )

  
        print("######################## T-test, var.equal=TRUE. ##############################################################################################")
        print("##################################################################################################################################")
        tTestB_1 <- t.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=TRUE,    var.equal=TRUE,  conf.level=0.95)
        print( tTestB_1  )
        cat( "\n\nExact p-value:", tTestB_1$p.value, "\n\n\n\n\n" )

        print("##################################################################################################################################")
        tTestB_3 <- t.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=TRUE,    var.equal=TRUE,  conf.level=0.95)
        print( tTestB_3  )
        cat( "\n\nExact p-value:", tTestB_3$p.value, "\n\n\n\n\n" )

        print("##################################################################################################################################")
        tTestB_5 <- t.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=TRUE,    var.equal=TRUE,  conf.level=0.95)
        print( tTestB_5  )
        cat( "\n\nExact p-value:", tTestB_5$p.value, "\n\n\n\n\n" )

        sink() 
}






MyHypothesisTest_3 <- function(vector1,  peakPosition_Local,  file2) {
        Length1<- length(peakPosition_Local)
        up1    <- vector1[ c( (peakPosition_Local[1]-Length1) : (peakPosition_Local[1]-1) ) ]
        peak1  <- vector1[ peakPosition_Local ]
        down1  <- vector1[ c( (peakPosition_Local[Length1-1]+1 )) : (peakPosition_Local[Length1-1]+Length1 )  ]
  
        print("######################################### up  regions  VS  peak regions #########################################################################################")
        MyHypothesisTest_2(vector1=up1, vector2=peak1, file1=paste( file2, "_up-vs-peak", sep = "") )
        cat("\n\n\n\n\n\n\n\n\n\n")  
  
        print("######################################### down  regions  VS  peak regions #########################################################################################")
        MyHypothesisTest_2(vector1=down1, vector2=peak1, file1=paste( file2, "_down-vs-peak",   sep = "") )
        cat("\n\n\n\n\n\n\n\n\n\n")  
  
        print("######################################### up  regions  VS  down regions #########################################################################################")
        MyHypothesisTest_2(vector1=up1, vector2=down1, file1=paste( file2, "_up-vs-down",   sep = "") )
        cat("\n\n\n\n\n\n\n\n\n\n")  
}






## scatter Diagram
MyScatterDiagram_1 <- function(vector2,  path2,   fileName2,  xScale2=0.001,  xLab2="peaks (x1000)",   yLab2,  title2,  height2=4,  width2=4,  yMin2=0, yMax2=2,  alpha2=0.5) {
        vector2[vector2>yMax2] <- yMax2
        vector2[vector2<yMin2] <- yMin2
        dataframeA <- data.frame( xAxis = c(1:length(vector2))*xScale2,  yAxis = vector2 ) 

        FigureTemp1 <- ggplot( data = dataframeA, aes(x = xAxis, y = yAxis) ) + 
                geom_point(size=0.1, colour="grey45", alpha=alpha2, fill="grey45", shape=20 ) +  
                ylim(yMin2, yMax2) +  xlab(xLab2) +   ylab(yLab2) +   ggtitle(title2) + 
                MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL)

        FigureTemp2 <- ggplot( data = dataframeA, aes(x = xAxis, y = yAxis) ) + 
                geom_point(size=0.1, colour="red", alpha=alpha2, fill="red", shape=20 ) +  
                ylim(yMin2, yMax2) +  xlab(xLab2) +   ylab(yLab2) +   ggtitle(title2) + 
                MyTheme_1( textSize=14, hjust1=NULL, vjust1=NULL,  angle1=NULL )

        MySaveGgplot2_2(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-scatterDiagram",        sep="",  collapse=NULL),  height1=height2,  width1=width2)
        MySaveGgplot2_2(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-scatterDiagram-Red",    sep="",  collapse=NULL),  height1=height2,  width1=width2)
}






## relative frequency histogram, number of sample type is 1.
MyHistogram_1 <- function(vector2,  path2,   fileName2,  title2,  xLab2, height2=4,  width2=4,   xMin2=0,  xMax2=1.2,    yMin2=0,  yMax2=0.4) {
        vector2[vector2>xMax2] <- xMax2
        vector2[vector2<xMin2] <- xMin2
        dataframeB  <- data.frame( xAxis = vector2 ) 
  
        FigureTemp1 <- ggplot(data=dataframeB, aes(x=xAxis) )  +  xlab(xLab2) + ylab("Relative Frequency") +  ggtitle(title2)  +  
                ylim(yMin2, yMax2) +  geom_histogram( binwidth = 0.05, aes( y = (..count..)/sum(..count..)), fill="grey45"  )  + 
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )  

        FigureTemp2 <- ggplot(data=dataframeB, aes(x=xAxis) )  +  xlab(xLab2) + ylab("Relative Frequency") +  ggtitle(title2)  + 
                geom_histogram( binwidth = 0.05, aes( y = (..count..)/sum(..count..)), fill="grey45"  )  +  
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) + 
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 ) 

        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-frequency-limitY",    sep="",  collapse=NULL),  height1=height2,  width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-frequency",           sep="",  collapse=NULL),  height1=height2,  width1=width2)
}





## absolute frequency histogram, number of sample type is 1.
MyHistogram_2 <- function(vector2, path2,   fileName2,  title2,  xLab2,  height2=4,  width2=4,   xMin2=0,  xMax2=1.2,   yMin2=0,  yMax2=100) {
        vector2[vector2>xMax2] <- xMax2
        vector2[vector2<xMin2] <- xMin2
        dataframeB  <- data.frame( xAxis = vector2 ) 
  
        FigureTemp1 <- ggplot(data=dataframeB, aes(x=xAxis) )  +  xlab(xLab2) + ylab("Absolute Frequency") +  ggtitle(title2)  +  
                ylim(yMin2, yMax2) +  geom_histogram( binwidth = 0.05, aes( y = (..count..) )  , fill="grey45" )  + 
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 ) 

        FigureTemp2 <- ggplot(data=dataframeB, aes(x=xAxis) )  +  xlab(xLab2) + ylab("Absolute Frequency") +  ggtitle(title2)  +  
                geom_histogram( binwidth = 0.05, aes( y = (..count..) ) , fill="grey45"  )  + 
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 ) 

        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-count-limitY",    sep="",  collapse=NULL),  height1=height2,  width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-count",           sep="",  collapse=NULL),  height1=height2,  width1=width2)
}






## density histogram, number of sample type is 1.   geom_histogram(x, alpha, colour, fill, linetype, size,  weight)
MyHistogram_3 <- function(vector2, path2,   fileName2,  title2,  xLab2,  height2=4,  width2=4,   xMin2=0,  xMax2=1.2,   yMin2=0,  yMax2=10) {
        vector2[vector2>xMax2] <- xMax2
        vector2[vector2<xMin2] <- xMin2
        dataframeB  <- data.frame( xAxis = vector2 ) 
  
        FigureTemp1 <- ggplot(data=dataframeB, aes(x=xAxis) )  +  xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
                ylim(yMin2, yMax2) +  geom_histogram( binwidth = 0.05, aes( y = (..density..) )  )  + 
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 ) 

        FigureTemp2 <- ggplot(data=dataframeB, aes(x=xAxis) )  +  xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
                geom_histogram( binwidth = 0.05, aes( y = (..density..) )  )  + 
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 ) 

        FigureTemp3 <- ggplot(data=dataframeB, aes(x=xAxis) )  +  xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
                ylim(yMin2, yMax2) +  geom_histogram( binwidth = 0.05, aes( y = (..density..), fill = ..count.. )  )  + scale_fill_gradient("Count", low = "green", high = "red") +
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 ) 

        FigureTemp4 <- ggplot(data=dataframeB, aes(x=xAxis) )  +  xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
                geom_histogram( binwidth = 0.05, aes( y = (..density..), fill = ..count.. )  )  + scale_fill_gradient("Count", low = "green", high = "red") +
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 ) 

        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-density-limitY",          sep="",  collapse=NULL),  height1=height2,  width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-density",                 sep="",  collapse=NULL),  height1=height2,  width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "-count-density-limitY",    sep="",  collapse=NULL),  height1=height2,  width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "-count-density",           sep="",  collapse=NULL),  height1=height2,  width1=width2)
}







## stacking
MyHistogram_4 <- function(vector2, sampleType2,  path2,   fileName2,  title2,  xLab2, height2=4,  width2=4,  xMin2=0,  xMax2=1.2,   yMin2=0,  yMax2=100, alpha2=0.5) {
        vector2[vector2>xMax2] <- xMax2
        vector2[vector2<xMin2] <- xMin2
        dataframeB  <- data.frame( xAxis = vector2,  sampleType=sampleType2 ) 
  
        FigureTemp1 <- ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType) )  + geom_bar(binwidth = 0.05 , alpha=alpha2) +  
                xlab(xLab2) + ylab("Absolute frequency") +  ggtitle(title2)  +  
                ylim(yMin2, yMax2) +  geom_histogram( binwidth = 0.05, aes( y = (..count..) ) , alpha=alpha2 )  + 
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   

        FigureTemp2 <- ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType) )  + geom_bar(binwidth = 0.05 , alpha=alpha2) +  
                xlab(xLab2) + ylab("Absolute frequency") +  ggtitle(title2)  +  
                geom_histogram( binwidth = 0.05, aes( y = (..count..) ) , alpha=alpha2 )  + 
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )  

        MySaveGgplot2_2(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-count-limitY",    sep="",  collapse=NULL),  height1=height2,  width1=width2)
        MySaveGgplot2_2(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-count",           sep="",  collapse=NULL),  height1=height2,  width1=width2)
}







## dodge
MyHistogram_5 <- function(vector2, sampleType2,  path2,   fileName2,  title2,  xLab2, height2=4,  width2=4,  xMin2=0,  xMax2=1.2,   yMin2=0,  yMax2=10) {
        vector2[vector2>xMax2] <- xMax2
        vector2[vector2<xMin2] <- xMin2
        dataframeB  <- data.frame( xAxis = vector2,  sampleType=sampleType2 ) 
  
        FigureTemp1 <- ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType) )  + geom_bar(binwidth = 0.05 , position="dodge") +  
                xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
                ylim(yMin2, yMax2)  +  geom_histogram( binwidth = 0.05, aes( y = (..density..) )  )  + 
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   

        FigureTemp2 <- ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType) )  + geom_bar(binwidth = 0.05 , position="dodge") +  
                xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
                geom_histogram( binwidth = 0.05, aes( y = (..density..) )  )  + 
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )  

        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-density-limitY",    sep="",  collapse=NULL),  height1=height2,  width1=width2)
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-density",           sep="",  collapse=NULL),  height1=height2,  width1=width2)
}







## compare  probability  density
MyHistogram_6 <- function(vector2, sampleType2,  colours2,  path2,   fileName2,  title2,  xLab2, height2=4,  width2=4,  xMin2=0,  xMax2=1.5,   yMin2=0,  yMax2=10) {
        vector2[vector2>xMax2] <- xMax2
        vector2[vector2<xMin2] <- xMin2
        dataframeB  <- data.frame( xAxis = vector2,  sampleType=sampleType2 )
  
        
        FigureTemp1 <- { ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType, colour=sampleType) )   +  xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
                geom_density(mapping = NULL, data = NULL, stat = "density", position = "identity", na.rm = FALSE, alpha=0.3  ) +
                scale_colour_manual( values=colours2   ) + scale_fill_manual( values = colours2) +  ylim(yMin2, yMax2) + 
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                #geom_hline(yintercept=-0.01, lty=1, col="white", size=0.6) +
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )    +  guides( colour = guide_legend(override.aes = list(size=0, shape=1)) ) }
        MySaveGgplot2_2(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-density-limitY",      sep="",  collapse=NULL),  height1=height2,  width1=width2)
        
        
        FigureTemp2 <- { ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType, colour=sampleType) )   +  xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
                geom_density(mapping = NULL, data = NULL, stat = "density", position = "identity", na.rm = FALSE, alpha=0.3  ) +
                scale_colour_manual( values=colours2   ) + scale_fill_manual( values = colours2) +  
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                #geom_hline(yintercept=-0.01, lty=1, col="white", size=0.6) +
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )    +  guides( colour = guide_legend(override.aes = list(size=0, shape=1)) ) }
        MySaveGgplot2_2(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-density",             sep="",  collapse=NULL),  height1=height2,  width1=width2)
        
        
        FigureTemp3 <- { ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType, colour=sampleType) )   +  xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
                geom_density(mapping = NULL, data = NULL, stat = "density", position = "identity", na.rm = FALSE, alpha=0.0  ) +
                scale_colour_manual( values=colours2   ) + scale_fill_manual( values = colours2) +  ylim(yMin2, yMax2) + 
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                #geom_hline(yintercept=-0.01, lty=1, col="white", size=0.6) +
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )    +  guides( colour = guide_legend(override.aes = list(size=1.5, shape=1)) ) }
        MySaveGgplot2_2(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "-density2-limitY",     sep="",  collapse=NULL),  height1=height2,  width1=width2)
        
        
        FigureTemp4 <- { ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType, colour=sampleType) )   +  xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
                geom_density(mapping = NULL, data = NULL, stat = "density", position = "identity", na.rm = FALSE, alpha=0.0  ) +
                scale_colour_manual( values=colours2   ) + scale_fill_manual( values = colours2) +  
                scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                #geom_hline(yintercept=-0.01, lty=1, col="white", size=0.6) +
                MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   +  guides( colour = guide_legend(override.aes = list(size=1.5, shape=1)) ) }
        MySaveGgplot2_2(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "-density2",            sep="",  collapse=NULL),  height1=height2,  width1=width2)
        
        
        FigureTemp5 <- { ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType, colour=sampleType) )   +  xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
                        geom_line(stat="density", alpha=1.0 ) +
                        scale_colour_manual( values=colours2   ) + scale_fill_manual( values = colours2) +  ylim(yMin2, yMax2) + 
                        scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                        #geom_hline(yintercept=-0.01, lty=1, col="white", size=0.6) +
                        MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )    +  guides( colour = guide_legend(override.aes = list(size=1.5, shape=1)) ) }
        MySaveGgplot2_2(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "-density3-limitY",     sep="",  collapse=NULL),  height1=height2,  width1=width2)
        
        
        FigureTemp6 <- { ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType, colour=sampleType) )   +  xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
                        geom_line(stat="density", alpha=1.0 ) +
                        scale_colour_manual( values=colours2   ) + scale_fill_manual( values = colours2) +  
                        scale_x_continuous(limits=c(xMin2, xMax2),  breaks=c(0, 0.3,  0.6, 0.9, 1.2), labels=c("0", "0.3",  "0.6", "0.9",   "1.2+") ) +  
                        #geom_hline(yintercept=-0.01, lty=1, col="white", size=0.6) +
                        MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   +  guides( colour = guide_legend(override.aes = list(size=1.5, shape=1)) ) }
        MySaveGgplot2_2(ggplot2Figure1=FigureTemp6,  path1=path2, fileName1=paste(fileName2, "-density3",            sep="",  collapse=NULL),  height1=height2,  width1=width2)
        
}










## reduce columns of a matrix by average the nearest columns.
reduceMatrixCol <- function(matrix_1, colNum_1=10 ) {     ## colNum_1: number of columns after reduced.
        allCol  <- ncol(matrix_1)
        allRow  <- nrow(matrix_1)
        matrix2 <- matrix(nrow = allRow, ncol = colNum_1)
        aveCol  <- allCol/colNum_1
        for(i  in  c(0:(colNum_1-1)) ) {
                index2 <- c( (aveCol*i+1):(aveCol*(i+1)) )
                matrix2[,i+1] <- rowMeans(matrix_1[, index2])
        }
        return(matrix2)
}






## reduce rows of a matrix by average the nearest rows.   
reduceMatrixRow <- function(matrix_1,    rowNum_1=10 ) {     ## rowNum_1: number of rows after reduced.
        allCol  <- ncol(matrix_1)
        allRow  <- nrow(matrix_1)
        matrix2 <- matrix(nrow = rowNum_1, ncol = allCol)
        aveRow  <- allRow/rowNum_1
        for(i  in  c(0:(rowNum_1-1)) )  {
                index2 <- c( (aveRow*i+1):(aveRow*(i+1)) )
                matrix2[i+1,] <- colMeans(matrix_1[index2, ])
        }
        return(matrix2)
}














sink( paste(Part1_g,  "/1-my_all_functions.txt",  sep = "") )
cat('\n
MyTheme_1(textSize1=14, hjust1=NULL, vjust1=NULL,  angle1=NULL)
MySaveGgplot2_1(ggplot2Figure1,  path1, fileName1,  height1, width1)
MySaveGgplot2_2(ggplot2Figure1,  path1, fileName1,  height1, width1)
MyAverageLines_1(vector2,  numSample2,  sampleType2, sampleRank2,   colours2,  path2,  fileName2,  title2,  xLab2,  yLab2,   Ymin2=0,   Ymax2=10,   height2=3, width2=4, center2="0")
MyAverageLines_2(vector2,  SEM2, numSample2,  sampleType2, sampleRank2,   colours2,  path2,  fileName2,  title2,  xLab2,  yLab2,   Ymin2=0,   Ymax2=10,   height2=3, width2=4, center2="0")
MyAverageLines_3(vector2,  numSample2,  sampleType2, sampleRank2,   colours2,  path2,  fileName2,  title2,  xLab2,  yLab2,   Ymin2=0,   Ymax2=10,   height2=3, width2=4, center2="0")
MyAverageLines_4(vector2,  SEM2, numSample2,  sampleType2, sampleRank2,   colours2,  path2,  fileName2,  title2,  xLab2,  yLab2,   Ymin2=0,   Ymax2=10,   height2=3, width2=4, center2="0") 
whisk_1(df, cond_col=1, val_col=2)
MyBoxViolinPlot_1(vector2,   sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,   Ymin2=0, Ymax2=3) 
MyBoxViolinPlot_2(vector2,   sampleType2,  sampleRank2,    path2,   fileName2,    title2,  xLab2,  yLab2,   height2=4,  width2=4,   Ymin2=0, Ymax2=3, alpha2=0.4)
MyBoxViolinPlot_3_s11(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3) 
MyBoxViolinPlot_3_s10(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3) 
MyBoxViolinPlot_3_s9(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3)
MyBoxViolinPlot_3_s8(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3)
MyBoxViolinPlot_3_s7(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3) 
MyBoxViolinPlot_3_s6(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3)
MyBoxViolinPlot_3_s5(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3)
MyBoxViolinPlot_3_s4(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3) 
MyBoxViolinPlot_3_s3(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3)
MyBoxViolinPlot_3_s2(vector2,  sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,  Ymin2=0,  Ymax2=3) 
MyNormalizeLinear_1(vector2, lower2=-1, upper2=1)
MyLmEquation_1(m)
MyTurnoverRate_1( xAxis2=c(0, 1, 2, 4, 6, 8),  yAxis2 )
MyComputeNCV_1(x1, x2,  b1, b2,  k1, k2,  n1, n2)
MyCor2Vars_1(vector1, vector2, file1)
MyHypothesisTest_1(vector1, vector2, file1) 
MyHypothesisTest_2(vector1, vector2, file1) 
MyHypothesisTest_3(vector1,  peakPosition_Local,  file2)
MyScatterDiagram_1(vector2,  path2,   fileName2,  xScale2=0.001,  xLab2="peaks (x1000)",   yLab2,  title2,  height2=4,  width2=4,  yMin2=0, yMax2=2,  alpha2=0.5) 
MyHistogram_1(vector2,  path2,   fileName2,  title2,  xLab2, height2=4,  width2=4,   xMin2=0,  xMax2=1.2,    yMin2=0,  yMax2=0.4) 
MyHistogram_2(vector2, path2,   fileName2,  title2,  xLab2,  height2=4,  width2=4,   xMin2=0,  xMax2=1.2,   yMin2=0,  yMax2=100)
MyHistogram_3(vector2, path2,   fileName2,  title2,  xLab2,  height2=4,  width2=4,   xMin2=0,  xMax2=1.2,   yMin2=0,  yMax2=10)
MyHistogram_4(vector2, sampleType2,  path2,   fileName2,  title2,  xLab2, height2=4,  width2=4,  xMin2=0,  xMax2=1.2,   yMin2=0,  yMax2=100, alpha2=0.5) 
MyHistogram_5(vector2, sampleType2,  path2,   fileName2,  title2,  xLab2, height2=4,  width2=4,  xMin2=0,  xMax2=1.2,   yMin2=0,  yMax2=10) 
MyHistogram_6(vector2, sampleType2,  colours2,  path2,   fileName2,  title2,  xLab2, height2=4,  width2=4,  xMin2=0,  xMax2=1.5,   yMin2=0,  yMax2=10)
reduceMatrixCol(matrix_1, colNum_1=10 )
reduceMatrixRow(matrix_1,    rowNum_1=10 ) 
\n')
sink()
































