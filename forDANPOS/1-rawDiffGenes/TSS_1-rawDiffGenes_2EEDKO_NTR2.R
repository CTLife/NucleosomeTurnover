#############################################################################################################################
## Part  2:  Read in the input matrix data, and reduce the dimensions (rows or columns) of matrix, and save the information of input files.
##              We should get a n*m matrix that means it contains n rows (n DNA fragment) and m columns (m bin, 1 bin=10bp).
##              The element of matrix is reads density: ChIP-input (If ChIP-input <0, we let it = 0)
#############################################################################################################################




## one1 is down,  two2 is retain, three3 is  up.
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
one1_H3_normal_rep1         <- read.table("downGenes/nor/H3_TSS_heatmap/a.xls",                   header=TRUE,   sep="",   quote = "",   comment.char = "")  
one1_week0_EEDheto_rep1     <- read.table("downGenes/nor/week0heto_rep1_TSS_heatmap/a.xls",       header=TRUE,   sep="",   quote = "",   comment.char = "")  
one1_week0_EEDheto_rep2     <- read.table("downGenes/nor/week0heto_rep2_TSS_heatmap/a.xls",       header=TRUE,   sep="",   quote = "",   comment.char = "")  
one1_week0_EEDko_rep1       <- read.table("downGenes/nor/week0ko_rep1_TSS_heatmap/a.xls",         header=TRUE,   sep="",   quote = "",   comment.char = "")  
one1_week4_EEDheto_rep1     <- read.table("downGenes/nor/week4heto_rep1_TSS_heatmap/a.xls",       header=TRUE,   sep="",   quote = "",   comment.char = "")  
one1_week4_EEDheto_rep2     <- read.table("downGenes/nor/week4heto_rep2_TSS_heatmap/a.xls",       header=TRUE,   sep="",   quote = "",   comment.char = "")  
one1_week4_EEDko_rep1       <- read.table("downGenes/nor/week4ko_rep1_TSS_heatmap/a.xls",         header=TRUE,   sep="",   quote = "",   comment.char = "")  
one1_week4_EEDko_rep2       <- read.table("downGenes/nor/week4ko_rep2_TSS_heatmap/a.xls",         header=TRUE,   sep="",   quote = "",   comment.char = "")  

two2_H3_normal_rep1       <- read.table("retainGenes/nor/H3_TSS_heatmap/a.xls",                 header=TRUE,   sep="",   quote = "",   comment.char = "")  
two2_week0_EEDheto_rep1   <- read.table("retainGenes/nor/week0heto_rep1_TSS_heatmap/a.xls",     header=TRUE,   sep="",   quote = "",   comment.char = "")  
two2_week0_EEDheto_rep2   <- read.table("retainGenes/nor/week0heto_rep2_TSS_heatmap/a.xls",     header=TRUE,   sep="",   quote = "",   comment.char = "")  
two2_week0_EEDko_rep1     <- read.table("retainGenes/nor/week0ko_rep1_TSS_heatmap/a.xls",       header=TRUE,   sep="",   quote = "",   comment.char = "")  
two2_week4_EEDheto_rep1   <- read.table("retainGenes/nor/week4heto_rep1_TSS_heatmap/a.xls",     header=TRUE,   sep="",   quote = "",   comment.char = "")  
two2_week4_EEDheto_rep2   <- read.table("retainGenes/nor/week4heto_rep2_TSS_heatmap/a.xls",     header=TRUE,   sep="",   quote = "",   comment.char = "")  
two2_week4_EEDko_rep1     <- read.table("retainGenes/nor/week4ko_rep1_TSS_heatmap/a.xls",       header=TRUE,   sep="",   quote = "",   comment.char = "")  
two2_week4_EEDko_rep2     <- read.table("retainGenes/nor/week4ko_rep2_TSS_heatmap/a.xls",       header=TRUE,   sep="",   quote = "",   comment.char = "")  

three3_H3_normal_rep1           <- read.table("upGenes/nor/H3_TSS_heatmap/a.xls",                     header=TRUE,   sep="",   quote = "",   comment.char = "")  
three3_week0_EEDheto_rep1       <- read.table("upGenes/nor/week0heto_rep1_TSS_heatmap/a.xls",         header=TRUE,   sep="",   quote = "",   comment.char = "")  
three3_week0_EEDheto_rep2       <- read.table("upGenes/nor/week0heto_rep2_TSS_heatmap/a.xls",         header=TRUE,   sep="",   quote = "",   comment.char = "")  
three3_week0_EEDko_rep1         <- read.table("upGenes/nor/week0ko_rep1_TSS_heatmap/a.xls",           header=TRUE,   sep="",   quote = "",   comment.char = "")  
three3_week4_EEDheto_rep1       <- read.table("upGenes/nor/week4heto_rep1_TSS_heatmap/a.xls",         header=TRUE,   sep="",   quote = "",   comment.char = "")  
three3_week4_EEDheto_rep2       <- read.table("upGenes/nor/week4heto_rep2_TSS_heatmap/a.xls",         header=TRUE,   sep="",   quote = "",   comment.char = "")  
three3_week4_EEDko_rep1         <- read.table("upGenes/nor/week4ko_rep1_TSS_heatmap/a.xls",           header=TRUE,   sep="",   quote = "",   comment.char = "")  
three3_week4_EEDko_rep2         <- read.table("upGenes/nor/week4ko_rep2_TSS_heatmap/a.xls",           header=TRUE,   sep="",   quote = "",   comment.char = "")  




dim(one1_H3_normal_rep1)
dim(one1_week0_EEDheto_rep1)
dim(one1_week0_EEDheto_rep2)  
dim(one1_week0_EEDko_rep1)  
dim(one1_week4_EEDheto_rep1)
dim(one1_week4_EEDheto_rep2)
dim(one1_week4_EEDko_rep1) 
dim(one1_week4_EEDko_rep2)

dim(two2_H3_normal_rep1)
dim(two2_week0_EEDheto_rep1)
dim(two2_week0_EEDheto_rep2)  
dim(two2_week0_EEDko_rep1)  
dim(two2_week4_EEDheto_rep1)
dim(two2_week4_EEDheto_rep2)
dim(two2_week4_EEDko_rep1) 
dim(two2_week4_EEDko_rep2)

dim(three3_H3_normal_rep1)
dim(three3_week0_EEDheto_rep1)
dim(three3_week0_EEDheto_rep2)  
dim(three3_week0_EEDko_rep1)  
dim(three3_week4_EEDheto_rep1)
dim(three3_week4_EEDheto_rep2)
dim(three3_week4_EEDko_rep1) 
dim(three3_week4_EEDko_rep2)



summary( unlist(one1_H3_normal_rep1) )
summary( unlist(one1_week0_EEDheto_rep1) )
summary( unlist(one1_week0_EEDheto_rep2) )  
summary( unlist(one1_week0_EEDko_rep1) )  
summary( unlist(one1_week4_EEDheto_rep1) )
summary( unlist(one1_week4_EEDheto_rep2) )
summary( unlist(one1_week4_EEDko_rep1) ) 
summary( unlist(one1_week4_EEDko_rep2) )

summary( unlist(two2_H3_normal_rep1) )
summary( unlist(two2_week0_EEDheto_rep1) )
summary( unlist(two2_week0_EEDheto_rep2) )  
summary( unlist(two2_week0_EEDko_rep1) )  
summary( unlist(two2_week4_EEDheto_rep1) )
summary( unlist(two2_week4_EEDheto_rep2) )
summary( unlist(two2_week4_EEDko_rep1) ) 
summary( unlist(two2_week4_EEDko_rep2) )

summary( unlist(three3_H3_normal_rep1) )
summary( unlist(three3_week0_EEDheto_rep1) )
summary( unlist(three3_week0_EEDheto_rep2) )  
summary( unlist(three3_week0_EEDko_rep1) )  
summary( unlist(three3_week4_EEDheto_rep1) )
summary( unlist(three3_week4_EEDheto_rep2) )
summary( unlist(three3_week4_EEDko_rep1) ) 
summary( unlist(three3_week4_EEDko_rep2) )


myTestVector1 <- three3_H3_normal_rep1
myTestVector1[15, 38]  <- NA
summary( unlist(myTestVector1) )
myTestVector1[12, 558]  <- Inf
summary( unlist(myTestVector1) )






#### remove 1st-4th columns  and only select 600 columns. (TSS-3kb ~ TSS+3kb)
one1_H3_normal_rep1      <- as.matrix(one1_H3_normal_rep1[, (205:804)])
one1_week0_EEDheto_rep1  <- as.matrix(one1_week0_EEDheto_rep1[, (205:804)])
one1_week0_EEDheto_rep2  <- as.matrix(one1_week0_EEDheto_rep2[, (205:804)])  
one1_week0_EEDko_rep1    <- as.matrix(one1_week0_EEDko_rep1[, (205:804)])  
one1_week4_EEDheto_rep1  <- as.matrix(one1_week4_EEDheto_rep1[, (205:804)])
one1_week4_EEDheto_rep2  <- as.matrix(one1_week4_EEDheto_rep2[, (205:804)])
one1_week4_EEDko_rep1    <- as.matrix(one1_week4_EEDko_rep1[, (205:804)]) 
one1_week4_EEDko_rep2    <- as.matrix(one1_week4_EEDko_rep2[, (205:804)])

two2_H3_normal_rep1      <- as.matrix(two2_H3_normal_rep1[, (205:804)])
two2_week0_EEDheto_rep1  <- as.matrix(two2_week0_EEDheto_rep1[, (205:804)])
two2_week0_EEDheto_rep2  <- as.matrix(two2_week0_EEDheto_rep2[, (205:804)])  
two2_week0_EEDko_rep1    <- as.matrix(two2_week0_EEDko_rep1[, (205:804)])  
two2_week4_EEDheto_rep1  <- as.matrix(two2_week4_EEDheto_rep1[, (205:804)])
two2_week4_EEDheto_rep2  <- as.matrix(two2_week4_EEDheto_rep2[, (205:804)])
two2_week4_EEDko_rep1    <- as.matrix(two2_week4_EEDko_rep1[, (205:804)]) 
two2_week4_EEDko_rep2    <- as.matrix(two2_week4_EEDko_rep2[, (205:804)])

three3_H3_normal_rep1      <- as.matrix(three3_H3_normal_rep1[, (205:804)])
three3_week0_EEDheto_rep1  <- as.matrix(three3_week0_EEDheto_rep1[, (205:804)])
three3_week0_EEDheto_rep2  <- as.matrix(three3_week0_EEDheto_rep2[, (205:804)])  
three3_week0_EEDko_rep1    <- as.matrix(three3_week0_EEDko_rep1[, (205:804)])  
three3_week4_EEDheto_rep1  <- as.matrix(three3_week4_EEDheto_rep1[, (205:804)])
three3_week4_EEDheto_rep2  <- as.matrix(three3_week4_EEDheto_rep2[, (205:804)])
three3_week4_EEDko_rep1    <- as.matrix(three3_week4_EEDko_rep1[, (205:804)]) 
three3_week4_EEDko_rep2    <- as.matrix(three3_week4_EEDko_rep2[, (205:804)])






dim(one1_H3_normal_rep1) 
dim(one1_week0_EEDheto_rep1)
dim(one1_week0_EEDheto_rep2)  
dim(one1_week0_EEDko_rep1)  
dim(one1_week4_EEDheto_rep1)
dim(one1_week4_EEDheto_rep2)
dim(one1_week4_EEDko_rep1) 
dim(one1_week4_EEDko_rep2)

dim(two2_H3_normal_rep1) 
dim(two2_week0_EEDheto_rep1)
dim(two2_week0_EEDheto_rep2)  
dim(two2_week0_EEDko_rep1)  
dim(two2_week4_EEDheto_rep1)
dim(two2_week4_EEDheto_rep2)
dim(two2_week4_EEDko_rep1) 
dim(two2_week4_EEDko_rep2)

dim(three3_H3_normal_rep1) 
dim(three3_week0_EEDheto_rep1)
dim(three3_week0_EEDheto_rep2)  
dim(three3_week0_EEDko_rep1)  
dim(three3_week4_EEDheto_rep1)
dim(three3_week4_EEDheto_rep2)
dim(three3_week4_EEDko_rep1) 
dim(three3_week4_EEDko_rep2)



summary( as.vector(one1_H3_normal_rep1) )
summary( as.vector(one1_week0_EEDheto_rep1) )
summary( as.vector(one1_week0_EEDheto_rep2) )  
summary( as.vector(one1_week0_EEDko_rep1) )  
summary( as.vector(one1_week4_EEDheto_rep1) )
summary( as.vector(one1_week4_EEDheto_rep2) )
summary( as.vector(one1_week4_EEDko_rep1) ) 
summary( as.vector(one1_week4_EEDko_rep2) )

summary( as.vector(two2_H3_normal_rep1) )
summary( as.vector(two2_week0_EEDheto_rep1) )
summary( as.vector(two2_week0_EEDheto_rep2) )  
summary( as.vector(two2_week0_EEDko_rep1) )  
summary( as.vector(two2_week4_EEDheto_rep1) )
summary( as.vector(two2_week4_EEDheto_rep2) )
summary( as.vector(two2_week4_EEDko_rep1) ) 
summary( as.vector(two2_week4_EEDko_rep2) )

summary( as.vector(three3_H3_normal_rep1) )
summary( as.vector(three3_week0_EEDheto_rep1) )
summary( as.vector(three3_week0_EEDheto_rep2) )  
summary( as.vector(three3_week0_EEDko_rep1) )  
summary( as.vector(three3_week4_EEDheto_rep1) )
summary( as.vector(three3_week4_EEDheto_rep2) )
summary( as.vector(three3_week4_EEDko_rep1) ) 
summary( as.vector(three3_week4_EEDko_rep2) )














## average biological replicates
Average_one1_H3_normal      <- (one1_H3_normal_rep1 +  one1_H3_normal_rep1)/2
Average_one1_week0_EEDheto  <- (one1_week0_EEDheto_rep1+ one1_week0_EEDheto_rep2)/2
Average_one1_week0_EEDko    <- (one1_week0_EEDko_rep1+ one1_week0_EEDko_rep1)/2  
Average_one1_week4_EEDheto  <- (one1_week4_EEDheto_rep1+ one1_week4_EEDheto_rep2)/2
Average_one1_week4_EEDko    <- (one1_week4_EEDko_rep1+ one1_week4_EEDko_rep2)/2 

Average_two2_H3_normal      <- (two2_H3_normal_rep1 +  two2_H3_normal_rep1)/2
Average_two2_week0_EEDheto  <- (two2_week0_EEDheto_rep1+ two2_week0_EEDheto_rep2)/2
Average_two2_week0_EEDko    <- (two2_week0_EEDko_rep1+ two2_week0_EEDko_rep1)/2  
Average_two2_week4_EEDheto  <- (two2_week4_EEDheto_rep1+ two2_week4_EEDheto_rep2)/2
Average_two2_week4_EEDko    <- (two2_week4_EEDko_rep1+ two2_week4_EEDko_rep2)/2 

Average_three3_H3_normal      <- (three3_H3_normal_rep1 +  three3_H3_normal_rep1)/2
Average_three3_week0_EEDheto  <- (three3_week0_EEDheto_rep1+ three3_week0_EEDheto_rep2)/2
Average_three3_week0_EEDko    <- (three3_week0_EEDko_rep1+ three3_week0_EEDko_rep1)/2  
Average_three3_week4_EEDheto  <- (three3_week4_EEDheto_rep1+ three3_week4_EEDheto_rep2)/2
Average_three3_week4_EEDko    <- (three3_week4_EEDko_rep1+ three3_week4_EEDko_rep2)/2 


dim(Average_one1_H3_normal) 
dim(Average_one1_week0_EEDheto) 
dim(Average_one1_week0_EEDko) 
dim(Average_one1_week4_EEDheto) 
dim(Average_one1_week4_EEDko) 

dim(Average_two2_H3_normal) 
dim(Average_two2_week0_EEDheto) 
dim(Average_two2_week0_EEDko) 
dim(Average_two2_week4_EEDheto) 
dim(Average_two2_week4_EEDko) 

dim(Average_three3_H3_normal) 
dim(Average_three3_week0_EEDheto) 
dim(Average_three3_week0_EEDko) 
dim(Average_three3_week4_EEDheto) 
dim(Average_three3_week4_EEDko) 



summary( as.vector(Average_one1_H3_normal) ) 
summary( as.vector(Average_one1_week0_EEDheto) ) 
summary( as.vector(Average_one1_week0_EEDko) ) 
summary( as.vector(Average_one1_week4_EEDheto) ) 
summary( as.vector(Average_one1_week4_EEDko) ) 

summary( as.vector(Average_two2_H3_normal) ) 
summary( as.vector(Average_two2_week0_EEDheto) ) 
summary( as.vector(Average_two2_week0_EEDko) ) 
summary( as.vector(Average_two2_week4_EEDheto) ) 
summary( as.vector(Average_two2_week4_EEDko) ) 

summary( as.vector(Average_three3_H3_normal) ) 
summary( as.vector(Average_three3_week0_EEDheto) ) 
summary( as.vector(Average_three3_week0_EEDko) ) 
summary( as.vector(Average_three3_week4_EEDheto) ) 
summary( as.vector(Average_three3_week4_EEDko) ) 









## average every row
row_one1_H3_normal_rep1      <- rowMeans(one1_H3_normal_rep1[, 201:400])
row_one1_week0_EEDheto_rep1  <- rowMeans(one1_week0_EEDheto_rep1[, 201:400])
row_one1_week0_EEDheto_rep2  <- rowMeans(one1_week0_EEDheto_rep2[, 201:400])  
row_one1_week0_EEDko_rep1    <- rowMeans(one1_week0_EEDko_rep1[, 201:400])  
row_one1_week4_EEDheto_rep1  <- rowMeans(one1_week4_EEDheto_rep1[, 201:400])
row_one1_week4_EEDheto_rep2  <- rowMeans(one1_week4_EEDheto_rep2[, 201:400])
row_one1_week4_EEDko_rep1    <- rowMeans(one1_week4_EEDko_rep1[, 201:400]) 
row_one1_week4_EEDko_rep2    <- rowMeans(one1_week4_EEDko_rep2[, 201:400])

row_Average_one1_H3_normal      <- rowMeans(Average_one1_H3_normal[, 201:400])
row_Average_one1_week0_EEDheto  <- rowMeans(Average_one1_week0_EEDheto[, 201:400])
row_Average_one1_week0_EEDko    <- rowMeans(Average_one1_week0_EEDko[, 201:400]) 
row_Average_one1_week4_EEDheto  <- rowMeans(Average_one1_week4_EEDheto[, 201:400])
row_Average_one1_week4_EEDko    <- rowMeans(Average_one1_week4_EEDko[, 201:400])

length(row_one1_H3_normal_rep1)
length(row_one1_week0_EEDheto_rep1)
length(row_one1_week0_EEDheto_rep2)  
length(row_one1_week0_EEDko_rep1) 
length(row_one1_week4_EEDheto_rep1)
length(row_one1_week4_EEDheto_rep2)
length(row_one1_week4_EEDko_rep1)
length(row_one1_week4_EEDko_rep2)

length(row_Average_one1_H3_normal)
length(row_Average_one1_week0_EEDheto)
length(row_Average_one1_week0_EEDko) 
length(row_Average_one1_week4_EEDheto)
length(row_Average_one1_week4_EEDko)




row_two2_H3_normal_rep1      <- rowMeans(two2_H3_normal_rep1[, 201:400])
row_two2_week0_EEDheto_rep1  <- rowMeans(two2_week0_EEDheto_rep1[, 201:400])
row_two2_week0_EEDheto_rep2  <- rowMeans(two2_week0_EEDheto_rep2[, 201:400])  
row_two2_week0_EEDko_rep1    <- rowMeans(two2_week0_EEDko_rep1[, 201:400])  
row_two2_week4_EEDheto_rep1  <- rowMeans(two2_week4_EEDheto_rep1[, 201:400])
row_two2_week4_EEDheto_rep2  <- rowMeans(two2_week4_EEDheto_rep2[, 201:400])
row_two2_week4_EEDko_rep1    <- rowMeans(two2_week4_EEDko_rep1[, 201:400]) 
row_two2_week4_EEDko_rep2    <- rowMeans(two2_week4_EEDko_rep2[, 201:400])

row_Average_two2_H3_normal      <- rowMeans(Average_two2_H3_normal[, 201:400])
row_Average_two2_week0_EEDheto  <- rowMeans(Average_two2_week0_EEDheto[, 201:400])
row_Average_two2_week0_EEDko    <- rowMeans(Average_two2_week0_EEDko[, 201:400]) 
row_Average_two2_week4_EEDheto  <- rowMeans(Average_two2_week4_EEDheto[, 201:400])
row_Average_two2_week4_EEDko    <- rowMeans(Average_two2_week4_EEDko[, 201:400])

length(row_two2_H3_normal_rep1)
length(row_two2_week0_EEDheto_rep1)
length(row_two2_week0_EEDheto_rep2)  
length(row_two2_week0_EEDko_rep1) 
length(row_two2_week4_EEDheto_rep1)
length(row_two2_week4_EEDheto_rep2)
length(row_two2_week4_EEDko_rep1)
length(row_two2_week4_EEDko_rep2)

length(row_Average_two2_H3_normal)
length(row_Average_two2_week0_EEDheto)
length(row_Average_two2_week0_EEDko) 
length(row_Average_two2_week4_EEDheto)
length(row_Average_two2_week4_EEDko)




row_three3_H3_normal_rep1      <- rowMeans(three3_H3_normal_rep1[, 201:400])
row_three3_week0_EEDheto_rep1  <- rowMeans(three3_week0_EEDheto_rep1[, 201:400])
row_three3_week0_EEDheto_rep2  <- rowMeans(three3_week0_EEDheto_rep2[, 201:400])  
row_three3_week0_EEDko_rep1    <- rowMeans(three3_week0_EEDko_rep1[, 201:400])  
row_three3_week4_EEDheto_rep1  <- rowMeans(three3_week4_EEDheto_rep1[, 201:400])
row_three3_week4_EEDheto_rep2  <- rowMeans(three3_week4_EEDheto_rep2[, 201:400])
row_three3_week4_EEDko_rep1    <- rowMeans(three3_week4_EEDko_rep1[, 201:400]) 
row_three3_week4_EEDko_rep2    <- rowMeans(three3_week4_EEDko_rep2[, 201:400])

row_Average_three3_H3_normal      <- rowMeans(Average_three3_H3_normal[, 201:400])
row_Average_three3_week0_EEDheto  <- rowMeans(Average_three3_week0_EEDheto[, 201:400])
row_Average_three3_week0_EEDko    <- rowMeans(Average_three3_week0_EEDko[, 201:400]) 
row_Average_three3_week4_EEDheto  <- rowMeans(Average_three3_week4_EEDheto[, 201:400])
row_Average_three3_week4_EEDko    <- rowMeans(Average_three3_week4_EEDko[, 201:400])

length(row_three3_H3_normal_rep1)
length(row_three3_week0_EEDheto_rep1)
length(row_three3_week0_EEDheto_rep2)  
length(row_three3_week0_EEDko_rep1) 
length(row_three3_week4_EEDheto_rep1)
length(row_three3_week4_EEDheto_rep2)
length(row_three3_week4_EEDko_rep1)
length(row_three3_week4_EEDko_rep2)

length(row_Average_three3_H3_normal)
length(row_Average_three3_week0_EEDheto)
length(row_Average_three3_week0_EEDko) 
length(row_Average_three3_week4_EEDheto)
length(row_Average_three3_week4_EEDko)







## average every column
column_one1_H3_normal_rep1      <- colMeans(one1_H3_normal_rep1)
column_one1_week0_EEDheto_rep1  <- colMeans(one1_week0_EEDheto_rep1)
column_one1_week0_EEDheto_rep2  <- colMeans(one1_week0_EEDheto_rep2)  
column_one1_week0_EEDko_rep1    <- colMeans(one1_week0_EEDko_rep1)  
column_one1_week4_EEDheto_rep1  <- colMeans(one1_week4_EEDheto_rep1)
column_one1_week4_EEDheto_rep2  <- colMeans(one1_week4_EEDheto_rep2)
column_one1_week4_EEDko_rep1    <- colMeans(one1_week4_EEDko_rep1) 
column_one1_week4_EEDko_rep2    <- colMeans(one1_week4_EEDko_rep2)

column_Average_one1_H3_normal      <- colMeans(Average_one1_H3_normal)
column_Average_one1_week0_EEDheto  <- colMeans(Average_one1_week0_EEDheto)
column_Average_one1_week0_EEDko    <- colMeans(Average_one1_week0_EEDko) 
column_Average_one1_week4_EEDheto  <- colMeans(Average_one1_week4_EEDheto)
column_Average_one1_week4_EEDko    <- colMeans(Average_one1_week4_EEDko)

length(column_one1_H3_normal_rep1)
length(column_one1_week0_EEDheto_rep1)
length(column_one1_week0_EEDheto_rep2)  
length(column_one1_week0_EEDko_rep1) 
length(column_one1_week4_EEDheto_rep1)
length(column_one1_week4_EEDheto_rep2)
length(column_one1_week4_EEDko_rep1)
length(column_one1_week4_EEDko_rep2)

length(column_Average_one1_H3_normal)
length(column_Average_one1_week0_EEDheto)
length(column_Average_one1_week0_EEDko) 
length(column_Average_one1_week4_EEDheto)
length(column_Average_one1_week4_EEDko)







column_two2_H3_normal_rep1      <- colMeans(two2_H3_normal_rep1)
column_two2_week0_EEDheto_rep1  <- colMeans(two2_week0_EEDheto_rep1)
column_two2_week0_EEDheto_rep2  <- colMeans(two2_week0_EEDheto_rep2)  
column_two2_week0_EEDko_rep1    <- colMeans(two2_week0_EEDko_rep1)  
column_two2_week4_EEDheto_rep1  <- colMeans(two2_week4_EEDheto_rep1)
column_two2_week4_EEDheto_rep2  <- colMeans(two2_week4_EEDheto_rep2)
column_two2_week4_EEDko_rep1    <- colMeans(two2_week4_EEDko_rep1) 
column_two2_week4_EEDko_rep2    <- colMeans(two2_week4_EEDko_rep2)

column_Average_two2_H3_normal      <- colMeans(Average_two2_H3_normal)
column_Average_two2_week0_EEDheto  <- colMeans(Average_two2_week0_EEDheto)
column_Average_two2_week0_EEDko    <- colMeans(Average_two2_week0_EEDko) 
column_Average_two2_week4_EEDheto  <- colMeans(Average_two2_week4_EEDheto)
column_Average_two2_week4_EEDko    <- colMeans(Average_two2_week4_EEDko)

length(column_two2_H3_normal_rep1)
length(column_two2_week0_EEDheto_rep1)
length(column_two2_week0_EEDheto_rep2)  
length(column_two2_week0_EEDko_rep1) 
length(column_two2_week4_EEDheto_rep1)
length(column_two2_week4_EEDheto_rep2)
length(column_two2_week4_EEDko_rep1)
length(column_two2_week4_EEDko_rep2)

length(column_Average_two2_H3_normal)
length(column_Average_two2_week0_EEDheto)
length(column_Average_two2_week0_EEDko) 
length(column_Average_two2_week4_EEDheto)
length(column_Average_two2_week4_EEDko)






column_three3_H3_normal_rep1      <- colMeans(three3_H3_normal_rep1)
column_three3_week0_EEDheto_rep1  <- colMeans(three3_week0_EEDheto_rep1)
column_three3_week0_EEDheto_rep2  <- colMeans(three3_week0_EEDheto_rep2)  
column_three3_week0_EEDko_rep1    <- colMeans(three3_week0_EEDko_rep1)  
column_three3_week4_EEDheto_rep1  <- colMeans(three3_week4_EEDheto_rep1)
column_three3_week4_EEDheto_rep2  <- colMeans(three3_week4_EEDheto_rep2)
column_three3_week4_EEDko_rep1    <- colMeans(three3_week4_EEDko_rep1) 
column_three3_week4_EEDko_rep2    <- colMeans(three3_week4_EEDko_rep2)

column_Average_three3_H3_normal      <- colMeans(Average_three3_H3_normal)
column_Average_three3_week0_EEDheto  <- colMeans(Average_three3_week0_EEDheto)
column_Average_three3_week0_EEDko    <- colMeans(Average_three3_week0_EEDko) 
column_Average_three3_week4_EEDheto  <- colMeans(Average_three3_week4_EEDheto)
column_Average_three3_week4_EEDko    <- colMeans(Average_three3_week4_EEDko)

length(column_three3_H3_normal_rep1)
length(column_three3_week0_EEDheto_rep1)
length(column_three3_week0_EEDheto_rep2)  
length(column_three3_week0_EEDko_rep1) 
length(column_three3_week4_EEDheto_rep1)
length(column_three3_week4_EEDheto_rep2)
length(column_three3_week4_EEDko_rep1)
length(column_three3_week4_EEDko_rep2)

length(column_Average_three3_H3_normal)
length(column_Average_three3_week0_EEDheto)
length(column_Average_three3_week0_EEDko) 
length(column_Average_three3_week4_EEDheto)
length(column_Average_three3_week4_EEDko)













numOfColumns1 <- ncol(Average_one1_H3_normal)
numOfRows1    <- nrow(Average_one1_H3_normal)  
numOfColumns1
numOfRows1

numOfColumns2 <- ncol(Average_two2_H3_normal)
numOfRows2    <- nrow(Average_two2_H3_normal)  
numOfColumns2
numOfRows2

numOfColumns3 <- ncol(Average_three3_H3_normal)
numOfRows3    <- nrow(Average_three3_H3_normal)  
numOfColumns3
numOfRows3


sink(paste(Part2_g,  "/1-numberOf-columns-rows.txt",           sep = ""))
cat("number Of Columns of one1:", numOfColumns1, "\n")
cat("number Of Rows of one1:",    numOfRows1, "\n\n\n")

cat("number Of Columns of two2:", numOfColumns2, "\n")
cat("number Of Rows of two2:",    numOfRows2, "\n\n\n")

cat("number Of Columns of three3:", numOfColumns3, "\n")
cat("number Of Rows of three3:",    numOfRows3, "\n\n\n")
sink() 











## The standard error of the mean (SEM) 
SEM_one1_week0_EEDheto <- apply(Average_one1_week0_EEDheto,    2, sd)/sqrt( nrow(Average_one1_week0_EEDheto)-1 )  
SEM_one1_week0_EEDko   <- apply(Average_one1_week0_EEDko,      2, sd)/sqrt( nrow(Average_one1_week0_EEDko)-1 )    
SEM_one1_week4_EEDheto <- apply(Average_one1_week4_EEDheto,    2, sd)/sqrt( nrow(Average_one1_week4_EEDheto)-1 )  
SEM_one1_week4_EEDko   <- apply(Average_one1_week4_EEDko,      2, sd)/sqrt( nrow(Average_one1_week4_EEDko)-1 )   
NonZero_one1_1   <- seq(from = 0,  by=20,  length.out=numOfColumns1/20)
NonZero_one1_1[1] <- 1
SEM_one1_week0_EEDheto[-NonZero_one1_1] <- NA 
SEM_one1_week0_EEDko[-NonZero_one1_1]   <- NA   
SEM_one1_week4_EEDheto[-NonZero_one1_1] <- NA  
SEM_one1_week4_EEDko[-NonZero_one1_1]   <- NA   

SEM_two2_week0_EEDheto <- apply(Average_two2_week0_EEDheto,    2, sd)/sqrt( nrow(Average_two2_week0_EEDheto)-1 )  
SEM_two2_week0_EEDko   <- apply(Average_two2_week0_EEDko,      2, sd)/sqrt( nrow(Average_two2_week0_EEDko)-1 )    
SEM_two2_week4_EEDheto <- apply(Average_two2_week4_EEDheto,    2, sd)/sqrt( nrow(Average_two2_week4_EEDheto)-1 )  
SEM_two2_week4_EEDko   <- apply(Average_two2_week4_EEDko,      2, sd)/sqrt( nrow(Average_two2_week4_EEDko)-1 )   
NonZero_two2_1   <- seq(from = 0,  by=20,  length.out=numOfColumns1/20)
NonZero_two2_1[1] <- 1
SEM_two2_week0_EEDheto[-NonZero_two2_1] <- NA 
SEM_two2_week0_EEDko[-NonZero_two2_1]   <- NA   
SEM_two2_week4_EEDheto[-NonZero_two2_1] <- NA  
SEM_two2_week4_EEDko[-NonZero_two2_1]   <- NA   

SEM_three3_week0_EEDheto <- apply(Average_three3_week0_EEDheto,    2, sd)/sqrt( nrow(Average_three3_week0_EEDheto)-1 )  
SEM_three3_week0_EEDko   <- apply(Average_three3_week0_EEDko,      2, sd)/sqrt( nrow(Average_three3_week0_EEDko)-1 )    
SEM_three3_week4_EEDheto <- apply(Average_three3_week4_EEDheto,    2, sd)/sqrt( nrow(Average_three3_week4_EEDheto)-1 )  
SEM_three3_week4_EEDko   <- apply(Average_three3_week4_EEDko,      2, sd)/sqrt( nrow(Average_three3_week4_EEDko)-1 )   
NonZero_three3_1   <- seq(from = 0,  by=20,  length.out=numOfColumns1/20)
NonZero_three3_1[1] <- 1
SEM_three3_week0_EEDheto[-NonZero_three3_1] <- NA 
SEM_three3_week0_EEDko[-NonZero_three3_1]   <- NA   
SEM_three3_week4_EEDheto[-NonZero_three3_1] <- NA  
SEM_three3_week4_EEDko[-NonZero_three3_1]   <- NA   








## reduce rows:
reduceRow1_Average_one1_week0_EEDheto 	<-	reduceMatrixRow(Average_one1_week0_EEDheto,    rowNum_1=5 ) 
reduceRow1_Average_one1_week0_EEDko   	<-	reduceMatrixRow(Average_one1_week0_EEDko,      rowNum_1=5 )   
reduceRow1_Average_one1_week4_EEDheto 	<-	reduceMatrixRow(Average_one1_week4_EEDheto,    rowNum_1=5 ) 
reduceRow1_Average_one1_week4_EEDko      <-	reduceMatrixRow(Average_one1_week4_EEDko,      rowNum_1=5 )    
dim(reduceRow1_Average_one1_week0_EEDheto) 
dim(reduceRow1_Average_one1_week0_EEDko)   
dim(reduceRow1_Average_one1_week4_EEDheto) 
dim(reduceRow1_Average_one1_week4_EEDko)    

reduceRow1_Average_two2_week0_EEDheto 	<-	reduceMatrixRow(Average_two2_week0_EEDheto,    rowNum_1=5 ) 
reduceRow1_Average_two2_week0_EEDko   	<-	reduceMatrixRow(Average_two2_week0_EEDko,      rowNum_1=5 )   
reduceRow1_Average_two2_week4_EEDheto 	<-	reduceMatrixRow(Average_two2_week4_EEDheto,    rowNum_1=5 ) 
reduceRow1_Average_two2_week4_EEDko    <-	reduceMatrixRow(Average_two2_week4_EEDko,      rowNum_1=5 )    
dim(reduceRow1_Average_two2_week0_EEDheto) 
dim(reduceRow1_Average_two2_week0_EEDko)   
dim(reduceRow1_Average_two2_week4_EEDheto) 
dim(reduceRow1_Average_two2_week4_EEDko)    

reduceRow1_Average_three3_week0_EEDheto 	<-	reduceMatrixRow(Average_three3_week0_EEDheto,    rowNum_1=5 ) 
reduceRow1_Average_three3_week0_EEDko   	<-	reduceMatrixRow(Average_three3_week0_EEDko,      rowNum_1=5 )   
reduceRow1_Average_three3_week4_EEDheto 	<-	reduceMatrixRow(Average_three3_week4_EEDheto,    rowNum_1=5 ) 
reduceRow1_Average_three3_week4_EEDko        <-	reduceMatrixRow(Average_three3_week4_EEDko,      rowNum_1=5 )    
dim(reduceRow1_Average_three3_week0_EEDheto) 
dim(reduceRow1_Average_three3_week0_EEDko)   
dim(reduceRow1_Average_three3_week4_EEDheto) 
dim(reduceRow1_Average_three3_week4_EEDko)    







## reduce rows:
reduceRow2_Average_one1_week0_EEDheto 	<-	reduceMatrixRow(Average_one1_week0_EEDheto,    rowNum_1=20 ) 
reduceRow2_Average_one1_week0_EEDko   	<-	reduceMatrixRow(Average_one1_week0_EEDko,      rowNum_1=20 )   
reduceRow2_Average_one1_week4_EEDheto 	<-	reduceMatrixRow(Average_one1_week4_EEDheto,    rowNum_1=20 ) 
reduceRow2_Average_one1_week4_EEDko      <-	reduceMatrixRow(Average_one1_week4_EEDko,      rowNum_1=20 )    
dim(reduceRow2_Average_one1_week0_EEDheto) 
dim(reduceRow2_Average_one1_week0_EEDko)   
dim(reduceRow2_Average_one1_week4_EEDheto) 
dim(reduceRow2_Average_one1_week4_EEDko)    

reduceRow2_Average_two2_week0_EEDheto 	<-	reduceMatrixRow(Average_two2_week0_EEDheto,    rowNum_1=20 ) 
reduceRow2_Average_two2_week0_EEDko   	<-	reduceMatrixRow(Average_two2_week0_EEDko,      rowNum_1=20 )   
reduceRow2_Average_two2_week4_EEDheto 	<-	reduceMatrixRow(Average_two2_week4_EEDheto,    rowNum_1=20 ) 
reduceRow2_Average_two2_week4_EEDko    <-	reduceMatrixRow(Average_two2_week4_EEDko,      rowNum_1=20 )    
dim(reduceRow2_Average_two2_week0_EEDheto) 
dim(reduceRow2_Average_two2_week0_EEDko)   
dim(reduceRow2_Average_two2_week4_EEDheto) 
dim(reduceRow2_Average_two2_week4_EEDko)    

reduceRow2_Average_three3_week0_EEDheto 	<-	reduceMatrixRow(Average_three3_week0_EEDheto,    rowNum_1=20 ) 
reduceRow2_Average_three3_week0_EEDko   	<-	reduceMatrixRow(Average_three3_week0_EEDko,      rowNum_1=20 )   
reduceRow2_Average_three3_week4_EEDheto 	<-	reduceMatrixRow(Average_three3_week4_EEDheto,    rowNum_1=20 ) 
reduceRow2_Average_three3_week4_EEDko        <-	reduceMatrixRow(Average_three3_week4_EEDko,      rowNum_1=20 )    
dim(reduceRow2_Average_three3_week0_EEDheto) 
dim(reduceRow2_Average_three3_week0_EEDko)   
dim(reduceRow2_Average_three3_week4_EEDheto) 
dim(reduceRow2_Average_three3_week4_EEDko)    
















## reduce columns, per bin contains 1000bp.
reduceColumn1_Average_one1_week0_EEDheto 	<-	reduceMatrixCol(Average_one1_week0_EEDheto,    colNum_1=6 ) 
reduceColumn1_Average_one1_week0_EEDko   	<-	reduceMatrixCol(Average_one1_week0_EEDko,      colNum_1=6 )   
reduceColumn1_Average_one1_week4_EEDheto 	<-	reduceMatrixCol(Average_one1_week4_EEDheto,    colNum_1=6 ) 
reduceColumn1_Average_one1_week4_EEDko   	<-	reduceMatrixCol(Average_one1_week4_EEDko,      colNum_1=6 )    
dim(reduceColumn1_Average_one1_week0_EEDheto) 
dim(reduceColumn1_Average_one1_week0_EEDko)   
dim(reduceColumn1_Average_one1_week4_EEDheto) 
dim(reduceColumn1_Average_one1_week4_EEDko)    

reduceColumn1_Average_two2_week0_EEDheto 	<-	reduceMatrixCol(Average_two2_week0_EEDheto,    colNum_1=6 ) 
reduceColumn1_Average_two2_week0_EEDko   	<-	reduceMatrixCol(Average_two2_week0_EEDko,      colNum_1=6 )   
reduceColumn1_Average_two2_week4_EEDheto 	<-	reduceMatrixCol(Average_two2_week4_EEDheto,    colNum_1=6 ) 
reduceColumn1_Average_two2_week4_EEDko  	        <-	reduceMatrixCol(Average_two2_week4_EEDko,      colNum_1=6 )    
dim(reduceColumn1_Average_two2_week0_EEDheto) 
dim(reduceColumn1_Average_two2_week0_EEDko)   
dim(reduceColumn1_Average_two2_week4_EEDheto) 
dim(reduceColumn1_Average_two2_week4_EEDko)    

reduceColumn1_Average_three3_week0_EEDheto 	<-	reduceMatrixCol(Average_three3_week0_EEDheto,    colNum_1=6 ) 
reduceColumn1_Average_three3_week0_EEDko   	<-	reduceMatrixCol(Average_three3_week0_EEDko,      colNum_1=6 )   
reduceColumn1_Average_three3_week4_EEDheto 	<-	reduceMatrixCol(Average_three3_week4_EEDheto,    colNum_1=6 ) 
reduceColumn1_Average_three3_week4_EEDko  	<-	reduceMatrixCol(Average_three3_week4_EEDko,      colNum_1=6 )    
dim(reduceColumn1_Average_three3_week0_EEDheto) 
dim(reduceColumn1_Average_three3_week0_EEDko)   
dim(reduceColumn1_Average_three3_week4_EEDheto) 
dim(reduceColumn1_Average_three3_week4_EEDko)    







## reduce columns, per bin contains 200bp.
reduceColumn2_Average_one1_week0_EEDheto 	<-	reduceMatrixCol(Average_one1_week0_EEDheto,    colNum_1=30 ) 
reduceColumn2_Average_one1_week0_EEDko   	<-	reduceMatrixCol(Average_one1_week0_EEDko,      colNum_1=30 )   
reduceColumn2_Average_one1_week4_EEDheto 	<-	reduceMatrixCol(Average_one1_week4_EEDheto,    colNum_1=30 ) 
reduceColumn2_Average_one1_week4_EEDko   	<-	reduceMatrixCol(Average_one1_week4_EEDko,      colNum_1=30 )    
dim(reduceColumn2_Average_one1_week0_EEDheto) 
dim(reduceColumn2_Average_one1_week0_EEDko)   
dim(reduceColumn2_Average_one1_week4_EEDheto) 
dim(reduceColumn2_Average_one1_week4_EEDko)    

reduceColumn2_Average_two2_week0_EEDheto 	<-	reduceMatrixCol(Average_two2_week0_EEDheto,    colNum_1=30 ) 
reduceColumn2_Average_two2_week0_EEDko   	<-	reduceMatrixCol(Average_two2_week0_EEDko,      colNum_1=30 )   
reduceColumn2_Average_two2_week4_EEDheto 	<-	reduceMatrixCol(Average_two2_week4_EEDheto,    colNum_1=30 ) 
reduceColumn2_Average_two2_week4_EEDko  	<-	reduceMatrixCol(Average_two2_week4_EEDko,      colNum_1=30 )    
dim(reduceColumn2_Average_two2_week0_EEDheto) 
dim(reduceColumn2_Average_two2_week0_EEDko)   
dim(reduceColumn2_Average_two2_week4_EEDheto) 
dim(reduceColumn2_Average_two2_week4_EEDko)    

reduceColumn2_Average_three3_week0_EEDheto 	<-	reduceMatrixCol(Average_three3_week0_EEDheto,    colNum_1=30 ) 
reduceColumn2_Average_three3_week0_EEDko   	<-	reduceMatrixCol(Average_three3_week0_EEDko,      colNum_1=30 )   
reduceColumn2_Average_three3_week4_EEDheto 	<-	reduceMatrixCol(Average_three3_week4_EEDheto,    colNum_1=30 ) 
reduceColumn2_Average_three3_week4_EEDko  	<-	reduceMatrixCol(Average_three3_week4_EEDko,      colNum_1=30 )    
dim(reduceColumn2_Average_three3_week0_EEDheto) 
dim(reduceColumn2_Average_three3_week0_EEDko)   
dim(reduceColumn2_Average_three3_week4_EEDheto) 
dim(reduceColumn2_Average_three3_week4_EEDko)    








col_NOL_max <- max( c(column_one1_week0_EEDheto_rep1, column_one1_week0_EEDheto_rep2, column_one1_week0_EEDko_rep1, 
                      column_one1_week4_EEDheto_rep1, column_one1_week4_EEDheto_rep2,  column_one1_week4_EEDko_rep1, column_one1_week4_EEDko_rep2,
                      column_two2_week0_EEDheto_rep1, column_two2_week0_EEDheto_rep2, column_two2_week0_EEDko_rep1, 
                      column_two2_week4_EEDheto_rep1, column_two2_week4_EEDheto_rep2,  column_two2_week4_EEDko_rep1, column_two2_week4_EEDko_rep2, 
                      column_three3_week0_EEDheto_rep1, column_three3_week0_EEDheto_rep2,  column_three3_week0_EEDko_rep1, 
                      column_three3_week4_EEDheto_rep1, column_three3_week4_EEDheto_rep2,  column_three3_week4_EEDko_rep1, column_three3_week4_EEDko_rep2 ) )

col_NOL_min <- min( c(column_one1_week0_EEDheto_rep1, column_one1_week0_EEDheto_rep2, column_one1_week0_EEDko_rep1, 
                      column_one1_week4_EEDheto_rep1, column_one1_week4_EEDheto_rep2,  column_one1_week4_EEDko_rep1, column_one1_week4_EEDko_rep2,
                      column_two2_week0_EEDheto_rep1, column_two2_week0_EEDheto_rep2, column_two2_week0_EEDko_rep1, 
                      column_two2_week4_EEDheto_rep1, column_two2_week4_EEDheto_rep2,  column_two2_week4_EEDko_rep1, column_two2_week4_EEDko_rep2, 
                      column_three3_week0_EEDheto_rep1, column_three3_week0_EEDheto_rep2,  column_three3_week0_EEDko_rep1, 
                      column_three3_week4_EEDheto_rep1, column_three3_week4_EEDheto_rep2,  column_three3_week4_EEDko_rep1, column_three3_week4_EEDko_rep2 ) )

col_NOL_max 
col_NOL_min






row_NOL_max <- max( c(row_one1_week0_EEDheto_rep1, row_one1_week0_EEDheto_rep2, row_one1_week0_EEDko_rep1, 
                      row_one1_week4_EEDheto_rep1, row_one1_week4_EEDheto_rep2,  row_one1_week4_EEDko_rep1, row_one1_week4_EEDko_rep2,
                      row_two2_week0_EEDheto_rep1, row_two2_week0_EEDheto_rep2, row_two2_week0_EEDko_rep1, 
                      row_two2_week4_EEDheto_rep1, row_two2_week4_EEDheto_rep2,  row_two2_week4_EEDko_rep1, row_two2_week4_EEDko_rep2, 
                      row_three3_week0_EEDheto_rep1, row_three3_week0_EEDheto_rep2,  row_three3_week0_EEDko_rep1, 
                      row_three3_week4_EEDheto_rep1, row_three3_week4_EEDheto_rep2,  row_three3_week4_EEDko_rep1, row_three3_week4_EEDko_rep2 ) )

row_NOL_min <- min( c(row_one1_week0_EEDheto_rep1, row_one1_week0_EEDheto_rep2, row_one1_week0_EEDko_rep1, 
                      row_one1_week4_EEDheto_rep1, row_one1_week4_EEDheto_rep2,  row_one1_week4_EEDko_rep1, row_one1_week4_EEDko_rep2,
                      row_two2_week0_EEDheto_rep1, row_two2_week0_EEDheto_rep2, row_two2_week0_EEDko_rep1, 
                      row_two2_week4_EEDheto_rep1, row_two2_week4_EEDheto_rep2,  row_two2_week4_EEDko_rep1, row_two2_week4_EEDko_rep2, 
                      row_three3_week0_EEDheto_rep1, row_three3_week0_EEDheto_rep2,  row_three3_week0_EEDko_rep1, 
                      row_three3_week4_EEDheto_rep1, row_three3_week4_EEDheto_rep2,  row_three3_week4_EEDko_rep1, row_three3_week4_EEDko_rep2 ) )

row_NOL_max 
row_NOL_min







sink(paste(Part2_g,  "/2-max-min.txt",           sep = ""))
cat("col_NOL_max:", col_NOL_max, "\n")
cat("col_NOL_min:", col_NOL_min, "\n\n\n")
cat("row_NOL_max:", row_NOL_max, "\n")
cat("row_NOL_min:", row_NOL_min, "\n\n\n")
sink() 











## SMOOTH
Position_1   <-  seq(from = -3,  by=0.01,  length.out=600)   ## unit is kb
column_one1_week0_EEDheto_rep1_smooth <-  ksmooth(x=Position_1,   y=column_one1_week0_EEDheto_rep1,     kernel = "normal",   bandwidth = 0.1)$y
column_one1_week0_EEDheto_rep2_smooth <-  ksmooth(x=Position_1,   y=column_one1_week0_EEDheto_rep2,     kernel = "normal",   bandwidth = 0.1)$y
column_one1_week0_EEDko_rep1_smooth   <-  ksmooth(x=Position_1,   y=column_one1_week0_EEDko_rep1,       kernel = "normal",   bandwidth = 0.1)$y
column_one1_week4_EEDheto_rep1_smooth <-  ksmooth(x=Position_1,   y=column_one1_week4_EEDheto_rep1,     kernel = "normal",   bandwidth = 0.1)$y
column_one1_week4_EEDheto_rep2_smooth <-  ksmooth(x=Position_1,   y=column_one1_week4_EEDheto_rep2,     kernel = "normal",   bandwidth = 0.1)$y
column_one1_week4_EEDko_rep1_smooth   <-  ksmooth(x=Position_1,   y=column_one1_week4_EEDko_rep1,       kernel = "normal",   bandwidth = 0.1)$y
column_one1_week4_EEDko_rep2_smooth   <-  ksmooth(x=Position_1,   y=column_one1_week4_EEDko_rep2,       kernel = "normal",   bandwidth = 0.1)$y
column_Average_one1_week0_EEDheto_smooth <-  ksmooth(x=Position_1,   y=column_Average_one1_week0_EEDheto,     kernel = "normal",   bandwidth = 0.1)$y
column_Average_one1_week0_EEDko_smooth   <-  ksmooth(x=Position_1,   y=column_Average_one1_week0_EEDko,       kernel = "normal",   bandwidth = 0.1)$y
column_Average_one1_week4_EEDheto_smooth <-  ksmooth(x=Position_1,   y=column_Average_one1_week4_EEDheto,     kernel = "normal",   bandwidth = 0.1)$y
column_Average_one1_week4_EEDko_smooth   <-  ksmooth(x=Position_1,   y=column_Average_one1_week4_EEDko,       kernel = "normal",   bandwidth = 0.1)$y

column_two2_week0_EEDheto_rep1_smooth <-  ksmooth(x=Position_1,   y=column_two2_week0_EEDheto_rep1,     kernel = "normal",   bandwidth = 0.1)$y
column_two2_week0_EEDheto_rep2_smooth <-  ksmooth(x=Position_1,   y=column_two2_week0_EEDheto_rep2,     kernel = "normal",   bandwidth = 0.1)$y
column_two2_week0_EEDko_rep1_smooth   <-  ksmooth(x=Position_1,   y=column_two2_week0_EEDko_rep1,       kernel = "normal",   bandwidth = 0.1)$y
column_two2_week4_EEDheto_rep1_smooth <-  ksmooth(x=Position_1,   y=column_two2_week4_EEDheto_rep1,     kernel = "normal",   bandwidth = 0.1)$y
column_two2_week4_EEDheto_rep2_smooth <-  ksmooth(x=Position_1,   y=column_two2_week4_EEDheto_rep2,     kernel = "normal",   bandwidth = 0.1)$y
column_two2_week4_EEDko_rep1_smooth   <-  ksmooth(x=Position_1,   y=column_two2_week4_EEDko_rep1,       kernel = "normal",   bandwidth = 0.1)$y
column_two2_week4_EEDko_rep2_smooth   <-  ksmooth(x=Position_1,   y=column_two2_week4_EEDko_rep2,       kernel = "normal",   bandwidth = 0.1)$y
column_Average_two2_week0_EEDheto_smooth <-  ksmooth(x=Position_1,   y=column_Average_two2_week0_EEDheto,     kernel = "normal",   bandwidth = 0.1)$y
column_Average_two2_week0_EEDko_smooth   <-  ksmooth(x=Position_1,   y=column_Average_two2_week0_EEDko,       kernel = "normal",   bandwidth = 0.1)$y
column_Average_two2_week4_EEDheto_smooth <-  ksmooth(x=Position_1,   y=column_Average_two2_week4_EEDheto,     kernel = "normal",   bandwidth = 0.1)$y
column_Average_two2_week4_EEDko_smooth   <-  ksmooth(x=Position_1,   y=column_Average_two2_week4_EEDko,       kernel = "normal",   bandwidth = 0.1)$y

column_three3_week0_EEDheto_rep1_smooth <-  ksmooth(x=Position_1,   y=column_three3_week0_EEDheto_rep1,     kernel = "normal",   bandwidth = 0.1)$y
column_three3_week0_EEDheto_rep2_smooth <-  ksmooth(x=Position_1,   y=column_three3_week0_EEDheto_rep2,     kernel = "normal",   bandwidth = 0.1)$y
column_three3_week0_EEDko_rep1_smooth   <-  ksmooth(x=Position_1,   y=column_three3_week0_EEDko_rep1,       kernel = "normal",   bandwidth = 0.1)$y
column_three3_week4_EEDheto_rep1_smooth <-  ksmooth(x=Position_1,   y=column_three3_week4_EEDheto_rep1,     kernel = "normal",   bandwidth = 0.1)$y
column_three3_week4_EEDheto_rep2_smooth <-  ksmooth(x=Position_1,   y=column_three3_week4_EEDheto_rep2,     kernel = "normal",   bandwidth = 0.1)$y
column_three3_week4_EEDko_rep1_smooth   <-  ksmooth(x=Position_1,   y=column_three3_week4_EEDko_rep1,       kernel = "normal",   bandwidth = 0.1)$y
column_three3_week4_EEDko_rep2_smooth   <-  ksmooth(x=Position_1,   y=column_three3_week4_EEDko_rep2,       kernel = "normal",   bandwidth = 0.1)$y
column_Average_three3_week0_EEDheto_smooth <-  ksmooth(x=Position_1,   y=column_Average_three3_week0_EEDheto,     kernel = "normal",   bandwidth = 0.1)$y
column_Average_three3_week0_EEDko_smooth   <-  ksmooth(x=Position_1,   y=column_Average_three3_week0_EEDko,       kernel = "normal",   bandwidth = 0.1)$y
column_Average_three3_week4_EEDheto_smooth <-  ksmooth(x=Position_1,   y=column_Average_three3_week4_EEDheto,     kernel = "normal",   bandwidth = 0.1)$y
column_Average_three3_week4_EEDko_smooth   <-  ksmooth(x=Position_1,   y=column_Average_three3_week4_EEDko,       kernel = "normal",   bandwidth = 0.1)$y


col_NOL_max2 <-max( c(column_one1_week0_EEDheto_rep1_smooth, column_one1_week0_EEDheto_rep2_smooth, column_one1_week0_EEDko_rep1_smooth, 
               column_one1_week4_EEDheto_rep1_smooth, column_one1_week4_EEDheto_rep2_smooth,  column_one1_week4_EEDko_rep1_smooth, column_one1_week4_EEDko_rep2_smooth,
               column_two2_week0_EEDheto_rep1_smooth, column_two2_week0_EEDheto_rep2_smooth, column_two2_week0_EEDko_rep1_smooth, 
               column_two2_week4_EEDheto_rep1_smooth, column_two2_week4_EEDheto_rep2_smooth,  column_two2_week4_EEDko_rep1_smooth, column_two2_week4_EEDko_rep2_smooth, 
               column_three3_week0_EEDheto_rep1_smooth, column_three3_week0_EEDheto_rep2_smooth, column_three3_week0_EEDko_rep1_smooth, 
               column_three3_week4_EEDheto_rep1_smooth, column_three3_week4_EEDheto_rep2_smooth,  column_three3_week4_EEDko_rep1_smooth, column_three3_week4_EEDko_rep2_smooth ) )

col_NOL_min2 <-min( c(column_one1_week0_EEDheto_rep1_smooth, column_one1_week0_EEDheto_rep2_smooth, column_one1_week0_EEDko_rep1_smooth, 
               column_one1_week4_EEDheto_rep1_smooth, column_one1_week4_EEDheto_rep2_smooth,  column_one1_week4_EEDko_rep1_smooth, column_one1_week4_EEDko_rep2_smooth,
               column_two2_week0_EEDheto_rep1_smooth, column_two2_week0_EEDheto_rep2_smooth, column_two2_week0_EEDko_rep1_smooth, 
               column_two2_week4_EEDheto_rep1_smooth, column_two2_week4_EEDheto_rep2_smooth,  column_two2_week4_EEDko_rep1_smooth, column_two2_week4_EEDko_rep2_smooth, 
               column_three3_week0_EEDheto_rep1_smooth, column_three3_week0_EEDheto_rep2_smooth, column_three3_week0_EEDko_rep1_smooth, 
               column_three3_week4_EEDheto_rep1_smooth, column_three3_week4_EEDheto_rep2_smooth,  column_three3_week4_EEDko_rep1_smooth, column_three3_week4_EEDko_rep2_smooth ) )

col_NOL_max2
col_NOL_min2










####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################




