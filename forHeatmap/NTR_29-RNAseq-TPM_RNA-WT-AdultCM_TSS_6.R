

####################################################################################################################################################################################################################################
A0_Gene_TPM_Rep1      <- read.table("Other_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM.24413genes.TPM.txt",          header=FALSE,   sep="\t",   quote = "",   comment.char = "")  

A1_H2BGFP_H3_Rep1     <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/1_H3_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.1_H3_Rep1.wig.heatmap.xls",                        header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week0_Rep1  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/2_H2bGFP-week0_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.2_H2bGFP-week0_Rep1.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week0_Rep2  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/2_H2bGFP-week0_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.2_H2bGFP-week0_Rep2.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week1_Rep1  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/1_H2bGFP-week1_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.1_H2bGFP-week1_Rep1.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week1_Rep2  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/1_H2bGFP-week1_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.1_H2bGFP-week1_Rep2.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week2_Rep1  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/1_H2bGFP-week2_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.1_H2bGFP-week2_Rep1.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week2_Rep2  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/1_H2bGFP-week2_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.1_H2bGFP-week2_Rep2.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week4_Rep1  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/1_H2bGFP-week4_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.1_H2bGFP-week4_Rep1.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week4_Rep2  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/2_H2BGFP-CM-week4-2015_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.2_H2BGFP-CM-week4-2015_Rep1.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week4_Rep3  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/2_H2BGFP-CM-week4-2015_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.2_H2BGFP-CM-week4-2015_Rep2.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week6_Rep1  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/1_H2bGFP-week6_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.1_H2bGFP-week6_Rep1.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week6_Rep2  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/1_H2bGFP-week6_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.1_H2bGFP-week6_Rep2.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week8_Rep1  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/1_H2bGFP-week8_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.1_H2bGFP-week8_Rep1.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week0_EEDheto_Rep1  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/3_H2BGFP-CM-week0-EEDheto_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.3_H2BGFP-CM-week0-EEDheto_Rep1.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week0_EEDheto_Rep2  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/3_H2BGFP-CM-week0-EEDheto_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.3_H2BGFP-CM-week0-EEDheto_Rep2.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week0_EEDko_Rep1    <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/3_H2BGFP-CM-week0-EEDko_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.3_H2BGFP-CM-week0-EEDko_Rep1.wig.heatmap.xls",        header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week4_EEDheto_Rep1  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/4_H2BGFP-CM-week4-EEDheto_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.4_H2BGFP-CM-week4-EEDheto_Rep1.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week4_EEDheto_Rep2  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/4_H2BGFP-CM-week4-EEDheto_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.4_H2BGFP-CM-week4-EEDheto_Rep2.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week4_EEDko_Rep1    <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/4_H2BGFP-CM-week4-EEDko_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.4_H2BGFP-CM-week4-EEDko_Rep1.wig.heatmap.xls",        header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_week4_EEDko_Rep2    <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/4_H2BGFP-CM-week4-EEDko_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.4_H2BGFP-CM-week4-EEDko_Rep2.wig.heatmap.xls",        header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_banding_Rep1  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/5_H2BGFP-CM-week2.5-banding_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.5_H2BGFP-CM-week2.5-banding_Rep1.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_banding_Rep2  <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/5_H2BGFP-CM-week2.5-banding_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.5_H2BGFP-CM-week2.5-banding_Rep2.wig.heatmap.xls",    header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_sham_Rep1     <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/5_H2BGFP-CM-week2.5-sham_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.5_H2BGFP-CM-week2.5-sham_Rep1.wig.heatmap.xls",          header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H2BGFP_sham_Rep2     <- read.table("H2BGFP_Matrix/29-RNAseq-TPM/RNA-WT-AdultCM/5_H2BGFP-CM-week2.5-sham_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.5_H2BGFP-CM-week2.5-sham_Rep2.wig.heatmap.xls",          header=FALSE,   sep="\t",   quote = "",   comment.char = "")  

A2_360R1_Fog2_adult_Rep1  <- read.table("/media/yp/helab-1/BillPuLab_360R1_TFs_SE50/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/Fog2-adult_TSS_heatmap/RNA-WT-AdultCM.tss.Fog2-adult.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A2_360R1_Fog2_E12_Rep1    <- read.table("/media/yp/helab-1/BillPuLab_360R1_TFs_SE50/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/Fog2-E12.5_TSS_heatmap/RNA-WT-AdultCM.tss.Fog2-E12.5.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A2_360R1_Mef2c_adult_Rep1 <- read.table("/media/yp/helab-1/BillPuLab_360R1_TFs_SE50/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/Mef2c-adult_TSS_heatmap/RNA-WT-AdultCM.tss.Mef2c-adult.wig.heatmap.xls",         header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A2_360R1_Mef2c_E12_Rep1   <- read.table("/media/yp/helab-1/BillPuLab_360R1_TFs_SE50/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/Mef2c-E12.5_TSS_heatmap/RNA-WT-AdultCM.tss.Mef2c-E12.5.wig.heatmap.xls",         header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A2_360R1_Srf_adult_Rep1   <- read.table("/media/yp/helab-1/BillPuLab_360R1_TFs_SE50/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/Srf-adult_TSS_heatmap/RNA-WT-AdultCM.tss.Srf-adult.wig.heatmap.xls",             header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A2_360R1_Tbx5_adult_Rep1  <- read.table("/media/yp/helab-1/BillPuLab_360R1_TFs_SE50/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/Tbx5-adult_TSS_heatmap/RNA-WT-AdultCM.tss.Tbx5-adult.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A2_360R1_Tbx5_E12_Rep1    <- read.table("/media/yp/helab-1/BillPuLab_360R1_TFs_SE50/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/Tbx5-E12.5_TSS_heatmap/RNA-WT-AdultCM.tss.Tbx5-E12.5.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  

A3_MNase_EEDheto_Rep1   <- read.table("/media/yp/RAID0/4-H2BGFP/2-MNase/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/EEDhetoM-lib1610_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.EEDhetoM-lib1610_Rep1.wig.heatmap.xls",     header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A3_MNase_EEDheto_Rep2   <- read.table("/media/yp/RAID0/4-H2BGFP/2-MNase/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/EEDhetoF-lib1612_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.EEDhetoF-lib1612_Rep2.wig.heatmap.xls",     header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A3_MNase_EEDheto_merge  <- read.table("/media/yp/RAID0/4-H2BGFP/2-MNase/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/EEDhetoMF-merge_TSS_heatmap/RNA-WT-AdultCM.tss.EEDhetoMF-merge.wig.heatmap.xls",                 header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A3_MNase_EEDko_Rep1     <- read.table("/media/yp/RAID0/4-H2BGFP/2-MNase/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/EEDkoM-lib1611_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.EEDkoM-lib1611_Rep1.wig.heatmap.xls",         header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A3_MNase_EEDko_Rep2     <- read.table("/media/yp/RAID0/4-H2BGFP/2-MNase/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/EEDkoF-lib1613_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.EEDkoF-lib1613_Rep2.wig.heatmap.xls",         header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A3_MNase_EEDko_merge    <- read.table("/media/yp/RAID0/4-H2BGFP/2-MNase/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/EEDkoMF-merge_TSS_heatmap/RNA-WT-AdultCM.tss.EEDkoMF-merge.wig.heatmap.xls",                     header=FALSE,   sep="\t",   quote = "",   comment.char = "")  

A4_EED_Adult_Rep1 <- read.table("/media/yp/RAID0/4-H2BGFP/4-EEDadult2/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/EED-lib1663-adult_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.EED-lib1663-adult_Rep1.wig.heatmap.xls",             header=FALSE,   sep="\t",   quote = "",   comment.char = "")                                          
A4_EED_Adult_Rep2 <- read.table("/media/yp/RAID0/4-H2BGFP/4-EEDadult2/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/EED-lib1664-adult_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.EED-lib1664-adult_Rep2.wig.heatmap.xls",             header=FALSE,   sep="\t",   quote = "",   comment.char = "")                                          
A4_EED_P5_Rep1    <- read.table("/media/yp/RAID0/5-EED-Rescue/4_EEDP5/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/EEDP5-lib1601_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.EEDP5-lib1601_Rep1.wig.heatmap.xls",                     header=FALSE,   sep="\t",   quote = "",   comment.char = "")                                          
A4_EED_P5_Rep2    <- read.table("/media/yp/RAID0/5-EED-Rescue/4_EEDP5/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/EEDP5-lib1602_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.EEDP5-lib1602_Rep2.wig.heatmap.xls",                     header=FALSE,   sep="\t",   quote = "",   comment.char = "")                                          

A5_HDAC_HetoHDAC1   <- read.table("/media/yp/RAID0/5-EED-Rescue/5-HDAC/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/heto-lib1643-HDAC1_TSS_heatmap/RNA-WT-AdultCM.tss.1-BED_1_heto-lib1643-HDAC1-adultCMs_Rep1.bgsub.Fnor.smooth.wig.heatmap.xls",                 header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A5_HDAC_HetoHDAC2   <- read.table("/media/yp/RAID0/5-EED-Rescue/5-HDAC/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/heto-lib1644-HDAC2_TSS_heatmap/RNA-WT-AdultCM.tss.1-BED_1_heto-lib1644-HDAC2-adultCMs_Rep1.bgsub.Fnor.smooth.wig.heatmap.xls",                 header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A5_HDAC_HomoHDAC1   <- read.table("/media/yp/RAID0/5-EED-Rescue/5-HDAC/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/homo-lib1641-HDAC1_TSS_heatmap/RNA-WT-AdultCM.tss.1-BED_1_homo-lib1641-HDAC1-adultCMs_Rep1.bgsub.Fnor.smooth.wig.heatmap.xls",                 header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A5_HDAC_HomoHDAC2   <- read.table("/media/yp/RAID0/5-EED-Rescue/5-HDAC/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/homo-lib1642-HDAC2_TSS_heatmap/RNA-WT-AdultCM.tss.1-BED_1_homo-lib1642-HDAC2-adultCMs_Rep1.bgsub.Fnor.smooth.wig.heatmap.xls",                 header=FALSE,   sep="\t",   quote = "",   comment.char = "")  

A6_H3K27ac_EEDko_rep1 <- read.table("/media/yp/RAID0/5-EED-Rescue/1-H3K27ac/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/1_H3K27ac-EEDko-AdultCM_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.1_H3K27ac-EEDko-AdultCM_Rep1.wig.heatmap.xls",             header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A6_H3K27ac_EEDko_rep2 <- read.table("/media/yp/RAID0/5-EED-Rescue/1-H3K27ac/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/1_H3K27ac-EEDko-AdultCM_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.1_H3K27ac-EEDko-AdultCM_Rep2.wig.heatmap.xls",             header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A6_H3K27ac_WT_rep1    <- read.table("/media/yp/RAID0/5-EED-Rescue/1-H3K27ac/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/1_H3K27ac-WT-AdultCM_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.1_H3K27ac-WT-AdultCM_Rep1.wig.heatmap.xls",                   header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A6_H3K27ac_WT_rep2    <- read.table("/media/yp/RAID0/5-EED-Rescue/1-H3K27ac/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/1_H3K27ac-WT-AdultCM_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.1_H3K27ac-WT-AdultCM_Rep2.wig.heatmap.xls",                   header=FALSE,   sep="\t",   quote = "",   comment.char = "") 

A7_H3K27me3_EEDko_rep1     <- read.table("/media/yp/RAID0/5-EED-Rescue/2-H3K27me3/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/1_H3K27me3-EEDko-AdultCM_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.1_H3K27me3-EEDko-AdultCM_Rep1.wig.heatmap.xls",     header=FALSE,   sep="\t",   quote = "",   comment.char = "")
A7_H3K27me3_EEDko_rep2     <- read.table("/media/yp/RAID0/5-EED-Rescue/2-H3K27me3/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/1_H3K27me3-EEDko-AdultCM_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.1_H3K27me3-EEDko-AdultCM_Rep2.wig.heatmap.xls",     header=FALSE,   sep="\t",   quote = "",   comment.char = "") 
A7_H3K27me3_WT_rep1        <- read.table("/media/yp/RAID0/5-EED-Rescue/2-H3K27me3/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/1_H3K27me3-WT-AdultCM_Rep1_TSS_heatmap/RNA-WT-AdultCM.tss.1_H3K27me3-WT-AdultCM_Rep1.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A7_H3K27me3_WT_rep2        <- read.table("/media/yp/RAID0/5-EED-Rescue/2-H3K27me3/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/1_H3K27me3-WT-AdultCM_Rep2_TSS_heatmap/RNA-WT-AdultCM.tss.1_H3K27me3-WT-AdultCM_Rep2.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "") 
A7_H3K27me3_WT_rep3        <- read.table("/media/yp/RAID0/5-EED-Rescue/2-H3K27me3/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/1_H3K27me3-WT-AdultCM_Rep3_TSS_heatmap/RNA-WT-AdultCM.tss.1_H3K27me3-WT-AdultCM_Rep3.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  

A8_H3K4me1_Adult_Rep1    <- read.table("/media/yp/YongPeng8/2014NC-DNAme/SRP033385-HistoneMod/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/H3K4me1-CM-Adult_TSS_heatmap/RNA-WT-AdultCM.tss.H3K4me1-CM-Adult.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A8_H3K4me3_Adult_Rep1    <- read.table("/media/yp/YongPeng8/2014NC-DNAme/SRP033385-HistoneMod/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/H3K4me3-CM-Adult_TSS_heatmap/RNA-WT-AdultCM.tss.H3K4me3-CM-Adult.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A8_H3K27ac_Adult_Rep1    <- read.table("/media/yp/YongPeng8/2014NC-DNAme/SRP033385-HistoneMod/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/H3K27ac-CM-Adult_TSS_heatmap/RNA-WT-AdultCM.tss.H3K27ac-CM-Adult.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "") 
A8_H3K27me3_Adult_Rep1   <- read.table("/media/yp/YongPeng8/2014NC-DNAme/SRP033385-HistoneMod/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/H3K27me3-CM-Adult_TSS_heatmap/RNA-WT-AdultCM.tss.H3K27me3-CM-Adult.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A8_MeCP2_Adult_Rep1      <- read.table("/media/yp/YongPeng8/2014NC-DNAme/SRP033385-HistoneMod/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/MeCP2-adult_TSS_heatmap/RNA-WT-AdultCM.tss.MeCP2-adult.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A8_MeCP2_TAC_Rep1        <- read.table("/media/yp/YongPeng8/2014NC-DNAme/SRP033385-HistoneMod/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/MeCP2-TAC_TSS_heatmap/RNA-WT-AdultCM.tss.MeCP2-TAC.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  

A9_GATA4_Ab_Adult_Rep1        <- read.table("/media/yp/YongPeng8/AHeLab/GSE52123-2014NC/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/1_GATA4-Ab-Adult_TSS_heatmap/RNA-WT-AdultCM.tss.1_GATA4-Ab-Adult.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A9_GATA4_fb_Adult_Rep1        <- read.table("/media/yp/YongPeng8/AHeLab/GSE52123-2014NC/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/1_GATA4-fb-Adult_TSS_heatmap/RNA-WT-AdultCM.tss.1_GATA4-fb-Adult.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A9_GATA4_fb_Banding_Rep1      <- read.table("/media/yp/YongPeng8/AHeLab/GSE52123-2014NC/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/2_GATA4-fb-Banding_TSS_heatmap/RNA-WT-AdultCM.tss.2_GATA4-fb-Banding.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A9_GATA4_fb_Sham_Rep1         <- read.table("/media/yp/YongPeng8/AHeLab/GSE52123-2014NC/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/2_GATA4-fb-Sham_TSS_heatmap/RNA-WT-AdultCM.tss.2_GATA4-fb-Sham.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A9_H3K4me1_Adult_Rep1         <- read.table("/media/yp/YongPeng8/AHeLab/GSE52123-2014NC/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/3_H3K4me1-WT-Adult_TSS_heatmap/RNA-WT-AdultCM.tss.3_H3K4me1-WT-Adult.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A9_H3K4me3_Adult_Rep1         <- read.table("/media/yp/YongPeng8/AHeLab/GSE52123-2014NC/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/3_H3K4me3-WT-Adult_TSS_heatmap/RNA-WT-AdultCM.tss.3_H3K4me3-WT-Adult.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A9_H3K27ac_Adult_Rep1         <- read.table("/media/yp/YongPeng8/AHeLab/GSE52123-2014NC/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/3_H3K27ac-WT-Adult_TSS_heatmap/RNA-WT-AdultCM.tss.3_H3K27ac-WT-Adult.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A9_Pol2_Adult_Rep1            <- read.table("/media/yp/YongPeng8/AHeLab/GSE52123-2014NC/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/3_Pol2-WT-Adult_TSS_heatmap/RNA-WT-AdultCM.tss.3_Pol2-WT-Adult.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "") 
A9_H3K27me3_Adult_Rep1        <- read.table("/media/yp/YongPeng8/AHeLab/GSE52123-2014NC/DANPOS2/4-Profile/29-RNAseq-TPM/RNA-WT-AdultCM/4_H3K27me3-WT-Adult_TSS_heatmap/RNA-WT-AdultCM.tss.4_H3K27me3-WT-Adult.wig.heatmap.xls",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  




dim(A0_Gene_TPM_Rep1)

dim(A1_H2BGFP_H3_Rep1)
dim(A1_H2BGFP_week0_Rep1)
dim(A1_H2BGFP_week0_Rep2) 
dim(A1_H2BGFP_week1_Rep1) 
dim(A1_H2BGFP_week1_Rep2) 
dim(A1_H2BGFP_week2_Rep1) 
dim(A1_H2BGFP_week2_Rep2)
dim(A1_H2BGFP_week4_Rep1)
dim(A1_H2BGFP_week4_Rep2) 
dim(A1_H2BGFP_week4_Rep3) 
dim(A1_H2BGFP_week6_Rep1)
dim(A1_H2BGFP_week6_Rep2) 
dim(A1_H2BGFP_week8_Rep1) 
dim(A1_H2BGFP_week0_EEDheto_Rep1)  
dim(A1_H2BGFP_week0_EEDheto_Rep2)  
dim(A1_H2BGFP_week0_EEDko_Rep1)  
dim(A1_H2BGFP_week4_EEDheto_Rep1)  
dim(A1_H2BGFP_week4_EEDheto_Rep2)  
dim(A1_H2BGFP_week4_EEDko_Rep1)  
dim(A1_H2BGFP_week4_EEDko_Rep2)  
dim(A1_H2BGFP_banding_Rep1)  
dim(A1_H2BGFP_banding_Rep2)  
dim(A1_H2BGFP_sham_Rep1)  
dim(A1_H2BGFP_sham_Rep2)  

dim(A2_360R1_Fog2_adult_Rep1)  
dim(A2_360R1_Fog2_E12_Rep1)               
dim(A2_360R1_Mef2c_adult_Rep1)  
dim(A2_360R1_Mef2c_E12_Rep1)  
dim(A2_360R1_Srf_adult_Rep1)   
dim(A2_360R1_Tbx5_adult_Rep1)  
dim(A2_360R1_Tbx5_E12_Rep1)   

dim(A3_MNase_EEDheto_Rep1)   
dim(A3_MNase_EEDheto_Rep2)    
dim(A3_MNase_EEDheto_merge)  
dim(A3_MNase_EEDko_Rep1)     
dim(A3_MNase_EEDko_Rep2)     
dim(A3_MNase_EEDko_merge)    

dim(A4_EED_Adult_Rep1)                                         
dim(A4_EED_Adult_Rep2)                                   
dim(A4_EED_P5_Rep1)                                
dim(A4_EED_P5_Rep2)                                          

dim(A5_HDAC_HetoHDAC1)   
dim(A5_HDAC_HetoHDAC2)   
dim(A5_HDAC_HomoHDAC1)    
dim(A5_HDAC_HomoHDAC2)   

dim(A6_H3K27ac_EEDko_rep1) 
dim(A6_H3K27ac_EEDko_rep2) 
dim(A6_H3K27ac_WT_rep1)     
dim(A6_H3K27ac_WT_rep2)  

dim(A7_H3K27me3_EEDko_rep1)
dim(A7_H3K27me3_EEDko_rep2)
dim(A7_H3K27me3_WT_rep1)
dim(A7_H3K27me3_WT_rep2)
dim(A7_H3K27me3_WT_rep3)


dim(A8_H3K4me1_Adult_Rep1)
dim(A8_H3K4me3_Adult_Rep1)
dim(A8_H3K27ac_Adult_Rep1)
dim(A8_H3K27me3_Adult_Rep1)
dim(A8_MeCP2_Adult_Rep1)
dim(A8_MeCP2_TAC_Rep1)

dim(A9_GATA4_Ab_Adult_Rep1)
dim(A9_GATA4_fb_Adult_Rep1)
dim(A9_GATA4_fb_Banding_Rep1)
dim(A9_GATA4_fb_Sham_Rep1)
dim(A9_H3K4me1_Adult_Rep1)
dim(A9_H3K4me3_Adult_Rep1)
dim(A9_H3K27ac_Adult_Rep1)
dim(A9_Pol2_Adult_Rep1)
dim(A9_H3K27me3_Adult_Rep1)




A1_H2BGFP_H3_Rep1[c(9000:9200), 1]
A0_Gene_TPM_Rep1[c(9000:9200), 1]   
A0_Gene_TPM_Rep1[1,]            





identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_H3_Rep1[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week0_Rep1[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week0_Rep2[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week1_Rep1[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week1_Rep2[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week2_Rep1[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week2_Rep2[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week4_Rep1[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week4_Rep2[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week4_Rep3[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week6_Rep1[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week6_Rep2[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week8_Rep1[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week0_EEDheto_Rep1[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week0_EEDheto_Rep2[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week0_EEDko_Rep1[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week4_EEDheto_Rep1[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week4_EEDheto_Rep2[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week4_EEDko_Rep1[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_week4_EEDko_Rep2[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_banding_Rep1[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_banding_Rep2[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_sham_Rep1[, 1] )
identical(A1_H2BGFP_H3_Rep1[, 1],  A1_H2BGFP_sham_Rep2[, 1] )

identical(A1_H2BGFP_H3_Rep1[, 1],  A2_360R1_Fog2_adult_Rep1[, 1])  
identical(A1_H2BGFP_H3_Rep1[, 1],  A2_360R1_Fog2_E12_Rep1[, 1])               
identical(A1_H2BGFP_H3_Rep1[, 1],  A2_360R1_Mef2c_adult_Rep1[, 1])  
identical(A1_H2BGFP_H3_Rep1[, 1],  A2_360R1_Mef2c_E12_Rep1[, 1])  
identical(A1_H2BGFP_H3_Rep1[, 1],  A2_360R1_Srf_adult_Rep1[, 1])   
identical(A1_H2BGFP_H3_Rep1[, 1],  A2_360R1_Tbx5_adult_Rep1[, 1])  
identical(A1_H2BGFP_H3_Rep1[, 1],  A2_360R1_Tbx5_E12_Rep1[, 1])   

identical(A1_H2BGFP_H3_Rep1[, 1],  A3_MNase_EEDheto_Rep1[, 1])   
identical(A1_H2BGFP_H3_Rep1[, 1],  A3_MNase_EEDheto_Rep2[, 1])    
identical(A1_H2BGFP_H3_Rep1[, 1],  A3_MNase_EEDheto_merge[, 1])  
identical(A1_H2BGFP_H3_Rep1[, 1],  A3_MNase_EEDko_Rep1[, 1])     
identical(A1_H2BGFP_H3_Rep1[, 1],  A3_MNase_EEDko_Rep2[, 1])     
identical(A1_H2BGFP_H3_Rep1[, 1],  A3_MNase_EEDko_merge[, 1])    

identical(A1_H2BGFP_H3_Rep1[, 1],  A4_EED_Adult_Rep1[, 1])                                         
identical(A1_H2BGFP_H3_Rep1[, 1],  A4_EED_Adult_Rep2[, 1])                                   
identical(A1_H2BGFP_H3_Rep1[, 1],  A4_EED_P5_Rep1[, 1])                                
identical(A1_H2BGFP_H3_Rep1[, 1],  A4_EED_P5_Rep2[, 1])                                          

identical(A1_H2BGFP_H3_Rep1[, 1],  A5_HDAC_HetoHDAC1[, 1])   
identical(A1_H2BGFP_H3_Rep1[, 1],  A5_HDAC_HetoHDAC2[, 1])   
identical(A1_H2BGFP_H3_Rep1[, 1],  A5_HDAC_HomoHDAC1[, 1])    
identical(A1_H2BGFP_H3_Rep1[, 1],  A5_HDAC_HomoHDAC2[, 1])   

identical(A1_H2BGFP_H3_Rep1[, 1],  A6_H3K27ac_EEDko_rep1[, 1]) 
identical(A1_H2BGFP_H3_Rep1[, 1],  A6_H3K27ac_EEDko_rep2[, 1]) 
identical(A1_H2BGFP_H3_Rep1[, 1],  A6_H3K27ac_WT_rep1[, 1])     
identical(A1_H2BGFP_H3_Rep1[, 1],  A6_H3K27ac_WT_rep2[, 1])  

identical(A1_H2BGFP_H3_Rep1[, 1],  A7_H3K27me3_EEDko_rep1[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A7_H3K27me3_EEDko_rep2[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A7_H3K27me3_WT_rep1[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A7_H3K27me3_WT_rep2[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A7_H3K27me3_WT_rep3[, 1])

identical(A1_H2BGFP_H3_Rep1[, 1],  A8_H3K4me1_Adult_Rep1[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A8_H3K4me3_Adult_Rep1[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A8_H3K27ac_Adult_Rep1[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A8_H3K27me3_Adult_Rep1[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A8_MeCP2_Adult_Rep1[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A8_MeCP2_TAC_Rep1[, 1])

identical(A1_H2BGFP_H3_Rep1[, 1],  A9_GATA4_Ab_Adult_Rep1[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A9_GATA4_fb_Adult_Rep1[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A9_GATA4_fb_Banding_Rep1[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A9_GATA4_fb_Sham_Rep1[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A9_H3K4me1_Adult_Rep1[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A9_H3K4me3_Adult_Rep1[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A9_H3K27ac_Adult_Rep1[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A9_Pol2_Adult_Rep1[, 1])
identical(A1_H2BGFP_H3_Rep1[, 1],  A9_H3K27me3_Adult_Rep1[, 1])




dim(A1_H2BGFP_H3_Rep1)
myRowNames <- A1_H2BGFP_H3_Rep1[-1, 1]   
length(myRowNames)
myRowNames[1:10]

####################################################################################################################################################################################################################################













####################################################################################################################################################################################################################################
B0_Gene_TPM_Rep1     <- A0_Gene_TPM_Rep1[-1, -c(1,2)]

B1_H2BGFP_H3_Rep1    <- as.matrix(A1_H2BGFP_H3_Rep1[-1, -c(1:4)]) 
B1_H2BGFP_week0_Rep1 <- as.matrix(A1_H2BGFP_week0_Rep1[-1, -c(1:4)]) * 1.08
B1_H2BGFP_week0_Rep2 <- as.matrix(A1_H2BGFP_week0_Rep2[-1, -c(1:4)]) * 1.08
B1_H2BGFP_week1_Rep1 <- as.matrix(A1_H2BGFP_week1_Rep1[-1, -c(1:4)]) 
B1_H2BGFP_week1_Rep2 <- as.matrix(A1_H2BGFP_week1_Rep2[-1, -c(1:4)]) 
B1_H2BGFP_week2_Rep1 <- as.matrix(A1_H2BGFP_week2_Rep1[-1, -c(1:4)]) 
B1_H2BGFP_week2_Rep2 <- as.matrix(A1_H2BGFP_week2_Rep2[-1, -c(1:4)]) 
B1_H2BGFP_week4_Rep1 <- as.matrix(A1_H2BGFP_week4_Rep1[-1, -c(1:4)]) 
B1_H2BGFP_week4_Rep2 <- as.matrix(A1_H2BGFP_week4_Rep2[-1, -c(1:4)])
B1_H2BGFP_week4_Rep3 <- as.matrix(A1_H2BGFP_week4_Rep3[-1, -c(1:4)]) 
B1_H2BGFP_week6_Rep1 <- as.matrix(A1_H2BGFP_week6_Rep1[-1, -c(1:4)]) 
B1_H2BGFP_week6_Rep2 <- as.matrix(A1_H2BGFP_week6_Rep2[-1, -c(1:4)]) 
B1_H2BGFP_week8_Rep1 <- as.matrix(A1_H2BGFP_week8_Rep1[-1, -c(1:4)]) 
B1_H2BGFP_week0_EEDheto_Rep1 <- as.matrix(A1_H2BGFP_week0_EEDheto_Rep1[-1, -c(1:4)]) 
B1_H2BGFP_week0_EEDheto_Rep2 <- as.matrix(A1_H2BGFP_week0_EEDheto_Rep2[-1, -c(1:4)]) 
B1_H2BGFP_week0_EEDko_Rep1   <- as.matrix(A1_H2BGFP_week0_EEDko_Rep1[-1, -c(1:4)]) 
B1_H2BGFP_week4_EEDheto_Rep1 <- as.matrix(A1_H2BGFP_week4_EEDheto_Rep1[-1, -c(1:4)]) 
B1_H2BGFP_week4_EEDheto_Rep2 <- as.matrix(A1_H2BGFP_week4_EEDheto_Rep2[-1, -c(1:4)]) 
B1_H2BGFP_week4_EEDko_Rep1   <- as.matrix(A1_H2BGFP_week4_EEDko_Rep1[-1, -c(1:4)]) 
B1_H2BGFP_week4_EEDko_Rep2   <- as.matrix(A1_H2BGFP_week4_EEDko_Rep2[-1, -c(1:4)]) 
B1_H2BGFP_banding_Rep1 <- as.matrix(A1_H2BGFP_banding_Rep1[-1, -c(1:4)]) 
B1_H2BGFP_banding_Rep2 <- as.matrix(A1_H2BGFP_banding_Rep2[-1, -c(1:4)]) 
B1_H2BGFP_sham_Rep1    <- as.matrix(A1_H2BGFP_sham_Rep1[-1, -c(1:4)]) 
B1_H2BGFP_sham_Rep2    <- as.matrix(A1_H2BGFP_sham_Rep2[-1, -c(1:4)]) 

B2_360R1_Fog2_adult_Rep1   <- as.matrix(A2_360R1_Fog2_adult_Rep1[-1, -c(1:4)])  
B2_360R1_Fog2_E12_Rep1     <- as.matrix(A2_360R1_Fog2_E12_Rep1[-1, -c(1:4)])               
B2_360R1_Mef2c_adult_Rep1  <- as.matrix(A2_360R1_Mef2c_adult_Rep1[-1, -c(1:4)])  
B2_360R1_Mef2c_E12_Rep1    <- as.matrix(A2_360R1_Mef2c_E12_Rep1[-1, -c(1:4)])  
B2_360R1_Srf_adult_Rep1    <- as.matrix(A2_360R1_Srf_adult_Rep1[-1, -c(1:4)])   
B2_360R1_Tbx5_adult_Rep1   <- as.matrix(A2_360R1_Tbx5_adult_Rep1[-1, -c(1:4)])  
B2_360R1_Tbx5_E12_Rep1     <- as.matrix(A2_360R1_Tbx5_E12_Rep1[-1, -c(1:4)])   

B3_MNase_EEDheto_Rep1  <- as.matrix(A3_MNase_EEDheto_Rep1[-1, -c(1:4)])   
B3_MNase_EEDheto_Rep2  <- as.matrix(A3_MNase_EEDheto_Rep2[-1, -c(1:4)])    
B3_MNase_EEDheto_merge <- as.matrix(A3_MNase_EEDheto_merge[-1, -c(1:4)])  
B3_MNase_EEDko_Rep1    <- as.matrix(A3_MNase_EEDko_Rep1[-1, -c(1:4)])     
B3_MNase_EEDko_Rep2    <- as.matrix(A3_MNase_EEDko_Rep2[-1, -c(1:4)])     
B3_MNase_EEDko_merge   <- as.matrix(A3_MNase_EEDko_merge[-1, -c(1:4)])    

B4_EED_Adult_Rep1 <- as.matrix(A4_EED_Adult_Rep1[-1, -c(1:4)])                                         
B4_EED_Adult_Rep2 <- as.matrix(A4_EED_Adult_Rep2[-1, -c(1:4)])                                   
B4_EED_P5_Rep1    <- as.matrix(A4_EED_P5_Rep1[-1, -c(1:4)])                                
B4_EED_P5_Rep2    <- as.matrix(A4_EED_P5_Rep2[-1, -c(1:4)])                                          

B5_HDAC_HetoHDAC1 <- as.matrix(A5_HDAC_HetoHDAC1[-1, -c(1:4)])   
B5_HDAC_HetoHDAC2 <- as.matrix(A5_HDAC_HetoHDAC2[-1, -c(1:4)])   
B5_HDAC_HomoHDAC1 <- as.matrix(A5_HDAC_HomoHDAC1[-1, -c(1:4)])    
B5_HDAC_HomoHDAC2 <- as.matrix(A5_HDAC_HomoHDAC2[-1, -c(1:4)])   

B6_H3K27ac_EEDko_rep1 <- as.matrix(A6_H3K27ac_EEDko_rep1[-1, -c(1:4)]) 
B6_H3K27ac_EEDko_rep2 <- as.matrix(A6_H3K27ac_EEDko_rep2[-1, -c(1:4)]) 
B6_H3K27ac_WT_rep1    <- as.matrix(A6_H3K27ac_WT_rep1[-1, -c(1:4)])     
B6_H3K27ac_WT_rep2    <- as.matrix(A6_H3K27ac_WT_rep2[-1, -c(1:4)])  

B7_H3K27me3_EEDko_rep1 <- as.matrix(A7_H3K27me3_EEDko_rep1[-1, -c(1:4)]) 
B7_H3K27me3_EEDko_rep2 <- as.matrix(A7_H3K27me3_EEDko_rep2[-1, -c(1:4)])  
B7_H3K27me3_WT_rep1    <- as.matrix(A7_H3K27me3_WT_rep1[-1, -c(1:4)])     
B7_H3K27me3_WT_rep2    <- as.matrix(A7_H3K27me3_WT_rep2[-1, -c(1:4)])  
B7_H3K27me3_WT_rep3    <- as.matrix(A7_H3K27me3_WT_rep3[-1, -c(1:4)])  


B8_H3K4me1_Adult_Rep1  <- as.matrix(A8_H3K4me1_Adult_Rep1[-1, -c(1:4)])
B8_H3K4me3_Adult_Rep1  <- as.matrix(A8_H3K4me3_Adult_Rep1[-1, -c(1:4)])
B8_H3K27ac_Adult_Rep1  <- as.matrix(A8_H3K27ac_Adult_Rep1[-1, -c(1:4)])
B8_H3K27me3_Adult_Rep1 <- as.matrix(A8_H3K27me3_Adult_Rep1[-1, -c(1:4)])
B8_MeCP2_Adult_Rep1    <- as.matrix(A8_MeCP2_Adult_Rep1[-1, -c(1:4)])
B8_MeCP2_TAC_Rep1      <- as.matrix(A8_MeCP2_TAC_Rep1[-1, -c(1:4)])

B9_GATA4_Ab_Adult_Rep1   <- as.matrix(A9_GATA4_Ab_Adult_Rep1[-1, -c(1:4)])
B9_GATA4_fb_Adult_Rep1   <- as.matrix(A9_GATA4_fb_Adult_Rep1[-1, -c(1:4)])
B9_GATA4_fb_Banding_Rep1 <- as.matrix(A9_GATA4_fb_Banding_Rep1[-1, -c(1:4)])
B9_GATA4_fb_Sham_Rep1    <- as.matrix(A9_GATA4_fb_Sham_Rep1[-1, -c(1:4)])
B9_H3K4me1_Adult_Rep1    <- as.matrix(A9_H3K4me1_Adult_Rep1[-1, -c(1:4)])
B9_H3K4me3_Adult_Rep1    <- as.matrix(A9_H3K4me3_Adult_Rep1[-1, -c(1:4)])
B9_H3K27ac_Adult_Rep1    <- as.matrix(A9_H3K27ac_Adult_Rep1[-1, -c(1:4)])
B9_Pol2_Adult_Rep1       <- as.matrix(A9_Pol2_Adult_Rep1[-1, -c(1:4)])
B9_H3K27me3_Adult_Rep1   <- as.matrix(A9_H3K27me3_Adult_Rep1[-1, -c(1:4)])





dim(B0_Gene_TPM_Rep1)

dim(B1_H2BGFP_H3_Rep1)
dim(B1_H2BGFP_week0_Rep1)
dim(B1_H2BGFP_week0_Rep2) 
dim(B1_H2BGFP_week1_Rep1) 
dim(B1_H2BGFP_week1_Rep2) 
dim(B1_H2BGFP_week2_Rep1) 
dim(B1_H2BGFP_week2_Rep2)
dim(B1_H2BGFP_week4_Rep1)
dim(B1_H2BGFP_week4_Rep2) 
dim(B1_H2BGFP_week4_Rep3) 
dim(B1_H2BGFP_week6_Rep1)
dim(B1_H2BGFP_week6_Rep2) 
dim(B1_H2BGFP_week8_Rep1) 
dim(B1_H2BGFP_week0_EEDheto_Rep1)  
dim(B1_H2BGFP_week0_EEDheto_Rep2)  
dim(B1_H2BGFP_week0_EEDko_Rep1)  
dim(B1_H2BGFP_week4_EEDheto_Rep1)  
dim(B1_H2BGFP_week4_EEDheto_Rep2)  
dim(B1_H2BGFP_week4_EEDko_Rep1)  
dim(B1_H2BGFP_week4_EEDko_Rep2)  
dim(B1_H2BGFP_banding_Rep1)  
dim(B1_H2BGFP_banding_Rep2)  
dim(B1_H2BGFP_sham_Rep1)  
dim(B1_H2BGFP_sham_Rep2)  

dim(B2_360R1_Fog2_adult_Rep1)  
dim(B2_360R1_Fog2_E12_Rep1)               
dim(B2_360R1_Mef2c_adult_Rep1)  
dim(B2_360R1_Mef2c_E12_Rep1)  
dim(B2_360R1_Srf_adult_Rep1)   
dim(B2_360R1_Tbx5_adult_Rep1)  
dim(B2_360R1_Tbx5_E12_Rep1)   

dim(B3_MNase_EEDheto_Rep1)   
dim(B3_MNase_EEDheto_Rep2)    
dim(B3_MNase_EEDheto_merge)  
dim(B3_MNase_EEDko_Rep1)     
dim(B3_MNase_EEDko_Rep2)     
dim(B3_MNase_EEDko_merge)    

dim(B4_EED_Adult_Rep1)                                         
dim(B4_EED_Adult_Rep2)                                   
dim(B4_EED_P5_Rep1)                                
dim(B4_EED_P5_Rep2)                                          

dim(B5_HDAC_HetoHDAC1)   
dim(B5_HDAC_HetoHDAC2)   
dim(B5_HDAC_HomoHDAC1)    
dim(B5_HDAC_HomoHDAC2)   

dim(B6_H3K27ac_EEDko_rep1) 
dim(B6_H3K27ac_EEDko_rep2) 
dim(B6_H3K27ac_WT_rep1)     
dim(B6_H3K27ac_WT_rep2)  

dim(B7_H3K27me3_EEDko_rep1)
dim(B7_H3K27me3_EEDko_rep2)
dim(B7_H3K27me3_WT_rep1)
dim(B7_H3K27me3_WT_rep2)
dim(B7_H3K27me3_WT_rep3)


dim(B8_H3K4me1_Adult_Rep1)
dim(B8_H3K4me3_Adult_Rep1)
dim(B8_H3K27ac_Adult_Rep1)
dim(B8_H3K27me3_Adult_Rep1)
dim(B8_MeCP2_Adult_Rep1)
dim(B8_MeCP2_TAC_Rep1)

dim(B9_GATA4_Ab_Adult_Rep1)
dim(B9_GATA4_fb_Adult_Rep1)
dim(B9_GATA4_fb_Banding_Rep1)
dim(B9_GATA4_fb_Sham_Rep1)
dim(B9_H3K4me1_Adult_Rep1)
dim(B9_H3K4me3_Adult_Rep1)
dim(B9_H3K27ac_Adult_Rep1)
dim(B9_Pol2_Adult_Rep1)
dim(B9_H3K27me3_Adult_Rep1)



B7_H3K27me3_EEDko_rep1 <- (B7_H3K27me3_EEDko_rep1 - 0.4623429)/0.304912 
B7_H3K27me3_EEDko_rep2 <- (B7_H3K27me3_EEDko_rep2 - 0.5356368)/0.340226 
B7_H3K27me3_WT_rep1    <- (B7_H3K27me3_WT_rep1 - 0.2293108)/0.2612208     
B7_H3K27me3_WT_rep2    <- (B7_H3K27me3_WT_rep2 - 0.2072042)/0.2477243  
B7_H3K27me3_WT_rep3    <- (B7_H3K27me3_WT_rep3 - 0.2226703)/0.2570662  




noMoreThan <- function(matrix1) {
  matrix1A =  matrix1[,1:500] 
  matrix1B =  matrix1[,501:1000]
  print(ncol(matrix1))
  N1 <- nrow(matrix1)
  matrix2 <- matrix1A
  for (i  in c(1:N1) ) {
      means1A <- mean( as.numeric( matrix1A[i,] ) )
      means1B <- mean( as.numeric( matrix1B[i,] ) )
      if(means1A > means1B+0.01) { matrix2[i,] <- means1B*matrix1A[i,]/means1A}
  }
  return(matrix2) 
}

B7_H3K27me3_EEDko_rep1   <- noMoreThan( cbind(B7_H3K27me3_EEDko_rep1, B7_H3K27me3_WT_rep1) ) 
B7_H3K27me3_EEDko_rep2   <- noMoreThan( cbind(B7_H3K27me3_EEDko_rep2, B7_H3K27me3_WT_rep1) ) 

dim(B7_H3K27me3_EEDko_rep1) 
dim(B7_H3K27me3_EEDko_rep2) 



####################################################################################################################################################################################################################################















####################################################################################################################################################################################################################################
C0_Average_Gene_TPM      <- B0_Gene_TPM_Rep1

C1_Average_H3     <- B1_H2BGFP_H3_Rep1
C1_Average_week0  <- (B1_H2BGFP_week0_Rep1+B1_H2BGFP_week0_Rep2)/2
C1_Average_week1  <- (B1_H2BGFP_week1_Rep1+B1_H2BGFP_week1_Rep2)/2 
C1_Average_week2  <- (B1_H2BGFP_week2_Rep1+B1_H2BGFP_week2_Rep2)/2 
C1_Average_week4  <- (B1_H2BGFP_week4_Rep1+B1_H2BGFP_week4_Rep2+B1_H2BGFP_week4_Rep3)/3
C1_Average_week6  <- B1_H2BGFP_week6_Rep2
C1_Average_week8  <- B1_H2BGFP_week8_Rep1 
C1_Average_week0_EEDheto  <- (B1_H2BGFP_week0_EEDheto_Rep1+B1_H2BGFP_week0_EEDheto_Rep2)/2  
C1_Average_week0_EEDko    <- B1_H2BGFP_week0_EEDko_Rep1  
C1_Average_week4_EEDheto  <- B1_H2BGFP_week4_EEDheto_Rep2 
C1_Average_week4_EEDko    <- (B1_H2BGFP_week4_EEDko_Rep1+B1_H2BGFP_week4_EEDko_Rep2)/2  
C1_Average_banding  <- (B1_H2BGFP_banding_Rep1+B1_H2BGFP_banding_Rep2)/2  
C1_Average_sham     <- (B1_H2BGFP_sham_Rep1+B1_H2BGFP_sham_Rep2)/2 
 
C2_Average_360R1_Fog2_adult  <- B2_360R1_Fog2_adult_Rep1  
C2_Average_360R1_Fog2_E12    <- B2_360R1_Fog2_E12_Rep1               
C2_Average_360R1_Mef2c_adult <- B2_360R1_Mef2c_adult_Rep1  
C2_Average_360R1_Mef2c_E12   <- B2_360R1_Mef2c_E12_Rep1  
C2_Average_360R1_Srf_adult   <- B2_360R1_Srf_adult_Rep1   
C2_Average_360R1_Tbx5_adult  <- B2_360R1_Tbx5_adult_Rep1  
C2_Average_360R1_Tbx5_E12    <- B2_360R1_Tbx5_E12_Rep1   

C3_Average_MNase_EEDheto        <- B3_MNase_EEDheto_Rep2   
C3_Average_MNase_EEDheto_merge  <- B3_MNase_EEDheto_merge   
C3_Average_MNase_EEDko          <- B3_MNase_EEDko_Rep2     
C3_Average_MNase_EEDko_merge    <- B3_MNase_EEDko_merge      

C4_Average_EED_Adult <- (B4_EED_Adult_Rep1 + B4_EED_Adult_Rep2)/2                                      
C4_Average_EED_P5    <- (B4_EED_P5_Rep1 + B4_EED_P5_Rep2)/2                                         

C5_Average_HetoHDAC1 <- B5_HDAC_HetoHDAC1  
C5_Average_HetoHDAC2 <- B5_HDAC_HetoHDAC2  
C5_Average_HomoHDAC1 <- B5_HDAC_HomoHDAC1 
C5_Average_HomoHDAC2 <- B5_HDAC_HomoHDAC2 

C6_Average_H3K27ac_EEDko <- (B6_H3K27ac_EEDko_rep1 + B6_H3K27ac_EEDko_rep2)/2 
C6_Average_H3K27ac_WT    <- (B6_H3K27ac_WT_rep1 + B6_H3K27ac_WT_rep2)/2    

C7_Average_H3K27me3_EEDko <- (B7_H3K27me3_EEDko_rep1 + B7_H3K27me3_EEDko_rep2)/2 
C7_Average_H3K27me3_WT    <- (B7_H3K27me3_WT_rep1 + B7_H3K27me3_WT_rep2 + B7_H3K27me3_WT_rep3)/3    

C8_Average_H3K4me1_Adult  <- B8_H3K4me1_Adult_Rep1
C8_Average_H3K4me3_Adult  <- B8_H3K4me3_Adult_Rep1
C8_Average_H3K27ac_Adult  <- B8_H3K27ac_Adult_Rep1
C8_Average_H3K27me3_Adult <- B8_H3K27me3_Adult_Rep1
C8_Average_MeCP2_Adult    <- B8_MeCP2_Adult_Rep1
C8_Average_MeCP2_TAC      <- B8_MeCP2_TAC_Rep1

C9_Average_GATA4_Ab_Adult   <- B9_GATA4_Ab_Adult_Rep1
C9_Average_GATA4_fb_Adult   <- B9_GATA4_fb_Adult_Rep1
C9_Average_GATA4_fb_Banding <- B9_GATA4_fb_Banding_Rep1
C9_Average_GATA4_fb_Sham    <- B9_GATA4_fb_Sham_Rep1
C9_Average_H3K4me1_Adult    <- B9_H3K4me1_Adult_Rep1
C9_Average_H3K4me3_Adult    <- B9_H3K4me3_Adult_Rep1
C9_Average_H3K27ac_Adult    <- B9_H3K27ac_Adult_Rep1
C9_Average_Pol2_Adult       <- B9_Pol2_Adult_Rep1
C9_Average_H3K27me3_Adult   <- B9_H3K27me3_Adult_Rep1






dim(C0_Average_Gene_TPM)

dim(C1_Average_H3)
dim(C1_Average_week0)
dim(C1_Average_week1)
dim(C1_Average_week2)
dim(C1_Average_week4)
dim(C1_Average_week6)
dim(C1_Average_week8)
dim(C1_Average_week0_EEDheto)  
dim(C1_Average_week0_EEDko) 
dim(C1_Average_week4_EEDheto)
dim(C1_Average_week4_EEDko) 
dim(C1_Average_banding) 
dim(C1_Average_sham)

dim(C2_Average_360R1_Fog2_adult) 
dim(C2_Average_360R1_Fog2_E12)              
dim(C2_Average_360R1_Mef2c_adult)
dim(C2_Average_360R1_Mef2c_E12)  
dim(C2_Average_360R1_Srf_adult)   
dim(C2_Average_360R1_Tbx5_adult) 
dim(C2_Average_360R1_Tbx5_E12)  

dim(C3_Average_MNase_EEDheto)   
dim(C3_Average_MNase_EEDheto_merge)   
dim(C3_Average_MNase_EEDko)     
dim(C3_Average_MNase_EEDko_merge)    

dim(C4_Average_EED_Adult)                                   
dim(C4_Average_EED_P5)                                        

dim(C5_Average_HetoHDAC1)  
dim(C5_Average_HetoHDAC2) 
dim(C5_Average_HomoHDAC1)
dim(C5_Average_HomoHDAC2) 

dim(C6_Average_H3K27ac_EEDko)
dim(C6_Average_H3K27ac_WT)   

dim(C7_Average_H3K27me3_EEDko)
dim(C7_Average_H3K27me3_WT)   

dim(C8_Average_H3K4me1_Adult)
dim(C8_Average_H3K4me3_Adult)
dim(C8_Average_H3K27ac_Adult)
dim(C8_Average_H3K27me3_Adult)
dim(C8_Average_MeCP2_Adult)
dim(C8_Average_MeCP2_TAC)

dim(C9_Average_GATA4_Ab_Adult)
dim(C9_Average_GATA4_fb_Adult)
dim(C9_Average_GATA4_fb_Banding)
dim(C9_Average_GATA4_fb_Sham)
dim(C9_Average_H3K4me1_Adult)
dim(C9_Average_H3K4me3_Adult)
dim(C9_Average_H3K27ac_Adult)
dim(C9_Average_Pol2_Adult)
dim(C9_Average_H3K27me3_Adult)



####################################################################################################################################################################################################################################













####################################################################################################################################################################################################################################
D0_Average_Gene_TPM      <- C0_Average_Gene_TPM

D1_Average_H3 <- C1_Average_H3[, c(225:275)]     
D1_Average_week0 <- C1_Average_week0[, c(225:275)]  
D1_Average_week1 <- C1_Average_week1[, c(225:275)]   
D1_Average_week2 <- C1_Average_week2[, c(225:275)]  
D1_Average_week4 <- C1_Average_week4[, c(225:275)]  
D1_Average_week6 <- C1_Average_week6[, c(225:275)]  
D1_Average_week8 <- C1_Average_week8[, c(225:275)]  
D1_Average_week0_EEDheto <- C1_Average_week0_EEDheto[, c(225:275)]  
D1_Average_week0_EEDko <- C1_Average_week0_EEDko[, c(225:275)]   
D1_Average_week4_EEDheto <- C1_Average_week4_EEDheto[, c(225:275)]  
D1_Average_week4_EEDko <- C1_Average_week4_EEDko[, c(225:275)]   
D1_Average_banding <- C1_Average_banding[, c(225:275)]  
D1_Average_sham <- C1_Average_sham[, c(225:275)]     
	
D2_Average_360R1_Fog2_adult  <- C2_Average_360R1_Fog2_adult[, c(225:275)]  
D2_Average_360R1_Fog2_E12    <- C2_Average_360R1_Fog2_E12[, c(225:275)]               
D2_Average_360R1_Mef2c_adult <- C2_Average_360R1_Mef2c_adult[, c(225:275)] 
D2_Average_360R1_Mef2c_E12   <- C2_Average_360R1_Mef2c_E12[, c(225:275)]   
D2_Average_360R1_Srf_adult   <- C2_Average_360R1_Srf_adult[, c(225:275)]    
D2_Average_360R1_Tbx5_adult  <- C2_Average_360R1_Tbx5_adult[, c(225:275)]  
D2_Average_360R1_Tbx5_E12    <- C2_Average_360R1_Tbx5_E12[, c(225:275)]   

D3_Average_MNase_EEDheto       <- C3_Average_MNase_EEDheto[, c(225:275)]        
D3_Average_MNase_EEDheto_merge <- C3_Average_MNase_EEDheto_merge[, c(225:275)]     
D3_Average_MNase_EEDko         <- C3_Average_MNase_EEDko[, c(225:275)]                
D3_Average_MNase_EEDko_merge   <- C3_Average_MNase_EEDko_merge[, c(225:275)]        

D4_Average_EED_Adult <- C4_Average_EED_Adult[, c(225:275)]                                     
D4_Average_EED_P5    <- C4_Average_EED_P5[, c(225:275)]                                          

D5_Average_HetoHDAC1 <- C5_Average_HetoHDAC1[, c(225:275)] 
D5_Average_HetoHDAC2 <- C5_Average_HetoHDAC2[, c(225:275)]  
D5_Average_HomoHDAC1 <- C5_Average_HomoHDAC1[, c(225:275)] 
D5_Average_HomoHDAC2 <- C5_Average_HomoHDAC2[, c(225:275)]  

D6_Average_H3K27ac_EEDko <- C6_Average_H3K27ac_EEDko[, c(225:275)] 
D6_Average_H3K27ac_WT    <- C6_Average_H3K27ac_WT[, c(225:275)]    

D7_Average_H3K27me3_EEDko <- C7_Average_H3K27me3_EEDko[, c(225:275)] 
D7_Average_H3K27me3_WT    <- C7_Average_H3K27me3_WT[, c(225:275)]    

D8_Average_H3K4me1_Adult  <- C8_Average_H3K4me1_Adult[, c(225:275)]
D8_Average_H3K4me3_Adult  <- C8_Average_H3K4me3_Adult[, c(225:275)]
D8_Average_H3K27ac_Adult  <- C8_Average_H3K27ac_Adult[, c(225:275)]
D8_Average_H3K27me3_Adult <- C8_Average_H3K27me3_Adult[, c(225:275)]
D8_Average_MeCP2_Adult    <- C8_Average_MeCP2_Adult[, c(225:275)]
D8_Average_MeCP2_TAC      <- C8_Average_MeCP2_TAC[, c(225:275)]

D9_Average_GATA4_Ab_Adult   <- C9_Average_GATA4_Ab_Adult[, c(225:275)]
D9_Average_GATA4_fb_Adult   <- C9_Average_GATA4_fb_Adult[, c(225:275)]
D9_Average_GATA4_fb_Banding <- C9_Average_GATA4_fb_Banding[, c(225:275)]
D9_Average_GATA4_fb_Sham    <- C9_Average_GATA4_fb_Sham[, c(225:275)]
D9_Average_H3K4me1_Adult    <- C9_Average_H3K4me1_Adult[, c(225:275)]
D9_Average_H3K4me3_Adult    <- C9_Average_H3K4me3_Adult[, c(225:275)]
D9_Average_H3K27ac_Adult    <- C9_Average_H3K27ac_Adult[, c(225:275)]
D9_Average_Pol2_Adult       <- C9_Average_Pol2_Adult[, c(225:275)]
D9_Average_H3K27me3_Adult   <- C9_Average_H3K27me3_Adult[, c(225:275)]










dim(D0_Average_Gene_TPM)

dim(D1_Average_H3)
dim(D1_Average_week0)
dim(D1_Average_week1)
dim(D1_Average_week2)
dim(D1_Average_week4)
dim(D1_Average_week6)
dim(D1_Average_week8)
dim(D1_Average_week0_EEDheto)  
dim(D1_Average_week0_EEDko) 
dim(D1_Average_week4_EEDheto)
dim(D1_Average_week4_EEDko) 
dim(D1_Average_banding) 
dim(D1_Average_sham)

dim(D2_Average_360R1_Fog2_adult) 
dim(D2_Average_360R1_Fog2_E12)              
dim(D2_Average_360R1_Mef2c_adult)
dim(D2_Average_360R1_Mef2c_E12)  
dim(D2_Average_360R1_Srf_adult)   
dim(D2_Average_360R1_Tbx5_adult) 
dim(D2_Average_360R1_Tbx5_E12)  

dim(D3_Average_MNase_EEDheto)   
dim(D3_Average_MNase_EEDheto_merge)   
dim(D3_Average_MNase_EEDko)     
dim(D3_Average_MNase_EEDko_merge)    

dim(D4_Average_EED_Adult)                                   
dim(D4_Average_EED_P5)                                        

dim(D5_Average_HetoHDAC1)  
dim(D5_Average_HetoHDAC2) 
dim(D5_Average_HomoHDAC1)
dim(D5_Average_HomoHDAC2) 

dim(D6_Average_H3K27ac_EEDko)
dim(D6_Average_H3K27ac_WT)   

dim(D7_Average_H3K27me3_EEDko)
dim(D7_Average_H3K27me3_WT)   

dim(D8_Average_H3K4me1_Adult)
dim(D8_Average_H3K4me3_Adult)
dim(D8_Average_H3K27ac_Adult)
dim(D8_Average_H3K27me3_Adult)
dim(D8_Average_MeCP2_Adult)
dim(D8_Average_MeCP2_TAC)

dim(D9_Average_GATA4_Ab_Adult)
dim(D9_Average_GATA4_fb_Adult)
dim(D9_Average_GATA4_fb_Banding)
dim(D9_Average_GATA4_fb_Sham)
dim(D9_Average_H3K4me1_Adult)
dim(D9_Average_H3K4me3_Adult)
dim(D9_Average_H3K27ac_Adult)
dim(D9_Average_Pol2_Adult)
dim(D9_Average_H3K27me3_Adult)


####################################################################################################################################################################################################################################










####################################################################################################################################################################################################################################
howToSort_1A <- function(vector1) {
  valueToSort = mean( as.numeric( vector1 ) )
  print(valueToSort)
  return(valueToSort)
}

howToSort_1B <- function(vector1) {
  valueToSort = mean( as.numeric( vector1[4:5] ) )
  print(valueToSort)
  return(valueToSort)
}

howToSort_2A <- function(vector1) {
  valueToSort = mean( as.numeric( vector1[75:125]) )/mean( as.numeric( vector1[275:325]) )
  print(valueToSort)
  return(valueToSort)
}

howToSort_2B <- function(vector1) {
  valueToSort = mean( as.numeric( vector1[8] ) + 0.001)/mean( as.numeric( vector1[9]) + 0.001)
  print(valueToSort)
  return(valueToSort)
}

howToSort_2C <- function(vector1) {
  valueToSort = mean( as.numeric( vector1[1:201] ) + 0.00001)/mean( as.numeric( vector1[202:402]) + 0.00001)
  print(valueToSort)
  return(valueToSort)
}

howToSort_2D <- function(vector1) {
  valueToSort = mean( as.numeric( vector1[1:500] ) ) - mean( as.numeric( vector1[501:1000]))
  print(valueToSort)
  return(valueToSort)
}

howToSort_2E <- function(vector1) {
  valueToSort = mean( as.numeric( vector1[4:5] ) ) - mean( as.numeric( vector1[6:7]))
  print(valueToSort)
  return(valueToSort)
}

howToSort_2F <- function(vector1) {
  valueToSort = mean( as.numeric( vector1[7:8] ) + 0.00001 ) / mean( as.numeric( vector1[9:10])+ 0.00001 )
  print(valueToSort)
  return(valueToSort)
}

NFRlength1 <- function(vector1) {
     n1 = length(vector1)
     n2 = n1/2
     ## n2 = which.min(vector1[245:255]) + 244
     m = -1  ## the length of NFR (nucleosome free region)
     thres1 = -0.001
     for(i in c(n2:5) ) {
       if( (vector1[i-1]<vector1[i] + thres1) & (vector1[i-2]<vector1[i] + thres1) & (vector1[i-3]<vector1[i] + thres1)   & (vector1[i-4]<vector1[i] + thres1)   ) { m = m+(n2-i); break; }
       if(vector1[i] > 0.4) {m = m+(n2-i); break; }
     }
     for(i in c(n2:(n1-4)) ) {
       if( (vector1[i+1]<vector1[i] + thres1) & (vector1[i+2]<vector1[i] + thres1) & (vector1[i+3]<vector1[i] + thres1)  & (vector1[i+4]<vector1[i] + thres1) ) { m = m+(i-n2); break; }
       if(vector1[i] > 0.4) {m = m+(i-n2); break; }
     }
     return(m)
}

NFRlength2 <- function(vector1) {
  n1 = length(vector1)
  n2 = which.min(vector1[245:255]) + 245
  ## n2 = n1/2
  print(n2)
  m = 0  ## the length of NFR (nucleosome free region)
  for(i in c(n2:5) ) {
    if( (vector1[i-1] < vector1[i]   )  & (vector1[i-2] < vector1[i]   )   ) { m = m+(n2-i); break; }
  }
  for(i in c(n2:(n1-4)) ) {
    if( (vector1[i+1] < vector1[i]   )  & (vector1[i+2] < vector1[i]   ) ) { m = m+(i-n2); break; }
  }
  #mean1 <- mean(vector1[245:255])
  #m = m + mean1/100
  return(m)
}



AllResults_g <- "NTR_Heatmap/29-RNAseq-TPM/RNA-WT-AdultCM/TSS/sort6_Pol2"
if( ! file.exists(AllResults_g) ) { dir.create(path=AllResults_g,   recursive = TRUE) }

toSort1   <- apply( C9_Average_Pol2_Adult,  1, howToSort_1A )   
print( as.numeric(sort(toSort1) ) )
 
length(toSort1)
index1     <- order(toSort1)        
length(index1)
toSort1[ index1[1:100] ]


myRowNames2 <- myRowNames[index1]
length(myRowNames2)



####################################################################################################################################################################################################################################









####################################################################################################################################################################################################################################
dir1 <- paste(AllResults_g,  "/1_merged_500bp",  sep="")
if( ! file.exists(dir1) ) { dir.create(path=dir1,   recursive = TRUE) }

numOfRows <- nrow(A1_H2BGFP_H3_Rep1)
zeroCol   <- cbind( rep("NA", numOfRows ), rep("NA", numOfRows ),  rep("NA", numOfRows ), rep("NA", numOfRows ),  rep("NA", numOfRows ), 
                    rep("NA", numOfRows ), rep("NA", numOfRows ),  rep("NA", numOfRows ), rep("NA", numOfRows ),  rep("NA", numOfRows ) )                               

mode(myRowNames2)  
myRowNames3 <- c("name1", as.character(myRowNames2) ) 
myRowNames4 <- c("name2", as.character(myRowNames2) ) 
length(myRowNames3)
length(myRowNames4)

myFirstLine <- A1_H2BGFP_H3_Rep1[1, c(229:279)]
names(myFirstLine)
rownames(D1_Average_H3)
colnames(D1_Average_H3)


E0_Gene_TPM_out <- A0_Gene_TPM_Rep1[c(1, index1+1),]

E1_H3_WT_adult_sort      <-  	rbind( myFirstLine ,    D1_Average_H3[index1,    ]    )
E1_H2BGFP_WT_week0_sort  <-  	rbind( myFirstLine ,    D1_Average_week0[index1, ]    )
E1_H2BGFP_WT_week1_sort  <-   rbind( myFirstLine ,    D1_Average_week1[index1, ]    ) 
E1_H2BGFP_WT_week2_sort  <-  	rbind( myFirstLine ,    D1_Average_week2[index1, ]    )
E1_H2BGFP_WT_week4_sort  <-  	rbind( myFirstLine ,    D1_Average_week4[index1, ]    )
E1_H2BGFP_WT_week6_sort  <-  	rbind( myFirstLine ,    D1_Average_week6[index1, ]    )
E1_H2BGFP_WT_week8_sort  <-   rbind( myFirstLine ,    D1_Average_week8[index1, ]    ) 
E1_H2BGFP_WT_out         <-   cbind(myRowNames3, myRowNames4, 
                                    E1_H3_WT_adult_sort,      zeroCol,   E1_H2BGFP_WT_week0_sort,  zeroCol, 
                                    E1_H2BGFP_WT_week1_sort,  zeroCol,   E1_H2BGFP_WT_week2_sort,  zeroCol,  
                                    E1_H2BGFP_WT_week4_sort,  zeroCol,   E1_H2BGFP_WT_week6_sort,  zeroCol,   E1_H2BGFP_WT_week8_sort )


E1_week0_EEDheto_sort      <-  	rbind( myFirstLine ,    D1_Average_week0_EEDheto[index1,  ]    )
E1_week0_EEDko_sort        <-  	rbind( myFirstLine ,    D1_Average_week0_EEDko[index1,    ]    )
E1_week4_EEDheto_sort      <-  	rbind( myFirstLine ,    D1_Average_week4_EEDheto[index1,  ]    )
E1_week4_EEDko_sort        <-  	rbind( myFirstLine ,    D1_Average_week4_EEDko[index1,    ]    )
E1_H2BGFP_CKO_out          <-   cbind(myRowNames3, myRowNames4, E1_week0_EEDheto_sort,  zeroCol,   E1_week0_EEDko_sort,  zeroCol, E1_week4_EEDheto_sort,  zeroCol,   E1_week4_EEDko_sort )


E1_banding_sort     <-  	rbind( myFirstLine ,    D1_Average_banding[index1,  ]    )
E1_sham_sort        <-  	rbind( myFirstLine ,    D1_Average_sham[index1,     ]    )
E1_H2BGFP_TAC_out   <-    cbind(myRowNames3, myRowNames4,  E1_sham_sort,  zeroCol,   E1_banding_sort )




E2_360R1_Fog2_adult_sort     <-  rbind( myFirstLine ,    D2_Average_360R1_Fog2_adult[index1,  ]    )
E2_360R1_Fog2_E12_sort       <-  rbind( myFirstLine ,    D2_Average_360R1_Fog2_E12[index1,  ]    )
E2_360R1_Mef2c_adult_sort    <-  rbind( myFirstLine ,    D2_Average_360R1_Mef2c_adult[index1, ]    )
E2_360R1_Mef2c_E12_sort      <-  rbind( myFirstLine ,    D2_Average_360R1_Mef2c_E12[index1, ]    )
E2_360R1_Srf_adult_sort      <-  rbind( myFirstLine ,    D2_Average_360R1_Srf_adult[index1, ]    )
E2_360R1_Tbx5_adult_sort     <-  rbind( myFirstLine ,    D2_Average_360R1_Tbx5_adult[index1,] )
E2_360R1_Tbx5_E12_sort       <-  rbind( myFirstLine ,    D2_Average_360R1_Tbx5_E12[index1,] )
E2_360R1_adult_out           <-  cbind(myRowNames3, myRowNames4,  E2_360R1_Fog2_adult_sort ,  zeroCol,   E2_360R1_Mef2c_adult_sort,   zeroCol,    E2_360R1_Srf_adult_sort,  zeroCol, E2_360R1_Tbx5_adult_sort )
E2_360R1_E12_out             <-  cbind(myRowNames3, myRowNames4,  E2_360R1_Fog2_E12_sort ,    zeroCol,   E2_360R1_Mef2c_E12_sort,     zeroCol,    E2_360R1_Tbx5_E12_sort )




E3_Average_MNase_EEDheto_sort     <-  	rbind( myFirstLine ,    D3_Average_MNase_EEDheto[index1,  ]    )
E3_Average_MNase_EEDko_sort       <-  	rbind( myFirstLine ,    D3_Average_MNase_EEDko[index1,    ]    )
E3_Average_MNase_out              <-    cbind(myRowNames3, myRowNames4,  E3_Average_MNase_EEDheto_sort,          zeroCol,   E3_Average_MNase_EEDko_sort )
E3_Average_MNase_out2             <-    cbind(myRowNames3, myRowNames4,  E3_Average_MNase_EEDheto_sort[,10:40],   zeroCol,   E3_Average_MNase_EEDko_sort[,10:40] )

E3_merge_MNase_EEDheto_sort     <-  	rbind( myFirstLine ,    D3_Average_MNase_EEDheto_merge[index1,  ]    )
E3_merge_MNase_EEDko_sort       <-  	rbind( myFirstLine ,    D3_Average_MNase_EEDko_merge[index1,    ]    )
E3_merge_MNase_out              <-    cbind(myRowNames3, myRowNames4,  E3_merge_MNase_EEDheto_sort,          zeroCol,   E3_merge_MNase_EEDko_sort )
E3_merge_MNase_out2             <-    cbind(myRowNames3, myRowNames4,  E3_merge_MNase_EEDheto_sort[,10:40],   zeroCol,   E3_merge_MNase_EEDko_sort[,10:40] )



E4_EED_Adult_sort   <-  	rbind( myFirstLine ,    D4_Average_EED_Adult[index1,  ]    )
E4_EED_P5_sort      <-  	rbind( myFirstLine ,    D4_Average_EED_P5[index1,    ]    )
E4_EED_out          <-    cbind(myRowNames3, myRowNames4,  E4_EED_P5_sort,  zeroCol,   E4_EED_Adult_sort )


E5_HetoHDAC1_sort   <-  	rbind( myFirstLine ,    D5_Average_HetoHDAC1[index1,  ]    )
E5_HetoHDAC2_sort   <-  	rbind( myFirstLine ,    D5_Average_HetoHDAC2[index1,  ]    )
E5_HomoHDAC1_sort   <-  	rbind( myFirstLine ,    D5_Average_HomoHDAC1[index1,  ]    )
E5_HomoHDAC2_sort   <-  	rbind( myFirstLine ,    D5_Average_HomoHDAC2[index1,  ]    )
E5_HDAC_out         <-   cbind(myRowNames3, myRowNames4, E5_HetoHDAC1_sort,  zeroCol,   E5_HomoHDAC1_sort,  zeroCol, E5_HetoHDAC2_sort,  zeroCol,   E5_HomoHDAC2_sort )


E6_H3K27ac_WT_sort         <-  	rbind( myFirstLine ,    D6_Average_H3K27ac_WT[index1,  ]    )
E6_H3K27ac_EEDko_sort      <-  	rbind( myFirstLine ,    D6_Average_H3K27ac_EEDko[index1,    ]    )
E6_H3K27ac_out             <-   cbind(myRowNames3, myRowNames4,  E6_H3K27ac_WT_sort,  zeroCol,   E6_H3K27ac_EEDko_sort )


E7_H3K27me3_WT_sort         <-  rbind( myFirstLine ,    D7_Average_H3K27me3_WT[index1,  ]    )
E7_H3K27me3_EEDko_sort      <-  rbind( myFirstLine ,    D7_Average_H3K27me3_EEDko[index1,    ]    )
E7_H3K27me3_out             <-  cbind(myRowNames3, myRowNames4,  E7_H3K27me3_WT_sort,  zeroCol,   E7_H3K27me3_EEDko_sort )


E8_Average_H3K4me1_Adult    <-  rbind( myFirstLine ,    D8_Average_H3K4me1_Adult[index1,  ]    )
E8_Average_H3K4me3_Adult    <-  rbind( myFirstLine ,    D8_Average_H3K4me3_Adult[index1,  ]    )
E8_Average_H3K27ac_Adult    <-  rbind( myFirstLine ,    D8_Average_H3K27ac_Adult[index1,  ]    )
E8_Average_H3K27me3_Adult   <-  rbind( myFirstLine ,    D8_Average_H3K27me3_Adult[index1,  ]    )
E8_Average_MeCP2_Adult      <-  rbind( myFirstLine ,    D8_Average_MeCP2_Adult[index1,  ]    )
E8_Average_MeCP2_TAC        <-  rbind( myFirstLine ,    D8_Average_MeCP2_TAC[index1,  ]    )
E8_DNAme_2014NC_out         <-  cbind(myRowNames3, myRowNames4,  
                                      E8_Average_H3K4me1_Adult,  zeroCol,   E8_Average_H3K4me3_Adult, zeroCol,   E8_Average_H3K27ac_Adult,   zeroCol,  
                                      E8_Average_H3K27me3_Adult, zeroCol,   E8_Average_MeCP2_Adult,   zeroCol,   E8_Average_MeCP2_TAC )



E9_Average_GATA4_Ab_Adult    <-  rbind( myFirstLine ,    D9_Average_GATA4_Ab_Adult[index1,  ]    )
E9_Average_GATA4_fb_Adult    <-  rbind( myFirstLine ,    D9_Average_GATA4_fb_Adult[index1,  ]    )
E9_Average_GATA4_fb_Banding  <-  rbind( myFirstLine ,    D9_Average_GATA4_fb_Banding[index1,  ]    )
E9_Average_GATA4_fb_Sham     <-  rbind( myFirstLine ,    D9_Average_GATA4_fb_Sham[index1,  ]    )
E9_Average_H3K4me1_Adult     <-  rbind( myFirstLine ,    D9_Average_H3K4me1_Adult[index1,  ]    )
E9_Average_H3K4me3_Adult     <-  rbind( myFirstLine ,    D9_Average_H3K4me3_Adult[index1,  ]    )
E9_Average_H3K27ac_Adult     <-  rbind( myFirstLine ,    D9_Average_H3K27ac_Adult[index1,  ]    )
E9_Average_Pol2_Adult        <-  rbind( myFirstLine ,    D9_Average_Pol2_Adult[index1,  ]    )
E9_Average_H3K27me3_Adult    <-  rbind( myFirstLine ,    D9_Average_H3K27me3_Adult[index1,  ]    )
E9_AHe_2014NC_GATA4_out      <-  cbind(myRowNames3, myRowNames4, E9_Average_GATA4_Ab_Adult,  zeroCol,  E9_Average_GATA4_fb_Adult, zeroCol,  E9_Average_GATA4_fb_Banding, zeroCol,  E9_Average_GATA4_fb_Sham )
E9_AHe_2014NC_HisMod_out     <-  cbind(myRowNames3, myRowNames4, 
                                       E9_Average_H3K4me1_Adult,    zeroCol,  E9_Average_H3K4me3_Adult,  zeroCol,  
                                       E9_Average_H3K27ac_Adult,    zeroCol,  E9_Average_H3K27me3_Adult, zeroCol, E9_Average_Pol2_Adult)





dim(E0_Gene_TPM_out)
dim(E1_H2BGFP_WT_out)
dim(E1_H2BGFP_CKO_out)
dim(E1_H2BGFP_TAC_out)
dim(E2_360R1_adult_out)
dim(E2_360R1_E12_out)
dim(E3_Average_MNase_out)
dim(E3_merge_MNase_out)
dim(E3_Average_MNase_out2)
dim(E3_merge_MNase_out2)
dim(E4_EED_out )
dim(E5_HDAC_out)
dim(E6_H3K27ac_out)
dim(E7_H3K27me3_out)
dim(E8_DNAme_2014NC_out)
dim(E9_AHe_2014NC_GATA4_out)
dim(E9_AHe_2014NC_HisMod_out)



write.table( E0_Gene_TPM_out,        file=paste(dir1, "/E0_Gene_TPM_out.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E1_H2BGFP_WT_out,       file=paste(dir1, "/E1_H2BGFP_WT_out.txt",    sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E1_H2BGFP_CKO_out,      file=paste(dir1, "/E1_H2BGFP_CKO_out.txt",   sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E1_H2BGFP_TAC_out,      file=paste(dir1, "/E1_H2BGFP_TAC_out.txt",   sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E2_360R1_adult_out,     file=paste(dir1, "/E2_360R1_adult_out.txt",  sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E2_360R1_E12_out,       file=paste(dir1, "/E2_360R1_E12_out.txt",    sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E3_Average_MNase_out,   file=paste(dir1, "/E3_Average_MNase_out.txt",sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E3_merge_MNase_out,     file=paste(dir1, "/E3_merge_MNase_out.txt",  sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E3_Average_MNase_out2,  file=paste(dir1, "/E3_Average_MNase_out2.txt",sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E3_merge_MNase_out2,    file=paste(dir1, "/E3_merge_MNase_out2.txt",  sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E4_EED_out,             file=paste(dir1, "/E4_EED_out.txt",          sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E5_HDAC_out,            file=paste(dir1, "/E5_HDAC_out.txt",         sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E6_H3K27ac_out,         file=paste(dir1, "/E6_H3K27ac_out.txt",      sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E7_H3K27me3_out,        file=paste(dir1, "/E7_H3K27me3_out.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E8_DNAme_2014NC_out,    file=paste(dir1, "/E8_DNAme_2014NC_out.txt", sep=""),      quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E9_AHe_2014NC_GATA4_out,        file=paste(dir1, "/E9_AHe_2014NC_GATA4_out.txt",sep=""),      quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( E9_AHe_2014NC_HisMod_out,       file=paste(dir1, "/E9_AHe_2014NC_HisMod_out.txt",sep=""),            quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)



### H2BGFP(EEDhete, KO),  H3K27me3(EEDhete, KO)
myForCluster1A  <- cbind(myRowNames3, myRowNames4, E1_week0_EEDheto_sort,  zeroCol,   E1_week0_EEDko_sort,  zeroCol, E1_week4_EEDheto_sort,  zeroCol,   E1_week4_EEDko_sort ,    zeroCol,   
                         E7_H3K27me3_WT_sort,      zeroCol,   E7_H3K27me3_EEDko_sort)
write.table( myForCluster1A,        file=paste(dir1, "/myForCluster1A.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)



### H2BGFP(EEDhete, KO),  H3K27ac(EEDhete, KO)
myForCluster1B  <- cbind(myRowNames3, myRowNames4, E1_week0_EEDheto_sort,  zeroCol,   E1_week0_EEDko_sort,    zeroCol, E1_week4_EEDheto_sort,  zeroCol,   E1_week4_EEDko_sort ,   zeroCol,   
                         E6_H3K27ac_WT_sort,  zeroCol, E6_H3K27ac_EEDko_sort )
write.table( myForCluster1B,        file=paste(dir1, "/myForCluster1B.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)



### H2BGFP(EEDhete, KO),  H3K27me3(EEDhete, KO),  H3K27ac(EEDhete, KO)
myForCluster1C  <- cbind(myRowNames3, myRowNames4, E1_week0_EEDheto_sort,  zeroCol,   E1_week0_EEDko_sort,    zeroCol, E1_week4_EEDheto_sort,  zeroCol,   E1_week4_EEDko_sort ,       zeroCol,   
                         E7_H3K27me3_WT_sort,      zeroCol,   E7_H3K27me3_EEDko_sort,  zeroCol, E6_H3K27ac_WT_sort,  zeroCol, E6_H3K27ac_EEDko_sort )
write.table( myForCluster1C,        file=paste(dir1, "/myForCluster1C.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)



### H2BGFP(EEDhete, KO),  H3K27me3(EEDhete, KO),  H3K27ac(EEDhete, KO),  MNase(EEDhete, KO), 
myForCluster1D  <- cbind(myRowNames3, myRowNames4, E1_week0_EEDheto_sort,  zeroCol,   E1_week0_EEDko_sort,  zeroCol, E1_week4_EEDheto_sort,  zeroCol,   E1_week4_EEDko_sort ,       zeroCol,   
                         E7_H3K27me3_WT_sort,      zeroCol,   E7_H3K27me3_EEDko_sort,  zeroCol, E6_H3K27ac_WT_sort,  zeroCol, 
                         E6_H3K27ac_EEDko_sort, zeroCol,  E3_Average_MNase_EEDheto_sort,   zeroCol,   E3_Average_MNase_EEDko_sort )
write.table( myForCluster1D,        file=paste(dir1, "/myForCluster1D.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)



### H2BGFP(EEDhete, KO),  H3K27me3(EEDhete, KO),  MNase(EEDhete, KO), 
myForCluster1E  <- cbind(myRowNames3, myRowNames4, E1_week0_EEDheto_sort,  zeroCol,   E1_week0_EEDko_sort,  zeroCol, E1_week4_EEDheto_sort,  zeroCol,   E1_week4_EEDko_sort ,       zeroCol,   
                         E7_H3K27me3_WT_sort,      zeroCol,   E7_H3K27me3_EEDko_sort,  zeroCol, E3_Average_MNase_EEDheto_sort,   zeroCol,   E3_Average_MNase_EEDko_sort )
write.table( myForCluster1E,        file=paste(dir1, "/myForCluster1E.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)



### H2BGFP(EEDhete, KO),   MNase(EEDhete, KO), 
myForCluster1F  <- cbind(myRowNames3, myRowNames4, E1_week0_EEDheto_sort,  zeroCol,   E1_week0_EEDko_sort,  zeroCol, E1_week4_EEDheto_sort,  zeroCol,   E1_week4_EEDko_sort ,       zeroCol,   
                         E3_Average_MNase_EEDheto_sort,   zeroCol,   E3_Average_MNase_EEDko_sort )
write.table( myForCluster1F,        file=paste(dir1, "/myForCluster1F.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)




### H2BGFP(EEDhete, KO),   Fog2, Mef2c, Srf, Tbx5
myForCluster1G  <- cbind(myRowNames3, myRowNames4, E1_week0_EEDheto_sort,  zeroCol,   E1_week0_EEDko_sort,  zeroCol, E1_week4_EEDheto_sort,  zeroCol,   E1_week4_EEDko_sort ,       zeroCol,   
                         E2_360R1_Fog2_adult_sort ,  zeroCol,   E2_360R1_Mef2c_adult_sort,   zeroCol,    E2_360R1_Srf_adult_sort,  zeroCol, E2_360R1_Tbx5_adult_sort )
write.table( myForCluster1G,        file=paste(dir1, "/myForCluster1G.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)




### H2BGFP(EEDhete, KO),   H3K4me1, H3K4me3, H3K27ac, H3K27me3, Pol2
myForCluster1H  <- cbind(myRowNames3, myRowNames4, E1_week0_EEDheto_sort,  zeroCol,   E1_week0_EEDko_sort,  zeroCol, E1_week4_EEDheto_sort,  zeroCol,   E1_week4_EEDko_sort ,       zeroCol,   
                         E9_Average_H3K4me1_Adult,    zeroCol,  E9_Average_H3K4me3_Adult,  zeroCol,  E9_Average_H3K27ac_Adult,    zeroCol,  E9_Average_H3K27me3_Adult, zeroCol, E9_Average_Pol2_Adult )
write.table( myForCluster1H,        file=paste(dir1, "/myForCluster1H.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)



### H2BGFP(EEDhete, KO),   Pol2
myForCluster1I  <- cbind(myRowNames3, myRowNames4, E1_week0_EEDheto_sort,  zeroCol,   E1_week0_EEDko_sort,  zeroCol, E1_week4_EEDheto_sort,  zeroCol,   E1_week4_EEDko_sort ,       zeroCol,   
                         E9_Average_Pol2_Adult )
write.table( myForCluster1I,        file=paste(dir1, "/myForCluster1I.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)



### H2BGFP(WT),   H3K4me1, H3K4me3, H3K27ac, H3K27me3, Pol2
myForCluster1J  <- cbind(myRowNames3, myRowNames4, 
                         E1_H3_WT_adult_sort,      zeroCol,   E1_H2BGFP_WT_week0_sort,  zeroCol,   E1_H2BGFP_WT_week1_sort,  zeroCol,   E1_H2BGFP_WT_week2_sort,  zeroCol,  
                         E1_H2BGFP_WT_week4_sort,  zeroCol,   E1_H2BGFP_WT_week6_sort,  zeroCol,   E1_H2BGFP_WT_week8_sort   ,  zeroCol,  
                         E9_Average_H3K4me1_Adult,    zeroCol,  E9_Average_H3K4me3_Adult,  zeroCol,  E9_Average_H3K27ac_Adult,    zeroCol,  E9_Average_H3K27me3_Adult, zeroCol, E9_Average_Pol2_Adult )
write.table( myForCluster1J,        file=paste(dir1, "/myForCluster1J.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)



### H2BGFP(WT),   Pol2
myForCluster1K  <- cbind(myRowNames3, myRowNames4, 
                         E1_H3_WT_adult_sort,      zeroCol,   E1_H2BGFP_WT_week0_sort,  zeroCol,   E1_H2BGFP_WT_week1_sort,  zeroCol,   E1_H2BGFP_WT_week2_sort,  zeroCol,  
                         E1_H2BGFP_WT_week4_sort,  zeroCol,   E1_H2BGFP_WT_week6_sort,  zeroCol,   E1_H2BGFP_WT_week8_sort   ,  zeroCol,  
                         E9_Average_Pol2_Adult )
write.table( myForCluster1K,        file=paste(dir1, "/myForCluster1K.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)



### H2BGFP(WT),   Fog2, Mef2c, Srf, Tbx5
myForCluster1L  <- cbind(myRowNames3, myRowNames4, 
                         E1_H3_WT_adult_sort,      zeroCol,   E1_H2BGFP_WT_week0_sort,  zeroCol,   E1_H2BGFP_WT_week1_sort,  zeroCol,   E1_H2BGFP_WT_week2_sort,  zeroCol,  
                         E1_H2BGFP_WT_week4_sort,  zeroCol,   E1_H2BGFP_WT_week6_sort,  zeroCol,   E1_H2BGFP_WT_week8_sort   ,  zeroCol,  
                         E2_360R1_Fog2_adult_sort ,  zeroCol,   E2_360R1_Mef2c_adult_sort,   zeroCol,    E2_360R1_Srf_adult_sort,  zeroCol, E2_360R1_Tbx5_adult_sort  )
write.table( myForCluster1L,        file=paste(dir1, "/myForCluster1L.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)



### H2BGFP(EEDhete, KO),    H3K27me3(EEDhete, KO), H3K4me3 (WT) 
myForCluster1M  <- cbind(myRowNames3, myRowNames4, E1_week0_EEDheto_sort,  zeroCol,   E1_week0_EEDko_sort,  zeroCol, E1_week4_EEDheto_sort,  zeroCol,   E1_week4_EEDko_sort ,       zeroCol,   
                         E7_H3K27me3_WT_sort,      zeroCol,   E7_H3K27me3_EEDko_sort,  zeroCol, E9_Average_H3K4me3_Adult)
write.table( myForCluster1M,        file=paste(dir1, "/myForCluster1M.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)



### H2BGFP(EEDhete, KO),    H3K27me3(EEDhete), H3K4me3 (WT) 
myForCluster1N  <- cbind(myRowNames3, myRowNames4, E1_week0_EEDheto_sort,  zeroCol,   E1_week0_EEDko_sort,  zeroCol, E1_week4_EEDheto_sort,  zeroCol,   E1_week4_EEDko_sort ,       zeroCol,   
                         E9_Average_H3K27me3_Adult,   zeroCol,   E9_Average_H3K4me3_Adult)
write.table( myForCluster1N,        file=paste(dir1, "/myForCluster1N.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)


### H3K27me3(EEDhete, KO), H3K4me3 (WT) 
myForCluster1O  <- cbind(myRowNames3, myRowNames4,   E9_Average_H3K27me3_Adult,   zeroCol,   E9_Average_H3K4me3_Adult)
write.table( myForCluster1O,        file=paste(dir1, "/myForCluster1O.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)



### H3K27me3(EEDhete, KO), H3K4me3 (WT) , H3K27ac (WT) 
myForCluster1P  <- cbind(myRowNames3, myRowNames4,   E9_Average_H3K27me3_Adult,   zeroCol,   E9_Average_H3K4me3_Adult,  zeroCol,   E9_Average_H3K27ac_Adult)
write.table( myForCluster1P,        file=paste(dir1, "/myForCluster1P.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)



### H2BGFP(WT),   H3K4me1, H3K4me3, H3K27ac, H3K27me3, Pol2
myForCluster1Q  <- cbind(myRowNames3, myRowNames4, 
                         E1_H3_WT_adult_sort,      zeroCol,   E1_H2BGFP_WT_week0_sort,  zeroCol,   E1_H2BGFP_WT_week1_sort,  zeroCol,   E1_H2BGFP_WT_week2_sort,  zeroCol,  
                         E1_H2BGFP_WT_week4_sort,  zeroCol,   E1_H2BGFP_WT_week6_sort,  zeroCol,   E1_H2BGFP_WT_week8_sort   ,  zeroCol,  
                         E9_Average_H3K27me3_Adult, zeroCol,   E9_Average_H3K4me3_Adult )
write.table( myForCluster1Q,        file=paste(dir1, "/myForCluster1Q.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)




####################################################################################################################################################################################################################################








####################################################################################################################################################################################################################################
dir2 <- paste(AllResults_g,  "/2_merged_5kb",  sep="")
if( ! file.exists(dir2) ) { dir.create(path=dir2,   recursive = TRUE) }

numOfRows <- nrow(A1_H2BGFP_H3_Rep1)
zeroCol   <- cbind( rep("NA", numOfRows ), rep("NA", numOfRows ),  rep("NA", numOfRows ), rep("NA", numOfRows ),  rep("NA", numOfRows ), rep("NA", numOfRows ), rep("NA", numOfRows ),  rep("NA", numOfRows ), rep("NA", numOfRows ),  rep("NA", numOfRows ) )                               

mode(myRowNames2)  
myRowNames3 <- c("name1", as.character(myRowNames2) ) 
myRowNames4 <- c("name2", as.character(myRowNames2) ) 
length(myRowNames3)
length(myRowNames4)

myFirstLine <- A1_H2BGFP_H3_Rep1[1, c(5:504)]
names(myFirstLine)
rownames(C1_Average_H3)
colnames(C1_Average_H3)

F0_Gene_TPM_out <- A0_Gene_TPM_Rep1[c(1, index1+1),]

F1_H3_WT_adult_sort      <-  	rbind( myFirstLine ,    C1_Average_H3[index1,    ]    )
F1_H2BGFP_WT_week0_sort  <-  	rbind( myFirstLine ,    C1_Average_week0[index1, ]    )
F1_H2BGFP_WT_week1_sort  <-   rbind( myFirstLine ,    C1_Average_week1[index1, ]    ) 
F1_H2BGFP_WT_week2_sort  <-  	rbind( myFirstLine ,    C1_Average_week2[index1, ]    )
F1_H2BGFP_WT_week4_sort  <-  	rbind( myFirstLine ,    C1_Average_week4[index1, ]    )
F1_H2BGFP_WT_week6_sort  <-  	rbind( myFirstLine ,    C1_Average_week6[index1, ]    )
F1_H2BGFP_WT_week8_sort  <-   rbind( myFirstLine ,    C1_Average_week8[index1, ]    ) 
F1_H2BGFP_WT_out         <-   cbind(myRowNames3, myRowNames4,  F1_H3_WT_adult_sort,      zeroCol,   F1_H2BGFP_WT_week0_sort,  zeroCol, 
                                    F1_H2BGFP_WT_week1_sort,  zeroCol,   F1_H2BGFP_WT_week2_sort,  zeroCol,  
                                    F1_H2BGFP_WT_week4_sort,  zeroCol,   F1_H2BGFP_WT_week6_sort,  zeroCol,   F1_H2BGFP_WT_week8_sort )

F1_week0_EEDheto_sort      <-  	rbind( myFirstLine ,    C1_Average_week0_EEDheto[index1,  ]    )
F1_week0_EEDko_sort        <-  	rbind( myFirstLine ,    C1_Average_week0_EEDko[index1,    ]    )
F1_week4_EEDheto_sort      <-  	rbind( myFirstLine ,    C1_Average_week4_EEDheto[index1,  ]    )
F1_week4_EEDko_sort        <-  	rbind( myFirstLine ,    C1_Average_week4_EEDko[index1,    ]    )
F1_H2BGFP_CKO_out         <-   cbind(myRowNames3, myRowNames4, F1_week0_EEDheto_sort,  zeroCol,   F1_week0_EEDko_sort,  zeroCol, 
                                     F1_week4_EEDheto_sort,  zeroCol,   F1_week4_EEDko_sort )

F1_banding_sort     <-  	rbind( myFirstLine ,    C1_Average_banding[index1,  ]    )
F1_sham_sort        <-  	rbind( myFirstLine ,    C1_Average_sham[index1,    ]    )
F1_H2BGFP_TAC_out   <-    cbind(myRowNames3, myRowNames4,  F1_sham_sort,  zeroCol,   F1_banding_sort )



F2_360R1_Fog2_adult_sort     <-  rbind( myFirstLine ,    C2_Average_360R1_Fog2_adult[index1,  ]    )
F2_360R1_Fog2_E12_sort       <-  rbind( myFirstLine ,    C2_Average_360R1_Fog2_E12[index1,  ]    )
F2_360R1_Mef2c_adult_sort    <-  rbind( myFirstLine ,    C2_Average_360R1_Mef2c_adult[index1, ]    )
F2_360R1_Mef2c_E12_sort      <-  rbind( myFirstLine ,    C2_Average_360R1_Mef2c_E12[index1, ]    )
F2_360R1_Srf_adult_sort      <-  rbind( myFirstLine ,    C2_Average_360R1_Srf_adult[index1, ]    )
F2_360R1_Tbx5_adult_sort     <-  rbind( myFirstLine ,    C2_Average_360R1_Tbx5_adult[index1,] )
F2_360R1_Tbx5_E12_sort       <-  rbind( myFirstLine ,    C2_Average_360R1_Tbx5_E12[index1,] )
F2_360R1_adult_out           <-  cbind(myRowNames3, myRowNames4,  F2_360R1_Fog2_adult_sort ,  zeroCol,   F2_360R1_Mef2c_adult_sort,   zeroCol,    F2_360R1_Srf_adult_sort,  zeroCol, F2_360R1_Tbx5_adult_sort )
F2_360R1_E12_out             <-  cbind(myRowNames3, myRowNames4,  F2_360R1_Fog2_E12_sort ,    zeroCol,   F2_360R1_Mef2c_E12_sort,     zeroCol,    F2_360R1_Tbx5_E12_sort )




F3_Average_MNase_EEDheto_sort     <-  	rbind( myFirstLine ,    C3_Average_MNase_EEDheto[index1,  ]    )
F3_Average_MNase_EEDko_sort       <-  	rbind( myFirstLine ,    C3_Average_MNase_EEDko[index1,    ]    )
F3_Average_MNase_out              <-    cbind(myRowNames3, myRowNames4,  F3_Average_MNase_EEDheto_sort,  zeroCol,   F3_Average_MNase_EEDko_sort )

F3_merge_MNase_EEDheto_sort     <-  	rbind( myFirstLine ,    C3_Average_MNase_EEDheto_merge[index1,  ]    )
F3_merge_MNase_EEDko_sort       <-  	rbind( myFirstLine ,    C3_Average_MNase_EEDko_merge[index1,    ]    )
F3_merge_MNase_out              <-    cbind(myRowNames3, myRowNames4,  F3_merge_MNase_EEDheto_sort,  zeroCol,   F3_merge_MNase_EEDko_sort )

F4_EED_Adult_sort   <-  	rbind( myFirstLine ,    C4_Average_EED_Adult[index1,  ]    )
F4_EED_P5_sort      <-  	rbind( myFirstLine ,    C4_Average_EED_P5[index1,    ]    )
F4_EED_out          <-    cbind(myRowNames3, myRowNames4,  F4_EED_P5_sort,  zeroCol,   F4_EED_Adult_sort )

F5_HetoHDAC1_sort   <-  	rbind( myFirstLine ,    C5_Average_HetoHDAC1[index1,  ]    )
F5_HetoHDAC2_sort   <-  	rbind( myFirstLine ,    C5_Average_HetoHDAC2[index1,  ]    )
F5_HomoHDAC1_sort   <-  	rbind( myFirstLine ,    C5_Average_HomoHDAC1[index1,  ]    )
F5_HomoHDAC2_sort   <-  	rbind( myFirstLine ,    C5_Average_HomoHDAC2[index1,  ]    )
F5_HDAC_out   <-   cbind(myRowNames3, myRowNames4, 
                         F5_HetoHDAC1_sort,  zeroCol,   F5_HomoHDAC1_sort,  zeroCol, 
                         F5_HetoHDAC2_sort,  zeroCol,   F5_HomoHDAC2_sort )

F6_H3K27ac_WT_sort         <-  	rbind( myFirstLine ,    C6_Average_H3K27ac_WT[index1,  ]    )
F6_H3K27ac_EEDko_sort      <-  	rbind( myFirstLine ,    C6_Average_H3K27ac_EEDko[index1,    ]    )
F6_H3K27ac_out             <-   cbind(myRowNames3, myRowNames4,  F6_H3K27ac_WT_sort,  zeroCol,   F6_H3K27ac_EEDko_sort )

F7_H3K27me3_WT_sort         <-  rbind( myFirstLine ,    C7_Average_H3K27me3_WT[index1,  ]    )
F7_H3K27me3_EEDko_sort      <-  rbind( myFirstLine ,    C7_Average_H3K27me3_EEDko[index1,    ]    )
F7_H3K27me3_out             <-  cbind(myRowNames3, myRowNames4,  F7_H3K27me3_WT_sort,  zeroCol,   F7_H3K27me3_EEDko_sort )




F8_Average_H3K4me1_Adult    <-  rbind( myFirstLine ,    C8_Average_H3K4me1_Adult[index1,  ]    )
F8_Average_H3K4me3_Adult    <-  rbind( myFirstLine ,    C8_Average_H3K4me3_Adult[index1,  ]    )
F8_Average_H3K27ac_Adult    <-  rbind( myFirstLine ,    C8_Average_H3K27ac_Adult[index1,  ]    )
F8_Average_H3K27me3_Adult   <-  rbind( myFirstLine ,    C8_Average_H3K27me3_Adult[index1,  ]    )
F8_Average_MeCP2_Adult      <-  rbind( myFirstLine ,    C8_Average_MeCP2_Adult[index1,  ]    )
F8_Average_MeCP2_TAC        <-  rbind( myFirstLine ,    C8_Average_MeCP2_TAC[index1,  ]    )
F8_DNAme_2014NC_out         <-  cbind(myRowNames3, myRowNames4,  
                                      F8_Average_H3K4me1_Adult,  zeroCol,   F8_Average_H3K4me3_Adult, zeroCol,   F8_Average_H3K27ac_Adult,   zeroCol,  
                                      F8_Average_H3K27me3_Adult, zeroCol,   F8_Average_MeCP2_Adult,   zeroCol,   F8_Average_MeCP2_TAC )



F9_Average_GATA4_Ab_Adult    <-  rbind( myFirstLine ,    C9_Average_GATA4_Ab_Adult[index1,  ]    )
F9_Average_GATA4_fb_Adult    <-  rbind( myFirstLine ,    C9_Average_GATA4_fb_Adult[index1,  ]    )
F9_Average_GATA4_fb_Banding  <-  rbind( myFirstLine ,    C9_Average_GATA4_fb_Banding[index1,  ]    )
F9_Average_GATA4_fb_Sham     <-  rbind( myFirstLine ,    C9_Average_GATA4_fb_Sham[index1,  ]    )
F9_Average_H3K4me1_Adult     <-  rbind( myFirstLine ,    C9_Average_H3K4me1_Adult[index1,  ]    )
F9_Average_H3K4me3_Adult     <-  rbind( myFirstLine ,    C9_Average_H3K4me3_Adult[index1,  ]    )
F9_Average_H3K27ac_Adult     <-  rbind( myFirstLine ,    C9_Average_H3K27ac_Adult[index1,  ]    )
F9_Average_Pol2_Adult        <-  rbind( myFirstLine ,    C9_Average_Pol2_Adult[index1,  ]    )
F9_Average_H3K27me3_Adult    <-  rbind( myFirstLine ,    C9_Average_H3K27me3_Adult[index1,  ]    )
F9_AHe_2014NC_GATA4_out      <-  cbind(myRowNames3, myRowNames4, F9_Average_GATA4_Ab_Adult,  zeroCol,  F9_Average_GATA4_fb_Adult, zeroCol,  F9_Average_GATA4_fb_Banding, zeroCol,  F9_Average_GATA4_fb_Sham )
F9_AHe_2014NC_HisMod_out     <-  cbind(myRowNames3, myRowNames4, 
                                       F9_Average_H3K4me1_Adult,    zeroCol,  F9_Average_H3K4me3_Adult,  zeroCol,  
                                       F9_Average_H3K27ac_Adult,    zeroCol,  F9_Average_H3K27me3_Adult, zeroCol, F9_Average_Pol2_Adult)





dim(F0_Gene_TPM_out)
dim(F1_H2BGFP_WT_out)
dim(F1_H2BGFP_CKO_out)
dim(F1_H2BGFP_TAC_out)
dim(F2_360R1_adult_out)
dim(F2_360R1_E12_out)
dim(F3_Average_MNase_out)
dim(F3_merge_MNase_out)
dim(F4_EED_out )
dim(F5_HDAC_out)
dim(F6_H3K27ac_out)
dim(F7_H3K27me3_out)
dim(F8_DNAme_2014NC_out)
dim(F9_AHe_2014NC_GATA4_out)
dim(F9_AHe_2014NC_HisMod_out)

write.table( F0_Gene_TPM_out,        file=paste(dir2, "/F0_Gene_TPM_out.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( F1_H2BGFP_WT_out,       file=paste(dir2, "/F1_H2BGFP_WT_out.txt",    sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( F1_H2BGFP_CKO_out,      file=paste(dir2, "/F1_H2BGFP_CKO_out.txt",   sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( F1_H2BGFP_TAC_out,      file=paste(dir2, "/F1_H2BGFP_TAC_out.txt",   sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( F2_360R1_adult_out,     file=paste(dir2, "/F2_360R1_adult_out.txt",  sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( F2_360R1_E12_out,       file=paste(dir2, "/F2_360R1_E12_out.txt",    sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( F3_Average_MNase_out,   file=paste(dir2, "/F3_Average_MNase_out.txt",sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( F3_merge_MNase_out,     file=paste(dir2, "/F3_merge_MNase_out.txt",  sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( F4_EED_out,             file=paste(dir2, "/F4_EED_out.txt",          sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( F5_HDAC_out,            file=paste(dir2, "/F5_HDAC_out.txt",         sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( F6_H3K27ac_out,         file=paste(dir2, "/F6_H3K27ac_out.txt",      sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( F7_H3K27me3_out,        file=paste(dir2, "/F7_H3K27me3_out.txt",     sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( F8_DNAme_2014NC_out,    file=paste(dir2, "/F8_DNAme_2014NC_out.txt", sep=""),      quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( F9_AHe_2014NC_GATA4_out, file=paste(dir2, "/F9_AHe_2014NC_GATA4_out.txt",sep=""),      quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( F9_AHe_2014NC_HisMod_out,     file=paste(dir2, "/F9_AHe_2014NC_HisMod_out.txt",sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)

####################################################################################################################################################################################################################################







####################################################################################################################################################################################################################################
dir3 <- paste(AllResults_g,  "/3_biogicalReplicates_5kb",  sep="")
if( ! file.exists(dir3) ) { dir.create(path=dir3,   recursive = TRUE) }

numOfRows <- nrow(A1_H2BGFP_H3_Rep1)
zeroCol   <- cbind( rep("NA", numOfRows ), rep("NA", numOfRows ),  rep("NA", numOfRows ), rep("NA", numOfRows ),  rep("NA", numOfRows ), rep("NA", numOfRows ), rep("NA", numOfRows ),  rep("NA", numOfRows ), rep("NA", numOfRows ),  rep("NA", numOfRows ) )                               

mode(myRowNames2)  
myRowNames3 <- c("name1", as.character(myRowNames2) ) 
myRowNames4 <- c("name2", as.character(myRowNames2) ) 
length(myRowNames3)
length(myRowNames4)

myFirstLine <- A1_H2BGFP_H3_Rep1[1, c(5:504)]
names(myFirstLine)
rownames(C1_Average_H3)
colnames(C1_Average_H3)

G0_Gene_TPM_out <- A0_Gene_TPM_Rep1[c(1, index1+1),]

G1_H2BGFP_H3_Rep1    <- rbind( myFirstLine ,    B1_H2BGFP_H3_Rep1[index1, ]    ) 
G1_H2BGFP_week0_Rep1 <- rbind( myFirstLine ,    B1_H2BGFP_week0_Rep1[index1, ]    ) 
G1_H2BGFP_week0_Rep2 <- rbind( myFirstLine ,    B1_H2BGFP_week0_Rep2[index1, ]    ) 
G1_H2BGFP_week1_Rep1 <- rbind( myFirstLine ,    B1_H2BGFP_week1_Rep1[index1, ]    ) 
G1_H2BGFP_week1_Rep2 <- rbind( myFirstLine ,    B1_H2BGFP_week1_Rep2[index1, ]    ) 
G1_H2BGFP_week2_Rep1 <- rbind( myFirstLine ,    B1_H2BGFP_week2_Rep1[index1, ]    ) 
G1_H2BGFP_week2_Rep2 <- rbind( myFirstLine ,    B1_H2BGFP_week2_Rep2[index1, ]    ) 
G1_H2BGFP_week4_Rep1 <- rbind( myFirstLine ,    B1_H2BGFP_week4_Rep1[index1, ]    ) 
G1_H2BGFP_week4_Rep2 <- rbind( myFirstLine ,    B1_H2BGFP_week4_Rep2[index1, ]    )
G1_H2BGFP_week4_Rep3 <- rbind( myFirstLine ,    B1_H2BGFP_week4_Rep3[index1, ]    ) 
G1_H2BGFP_week6_Rep1 <- rbind( myFirstLine ,    B1_H2BGFP_week6_Rep1[index1, ]    ) 
G1_H2BGFP_week6_Rep2 <- rbind( myFirstLine ,    B1_H2BGFP_week6_Rep2[index1, ]    ) 
G1_H2BGFP_week8_Rep1 <- rbind( myFirstLine ,    B1_H2BGFP_week8_Rep1[index1, ]    ) 
G1_H2BGFP_week0_EEDheto_Rep1 <- rbind( myFirstLine ,    B1_H2BGFP_week0_EEDheto_Rep1[index1, ]    ) 
G1_H2BGFP_week0_EEDheto_Rep2 <- rbind( myFirstLine ,    B1_H2BGFP_week0_EEDheto_Rep2[index1, ]    ) 
G1_H2BGFP_week0_EEDko_Rep1   <- rbind( myFirstLine ,    B1_H2BGFP_week0_EEDko_Rep1[index1, ]    ) 
G1_H2BGFP_week4_EEDheto_Rep1 <- rbind( myFirstLine ,    B1_H2BGFP_week4_EEDheto_Rep1[index1, ]    ) 
G1_H2BGFP_week4_EEDheto_Rep2 <- rbind( myFirstLine ,    B1_H2BGFP_week4_EEDheto_Rep2[index1, ]    ) 
G1_H2BGFP_week4_EEDko_Rep1   <- rbind( myFirstLine ,    B1_H2BGFP_week4_EEDko_Rep1[index1, ]    ) 
G1_H2BGFP_week4_EEDko_Rep2   <- rbind( myFirstLine ,    B1_H2BGFP_week4_EEDko_Rep2[index1, ]    ) 
G1_H2BGFP_banding_Rep1 <- rbind( myFirstLine ,    B1_H2BGFP_banding_Rep1[index1, ]    ) 
G1_H2BGFP_banding_Rep2 <- rbind( myFirstLine ,    B1_H2BGFP_banding_Rep2[index1, ]    ) 
G1_H2BGFP_sham_Rep1    <- rbind( myFirstLine ,    B1_H2BGFP_sham_Rep1[index1, ]    ) 
G1_H2BGFP_sham_Rep2    <- rbind( myFirstLine ,    B1_H2BGFP_sham_Rep2[index1, ]    ) 

G1_H2BGFP_WT_out         <-   cbind(myRowNames3, myRowNames4, 
                                    G1_H2BGFP_H3_Rep1,      zeroCol,   G1_H2BGFP_week0_Rep1,  zeroCol, G1_H2BGFP_week0_Rep2, zeroCol,  G1_H2BGFP_week1_Rep1, zeroCol, G1_H2BGFP_week1_Rep2,  zeroCol,   
                                    G1_H2BGFP_week2_Rep1,   zeroCol,   G1_H2BGFP_week2_Rep2,  zeroCol, G1_H2BGFP_week4_Rep1, zeroCol,  G1_H2BGFP_week4_Rep2, zeroCol, G1_H2BGFP_week4_Rep3,  zeroCol,    
                                    G1_H2BGFP_week6_Rep2,   zeroCol,   G1_H2BGFP_week8_Rep1  )

G1_H2BGFP_EEDcKO_out    <-   cbind(myRowNames3, myRowNames4, 
                                    G1_H2BGFP_week0_EEDheto_Rep1,   zeroCol,   G1_H2BGFP_week0_EEDheto_Rep2,  zeroCol, G1_H2BGFP_week0_EEDko_Rep1,  zeroCol, 
                                    G1_H2BGFP_week4_EEDheto_Rep2,   zeroCol,   G1_H2BGFP_week4_EEDko_Rep1,    zeroCol, G1_H2BGFP_week4_EEDko_Rep2 )

G1_H2BGFP_TAC_out         <-   cbind(myRowNames3, myRowNames4, G1_H2BGFP_sham_Rep1,      zeroCol,   G1_H2BGFP_sham_Rep2,  zeroCol,  G1_H2BGFP_banding_Rep1, zeroCol, G1_H2BGFP_banding_Rep2  )




G3_MNase_EEDheto_Rep1  <- rbind( myFirstLine ,    B3_MNase_EEDheto_Rep1[index1, ]    )   
G3_MNase_EEDheto_Rep2  <- rbind( myFirstLine ,    B3_MNase_EEDheto_Rep2[index1, ]    )    
G3_MNase_EEDko_Rep1    <- rbind( myFirstLine ,    B3_MNase_EEDko_Rep1[index1, ]    )     
G3_MNase_EEDko_Rep2    <- rbind( myFirstLine ,    B3_MNase_EEDko_Rep2[index1, ]    )     
G3_MNase_out           <- cbind(myRowNames3, myRowNames4, G3_MNase_EEDheto_Rep1,      zeroCol,   G3_MNase_EEDheto_Rep2,  zeroCol, G3_MNase_EEDko_Rep1, zeroCol, G3_MNase_EEDko_Rep2  )


G4_EED_Adult_Rep1 <- rbind( myFirstLine ,    B4_EED_Adult_Rep1[index1, ]    )                                          
G4_EED_Adult_Rep2 <- rbind( myFirstLine ,    B4_EED_Adult_Rep2[index1, ]    )                                    
G4_EED_P5_Rep1    <- rbind( myFirstLine ,    B4_EED_P5_Rep1[index1, ]    )                                 
G4_EED_P5_Rep2    <- rbind( myFirstLine ,    B4_EED_P5_Rep2[index1, ]    )                                           
G4_EED_out      <- cbind(myRowNames3, myRowNames4, G4_EED_P5_Rep1,      zeroCol,   G4_EED_P5_Rep2,  zeroCol, G4_EED_Adult_Rep1, zeroCol, G4_EED_Adult_Rep2  )


G5_HDAC_HetoHDAC1 <- rbind( myFirstLine ,    B5_HDAC_HetoHDAC1[index1, ]    )      
G5_HDAC_HetoHDAC2 <- rbind( myFirstLine ,    B5_HDAC_HetoHDAC2[index1, ]    )      
G5_HDAC_HomoHDAC1 <- rbind( myFirstLine ,    B5_HDAC_HomoHDAC1[index1, ]    )     
G5_HDAC_HomoHDAC2 <- rbind( myFirstLine ,    B5_HDAC_HomoHDAC2[index1, ]    )    
G5_HDAC_out      <- cbind(myRowNames3, myRowNames4, G5_HDAC_HetoHDAC1,  zeroCol,   G5_HDAC_HomoHDAC1, zeroCol, 
                          G5_HDAC_HetoHDAC2,  zeroCol, G5_HDAC_HomoHDAC2  )


G6_H3K27ac_EEDko_rep1 <- rbind( myFirstLine ,    B6_H3K27ac_EEDko_rep1[index1, ]    )  
G6_H3K27ac_EEDko_rep2 <- rbind( myFirstLine ,    B6_H3K27ac_EEDko_rep2[index1, ]    )  
G6_H3K27ac_WT_rep1    <- rbind( myFirstLine ,    B6_H3K27ac_WT_rep1[index1, ]    )      
G6_H3K27ac_WT_rep2    <- rbind( myFirstLine ,    B6_H3K27ac_WT_rep2[index1, ]    ) 
G6_H3K27ac_out        <- cbind(myRowNames3, myRowNames4, G6_H3K27ac_WT_rep1,      zeroCol,   G6_H3K27ac_WT_rep2,  zeroCol, G6_H3K27ac_EEDko_rep1, zeroCol, G6_H3K27ac_EEDko_rep2  )


G7_H3K27me3_EEDko_rep1 <- rbind( myFirstLine ,    B7_H3K27me3_EEDko_rep1[index1, ]    )  
G7_H3K27me3_EEDko_rep2 <- rbind( myFirstLine ,    B7_H3K27me3_EEDko_rep2[index1, ]    )   
G7_H3K27me3_WT_rep1    <- rbind( myFirstLine ,    B7_H3K27me3_WT_rep1[index1, ]    )      
G7_H3K27me3_WT_rep2    <- rbind( myFirstLine ,    B7_H3K27me3_WT_rep2[index1, ]    )   
G7_H3K27me3_WT_rep3    <- rbind( myFirstLine ,    B7_H3K27me3_WT_rep3[index1, ]    ) 
G7_H3K27me3_out      <- cbind(myRowNames3, myRowNames4, G7_H3K27me3_WT_rep1, zeroCol,   G7_H3K27me3_WT_rep2,  zeroCol, G7_H3K27me3_WT_rep3,
                           zeroCol, G7_H3K27me3_EEDko_rep1,   zeroCol, G7_H3K27me3_EEDko_rep2)


dim(G0_Gene_TPM_out)
dim(G1_H2BGFP_WT_out)
dim(G1_H2BGFP_EEDcKO_out)
dim(G1_H2BGFP_TAC_out)
dim(G3_MNase_out)
dim(G4_EED_out)
dim(G5_HDAC_out)
dim(G6_H3K27ac_out)
dim(G7_H3K27me3_out)

write.table( G0_Gene_TPM_out,                     file=paste(dir3, "/G0_Gene_TPM_out.txt",    sep=""),         quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( G1_H2BGFP_WT_out,                    file=paste(dir3, "/G1_H2BGFP_WT_out.txt",    sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( G1_H2BGFP_EEDcKO_out,                file=paste(dir3, "/G1_H2BGFP_EEDcKO_out.txt",   sep=""),     quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( G1_H2BGFP_TAC_out,                   file=paste(dir3, "/G1_H2BGFP_TAC_out.txt",   sep=""),        quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( G3_MNase_out,                        file=paste(dir3, "/G3_MNase_out.txt",sep=""),                quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( G4_EED_out,                          file=paste(dir3, "/G4_EED_out.txt",  sep=""),                quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( G5_HDAC_out,                         file=paste(dir3, "/G5_HDAC_out.txt",  sep=""),               quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( G6_H3K27ac_out,                      file=paste(dir3, "/G6_H3K27ac_out.txt",          sep=""),    quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)
write.table( G7_H3K27me3_out,                     file=paste(dir3, "/G7_H3K27me3_out.txt",         sep=""),    quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = FALSE)


####################################################################################################################################################################################################################################
















