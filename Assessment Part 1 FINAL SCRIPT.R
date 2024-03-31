#Assessment Part 1#

####Loading tables####
em_counts=read.table("C:\\Users\\Samarth Islam Monaz\\OneDrive\\Desktop\\R2 Assessment\\counts.csv",header=TRUE,row.names=1,sep="\t")
ss = read.table("C:\\Users\\Samarth Islam Monaz\\OneDrive\\Desktop\\R2 Assessment\\ss.csv",header=TRUE,row.names=1,sep="\t")
ss_patient= read.table("C:\\Users\\Samarth Islam Monaz\\OneDrive\\Desktop\\R2 Assessment\\ss_per_patient.csv",header=TRUE,row.names=1,sep="\t")

####Summary Stats####

#Calculating the means of each quantitative columns 
Age_mean = mean(ss$age)
#BMI
BMI_mean = mean(ss$bmi)
#slide area
slidearea_mean = mean(ss$slide_area)
#days survived
survival_mean=mean(ss$days_survived)


####Turning continuous variables discrete####

#age
age_discrete = cut(labels=c("low","high"),ss$age,2) 
age_discrete 
ss$age_discrete  = age_discrete 

#slide_area
slide_area_discrete=cut(labels=c("low","high"),ss$slide_area,2)
slide_area_discrete
ss$slide_area_discrete = slide_area_discrete

#bmi
bmi_discrete=cut(labels=c("low","high"),ss$bmi,2)
bmi_discrete
ss$bmi_discrete = bmi_discrete

#Converting the discrete columns to a factor
ss$sample_group=as.factor(ss$sample_group)
ss$batch=as.factor(ss$batch)
ss$patient=as.factor(ss$patient)

#Modifying the order of which our groups are handled by specifying levels
ss$sample_group= factor(ss$sample_group, levels = c("HC", "GB_2", "GB_1"))
summary_of_ss=summary(ss)

#Splitting our data up by the SAMPLE_GROUP 
ss_HEALTHY = ss[ss$sample_group == "HC",]
ss_HEALTHY

ss_GB_2=ss[ss$sample_group == "GB_2",]
ss_GB_2

ss_GB_1=ss[ss$sample_group == "GB_1",]
ss_GB_1

#calculate summary stats for each group and save as a table.
summary(ss_HEALTHY)
summary(ss_GB_2)
summary(ss_GB_1)


#Looking at dispersion of the data by calculating the S. deviation for each of the variables for each sample group

#For HEALTHY group
sd(ss_HEALTHY$age)
sd(ss_HEALTHY$bmi)
sd(ss_HEALTHY$slide_area)
sd(ss_HEALTHY$days_survived)

#For GB_1 group
sd(ss_GB_1$age)
sd(ss_GB_1$bmi)
sd(ss_GB_1$slide_area)
sd(ss_GB_1$days_survived)

#For GB_2 group
sd(ss_GB_2$age)
sd(ss_GB_2$bmi)
sd(ss_GB_2$slide_area)
sd(ss_GB_2$days_survived)


####Observing differences in our CONTINUOUS VARIABLES using ANOVA####

#For AGE
AGE_Anova <- aov(age ~ sample_group, data = ss)
summary(AGE_Anova)

#For BMI
BMI_Anova <- aov(bmi~ sample_group, data = ss)
summary(BMI_Anova)

#For SURVIVAL DAYS
Surv_Anova <- aov(days_survived~ sample_group, data = ss)
summary(Surv_Anova)

#FOR SLIDE AREA
Slide_Anova <- aov(slide_area~ sample_group, data = ss)
summary(Slide_Anova)


####Getting adjusted p value for each comparison between the different sample groups using TUKEY####

#For AGE
Age_Tukey=TukeyHSD(AGE_Anova)
Age_Tukey

#For BMI
BMI_Tukey=TukeyHSD(BMI_Anova )
BMI_Tukey

#For SURVIVAL DAYS
Surv_Tukey=TukeyHSD(Surv_Anova)
Surv_Tukey

#FOR SLIDE AREA
Slide_Tukey=TukeyHSD(Slide_Anova) 
Slide_Tukey


####Normalizing raw counts data using DSEQ####

#Installing and loading the library
BiocManager::install("DESeq2")
library(DESeq2)

#Preparing the counts by removing non-expressed genes and convert to matrix
em_countsnorm = subset(em_counts,apply(em_counts, 1, mean) >= 1)
em_countsnorm = as.matrix(em_countsnorm)

#Now ready to run DESeq2 to normalize

#sample sheet sample_group column must be a factor. 
ss$sample_group = factor(ss$sample_group)
dds = DESeqDataSetFromMatrix(countData=em_countsnorm, colData=ss, design=~sample_group)
dds = DESeq(dds)

#extracting from dds
normalized_counts = as.data.frame(counts(dds, normalized=TRUE))

#rounding the data (makes it easier to look at)
normalized_counts = round(normalized_counts,2)


####Making PCA using normalized counts####
#PC1 vs PC2 function
make_pc1_pc2 = function(colour_groups, e_data)
{
  # scale data
  e_data_scaled = na.omit(data.frame(t(scale(t(e_data)))))
  
  # run PCA
  xx = prcomp(t(e_data_scaled))
  pca_coordinates = data.frame(xx$x)
  
  # get % variation
  vars = apply(xx$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  x_axis_label = paste("PC1 (" ,prop_x, "%)", sep="")
  y_axis_label = paste("PC2 (" ,prop_y, "%)", sep="")
  
  # plot  
  ggpFNC = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = colour_groups)) +
    geom_point(size = 3) +
    labs(title = "PCA", x= x_axis_label, y= y_axis_label) +
    theme_bw()
  
  return(ggpFNC)
}


#Plotting PC1vsPC2 for variables
library(ggplot2)
pca_counts_batch= make_pc1_pc2(ss$batch, normalized_counts)
pca_counts_batch

pca_counts_slide= make_pc1_pc2(ss$slide_area_discrete, normalized_counts)
pca_counts_slide

pca_counts_age= make_pc1_pc2(ss$age_discrete, normalized_counts)
pca_counts_age

pca_counts_bmi= make_pc1_pc2(ss$bmi_discrete, normalized_counts)
pca_counts_bmi


####Correcting for SLIDE AREA####

#Loading library
library(sva)

#Corrected code makes sure the batch and group information is a factor
sample_group = factor(ss$sample_group)
batch = factor(ss$slide_area_discrete)

#Correcting the counts table using CombatSeq
counts_corrected_slidearea = ComBat_seq(as.matrix(em_counts), batch=slide_area_discrete, group=sample_group)
counts_corrected_slidearea= data.frame(counts_corrected)

#Running DESEQ2 with corrected data to normalize
dds_corrected_slide = DESeqDataSetFromMatrix(countData=counts_corrected_slidearea, colData=ss, design=~sample_group)
dds_corrected_slide = DESeq(dds_corrected_slide)

#extracting from dds
normalized_correctedcounts_slide = as.data.frame(counts(dds_corrected_slide, normalized=TRUE))

#rounding the data
normalized_correctedcounts_slide = round(normalized_correctedcounts,2)
normalized_correctedcounts_slide

#Corrected PC1vsPC2 for slide area
pca_counts_slide_corr= make_pc1_pc2(ss$slide_area_discrete, normalized_correctedcounts_slide)
pca_counts_slide_corr


####Getting corrected DE tables from dds####

#GB1vsHC
de_GB1vsHC = results(dds_corrected_slide , c("sample_group","GB_1","HC"))
de_GB1vsHC = as.data.frame(de_GB1vsHC)
de_GB1vsHC = de_GB1vsHC[order(de_GB1vsHC$padj), ]
#parsing
de_GB1vsHC$ID = row.names(de_GB1vsHC)
de_GB1vsHC = de_GB1vsHC[,c(7,2,5,6)]
colnames(de_GB1vsHC) = c("ID","log2fold","p","p.adj")

#GB2vsHC
de_GB2vsHC = results(dds_corrected_slide , c("sample_group","GB_2","HC"))
de_GB2vsHC = as.data.frame(de_GB2vsHC)
de_GB2vsHC = de_GB2vsHC[order(de_GB2vsHC$padj), ]
#parsing
de_GB2vsHC$ID = row.names(de_GB2vsHC)
de_GB2vsHC = de_GB2vsHC[,c(7,2,5,6)]
colnames(de_GB2vsHC) = c("ID","log2fold","p","p.adj")

#GB1vsGB2
de_GB1vsGB2 = results(dds_corrected_slide , c("sample_group","GB_1","GB_2"))
de_GB1vsGB2 = as.data.frame(de_GB1vsGB2)
de_GB1vsGB2 = de_GB1vsGB2[order(de_GB1vsGB2$padj), ]
#parsing
de_GB1vsGB2$ID = row.names(de_GB1vsGB2)
de_GB1vsGB2 = de_GB1vsGB2[,c(7,2,5,6)]
colnames(de_GB1vsGB2) = c("ID","log2fold","p","p.adj")


#Subsetting for significant genes

#GB1vsHC
de_sig_GB1HC =subset(de_GB1vsHC, p.adj < 0.05 & abs(log2fold) > 1)
nrow(de_sig_GB1HC)

#GB2vsHC
de_sig_GB2HC=subset(de_GB2vsHC, p.adj < 0.05 & abs(log2fold) > 1)
nrow(de_sig_GB2HC)

#GB1vsGB2
de_sig_GB1GB2=subset(de_GB1vsGB2, p.adj < 0.05 & abs(log2fold) > 1)
nrow(de_sig_GB1GB2)


