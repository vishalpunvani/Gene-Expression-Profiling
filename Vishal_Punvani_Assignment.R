#**Note : ALL important points are written as comments between * and *""

## Install and load all dependencies 
install.packages("devtools")
install.packages("dplyr")
install.packages("BiocManager")

library(devtools)
library(BiocManager)
 
BiocManager::install("GSVA")
install_github("cmap/cmapR")
install_github("vqv/ggbiplot")
install_github("sinhrks/ggfortify")

library(cmapR)
library(ggplot2)
library(ggfortify)
library(ggbiplot)
library(dplyr)




#1. Read the gene expression file 
#2. Store column description containing description of each sample in a dataframe
#3. Create a dataframe containing the gene matrix 



#setwd("C:/Users/vishal.punvani/Documents/Finalised Assignment_01.01.2020-20200121T152954Z-001/Finalised Assignment_01.01.2020")
paad_file <- parse_gctx("Copy of PAAD.gct")

column_metadata <- paad_file@cdesc

gene_matrix <- as.data.frame(paad_file@mat)



##Question 1
#remove rows/genes containing any NaN values 
#**4367** genes had NaN values for atleast one sample value


gene_matrix_cleaned <-  gene_matrix[which(apply(gene_matrix,1,function(x) {sum(is.nan(x))})==0),]



##Create all boxplots 
#From the boxplot we can see the spread and center of data
#We can see that gene expressions are normally distributed without with values in the *range -1.8 to 21* and *median at ~9* 


sapply(seq_along(gene_matrix_cleaned), function(i) {
  y <- gene_matrix_cleaned[, i]
  png(filename = paste0(colnames(gene_matrix_cleaned)[i],"_boxplot.png") )
  boxplot(y, outline=FALSE, ylab="gene expression",xlab=(colnames(gene_matrix_cleaned)[i])
  )
  dev.off()
  
})




##Question 2
#For creating and labeling the PCA plot create a copy of paad file named pca_labels 
#Add type_of_tumor column to the column description of pca_labels file to show whether the "histological_type_other column contains neuroendocrine"
#This is done by searching for "nueroendocrine" in the column "histological_type_other" using grep utility
#Also, use pca_data to store the gene matrix which was cleaned in the last step 



pca_labels <- paad_file
pca_labels@cdesc$type_of_tumor <- ifelse(grepl("neuroendocrine",pca_labels@cdesc$histological_type_other,ignore.case = T),"Neuroendocrine","Exocrine")
head(pca_labels@cdesc)
pca_data <- as.data.frame(t(gene_matrix_cleaned))


##Plot of PCA1 against PCA2 
#The position of a point on the plot shows the *direction in which each sample varies the most* 
#The plot shows that the first principle component is not able to explain the variance of the nueroendocrine sample gene expressions
#The high values of second principle component scores for neuroendocrine samples shows that the nueroendocrine samples vary in completely different direction as compared to Exocrine samples
#As can be seen from the blue color dots, the plot shows that *Nueroendocrine tumors are clearly seperable from the Adenocarcinoma Tumor*

pca.out <- prcomp(pca_data,scale. = T)
autoplot(pca.out,data=pca_labels@cdesc,colour='type_of_tumor')+ ggtitle("PCA plot showing difference in direction of variance of \n Neuroendocrine and exocrine samples") 

ggsave(filename = "pca_plot.png")





##Plotting % of variance explained by the principal components
#This plot shows that the first two principal components explain only about 30% of the variance in the data and therefore inorder to have a better explaination of data we need to look at more number of principal components

pca.var <- pca.out$sdev^2
pve <- pca.var/sum(pca.var) *100
png(filename = "percentage_of_variance_explined.png")
plot(cumsum(pve),xlab = 'Principal Component',ylab = 'Percentage of variance explained',ylim=c(0,100),type='b',main='Percentage of variace explained \n as a function of number of principal components')
dev.off()


###Removing neuroendocrine samples from the cleaned gene matrix

samples_to_remove <- rownames(pca_labels@cdesc)[which(grepl("neuroendocrine",pca_labels@cdesc$histological_type_other))]
gene_matrix_without_neuroendocrine<-gene_matrix_cleaned[,-which(colnames(gene_matrix_cleaned) %in% samples_to_remove)]




##Question 3

#Read the file containing the 25 genes defining the IFN signature

setwd("C:/Users/vishal.punvani/Documents/Finalised Assignment_01.01.2020-20200121T152954Z-001/Finalised Assignment_01.01.2020")
IFN_t1 <- read.table("Copy of type1_IFN.txt", header = FALSE, sep = "", dec = ".")
colnames(IFN_t1) <- c("gene")


#GSVA is an unsupervised GSE algorithm 
#1. First The cumulative density function is estimated for every gene using all samples 
#**Since the gene expression is given in the form of integer values therefore we will use poisson distribution function**
#2. Next the GSVA Ranks each gene for all samples 
#3. The rank says how positively or negatively genes are expressed in each sample
#4. Once ranks are calculated then the GSVA scores are calculated using the ranked genes and input gene signature
#5. The GSVA score is calculated by taking maximum deviation between the distribution for gene being present in the IFN signature and gene not being present in the IFN signature in the positive and negative directions and subracting these deviations

setwd("C:/Users/vishal.punvani/Documents/Finalised Assignment_01.01.2020-20200121T152954Z-001/Finalised Assignment_01.01.2020")
gene_set_scores <- GSVA::gsva(as.matrix(gene_matrix_without_neuroendocrine),IFN_t1,kcdf="Poisson")


#Merge sample description file with GSVA scores
#After merging, We have divided samples into three parts using percentile (ntile function) gives us 3 subtypes and created a column named percentile_dev for it 

#**SUBTYPE3**
#Samples having high +ve GSVA(>0.3) score indicates that IFN genes in the gene sample are positively enriched as compared to genes not in the gene set.
#**SUBTYPE1**
#Samples having high -ve GSVA(<-0.3) score indicates that IFN genes in the gene sample are negatively enriched as compared to genes not in the gene set which means sample is deprived of genes in the IFN gene signature 
#**SUBTYPE2**
#Samples having 0.3>0>-0.3 GSVA score indicates that IFN genes in the gene sample are neither positive nor negatively enriched as compared to genes not in the gene set.


paad_file_with_IFN_scores <- merge(paad_file@cdesc,t(gene_set_scores),by = 'row.names')
paad_file_with_IFN_scores$percentile_div <- ntile(paad_file_with_IFN_scores$gene,3)



## Bonus Question 
#Read data containing expression matrix for normal tissue

normal_tissue_data <- read.csv("Copy of Pancreas_log_tpm_RNAseq_mat.csv")
rownames(normal_tissue_data) <- normal_tissue_data$X
normal_tissue_data$X <- NULL

#checking structure of data

head(normal_tissue_data)

#checking number of samples and genes

dim(normal_tissue_data)

###Get the GSVA scores for the normal gene expression dataset
#**Since the gene expression is given in the form of log normalized values therefore we will use Gaussian distribution function**

normal_gene_set_scores <- GSVA::gsva(as.matrix(normal_tissue_data),IFN_t1,kcdf="Gaussian")


plot_dataframe <- as.data.frame(t(gene_set_scores))
plot_dataframe$Type_of_tissue <- "Tumor"
plot_dataframe <- bind_rows(plot_dataframe,as.data.frame(t(normal_gene_set_scores)))
plot_dataframe[is.na(plot_dataframe)] <- "Normal_Tissue"



p <- ggplot(plot_dataframe, aes(x=gene, fill=Type_of_tissue)) +
  geom_density(alpha=0.4) +scale_color_brewer(palette="Paired") + theme_classic() + labs(x="GSVA Score",y="Density") + ggtitle(label = "Comparision of GSVA scores for Normal \n and Tumor Tissues")
p
ggsave(filename = "Density based comparision of Normal Tissue and Cancer Tissue IFN scores.png")

###As seen in the Density plot for GSVA scores of pancreatic noramal tissues and pancreatic cancer tissues we can get the following insights :
#**Insight 1**
#The blue color plot indicating Tumor pancreatic tissues have higher density for large +ve GSVA scores whereas the red plot for normal tissues has low density in this region 
#This is due to the fact that when body is infected with pancreatic tumor, the immune activation system will release IFN proteins resulting in large amount of IFN genes being found in the gene matrix as indicated by the high density of large +ve GSVA score in the plot
#**Insight 2**
#As can be seen from the plot, the density of GSVA score close to 0 is more for normal tissues as compare to Tumor Tissue which indicates that in the normal state the IFN proteins are not released by host cells and therefore are not found in the gene sample resulting in high density of GSVA scores close to 0 for normal tissues as compared to tumor tissues