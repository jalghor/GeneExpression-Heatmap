# Tutorial Referenced: Analyzing data from GEO - Work in Progress. Retrieved from https://sbc.shef.ac.uk/geo_tutorial/tutorial.nb.html#Further_visualisation 
# Code here is following the examples from the tutorial above using a different dataset related to type 2 diabetes


# Install Packages (No need to install if you already have them - if not - install once)

install.packages("ggplot2")
install.packages("readr")
install.packages("BiocManager")
install.packages("ggrepel")
install.packages("tidyr")

BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("org.Hs.eg.db")


# Loading packages 

library(dplyr)
library(GEOquery)
library(pheatmap)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(limma)



Dataset1 <- "GSE25724"  # Assign the ID to the data accession GEO number 
filelist <- getGEO(Dataset1)  # Download the data and store files found in a list 
length(filelist)           # Find how many files there in the list 
filelist <- filelist[[1]]      # Use the first file in the list 
filelist

summary(exprs(filelist)) # To get and print the distribution of the data 

# Value are not normalized so need to perform log2 transformation on the data then check it using a boxplot that should show similar distributions for each sample 
 
exprs(filelist) <- log2(exprs(filelist))
boxplot(exprs(filelist),outline=FALSE, main = 'Normlized Distribution of Data Samples from GSE25724 Dataset', xlab="Samples",
        ylab="Distribution Value") 

# Extracting clinical variables for analysis 

sampleInfo <- pData(filelist)
sampleInfo <- select(sampleInfo, source_name_ch1,characteristics_ch1.1) # Select variables for analysis
sampleInfo <- rename(sampleInfo, cell_type = source_name_ch1, Disease_Status= characteristics_ch1.1) # Rename variables 
#sampleInfo # To display selected columns 

# Calculate correlation between samples to identify variations in the data 
corMatrix <- cor(exprs(filelist), use="c")
corMatrix

## To creat heatmap of correlations and add descriptions to rows and columns in first heatmap 
rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo, main = "Clustered Distrubtion (Distance) of Expressed Genes in Dataset GSE25724") ### to create heatmap   

################################## Further analysis of gene expression 

design <- model.matrix(~0+sampleInfo$cell_type)
colnames(design) <- c("Diabetic","Non_Diabetic") # Rename coulmns 
design

summary(exprs(filelist))

cutoff <- median(exprs(filelist)) # Find the cutoff which is median expression level
is_expressed <- exprs(filelist) > cutoff ##  Finding the expressed genes above the cutoff

keep <- rowSums(is_expressed) > 2 # Keep the genes that are expressed in more than one samples 

table(keep) # Show the number of genes that were above and bolow cutoffs 
filelist <- filelist[keep,] # Save expressed genes only

# To estimate expression levels 
fit <- lmFit(exprs(filelist), design) 
head(fit$coefficients)

# Find contrast between diabetic and non-diabetic groups for differnical analysis  
contrasts <- makeContrasts(Non_Diabetic - Diabetic, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
topTable(fit2)

# To calculate the overall differentially expressed genes
decideTests(fit2)
table(decideTests(fit2))

# To clean and annotated results of identified genes 
anno <- select(anno, Gene = "Gene", ENTREZ_GENE_ID)
fit2$genes <- anno
topTable(fit2)

# Create table with full gene results
full_results <- topTable(fit2, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

# Use to top 20 genes for illustration by filtering genes 

p_cutoff <- 0.05
fc_cutoff <- 1
topN <- 20


ids_of_interest <- mutate(full_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(ID)

gene_names <- mutate(full_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(Gene) 

# Pull rows and coulmns of ids_of_interest 
gene_matrix <- exprs(filelist)[ids_of_interest,]

pheatmap(gene_matrix,
         labels_row = gene_names, scale="row", annotation_col=sampleInfo, main = "Top 20 Differentially Expressed Genes of Diabetic and Non_diabetic Samples in Dataset GSE25724") ### to create heatmap   

