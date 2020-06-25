#Packages ----
#These 3 packages (readr, dplyr and magrittr) are
#part of the tidyverse package
#library(readr)
#library(dplyr) # Manipulate tibbles
#library(magrittr) #Pipelines?
library(tximport) #Importer of Salmon-quantified data
library(DESeq2)
library(tidyverse)
library(dplyr)
library(biomaRt) #API access for gene annotation
library(gridExtra) #Multiple ggplots in one figure
library(plotly) #External plotter
library(pheatmap) #Pretty heatmaps
library(RColorBrewer) #Pretty colors
library(clusterProfiler) #gene ontology stuff
library(org.Hs.eg.db) #Annotation DB for Homo sapiens + ensembl gene
library(grid)
library(lattice)
library(cowplot)

#Set working directory, where the data is. This is purely for example
#setwd("e:/temp/")
setwd(getwd())

#The readr library contains all the read_ functions
#Tab complete will also work with filenames
#read_csv will not produce a data.frame but a tibble
#Pipes are %>%
#View is a function of Rstudio to visualize tables

#Open the data ----
sraruntable = read_csv(file = "SraRunTable.txt")
#read_csv produces a tibble

#Show column and row names:
colnames(sraruntable)
rownames(sraruntable)

#Select the columns (fields) that you want:
#Select() is a function for extracting data out of CSV files because it is
#(again) some weird special-format data table instead of a regular df.
sample_table = dplyr::select(sraruntable, `Sample Name`,
                      source_name, treatment,
                      Cell_Line, Cell_type, time_point)

#sample_table contains 4 rows of each sample, make unique:
unique_sample_table = unique(sample_table)
#alt_unique_sample_table = slice(sample_table, seq(1,48,by=4))

#Extract the sample IDs from the table
sample_ids = as.vector(unique_sample_table[,1])
#as.vector does not work on a tibble
#Alternatively, pull data out with pull
sample_ids = dplyr::pull(.data=unique_sample_table,var = `Sample Name`)

#The paths are the same as the IDs + /quant.sf
sample_paths = paste0(sample_ids, '/quant.sf')

#Now need to get the tx2gene table
#There are no column names in it, add them
tx2gene_data = read_csv(file = "tx2gene.txt", col_names = c("Ensembl Transcripts",
                                                            "GEnsemble Genes"))

#Import the reads from the salmon-analysed files into count_data
count_data = tximport::tximport(files = sample_paths ,
         type = "salmon",
         tx2gene = tx2gene_data ,
         ignoreTxVersion = TRUE)

#Now contains the following sets of information:
#count_data$abundance[1:6,]
#count_data$counts[1:6,]
#count_data$length[1:6,]
#count_data$countsFromAbundance

#Use unique_data_table, which is still a tibble
#Convert tibble to data.frame
unique_sample_table = as.data.frame(unique_sample_table)
#Change name of first column to 'Sample'
colnames(unique_sample_table)[1] = 'Sample'
#For DESeq, we need a list of conditions
#Add those conditions to the unique_sample_table
unique_sample_table$conditions = factor(rep(c('mock_nhbe', 
      'infected_nhbe', 'mock_a549', 'infected_a549'), each = 3))
#Data.frames can contain factor columns
is.factor(unique_sample_table$conditions)
#You can convert a column back into a vector by reading it:
conditions = unique_sample_table$conditions

#Import the data from tximport into deseq
deseq_data = DESeqDataSetFromTximport(txi = count_data,
                                      design =~ conditions,
                                      colData = unique_sample_table)

#From colData, these variables are also now stored in deseq_data:
base::unique(deseq_data$Sample)
base::unique(deseq_data$source_name)
base::unique(deseq_data$treatment)
base::unique(deseq_data$Cell_Line)
base::unique(deseq_data$Cell_type)
base::unique(deseq_data$time_point)
base::unique(deseq_data$conditions)

#To extract read counts from deseq_data, use counts()
counts(deseq_data)[1:6,1:3]
#The DESeq data is still the same as what we imported with tximport
#However, DESeqDataSetFromTximport does turn numbers to integers
count_data$counts[1:6,1:3]

#For normalization of RNA_Seq data, we need to estimate the
#size factors, which ARE THE NORMALIZATION FACTORS
#Notice how the data in deseq_data grows by ~50%
deseq_data = DESeq2::estimateSizeFactors(deseq_data)
DESeq2::sizeFactors(deseq_data)

#Show the normalization factors
DESeq2::normalizationFactors(deseq_data)[1:6,1:3]

#What count now does is multiply the counts with the 
#normalization factor
DESeq2::counts(deseq_data, normalized=TRUE)[1:6, 1:3]

#Which is the same as doing it manually by dividing 
#the data with the normalization factor:
#(why is it divided by, anyway?)
DESeq2::counts(deseq_data, normalized=FALSE)[1:6, 1:3]/DESeq2::normalizationFactors(deseq_data)[1:6,1:3]

#Plot the number of counts of each sample. BIG figure
boxplot(counts(deseq_data, normalized=TRUE))

#Transform the (normalized) count data into one which has a constant
#(approximate) amount of variance across the means (homoskedastic)
#This transformation is similar to a log2 transformation, but is 
#referred to as an 'rlog'
log2(DESeq2::counts(deseq_data, normalized=TRUE)[1:6, 1:3])
vst = DESeq2::varianceStabilizingTransformation(deseq_data)

#We can read the data from this transformed dataset with 'assay'
assay(vst)[1:6,1:3]

#Plot the number of (transformed) counts of each sample, BIG figure
boxplot(assay(vst))

#PCA plot ----
#We can plot the data as a PCA, this uses the top 500 (?) most
#differentially expressed genes (DEGs)
unique(conditions)
base::length(unique(conditions))
plotPCA(vst, intgroup='conditions') +
  theme_bw() + 
  theme(aspect.ratio = 0.75)
#observe the cell line effect

# split experiment in two based on the cell lines nhbe and a540
dds_nhbe = deseq_data[,1:6] #This is nhbe
dds_a549 = deseq_data[,7:12] #This is A549
dds_nhbe = DESeq2::estimateSizeFactors(dds_nhbe)
dds_a549 = DESeq2::estimateSizeFactors(dds_a549)

#The function plotPCA only accepts DESeqTransform objects, which are
#produced by the function rlog or varianceStabilizingTransformation
vst_nhbe = DESeq2::varianceStabilizingTransformation(dds_nhbe)
vst_a549 = DESeq2::varianceStabilizingTransformation(dds_a549)

#Plot a PCA analysis, plotPCA is based on ggplot, so the function
#par(mfrow=) does not work on it. Use gridExtra instead:
p1 = DESeq2::plotPCA(object = vst_nhbe, intgroup='conditions',
                     ntop = 500, returnData = FALSE) + 
  theme_bw() + 
  theme(aspect.ratio = 0.75) +
  ggtitle("nhbe") +
  theme(plot.title = element_text(hjust = 0.5)) #allign center
p2 = DESeq2::plotPCA(object = vst_a549, intgroup='conditions',
                     ntop = 500, returnData = FALSE) + 
  theme_bw() + 
  theme(aspect.ratio = 0.75) +
  ggtitle("a549") +
  theme(plot.title = element_text(hjust = 0.5)) #allign center
gridExtra::grid.arrange(p1, p2, nrow = 1, ncol = 2, top = "Principal component analysis")

#Cluster analysis----
#Cluster analysis using hclust (Hierarchical Clustering)
#hclust expects the SAMPLES ON THE ROWS, AND PARAMS ON COLS
#MEANING THAT YOU MUST TRANSPOSE PRIOR TO CALCULATING THE
#DISTANCE MATRIX
d_nhbe = SummarizedExperiment::assay(vst_nhbe) #Extract normalized + transformed data out
d_a549 = SummarizedExperiment::assay(vst_a549)

#Transpose the data; needed for dist()
d_nhbe = t(d_nhbe)
d_a549 = t(d_a549)

#Calculate the distance matrix
d_nhbe = dist(d_nhbe)
d_a549 = dist(d_a549)

h_nhbe = hclust(d_nhbe)
h_a549 = hclust(d_a549)

#Plot the hclust figures
par(mfrow=c(1,2))
plot(h_nhbe, labels = conditions[1:6])
plot(h_a549, labels = conditions[7:12])
par(mfrow=c(1,1))
remove(list=c("d_nhbe", "d_a549", "h_nhbe", "h_a549"))
gc()

#You can compact this into the following:
#This plot is an unsupervised cluster analysis, meaning it will
#not limit itself to a predefined number of subgroups
#hang defines the length of the lowest "sticks"
par(mfrow=(c(1, 2)))
plot(hclust(dist(t(assay(vst_nhbe)))), hang = 0.1, labels = conditions[1:6])
plot(hclust(dist(t(assay(vst_a549)))), hang = 0.1, labels = conditions[7:12])
par(mfrow=(c(1, 1)))

#For supervised clustering, we specify a number of clusters
#that the data should be fitted to. This can be done using K-
#means clustering
k = kmeans(t(assay(vst_nhbe)), centers=2)
k$cluster

# 3 steps to DESeq2 analysis
# 1) estimate size factors (normalisation)
# 2) estimate dispersions
# 3) apply statistics (Wald Test)
#So far we've only done the first step. 

##Split the data up again ----

#Note that what we did above does not entirely work as the data
#was not split early enough. So import the data again and create
#new variables:
sample_files_nhbe = sample_paths[1:6]
sample_files_a549 = sample_paths[7:12]
sample_table_nhbe = unique(sample_table)[1:6,]
sample_table_a549 = unique(sample_table)[7:12,]

#Alternatively we can use the dplyr function called "filter"
sample_table_nhbe = dplyr::filter(unique(sample_table), Cell_Line=="NHBE")

#Overwrite the conditions column in sample_table and change the factor
#to only contain mock and infected. If column conditions contains
#factors which are not present, then the function estimateDispersions
#will throw errors.
sample_table_nhbe$conditions = factor(rep(c('mock', 'infected'), each=3), 
                                      levels=c('mock', 'infected'))
sample_table_a549$conditions = factor(rep(c('mock', 'infected'), each=3), 
                                      levels=c('mock', 'infected'))

#Import the count data again, all the way from the salmon files:
count_data_nhbe = tximport::tximport(files = sample_files_nhbe,
                                    type = "salmon",
                                    tx2gene = tx2gene_data,
                                    ignoreTxVersion = TRUE)
count_data_a549 = tximport::tximport(files = sample_files_a549,
                                    type = "salmon",
                                    tx2gene = tx2gene_data,
                                    ignoreTxVersion = TRUE)

#Import the txi data into DESeq
dds_nhbe = DESeq2::DESeqDataSetFromTximport(txi=count_data_nhbe, 
                                    colData=sample_table_nhbe,
                                    design=~conditions)
dds_a549 = DESeq2::DESeqDataSetFromTximport(txi=count_data_a549, 
                                    colData=sample_table_a549,
                                    design=~conditions)

#Now calculate the size factors (normalization)
dds_nhbe = DESeq2::estimateSizeFactors(dds_nhbe)
dds_a549 = DESeq2::estimateSizeFactors(dds_a549)

#Calculate the dispersion factors and add them to dds_
dds_nhbe = DESeq2::estimateDispersions(dds_nhbe)
dds_a549 = DESeq2::estimateDispersions(dds_a549)

#The dispersion estimates can be plotted using plotDispEsts
par(mfrow=c(1,2))
DESeq2::plotDispEsts(object = dds_nhbe, legend = TRUE, 
                     xlab = 'mean of normalized counts (nhbe)')
DESeq2::plotDispEsts(object = dds_a549, legend = TRUE, 
                     xlab = 'mean of normalized counts (a549)')
par(mfrow=c(1,1))

#Perform the Wald tests for both the nhbe and the a549 dataset
dds_nhbe = nbinomWaldTest(dds_nhbe); dds_a549 = nbinomWaldTest(dds_a549);

#The function estimateSizeFactors, estimateDispersions, and Wald test
#is summarized in the function DESeq()
#dds_nhbe = DESeq(dds_nhbe)
#dds_a549 = DESeq(dds_a549)

#Simple results can be obtained per gene, with plotCounts:
DESeq2::plotCounts(dds_nhbe, gene='ENSG00000265794',
           intgroup='conditions')

#The results can be extracted from DESeq2 with results:
results_df_nhbe = as.data.frame(DESeq2::results(dds_nhbe))
results_df_a549 = as.data.frame(DESeq2::results(dds_a549))

#Alternatively results can be extracted from the DESeq table by ways of 'contrasting'
#a group to another group. The parameter 'contrast' needs 3 arguments: the condition 
#that you're comparing, and a/b ie infected/control
#Remember that in the experimental design we set up experimental conditions 'mock' and
#the 'infected'. We can use these factor levels
dds_nhbe$conditions
DESeq2::results(dds_nhbe, contrast = c('conditions', 'infected', 'mock'))
#You can see the total number of genes with a 'nonzero' read count in the summary:
summary(DESeq2::results(dds_nhbe, contrast = c('conditions', 'infected', 'mock')))
#Keep in mind that results() does not filter the output. You will get the same number
#of rows as what was present in the dds dataset to begin with, but now it will have
#some values which are NA, as they fell below the cutoffs. 

#If we request the sum of all complete cases, that is,
#cases which are non-zero and which exceed the "low counts"
#filter, we get, out of the 60233 reads, 28451 non-zero ones,
#and 15983 results which are complete.
sum(complete.cases(results_df_nhbe))
sum(complete.cases(results_df_a549))
#which is the same as listing the dimensions of a filtered table:
dim(results_df_nhbe[complete.cases(results_df_nhbe),])
dim(results_df_a549[complete.cases(results_df_a549),])
summary(DESeq2::results(dds_nhbe))
summary(DESeq2::results(dds_a549))

#Now we need to filter the data for valid and interesting cases
filter_df1_nhbe = results_df_nhbe[complete.cases(results_df_nhbe),]
filter_df2_nhbe = filter_df1_nhbe[filter_df1_nhbe$padj < 0.05,]
filter_df3_nhbe = filter_df2_nhbe[abs(filter_df2_nhbe$log2FoldChange)>=1,]

#Also filter the A549 data
filter_df1_a549 = results_df_a549[complete.cases(results_df_a549),]
filter_df2_a549 = filter_df1_a549[filter_df1_a549$padj < 0.05,]
filter_df3_a549 = filter_df2_a549[abs(filter_df2_a549$log2FoldChange)>=1,]

par(mfrow=c(1, 2))
DESeq2::plotMA(DESeq2::results(dds_nhbe), alpha = 0.05,
               main = "MA plot of nhbe, red is p < 0.05", 
               xlab = "mean of normalized counts (nhbe)")
DESeq2::plotMA(DESeq2::results(dds_a549), alpha = 0.05,
               main = "MA plot of a549, red is p < 0.05", 
               xlab = "mean of normalized counts (a549)")
par(mfrow=c(1, 1))

##Volcano plots ----
#We need to add a column to filter_df1 which is a TRUE or FALSE based on
#whether the row has a padj < 0.05 and a fc greater than 1 
filter_df1_nhbe$test = filter_df1_nhbe$padj < 0.05 & abs(filter_df1_nhbe$log2FoldChange) > 1
filter_df1_a549$test = filter_df1_a549$padj < 0.05 & abs(filter_df1_a549$log2FoldChange) > 1

#Convert rownames to a column. var is the name of the new column
if (tibble::has_rownames(filter_df1_nhbe)) {
  filter_df1_nhbe = tibble::rownames_to_column(filter_df1_nhbe, var='ensgene')
}
if (tibble::has_rownames(filter_df1_a549)) {
  filter_df1_a549 = tibble::rownames_to_column(filter_df1_a549, var='ensgene')
}

p1 = ggplot(filter_df1_nhbe, aes(x=log2FoldChange, 
                           y=-log10(padj), 
                           name=ensgene)) +
  geom_point(aes(colour=test), size=1, alpha=0.3) +
  scale_colour_manual(values=c('black', 'red')) +
  geom_vline(xintercept=1, colour='green', linetype=3) +
  geom_vline(xintercept=-1, colour='green', linetype=3) +
  geom_hline(yintercept=-log10(0.05), colour='blue', linetype=3) +
  theme_bw() +
  theme(legend.position = 'none')

p2 = ggplot(filter_df1_a549, aes(x=log2FoldChange, 
                           y=-log10(padj), 
                           name=ensgene)) +
  geom_point(aes(colour=test), size=1, alpha=0.3) +
  scale_colour_manual(values=c('black', 'red')) +
  geom_vline(xintercept=1, colour='green', linetype=3) +
  geom_vline(xintercept=-1, colour='green', linetype=3) +
  geom_hline(yintercept=-log10(0.05), colour='blue', linetype=3) +
  theme_bw() +
  theme(legend.position = 'none')
#Output both plots in one figure:
gridExtra::grid.arrange(p1, p2, nrow = 1)

#We can use plotly, which is an interactive plot viewer, so we can
#hoover over points to see what they are
plotly::ggplotly(p1)

##Biomart part ----
biomaRt::listMarts()
#View(listDatasets(ensembl99))

#Load in the homosapiens ensembl version 99 mart
ensembl99 = biomaRt::useEnsembl(biomart="ensembl", version=99)
ensembl99 = biomaRt::useDataset("hsapiens_gene_ensembl", mart=ensembl99)

#The biomaRt api has access to the following long list
#of attributes, including ensembl_gene_id and description:
#View(biomaRt::listAttributes(ensembl99))
#And you can use the folllowing filters on datasets:
#View(listFilters(ensembl99))

#Data can be obtained from ensembl using the function getBM()
#If you are requesting data which contains several lines of data
#per gene, then this will be displayed on several lines, which 
#also results in (some) duplicate data.
biomaRt::getBM(attributes=c('ensembl_gene_id', 'ensembl_gene_id_version',
                   'ensembl_transcript_id', 'ensembl_transcript_id_version',
                   'external_gene_name'), 
      filters = c('ensembl_gene_id'), 
      values = filter_df1_nhbe$ensgene[1:6],
      mart = ensembl99)

#Now to request the data for all complete data in the nhbe dataset
#we will use the list of ensgenes in filter_df1_nhbe to make the request
annotation_nhbe = biomaRt::getBM(attributes=c('ensembl_gene_id',
                                'chromosome_name',
                                'start_position',
                                'end_position',
                                'strand',
                                'gene_biotype',
                                'external_gene_name',
                                'description'),
                   filters = c('ensembl_gene_id'),
                   values = filter_df1_nhbe$ensgene,
                   mart = ensembl99)

annotation_a549 = biomaRt::getBM(attributes=c('ensembl_gene_id',
                                              'chromosome_name',
                                              'start_position',
                                              'end_position',
                                              'strand',
                                              'gene_biotype',
                                              'external_gene_name',
                                              'description'),
                                 filters = c('ensembl_gene_id'),
                                 values = filter_df1_a549$ensgene,
                                 mart = ensembl99)

#We have now gathered the data, we can compare the obtained data
#against the original query using dims:
dim(filter_df1_nhbe)
dim(annotation_nhbe)
if (length(filter_df1_nhbe$ensgene) == length(annotation_nhbe$ensembl_gene_id)) {
  print("Input and output of biomaRt query were the same")
} else {
  print("biomaRt query not the same!!")
}

#We can add the obtained data to an existing dataframe (such as filter_df1_nbhe)
#for including such data into ggplots
annotated_filter_df1_nhbe = left_join(filter_df1_nhbe, annotation_nhbe,
                         by=c('ensgene'='ensembl_gene_id'))
annotated_filter_df1_a549 = left_join(filter_df1_a549, annotation_a549,
                          by=c('ensgene'='ensembl_gene_id'))

#Make a ggplot for plotly with the external_gene_name in it ---- 
#Now feed everything into a ggplot
p1 = ggplot(annotated_filter_df1_nhbe, aes(x=log2FoldChange, 
                             y=-log10(padj), 
                             name=external_gene_name,
                             desc=description)) +
  geom_point(aes(colour=test), size=1, alpha=0.3) +
  scale_colour_manual(values=c('black', 'red')) +
  geom_vline(xintercept=1, colour='green', linetype=3) +
  geom_vline(xintercept=-1, colour='green', linetype=3) +
  geom_hline(yintercept=-log10(0.05), colour='blue', linetype=3) +
  theme_bw() +
  theme(legend.position = 'none')

p2 = ggplot(annotated_filter_df1_a549, aes(x=log2FoldChange, 
                                           y=-log10(padj), 
                                           name=external_gene_name,
                                           desc=description)) +
  geom_point(aes(colour=test), size=1, alpha=0.3) +
  scale_colour_manual(values=c('black', 'red')) +
  geom_vline(xintercept=1, colour='green', linetype=3) +
  geom_vline(xintercept=-1, colour='green', linetype=3) +
  geom_hline(yintercept=-log10(0.05), colour='blue', linetype=3) +
  theme_bw() +
  theme(legend.position = 'none')
gridExtra::grid.arrange(p1, p2, nrow = 1)

#We can now view the graph using ggplotly. When hoovering over a point
#it will display the gene name and its description, as we added this data
#to the ggplot
ggplotly(p1)

##Methods of selecting data ----
#We need the differentially expressed genes, and the data filtered
#using varianceStabilizingTransformation. The function assay() pulls
#data out of the type "DESeqTransform and exports it into a simple
#matrix. We can test the data type with the function class().
vst_nhbe_matrix = SummarizedExperiment::assay(vst_nhbe)
vst_a549_matrix = SummarizedExperiment::assay(vst_a549)
class(vst_nhbe) #vst_nhbe is produced from dds_nhbe, which itself is made from the imported count data 
class(vst_a549)
class(vst_nhbe_matrix) #These are the matrices produced by the assay() function
class(vst_a549_matrix)

#Right now the vst_nhbe_matrix (and the a549 one) contain the full dataset, so that
#contains ~60k different genes, each with (or indeed without) their individual count
#data per sample. 
dim(vst_nhbe_matrix);
dim(vst_a549_matrix);
vst_nhbe_matrix = vst_nhbe_matrix[complete.cases(results_df_nhbe),]
vst_a549_matrix = vst_a549_matrix[complete.cases(results_df_a549),]

#If we run the dim() function again, we'll see that the total number of rows
#has decreased from 60,233 rows to 15,983 for NHBE and 14,995 for A549.
#THESE ARE THE COMPLETE CASES
dim(vst_nhbe_matrix);
dim(vst_a549_matrix);

#We can read and write to colnames and rownames using the same function
#but in different order. Here we write the name of the samples to colnames
rownames(vst_nhbe_matrix) = annotation_nhbe$external_gene_name
rownames(vst_a549_matrix) = annotation_a549$external_gene_name
colnames(vst_nhbe_matrix) = conditions[1:6]
colnames(vst_a549_matrix) = conditions[7:12]

#To specifically select the genes which are differentially expressed and have a fold-change
#greater than -2 or +2, we select for them:
vst_nhbe_matrix = vst_nhbe_matrix[filter_df1_nhbe$padj < 0.05 & abs(filter_df1_nhbe$log2FoldChange) > 1,]
vst_a549_matrix = vst_a549_matrix[filter_df1_a549$padj < 0.05 & abs(filter_df1_a549$log2FoldChange) > 1,]

#Let's see how many of them are left now:
dim(vst_nhbe_matrix); dim(vst_a549_matrix);

#Alternatively we could have done this in another way, by selecting the ensembl geneIDs out
#of the original 60k row tables
vst_nhbe_matrix = SummarizedExperiment::assay(vst_nhbe)
vst_a549_matrix = SummarizedExperiment::assay(vst_a549)
dim(vst_nhbe_matrix); dim(vst_a549_matrix);
#The differentially expressed genes are as follows:
row.names(results_df_nhbe[complete.cases(results_df_nhbe)
                          & results_df_nhbe$padj < 0.05
                          & abs(results_df_nhbe$log2FoldChange) > 1,])
row.names(results_df_a549[complete.cases(results_df_a549)
                          & results_df_a549$padj < 0.05
                          & abs(results_df_a549$log2FoldChange) > 1,])
#The ensembl gene IDs inside are vst_nhbe_matrix are stored as rownames, thus we can filter for them
vst_nhbe_matrix = vst_nhbe_matrix[complete.cases(results_df_nhbe)
                    & results_df_nhbe$padj < 0.05
                    & abs(results_df_nhbe$log2FoldChange) > 1,]
vst_a549_matrix = vst_a549_matrix[complete.cases(results_df_a549)
                    & results_df_a549$padj < 0.05
                    & abs(results_df_a549$log2FoldChange) > 1,]
dim(vst_nhbe_matrix); dim(vst_a549_matrix)

#We need to add the gene names back to the rownames.
#The data.frame annotation_nhbe and annotation_a549 contain the ensembl IDs
#and the external_gene_name of all the ones in our data. In this case we only want
#the gene names of those which we selected out of the data as being p<0.05 and
#with log2fc greater than 1. 
row.names(vst_nhbe_matrix)[1:3]
row.names(vst_a549_matrix)[1:3]
row.names(vst_nhbe_matrix) <- annotation_nhbe$external_gene_name[annotation_nhbe$ensembl_gene_id
                                                                 %in% row.names(vst_nhbe_matrix)]
row.names(vst_a549_matrix) <- annotation_a549$external_gene_name[annotation_a549$ensembl_gene_id
                                                                 %in% row.names(vst_a549_matrix)]

#I'd also like to have sample names in my heatmaps.
colnames(vst_nhbe_matrix) <- paste(1:6, conditions[1:6]);
colnames(vst_a549_matrix) <- paste(7:12, conditions[7:12]);

#Heatmap ----
#Produce a heatmap of the differentially expressed genes:
#The heatmap requires a list (here we use VST-corrected gene expression)
#With the gene names (or at least ID numbers) on the rownames and the sample
#names on the colnames.
#VST correction was done to transform the dataset into one which is approximately
#homoskedastic. Meaning that the relative variance is equal between low and high counts
row.names(vst_nhbe_matrix)[1:3]
row.names(vst_a549_matrix)[1:3]
colnames(vst_nhbe_matrix)[1:6]
colnames(vst_a549_matrix)[1:6]

#Now produce the actual heatmap, first let's try it out with the base
#heatmap function
stats::heatmap(vst_nhbe_matrix)
stats::heatmap(vst_a549_matrix)

#Produce a nicer heatmap using the pheatmap package
pheatmap::pheatmap(vst_nhbe_matrix, fontsize_row=4, scale='row')
pheatmap::pheatmap(vst_a549_matrix, fontsize_row=4, scale='row')
#cowplot::plot_grid(p1, p2, ncol = 2, labels = c('A', 'B'))
#The package heatmap does not permit writing its output to a variable
#thus I can't use multiplexing functions such as cowplot::plotgrid

#Heatmap using pheatmap with different colours. We will use the package
#RColorBrewer as a source for the colors. In this package there is also
#a function called colorRampPalette which smoothes out the palette
RColorBrewer::display.brewer.all();

#Produce the heatmap in black&white
RColorBrewer::display.brewer.pal(n = 9, name = "Greys")
pheatmap::pheatmap(vst_nhbe_matrix, fontsize_row=4, scale='row',
         color=colorRampPalette(brewer.pal(n = 9, name = "Greys"))(100))

#Use a stupid set of colors called "paried", just why would anyone use this
#for continues data? Maybe useful for categorical data? Who knows
RColorBrewer::display.brewer.pal(n = 12, name = "Paired")
pheatmap::pheatmap(vst_nhbe_matrix, fontsize_row=4, scale='row',
         color=colorRampPalette(brewer.pal(n = 12, name = "Paired"))(100)) #just why

#Add additional elements to the heatmap, such as dendrogram tree-cuts
RColorBrewer::display.brewer.pal(n = 7, name = "Blues")
#We can interpolate the colours using the function colorRampPalette
#At the end of this function, we tell it how many values we want to
#have returned to us (in a vector)
grDevices::colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)

#construct the pretty pheatmap with (automatically determined) treecuts 
pheatmap(mat = vst_nhbe_matrix,
         kmeans_k = NA,
         border_color = "#000000",
         fontsize_row = 4, 
         scale = 'row',
         color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100), 
         cutree_cols = 2,
         cutree_rows = 2,
         angle_col = 45,
         main = "Heatmap of nhbe DEGs",
         legend = TRUE)

#Let's clean up the working memory:
gc();

#GO enrichment analysis ----
#For the GO analysis, it's better to use entrez gene IDs instead of ensembl gene IDs
#Load in the homosapiens ensembl version 99 mart
ensembl99 = biomaRt::useEnsembl(biomart="ensembl", version=99)
ensembl99 = biomaRt::useDataset("hsapiens_gene_ensembl", mart=ensembl99)

#We can get those using the biomart api. The code that I'm putting after values
#represents a vector of 141 values for NHBE and 58 values for A549 which were filtered
#from the original data for complete.cases, adjusted P value and a fold-change > 1.
#These 141 & 58 ensembl IDs are in essence our "genes of interest.
entrez_nhbe = biomaRt::getBM(attributes=c('ensembl_gene_id',
                                          'entrezgene_id'),
                             filters = c('ensembl_gene_id'),
                             values = row.names(results_df_nhbe[complete.cases(results_df_nhbe) & results_df_nhbe$padj < 0.05 & abs(results_df_nhbe$log2FoldChange) > 1,]),
                             mart = ensembl99)
entrez_a549 = biomaRt::getBM(attributes=c('ensembl_gene_id',
                                          'entrezgene_id'),
                             filters = c('ensembl_gene_id'),
                             values = row.names(results_df_a549[complete.cases(results_df_a549) & results_df_a549$padj < 0.05 & abs(results_df_a549$log2FoldChange) > 1,]),
                             mart = ensembl99)
#However not every ensembl ID has an entrez ID. Check if this is the case:
length(entrez_nhbe$ensembl_gene_id[!is.na(entrez_nhbe$ensembl_gene_id)]);
length(entrez_nhbe$entrezgene_id[!is.na(entrez_nhbe$entrezgene_id)]);
print('the difference is: '); 
length(entrez_nhbe$ensembl_gene_id[!is.na(entrez_nhbe$ensembl_gene_id)]) - length(entrez_nhbe$entrezgene_id[!is.na(entrez_nhbe$entrezgene_id)]);
#And for A549:
length(entrez_a549$ensembl_gene_id[!is.na(entrez_a549$ensembl_gene_id)]);
length(entrez_a549$entrezgene_id[!is.na(entrez_a549$entrezgene_id)]);
print('the difference is: '); 
length(entrez_a549$ensembl_gene_id[!is.na(entrez_a549$ensembl_gene_id)]) - length(entrez_a549$entrezgene_id[!is.na(entrez_a549$entrezgene_id)]);
#We can get the list of valid entrez IDs as follows:
entrez_nhbe$entrezgene_id[!is.na(entrez_nhbe$entrezgene_id)][1:3];
#And for A549:
entrez_a549$entrezgene_id[!is.na(entrez_a549$entrezgene_id)][1:3];

#For the GO enrichment analysis, we need to have a list of "universe" genes which
#represent all the genes in the pool. In this case it would be fair to consider 
#that the tested genes, which are the ones after the complete.cases() function
#represent all "universe" genes.

entrez_nhbe_uni = biomaRt::getBM(attributes=c('ensembl_gene_id',
                                              'entrezgene_id'),
                                 filters = c('ensembl_gene_id'),
                                 values = row.names(results_df_nhbe[complete.cases(results_df_nhbe),]),
                                 mart = ensembl99)
entrez_a549_uni = biomaRt::getBM(attributes=c('ensembl_gene_id',
                                              'entrezgene_id'),
                                 filters = c('ensembl_gene_id'),
                                 values = row.names(results_df_a549[complete.cases(results_df_a549),]),
                                 mart = ensembl99)
#Check the lengths and see how many ensembl IDs do not have an entrez ID:
length(entrez_nhbe_uni$ensembl_gene_id[!is.na(entrez_nhbe_uni$ensembl_gene_id)]); length(entrez_nhbe_uni$entrezgene_id[!is.na(entrez_nhbe_uni$entrezgene_id)]);
#The differfnce between the number of ensembl IDs and the number of entrez IDs is:
length(entrez_nhbe_uni$ensembl_gene_id[!is.na(entrez_nhbe_uni$ensembl_gene_id)]) - length(entrez_nhbe_uni$entrezgene_id[!is.na(entrez_nhbe_uni$entrezgene_id)]);

#Also do this for the data of A549:
length(entrez_a549_uni$ensembl_gene_id[!is.na(entrez_a549_uni$ensembl_gene_id)]); length(entrez_a549_uni$entrezgene_id[!is.na(entrez_a549_uni$entrezgene_id)]);
#The differfnce between the number of ensembl IDs and the number of entrez IDs is:
length(entrez_a549_uni$ensembl_gene_id[!is.na(entrez_a549_uni$ensembl_gene_id)]) - length(entrez_a549_uni$entrezgene_id[!is.na(entrez_a549_uni$entrezgene_id)]);

#Now to actually do the GO analysis. ont stands for "ontology", the type of process
#that you wish to test for. BP for biological process, MF for molecular function
ego_nhbe = clusterProfiler::enrichGO(gene = entrez_nhbe$entrezgene_id[!is.na(entrez_nhbe$entrezgene_id)],
                                OrgDb = org.Hs.eg.db,
                                ont = "BP",
                                universe = as.character(entrez_nhbe_uni$entrezgene_id))

ego_a549 = clusterProfiler::enrichGO(gene = entrez_a549$entrezgene_id[!is.na(entrez_a549$entrezgene_id)],
                                     OrgDb = org.Hs.eg.db,
                                     ont = "BP",
                                     universe = as.character(entrez_a549_uni$entrezgene_id))
#The result table can be shown using summary()
#summary(ego_nhbe);
#summary(ego_a549);
#The output says that summary() method is deprecated and suggests to use as.data.frame instead
#so sure we can do that with the data:
as.data.frame(ego_nhbe)
as.data.frame(ego_a549)
#You can use view() to output it (in Rstudio) into more conveniently-visualized tables
#view(as.data.frame(ego_nhbe))
#view(as.data.frame(ego_a549))
#In the output table, we are most intrested in the column "qvalue", which represents 
#the final adjusted p-value. (Not sure what Padj means in this table)

#Visualize Gene Ontology results ----
p1 <- graphics::barplot(ego_nhbe, showCategory = 20)
p2 <- graphics::barplot(ego_a549, showCategory = 20)
#We can plot two graphs into one figure
gridExtra::grid.arrange(p1, p2, ncol=2)

#Visualize with dotplots. This function dotplot uses ggplot as its base
p1 <- clusterProfiler::dotplot(ego_nhbe, showCategory = 15)
p2 <- clusterProfiler::dotplot(ego_a549, showCategory = 15)
gridExtra::grid.arrange(p1, p2, ncol=2)

#cnetplot: connectivitynet plot
p1 <- clusterProfiler::cnetplot(ego_nhbe, circular = TRUE, colorEdge = TRUE)
p2 <- clusterProfiler::cnetplot(ego_a549, circular = TRUE, colorEdge = TRUE)
gridExtra::grid.arrange(p1, ncol=1)
gridExtra::grid.arrange(p2, ncol=1)
cowplot::plot_grid(p1, p2, ncol = 2, labels = c('A', 'B'))

p1 <- clusterProfiler::goplot(ego_nhbe)
p2 <- clusterProfiler::goplot(ego_a549)
cowplot::plot_grid(p1, p2, ncol = 2, labels = c('A', 'B'))

#Gene ontology analysis with KEGG ----
#This function actually relies on an input of entrez gene IDs and will not accept
#ensembl IDs at all.
kegg_nhbe = clusterProfiler::enrichKEGG(gene = entrez_nhbe$entrezgene_id[!is.na(entrez_nhbe$entrezgene_id)],
                                        organism = 'hsa',
                                        pvalueCutoff = 0.05,
                                        universe = as.character(entrez_nhbe_uni$entrezgene_id))

kegg_a549 = clusterProfiler::enrichKEGG(gene = entrez_a549$entrezgene_id[!is.na(entrez_a549$entrezgene_id)],
                                        organism = 'hsa',
                                        pvalueCutoff = 0.05,
                                        universe = as.character(entrez_a549_uni$entrezgene_id))
#As with enrichGO, the function summary() is deprecated and should not be used
#summary(kegg_nhbe)
#summary(kegg_a549)
as.data.frame(kegg_nhbe)
as.data.frame((kegg_a549))


#Writing stuff out to a table ----
getwd()
readr::write_tsv(x = annotated_filter_df1_nhbe[annotated_filter_df1_nhbe$test==TRUE,],
                 path = "NHBE_differntially_ex_genes.txt",
                 na = "NA",
                 append = FALSE,
                 col_names = TRUE)
#When importing into excel, make sure that the table file is imported with gene_names as
#string (text) and not as general, as certain gene names will get converted into dates
#by excel.

#The current workspace can be saved with save.image(). Note that the parameter
#"compression_level" doesn't work in this function. 
save.image(file = ".RData", compress = "bzip2")
#And we can save individual objects with save()
save(list = "annotated_filter_df1_nhbe", file = "annotated_filter_df1_nhbe",
     compress = "bzip2", compression_level = 9, precheck = TRUE)


#Manual PCA plot----
#Calculate a PCA plot manually using prcomp(), which is an internal function of R
#stats module.
#prcomp() calculates the principal components accross the rows and not across the
#columns. So we (optionally) need to transpose the data so that the samples are on
#the rows and the parameters measured are the columns
#We'll use the variance stabilized data (the VST)
t(SummarizedExperiment::assay(vst_nhbe))[,1:3]
dim(t(SummarizedExperiment::assay(vst_nhbe)))
pca_nhbe = prcomp(t(SummarizedExperiment::assay(vst_nhbe)))
pca_a549 = prcomp(t(SummarizedExperiment::assay(vst_a549)))
#The number of principal components is dependent on the number of samples?
#Set some rownames, to sample names
rownames(pca_nhbe$x) = sample_table_nhbe$`Sample Name`
rownames(pca_a549$x) = sample_table_a549$`Sample Name`

#Let's turn it into a figure:
#We can add formulas inside the ggplot to dynamically add data to (for instance) the axes
#Round up numbers with round(number, digits)
#The percent variance explained in a PCA is calculated by sdev^2/sum(sdev^2)*100
p1 <- ggplot2::ggplot(data = as.data.frame(pca_nhbe$x), aes(x = PC1, y = PC2, colour=factor(sample_table_nhbe$conditions, labels = c("\nmock\ninfected\n","\nSarscov2\ninfected\n")))) +
  geom_point(size = 2) +
  scale_colour_discrete("Conditions") +
  theme_bw() +
  xlab(label = paste0("PC1 (", round((pca_nhbe$sdev^2/sum(pca_nhbe$sdev^2)*100)[1],0), "%)")) +
  ylab(label = paste0("PC2 (", round((pca_nhbe$sdev^2/sum(pca_nhbe$sdev^2)*100)[2],0), "%)"))
p2 <- ggplot2::ggplot(data = as.data.frame(pca_a549$x), aes(x = PC1, y = PC2, colour=factor(sample_table_a549$conditions, labels = c("\nmock\ninfected\n","\nSarscov2\ninfected\n")))) + 
  geom_point(size = 2) +
  scale_colour_discrete("Conditions") +
  theme_bw() +
  xlab(label = paste0("PC1 (", round((pca_a549$sdev^2/sum(pca_a549$sdev^2)*100)[1],0), "%)")) +
  ylab(label = paste0("PC2 (", round((pca_a549$sdev^2/sum(pca_a549$sdev^2)*100)[2],0), "%)"))
cowplot::plot_grid(p1, p2, ncol = 2, labels = c('A', 'B'))

#Now let's try to do the same but specifically for the top500 most differentially expressed genes.
#We can probably do this with a sorted filter_df1_nhbe table. We can sort columns with the
#dplyr function "arrange"
#The following function will extract the top 500 VST data:
assay(vst_nhbe)[dplyr::arrange(.data = filter_df1_nhbe, desc(abs(log2FoldChange)))[1:500,]$ensgene,]
#And we can get the ensembl IDs from it easily:
dplyr::arrange(.data = filter_df1_nhbe, desc(abs(log2FoldChange)))[1:500,]$ensgene

pca_nhbe_top500 = prcomp(t(SummarizedExperiment::assay(vst_nhbe)[dplyr::arrange(.data = filter_df1_nhbe, desc(abs(log2FoldChange)))[1:500,]$ensgene,]))
pca_a549_top500 = prcomp(t(SummarizedExperiment::assay(vst_a549)[dplyr::arrange(.data = filter_df1_a549, desc(abs(log2FoldChange)))[1:500,]$ensgene,]))

#add rownames. This isn't needed for the ggplot per se. 
rownames(pca_nhbe_top500$x) = sample_table_nhbe$`Sample Name`
rownames(pca_a549_top500$x) = sample_table_a549$`Sample Name`

#Plot the top500 most differentially expressed genes PCA on an XY scatter
#The axes labels will include the percentage variance explained of the 
#principal components 1 and 2.
p1 <- ggplot2::ggplot(data = as.data.frame(pca_nhbe_top500$x), aes(x = PC1, y = PC2, colour=factor(sample_table_nhbe$conditions, labels = c("\nmock\ninfected\n","\nSarscov2\ninfected\n")))) +
  geom_point(size = 2) +
  scale_colour_discrete("Conditions") +
  theme_bw() +
  xlab(label = paste0("PC1 (", round((pca_nhbe_top500$sdev^2/sum(pca_nhbe_top500$sdev^2)*100)[1],0), "%)")) +
  ylab(label = paste0("PC2 (", round((pca_nhbe_top500$sdev^2/sum(pca_nhbe_top500$sdev^2)*100)[2],0), "%)"))
p2 <- ggplot2::ggplot(data = as.data.frame(pca_a549_top500$x), aes(x = PC1, y = PC2, colour=factor(sample_table_a549$conditions, labels = c("\nmock\ninfected\n","\nSarscov2\ninfected\n")))) + 
  geom_point(size = 2) +
  scale_colour_discrete("Conditions") +
  theme_bw() +
  xlab(label = paste0("PC1 (", round((pca_a549_top500$sdev^2/sum(pca_a549_top500$sdev^2)*100)[1],0), "%)")) +
  ylab(label = paste0("PC2 (", round((pca_a549_top500$sdev^2/sum(pca_a549_top500$sdev^2)*100)[2],0), "%)"))
cowplot::plot_grid(p1, p2, ncol = 2, labels = c('A', 'B'))















