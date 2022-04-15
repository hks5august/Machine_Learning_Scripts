#Load packages
library(Rtsne)
library(ggplot2)

#Set seed
set.seed(2)

#load data
data<- read.table('https://raw.githubusercontent.com/PineBiotech/bioinformatics/master/15gene_transposed.txt', sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=TRUE, row.names=1)

#Take data without groups
df <- data[2:16]

#Create t-SNE object
tsne_realData <- Rtsne(df, perplexity=15, verbose=FALSE, max_iter = 550, check_duplicates = FALSE)

#Add colors
colors = rainbow(length(unique(data$Group)))

#Make tSNE object as data frame
tsne_df <- as.data.frame(tsne_realData$Y)

#Extract labels of samples
lab <- data$Group

#Cobine samples with t-SNE plot
tsne_df1 <- cbind(tsne_df,lab)

#Add column names
colnames(tsne_df1) <- c("X","Y","Labels")

#Create colored plot with labels of sample groups
plot<- ggplot(tsne_df1, aes(X, Y, colour = Labels)) +geom_point()

#Check dimensions of t-SNE input data
dimensions <- dim(tsne_df1)

#Print dimensions of t-SNE object dataframe
print(dimensions)

#Show plot
plot
