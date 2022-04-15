#Load libraries
library(ggfortify)
library(cluster)

#Load data
ExpressionTable <- read.table('https://raw.githubusercontent.com/PineBiotech/bioinformatics/master/15gene_transposed.txt', sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=TRUE, row.names=1)

#Take data without groups
df <- ExpressionTable[2:16]

#Generate PCA components
pca_res <- prcomp(df, scale.=TRUE)

#Plot PCA plot
plot <- autoplot(pca_res, data=ExpressionTable, frame=F, colour='Group', frame.type='norm')
plot(plot)

#Create Frame
plot2 <- autoplot(pca_res, data=ExpressionTable, frame=TRUE, colour='Group', frame.type='norm')
plot(plot2)

#Add labels
plot3 <- autoplot(pca_res, data=ExpressionTable, frame=TRUE, colour='Group', frame.type='norm', label=TRUE, label.size=3)
plot(plot3)


#Print top scaled values of pca
head(pca_res$scale)