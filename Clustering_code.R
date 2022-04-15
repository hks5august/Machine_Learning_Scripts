#Load libraries
library(stats)
library(tools)

#Load data
data <- read.table ("https://raw.githubusercontent.com/pine-bio-support/Unsupervised-codeomics/main/CellLines_7000Genes_marked_UI.txt", header = TRUE, row.names = 1)

#Transpose data
data1 <- t(data)

#cutree to get the required number of clusters and output cluster IDs for all the samples
clust_ids <- cutree(cl, k=5)
head(clust_ids)

#Output table with cluster IDs 
txt_file_c <- paste("file_name", "_", "cRes_hclust.ClustersOnly", ".txt", sep="")
cat ("\"Object\"\t\"ClusterID\"\n", file=txt_file_c, append=FALSE)

#Write into a file
write.table (clust_ids, file=txt_file_c, append=TRUE, sep="\t", col.names=FALSE, row.names=TRUE)

#Output full table with cluster IDs 
txt_file <- paste("file_name", "_", "cRes_hclust.FullTable", ".txt", sep="")	
out_table <- cbind (clust_ids, data)
colnames (out_table) [1] <- "ClusterID"
cat ("\"Object\"\t", file="Full_table.txt", append=FALSE)

#Write into a file
write.table (out_table, file="Full_table.txt", append=TRUE, sep="\t", col.names=TRUE, row.names=TRUE)

dim(out_table)





#kmean clustering
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

#Load data
data <- read.table ("https://raw.githubusercontent.com/pine-bio-support/Unsupervised-codeomics/main/CellLines_7000Genes_marked_UI.txt", header = TRUE, row.names = 1)

#Transpose data
data1 <- t(data)


#Perform scaling of data required for k-mean clustering
df <- scale(data1)

#Perform k-mean clustering and define value of centers or the clusters that are required
kmean_res <- kmeans(df, centers = 5, nstart = 25)

#visualize clustering
fviz_cluster(kmean_res, data = df, ellipse = TRUE, ellipse.type = "convex", pointsize = 1.5, labelsize = 7)


#Extract cluster ID with sample
clusters <- kmean_res$cluster


#Write into a file

cat ("\"Object\"\t Cluster_ID\n", file="kmean_results.txt", append=TRUE)

#Write into a file
write.table (clusters, file="kmean_results.txt", append=TRUE, sep="\t", col.names=F, row.names=T)



