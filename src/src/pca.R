start <- Sys.time()

file.path <- "~\\projects\\transcriptomics\\data\\"
file.name <- "expressionMatrix_FPKM.csv"
ribo.name <- "riboGenes.csv"

file <- paste(file.path, file.name, sep = '')
pre.data <- read.csv(file = file, row.names = 1)

post.data <- pre.data
for (i in 1:24) {
  post.data <- post.data[post.data[,i] != 0,]
}
post.data <- log(post.data)

par(mfrow = c(1,2))

pca <- prcomp(post.data)
loadings <- pca$rotation[,1:2]
PC1 <- loadings[,1]
PC2 <- loadings[,2]
plot(PC1, PC2)

ribo.file <- paste(file.path, ribo.name, sep = '')
ribo.genes <- unname(unlist(read.csv(file = ribo.file, header = FALSE)[1,]))

ribo.data <- post.data[rownames(post.data) %in% ribo.genes,]
ribo.pca <- prcomp(ribo.data)
ribo.loadings <- ribo.pca$rotation[,1:2]
ribo.PC1 <- ribo.loadings[,1]
ribo.PC2 <- ribo.loadings[,2]
plot(ribo.PC1, ribo.PC2)

end <- Sys.time()
print(end - start)
