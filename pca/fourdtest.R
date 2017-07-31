start <- Sys.time()

#factoextra may become more useful later
# install.packages("devtools")
#devtools::install_github("kassambara/factoextra")
library("ggfortify")
library("factoextra")
require("factoextra")
require("ggfortify")


#CHANGE DIRECTORY
fourdset = read.csv("~/_ISB/_workspace/pca/fourdset.csv", header=TRUE)

#not used right now
sdx = sd(fourdset[,1])
sdy = sd(fourdset[,2])
sdz = sd(fourdset[,3])
sdt = sd(fourdset[,4])

fourdset.pca <- prcomp(fourdset,scale.=TRUE)
print(summary(fourdset.pca))

print(fourdset.pca)

#IMPORTANT PLOT OF POINTS WITH EIGENS
#biplot(fourdset.pca,scale=0)

autoplot(fourdset.pca)

#variability confirmation tests

sdev <- fourdset.pca$sdev

fourdset.var <- sdev^2

fourdset.varexp <- fourdset.var/sum(fourdset.var)

#print(fourdset.varexp)
#IMPORTANT PLOT OF VARIABILITIES' EXPLANATIONS
#plot(fourdset.varexp,xlab="pc", ylab="proportion of var explained", type = "b")

#everything below is for testing or not used yet

x <- data.frame(fourdset[,1])
y <- data.frame(fourdset[,2])
z <- data.frame(fourdset[,3])

#print(x)
#print(y)

fourdtable <- table(x$y)

#barplot(fourdtable, main="xd", xlim = c(0,100) , ylim = c(0,100), xlab="x", ylab="y", col="steelblue")

end <- Sys.time()
print(end - start)