start <- Sys.time()

file.path <- "~\\projects\\transcriptomics\\data\\"
file.name <- "expressionMatrix_FPKM.csv"
ribo.name <- "riboGenes.csv"

q = 0.970 # 95% confidence

Col <- function(x) {
  apply(cbind(pre.data[,x * 3 - 2], pre.data[,x * 3 - 1], pre.data[,x * 3]),
    1, Outlier)
}

# Dixon's Q Test
Outlier <- function(x) {
  # x <- sort(x)
  # q.max <- (x[3] - x[2]) / (x[3] - x[1])
  # q.min <- (x[2] - x[1]) / (x[3] - x[1])
  # if (!is.na(q.max) && q.max > q.min && q.max > q) {
  #   x[3] <- NA
  #   return(mean(x, na.rm = TRUE))
  # }
  # if (!is.na(q.min) && q.min > q) {
  #   x[1] <- NA
  #   return(mean(x, na.rm = TRUE))
  # }
  # return(mean(x))
}

Graph <- function(data, text) {

  fit <- lm(log10(RPF) ~ log10(mRNA), data = data)

  outlier.test <- outlierTest(fit)
  outliers <- names(outlier.test[[1]])
  print(outlier.test)

  data$Color <- 'black'
  data$Color[rownames(data) %in% outliers] <- 'red'

  xlim <- c(1e0, 1e6)
  ylim <- xlim
  plot(data$mRNA, data$RPF, log = 'xy', xlab = '', ylab = '', col = data$Color,
  xlim = xlim, ylim = ylim, main = paste(text, 'Time Point', i))

  title(xlab = 'mRNA (FPKM)', ylab = 'RPF (FPKM)')
  abline(fit, col = 'blue')

  r2 <- round(summary(fit)$r.squared, 2)
  text(2e1, 5e5, paste('r2 =', r2))

  m <- round(coef(fit)[2], 3)
  b <- round(coef(fit)[1], 3)
  text(2e1, 1.5e5, paste('y = ', m, 'x + ', b, sep = ''))

  f <- summary(fit)$fstatistic
  p <- signif(pf(f[1], f[2], f[3], lower.tail = FALSE), 3)
  text(2e1, 5e4, paste('p =', p))
}

file <- paste(file.path, file.name, sep = '')
pre.data <- read.csv(file = file, row.names = 1)

post.data <- pre.data
for (i in 1:24) {
  post.data <- post.data[post.data[,i] != 0,]
}

par(mfrow = c(1,2))

pca <- prcomp(t(post.data))
loadings <- pca$rotation[,1:2]
PC1 <- loadings[,1]
PC2 <- loadings[,2]
plot(PC1, PC2)

ribo.file <- paste(file.path, ribo.name, sep = '')
ribo.genes <- unname(unlist(read.csv(file = ribo.file, header = FALSE)[1,]))

ribo.data <- post.data[rownames(post.data) %in% ribo.genes,]
ribo.pca <- prcomp(t(ribo.data))
ribo.loadings <- ribo.pca$rotation[,1:2]
ribo.PC1 <- ribo.loadings[,1]
ribo.PC2 <- ribo.loadings[,2]
plot(ribo.PC1, ribo.PC2)

# total.data <- list()
# ribo.data <- list()
# for (i in 1:4) {
#   total.data[[i]] <- data.frame(Col(i), Col(i + 4))
#   rownames(total.data[[i]]) <- rownames(pre.data)
#   colnames(total.data[[i]]) <- c('mRNA', 'RPF')
#   total.data[[i]] <- total.data[[i]][total.data[[i]][,1] != 0,]
#   total.data[[i]] <- total.data[[i]][total.data[[i]][,2] != 0,]
#   ribo.data[[i]] <- total.data[[i]][rownames(total.data[[i]]) %in% ribo.genes,]
# }
#
# par(mfrow = c(2, 4))
#
# for (i in 1:4) {
#   Graph(total.data[[i]], 'Total')
# }
#
# for (i in 1:4) {
#  Graph(ribo.data[[i]], 'Ribosomal')
# }

end <- Sys.time()
print(end - start)
