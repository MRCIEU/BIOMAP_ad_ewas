#Function to remove outliers using the Tukey IQR*3 method

IQR.removal <- function(meth) {
  rowIQR <- rowIQRs(meth, na.rm = T)
  row2575 <- rowQuantiles(meth, probs = c(0.25, 0.75), na.rm = T)
  maskL <- meth < row2575[, 1] - 3 * rowIQR
  maskU <- meth > row2575[, 2] + 3 * rowIQR
  meth[maskL] <- NA
  meth[maskU] <- NA
  return(meth)
}

#load data
load("norm.beta.Robj")
dim(norm.beta)

#remove outliers
log.iqr <- data.frame(cpgs = row.names(norm.beta), NAs.before.IQR3 = rowSums(is.na(norm.beta)))
norm.beta <- IQR.removal(norm.beta)
log.iqr$NAs.after.IQR3 <- rowSums(is.na(norm.beta))
table(log.iqr$NAs.after.IQR3)