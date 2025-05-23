library(parallel)
library(MASS)
library(tidyr)
library(dplyr)
suppressMessages(library(meffil))
suppressMessages(library(matrixStats))

main <- function()
{
  arguments <- commandArgs(T)
  
  filepath <- arguments[1]
  methylationfile <- arguments[2]
  grmfile <- arguments[3]
  cov_file <- arguments[4]
  out_file <- arguments[5]
  # removed transform option as we're only running untransformed
  #transform <- arguments[6]
  #print(transform)
  # removed nthreads option because it was causing issues
  #nthreads <- as.numeric(arguments[7])
  chunks <- as.numeric(arguments[8])
  jid <- as.numeric(arguments[9])
  meth_array <- arguments[10]
  
  source(paste0(filepath, "/resources/methylation/polygenic_genable.R"))
  source(paste0(filepath, "/resources/methylation/polylik_genable.R"))
  
  message("Reading methylation data...")
  load(methylationfile)
  
  if(!is.na(jid))
  {
    chunksize <- ceiling(nrow(norm.beta) / chunks)
    i1 <- chunksize * (jid-1) + 1
    i2 <- min(nrow(norm.beta), chunksize * jid)
    norm.beta <- norm.beta[i1:i2,]
  } 
  
  # Remove all IDs that have any NAs in the covariate file
  covs <- read.table(cov_file, he=T,stringsAsFactors=F,colClasse=c("Sex_factor"="character"))
  covs <- covs[,!colnames(covs)=="Treg"]
  
  rownames(covs) <- covs$IID
  covs <- subset(covs, IID %in% colnames(norm.beta), select=-c(IID))
  norm.beta <- norm.beta[, colnames(norm.beta) %in% rownames(covs)]
  covs <- covs[match(colnames(norm.beta), rownames(covs)), ]
  stopifnot(all(rownames(covs) == colnames(norm.beta)))
  
  ## process CpGs on sex chromosomes
  annots <- meffil.get.features(meth_array)
  annots <- annots[!is.na(annots$chromosome),]
  x_probes <- annots$name[annots$chromosome == "chrX"]
  y_probes <- annots$name[annots$chromosome == "chrY"]
  
  beta.x.female <- norm.beta[rownames(norm.beta) %in% x_probes, colnames(norm.beta) %in% rownames(covs)[covs$Sex_factor=="F"]]
  beta.x.male <- norm.beta[rownames(norm.beta) %in% x_probes, colnames(norm.beta) %in% rownames(covs)[covs$Sex_factor=="M"]]
  beta.y.male <- norm.beta[rownames(norm.beta) %in% y_probes, colnames(norm.beta) %in% rownames(covs)[covs$Sex_factor=="M"]]
  
  if(is.matrix(beta.x.female) == F){
    beta.x.female <- matrix(beta.x.female, nrow=1)
  }
  rownames(beta.x.female) <- rownames(norm.beta)[rownames(norm.beta) %in% x_probes]
  colnames(beta.x.female) <- colnames(norm.beta)[colnames(norm.beta) %in% rownames(covs)[covs$Sex_factor=="F"]]
  
  
  if(is.matrix(beta.x.male) == F){
    beta.x.male <- matrix(beta.x.male, nrow=1)
  }
  rownames(beta.x.male) <- rownames(norm.beta)[rownames(norm.beta) %in% x_probes]
  colnames(beta.x.male) <- colnames(norm.beta)[colnames(norm.beta) %in% rownames(covs)[covs$Sex_factor=="M"]]
  
  
  if(is.matrix(beta.y.male) == F){
    beta.y.male <- matrix(beta.y.male, nrow=1)
  }
  rownames(beta.y.male) <- rownames(norm.beta)[rownames(norm.beta) %in% y_probes]
  colnames(beta.y.male) <- colnames(norm.beta)[colnames(norm.beta) %in% rownames(covs)[covs$Sex_factor=="M"]]
  
  ## process CpGs on autosomes
  autosomal_probes <- annots$name[!annots$chromosome %in% c("chrX","chrY")]
  norm.beta <- norm.beta[rownames(norm.beta) %in% autosomal_probes,]
  
  
  g <- grep("_factor",names(covs))
  if(length(g) > 0)
  {
    for (i in 1:length(g))
    {
      covs[,g[i]]<-as.factor(covs[,g[i]])
    }
  }
  
  covs.female <- covs[covs$Sex_factor=="F",]
  covs.male <- covs[covs$Sex_factor=="M",]
  
  grm <- readGRM(grmfile)
  kin <- makeGRMmatrix(grm)
  kin <- kin[rownames(kin) %in% colnames(norm.beta), colnames(kin) %in% colnames(norm.beta)]
  
  message("Kinship matrix is ", nrow(kin), " by ", nrow(kin))
  
  index <- match(rownames(kin), colnames(norm.beta))
  norm.beta <- norm.beta[,index]
  
  index_female <- match(rownames(kin), colnames(beta.x.female))
  beta.x.female <- beta.x.female[,index_female]
  
  index_male_x <- match(rownames(kin), colnames(beta.x.male))
  beta.x.male <- beta.x.male[,index_male_x]
  
  index_male_y <- match(rownames(kin), colnames(beta.y.male))
  beta.y.male <- beta.x.male[,index_male_y]
  
  covs <- covs[match(colnames(norm.beta), rownames(covs)), ]
  stopifnot(all(rownames(kin) == colnames(norm.beta)))
  stopifnot(all(rownames(covs) == colnames(norm.beta)))
  
  message("Calculating eigenvectors")
  
  relmat <- kin * 2
  tmp <- t(relmat)
  relmat[upper.tri(relmat)] <- tmp[upper.tri(tmp)]
  eig <- eigen(relmat, symmetric=TRUE)
  class(eig) <- "list"
  message(class(eig))
  print(str(eig))
  
  rm(tmp, relmat)
  
  if(nrow(norm.beta)>0 & ncol(norm.beta)>0){
    index<-sapply(covs,function(.col){all(is.na(.col) | .col[1L] == .col)})
    index[is.na(index)] <- FALSE
    covs <- covs[match(colnames(norm.beta), rownames(covs)),!index]
    run.adjust.cov(norm.beta, covs, kin, eig, paste0(out_file,".",jid,".rds"))
  }
  
  if(nrow(beta.x.female)>0 & ncol(beta.x.female)>0){
    index_female <-sapply(covs.female,function(.col){all(is.na(.col) | .col[1L] == .col)})
    index_female[is.na(index_female)] <- FALSE
    covs.female <- covs.female[match(colnames(beta.x.female), rownames(covs.female)),!index_female]
    run.adjust.cov(beta.x.female, covs.female, kin, eig, paste0(out_file, ".Female.chrX.", jid, ".RData"))
  }
  
  if(nrow(beta.x.male)>0 & ncol(beta.x.male)>0){
    index_male <-sapply(covs.male,function(.col){all(is.na(.col) | .col[1L] == .col)})
    index_male[is.na(index_male)] <- FALSE
    covs.male <- covs.male[match(colnames(beta.x.male), rownames(covs.male)),!index_male]
    run.adjust.cov(beta.x.male, covs.male, kin, eig, paste0(out_file,".Male.chrX.", jid, ".RData"))
  }
  
  if(nrow(beta.y.male)>0 & ncol(beta.y.male)>0){
    index_male <-sapply(covs.male,function(.col){all(is.na(.col) | .col[1L] == .col)})
    index_male[is.na(index_male)] <- FALSE
    covs.male <- covs.male[match(colnames(beta.x.male), rownames(covs.male)),!index_male]
    run.adjust.cov(beta.y.male, covs.male, kin, eig, paste0(out_file,".Male.chrY.", jid, ".RData"))
  }
}

run.adjust.cov <- function(betas, covs, kin, eig, out_file)
{
  #print(transform)
  if(nrow(betas) > 0){
    betas.copy <- is.na(betas)
    
    # we don't want to run in serial
    #if(is.na(nthreads) | nthreads == 1){
    #  out <- adjust.relatedness.serial(betas, covs, kin, eig)
    #} else {
      #message("Running with ", nthreads, " threads")
      out <- adjust.relatedness(betas, covs, kin, eig)
    #}
    
    betas <- out$x
    classes <- data.frame(cpg=rownames(betas), cl=out$cl)
    betas[betas.copy] <- NA
    
    index <- which(is.na(betas), arr.ind = TRUE) 
    if (length(index)>0){
      message("Replace ",length(index)," missing values with rowmeans")
      betas[index] <- rowMeans(betas, na.rm = TRUE)[index[, "row"]] 
    }
    saveRDS(betas, file=out_file)
    save(classes, file=paste0(out_file, "_classes"))
  }
}

readGRM <- function(rootname)
{
  bin.file.name <- paste(rootname, ".grm.bin", sep="")
  n.file.name <- paste(rootname, ".grm.N.bin", sep="")
  id.file.name <- paste(rootname, ".grm.id", sep="")
  
  cat("Reading IDs\n")
  id <- read.table(id.file.name)
  n <- dim(id)[1]
  cat("Reading GRM\n")
  bin.file <- file(bin.file.name, "rb")
  grm <- readBin(bin.file, n=n*(n+1)/2, what=numeric(0), size=4)
  close(bin.file)
  cat("Reading N\n")
  n.file <- file(n.file.name, "rb")
  N <- readBin(n.file, n=n*(n+1)/2, what=numeric(0), size=4)
  close(n.file)
  
  cat("Creating data frame\n")
  l <- list()
  for(i in 1:n)
  {
    l[[i]] <- 1:i
  }
  col1 <- rep(1:n, 1:n)
  col2 <- unlist(l)
  grm <- data.frame(id1=col1, id2=col2, N=N, grm=grm)	
  
  ret <- list()
  ret$grm <- grm
  ret$id <- id
  return(ret)
}

makeGRMmatrix <- function(grm)
{
  mat <- diag(nrow(grm$id))
  mat[upper.tri(mat, diag=TRUE)] <- grm$grm$grm
  mat <- t(mat)
  nsnpvec <- subset(grm$grm, id1 != id2)$N
  mat[upper.tri(mat, diag=FALSE)] <- nsnpvec
  rownames(mat) <- grm$id$V2
  colnames(mat) <- grm$id$V2
  return(mat)
}

adjust.relatedness.fast.1 <- function(x, covs, kin, eig, quiet=TRUE)
{
  x[!is.finite(x)] <- mean(x, na.rm=T)
  print("x size:")
  print(dim(x))
  print(length(x))

  d <- data.frame(X=x, covs)
  print(paste("names (d) = =",names(d)))
  
  # we aren't adjusting for covariates at this stage because we will do this in the EWAS
  #form <- as.formula(paste0("X ~ ", paste(names(d)[-1], collapse=" + ")))
  #d$X <- residuals(lm(form, d, na.action = na.omit))
  
  p_out <- try(polygenic(X, data=d, kinship.matrix=kin, eigenOfRel=eig, quiet=TRUE))
  
  iter <- 1
  while(class(p_out) == "try-error" & iter < 5)
  {
    message("trying again...")
    iter <- iter + 1
    p_out <- try(polygenic(X, data=d, kinship.matrix=kin, eigenOfRel=eig, quiet=TRUE))
  }
  if(class(p_out) == "try-error")
  {
    message("giving up, just using fixed effects model")
    #if(transform == "Yes"){
    #  a <- as.numeric(rntransform(d$X))
    #}
    #else{
      a <- d$X
    #}
  } else {
    #if(transform == "Yes"){
     # a <- as.numeric(rntransform(p_out$grresidualY))
    #}
    #else{
      a <- as.numeric(p_out$grresidualY)
    #}
  }
  return(list(x=a, cl=class(p_out)))
} 

adjust.relatedness.serial <- function(B, covs, kin, eig)
{
  cl <- array(0, nrow(B))
  #print(transform)
  for(i in 1:nrow(B))
  {
    message("Probe ",i, " of ", nrow(B))
    out <- adjust.relatedness.fast.1(B[i,], covs, kin, eig)
    B[i, ] <- out$x
    cl[i] <- out$cl
  }
  return(list(x=B, cl=cl))
}

get.index.list <- function(n, mc.cores)
{
  #print(mc.cores)
  B <	
    mc.cores <- ifelse(mc.cores < 1, 1, mc.cores)
  div <- floor(n / mc.cores)
  rem <- n %% mc.cores
  l1 <- lapply(1:div, function(x) (x-1) * mc.cores + 1:mc.cores)
  if(rem != 0) l1[[div+1]] <- l1[[div]][mc.cores] + 1:rem
  return(l1)
}

rntransform <- function(x)
{
  out <- rank(x) - 0.5
  out[is.na(x)] <- NA
  mP <- 0.5/max(out, na.rm = T)
  out <- out/(max(out, na.rm = T) + 0.5)
  out <- scale(qnorm(out))
  out
}

adjust.relatedness <- function(B, covs, kin, mc.cores=mc.cores)
{
  print(mc.cores)
  l1 <- get.index.list(nrow(B), mc.cores)
  #print(transform)
  
  l <- lapply(l1, function(ii)
  {
    res <- mclapply(ii, function(i)
    {
      message("Probe ", i, " of ", nrow(B))
      out <- adjust.relatedness.fast.1(B[i,], covs, kin, eig)
    }, mc.cores=mc.cores, mc.preschedule=FALSE)
    a <- do.call(rbind, lapply(res, function(x) x$x))
    b <- sapply(res, function(x) x$cl)
    return(list(x=a, cl=b))
  })
  x <- do.call(rbind, lapply(l, function(x) x$x))
  cl <- unlist(lapply(l, function(x) x$cl))
  rownames(x) <- rownames(B)
  colnames(x) <- colnames(B)
  return(list(x=x, cl=cl))
}

main()
