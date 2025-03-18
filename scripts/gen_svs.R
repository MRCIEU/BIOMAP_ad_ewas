# -------------------------------------------------------
# Script to generate surrogate variables 
#
# Please look through the script to change:
# file names and paths
# batch variable names
# phenotype name
# -------------------------------------------------------

## pkgs
library(tidyverse) # tidy data and code
library(sva) # calculating surrogate variables
library(SmartSVA) # calculating SVs
library(matrixStats) # imputing DNAm data
# devtools::install_github("thomasbattram/usefunc")
library(usefunc) # personal package of useful functions
library(data.table)

# loading and saving files
# please use either the command args or file name options below -
# whichever you're more comfortable with

# Command args loading/saving option
# args <- commandArgs(trailingOnly = TRUE)
# phen_file <- args[1]
# meth_file <- args[2]
# aries_dir <- args[3]
# out_file <- args[4]
# removed_out_file <- args[5]
# covars_outfile <- args[6]
# cc_covars_outfile <- args[7]
# heatmap_outfile <- args[8]

# File names loading/saving option
# Please change file names and add paths as appropriate!
phen_file <- "ad-data-cleaned.tsv"
phen_outfile <- "pheno_and_samplesheet-ahrr.tsv"
meth_file <- "clean-meth.RData" 
out_file <- "ad-svs.tsv"
removed_out_file <- "ad-removed-svs.RData"
covars_outfile <- "covars-no-cc-svs.txt"
cc_covars_outfile <- "covars-cc-svs.txt"
heatmap_outfile <- "svs-heatmap.png"

## read in data
phen_dat <- read_tsv(phen_file)
phen_dat <- as.data.frame(phen_dat)
phen_dat$IID <- as.character(phen_dat$IID)

load(meth_file)

# names of batch variables - please change!
batch_vars <- c("set", "batch")

# -------------------------------------------------------
# functions
# -------------------------------------------------------

#' Impute missing values in DNAm matrix
#' 
#' @param x DNAm matrix
#' @param FUN function to apply to "x" to get values to impute missing values
#' 
#' @return imputed DNAm matrix
impute_matrix <- function(x, FUN = function(x) rowMedians(x, na.rm = T)) {
  idx <- which(is.na(x), arr.ind = T)
  if (length(idx) > 0) {
    v <- FUN(x)
    v[which(is.na(v))] <- FUN(matrix(v, nrow = 1))
    x[idx] <- v[idx[, "row"]]
  }
  return(x)
}

# function to add quotes for weird trait names
addq <- function(x) paste0("`", x, "`")

#' Generate surrogate variables
#' 
#' @param trait trait of interest
#' @param phen_data phenotype data (data.frame or tibble) containing trait of interest and covariates
#' @param meth_data DNAm matrix
#' @param covariates character vector of covariates
#' @param nsv number of surrogate variables to use
#' @param IID name of the identifier used for DNAm samples. Default = "Sample_Name"
#' 
#' @return data.frame of surrogate variable values for each individual
generate_svs <- function(trait, phen_data, mdat, covariates = "", nsv, 
                         IID = "IID") {
  print("Starting SV generation")
  if (any(grepl("^sv", paste(covariates, collapse = "|")))) {
    covs <- covariates[-grep("^sv", covariates)]
  } else {
    covs <- covariates
  }
  print("SV covariates=")
  print(covs)
  
  print("original dim of pheno data is:")
  print(dim(phen_data))
  print(head(phen_data))
  phen <- phen_data %>%
    dplyr::select(one_of(IID), one_of(trait, covs)) %>%
    .[complete.cases(.), ]
  print("complete case dim of pheno data is:")
  print(dim(phen))
  print(head(phen))
  joint_samples <- intersect(as.character(colnames(mdat)),as.character(phen[[IID]]))
  print("length joint samples:")
  print(length(joint_samples))
  
  mdat <- mdat[, joint_samples]
  print("mdat dimensions after IID subset:")
  print(dim(mdat))
  phen <- phen %>%
    dplyr::filter(!!as.symbol(IID) %in% colnames(mdat))
  print("phen dimensions after mdat colnames subset:")
  print(dim(phen))
  print("all equal?")
  print(all.equal(rownames(phen),colnames(mdat)))
  
  participants.temp <- as.character(colnames(mdat))
  
  # models 
  trait_mod <- paste0("~ ", addq(trait))
  cov_mod <- paste(covs, collapse = " + ")
  if (length(covs)>0) {
    full_mod <- paste(trait_mod, cov_mod, sep = " + ")
    fom <- as.formula(full_mod)
    # null model
    fom0 <- as.formula(paste0("~ ", cov_mod))
    mod0 <- model.matrix(fom0, data = phen)
  } else {
    fom <- as.formula(trait_mod)
    mod0 <- NULL
  }
  
  print(paste("Svs cov mod =",fom))
  
  # full model - with variables of interest 
  mod <- model.matrix(fom, data = phen)
  
  message("mdat dimensions: ", paste0(nrow(mdat), ", ", ncol(mdat)))
  message("pheno dimensions: ", paste0(nrow(phen), ", ", ncol(phen)))
  
  # Estimate the surrogate variables
  tryCatch({
    svobj <- smartsva.cpp(mdat, mod, mod0, n.sv = nsv, VERBOSE = T)
    svs <- as.data.frame(svobj$sv, stringsAsFactors = F)
    svs[[IID]] <- phen[[IID]]
    # head(svs)
    colnames(svs)[1:nsv] <- paste0("sv", 1:nsv)
    rownames(svs) <- participants.temp
    
    return(svs)
    print(paste("N of SVs generated is",ncol(svs)))
  }, error = function(e) {err_msg(e, r_msg = TRUE, user_msg = paste("SV fail"))})
}

# -------------------------------------------------------
# sort data for generating SVs
# -------------------------------------------------------

## methylation data
print("methylation data:")
print(dim(meth))
mdata <- impute_matrix(meth)
rm(meth)
print("dim mdata:")
print(dim(mdata))

# 26/6/24 - added in as we want to adjust for smoking
meth_ahrr <- as.data.frame(t(mdata))
# phenotype data
# 26/6/24 moved down so that we can add in ahrr as a covariate from the imputed DNAm matrix
pheno <- phen_dat
pheno <- as.data.frame(pheno)
print("mean age =")
print(mean(pheno$age,na.rm = T))
rm(phen_dat)
# please change "IID" to your identifier name!
rownames(pheno) <- pheno$IID

pheno$ahrr <- meth_ahrr$cg05575921[match(rownames(pheno), rownames(meth_ahrr))]
rm(meth_ahrr)
print("pheno after ahrr:")
print(dim(pheno))
print("this is the pheno dataframe: are all the variables you want here? (check for sentrix ID??)")
write.table(pheno, file = phen_outfile, quote = F, col.names = T, row.names = F, sep = "\t")

# we want covs to be age, sex, ahrr (smoking), and cell counts
# please change to your own variable names!
# IID = sample study ID, Sample_Sentrix_ID = sample array ID, ad = phenotype
covs <- colnames(pheno)[!colnames(pheno) %in% c("IID", "Sample_Sentrix_ID", "ad", batch_vars)]

print("covs =")
print(covs)

# please change to your own variable name!
phen <- "ad"

## adjust methylation data for cell counts and take residuals

# make sure to edit the below line if any other covariates are added in your study
# and please change to your own variable names!
celltypes <- covs[!covs %in% c("age", "sex", "ahrr")]
form <- paste(paste0("pheno$", celltypes), collapse = " + ")

ori_time <- proc.time()
class(mdata)

participants <- colnames(mdata)
meth_resid <- lapply(1:nrow(mdata), function(x) {
  dnam <- as.numeric(mdata[x,, drop=T])
  resid(lm(as.formula(paste0("dnam ~ ", form))))
})
names(meth_resid) <- rownames(mdata)
mdata <- do.call(rbind, meth_resid)
colnames(mdata) <- participants

message("mdata dimensions: ", paste0(nrow(mdata), ", ", ncol(mdata)))
time_taken <- proc.time() - ori_time
time_taken 
rm(meth_resid)

# -------------------------------------------------------
# Generate SVs
# -------------------------------------------------------

# First we estimate the number of SVs using random matrix theory
# option 1 for this will work for most datasets and should be the default option
# if this throws an error relating to maximum vector size being exceeded,
# use option 2.

# option 1
n.sv <- EstDimRMT(mdata, FALSE)$dim + 1
print(paste("n.sv = ",n.sv))

# option 2
# n.sv.1 <- EstDimRMT(mdata[,1:2000], FALSE)$dim + 1
# print(paste("n.sv.1 = ",n.sv.1))
# n.sv.2 <- EstDimRMT(mdata[,2001:4000], FALSE)$dim + 1
# print(paste("n.sv.2 = ",n.sv.2))
# n.sv.3 <- EstDimRMT(mdata[,4001:ncol(mdata)], FALSE)$dim + 1
# print(paste("n.sv.3 = ",n.sv.3))
# sv_combo <- c(n.sv.1,n.sv.2,n.sv.3)
# n.sv <- mean(sv_combo)
# print("n.sv =")
# print(n.sv)

# 19/06/24 - SW - changed nsv so we dont squish all the variance into 10 SVs
svs <- generate_svs(trait = phen, 
                    phen_data = pheno, 
                    mdat = mdata, 
                    covariates = covs[!covs %in% celltypes], 
                    #nsv = 10, 
                    nsv=n.sv,
                    IID = "IID")

# svs <- read_tsv(out_file)

n_of_svs <- paste("n of SVs generated for analysis was:",ncol(svs))
print(paste("n of SVs generated for analysis was:",ncol(svs)))
print(paste("n of participants SVs generated for analysis was:",nrow(svs)))

# -------------------------------------------------------
# Check association between SVs and phenotype of interest
# -------------------------------------------------------

sv_nam <- grep("sv", colnames(svs), value=T)
print(length(sv_nam))
svs <- svs[,1:10]
svs$IID <- rownames(svs)
sv_nam <- grep("sv", colnames(svs), value=T)
print(head(svs))
print(length(sv_nam))

#' Assess association with SVs
#' 
#' @param svs table of SV values
#' @param pheno phenotype data
#' @param trait trait of interest (cov or phenotype)
regress_sv <- function(svs, pheno, trait)
{
  sv_check_dat <- pheno %>%
    left_join(svs) %>%
    dplyr::select(one_of(trait, sv_nam)) %>%
    na.omit
  
  sv_assoc <- map_dfr(sv_nam, function(x) {
    print(x)
    form <- as.formula(paste(x, trait, sep = " ~ "))
    fit <- lm(form, data = sv_check_dat)
    out_nums <- summary(fit)$coef[2, ]
    out_r2 <- summary(fit)$adj.r.squared
    out <- as_tibble(t(as.matrix(out_nums))) %>%
      mutate(sv = x, r2 = out_r2) %>%
      dplyr::select(sv, beta = Estimate, SE = `Std. Error`, t_val = `t value`, P = `Pr(>|t|)`, r2)
    return(out)
  })
  
  ## remove associations at P<0.01 (changing from P<0.05 as doing 10 tests here...)
  sv_to_rm <- sv_assoc %>% 
    dplyr::filter(P < 0.01) %>%
    pull(sv)
  
  out <- list(assoc = sv_assoc, to_rm = sv_to_rm)
  return(out)
}

# change "IID" and "Sample_Sentrix_ID" to your own variable names!
sv_vars <- colnames(pheno)[!colnames(pheno) %in% c("IID", "Sample_Sentrix_ID")]
print("testing sv 1-10 association with:")
print(sv_vars)
sv_assoc_list <- lapply(sv_vars, regress_sv, svs = svs, pheno = pheno)
names(sv_assoc_list) <- sv_vars

sv_assoc <- bind_rows(map(sv_assoc_list, "assoc"), .id = "variable")

## reorder the variables for the plot
sv_assoc$sv <- factor(sv_assoc$sv, levels = sv_nam)
var_levels <- c(sv_vars[!sv_vars %in% celltypes], celltypes)
sv_assoc$variable <- factor(sv_assoc$variable, levels = var_levels)

## Table
# SV | trait | beta

heatmap_r2 <- ggplot(sv_assoc, aes(sv, variable, fill = r2)) +
  geom_tile(color = "white") + 
  geom_text(aes(label = round(r2, 2))) + 
  scale_fill_gradient2(low = "white", high = "red", space = "Lab", name="r2")

ggsave(heatmap_outfile, plot = heatmap_r2)


## REMOVE SVS ASSOC WITH TRAIT HERE
# change "ad" to your own variable name!
sv_to_rm <- sv_assoc %>%
  dplyr::filter(variable == "ad") %>%
  dplyr::filter(P < 0.01) %>%
  pull(sv)

svs_out <- svs %>%
  dplyr::select(-one_of(sv_to_rm))

# -------------------------------------------------------
# Write out results and info about removed SVs
# -------------------------------------------------------

write.table(svs_out, file = out_file,
            sep = "\t", quote = F, col.names = T, row.names = F)

sv_rem_info <- lapply(sv_to_rm, function(x) {svs[[x]]})
names(sv_rem_info) <- sv_to_rm
# change "ad" to your own variable name!
sv_rem_info$sv_assoc <- sv_assoc %>%
  dplyr::filter(variable == "ad")

save(sv_rem_info, file = removed_out_file)

## Add SVs to covariate file
sv_covs <- grep("sv", colnames(svs_out), value=T)
message("Appending the following variables to the covariates file: ", paste(sv_covs, collapse = ", "))
covars <- readLines(covars_outfile)
covars <- unique(c(covars, sv_covs))
covars <- c(covars, "ahrr")
write(covars, file = covars_outfile)
covars_cc <- readLines(cc_covars_outfile)
covars_cc <- unique(c(covars_cc, sv_covs))
covars_cc <- c(covars_cc, "ahrr")
write(covars_cc, file = cc_covars_outfile)
