# ----------------------------------------
# Variance EWAS script - version 1 of 2 - with no cell count adjustment
# Generation Scotland
# initially written by Tom Battram
# Edited by Sarah Watkins, July 2024
#
# Please look through the script to change:
# file names and paths
# batch variable names
# phenotype name
# ----------------------------------------

## This script uses the "jlst" R package to test whether there are differences in
## DNAm variance, at each CpG site, between cases and controls

## pkgs
library(tidyverse) # tidy code and data
## devtools::install_github("jrs95/jlst")
library(jlst) # functions for variance tests
library(usefunc) # own package of useful functions

# loading and saving files
# please use either the command args or file name options below -
# whichever you're more comfortable with

# args <- commandArgs(trailingOnly = TRUE)
# phen_file <- args[1]
# meth_file <- args[2]
# svs_file <- args[3] # make this "" if not using
# covar_file <- args[4]
# pheno_summary_table_outfile_1 <- args[5]
# pheno_summary_table_outfile_2 <- args[6]
# pheno_summary_table_outfile_3 <- args[7]
# out_file <- args[8]

# message("Arguments are: ", args)

phen_file <- "pheno_and_samplesheet-ahrr.tsv"
meth_file <- "clean-meth.RData"
svs_file <- "ad-svs.tsv"
covar_file <- "covars-cc.txt"
pheno_summary_table_outfile_1 <- "[cohort_name]_var-ewas-pheno_summary_numeric_vars_meanSD-cc.csv"
pheno_summary_table_outfile_2 <- "[cohort_name]_var-ewas-pheno_summary_numeric_vars_ttest-cc.csv"
pheno_summary_table_outfile_3 <- "[cohort_name]_var-ewas-pheno_summary_categorical_vars_counts-cc.csv"
out_file <- "[cohort_name]_varewas-res-cc.tsv"

use_covs <- ifelse(grepl("cc", out_file), TRUE, FALSE)

message("It is ", use_covs, " that covariates will be used in this analysis")

## read in data
phen_dat <- read_tsv(phen_file)
# please change column IID to your own variable name!
phen_dat$IID <- as.character(phen_dat$IID)
meth <- new_load(meth_file)
covariates <- readLines(covar_file)

# ----------------------------------------------------------------
# ewas functions
# ----------------------------------------------------------------

prep_pheno_data <- function(phen, pheno_dat, svs_file, IID, covs)
{
  message("Preparing phenotype data")
  if (svs_file == "") return(na.omit(pheno_dat[, c(IID, phen, covs)]))
  ## read in svs and cell counts
  svs <- read_tsv(svs_file)
  # please change column IID to your own variable name!
  svs$IID <- as.character(svs$IID)
  
  sv_nam <- grep("sv", colnames(svs), value = T)
  # we don't need this next line as SVs are already in the covs file
  #covs <- c(covs, sv_nam)
  # Prepare phenotype data
  temp_phen <- pheno_dat %>%
    left_join(svs, by = setNames(IID, IID)) %>%
    dplyr::select(one_of(IID), one_of(phen), one_of(covs)) %>%
    na.omit(.)
  
  message("Phenotype data prepared")
  
  return(temp_phen)
}

varewas <- function(use_covs, cpg, temp_phen, temp_meth, covs) {
  covar_dat <- as.data.frame(temp_phen[, covs])
  if(!use_covs) covar_dat <- NULL
  var_out <- tryCatch({
    vartest(y = as.numeric(temp_meth[rownames(temp_meth) == cpg, , drop = TRUE]), 
            x = as.factor(temp_phen[[phen]]), 
            covar = covar_dat, 
            covar.var = use_covs, 
            type = 2
    )
  }, error = function(e) {
    usr_m <- paste0("Error in variance test of ", cpg)
    err_msg(e, r_msg = TRUE, user_msg = usr_m, to_return = cpg)
  })
  return(var_out)
}

tidy_var_res <- function(var_res, cpg)
{
  CO <- var_res$coef
  CO <- CO[which(rownames(CO) == "x1"),]
  CO[names(CO) == "Estimate"]
  out <- data.frame(probeID = cpg, 
                    BETA = CO[names(CO) == "Estimate"], 
                    SE = CO[names(CO) == "Std. Error"], 
                    P = CO[names(CO) == "Pr(>|t|)"]
  )
  rownames(out) <- NULL
  return(out)
}

calc_increments <- function(length_task, percent_increments = 10) 
{
  ### Calculate increments to output message
  increments <- round(length_task / percent_increments)
  message_at <- round(seq(increments, length_task, length.out = percent_increments))
  names(message_at) <- seq(percent_increments, 100, percent_increments)
  return(message_at)
}

output_percent_complete <- function(n_task, increments) 
{
  ### output percentage complete
  if (n_task %in% increments) {
    percent <- increments[increments == n_task]
    message(names(percent), "% complete.")
  }
}


run_ewas <- function(phen, pheno_dat, svs_file, meth_dat, IID, out_file, covs) 
{
  # prep pheno data
  temp_phen <- prep_pheno_data(phen, pheno_dat, svs_file, IID, covs)
  # sv_nam <- grep("sv", colnames(temp_phen), value = T)
  # covs <- c(covs, sv_nam)
  covs <- covs[covs %in% colnames(temp_phen)]
  print("covs are:")
  print(covs)
  
  # Match meth to Pheno
  temp_meth <- meth_dat[, na.omit(match(temp_phen[[IID]], as.character(colnames(meth_dat))))]
  temp_phen <- temp_phen[match(as.character(colnames(temp_meth)), temp_phen[[IID]]), ]
  
  # Get cases and controls per probe
  # Please edit if your phenotype is not 0=control, 1=case
  cases <- temp_phen[temp_phen[[phen]] == 1, IID, drop = TRUE]
  print("cases:")
  print(head(cases))
  
  n_cases <- rowSums(!is.na(temp_meth[, cases]))
  n_controls <- rowSums(!is.na(temp_meth[, !colnames(temp_meth) %in% cases]))
  probe_n_cc <- as_tibble(cbind(probeID = names(n_cases), N = n_cases + n_controls, N_cases = n_cases, N_controls = n_controls))
  
  ################
  
  ## new section - 11/07/24
  ## output descriptive stats for table 1
  
  # report on overall n of cases and controls
  total_n_cases <- ncol(temp_meth[, cases])
  total_n_controls <- ncol(temp_meth[, !colnames(temp_meth) %in% cases])
  
  print(paste("n cases:",total_n_cases))
  print(paste("n controls:",total_n_controls))
  
  # report on other variables
  cases_temp_phen <- temp_phen[cases,]
  control_temp_phen <- temp_phen[!rownames(temp_phen) %in% cases,]
  
  phen_summary_numeric_mean <- aggregate(temp_phen,by=list(temp_phen$ad),FUN="mean")
  print(phen_summary_numeric_mean)
  phen_summary_numeric_mean$stat <- "mean"
  phen_summary_numeric_sd <- aggregate(temp_phen,by=list(temp_phen$ad),FUN="sd")
  print(phen_summary_numeric_sd)
  phen_summary_numeric_sd$stat <- "sd"
  phen_summary_numeric <- rbind(phen_summary_numeric_mean,phen_summary_numeric_sd)
  write.csv(phen_summary_numeric,file=pheno_summary_table_outfile_1)
  
  # please edit the next two lines with your variable names
  ttest_vars <- temp_phen[,!colnames(temp_phen)%in%c(IID,"sex")]
  phen_summary_numeric_ttest <- sapply(ttest_vars[,!colnames(ttest_vars)=="ad"],function(x) t.test(x~ad,data=ttest_vars)$p.value)
  print(phen_summary_numeric_ttest)
  write.csv(phen_summary_numeric_ttest,file=pheno_summary_table_outfile_2)
  
  categorical_vars <- colnames(temp_phen[,sapply(temp_phen,class)%in%c('character','factor')])
  categorical_vars <- categorical_vars[!categorical_vars%in%c(IID)]
  # please edit with your variable name ($ad)
  phen_summary_count <- aggregate(temp_phen[,categorical_vars],by=list(temp_phen$ad),FUN="table")
  write.csv(phen_summary_count,file=pheno_summary_table_outfile_3)
  
  #################
  
  if (!all(temp_phen[[IID]] == colnames(temp_meth))) stop("phenotype and DNAm data not matched.")
  
  print(temp_phen)
  
  increments <- calc_increments(length_task = nrow(temp_meth))
  
  # Run variance EWAS using jlst
  message("Running EWAS")
  # please check this includes all expected covariates:
  # age, sex, SVs and cell counts
  print("about to run var EWAS. Covariates are:")
  print(covs)
  
  list_res <- lapply(1:nrow(temp_meth), function(x) {
    output_percent_complete(x, increments)
    cpg <- rownames(temp_meth)[x]
    
    res <- varewas(use_covs, cpg, temp_phen, temp_meth, covs)
    if (class(res) == "character") return(res)
    out <- tidy_var_res(res, cpg)
    return(out)
  })
  
  # free up some space
  rm(temp_meth)
  
  succ_res <- list_res[sapply(list_res, class) != "character"]
  
  res <- bind_rows(succ_res) %>%
    left_join(probe_n_cc, by = "probeID")
  
  write.table(res, file = out_file, sep = "\t", col.names = T, row.names = F, quote = F)
  res <- as.data.frame(res)
  res <- res[order(res$P),]
  print("Results!:")
  print(head(res))
}


# ----------------------------------------
# Run the EWAS
# ----------------------------------------
phen <- "ad"

run_ewas(phen = phen, 
         pheno_dat = phen_dat, 
         svs_file = svs_file,
         meth_dat = meth,
         IID = "IID", 
         out_file = out_file, 
         covs = covariates)

print("FIN")