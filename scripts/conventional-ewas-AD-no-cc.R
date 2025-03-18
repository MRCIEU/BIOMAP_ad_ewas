# ----------------------------------------
# ewaff EWAS script - 1 of 2 - no cell count adjustment
# Sarah changed to meffil 19/06/24
# because the model has changed to dnam ~ AD
# and the meffil summary report will be useful for analysts to easily check inflation etc
#
# Please look through the script to change:
# file names and paths
# batch variable names
# phenotype name
#
# ----------------------------------------

## pkgs
library(tidyverse) # tidy code and data
library(usefunc) # own package of useful functions
library(meffil)
ncores <- detectCores()
ncores <- ncores-4
options(mc.cores=ncores)
print("n cores =")
print(ncores)

# loading and saving files
# please use either the command args or file name options below -
# whichever you're more comfortable with

# args <- commandArgs(trailingOnly = TRUE)
# phen_file <- args[1]
# meth_file <- args[2]
# svs_file <- args[3] # make this "" if not using
# covar_file <- args[4]
# out_file <- args[5]
# pheno_summary_table_outfile <- args[6]
# pheno_summary_table_outfile_2 <- args[7]
# pheno_summary_table_outfile_3 <- args[8]
# ewas_report_file <- args[9]

phen_file <- "pheno_and_samplesheet-ahrr.tsv"
meth_file <- "clean-meth.RData"
svs_file <- "ad-svs.tsv"
covar_file <- "covars-no-cc-svs.txt"
out_file <- "AD-ewas-no-cc.Robj"
# output files - please replace [cohort_name] with your cohort name!
summary_stats_file <- "[cohort_name]_AD_ewas_summary_stats_no_cc.csv"
pheno_summary_table_outfile_1 <- "[cohort_name]_AD-ewas-pheno_summary_numeric_vars_meanSD_no_cc.csv"
pheno_summary_table_outfile_2 <- "[cohort_name]_AD-ewas-pheno_summary_numeric_vars_ttest_no_cc.csv"
pheno_summary_table_outfile_3 <- "[cohort_name]_AD-ewas-pheno_summary_categorical_vars_counts_no_cc.csv"
ewas_report_file <- "[cohort_name]_AD_ewas_report_no_cc.html"

## read in data
# change column "IID" to your own identifier
phen_dat <- read_tsv(phen_file)
phen_dat$IID <- as.character(phen_dat$IID)
meth <- new_load(meth_file)

# ----------------------------------------
# EWAS functions
# ----------------------------------------

prep_pheno_data <- function(phen, pheno_dat, svs_file, IID, covs)
{
  if (svs_file == "") return(na.omit(pheno_dat[, c(IID, phen, covs)]))
  ## read in svs and cell counts
  svs <- read_tsv(svs_file)
  svs$IID <- as.character(svs$IID)
  sv_nam <- grep("sv", colnames(svs), value = T)
  # we don't need this next line as SVs are already in the covs file
  #covs <- c(covs, sv_nam)
  # Prepare phenotype data
  temp_phen <- pheno_dat %>%
    left_join(svs, by = setNames(IID, IID)) %>%
    dplyr::select(one_of(IID), one_of(phen), one_of(covs)) %>%
    na.omit(.)
  
  return(temp_phen)
}

run_ewas <- function(phen, pheno_dat, svs_file, meth_dat, IID, out_file, covs) 
{
  # prep pheno data
  temp_phen <- prep_pheno_data(phen, pheno_dat, svs_file, IID, covs)

  sv_nam <- grep("sv", colnames(temp_phen), value = T)
  # we don't need this next line as SVs are already in the covs file
  #covs <- c(covs, sv_nam)
  covs <- covs[covs %in% colnames(temp_phen)]
  print("covs:")
  print(covs)
  
  # Match meth to Pheno
  temp_meth <- meth_dat[, na.omit(match(temp_phen[[IID]], colnames(meth_dat)))]
  temp_phen <- temp_phen[match(colnames(temp_meth), temp_phen[[IID]]), ]
  
  # Get cases and controls per probe
  # Please edit if your phenotype is not 0=control, 1=case
  cases <- temp_phen[temp_phen[[phen]] == 1, IID, drop = TRUE]
  cases <- as.character(cases)
  
  n_cases <- rowSums(!is.na(temp_meth[, cases]))
  n_controls <- rowSums(!is.na(temp_meth[, !colnames(temp_meth) %in% cases]))
  probe_n_cc <- as_tibble(cbind(probeID = names(n_cases), N_cases = n_cases, N_controls = n_controls))
  
  ################
  
  ## new section - 19/06/24
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
  
  model <- as.formula(paste0(phen, " ~ ", paste(c("methylation", covs), collapse = " + ")))
  
  print(head(temp_phen))
  print(dim(temp_phen))
  print(dim(temp_meth))
  temp_phen <- as.data.frame(temp_phen)
  print(temp_phen[,phen][1:5])
  print(length(temp_phen[,phen]))
  
  print(head(temp_phen[,covs]))
  print(dim(temp_phen[,covs]))
  temp_meth <- as.matrix(temp_meth)
  
  # Run EWAS using meffil
  # changed to meffil by Sarah 19/06/24
  message("Running EWAS")
  # please change sig.threshold to 2.4e-7 if your data is 450k
  ewas.parameters <- meffil.ewas.parameters(sig.threshold=9e-8,  ## EWAS p-value threshold
                                            max.plots=10, ## plot at most 100 CpG sites
                                            qq.inflation.method="median",  ## measure inflation using median
                                            model="all") ## select default EWAS model; 
  
  obj <- meffil.ewas(temp_meth,variable=temp_phen[,phen],covariates=temp_phen[,covs], isva=F, random.seed=23)
  print("ewas has run")
  save(obj, file=out_file)
  ewas.summary<-meffil.ewas.summary(obj,temp_meth,parameters=ewas.parameters)                              
  meffil.ewas.report(ewas.summary, output.file=ewas_report_file)
  
  print("ewas is saved")
}

# ----------------------------------------
# Run the EWAS
# ----------------------------------------

# please change to your variable name
phen <- "ad"

covariates <- readLines(covar_file)

# please change "IID" to your ID name
run_ewas(phen = phen, 
         pheno_dat = phen_dat, 
         svs_file = svs_file,
         meth_dat = meth,
         IID = "IID", 
         out_file = out_file, 
         covs = covariates)

ewas_results <- obj$analyses$all$table
ewas_results <- ewas_results[order(ewas_results$p.value),]
print(head(ewas_results))
write.csv(ewas_results, file=summary_stats_file,quote = F)
# please change threshold value to 2.4e-7 if your data is 450k
ewas_results <- ewas_results[ewas_results$p.value<9e-08,]
print("N of associations:")
print(dim(ewas_results))
print("Total N of participants:")
print(obj$analyses$all$table$n[1])


print("FIN")