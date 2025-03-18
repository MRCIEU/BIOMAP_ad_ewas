# BIOMAP_ad_ewas
V2 of BIOMAP atopic dermatitis EWAS pipeline. Builds on the original (V1) pipeline by Tom Battram (https://github.com/thomasbattram/ad-cse)

# *Please read these instructions before running the pipeline*

Before running this pipeline, you will need to have normalised and cleaned your methylation data (in these scripts the file containing the methylation matrix is named clean-meth.RData), and prepared a file with your phenotype data plus methylation covariates such as batch information (in these scripts this is named ad-data-cleaned.tsv). 

For each R script, you will need to look through and check where you need to edit/change/add variable names. In V2 scripts we have indicated within the script exactly where you need to do so.

# Please follow the pipeline as follows:

# 1. Estimating cell counts
Please follow the *original pipeline* to estimate 12-cell counts. Theer are two versions: one uses IDAT files, and the other uses a methylatin matrix (if you do not have IDAT files). The scripts and instructions can be found here:
https://github.com/thomasbattram/ad-cse/tree/main/scripts-for-external-ewas/estimating-cellcounts

The 'Estimating cell counts' R script will produce a file that is either named extended-blood-celltypes-idats.tsv or extended-blood-celltypes-epidish.tsv, depending on which version you create.

These cell counts should now be added to the phenotype file (ad-data-cleaned.tsv) - you need to write your own script to do this.

Finally, once you have merged the cell count and phenotype data, you need to create two vectors of covariate names - one without cell counts (this will probably just be age and sex) and one with those covariates plus cell count names. Write these to two text files - covars-no-cc.txt and covars-cc.txt

# 2. Generate surrogate variables 
Surrogate variables (SVs) will be used to control for batch and any unmeasured variation in the methylation data that is not related to our variable of interest or covariates. Because we are going to run the EWAS twice, once with and once without adjusting for cell counts, this script regresses out cell counts from the methylation data before generating the SVs.

Please run the script gen_svs.R in this repository (https://github.com/MRCIEU/BIOMAP_ad_ewas/blob/main/scripts/gen_svs.R)

This script will produce a file with the SVs (currently named ad-svs.tsv), and will add the names of the SVs to the text files containing the covariate names. 

It will also create a smoking variable from methylation level at the AHRR locus (cg05575921), as we also want to control for smoking in the EWAS. We are doing this using DNAm to avoid the issue of missing self-report smoking data, and to standardise it across cohorts. The script will add the AHRR variable to the phenotype dataframe, and write out a new phenotype file (pheno_and_samplesheet-ahrr.tsv).

Please check the heatmap produced by this script for association between SVs and other covariates.

Before running the EWAS, please double check these covariate files to make sure they contain the correct covariates!

# 3. Run conventional EWAS
The EWAS will be run twice - once not adjusting for cell counts (3a), and once adjusting for cell counts (3b).

# 3a. EWAS not adjusting for cell counts
Please run script conventional-ewas-AD-no-cc.R in this repository (https://github.com/MRCIEU/BIOMAP_ad_ewas/blob/main/scripts/conventional-ewas-AD-no-cc.R).

In the output file names at the start of the script, please add you cohort name in place of [cohort_name]

This script will produce the following 5 output files that you need to send to us:
- [cohort_name]_AD-ewas-pheno_summary_numeric_vars_meanSD_no_cc.csv
- [cohort_name]_AD-ewas-pheno_summary_numeric_vars_ttest_no_cc.csv
- [cohort_name]_AD-ewas-pheno_summary_categorical_vars_counts_no_cc.csv
- [cohort_name]_AD_ewas_report_no_cc.html
- [cohort_name]_AD_ewas_summary_stats_no_cc.csv

# 3b. EWAS adjusting for cell counts
Please run script conventional-ewas-AD-cc.R in this repository (https://github.com/MRCIEU/BIOMAP_ad_ewas/blob/main/scripts/conventional-ewas-AD-cc.R)

In the output file names at the start of the script, please add you cohort name in place of [cohort_name]

This script will produce the following 5 output files that you need to send to us:
- [cohort_name]_AD-ewas-pheno_summary_numeric_vars_meanSD_cc.csv
- [cohort_name]_AD-ewas-pheno_summary_numeric_vars_ttest_cc.csv
- [cohort_name]_AD-ewas-pheno_summary_categorical_vars_counts_cc.csv
- [cohort_name]_AD_ewas_report_cc.html
- [cohort_name]_AD_ewas_summary_stats_cc.csv

# 4. Run variance EWAS
The EWAS will be run twice - once not adjusting for cell counts (4a), and once adjusting for cell counts (4b).

# 3a. Variance EWAS not adjusting for cell counts
Please run script var-ewas-AD-no-cc.R in this repository (https://github.com/MRCIEU/BIOMAP_ad_ewas/blob/main/scripts/var-ewas-AD-no-cc.R).

This script will produce the following 4 output files that you need to send to us:
- [cohort_name]_var-ewas-pheno_summary_numeric_vars_meanSD_no-cc.csv
- [cohort_name]_var-ewas-pheno_summary_numeric_vars_ttest_no-cc.csv
- [cohort_name]_var-ewas-pheno_summary_categorical_vars_counts_no-cc.csv
- [cohort_name]_varewas-res-no-cc.tsv


# 3b. Variance EWAS adjusting for cell counts
Please run script var-ewas-AD-cc.R in this repository (xxx).

This script will produce the following 4 output files that you need to send to us:
- [cohort_name]_var-ewas-pheno_summary_numeric_vars_meanSD-cc.csv
- [cohort_name]_var-ewas-pheno_summary_numeric_vars_ttest-cc.csv
- [cohort_name]_var-ewas-pheno_summary_categorical_vars_counts-cc.csv
- [cohort_name]_varewas-res-cc.tsv
