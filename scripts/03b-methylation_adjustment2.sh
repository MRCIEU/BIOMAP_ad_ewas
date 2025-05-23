#!/bin/bash

set -e

# Initialize variables
config_file="./config"

# Parse options using getopts
while getopts "c:" opt; do
    case $opt in
        c) config_file=$OPTARG ;;
        *) echo "Usage: $0 -c <config_file>"
           exit 1 ;;
    esac
done

# Shift option arguments, so $1 becomes the first positional argument
shift $((OPTIND - 1))

set -e
echo "-----------------------------------------------"
echo ""
echo "Using config located at:" ${config_file}
echo ""
echo "-----------------------------------------------"
	
source ${config_file}
batch_number=${1}
exec &> >(tee ${section_03b_logfile}${batch_number})
print_version

if [ -n "${1}" ]
then
	re='^[0-9]+$'
	if ! [[ $batch_number =~ $re ]] ; then
		echo "error: Batch variable is not a number"
		exit 1
	fi
	i=${1}
	echo "Running batch ${i} of ${meth_chunks}"
else
	i="NA"
	echo "Running entire set on a single node using ${nthreads} threads."
fi


if [ "${related}" = "no" ]
then
	echo "You have specified that the data is not family data. Adjusting only for covariates."
	echo "Generating normally transformed residuals"
	${R_directory}Rscript resources/methylation/adjust_covs.R \
		${methylation_no_outliers} \
		${covariates_combined}.txt \
		${transformed_methylation_adjusted} \
		transformed \
		${nthreads} \
		${meth_chunks} \
		${i} \
        	${methylation_array}

	echo "Generating untransformed residuals"
	${R_directory}Rscript resources/methylation/adjust_covs.R \
		${methylation_no_outliers} \
		${covariates_combined}.txt \
		${untransformed_methylation_adjusted} \
		untransformed \
		${nthreads} \
		${meth_chunks} \
		${i} \
        	${methylation_array}

elif [ "${related}" = "yes" ]
then
	# For family data adjust methylation data for relatedness (take residuals after fitting pedigree matrix, i.e. GRAMMAR method)
	echo "You have specified that the data is family data. Adjusting for pedigree and covariates..."
	#echo "Generating normally transformed residuals"
	#${R_directory}Rscript resources/methylation/adjust_pedigree.R \
	#	${home_directory} \
	#	${methylation_no_outliers} \
	#	${grmfile_relateds} \
	#	${covariates_combined}.txt \
	#	${transformed_methylation_adjusted} \
	#	transformed \
	#	${nthreads} \
	#	${meth_chunks} \
	#	${i} \
        #	${methylation_array}
    
	echo "Generating untransformed residuals"
	${R_directory}Rscript resources/methylation/adjust_pedigree_notransform.R \
		${home_directory} \
		${methylation_no_outliers} \
		${grmfile_relateds} \
		${covariates_combined}.txt \
		${untransformed_methylation_adjusted} \
		untransformed \
		${nthreads} \
		${meth_chunks} \
		${i} \
       		${methylation_array}
fi

echo "Successfully adjusting covariates"
