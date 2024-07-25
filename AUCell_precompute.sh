#!/bin/bash

module load R/4.2

working_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"

while getopts ":c:a:g:o:" option; do
	case $option in 
		c)
			cellset_file="$OPTARG"
			;;
		a)
			assay_to_use="$OPTARG"
			;;
		g)
			genesets_file="$OPTARG"
			;;
		o)
			output_file="$OPTARG"
			;;
		*)
			echo "Usage $0 [-c cellset_file] [-a assay_to_use] [-g genesets_file] [-o output_file]"
			exit 1
			;;
		esac
	done

Rscript --vanilla /data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_scripts/AUCell_precompute.R $working_directory $cellset_file $assay_to_use $genesets_file $SLURM_CPUS_PER_TASK $output_file
