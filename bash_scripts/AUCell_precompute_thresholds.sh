#!/bin/bash

module load R/4.3

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

Rscript --vanilla /data/ruoffcj/projects/drug_treatment/source/AUCell_precompute_thresholds.R $cellset_file $assay_to_use $genesets_file $SLURM_CPUS_PER_TASK $output_file
