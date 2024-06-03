#!/bin/bash

module load R/4.3

while getopts ":c:" option; do
	case $option in 
		c)
			cell_line="$OPTARG"
			;;
		*)
			echo "Usage $0 [-c cell_line]"
			exit 1
			;;
		esac
	done

Rscript --vanilla /data/ruoffcj/projects/drug_treatment/source/find_cluster_markers_MAST.R $cell_line
