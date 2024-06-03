#!/bin/bash

while getopts ":c:" option; do
	case $option in 
		c) 
			cell_line="$OPTARG"
			;;

		esac
	done

module load R/4.2

Rscript --vanilla /data/ruoffcj/projects/drug_treatment/source/cell_cycle_regression.R $cell_line
