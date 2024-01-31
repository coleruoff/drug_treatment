#!/bin/bash

module load R/4.2

while getopts ":c:s:" option; do
	case $option in
		c)
			cell_line="$OPTARG"
			;;

		s)
			setting_to_use="$OPTARG"

		esac
	done

Rscript --vanilla /data/ruoffcj/projects/drug_treatment/source/find_resistant_cluster_markers.R $cell_line $setting_to_use
