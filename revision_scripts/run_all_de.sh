#!/bin/bash

module load R/4.4

working_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
source_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_scripts/"

#AUCell scoring
# Rscript --vanilla ${source_directory}/section2/cluster_de.R $working_directory

Rscript --vanilla ${source_directory}section2/global_rac_de.R $working_directory

# Rscript --vanilla ${source_directory}/section1/oren_cell_lines_de.R $working_directory
