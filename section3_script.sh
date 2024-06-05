#!/bin/bash

module load R/4.4

working_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
source_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_scripts/section3/"

echo "figure 3a"
Rscript --vanilla ${source_directory}figure_3a.R $working_directory

echo "create supercluster signatures"
#Rscript --vanilla ${source_directory}create_supercluster_signatures.R $working_directory

echo "figure 3b"
#Rscript --vanilla ${source_directory}figure_3b.R $working_directory

echo "figure 3c"
#Rscript --vanilla ${source_directory}figure_3c.R $working_directory

echo "figure 3d"
#Rscript --vanilla ${source_directory}figure_3d.R $working_directory

echo "figure 3e"
#Rscript --vanilla ${source_directory}figure_3e.R $working_directory