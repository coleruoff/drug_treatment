#!/bin/bash

module load R/4.4

working_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
source_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_scripts/section4/"

echo "figure 4a"
Rscript --vanilla ${source_directory}figure_4a.R $working_directory

echo "figure 4b"
Rscript --vanilla ${source_directory}figure_4b.R $working_directory

echo "figure 4c"
Rscript --vanilla ${source_directory}figure_4c.R $working_directory

echo "figure 4d"
Rscript --vanilla ${source_directory}figure_4d.R $working_directory

echo "figure 4e"
Rscript --vanilla ${source_directory}figure_4e.R $working_directory