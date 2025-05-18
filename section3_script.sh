#!/bin/bash

module load R/4.4

working_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
source_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_scripts/section3/"


echo "Global RAC Survival Analysis"
Rscript --vanilla ${source_directory}figure_3a.R $working_directory

echo "pan cancer cox regression"
Rscript --vanilla ${source_directory}pan_cancer_cox_regression.R $working_directory

echo "figure 3b"
Rscript --vanilla ${source_directory}figure_3b.R $working_directory

echo "figure 3c"
Rscript --vanilla ${source_directory}figure_3c.R $working_directory

echo "figure 3d"
Rscript --vanilla ${source_directory}figure_3d.R $working_directory

echo "create NR and R score spreadsheet"
Rscript --vanilla ${source_directory}table_S3.R $working_directory

