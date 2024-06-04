#!/bin/bash

module load R/4.4

working_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
source_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_scripts/section1/"

echo "yogesh_resistant_cells_processing"
#Rscript --vanilla ${source_directory}process_yogesh_data.R $working_directory

echo "memory_seq_processing"
#Rscript --vanilla ${source_directory}memory_seq_processing.R $working_directory

echo "yogesh_resistant_vs_control_de"
#Rscript --vanilla ${source_directory}yogesh_resistant_vs_control_de.R $working_directory

echo "watermelon_pc9_processing"
#Rscript --vanilla ${source_directory}watermelon_pc9_processing.R $working_directory

echo "create_resistance_signature"
#Rscript --vanilla ${source_directory}create_resistance_signature.R $working_directory

echo "create_supp_table1"
#Rscript --vanilla ${source_directory}create_supp_table1.R $working_directory

echo "figure_1a"
#Rscript --vanilla ${source_directory}figure_1a.R $working_directory 

echo "trapnell_processing_filtering"
Rscript --vanilla ${source_directory}trapnell_processing_filtering.R $working_directory

#AUCell scoring
echo "A549 scoring"
#Rscript --vanilla ${source_directory}../AUCell_precompute_thresholds.R $working_directory "sciPlex_data/A549_processed_filtered" "RNA" "raj_watermelon_resistance_signature" $SLURM_CPUS_PER_TASK &

echo "K562 scoring"
#Rscript --vanilla ${source_directory}../AUCell_precompute_thresholds.R $working_directory "sciPlex_data/K562_processed_filtered" "RNA" "raj_watermelon_resistance_signature" $SLURM_CPUS_PER_TASK &

echo "MCF7 scoring"
#Rscript --vanilla ${source_directory}../AUCell_precompute_thresholds.R $working_directory "sciPlex_data/MCF7_processed_filtered" "RNA" "raj_watermelon_resistance_signature" $SLURM_CPUS_PER_TASK 

echo "figure_1b"
Rscript --vanilla ${source_directory}figure_1b.R $working_directory

echo "figure_1c"
Rscript --vanilla ${source_directory}figure_1c.R $working_directory

echo "figure_supp1c"
Rscript --vanilla ${source_directory}figure_supp1c.R $working_directory 

echo "creating all_data.rds"
Rscript --vanilla ${source_directory}create_all_data_rac_object.R $working_directory

echo "figure_1d"
Rscript --vanilla ${source_directory}figure_1d.R $working_directory 