#!/bin/bash

module load R/4.4

working_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
source_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_scripts/section1/"

#process oren data
Rscript --vanilla ${source_directory}oren_pooled_processing.R $working_directory

# run oren DE
Rscript --vanilla ${source_directory}oren_cell_lines_de.R $working_directory

# create resistance signature
Rscript --vanilla ${source_directory}create_resistance_signature_rra.R $working_directory

# validate in oren watermelon. and plot
echo "figure_1a"
Rscript --vanilla ${source_directory}figure_1a.R $working_directory

echo "processing trapnell data"
# process trapnell data
Rscript --vanilla ${source_directory}trapnell_processing_filtering.R $working_directory

# score trapnell data
#AUCell scoring
echo "A549 scoring"
Rscript --vanilla ${source_directory}../aucell_scoring_and_threshold.R $working_directory "sciplex_data/A549_processed_filtered" "RNA" "resistance_signature" $SLURM_CPUS_PER_TASK &

echo "K562 scoring"
Rscript --vanilla ${source_directory}../aucell_scoring_and_threshold.R $working_directory "sciplex_data/K562_processed_filtered" "RNA" "resistance_signature" $SLURM_CPUS_PER_TASK &

echo "MCF7 scoring"
Rscript --vanilla ${source_directory}../aucell_scoring_and_threshold.R $working_directory "sciplex_data/MCF7_processed_filtered" "RNA" "resistance_signature" $SLURM_CPUS_PER_TASK

# plot trapnell boxplots
echo "figure_1b"
Rscript --vanilla ${source_directory}figure_1b.R $working_directory

echo "figure_1c"
# plot RAC barplots
Rscript --vanilla ${source_directory}figure_1c.R $working_directory

echo "figure_S1b"
# plot pre treatment RAC barplots
Rscript --vanilla ${source_directory}figure_S1b.R $working_directory

echo "creating all data object"
# Create all data object
Rscript --vanilla ${source_directory}create_all_data_rac_object.R $working_directory

#AUCell scoring
echo "A549 scoring"
Rscript --vanilla ${source_directory}../AUCell_precompute.R $working_directory "sciPlex_data/A549_processed_filtered" "RNA" "hallmarks" $SLURM_CPUS_PER_TASK
Rscript --vanilla ${source_directory}../AUCell_precompute.R $working_directory "sciPlex_data/A549_processed_filtered" "RNA" "ITH_meta_programs" $SLURM_CPUS_PER_TASK

echo "K562 scoring"
Rscript --vanilla ${source_directory}../AUCell_precompute.R $working_directory "sciPlex_data/K562_processed_filtered" "RNA" "hallmarks" $SLURM_CPUS_PER_TASK
Rscript --vanilla ${source_directory}../AUCell_precompute.R $working_directory "sciPlex_data/K562_processed_filtered" "RNA" "ITH_meta_programs" $SLURM_CPUS_PER_TASK

echo "MCF7 scoring"
Rscript --vanilla ${source_directory}../AUCell_precompute.R $working_directory "sciPlex_data/MCF7_processed_filtered" "RNA" "hallmarks" $SLURM_CPUS_PER_TASK
Rscript --vanilla ${source_directory}../AUCell_precompute.R $working_directory "sciPlex_data/MCF7_processed_filtered" "RNA" "ITH_meta_programs" $SLURM_CPUS_PER_TASK


echo "figure_1d"
# plot RAC hallmark and MP heatmap
Rscript --vanilla ${source_directory}figure_1d.R $working_directory

echo "figure_S1"
# Create supplementary table with genesets
Rscript --vanilla ${source_directory}table_S1.R $working_directory


