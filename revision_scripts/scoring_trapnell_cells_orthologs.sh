#!/bin/bash

module load R/4.4

working_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
source_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_scripts/"

#AUCell scoring
echo "A549 scoring"
Rscript --vanilla ${source_directory}aucell_scoring_and_threshold.R $working_directory "sciplex_data/A549_processed_filtered" "RNA" "yeast_human_orthologs_up" $SLURM_CPUS_PER_TASK &
# Rscript --vanilla ${source_directory}aucell_scoring_and_threshold.R $working_directory "sciplex_data/A549_processed_filtered" "RNA" "ecoli_human_orthologs_up" $SLURM_CPUS_PER_TASK &

echo "K562 scoring"
Rscript --vanilla ${source_directory}aucell_scoring_and_threshold.R $working_directory "sciplex_data/K562_processed_filtered" "RNA" "yeast_human_orthologs_up" $SLURM_CPUS_PER_TASK &
# Rscript --vanilla ${source_directory}aucell_scoring_and_threshold.R $working_directory "sciplex_data/K562_processed_filtered" "RNA" "ecoli_human_orthologs_up" $SLURM_CPUS_PER_TASK &

echo "MCF7 scoring"
Rscript --vanilla ${source_directory}aucell_scoring_and_threshold.R $working_directory "sciplex_data/MCF7_processed_filtered" "RNA" "yeast_human_orthologs_up" $SLURM_CPUS_PER_TASK 
# Rscript --vanilla ${source_directory}aucell_scoring_and_threshold.R $working_directory "sciplex_data/MCF7_processed_filtered" "RNA" "ecoli_human_orthologs_up" $SLURM_CPUS_PER_TASK 