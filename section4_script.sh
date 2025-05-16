#!/bin/bash

module load R/4.4

working_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
source_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_scripts/section4/"

echo "create_candida_auris_orthologs"
Rscript --vanilla ${source_directory}create_candida_auris_orthologs.R $working_directory

echo "create_ecoli_orthologs"
Rscript --vanilla ${source_directory}create_ecoli_orthologs.R $working_directory

#AUCell scoring
echo "A549 scoring"
Rscript --vanilla ${source_directory}../aucell_scoring_and_threshold.R $working_directory "sciplex_data/A549_processed_filtered" "RNA" "yeast_human_orthologs_up" $SLURM_CPUS_PER_TASK &
Rscript --vanilla ${source_directory}../aucell_scoring_and_threshold.R $working_directory "sciplex_data/A549_processed_filtered" "RNA" "ecoli_human_orthologs_up" $SLURM_CPUS_PER_TASK &

echo "K562 scoring"
Rscript --vanilla ${source_directory}../aucell_scoring_and_threshold.R $working_directory "sciplex_data/K562_processed_filtered" "RNA" "yeast_human_orthologs_up" $SLURM_CPUS_PER_TASK &
Rscript --vanilla ${source_directory}../aucell_scoring_and_threshold.R $working_directory "sciplex_data/K562_processed_filtered" "RNA" "ecoli_human_orthologs_up" $SLURM_CPUS_PER_TASK &

echo "MCF7 scoring"
Rscript --vanilla ${source_directory}../aucell_scoring_and_threshold.R $working_directory "sciplex_data/MCF7_processed_filtered" "RNA" "yeast_human_orthologs_up" $SLURM_CPUS_PER_TASK &
Rscript --vanilla ${source_directory}../aucell_scoring_and_threshold.R $working_directory "sciplex_data/MCF7_processed_filtered" "RNA" "ecoli_human_orthologs_up" $SLURM_CPUS_PER_TASK

echo "figure 4a"
Rscript --vanilla ${source_directory}figure_4a.R $working_directory

echo "figure 4b"
Rscript --vanilla ${source_directory}figure_4b.R $working_directory

echo "figure 4c"
Rscript --vanilla ${source_directory}figure_4c.R $working_directory

echo "figure 4d"
Rscript --vanilla ${source_directory}figure_4d.R $working_directory

echo table S4
Rscript --vanilla ${source_directory}table_S4.R $working_directory
