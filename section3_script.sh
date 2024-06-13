#!/bin/bash

module load R/4.4

working_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
source_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_scripts/section3/"

#CRISPR validation

echo "create regulon genesets"
#Rscript --vanilla ${source_directory}create_tf_regulons_genesets.R $working_directory

echo "score all cells for regulon activity"
#Rscript --vanilla ${source_directory}../AUCell_precompute.R $working_directory "sciPlex_data/A549_processed_filtered" "RNA" "tf_regulons" $SLURM_CPUS_PER_TASK
#Rscript --vanilla ${source_directory}../AUCell_precompute.R $working_directory "sciPlex_data/K562_processed_filtered" "RNA" "tf_regulons" $SLURM_CPUS_PER_TASK
#Rscript --vanilla ${source_directory}../AUCell_precompute.R $working_directory "sciPlex_data/MCF7_processed_filtered" "RNA" "tf_regulons" $SLURM_CPUS_PER_TASK

echo "calculate cluster TF activity"
#Rscript --vanilla ${source_directory}calculate_cluster_tf_activity.R $working_directory

echo "create top cluster TF lists"
#Rscript --vanilla ${source_directory}create_all_cluster_tf_lists.R $working_directory

echo "create supercluster consensus TFs"
#Rscript --vanilla ${source_directory}create_supercluster_consensus_tf_lists.R $working_directory

echo "figure 3a"
Rscript --vanilla ${source_directory}figure_3a.R $working_directory

echo "create TF gene spreadsheet"
Rscript --vanilla ${source_directory}table_S3.R $working_directory
