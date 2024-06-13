#!/bin/bash

module load R/4.4

working_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
source_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_scripts/section2/"

#AUCell scoring
echo "A549 scoring"
#Rscript --vanilla ${source_directory}../AUCell_precompute.R $working_directory "sciPlex_data/A549_processed_filtered" "RNA" "hallmarks" $SLURM_CPUS_PER_TASK
#Rscript --vanilla ${source_directory}../AUCell_precompute.R $working_directory "sciPlex_data/A549_processed_filtered" "RNA" "ITH_meta_programs" $SLURM_CPUS_PER_TASK

echo "K562 scoring"
#Rscript --vanilla ${source_directory}../AUCell_precompute.R $working_directory "sciPlex_data/K562_processed_filtered" "RNA" "hallmarks" $SLURM_CPUS_PER_TASK
#Rscript --vanilla ${source_directory}../AUCell_precompute.R $working_directory "sciPlex_data/K562_processed_filtered" "RNA" "ITH_meta_programs" $SLURM_CPUS_PER_TASK

echo "MCF7 scoring"
#Rscript --vanilla ${source_directory}../AUCell_precompute.R $working_directory "sciPlex_data/MCF7_processed_filtered" "RNA" "hallmarks" $SLURM_CPUS_PER_TASK
#Rscript --vanilla ${source_directory}../AUCell_precompute.R $working_directory "sciPlex_data/MCF7_processed_filtered" "RNA" "ITH_meta_programs" $SLURM_CPUS_PER_TASK

echo "figure 2a and S2a"
Rscript --vanilla ${source_directory}figure_2a.R $working_directory

echo "create gene universes"
#Rscript --vanilla ${source_directory}create_gene_universes.R $working_directory

echo "find global RAC markers"
#Rscript --vanilla ${source_directory}find_global_rac_markers.R $working_directory

echo "create ITH MP term2gene"
#Rscript --vanilla ${source_directory}create_mp_t2g.R $working_directory

echo "create global RAC signatures and ranks"
#Rscript --vanilla ${source_directory}create_global_rac_signatures_and_ranks.R $working_directory

echo "figure S2b"
Rscript --vanilla ${source_directory}figure_S2b.R $working_directory

echo "find all cluster markers"
#Rscript --vanilla ${source_directory}find_cluster_markers.R $working_directory

echo "create all clusters signatures and ranks"
#Rscript --vanilla ${source_directory}create_cluster_signatures_and_ranks.R $working_directory

echo "figure 2b"
Rscript --vanilla ${source_directory}figure_2b.R $working_directory

echo "create supercluster signatures"
#Rscript --vanilla ${source_directory}create_supercluster_signatures.R $working_directory

echo "figure 2c"
Rscript --vanilla ${source_directory}figure_2c.R $working_directory

echo "figure 2d"
Rscript --vanilla ${source_directory}figure_2d.R $working_directory

echo "figure S2c"
Rscript --vanilla ${source_directory}figure_S2c.R $working_directory

echo "figure 2e"
Rscript --vanilla ${source_directory}figure_2e.R $working_directory

echo "figure S2d"
Rscript --vanilla ${source_directory}figure_S2d.R $working_directory