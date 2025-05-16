#!/bin/bash

module load R/4.4

working_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
source_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_scripts/section2/"


echo "create ITH MP term2gene"
Rscript --vanilla ${source_directory}create_mp_t2g.R $working_directory

echo "global RAC DE"
Rscript --vanilla ${source_directory}global_rac_de.R $working_directory

echo "create global RAC signatures and ranks"
Rscript --vanilla ${source_directory}create_global_rac_signatures_and_ranks.R $working_directory

echo "figure S2a Global RAC enrichment"
Rscript --vanilla ${source_directory}figure_S2a.R $working_directory

echo "find all cluster markers"
# Rscript --vanilla ${source_directory}cluster_de.R $working_directory

echo "create all clusters signatures and ranks"
# Rscript --vanilla ${source_directory}create_cluster_signatures_and_ranks.R $working_directory

echo "figure 2a"
Rscript --vanilla ${source_directory}figure_2a.R $working_directory

echo "create supercluster signatures"
Rscript --vanilla ${source_directory}create_supercluster_signatures_rra.R $working_directory

echo "figure 2b"
Rscript --vanilla ${source_directory}figure_2b.R $working_directory

echo "figure 2c"
Rscript --vanilla ${source_directory}figure_2c.R $working_directory

echo "figure 2d"
Rscript --vanilla ${source_directory}figure_2d.R $working_directory

echo "figure S2b"
Rscript --vanilla ${source_directory}figure_S2b.R $working_directory

echo "figure 2e"
Rscript --vanilla ${source_directory}figure_2e.R $working_directory

echo "create rac and supercluster signature spreadsheet"
Rscript --vanilla ${source_directory}table_S2.R $working_directory



