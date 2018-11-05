#!/bin/bash

source activate galaxy
conda install bioconductor-sc3-scripts

## sc3-sc3-prepare
planemo tool_init --force \
					--macros \
                    --id 'sc3-sc3-prepare' \
                    --description 'Reads a normalised scater object for SC3' \
                    --name 'Prepare SC3 object' \
                    --requirement bioconductor-sc3-scripts@0.0.1 \
                    --example_command 'sc3-sc3-prepare.R -i DATA-DIR -t N_CORES -s 1712 -o OUTPUT-OBJECT-FILE' \
                    --example_input scater_log.rds \
                    --example_output R_sc3_object.rds \
                    --test_case \
                    --cite_url 'https://github.com/ebi-gene-expression-group/bioconductor-sc3-scripts.git' \
                    --help_from_command 'sc3-sc3-prepare.R -h'
                    

# sc3-sc3-calc-transfs 
    
planemo tool_init --force \
					--macros \
                    --id 'sc3-sc3-calc-transfs' \
                    --description 'Calculate CPM from serialized scater R object' \
                    --name 'SC3 calculate transformations' \
                    --requirement bioconductor-sc3-scripts@0.0.1 \
                    --example_command 'sc3-sc3-calc-transfs.R -i R_scater_serialized.rds -o R_scater_cpm.rds' \
                    --example_input R_sc3_object.rds \
                    --example_output R_sc3_transformed.rds \
                    --test_case \
                    --cite_url 'https://github.com/ebi-gene-expression-group/bioconductor-sc3-scripts.git' \
                    --help_from_command 'sc3-sc3-calc-transfs.R -h'
                                    
   

 # sc3-sc3-estimate-k.R
 
    
 planemo tool_init --force \
					--macros \
                    --id 'sc3-sc3-estimate-k' \
                    --description 'Estimate the best K for clustering' \
                    --name 'SC3 estimate k' \
                    --requirement bioconductor-sc3-scripts@0.0.1 \
                    --example_command 'sc3-sc3-estimate-k.R -i R_sc3_transformed.rds -t sc3-estimated-k.txt -o R_scater_qc.rds ' \
                    --example_input R_sc3_transformed.rds \
                    --example_output sc3-estimated-k.txt \
                    --test_case \
                    --help_from_command 'sc3-sc3-estimate-k.R -h' \
                    --cite_url 'https://github.com/ebi-gene-expression-group/bioconductor-sc3-scripts.git'  
                    --help_from_command 'sc3-sc3-estimate-k.R -h'
                                     
                    
# sc3-sc3-calc-dists.R

planemo tool_init --force \
					--macros \
                    --id 'sc3-sc3-calc-dists' \
                    --description 'Calculate distances between cells' \
                    --name 'SC3 calculate distances' \
                    --requirement bioconductor-sc3-scripts@0.0.1 \
                    --example_command 'sc3-sc3-calc-dists.R -i R_sc3_transformed.rds -o R_sc3_distance.rds' \
                    --example_input R_sc3_transformed.rds \
                    --example_output R_sc3_distance.rds \
                    --test_case \
                    --help_from_command 'sc3-sc3-calc-dists.R'
                    --cite_url 'https://github.com/ebi-gene-expression-group/bioconductor-sc3-scripts'                    
                    --help_from_command 'sc3-sc3-calc-dists.R -h'
 
 
# sc3-sc3-kmeans.R
 
    
 planemo tool_init --force \
					--macros \
                    --id 'sc3-sc3-kmeans' \
                    --description 'Perform k-Means clustering on cells' \
                    --name 'SC3 kmeans' \
                    --requirement bioconductor-sc3-scripts@0.0.1 \
                    --example_command 'sc3-sc3-kmeans.R -i R_sc3_distance.rds -k K_size -o R_sc3_kmeans.rds ' \
                    --example_input R_sc3_distance \
                    --example_input K_size \
                    --example_output R_sc3_kmeans.rds \
                    --test_case \
                    --help_from_command 'sc3-sc3-kmeans.R' \
                    --cite_url 'https://github.com/ebi-gene-expression-group/bioconductor-sc3-scripts' 
                    
# sc3-sc3-calc-consens.R 
 
 planemo tool_init --force \
					--macros \
                    --id 'sc3-sc3-calc-consens' \
                    --description 'Calculate consensus K-means clustering' \
                    --name 'SC3 calculate consensus' \
                    --requirement bioconductor-sc3-scripts@0.0.1 \
                    --example_command 'sc3-sc3-calc-consens.R -i R_sc3_kmeans.rds -d directory_store_data -o R_sc3_consensus.rds' \
                    --example_input R_sc3_kmeans.rds \
                    --example_input directory_store_data \
                    --example_output R_sc3_consensus.rds \
                    --test_case \
                    --help_from_command 'sc3-sc3-calc-consens.R' \
                    --cite_url 'https://github.com/ebi-gene-expression-group/bioconductor-sc3-scripts' 
                    
                    
                    
# sc3-sc3-calc-biology.R

 planemo tool_init --force \
					--macros \
                    --id 'sc3-sc3-calc-biology' \
                    --description 'SC3 calculates biologically relevant clusters' \
                    --name 'SC3 calculate biology' \
                    --requirement bioconductor-sc3-scripts@0.0.1 \
                    --example_command 'sc3-sc3-calc-biology.R -i R_sc3_kmeans.rds -o R_sc3_kmeans_biology.rds' \
                    --example_input R_sc3_kmeans.rds \
                    --example_output R_sc3_kmeans_biology.rds \
                    --test_case \
                    --help_from_command 'sc3-sc3-calc-biology.R' \
                    --cite_url 'https://github.com/ebi-gene-expression-group/bioconductor-sc3-scripts' 

