#!/usr/bin/env bash

write_rows_for_modules() {
  if [[ $cli != "None" ]]; then
    local cli_link="[$cli](https://github.com/ebi-gene-expression-group/$cli/tree/master)"
  else
    local cli_link="None"
  fi
  for mod in $(ls $path/*.xml | grep -v macro); do
    local desc=$(xpath -q -e /tool/description $mod | sed -E 's+</?description>++g')
    local mod_id=$(xpath -q -e /tool/@id $mod | sed -E 's+(id=)?"++g' | sed 's/ *//g' | tr '[:upper:]' '[:lower:]')
    local mod_instance_link=https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2F$mod_id%2F$mod_id
    local mod_ts_link=https://toolshed.g2.bx.psu.edu/view/ebi-gxa/$mod_id
    echo "| [$mod_id]($mod_instance_link)<sup>[TS]($mod_ts_link) | $desc | $cli_link | ${mappings[$mod_id]} |" >> S2_modules_aut.md
  done
}

# Read mappings to analysis areas
declare -A mappings
while IFS= read -r line; do
    mappings[${line%%:*}]=${line#*:}
done < module_analysis_mapping.txt

echo "**Table S2**: Decomposed modules contributed from the different tools. The name of the module links to an active Galaxy instance where that can be used, and TS provides a second link to the module in the Galaxy Toolshed." > S2_modules_aut.md
echo "Each module is linked to one of the cli-layers mentioned in [Table S1](S1_cli-layer.md) and assigned to one or more of the relevant analysis areas: Clustering (**C**), Differential expression/Marker detection (**DE-MD**), Trajectories (**T**), Cell type alignment (**CT**) and Dimensionality reduction (**DR**)." >> S2_modules_aut.md
# Requires perl-xml-xpath conda package for xpath
echo "" >> S2_modules_aut.md
echo "| Module | Description | cli-layer | Analysis areas |" >> S2_modules_aut.md
echo "|--------|-------------|-----------|----------------|" >> S2_modules_aut.md

tertiary_path=../tools/tertiary-analysis
# Scanpy
path=$tertiary_path/scanpy
cli=scanpy-scripts
write_rows_for_modules

# Seurat
path=$tertiary_path/seurat
cli=r-seurat-scripts
write_rows_for_modules

# SC3
path=$tertiary_path/sc3
cli=bioconductor-sc3-scripts
write_rows_for_modules

# scMap
path=$tertiary_path/scmap
cli=scmap-cli
write_rows_for_modules

# scater
path=$tertiary_path/scater
cli=bioconductor-scater-scripts
write_rows_for_modules

# Monocle3
path=$tertiary_path/monocle3
cli=monocle-scripts
write_rows_for_modules

# SCCAF
path=$tertiary_path/sccaf
cli=None
write_rows_for_modules

# SCEASY
path=$tertiary_path/sceasy
cli=None
write_rows_for_modules

# ucsc cell browser
path=$tertiary_path/ucsc-cell-browser
cli=None
write_rows_for_modules

# hca-client
path=$tertiary_path/data-hca
cli=hca-matrix-downloader
write_rows_for_modules

# atlas-client
path=$tertiary_path/data-scxa
cli=None
write_rows_for_modules

# droplet-utils
path=$tertiary_path/dropletutils
cli=dropletutils-scripts
write_rows_for_modules

# garnett
path=$tertiary_path/garnett
cli=garnett-cli
write_rows_for_modules

# scpred
path=$tertiary_path/scpred
cli=scpred-cli
write_rows_for_modules

# cell-types-analysis
path=$tertiary_path/cell-types-analysis
cli=cell-types-analysis
write_rows_for_modules
