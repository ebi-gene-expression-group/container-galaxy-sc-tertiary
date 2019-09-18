#!/usr/bin/env bash

write_rows_for_modules() {
  local cli_link="[$cli](https://github.com/ebi-gene-expression-group/$cli/tree/master)"
  for mod in $(ls $path/*.xml | grep -v macro); do
    local desc=$(xpath -q -e /tool/description $mod | sed -E 's+</?description>++g')
    local mod_id=$(xpath -q -e /tool/@id $mod | sed -E 's+(id=)?"++g' | sed 's/ *//g')
    local mod_instance_link=https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2F$mod_id%2F$mod_id
    local mod_ts_link=https://toolshed.g2.bx.psu.edu/view/ebi-gxa/$mod_id
    echo "| [$mod_id]($mod_instance_link)<sup>[TS]($mod_ts_link) | $desc | [$cli]($cli_link) | ${mappings[$mod_id]} |" >> O2_modules_aut.md
  done
}

# Read mappings to analysis areas
declare -A mappings
while IFS= read -r line; do
    mappings[${line%%:*}]=${line#*:}
done < module_analysis_mapping.txt
# Requires perl-xml-xpath conda package for xpath

echo "| Module | Description | cli-layer | Analysis areas |" > O2_modules_aut.md
echo "|--------|-------------|-----------|----------------|" >> O2_modules_aut.md

tertiary_path=../tools/tertiary-analysis
# Scanpy
path=$tertiary_path/scanpy
cli=scanpy-scripts
write_rows_for_modules

# Seurat
path=$tertiary_path/seurat
cli=seurat-scripts
write_rows_for_modules

# SC3
path=$tertiary_path/sc3
cli=sc3-scripts
write_rows_for_modules

# scMap
path=$tertiary_path/scmap
cli=scmap-cli
write_rows_for_modules

# scater
path=$tertiary_path/scater
cli=scater-scripts
write_rows_for_modules

# Monocle3
path=$tertiary_path/monocle3
cli=monocle-scripts
write_rows_for_modules

# SCCAF
#path=$tertiary_path/sccaf
#cli=None
#write_rows_for_modules

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
