#!/bin/env bash

docker_file_dev=Dockerfile_init_dev
docker_file_prod=Dockerfile_init_prod

from=quay.io/bgruening/galaxy-init:dev
version=1.1.0
software_version=19.05-dev

tee $docker_file_dev $docker_file_prod 1> /dev/null <<EOF
FROM $from

MAINTAINER EBI Expression Atlas Team <gene-expression@ebi.ac.uk>

LABEL Description="Hinxton Single Cell Interactive Analysis Environment"
LABEL software="Galaxy Init"
LABEL software.version="$software_version"
LABEL version="$version"


RUN virtualenv venv-workflow && \
    /bin/bash -c "source venv-workflow/bin/activate && \
                  pip install 'pip>=8.1' && \
                  pip install urllib3[secure] && \
                  pip install ephemeris==0.9.0 && \
                  deactivate && \
                  virtualenv --relocatable venv-workflow"

COPY config/job_conf.xml config/job_conf.xml
COPY config/container_resolvers_conf.xml config/container_resolvers_conf.xml
COPY config/sanitize_whitelist.txt config/sanitize_whitelist.txt

# Reqs/limits
COPY config/job_resource_params_conf.xml config/job_resource_params_conf.xml
COPY config/phenomenal_tools2container.yaml config/phenomenal_tools2container.yaml
COPY rules/k8s_destinations.py /galaxy-central/lib/galaxy/jobs/rules/k8s_destination.py

COPY html/welcome_embl_ebi_rgb_2018.png welcome/welcome_embl_ebi_rgb_2018.png
COPY html/welcome_sanger_RGB_Full_Colour_landscape.png welcome/welcome_sanger_RGB_Full_Colour_landscape.png
COPY html/welcome.html welcome/welcome.html
EOF


cat >> $docker_file_dev <<EOF
COPY config/tool_conf.xml config/tool_conf.xml
COPY tools tools

COPY workflows workflows_to_load
COPY post-start-actions.sh post-start-actions.sh
EOF

cat post-start-actions.sh post-start-actions-prod-part.sh > post-start-actions-prod.sh

cat >> $docker_file_prod <<EOF
COPY prod-workflows workflows_to_load
COPY post-start-actions-prod.sh post-start-actions.sh
COPY prod-tools.yaml prod-tools.yaml
EOF
