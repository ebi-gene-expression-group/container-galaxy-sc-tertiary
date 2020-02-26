FROM quay.io/ebigxa/galaxy-k8s:20.01_200206
#FROM galaxy/galaxy-k8s:20.01-feb26

#RUN virtualenv venv-workflow && /bin/bash -c "source venv-workflow/bin/activate \
#    && pip install 'pip>=8.1' && pip install urllib3[secure] && pip install ephemeris==0.9.0 \
#    && deactivate && virtualenv --relocatable venv-workflow"

#COPY config/job_conf.xml /galaxy/server/config/cont/job_conf.xml
#COPY config/tool_conf.xml /galaxy/server/config/cont/tool_conf.xml
#COPY config/container_resolvers_conf.xml /galaxy/server/config/cont/container_resolvers_conf.xml
#COPY config/sanitize_whitelist.txt /galaxy/server/config/sanitize_whitelist.txt
#COPY tools/tertiary-analysis /galaxy/server/tools/tertiary-analysis
#COPY tools/util /galaxy/server/tools/util

COPY html/welcome_embl_ebi_rgb_2018.png /galaxy/server/static/welcome_embl_ebi_rgb_2018.png
COPY html/welcome_sanger_RGB_Full_Colour_landscape.png /galaxy/server/static/welcome_sanger_RGB_Full_Colour_landscape.png
COPY html/welcome.html /galaxy/server/static/welcome.html

# Reqs/limits
# COPY config/job_resource_params_conf.xml /galaxy/server/config/job_resource_params_conf.xml
# COPY files/phenomtools2container.yaml /galaxy/server/config/phenomenal_tools2container.yaml
COPY rules/k8s_destinations.py /galaxy/server/lib/galaxy/jobs/rules/k8s_destinations.py
COPY rules/resource_bins.yaml /galaxy/server/lib/galaxy/jobs/rules/resource_bins.yaml

# COPY prod-workflows workflows_to_load
# COPY post-start-actions-prod.sh post-start-actions.sh
