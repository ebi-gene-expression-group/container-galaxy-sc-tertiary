FROM quay.io/ebigxa/galaxy-k8s:20.01_200206
#FROM galaxy/galaxy-k8s:20.01-feb26

# All config files are now injected at the helm level.

# Branding
COPY html/welcome_embl_ebi_rgb_2018.png /galaxy/server/static/welcome_embl_ebi_rgb_2018.png
COPY html/welcome_sanger_RGB_Full_Colour_landscape.png /galaxy/server/static/welcome_sanger_RGB_Full_Colour_landscape.png
COPY html/welcome.html /galaxy/server/static/welcome.html

# Reqs/limits
COPY rules/k8s_destinations.py /galaxy/server/lib/galaxy/jobs/rules/k8s_destinations.py
COPY rules/resource_bins.yaml /galaxy/server/lib/galaxy/jobs/rules/resource_bins.yaml
