# HSCG direct tools install

This part of the setup explains how to deploy all of our latest tools and workflows on a living Galaxy instance. This could have been achieved in different ways, but once the Galaxy instance is up, it becomes irrelevant how it was provisioned for the sake of adding our tools and workflows.

This instructions also adds third party tools that we have considered useful through the time that we have used this setup for production, exploratory analysis and trainings.

## Tools installation

Tools to be provisioned are mostly available through the main Galaxy Toolshed within our [ebi-gxa user](https://toolshed.g2.bx.psu.edu/view/ebi-gxa). There are additional third party tools that we provision as part of this process. All versions of each tool are installed, which means that workflows created at different points in time will work with the tools that they were originally created with. Normally, they user should be able to upgrade the tools in a workflow through the interface without much issues, but this design guarantees reproducibility.

Requirements:
- Conda
- Conda environment with Galaxy ephemeris. This can be achieved via:

```
conda create -n galaxy-ephemeris -c bioconda ephemeris
```

- A Galaxy API key for an administrator user in the instance where you would like to install the tools.
- Place `tool_conf_ebi-gxa.xml` inside the Galaxy config files.

To be executed in bash
```
# clone our tool definition file repository
git clone https://github.com/ebi-gene-expression-group/galaxy-gxa-tools-setup
# activate the ephemeris conda env, provided it was created as shown in the requirements
conda activate galaxy-ephemeris
cd galaxy-gxa-tools-setup
./scripts/install_tools.sh http://you-galaxy-url:port/ <GALAXY_ADMIN_API_KEY>
```

The above github repository is routinely updated with newer revisions of the tools, so if you want to keep
up-to-date, you might want to run this routinely.
