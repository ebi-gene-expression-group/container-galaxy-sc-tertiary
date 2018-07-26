# Galaxy HCA Single Cell Tertiary Analysis

This is a Galaxy init container for single cell RNA-Seq tertiary analysis tools, 
to be used as part of a larger orchestration of containers within Kubernetes. The init
container includes tools, workflows and defined settings on how to give certains resources
for the different tools.

## Building this container

To build this container simply execute in the root directory of this repo:

```
docker build -t <desired-docker-owner>/galaxy-sc-init:<desired-tag> -f Dockerfile_init .
```

Then push your container to docker registry with `docker push`. 

## Using this container

Follow these [instructions](https://tertiary-workflows-docs.readthedocs.io/en/latest/running_galaxy_sc_locally.html) 
to setup your environment if you haven't (only follow until before "Normal Run". In that setup, the 
helm configuration for the galaxy-stable chart, available on this repo at `helm-configs/`, is used.
To make use of the container that you just built, in your local copy of the helm config file, set the
init part to:

```yaml
init:
  image:
    repository: <desired-docker-owner>/galaxy-sc-init
    tag: <desired-tag>
    pullPolicy: Always
  force_copy: "__venv__,__config__,__galaxy-central__,__tools__"
```

Then proceed with the instructions for the 
[normal run](https://tertiary-workflows-docs.readthedocs.io/en/latest/running_galaxy_sc_locally.html#normal-run) 
in the previously mentioned page.
