[![Docker Repository on Quay](https://quay.io/repository/ebigxa/galaxy-sc-tertiary/status "Docker Repository on Quay")](https://quay.io/repository/ebigxa/galaxy-sc-tertiary)

# Galaxy HCA Single Cell Tertiary Analysis

This is a Galaxy init container for single cell RNA-Seq tertiary analysis tools,
to be used as part of a larger orchestration of containers within Kubernetes. The init
container includes tools, workflows and defined settings on how to give certains resources
for the different tools.

This Galaxy flavour aims to support computational scRNA-Seq data analysis in the
context of the Human Cell Atlas (HCA) Project.

## Simple run

To run this Galaxy instance through minikube follow these [instructions](https://tertiary-workflows-docs.readthedocs.io/en/latest/running_galaxy_sc_locally.html).

## Advanced: used your own flavour of the container

### Building this container

To build this container simply execute in the root directory of this repo:

```
docker build -t <desired-docker-owner>/galaxy-sc-init:<desired-tag> -f Dockerfile_init .
```

Then push your container to docker registry with `docker push`.

### Running in minikube

Follow these [instructions](https://tertiary-workflows-docs.readthedocs.io/en/latest/running_galaxy_sc_locally.html)
to setup your environment if you haven't (only follow until before "Normal Run"). In that setup, the
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

# Building community images

If for whatever reason you are more involved with Galaxy development (or ansible-galaxy-extras roles) itself and you need some of those changes to be used within a deployment of the type described here, then you will need to re-build the community images and use those within the settings configured in the helm config described above. For that, you need to execute the script `compose/build-orchestration-images.sh` with adequate arguments (see usage of that script) within https://github.com/bgruening/docker-galaxy-stable/ repo for community images. If you are happy to use Galaxy as it is from its release versions, you don't need to use this section at all and simply stay with directions up to the previous section.

Note: currently, `build-orchestration-images.sh` is part of an ongoing [pull request](https://github.com/bgruening/docker-galaxy-stable/pull/446)
