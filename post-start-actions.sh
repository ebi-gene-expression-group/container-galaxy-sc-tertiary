#!/bin/bash

echo "Adding workflows to Galaxy for user $GALAXY_DEFAULT_ADMIN_USER from dir $WORKFLOWS_DIR"


# Waits for the instance to be up and adds the workflows in a non-blocking fashion, so that tail below proceeds
# uses ephemeris. Virtualenv needed while pre 0.8 versions of ephemeris are being used in other places
source /export/venv-workflow/bin/activate && \
  /export/venv-workflow/bin/workflow-install -g http://127.0.0.1 -v -a $GALAXY_DEFAULT_ADMIN_KEY \
      -w /export/workflows_to_load --publish_workflows && \
  deactivate \
