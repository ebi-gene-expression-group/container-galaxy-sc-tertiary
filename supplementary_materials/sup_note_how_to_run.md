# Running Hinxton Singe Cell Galaxy Setup

Running the galaxy setup with our tools and workflows can be done in various ways, from easier to more complex (but more powerful):

- Directly through https://humancellatlas.usegalaxy.eu/. This is a Galaxy instance that is ready and hosted by the Frieburg Galaxy group on a DEN.BI funded HPC system with access to x cores and y GB of RAM through m nodes as of the time of this writing.
- On an existing Galaxy instance where you have administrative rights or where you can ask administrators to install tools that you need. Ask the administrator to follow (these instructions)[] to install our tools and workflows there.
- On your own machine through Kubernetes. This entails having Kubernetes installed on your laptop/desktop machine. This will allow you to run most tools, but of course performance is degraded compared to deployment in an HPC or cloud solution due to having less hardware available (less CPUs, RAM, disk, etc).
- Deploying on a cloud provider through Kubernetes. This will vary from cloud provider to cloud provider, but the essential steps are:
  - Deploy a Kubernetes cluster
  - Deploy a shared file system that is accessible from that Kubernetes cluster
  - Deploy the Galaxy helm chart against that cluster.
  - Deploy the tools as per instructions.
  As an example of the above, these instructions allow to deploy the setup in AWS.
