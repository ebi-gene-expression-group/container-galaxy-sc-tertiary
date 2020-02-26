from galaxy.jobs import JobDestination
from galaxy.util.logging import get_logger
import os
import yaml


# Sizes for the different set of resources to be used in Kubernetes
# CPU values are in CPU units (0.1 means 100 mili CPUs)
# Memory values are in Gigabytes.

# TODO: Read YAML file with mappings from here.
__unit_suffix = {'memory': "Gi"}

r_bins = yaml.load(open( os.path.join(
    os.path.dirname(__file__), "resource_bins.yaml" ), mode='r'))

__path_tool2container = os.path.join(
    os.path.dirname(__file__), "tools2container.yaml" )

log = get_logger(__name__)

def k8s_dispatcher(resource_params, rule_helper, no_docker_default_destination_id, tool):
    # Allocate extra time
    if not rule_helper.supports_docker(tool):
        #return JobDestination(id=no_docker_default_destination_id, params=resource_params)
        return no_docker_default_destination_id
    resource_params['docker_enabled'] = True
    return JobDestination(runner="k8s", params=resource_params)


def k8s_wrapper_tiny(resource_params, tool_id, job, rule_helper, no_docker_default_destination_id, tool):
    if not rule_helper.supports_docker(tool):
        #return JobDestination(id=no_docker_default_destination_id, params=resource_params)
        return no_docker_default_destination_id
    __read_assignments(resource_params, tool_id, job)
    # return __setup_resources(resource_params, settings=__tiny, job=job, follow_up_destination="small")
    return __setup_resources(resource_params, settings=r_bins['resource_bins']['tiny'], job=job, follow_up_destination="small")


def k8s_wrapper_small(resource_params, tool_id, job, rule_helper, no_docker_default_destination_id, tool):
    if not rule_helper.supports_docker(tool):
        #return JobDestination(id=no_docker_default_destination_id, params=resource_params)
        return no_docker_default_destination_id
    __read_assignments(resource_params, tool_id, job)
    return __setup_resources(resource_params, settings=r_bins['resource_bins']['small'], job=job, follow_up_destination="medium")


def k8s_wrapper_medium(resource_params, tool_id, job, rule_helper, no_docker_default_destination_id, tool):
    if not rule_helper.supports_docker(tool):
        #return JobDestination(id=no_docker_default_destination_id, params=resource_params)
        return no_docker_default_destination_id
    __read_assignments(resource_params, tool_id, job)
    return __setup_resources(resource_params, settings=r_bins['resource_bins']['medium'], job=job, follow_up_destination="large")


def k8s_wrapper_large(resource_params, tool_id, job, rule_helper, no_docker_default_destination_id, tool):
    if not rule_helper.supports_docker(tool):
        #return JobDestination(id=no_docker_default_destination_id, params=resource_params)
        return no_docker_default_destination_id
    __read_assignments(resource_params, tool_id, job)
    return __setup_resources(resource_params, settings=r_bins['resource_bins']['large'], job=job, follow_up_destination="xlarge")


def k8s_wrapper_xlarge(resource_params, tool_id, job, rule_helper, no_docker_default_destination_id, tool):
    if not rule_helper.supports_docker(tool):
        #return JobDestination(id=no_docker_default_destination_id, params=resource_params)
        return no_docker_default_destination_id
    __read_assignments(resource_params, tool_id, job)
    return __setup_resources(resource_params, settings=r_bins['resource_bins']['xlarge'], job=job, follow_up_destination="xxlarge")

def k8s_wrapper_xxlarge(resource_params, tool_id, job, rule_helper, no_docker_default_destination_id, tool):
    if not rule_helper.supports_docker(tool):
        #return JobDestination(id=no_docker_default_destination_id, params=resource_params)
        return no_docker_default_destination_id
    __read_assignments(resource_params, tool_id, job)
    return __setup_resources(resource_params, settings=r_bins['resource_bins']['xxlarge'], job=job)

def __read_assignments(resource_params, tool_id, job):
    stream = open(__path_tool2container, mode='r')
    mappings = yaml.load(stream)
    for mapping in mappings['assignment']:
        if tool_id in mapping['tools_id']:
            resource_params.update(mapping)
            break


def __setup_resources(resource_params, settings, job, follow_up_destination=None):
    resource_params['docker_enabled'] = True
    resource_params['resources'] = "all"
    resource_params['shared_home_dir'] = "$_GALAXY_JOB_HOME_DIR"
    __check_resource_params(resource_params, resource_type='cpu')
    __check_resource_params(resource_params, resource_type='memory')
    __merge_into_res_params(resource_params, settings, resource_type='cpu')
    __merge_into_res_params(resource_params, settings, resource_type='memory')
    job_destination=JobDestination(runner="k8s", params=resource_params, id=settings['dest_id'])
    if follow_up_destination is not None:
        job_destination['resubmit'] = [dict(
            condition="memory_limit_reached",
            destination="dynamic-k8s-"+follow_up_destination,
        )]
    job.destination_id = settings['dest_id']
    log.warning("Memory requests: "+resource_params['requests_memory'])
    log.warning("Memory limits: "+resource_params['limits_memory'])
    return job_destination


def __check_resource_params(resource_params, resource_type):
    for resource_key in ['requests_' + resource_type, 'limits_' + resource_type]:
        if resource_key not in resource_params:
            resource_params[resource_key] = 0


def __merge_into_res_params(resource_params, settings, resource_type):
    if resource_params['requests_' + resource_type] == 0 and resource_params['limits_' + resource_type] == 0:
        resource_params['requests_' + resource_type] = settings['requests_' + resource_type]
        resource_params['limits_' + resource_type] = settings['limits_' + resource_type]
    elif resource_params['requests_' + resource_type] != 0 and resource_params['limits_' + resource_type] == 0:
        resource_params['limits_' + resource_type] = max(resource_params['requests_' + resource_type],
                                                         settings['limits_' + resource_type])
    elif resource_params['requests_' + resource_type] == 0 and resource_params['limits_' + resource_type] != 0:
        resource_params['requests_' + resource_type] = min(settings['limits_' + resource_type],
                                                           resource_params['requests_' + resource_type])
    if resource_type in r_bins['unit_suffix']:
        resource_params['requests_' + resource_type] = str(resource_params['requests_' + resource_type]) + r_bins['unit_suffix'][resource_type]
        resource_params['limits_' + resource_type] = str(resource_params['limits_' + resource_type]) + r_bins['unit_suffix'][resource_type]
