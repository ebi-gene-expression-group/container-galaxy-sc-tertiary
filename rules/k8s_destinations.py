from galaxy.jobs import JobDestination
import yaml


# Sizes for the different set of resources to be used in Kubernetes
# CPU values are in CPU units (0.1 means 100 mili CPUs)
# Memory values are in Gigabytes.

# TODO: Read YAML file with mappings from here.

__tiny = {'requests_cpu': 0.1,
          'limits_cpu': 0.5,
          'requests_memory': 0.3,
          'limits_memory': 0.6,
          'dest_id': "dynamic-k8s-tiny"
          }

__small = {'requests_cpu': 0.4,
           'limits_cpu': 0.8,
           'requests_memory': 0.5,
           'limits_memory': 0.9,
           'dest_id': "dynamic-k8s-small"
           }

__medium = {'requests_cpu': 0.7,
            'limits_cpu': 2,
            'requests_memory': 0.8,
            'limits_memory': 2,
            'dest_id': "dynamic-k8s-medium"
            }

__large = {'requests_cpu': 1.5,
           'limits_cpu': 4,
           'requests_memory': 1.8,
           'limits_memory': 5,
           'dest_id': "dynamic-k8s-large"
           }

__xlarge = {'requests_cpu': 4,
            'limits_cpu': 8,
            'requests_memory': 8,
            'limits_memory': 16,
            'dest_id': "dynamic-k8s-xlarge"
            }

__xxlarge = {'requests_cpu': 4,
             'limits_cpu': 8,
             'requests_memory': 16,
             'limits_memory': 30,
             'dest_id': "dynamic-k8s-xxlarge"
              }

__path_tool2container = "config/phenomenal_tools2container.yaml"

def k8s_dispatcher(resource_params, rule_helper, no_docker_default_destination_id, tool):
    # Allocate extra time
    if not rule_helper.supports_docker(tool):
        #return JobDestination(id=no_docker_default_destination_id, params=resource_params)
        return no_docker_default_destination_id
    resource_params['docker_enabled'] = True
    return JobDestination(runner="k8s", params=resource_params)


def k8s_wrapper_tiny(resource_params, tool_id, job):
    __read_assignments(resource_params, tool_id, job)
    return __setup_resources(resource_params, settings=__tiny, job=job, follow_up_destination="small")


def k8s_wrapper_small(resource_params, tool_id, job):
    __read_assignments(resource_params, tool_id, job)
    return __setup_resources(resource_params, settings=__small, job=job, follow_up_destination="medium")


def k8s_wrapper_medium(resource_params, tool_id, job):
    __read_assignments(resource_params, tool_id, job)
    return __setup_resources(resource_params, settings=__medium, job=job, follow_up_destination="large")


def k8s_wrapper_large(resource_params, tool_id, job):
    __read_assignments(resource_params, tool_id, job)
    return __setup_resources(resource_params, settings=__large, job=job, follow_up_destination="xlarge")


def k8s_wrapper_xlarge(resource_params, tool_id, job):
    __read_assignments(resource_params, tool_id, job)
    return __setup_resources(resource_params, settings=__xlarge, job=job, follow_up_destination="xxlarge")

def k8s_wrapper_xxlarge(resource_params, tool_id, job):
    __read_assignments(resource_params, tool_id, job)
    return __setup_resources(resource_params, settings=__xlarge, job=job)

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
