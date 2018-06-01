from galaxy.jobs import JobDestination
import yaml


# Sizes for the different set of resources to be used in Kubernetes
# CPU values are in CPU units (0.1 means 100 mili CPUs)
# Memory values are in Gigabytes.

# TODO: Read YAML file with mappings from here.

__tiny = {'requests_cpu': 0.1,
          'limits_cpu': 0.5,
          'requests_memory': 0.3,
          'limits_memory': 0.6
          }

__small = {'requests_cpu': 0.4,
           'limits_cpu': 0.8,
           'requests_memory': 0.5,
           'limits_memory': 0.9
           }

__medium = {'requests_cpu': 0.7,
            'limits_cpu': 2,
            'requests_memory': 0.8,
            'limits_memory': 2
            }

__large = {'requests_cpu': 1.5,
           'limits_cpu': 4,
           'requests_memory': 1.8,
           'limits_memory': 5
           }

__xlarge = {'requests_cpu': 4,
            'limits_cpu': 8,
            'requests_memory': 8,
            'limits_memory': 16
            }

__path_tool2container = "config/phenomenal_tools2container.yaml"

def k8s_wrapper(resource_params):
    # Allocate extra time
    resource_params['docker_enabled'] = True
    return JobDestination(runner="k8s", params=resource_params)


def k8s_wrapper_tiny(resource_params, tool_id):
    __read_assignments(resource_params,tool_id)
    return __setup_resources(resource_params, settings=__tiny)


def k8s_wrapper_small(resource_params, tool_id):
    __read_assignments(resource_params, tool_id)
    return __setup_resources(resource_params, settings=__small)


def k8s_wrapper_medium(resource_params, tool_id):
    __read_assignments(resource_params, tool_id)
    return __setup_resources(resource_params, settings=__medium)


def k8s_wrapper_large(resource_params, tool_id):
    __read_assignments(resource_params, tool_id)
    return __setup_resources(resource_params, settings=__large)


def k8s_wrapper_xlarge(resource_params, tool_id):
    __read_assignments(resource_params, tool_id)
    return __setup_resources(resource_params, settings=__xlarge)


def __read_assignments(resource_params, tool_id):
    stream = open(__path_tool2container, mode='r')
    mappings = yaml.load(stream)
    for mapping in mappings['assignment']:
        if tool_id in mapping['tools_id']:
            resource_params.update(mapping)
            break


def __setup_resources(resource_params, settings):
    resource_params['docker_enabled'] = True
    __check_resource_params(resource_params, resource_type='cpu')
    __check_resource_params(resource_params, resource_type='memory')
    __merge_into_res_params(resource_params, settings, resource_type='cpu')
    __merge_into_res_params(resource_params, settings, resource_type='memory')
    return JobDestination(runner="k8s", params=resource_params)


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
