# Description

Please include a summary of the change and which issue is fixed. Please also include relevant motivation and context. List any dependencies that are required for this change.

Fixes # (issue)

## Type of change

- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] This change requires a documentation update

## Checklist

- [ ] I have made any required changes to upstream dependencies for a tool wrapper, and they are available in distribution channels (e.g. Pip, Conda).
- [ ] If I have updated the underlying software for a tool wrapper (e.g. scanpy-scripts by changing the value of `@TOOL_VERSION@`), then I have reset all 'build' values to 0 (e.g. `@TOOL_VERSION@+galaxy0`)
- [ ] If I have updated a tool wrapper without a software change, then I have bumped the associated 'build' values (e.g. `@TOOL_VERSION@+galaxy0` `@TOOL_VERSION@+galaxy1`)  
