#!/usr/bin/env python

import yaml

# Load input YAML
with open("software_versions_input.yml", "r") as f:
    input_versions = yaml.safe_load(f)

# Save plain version (optional)
with open("software_versions.yml", "w") as out:
    yaml.dump(input_versions, out, sort_keys=True)

# Create MultiQC-compatible YAML
output_mqc = {
    "id": "software_versions",
    "section_name": "Software Versions",
    "section_href": "https://github.com/nf-core/tools",
    "plot_type": "html",
    "description": "Software versions collected at the end of the workflow execution.",
    "data": {}
}

# Flatten tool versions correctly
for step_name, tools in input_versions.items():
    for tool_name, version in tools.items():
        output_mqc["data"][tool_name] = str(version)

# Save MultiQC-compatible YAML
with open("software_versions_mqc.yml", "w") as out:
    yaml.dump(output_mqc, out, sort_keys=False, default_flow_style=False)
