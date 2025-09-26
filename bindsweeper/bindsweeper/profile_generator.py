#!/usr/bin/env python3
"""Nextflow profile generation for parameter sweeps."""

import re
from typing import Any

from .parameter_converters import get_converter


def generate_profile_name(
    mode: str, swept_params: dict[str, Any], prefix: str = ""
) -> str:
    """Generate a profile name based on mode and swept parameters."""
    components = ["bindsweeper"]
    if prefix:
        components.append(prefix.rstrip("_"))
    components.append(mode)

    # Add swept parameter values to name
    for param_name, value in sorted(swept_params.items()):
        converter = get_converter(param_name)
        formatted_value = converter.format_value_for_name(value)

        # For parameters with custom formatters that include the parameter name, use formatted value directly
        if param_name in ["hotspots", "rfd_hotspots"] and formatted_value.startswith("hotspots"):
            components.append(formatted_value)
        else:
            # Shorten parameter names for readability
            param_short = param_name.replace("rfd_", "r").replace("fampnn_", "f").replace("mpnn_", "m").replace("af2_", "a").replace("boltz_", "b").replace("_", "")
            components.append(f"{param_short}{formatted_value}")

    # Join with underscores and ensure valid profile name
    profile_name = "_".join(components)
    # Replace any remaining invalid characters
    profile_name = re.sub(r"[^a-zA-Z0-9_]", "_", profile_name)

    return profile_name


def generate_profile_content(
    mode: str, all_params: dict[str, Any], swept_params: dict[str, Any], profile_name: str = None
) -> str:
    """Generate nextflow.config profile configuration content."""
    if profile_name is None:
        profile_name = generate_profile_name(mode, swept_params)

    # Convert all parameters to profile format
    profile_params = {}
    for param_name, value in all_params.items():
        converter = get_converter(param_name)
        profile_params.update(converter.to_profile_param(param_name, value))

    # Build profile content
    lines = [f"    {profile_name} {{", "        params {"]

    # Add mode first
    lines.append(f"            rfd_mode = '{mode}'")

    # Add all other parameters
    for key, value in sorted(profile_params.items()):
        if key == "rfd_mode":
            continue  # Already added

        # Format value appropriately
        if isinstance(value, str):
            if value.startswith("[") and value.endswith("]"):
                # Already formatted list
                lines.append(f'            {key} = "{value}"')
            else:
                lines.append(f'            {key} = "{value}"')
        elif isinstance(value, bool):
            lines.append(f"            {key} = {str(value).lower()}")
        elif isinstance(value, (int, float)):
            lines.append(f"            {key} = {value}")
        elif isinstance(value, list):
            # Format list
            list_str = "[" + ", ".join(f"'{item}'" for item in value) + "]"
            lines.append(f'            {key} = "{list_str}"')
        elif value is None:
            lines.append(f"            {key} = null")
        else:
            lines.append(f"            {key} = '{value}'")

    lines.extend(["        }", "    }"])

    return "\n".join(lines) + "\n"




def append_profiles_to_bindsweeper_config(
    profiles: list[str], config_path: str = "bindsweeper.config", dry_run: bool = False
) -> None:
    """Append generated profiles to bindsweeper.config file."""
    if dry_run:
        print("\nWould append profiles to bindsweeper.config:")
        for profile in profiles:
            print(profile)
        return

    try:
        # Check if config file exists, create if it doesn't
        import os

        if not os.path.exists(config_path):
            # Create new config file with profiles section
            content = "profiles {\n" + "".join(profiles) + "}\n"
        else:
            # Read existing file
            with open(config_path) as f:
                content = f.read()

            # Check if profiles section exists
            match = re.search(r"profiles\s*{\s*", content)
            if not match:
                # No profiles section, add it
                content = content.rstrip() + "\n\nprofiles {\n" + "".join(profiles) + "}\n"
            else:
                # Find the closing brace of the profiles section
                brace_count = 1
                search_pos = match.end()

                while search_pos < len(content) and brace_count > 0:
                    if content[search_pos] == '{':
                        brace_count += 1
                    elif content[search_pos] == '}':
                        brace_count -= 1
                    search_pos += 1

                # Insert profiles before the closing brace
                insert_pos = search_pos - 1  # Position of closing brace
                content = content[:insert_pos] + "".join(profiles) + content[insert_pos:]

        # Write to file
        with open(config_path, "w") as f:
            f.write(content)

        print(f"Successfully appended {len(profiles)} profiles to {config_path}")

    except Exception as e:
        raise ValueError(f"Failed to append profiles to {config_path}: {str(e)}")


def insert_profiles_into_config(
    config_path: str, profiles: list[str], dry_run: bool = False
) -> None:
    """Insert generated profiles into the nextflow.config file."""
    if dry_run:
        print("\nWould insert profiles into nextflow.config:")
        for profile in profiles:
            print(profile)
        return

    try:

        with open(config_path) as f:
            content = f.read()

        # Find the profiles section
        match = re.search(r"profiles\s*{\s*", content)
        if not match:
            raise ValueError("Profiles section not found in nextflow.config")

        # Calculate insertion position (after the opening brace)
        insert_pos = match.end()

        # Insert profiles
        new_content = (
            content[:insert_pos] + "\n" + "".join(profiles) + content[insert_pos:]
        )

        # Write back to file
        with open(config_path, "w") as f:
            f.write(new_content)

        print(f"Successfully inserted {len(profiles)} profiles into nextflow.config")

    except Exception as e:
        raise ValueError(f"Failed to insert profiles into nextflow.config: {str(e)}")
