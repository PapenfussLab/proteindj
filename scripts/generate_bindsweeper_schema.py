#!/usr/bin/env python3
"""Generate bindsweeper's binder_schema.json from the generated per-mode nextflow schemas.

Includes every pipeline parameter valid for each mode (excluding design_mode, which is
fixed per mode), so users can sweep any parameter the pipeline accepts rather than a
hand-curated subset. Parameters not required for a mode are marked nullable.
"""
import argparse
import json
from pathlib import Path

MODES = [
    "rfd_denovo",
    "rfd_foldcond",
    "rfd_motifscaff",
    "rfd_partialdiff",
    "bindcraft_denovo",
    "boltzgen_denovo",
    "boltzgen_motifscaff",
    "rfd4_denovo",
]

FIELDS_TO_COPY = ["description", "pattern", "enum", "minimum", "maximum"]


def build_mode_entry(mode_schema):
    required = set()
    for defn in mode_schema["definitions"].values():
        required.update(defn.get("required", []))

    entry = {}
    for defn in mode_schema["definitions"].values():
        for param, prop in defn.get("properties", {}).items():
            if param == "design_mode":
                continue
            param_type = prop["type"] if param in required else [prop["type"], "null"]
            param_schema = {"type": param_type}
            for field in FIELDS_TO_COPY:
                if field in prop:
                    param_schema[field] = prop[field]
            if prop.get("default") is not None:
                param_schema["example"] = prop["default"]
            entry[param] = param_schema
    return entry


def main():
    parser = argparse.ArgumentParser(
        description="Generate bindsweeper/bindsweeper/binder_schema.json from the generated per-mode nextflow schemas."
    )
    parser.add_argument(
        "-s", "--schemas-dir", type=str, default="schemas",
        help="Directory containing nextflow_schema_<mode>.json files (default: schemas)"
    )
    parser.add_argument(
        "-o", "--output", type=str,
        default="bindsweeper/bindsweeper/binder_schema.json",
        help="Output path for binder_schema.json"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output")
    args = parser.parse_args()

    schemas_dir = Path(args.schemas_dir)
    binder_schema = {}
    for mode in MODES:
        mode_file = schemas_dir / f"nextflow_schema_{mode}.json"
        with open(mode_file) as f:
            mode_schema = json.load(f)
        binder_schema[mode] = build_mode_entry(mode_schema)
        if args.verbose:
            print(f"{mode}: {len(binder_schema[mode])} parameters")

    out_path = Path(args.output)
    with open(out_path, "w") as f:
        json.dump(binder_schema, f, indent=2)
        f.write("\n")
    print(f"Wrote {out_path} with {len(binder_schema)} modes.")


if __name__ == "__main__":
    main()
