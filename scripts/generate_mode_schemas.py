import argparse
import csv
import json
import copy
from pathlib import Path

def parse_csv(csv_file):
    """
    Parse the CSV into a dict of {mode: {param: value or None}}
    Only parameters with a value in the <mode>_values column will override the schema default.
    """
    with open(csv_file, newline='', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        headers = next(reader)
        # Map: mode -> (param_col, value_col)
        mode_cols = {}
        for i, h in enumerate(headers):
            if h.endswith('_parameters'):
                mode = h.replace('_parameters', '')
                value_col = i + 1
                mode_cols[mode] = (i, value_col)
        param_overrides = {mode: {} for mode in mode_cols}
        for row in reader:
            for mode, (p_col, v_col) in mode_cols.items():
                param = row[p_col].strip()
                value = row[v_col].strip() if v_col < len(row) else ''
                if param:
                    # Only set override if value is non-empty
                    param_overrides[mode][param] = value if value else None
    return param_overrides

def convert_value(value, schema_param):
    """
    Convert string value from CSV to appropriate type based on schema definition.
    """
    if value is None:
        return None
    if schema_param.get('type') == 'boolean':
        return value.lower() == 'true'
    if schema_param.get('type') == 'integer':
        try:
            return int(value)
        except Exception:
            return value
    if schema_param.get('type') == 'number':
        try:
            return float(value)
        except Exception:
            return value
    if value.lower() == 'null':
        return None
    # Try to parse as JSON (for arrays etc.)
    try:
        return json.loads(value)
    except Exception:
        return value

def build_mode_schema(main_schema, mode, overrides):
    """
    Build a mode-specific schema, applying parameter overrides only where specified.
    """
    schema = copy.deepcopy(main_schema)
    schema['title'] = f"{mode} pipeline parameters"
    schema['description'] = f"Parameters for {mode} mode"
    schema['$id'] = main_schema['$id'].rsplit('/', 1)[0] + f"/nextflow_schema_{mode}.json"

    # For each definition section, filter and override defaults as needed
    for defn in schema['definitions'].values():
        if 'properties' not in defn:
            continue
        filtered = {}
        for param, prop in defn['properties'].items():
            if param in overrides:
                prop = copy.deepcopy(prop)
                override_val = overrides[param]
                if override_val is not None:
                    prop['default'] = convert_value(override_val, prop)
                filtered[param] = prop
        defn['properties'] = filtered
        # Update required list
        if 'required' in defn:
            defn['required'] = [p for p in defn['required'] if p in filtered]
            if not defn['required']:
                del defn['required']

    # Special handling for rfd_mode
    for defn in schema['definitions'].values():
        if 'rfd_mode' in defn.get('properties', {}):
            prop = defn['properties']['rfd_mode']
            prop['description'] = "Pipeline mode."
            
            if mode == 'custom':
                # Set full enum list for custom mode
                prop['enum'] = [
                    "binder_denovo",
                    "binder_foldcond",
                    "binder_motifscaff",
                    "binder_partialdiff",
                    "monomer_denovo",
                    "monomer_foldcond",
                    "monomer_motifscaff",
                    "monomer_partialdiff"
                ]
                # Preserve default if specified in CSV
                if 'rfd_mode' in overrides and overrides['rfd_mode'] is not None:
                    prop['default'] = convert_value(overrides['rfd_mode'], prop)
            else:
                # For non-custom modes, set single enum value
                val = convert_value(overrides.get('rfd_mode', mode), prop)
                prop['enum'] = [val]
                prop['default'] = val

    return schema


def main():
    parser = argparse.ArgumentParser(
        description="Generate mode-specific Nextflow schemas from a main schema and a CSV file."
    )
    parser.add_argument(
        "-s", "--schema", type=str, default="nextflow_schema.json",
        help="Path to the main schema JSON file (default: nextflow_schema.json)"
    )
    parser.add_argument(
        "-c", "--csv", type=str, default="schemas/mode_parameters.csv",
        help="Path to the CSV file with parameter overrides (default: schemas/mode_parameters.csv)"
    )
    parser.add_argument(
        "-o", "--output", type=str, default="schemas",
        help="Directory to write generated schemas (default: schemas)"
    )
    parser.add_argument(
        "-m", "--modes", nargs="*", default=None,
        help="Optional: Only generate schemas for these modes (space separated)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable verbose output"
    )
    args = parser.parse_args()

    # Load files
    with open(args.schema) as f:
        main_schema = json.load(f)
    param_overrides = parse_csv(args.csv)

    # Filter modes if requested
    modes = args.modes if args.modes else list(param_overrides.keys())
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)

    for mode in modes:
        if mode not in param_overrides:
            print(f"Warning: mode '{mode}' not found in CSV, skipping.")
            continue
        schema = build_mode_schema(main_schema, mode, param_overrides[mode])
        out_file = output_dir / f"nextflow_schema_{mode}.json"
        with open(out_file, "w") as f:
            json.dump(schema, f, indent=2)
        if args.verbose:
            print(f"Wrote {out_file} with {len(param_overrides[mode])} parameters.")

if __name__ == "__main__":
    main()
