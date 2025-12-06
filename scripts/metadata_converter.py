"""
Metadata Converter Module

This module provides functionality for converting and merging metadata files
into JSONL or CSV formats. It includes a base class `MetadataConverter` and
specific methods for handling metadata operations such as checking for the
presence of `fold_id` and merging multiple JSONL files into a single CSV file.

Classes:
    MetadataConverter:
        A base class for converting metadata files to JSONL or CSV formats.

Functions:
    None (all functionality is encapsulated within the `MetadataConverter` class).

Usage:
    This module is designed to be used as a command-line tool. It accepts input
    files or directories, processes metadata, and outputs the results in the
    desired format.

Example:
    python metadata_converter.py --input_files file1.jsonl file2.jsonl --output_file output.csv
"""

import csv
from typing import Dict, Any, Iterator
from pathlib import Path
import os
import argparse
import json
import re
import logging
import pickle
import numpy as np

logging.basicConfig(level=logging.INFO)


class MetadataConverter:
    """Base class for converting metadata files to JSONL."""
    def _is_fold_id_present(self, combined_entries, fold_id_to_check):
        # Iterate through the keys of combined_entries
        for key in combined_entries.keys():
            if key[0] == fold_id_to_check:  # Check if the first element of the tuple matches the fold_id
                return True
        return False
    
    def merge_all(self, metadata_fold_file: str, metadata_fold_seq_file: str, output_file: str):
        """
        Converts two JSONL files to a combined CSV file and removes empty columns only if they exist.

        Args:
            metadata_fold_file (str): Path to a JSONL file with metadata containing only fold_id's.
            metadata_fold_seq_file (str): Path to a JSONL file with metadata containing fold_id's and seq_id's.
            output_file (str): Path to the output CSV file.
        """
        import pandas as pd
        metadata_field_names = [
            # List of all possible fields
            # ID fields
            'description','fold_id', 'seq_id',
            # RFdiffusion fields
            'rfd_sampled_mask',
            # BindCraft design fields
            'bc_length','bc_plddt','bc_rmsd_target',
            # Fold design secondary structure and RoG
            'fold_helices', 'fold_strands', 'fold_total_ss', 'fold_RoG',
            # MPNN/FAMPNN fields
            'fampnn_avg_psce','mpnn_score',
            # AF2 fields
            'af2_pae_interaction','af2_pae_overall', 'af2_pae_binder', 'af2_pae_target',
            'af2_plddt_overall', 'af2_plddt_binder', 'af2_plddt_target',
            'af2_rmsd_overall','af2_rmsd_binder_bndaln','af2_rmsd_binder_tgtaln', 'af2_rmsd_target',
            # Boltz fields
            'boltz_rmsd_overall', 'boltz_rmsd_binder','boltz_rmsd_target',
            'boltz_conf_score', 'boltz_ptm', 'boltz_ptm_interface',
            'boltz_plddt', 'boltz_plddt_interface',
            'boltz_pde', 'boltz_pde_interface',
            # PyRosetta Analysis fields
            'pr_helices','pr_strands', 'pr_total_ss','pr_RoG',
            'pr_intface_BSA','pr_intface_shpcomp',
            'pr_intface_deltaG','pr_intface_deltaGtoBSA',
            'pr_intface_hbonds','pr_intface_unsat_hbonds',
            'pr_intface_packstat','pr_TEM','pr_surfhphobics',
            'seq_ext_coef','seq_length','seq_MW','seq_pI',
            # Sequence at the end for readability, followed by time stats
            'sequence','rfd_time','bc_time','af2_time'
        ]

        try:
            combined_entries = {}
            metadata_fold_data = {}
            metadata_fold_ids = set()

            # Fold ID only metadata processing
            logging.info(f"Processing fold_id-only metadata file: {metadata_fold_file}")
            with open(metadata_fold_file, 'r', encoding='utf-8') as file:
                metadata_fold_count = 0
                for line in file:
                    if not line.strip():
                        continue
                    data = json.loads(line.strip())
                    fold_id = data.get('fold_id')
                    if fold_id is None:
                        logging.warning("Found metadata entry without fold_id, skipping")
                        continue
                    
                    if fold_id not in metadata_fold_data:
                        metadata_fold_data[fold_id] = {}
                    metadata_fold_data[fold_id].update({k: v for k, v in data.items() if k != 'fold_id'})
                    metadata_fold_count += 1
                logging.info(f"Processed {metadata_fold_count} metadata fold-only entries for {len(metadata_fold_data)} fold_ids")

            # Fold/Seq ID metadata processing
            logging.info(f"Processing fold_id + seq_id metadata file: {metadata_fold_seq_file}")
            if os.path.getsize(metadata_fold_seq_file) > 0:
                with open(metadata_fold_seq_file, 'r', encoding='utf-8') as file:
                    metadata_fold_seq_count = 0
                    merge_count = 0
                    for line in file:
                        if not line.strip():
                            continue
                        data = json.loads(line.strip())
                        fold_id = data.get('fold_id')
                        seq_id = data.get('seq_id')
                        
                        if fold_id is None:
                            logging.warning("Found metadata entry without fold_id, skipping")
                            continue

                        metadata_fold_ids.add(fold_id)
                        key = (fold_id, seq_id)
                        
                        if key not in combined_entries:
                            combined_entries[key] = {}
                            metadata_fold_seq_count += 1
                        combined_entries[key].update(data)
                        
                        # Fold-only Metadata Merging
                        if fold_id in metadata_fold_data:
                            combined_entries[key].update(metadata_fold_data[fold_id])
                            merge_count += 1

                    logging.info(f"Processed {metadata_fold_seq_count} fold_id + seq_id metadata entries")
                    logging.info(f"Merged fold-only metadata into {merge_count} entries")

            # Fold-Only Entries
            fold_only_count = 0
            for fold_id, fold_only_entry in metadata_fold_data.items():
                if fold_id not in metadata_fold_ids:
                    key = (fold_id, None)
                    # Add a description field with fold_id format
                    combined_entries[key] = {
                        'description': f'fold_{fold_id}',
                        'fold_id': fold_id, 
                        **fold_only_entry
                    }
            logging.info(f"Added {fold_only_count} fold-only entries")

            # DataFrame Creation
            logging.info("Creating DataFrame from combined entries")
            df = pd.DataFrame(list(combined_entries.values()))
            logging.debug(f"Initial DataFrame shape: {df.shape}")

            # Column Cleaning
            empty_columns = df.columns[df.isna().all()].tolist()
            if empty_columns:
                logging.info(f"Removing {len(empty_columns)} empty columns: {empty_columns}")
                df = df.drop(columns=empty_columns)
            logging.debug(f"Cleaned DataFrame shape: {df.shape}")

            # Column Ordering
            ordered_cols = [col for col in metadata_field_names if col in df.columns]
            missing_cols = set(metadata_field_names) - set(ordered_cols)
            if missing_cols:
                logging.info(f"Columns not found in metadata: {missing_cols}")
            
            df = df.reindex(columns=ordered_cols)
            logging.info(f"Final output columns: {df.columns.tolist()}")

            # File Output
            # Convert Int64 to regular int for CSV output (preserving empty cells for NA values)
            # This ensures seq_id is written as integer in CSV, not float
            if 'seq_id' in df.columns:
                # Replace NA with empty string for CSV, will be empty cell
                df['seq_id'] = df['seq_id'].apply(lambda x: '' if pd.isna(x) else int(x))
                logging.info("Converted seq_id to integer strings for CSV output")
            
            df.to_csv(output_file, index=False)
            if os.path.exists(output_file):
                logging.info(f"Successfully wrote output to {output_file} ({df.shape[0]} rows, {df.shape[1]} columns)")
            else:
                logging.error("Output file creation failed")

            return True

        except (FileNotFoundError, json.JSONDecodeError, IOError) as ex:
            logging.error(f"Fatal error during processing: {str(ex)}", exc_info=True)
            return False
        except Exception as ex:
            logging.critical(f"Fatal error during processing: {str(ex)}", exc_info=True)
            return False
       

    def convert(self, input_file: Path, output_file: Path) -> bool:
        """Convert metadata file to JSONL format.
        
        Args:
            input_file: Path to metadata file
            output_file: Path where JSONL will be written
            
        Returns:
            True if conversion succeeded, False otherwise
        """

        try:
            
            with open(output_file, 'a') as f:
                
                for record in self._parse_metadata(input_file):
                    f.write(json.dumps(record) + '\n')
            return True
        except Exception as ex:
            print(ex)
            return False
    
    def _parse_metadata(self, input_file: Path) -> Iterator[Dict[str, Any]]:
        """Parse metadata file and yield records.
        
        Args:
            input_file: Path to metadata file
            
        Yields:
            Dictionary records to be written as JSONL
        """
        # This would be implemented by subclasses
        raise NotImplementedError("Subclasses must implement this method")
    
    def convert_split_by_description(self, input_file: Path) -> bool:
        """Convert metadata file to multiple JSON files, one per description.
        
        Args:
            input_file: Path to metadata file
            
        Returns:
            True if conversion succeeded, False otherwise
        """
        try:
            for record in self._parse_metadata(input_file):
                if 'description' in record:
                    # Use the description field for the filename
                    output_file = f"{record['description']}.json"
                    with open(output_file, 'w') as f:
                        f.write(json.dumps(record) + '\n')
                else:
                    print(f"Warning: Record missing description field: {record}")
            return True
        except Exception as ex:
            print(ex)
            return False

class AF2MetadataConverter(MetadataConverter):
    def _parse_metadata(self, input_file: Path) -> Iterator[Dict[str, Any]]:
        """
        Parse AlphaFold2 score.sc file and extract metadata with 'af2_' prefix for all fields.
        
        Args:
            input_file: Path to AlphaFold2 score.sc file
            
        Yields:
            Dictionary records containing AF2 metrics and identifiers with prefixed keys.
        """
        # Define the expected format and conversions
        float_fields = {
            'pae_interaction','pae_overall', 'pae_binder', 'pae_target',
            'plddt_overall', 'plddt_binder', 'plddt_target',
            'rmsd_overall','rmsd_binder_bndaln','rmsd_binder_tgtaln', 'rmsd_target',
            'time'
        }

        with open(input_file, 'r', encoding='utf-8') as f:
            # Read header
            header = None
            for line in f:
                if line.startswith("SCORE:"):
                    header = [col.strip() for col in line.split()[1:]]
                    break
            
            if not header:
                raise ValueError("No valid header found in score.sc file")
            
            # Process data lines
            for line in f:
                if not line.startswith("SCORE:"):
                    continue
                    
                parts = line.split()
                if len(parts) != len(header) + 1:  # Account for "SCORE:" prefix
                    continue
                
                # Extract values and convert to appropriate types
                values = parts[1:]
                record = {}
                for key, value in zip(header, values):
                    if key == "description":
                        record[key] = value  # Keep 'description' unchanged
                    else:
                        prefixed_key = f"af2_{key}"
                        
                        if key in float_fields:
                            try:
                                record[prefixed_key] = round(float(value), 2)
                            except ValueError:
                                record[prefixed_key] = None
                        else:
                            record[prefixed_key] = value
                
                # Round time to seconds
                record['af2_time'] = round(float(record['af2_time'])) 
                
                # Extract fold_id and seq_id from description
                description = record.get('description', '')
                match = re.match(r'fold_(\d+)_seq_(\d+)_af2pred', description)
                if match:
                    record['fold_id'] = int(match.group(1))
                    record['seq_id'] = int(match.group(2))
                else:
                    record['fold_id'] = None
                    record['seq_id'] = None
                
                yield record

class BCMetadataConverter(MetadataConverter):
    def _parse_metadata(self, input_file: Path) -> Iterator[Dict[str, Any]]:
        """
        Parse BindCraft analysis JSON files and yield their contents directly.
        
        Args:
            input_file: Path to BindCraft analysis JSON file
            
        Yields:
            Dictionary record with pre-formatted fields from analysis script
        """
        try:
            with open(input_file, 'r') as f:
                data = json.load(f)
                yield data
        except json.JSONDecodeError as e:
            logging.error(f"Invalid JSON in {input_file}: {e}")
        except Exception as e:
            logging.error(f"Error processing {input_file}: {e}")

    def save_jsonl_file(self, input_files: list, output_file: Path) -> bool:
        """
        Convert multiple JSON files to a single JSONL file with selected metadata.
        
        Args:
            input_files: List of JSON file paths
            output_file: Path to save the JSONL file
            
        Returns:
            True if conversion succeeded, False otherwise
        """
        try:
            with open(output_file, 'w') as out_file:
                for input_file in input_files:
                    for record in self._parse_metadata(input_file):
                        json.dump(record, out_file)
                        out_file.write('\n')
            return True
        except Exception as e:
            logging.error(f"Failed to create JSONL file {output_file}: {e}")
            return False

class BoltzMetadataConverter(MetadataConverter):
    def _parse_metadata(self, input_file: Path) -> Iterator[Dict[str, Any]]:
        """
        Parse Boltz-2 prediction JSON files and yield their contents directly.
        
        Args:
            input_file: Path to Boltz-2 JSON file
            
        Yields:
            Dictionary record with pre-formatted fields from alignment script
        """
        try:
            with open(input_file, 'r') as f:
                data = json.load(f)
                yield data
        except json.JSONDecodeError as e:
            logging.error(f"Invalid JSON in {input_file}: {e}")
        except Exception as e:
            logging.error(f"Error processing {input_file}: {e}")

class FAMPNNMetadataConverter(MetadataConverter):
    def _parse_metadata(self, input_file: Path) -> Iterator[Dict[str, Any]]:
        """
        Parse a Full-Atom MPNN JSON file and extract relevant metadata.
        
        Args:
            input_file: Path to the Full-Atom MPNN JSON file
            
        Yields:
            Dictionary records containing metadata for each design
        """
        with open(input_file, 'r') as f:
            data = json.load(f)
            design = data['design']
            sequence = data['sequence']
            fampnn_avg_psce = data['fampnn_avg_psce']
            # Extract fold_id and seq_id from the design name
            match = re.match(r'fold_(\d+)_seq_(\d+)', design)
            if match:
                fold_id = int(match.group(1))
                seq_id = int(match.group(2))
            else:
                # If the pattern doesn't match, set to None or handle as needed
                fold_id = None
                seq_id = None
            
            yield {
                "fold_id": fold_id,
                "seq_id": seq_id,
                "sequence": sequence,
                "fampnn_avg_psce": float(fampnn_avg_psce)
            }

class MPNNMetadataConverter(MetadataConverter):
    def _parse_metadata(self, input_file: Path) -> Iterator[Dict[str, Any]]:
        """
        Parse a ProteinMPNN JSON file and extract relevant metadata.
        
        Args:
            input_file: Path to the ProteinMPNN JSON file
            
        Yields:
            Dictionary records containing metadata for each design
        """
        with open(input_file, 'r') as f:
            data = json.load(f)
            design = data['design']
            sequence = data['sequence']
            score = data['score']
            # Extract fold_id and seq_id from the design name
            match = re.match(r'fold_(\d+)_seq_(\d+)', design)
            if match:
                fold_id = int(match.group(1))
                seq_id = int(match.group(2))
            else:
                # If the pattern doesn't match, set to None or handle as needed
                fold_id = None
                seq_id = None
            
            yield {
                "fold_id": fold_id,
                "seq_id": seq_id,
                "sequence": sequence,
                "mpnn_score": float(score)
            }

class RFDMetadataConverter(MetadataConverter):
    def _parse_metadata(self, input_file: Path) -> Iterator[Dict[str, Any]]:
        """
        Parse a .trb file (pickle format) and extract metadata.
        
        Args:
            input_file: Path to .trb file
            
        Yields:
            Dictionary record containing metadata with keys prefixed with 'rfd_'
        """
        try:
            # Extract fold_id from filename
            filename = os.path.basename(input_file)
            match = re.search(r'fold_(\d+)', filename)
            fold_id = int(match.group(1)) if match else None
            
            with open(input_file, 'rb') as f:
                trb_data = pickle.load(f)
            
            # Create JSON-serializable record with prefixed keys
            record = {"fold_id": fold_id}
            
            for key, value in trb_data.items():
                # Prefix key with 'rfd_'
                prefixed_key = f"rfd_{key}"
                
                # Round time to whole number
                if key == 'time':
                    record[prefixed_key] = round(value)
                # Convert numpy arrays to lists
                elif isinstance(value, np.ndarray):
                    record[prefixed_key] = value.tolist()
                # Handle other types
                else:
                    # Try to make it JSON serializable
                    try:
                        json.dumps({prefixed_key: value})
                        record[prefixed_key] = value
                    except (TypeError, OverflowError):
                        # If not serializable, convert to string
                        record[prefixed_key] = str(value)
            
            yield record
            
        except Exception as e:
            logging.error(f"Error processing file {input_file}: {e}")
            raise

    def save_json_file(self, input_file: Path, output_dir: Path) -> bool:
        """
        Convert a .trb file to a JSON file with the same base name.
        
        Args:
            input_file: Path to .trb file
            output_dir: Directory to save the JSON file
            
        Returns:
            True if conversion succeeded, False otherwise
        """
        try:
            output_file = output_dir / f"{input_file.stem}.json"
            with open(output_file, 'w') as out_file:
                for record in self._parse_metadata(input_file):
                    json.dump(record, out_file, indent=2)
            return True
        except Exception as e:
            logging.error(f"Failed to convert {input_file} to JSON: {e}")
            return False

    def save_jsonl_file(self, input_files: list, output_file: Path) -> bool:
        """
        Convert multiple .trb files to a single JSONL file with selected metadata.
        
        Args:
            input_files: List of .trb file paths
            output_file: Path to save the JSONL file
            
        Returns:
            True if conversion succeeded, False otherwise
        """
        try:
            with open(output_file, 'w') as out_file:
                for input_file in input_files:
                    for record in self._parse_metadata(input_file):
                        # Extract only the fields we want for the JSONL file
                        selected_data = {
                            "fold_id": record.get("fold_id"),
                            "rfd_time": record.get("rfd_time"),
                            "rfd_sampled_mask": record.get("rfd_sampled_mask")
                        }
                        json.dump(selected_data, out_file)
                        out_file.write('\n')
            return True
        except Exception as e:
            logging.error(f"Failed to create JSONL file {output_file}: {e}")
            return False

def main():
    parser = argparse.ArgumentParser(description='Extracts metadata and converts it to json line format')
    parser.add_argument('--input_files', '-i', nargs='+',
                        help='Path to input files e.g. log file')
    parser.add_argument('--input_dir', '-d', 
                        help='Path to input directory containing a set of input files')
    parser.add_argument('--input_ext', '-e', 
                        help='Input filename extension. Only applies if using an input directory e.g. ".json"')
    parser.add_argument('--converter', '-c', default="rfd",
                        choices=['af2','bc','boltz','fampnn','mpnn','rfd'], 
                        help='Converter to use. e.g. af2, bc, boltz, fampnn, mpnn,rfd')
    parser.add_argument('--output_dir', 
                        help='Output directory path')
    parser.add_argument('--output_file', '-o', default='metadata.jsonl', 
                        help='Output jsonl file path (default: metadata.jsonl)')
    parser.add_argument('--type', '-t', default='jsonl',
                        help='format of conversion whether csv or json (default: jsonl)')
    parser.add_argument('--split_by_description', action='store_true',
                    help='Create individual JSON files for each entry using the description field')
    args = parser.parse_args()
    
    converters = {
        "af2": AF2MetadataConverter,
        "bc": BCMetadataConverter,
        "boltz": BoltzMetadataConverter,
        "fampnn": FAMPNNMetadataConverter,
        "mpnn": MPNNMetadataConverter,
        "rfd": RFDMetadataConverter
    }
    
    if args.converter not in converters:
        raise ValueError(f"Unknown converter: {args.converter}")
    selected_converter = converters[args.converter]()

    if args.converter == 'rfd' and args.input_dir:
        input_dir = Path(args.input_dir)
        output_dir = Path(args.output_dir) if args.output_dir else input_dir
        output_dir.mkdir(exist_ok=True, parents=True)
        
        extension = args.input_ext if args.input_ext.startswith('.') else '.' + args.input_ext
        input_files = list(input_dir.glob(f'*{extension}'))
        
        if not input_files:
            print(f"No {extension} files found in {input_dir}")
            return
        
        # Process individual JSON files
        success_count = 0
        for input_file in input_files:
            if selected_converter.save_json_file(input_file, output_dir):
                success_count += 1
                print(f"Converted {input_file} to {output_dir}/{input_file.stem}.json")
        
        print(f"Successfully converted {success_count} out of {len(input_files)} files to individual JSON files.")
        
        # Also create a JSONL file with selected metadata
        output_jsonl = Path(args.output_file)
        if selected_converter.save_jsonl_file(input_files, output_jsonl):
            print(f"Successfully created JSONL file with selected metadata at {output_jsonl}")
        else:
            print(f"Failed to create JSONL file at {output_jsonl}")
        
        return
    if args.converter == 'bc' and args.input_dir:
        input_dir = Path(args.input_dir)
        output_dir = Path(args.output_dir) if args.output_dir else input_dir
        output_dir.mkdir(exist_ok=True, parents=True)
        
        extension = args.input_ext if args.input_ext.startswith('.') else '.' + args.input_ext
        input_files = list(input_dir.glob(f'*{extension}'))
        
        if not input_files:
            print(f"No {extension} files found in {input_dir}")
            return
        
        # Create a JSONL file with selected metadata
        output_jsonl = Path(args.output_file)
        if selected_converter.save_jsonl_file(input_files, output_jsonl):
            print(f"Successfully created JSONL file with selected metadata at {output_jsonl}")
        else:
            print(f"Failed to create JSONL file at {output_jsonl}")
        return
    if args.split_by_description:
        logging.info(f"Converting metadata to individual JSON files using {args.converter} converter")
        if isinstance(args.input_files, list):
            input_path = Path(args.input_files[0])
        else:
            input_path = Path(args.input_files)
        success = selected_converter.convert_split_by_description(input_path)
        if success:
            print(f"Conversion successful. Individual JSON files created.")
        else:
            print("Conversion failed.")
    elif args.type == 'jsonl' and args.input_files:
        if isinstance(args.input_files, list):
            input_path = Path(args.input_files[0])
        else:
            input_path = Path(args.input_files)
        logging.info(f"Converting metadata to JSONL using {args.converter} converter")
        logging.info("Only one input file accepted for jsonl conversion: %s", input_path)
        output_path = Path(args.output_file)
        success = selected_converter.convert(input_path, output_path)
        if success:
            print(f"Conversion successful. Output written to {output_path}")
        else:
            print("Conversion failed.")
    elif args.type == 'jsonl' and args.input_dir:
        input_dir = Path(args.input_dir)
        extension = '*' + args.input_ext
        input_files = list(input_dir.glob(extension))
        if not input_files:
            print(f"No files found in {input_dir}")
            return
        for input_file in input_files:
            success = selected_converter.convert(input_file, args.output_file)
            if success:
                print(f"Converted {input_file}")
            else:
                print(f"Failed to convert {input_file}")
    elif args.type == 'csv':
        logging.info("Combining metadata to CSV")
        logging.info("Two input files accepted: metadata_fold: %s, metadata_fold_seq: %s", args.input_files[0], args.input_files[1])
        
        if args.input_files[0] and args.input_files[1]:
            success = selected_converter.merge_all(Path(args.input_files[0]),Path(args.input_files[1]), Path(args.output_file))
            if success:
                print(f"Conversion to csv successful. Output written to {args.output_file}")
            else:
                print("Conversion to csv failed.")

if __name__ == "__main__":
    main()