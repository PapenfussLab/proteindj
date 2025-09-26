# BindSweeper Developer Documentation

BindSweeper is a parameter sweeping tool for the ProteinDJ pipeline. It allows systematic exploration of parameter combinations to optimize protein design outcomes.

## Architecture Overview

### Core Components

- **CLI (`cli.py`)**: Command-line interface and main orchestration
- **SweepConfig (`sweep_config.py`)**: Configuration parsing and validation
- **SweepEngine (`sweep_engine.py`)**: Parameter combination generation and execution
- **ProfileGenerator (`profile_generator.py`)**: Nextflow profile generation
- **ResultsProcessor (`results_processor.py`)**: Results aggregation and analysis

### Key Classes

#### SweepConfig
Loads and validates sweep configuration from YAML files. Supports:
- Mode-specific validation for binder and non-binder designs
- Parameter type validation (strings, numbers, enums, patterns)
- Schema validation against Nextflow parameter definitions

#### SweepEngine
Orchestrates the parameter sweep process:
- Generates cartesian product of sweep parameters
- Creates Nextflow profiles for each combination
- Executes pipeline runs in parallel or sequentially
- Handles quick test functionality

#### ProfileGenerator
Converts parameter combinations into Nextflow profiles:
- Generates unique profile names based on parameters
- Formats parameters according to Nextflow syntax
- Inserts profiles into `nextflow.config`

## Development Setup

### Prerequisites
- Python 3.9+
- uv package manager (recommended) or pip
- Access to ProteinDJ pipeline

### Installation

Development mode with uv:
```bash
cd bindsweeper/
uv sync --extra test
uv run bindsweeper --help
```

Development mode with pip:
```bash
cd bindsweeper/
pip install -e .
bindsweeper --help
```

### Running Tests

```bash
# Full test suite with coverage
uv run pytest
make test

# Fast tests without coverage
uv run pytest --no-cov
make test-fast

# Specific test module
uv run pytest tests/test_sweep_config.py -v
```

### Code Quality

```bash
# Format code
uv run ruff format
make format

# Lint code
uv run ruff check
make lint

# Fix auto-fixable lint issues
uv run ruff check --fix
```

## Configuration Format

### Basic Structure
```yaml
mode: binder_denovo  # Pipeline mode

# Optional: additional profile to include (e.g., cluster profile)
profile: milton

# Fixed parameters (override nextflow.config)
fixed_params:
  rfd_num_designs: 8
  seqs_per_design: 8

# Parameters to sweep
sweep_params:
  rfd_noise_scale:
    type: range
    min: 0.0
    max: 0.5
    step: 0.25
  
  rfd_hotspots:
    values:
      - "[B208,B232,B239]"
      - "[B210,B234]"
      - "[]"

# Results processing configuration
results_config:
  zip_results: true
  output_csv: merged_best.csv
```

### Sweep Types

#### Values List
```yaml
parameter_name:
  values:
    - value1
    - value2
    - value3
```

#### Range
```yaml
parameter_name:
  type: range
  min: 0.0
  max: 1.0
  step: 0.1
```

#### Linspace
```yaml
parameter_name:
  type: linspace
  start: 0.0
  stop: 1.0
  num: 11
```

### Profile Parameter

The optional `profile` parameter allows you to include additional Nextflow profiles in the pipeline command. This is useful for specifying cluster configurations, resource profiles, or other environment-specific settings.

**Examples:**

With profile:
```yaml
profile: milton
```
Generated command: `nextflow run main.nf -profile milton,bindsweeper_binder_denovo_param1_value1`

Without profile:
```yaml
# profile parameter omitted or null
```
Generated command: `nextflow run main.nf -profile bindsweeper_binder_denovo_param1_value1`

**Common use cases:**
- Cluster-specific profiles (e.g., `mitlon`, `slurm`, `pbs`)
- Resource profiles (e.g., `large_memory`, `gpu`)
- Environment profiles (e.g., `docker`, `singularity`)

## Adding New Features

### Adding New Sweep Types

1. Define the sweep type in `sweep_types.py`:
```python
@dataclass
class CustomSweep(SweepType):
    custom_param: float
    
    def generate_values(self) -> list[Any]:
        # Implementation
        return values
```

2. Register in `create_sweep()` function:
```python
def create_sweep(definition: dict) -> SweepType:
    sweep_type = definition.get("type", "values")
    if sweep_type == "custom":
        return CustomSweep(**definition)
    # ... existing types
```

### Adding Parameter Converters

For complex parameter formatting, add converters in `parameter_converters.py`:

```python
class CustomParameterConverter(ParameterConverter):
    def to_profile_param(self, name: str, value: Any) -> dict[str, Any]:
        # Convert value to Nextflow format
        return {name: formatted_value}
    
    def format_value_for_name(self, value: Any) -> str:
        # Format for profile naming
        return str(value)

# Register converter
PARAMETER_CONVERTERS["custom_param"] = CustomParameterConverter()
```

### Extending Results Processing

Add new result processors in `results_processor.py`:

```python
def process_custom_results(self, results_dir: str) -> pd.DataFrame:
    # Custom results processing logic
    return processed_data
```

## Testing Strategy

### Unit Tests
- Test individual components in isolation
- Mock external dependencies (file system, subprocess calls)
- Cover edge cases and error conditions

### Integration Tests
- Test full workflow with sample configurations
- Use test data in `scaffolds/` directory
- Validate generated profiles and commands

### Test Structure
```
tests/
├── __init__.py
├── conftest.py          # Shared fixtures
├── test_cli.py          # CLI interface tests
├── test_sweep_config.py # Configuration parsing tests
├── test_sweep_engine.py # Engine logic tests
└── test_sweep_types.py  # Sweep type tests
```

## Common Development Tasks

### Adding a New Mode
1. Update mode validation in `binder_validator.py` or `sweep_config.py`
2. Add mode-specific parameter validation
3. Test with sample configuration
4. Update documentation

### Debugging Parameter Issues
1. Use `--debug` flag for detailed logging
2. Check parameter validation in `sweep_config.py:validate_param_value()`
3. Verify parameter conversion in `parameter_converters.py`
4. Test profile generation with `--dry-run`

### Performance Optimization
- Use parallel execution where possible
- Cache expensive operations
- Profile with larger parameter sets
- Monitor memory usage during long sweeps

## File Structure

```
bindsweeper/
├── bindsweeper/
│   ├── __init__.py
│   ├── __main__.py
│   ├── cli.py                   # Main CLI interface
│   ├── sweep_config.py          # Configuration handling
│   ├── sweep_engine.py          # Core sweep logic
│   ├── sweep_types.py           # Parameter sweep types
│   ├── profile_generator.py     # Nextflow profile generation
│   ├── parameter_converters.py  # Parameter formatting
│   ├── results_processor.py     # Results analysis
│   ├── binder_validator.py      # Binder-specific validation
│   ├── binder_schema.json       # Binder parameter schema
│   └── utils.py                 # Utility functions
├── tests/                       # Test suite
├── scaffolds/                   # Test scaffold data
├── pyproject.toml              # Project configuration
├── Makefile                    # Development commands
└── README.md                   # User documentation
```

## Contribution Guidelines

1. **Code Style**: Follow PEP 8, use ruff for formatting
2. **Testing**: Add tests for new functionality
3. **Documentation**: Update docstrings and user documentation
4. **Type Hints**: Use type hints for all new code
5. **Error Handling**: Provide clear error messages
6. **Backward Compatibility**: Avoid breaking changes to configuration format

## Debugging Tips

### Common Issues

1. **Profile Generation Failures**
   - Check parameter formatting in `parameter_converters.py`
   - Verify Nextflow config syntax
   - Test with `--dry-run` flag

2. **Configuration Validation Errors**
   - Review parameter schemas
   - Check YAML syntax
   - Validate against expected types

3. **Pipeline Execution Failures**
   - Check Nextflow logs in output directory
   - Verify container accessibility
   - Test with smaller parameter sets

### Debugging Commands

```bash
# Dry run to see generated profiles
uv run bindsweeper --dry-run --config test_sweep.yaml --output-dir test_output

# Debug mode with detailed logging
uv run bindsweeper --debug --config test_sweep.yaml --output-dir test_output

# Quick test before full sweep
uv run bindsweeper --quick-test --config test_sweep.yaml --output-dir test_output
```