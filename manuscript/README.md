# Manuscript Generation

This directory contains the manuscript source files and generated outputs.

## Directory Structure

- `manuscript_files/` - Generated plots and tables for the main manuscript
- `supplemental_files/` - Generated plots and tables for the supplemental appendix
- `*.typ` - Typst manuscript source files
- `OD.bib` - Bibliography

## Generating Manuscript Outputs

All scripts for generating manuscript figures and tables are located in
`../scripts/manuscript/`.

### Prerequisites

1. Ensure all dependencies are installed:
   ```bash
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```

2. Run the test suite to verify everything works:
   ```bash
   just tests
   ```

### Running the Analysis

#### Option 1: Generate all outputs (recommended)

Run the main manuscript script which orchestrates all analysis:

```bash
julia --project=. scripts/manuscript/optimal-thresholds.jl
```

This will:
1. Run threshold optimizations (if not already cached)
2. Generate all main manuscript plots
3. Generate all supplemental plots
4. Create all tables
5. Run validation checks

**Note:** The optimization step can take several hours on first run. Results
are cached in `_research/` directory.

#### Option 2: Generate specific outputs

You can run individual scripts for specific outputs:

**Schematic plot:**
```bash
julia --project=. scripts/manuscript/schematic-plot.jl
```

**Main manuscript plots:**
```bash
julia --project=. scripts/manuscript/optimal-thresholds_plots.jl
```

**Supplemental plots:**
```bash
julia --project=. scripts/manuscript/supplemental_plots.jl
```

**Tables:**
```bash
julia --project=. scripts/manuscript/supplemental_tables.jl
```

**Validation checks:**
```bash
julia --project=. scripts/manuscript/optimal-thresholds_checks.jl
```

### Building the Manuscript

After generating all outputs, build the manuscript PDF:

```bash
just manuscript
```

Or manually:
```bash
typst compile manuscript/combined-manuscript.typ
```

## Output Locations

- **Main manuscript plots:** `manuscript_files/plots/`
- **Main manuscript tables:** `manuscript_files/tables/`
- **Supplemental plots:** `supplemental_files/plots/`
- **Supplemental tables:** `supplemental_files/tables/`

## Troubleshooting

**Issue:** Optimization takes too long
**Solution:** Results are cached. If you need to re-run, delete the relevant
files in `_research/ensemble/scenario-optimizations/`

**Issue:** Plots don't appear in manuscript
**Solution:** Ensure all scripts completed successfully and check that output
files exist in the expected locations

**Issue:** Path errors when running scripts
**Solution:** Always run scripts from the project root directory using the
paths shown above

## Script Descriptions

- **`optimal-thresholds.jl`** - Main orchestrator script
- **`optimal-thresholds_loading.jl`** - Runs optimizations and loads results
- **`optimal-thresholds_plots.jl`** - Generates main manuscript plots
- **`optimal-thresholds_checks.jl`** - Validation and sanity checks
- **`optimal-thresholds_tables.jl`** - Generates main manuscript tables
- **`supplemental_plots.jl`** - Generates supplemental plots
- **`supplemental_tables.jl`** - Generates supplemental tables
- **`plotting-setup.jl`** - Makie theme configuration
- **`schematic-plot.jl`** - Generates schematic figure
