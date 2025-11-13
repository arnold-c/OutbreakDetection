# Deletion Plan for Unused Scripts and Functions

**Project:** OutbreakDetection
**Date:** 2025-11-13
**Goal:** Remove all unused scripts from `/scripts` directory and reorganize manuscript scripts

---

## Pre-Deletion Checklist

- [ ] Current working directory is clean (or changes are committed)
- [ ] All tests pass: `just tests`
- [ ] Manuscript builds successfully: `just manuscript`
- [ ] Created backup branch: `jj new -m "backup: before script cleanup"`

---

## Deletion Plan Overview

**Total files to delete:** 17 scripts
**Total files to move:** 8 manuscript scripts
**Estimated time:** ~45 minutes
**Risk level:** Low (all files are unreferenced or being reorganized)

---

## Phase 1: Delete High-Confidence Unused Scripts (12 files)

These scripts have NO references in the codebase.

### Step 1.1: Delete exploratory ensemble simulation scripts
**Files:** 4 files
**Rationale:** Exploratory scripts for ensemble simulations, superseded by manuscript workflow

```bash
rm scripts/ensemble-sim.jl
rm scripts/ensemble-sim_single-scenario.jl
rm scripts/ensemble-sim_noise-visualizations.jl
rm scripts/ensemble-sim_optimal-accuracy-isocline.jl
```

**JJ Revision:**
```bash
jj describe -m "refactor: remove unused ensemble simulation scripts

Remove exploratory ensemble simulation scripts that are not referenced
anywhere in the codebase. These were superseded by the manuscript
workflow in manuscript/scripts/.

Deleted files:
- scripts/ensemble-sim.jl
- scripts/ensemble-sim_single-scenario.jl
- scripts/ensemble-sim_noise-visualizations.jl
- scripts/ensemble-sim_optimal-accuracy-isocline.jl"
```

**Verification:**
```bash
# Verification will be done at end of phase
```

---

### Step 1.2: Delete exploratory diagnostic testing scripts
**Files:** 3 files
**Rationale:** Exploratory scripts for diagnostic testing, superseded by manuscript workflow

```bash
rm scripts/ensemble-diag-testing_constant-thresholds.jl
rm scripts/ensemble-diag-testing_optimal-thresholds.jl
rm scripts/ensemble-diag-testing_optimal-thresholds_single-timeseries.jl
```

**JJ Revision:**
```bash
jj describe -m "refactor: remove unused diagnostic testing scripts

Remove exploratory diagnostic testing scripts that are not referenced
anywhere in the codebase. These were superseded by the manuscript
workflow in manuscript/scripts/.

Deleted files:
- scripts/ensemble-diag-testing_constant-thresholds.jl
- scripts/ensemble-diag-testing_optimal-thresholds.jl
- scripts/ensemble-diag-testing_optimal-thresholds_single-timeseries.jl"
```

**Verification:**
```bash
# Verification will be done at end of phase
```

---

### Step 1.3: Delete single simulation scripts
**Files:** 2 files
**Rationale:** Exploratory single simulation scripts, not used by manuscript

```bash
rm scripts/single-sim.jl
rm scripts/single-sim_setup.jl
```

**JJ Revision:**
```bash
jj describe -m "refactor: remove unused single simulation scripts

Remove single simulation scripts that are not referenced anywhere in
the codebase. These were exploratory scripts not used by the manuscript.

Deleted files:
- scripts/single-sim.jl
- scripts/single-sim_setup.jl"
```

**Verification:**
```bash
# Verification will be done at end of phase
```

---

### Step 1.4: Delete superseded plotting scripts
**Files:** 2 files
**Rationale:** Superseded by manuscript/scripts/schematic-plot.jl and src/line_plots.jl

```bash
rm scripts/outbreak-detection-schematic.jl
rm scripts/schematic-plots.jl
```

**JJ Revision:**
```bash
jj describe -m "refactor: remove superseded plotting scripts

Remove plotting scripts that have been superseded by:
- manuscript/scripts/schematic-plot.jl
- src/line_plots.jl

Deleted files:
- scripts/outbreak-detection-schematic.jl
- scripts/schematic-plots.jl"
```

**Verification:**
```bash
# Verification will be done at end of phase
# Verification will be done at end of phase
```

---

### Step 1.5: Delete duplicate line-plots script
**Files:** 1 file
**Rationale:** Functionality moved to src/line_plots.jl

```bash
rm scripts/line-plots.jl
```

**JJ Revision:**
```bash
jj describe -m "refactor: remove duplicate line-plots script

Remove scripts/line-plots.jl as functionality has been moved to
src/line_plots.jl which is actively used by the manuscript.

Deleted files:
- scripts/line-plots.jl"
```

**Verification:**
```bash
# Verification will be done at end of phase
# Verification will be done at end of phase
```

---

### Step 1.6: Delete utility calculation script
**Files:** 1 file
**Rationale:** Not referenced anywhere, likely one-off calculation

```bash
rm scripts/calculate-dynamical-noise-vaccintion-rates.jl
```

**JJ Revision:**
```bash
jj describe -m "refactor: remove unused vaccination rate calculation script

Remove one-off calculation script that is not referenced anywhere in
the codebase.

Deleted files:
- scripts/calculate-dynamical-noise-vaccintion-rates.jl"
```

**Verification:**
```bash
# Verification will be done at end of phase
```

---

### Step 1.7: Delete scratch file
**Files:** 1 file
**Rationale:** Scratch/temporary file

```bash
rm scripts/scratch.jl
```

**JJ Revision:**
```bash
jj describe -m "refactor: remove scratch file

Remove temporary scratch file from scripts directory.

Deleted files:
- scripts/scratch.jl"
```

**Verification:**
```bash
# Verification will be done at end of phase
```

---

## Phase 2: Delete Medium-Confidence Scripts (3 files)

These scripts are referenced in `@static if false` blocks but are disabled.

### Step 2.1: Delete disabled development scripts
**Files:** 3 files
**Rationale:** Scripts were disabled and are no longer needed

```bash
rm scripts/single-sim_plots.jl
rm scripts/ensemble-diag-testing_scenarios_plots.jl
rm scripts/ensemble-sim_optimal-accuracy-lineplot.jl
```

**Note:** `scripts/debugging.jl` doesn't exist, so skip it.

**JJ Revision:**
```bash
jj describe -m "refactor: remove disabled development scripts

Remove development scripts that were previously disabled via
@static if false. These scripts are no longer needed as the
manuscript workflow is complete.

Deleted files:
- scripts/single-sim_plots.jl
- scripts/ensemble-diag-testing_scenarios_plots.jl
- scripts/ensemble-sim_optimal-accuracy-lineplot.jl"
```

**Verification:**
```bash
# Verification will be done at end of phase
```

---

## Phase 3: Delete threshold-optimization.jl

### Step 3.1: Delete threshold-optimization.jl
**Files:** 1 file
**Rationale:** Functionality is fully covered by manuscript scripts

```bash
rm scripts/threshold-optimization.jl
```

**JJ Revision:**
```bash
jj describe -m "refactor: remove threshold optimization script

Remove threshold-optimization.jl as functionality is fully covered
by the manuscript workflow in manuscript/scripts/.

Deleted files:
- scripts/threshold-optimization.jl"
```

**Verification:**
```bash
# Verification will be done at end of phase
```

---

## Phase 4: Reorganize Manuscript Scripts

### Step 4.1: Move manuscript scripts to main scripts directory
**Files:** 8 files moved
**Rationale:** Consolidate all analysis scripts in one location

```bash
# Move all manuscript scripts to scripts/manuscript/
mkdir -p scripts/manuscript
mv manuscript/scripts/optimal-thresholds.jl scripts/manuscript/
mv manuscript/scripts/optimal-thresholds_loading.jl scripts/manuscript/
mv manuscript/scripts/optimal-thresholds_plots.jl scripts/manuscript/
mv manuscript/scripts/optimal-thresholds_checks.jl scripts/manuscript/
mv manuscript/scripts/optimal-thresholds_tables.jl scripts/manuscript/
mv manuscript/scripts/supplemental_plots.jl scripts/manuscript/
mv manuscript/scripts/supplemental_tables.jl scripts/manuscript/
mv manuscript/scripts/plotting-setup.jl scripts/manuscript/
mv manuscript/scripts/schematic-plot.jl scripts/manuscript/

# Remove now-empty scripts directory
rmdir manuscript/scripts
```

**JJ Revision:**
```bash
jj describe -m "refactor: consolidate manuscript scripts into scripts/manuscript

Move all manuscript generation scripts from manuscript/scripts/ to
scripts/manuscript/ to consolidate all analysis scripts in one location.

This improves project organization by:
- Keeping all executable scripts in the scripts/ directory
- Maintaining manuscript/ for manuscript source files only
- Making it clearer which scripts generate manuscript outputs

Moved files:
- manuscript/scripts/optimal-thresholds.jl -> scripts/manuscript/
- manuscript/scripts/optimal-thresholds_loading.jl -> scripts/manuscript/
- manuscript/scripts/optimal-thresholds_plots.jl -> scripts/manuscript/
- manuscript/scripts/optimal-thresholds_checks.jl -> scripts/manuscript/
- manuscript/scripts/optimal-thresholds_tables.jl -> scripts/manuscript/
- manuscript/scripts/supplemental_plots.jl -> scripts/manuscript/
- manuscript/scripts/supplemental_tables.jl -> scripts/manuscript/
- manuscript/scripts/plotting-setup.jl -> scripts/manuscript/
- manuscript/scripts/schematic-plot.jl -> scripts/manuscript/"
```

**Verification:**
```bash
# Update any references to the old paths
rg "manuscript/scripts/" --type julia
```

---

### Step 4.2: Update script references in moved files
**Files:** Multiple files modified
**Rationale:** Fix internal references after moving scripts

**Edit the following files to update paths:**

1. **`scripts/manuscript/optimal-thresholds.jl`**
   - Change `manuscriptdir("scripts", args...)` to `projectdir("scripts", "manuscript", args...)`
   - Update all `manuscript_scripts()` calls to use new path

2. **`scripts/manuscript/optimal-thresholds_loading.jl`**
   - Update `include(srcdir("ensemble-parameters.jl"))` if needed

3. **Any other files with cross-references**

**JJ Revision:**
```bash
jj describe -m "refactor: update paths in manuscript scripts

Update internal path references in manuscript scripts to reflect
their new location in scripts/manuscript/.

Modified files:
- scripts/manuscript/optimal-thresholds.jl
- scripts/manuscript/optimal-thresholds_loading.jl
- [any other files with path updates]"
```

**Verification:**
```bash
# Verification will be done at end of phase
# Verification will be done at end of phase
```

---

### Step 4.3: Update OutbreakDetection.jl to include manuscript scripts
**Files:** 1 file (modified)
**Rationale:** Enable LSP support for manuscript scripts by including them in the module

**Edit:** `src/OutbreakDetection.jl`

Replace the `@static if false` block (lines 59-64) with:
```julia
@static if false
    # Include manuscript scripts for LSP support
    # These scripts are not loaded at runtime but help the LSP
    # recognize plotting functions and provide autocomplete
    include("../scripts/manuscript/optimal-thresholds.jl")
    include("../scripts/manuscript/optimal-thresholds_loading.jl")
    include("../scripts/manuscript/optimal-thresholds_plots.jl")
    include("../scripts/manuscript/optimal-thresholds_checks.jl")
    include("../scripts/manuscript/optimal-thresholds_tables.jl")
    include("../scripts/manuscript/supplemental_plots.jl")
    include("../scripts/manuscript/supplemental_tables.jl")
    include("../scripts/manuscript/plotting-setup.jl")
    include("../scripts/manuscript/schematic-plot.jl")
end
```

**JJ Revision:**
```bash
jj describe -m "refactor: update LSP includes for manuscript scripts

Update @static if false block to include all manuscript scripts
in their new location. This enables LSP support for these scripts
by making the language server aware of the plotting functions and
other exports from the OutbreakDetection module.

The scripts are not loaded at runtime (due to @static if false)
but provide autocomplete and type checking support during development.

Modified files:
- src/OutbreakDetection.jl"
```

**Verification:**
```bash
# Verification will be done at end of phase
```

---

### Step 4.4: Create manuscript README
**Files:** 1 file created
**Rationale:** Document how to generate manuscript outputs

**Create:** `manuscript/README.md`

```markdown
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
julia --project=. ../scripts/manuscript/optimal-thresholds.jl
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
julia --project=. ../scripts/manuscript/schematic-plot.jl
```

**Main manuscript plots:**
```bash
julia --project=. ../scripts/manuscript/optimal-thresholds_plots.jl
```

**Supplemental plots:**
```bash
julia --project=. ../scripts/manuscript/supplemental_plots.jl
```

**Tables:**
```bash
julia --project=. ../scripts/manuscript/supplemental_tables.jl
```

**Validation checks:**
```bash
julia --project=. ../scripts/manuscript/optimal-thresholds_checks.jl
```

### Building the Manuscript

After generating all outputs, build the manuscript PDF:

```bash
# Verification will be done at end of phase
```

Or manually:
```bash
typst compile manuscript.typ
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
- **`supplemental_plots.jl`** - Generates supplemental plots
- **`supplemental_tables.jl`** - Generates all tables
- **`plotting-setup.jl`** - Makie theme configuration
- **`schematic-plot.jl`** - Generates schematic figure
```

**JJ Revision:**
```bash
jj describe -m "docs: add manuscript generation README

Add comprehensive README to manuscript/ directory explaining:
- How to run the optimization and plotting scripts
- Directory structure and output locations
- Step-by-step instructions for reproducing outputs
- Troubleshooting common issues
- Description of each script's purpose

This makes it clear that analysis scripts are now in
scripts/manuscript/ and provides clear instructions for
reproducing manuscript outputs.

Added files:
- manuscript/README.md"
```

**Verification:**
```bash
# Verify README is clear and accurate
cat manuscript/README.md
```

---

## Phase 5: Final Verification

### Step 5.1: Run full test suite
```bash
# Verification will be done at end of phase
```

**Expected result:** All tests pass

---

### Step 5.2: Build manuscript
```bash
# Verification will be done at end of phase
```

**Expected result:** Manuscript builds successfully with all figures

---

### Step 5.3: Verify no broken references
```bash
rg "manuscript/scripts/" --type julia
```

**Expected result:** No matches except in `@static if false` block

```bash
rg "scripts/(ensemble-sim|ensemble-diag|single-sim|outbreak-detection|schematic-plots|line-plots|calculate-dynamical)" --type julia
```

**Expected result:** No matches (or only in this deletion plan)

---

### Step 5.4: Check directory structure
```bash
ls -la scripts/
ls -la scripts/manuscript/
ls -la manuscript/
```

**Expected contents:**

`scripts/`:
- `manuscript/` directory with all manuscript generation scripts

`scripts/manuscript/`:
- `optimal-thresholds.jl`
- `optimal-thresholds_loading.jl`
- `optimal-thresholds_plots.jl`
- `optimal-thresholds_checks.jl`
- `optimal-thresholds_tables.jl`
- `supplemental_plots.jl`
- `supplemental_tables.jl`
- `plotting-setup.jl`
- `schematic-plot.jl`

`manuscript/`:
- `README.md` (new)
- `manuscript_files/`
- `supplemental_files/`
- `*.typ` files
- `OD.bib`
- No `scripts/` directory

---

### Step 5.5: Verify LSP support
**Rationale:** Ensure the LSP can provide autocomplete for manuscript scripts

1. Open `scripts/manuscript/optimal-thresholds_plots.jl` in your editor
2. Type `line_plot(` and verify autocomplete shows the function signature
3. Hover over `line_plot` and verify documentation appears
4. Check that no errors are shown for OutbreakDetection module functions

**Expected result:** Full LSP support with autocomplete and documentation

---

## Phase 6: Documentation Updates

### Step 6.1: Update main README if needed
**File:** `README.md`

Check if README references any deleted scripts or old manuscript/scripts/ path.
Update if necessary.

**JJ Revision (if changes made):**
```bash
jj describe -m "docs: update README after script reorganization

Update README to reflect:
- Removal of unused scripts from scripts/ directory
- New location of manuscript scripts in scripts/manuscript/
- Reference to manuscript/README.md for generation instructions"
```

---

### Step 6.2: Update DEPENDENCY_ANALYSIS.md
**File:** `DEPENDENCY_ANALYSIS.md`

Add a note at the top indicating cleanup was completed.

**JJ Revision:**
```bash
jj describe -m "docs: mark dependency analysis as completed

Add completion note to DEPENDENCY_ANALYSIS.md indicating that
the cleanup was performed on $(date +%Y-%m-%d).

Also note the reorganization of manuscript scripts to
scripts/manuscript/ directory."
```

---

## Rollback Plan

If anything breaks:

```bash
# View recent changes
jj log -r 'trunk()..@' --limit 15

# Abandon current change and go back to before cleanup
jj abandon @

# Or restore specific files
jj restore --from <revision-id> scripts/
jj restore --from <revision-id> manuscript/scripts/
```

---

## Summary

**Total revisions:** 14-15 (depending on documentation updates)
**Files deleted:** 17 scripts
**Files moved:** 9 (8 scripts + directory)
**Files created:** 1 (manuscript/README.md)
**Files modified:** 2-3 (OutbreakDetection.jl, possibly main README, path updates)

**Estimated time breakdown:**
- Phase 1: 15 minutes (7 revisions)
- Phase 2: 3 minutes (1 revision)
- Phase 3: 2 minutes (1 revision)
- Phase 4: 12 minutes (4 revisions - move, update paths, update LSP includes, create README)
- Phase 5: 5 minutes (verification)
- Phase 6: 3 minutes (documentation)
- **Total: ~40 minutes**

---

## Post-Cleanup Benefits

1. **Cleaner repository:** Only actively used scripts remain
2. **Better organization:** All analysis scripts in scripts/, manuscript sources in manuscript/
3. **Reduced confusion:** Clear separation between code and manuscript
4. **Easier maintenance:** Fewer files to maintain and update
5. **Better onboarding:** New contributors see only relevant code
6. **Clear documentation:** manuscript/README.md explains how to reproduce outputs
7. **Smaller repository:** Reduced disk space and clone time
8. **LSP support:** Full autocomplete and type checking for manuscript scripts

---

## Notes

- All deletions and moves are tracked in jj history and can be recovered if needed
- Each phase can be done independently
- Tests should pass after each step
- The manuscript should build successfully after all changes
- The new manuscript/README.md provides clear instructions for future use
- The `@static if false` block in OutbreakDetection.jl enables LSP support without runtime overhead
