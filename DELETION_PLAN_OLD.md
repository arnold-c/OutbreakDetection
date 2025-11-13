# Deletion Plan for Unused Scripts and Functions

**Project:** OutbreakDetection
**Date:** 2025-11-13
**Goal:** Remove all unused scripts from `/scripts` directory including medium-confidence files

---

## Pre-Deletion Checklist

- [ ] Current working directory is clean (or changes are committed)
- [ ] All tests pass: `just tests`
- [ ] Manuscript builds successfully: `just manuscript`
- [ ] Created backup branch: `jj new -m "backup: before script cleanup"`

---

## Deletion Plan Overview

**Total files to delete:** 17 scripts
**Estimated time:** ~30 minutes
**Risk level:** Low (all files are unreferenced)

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
just tests
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
just tests
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
just tests
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
just tests
just manuscript
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
just tests
just manuscript
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
just tests
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
just tests
```

---

## Phase 2: Delete Medium-Confidence Scripts (5 files)

These scripts are referenced in `@static if false` blocks but are disabled.

### Step 2.1: Remove disabled script includes from OutbreakDetection.jl
**Files:** 1 file (modified)
**Rationale:** Clean up disabled code block before deleting the scripts

**Edit:** `src/OutbreakDetection.jl`

Remove lines 59-64:
```julia
@static if false
    include("../scripts/single-sim_plots.jl")
    include("../scripts/ensemble-diag-testing_scenarios_plots.jl")
    include("../scripts/ensemble-sim_optimal-accuracy-lineplot.jl")
    include("../scripts/debugging.jl")
end
```

**JJ Revision:**
```bash
jj describe -m "refactor: remove disabled script includes

Remove @static if false block that included disabled development
scripts. These scripts will be deleted in subsequent commits.

Modified files:
- src/OutbreakDetection.jl"
```

**Verification:**
```bash
just tests
```

---

### Step 2.2: Delete disabled development scripts
**Files:** 4 files
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
just tests
```

---

## Phase 3: Delete threshold-optimization.jl

### Step 3.1: Delete threshold-optimization.jl
**Files:** 1 file

**Options:**
1. **Delete:** If functionality as fully covered by manuscript scripts

```bash
rm scripts/threshold-optimization.jl
```

**JJ Revision:**
```bash
jj describe -m "refactor: delete threshold optimization script as superseded

Delete threshold-optimization.jl as it has been superseded by the manuscript
scripts

Deleted files:
- scripts/threshold-optimization.jl
```

**Verification:**
```bash
just tests
```

***Update plan here***

---

## Phase 4: Final Verification

### Step 4.3: Verify no broken references
```bash
rg "scripts/(ensemble-sim|ensemble-diag|single-sim|outbreak-detection|schematic-plots|line-plots|calculate-dynamical)" --type julia
```

**Expected result:** No matches (or only in this deletion plan)

---

### Step 4.4: Check scripts directory
```bash
ls -la scripts/
```

**Expected contents:**
- `examples/` directory (if kept threshold-optimization)
- Possibly other files not in our deletion list

---

## Phase 5: Documentation Updates

### Step 5.1: Update README if needed
**File:** `README.md`

Check if README references any deleted scripts. Update if necessary.

**JJ Revision (if changes made):**
```bash
jj describe -m "docs: update README after script cleanup

Update README to reflect removal of unused scripts from scripts/
directory."
```

---

### Step 5.2: Update DEPENDENCY_ANALYSIS.md
**File:** `DEPENDENCY_ANALYSIS.md`

Add a note at the top indicating cleanup was completed.

**JJ Revision:**
```bash
jj describe -m "docs: mark dependency analysis as completed

Add completion note to DEPENDENCY_ANALYSIS.md indicating that
the cleanup was performed on [DATE]."
```

---

## Rollback Plan

If anything breaks:

```bash
# View recent changes
jj log -r 'trunk()..@' --limit 10

# Abandon current change and go back to before cleanup
jj abandon @

# Or restore specific files
jj restore --from <revision-id> scripts/
```

---

## Summary

**Total revisions:** 10-11 (depending on documentation updates)
**Files deleted:** 18 scripts
**Files modified:** 1-2 (OutbreakDetection.jl, possibly README)

**Estimated time breakdown:**
- Phase 1: 15 minutes (7 revisions)
- Phase 2: 5 minutes (2 revisions)
- Phase 3: 3 minutes (1 revision, optional)
- Phase 4: 5 minutes (verification)
- Phase 5: 2 minutes (documentation)
- **Total: ~30 minutes**

---

## Post-Cleanup Benefits

1. **Cleaner repository:** Only actively used scripts remain
2. **Reduced confusion:** No ambiguity about which scripts to use
3. **Easier maintenance:** Fewer files to maintain and update
4. **Better onboarding:** New contributors see only relevant code
5. **Smaller repository:** Reduced disk space and clone time

---

## Notes

- All deletions are tracked in jj history and can be recovered if needed
- Each phase can be done independently
- Tests should pass after each step
- The manuscript should build successfully after plotting script deletions
