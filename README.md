# OutbreakDetection

This code base is using the [Julia Language](https://julialang.org/) and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) to make a reproducible scientific project named
> OutbreakDetection

## Getting Started
### Installing dependencies

To (locally) reproduce this project, do the following:

0. Download this code base.
You can do this by cloning the repository or by downloading the zip file.
Notice that raw data are not included in the git-history and may need to be downloaded independently.
1. Install Julia. I would recommend you use the [`juliaup` installer](https://github.com/JuliaLang/juliaup) as it makes it much easier to deal with multiple versions of Julia, as well as keep them up to date.
2. Open a Julia console and do:

   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> using DrWatson
   julia> @quickactivate "OutbreakDetection"
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:

```julia
using DrWatson
@quickactivate "OutbreakDetection"
```
which auto-activate the project and enable local path handling from DrWatson.

### Running the project

To run the project, you can use the Justfile.
First, you should make sure you have Just installed on your computer.
If you have macOS you can  install it using homebrew.

```bash
brew install just
```

If you are on Linux, you can install Just with your package manager.
If you are on Windows, you can use the [Chocolatey](https://chocolatey.org/) package manager, or use the ubuntu package manager on [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) (recommended as it will give you access to a linux terminal).

A full list of download options for each OS can be found at the [Just GitHub repository](https://github.com/casey/just).

Once you have a version of Just installed you can use the following terminal command to run the Justfile:

```bash
just
```

If you would like to run a specific target, you can specify it after the `just` command e.g. `just manuscript`.
To list all available tasks, you can use the command `just --list` (or open the Justfile).


## Project Structure

```bash
OutbreakDetection  main  v1.11.3
❯ tree -L 3
.
├── _research
├── archived
├── data
│   ├── cases_and_deaths_estimates_2022.csv
│   ├── CFR_2022.csv
│   ├── input-gbd_region.csv
│   └── input-populations.csv
├── Justfile
├── Manifest.toml
├── manuscript
│   ├── ARCHIVE
│   │   ├── _extensions/
│   │   ├── manuscript.qmd
│   │   └── supplemental-appendix.qmd
│   ├── combined-manuscript.pdf
│   ├── combined-manuscript.typ
│   ├── manuscript_files
│   │   ├── plots/
│   │   └── tables/
│   ├── manuscript.typ
│   ├── OD.bib
│   ├── scripts
│   │   ├── optimal-thresholds_checks.jl
│   │   ├── optimal-thresholds_loading.jl
│   │   ├── optimal-thresholds_plots.jl
│   │   ├── optimal-thresholds_tables.jl
│   │   ├── optimal-thresholds.jl
│   │   ├── plotting-setup.jl
│   │   ├── schematic-plot.jl
│   │   ├── supplemental_plots.jl
│   │   └── supplemental_tables.jl
│   ├── supplemental_files
│   │   ├── plots
│   │   └── tables
│   ├── supplemental-appendix.typ
│   └── template.typ
├── notebooks
│   ├── Julia
│   └── R
│       ├── RDT-equivalence_files/
│       ├── RDT-equivalence.html
│       └── RDT-equivalence.rmd
├── out
│   ├── 2025-03-13_11:00:00_optimization-df.jld2
│   ├── ensemble
│   │   ├── optimal-threshold-results
│   │   ├── scenario-optimization-summaries
│   │   ├── scenario-optimizations
│   │   └── seasonal-infectivity-import
│   ├── optimization-df.jld2
│   ├── singlesim
│   │   ├── single-sim_arrays.jld2
│   │   └── single-sim_setup.jld2
│   └── TEST_optimization-df.jld2
├── OutbreakDetectionUtils
│   ├── LICENSE
│   ├── Manifest.toml
│   ├── Project.toml
│   ├── README.md
│   ├── src
│   │   ├── cleaning-functions.jl
│   │   ├── collect-thresholds-vec_functions.jl
│   │   ├── detection-thresholds.jl
│   │   ├── diag-testing-functions.jl
│   │   ├── DrWatson-helpers.jl
│   │   ├── dynamics-constants.jl
│   │   ├── ensemble-functions.jl
│   │   ├── noise-functions.jl
│   │   ├── optimal-threshold-functions.jl
│   │   ├── OutbreakDetectionUtils.jl
│   │   ├── scenario-optimizations.jl
│   │   ├── SEIR-model.jl
│   │   ├── structs.jl
│   │   ├── test-constants.jl
│   │   ├── threshold-optimization-functions.jl
│   │   └── transmission-functions.jl
│   └── test
│       ├── cleaning-functions.jl
│       ├── collect-thresholds-vec_functions.jl
│       ├── detection-thresholds.jl
│       ├── diag-testing-functions.jl
│       ├── ensemble-functions.jl
│       ├── Manifest.toml
│       ├── noise-functions.jl
│       ├── optimal-threshold-functions.jl
│       ├── Project.toml
│       ├── runtests.jl
│       └── SEIR-model.jl
├── plots
│   ├── dynamical-noise-schematic.svg
│   ├── ensemble
│   │   ├── optimal-thresholds
│   │   ├── scenario-optimizations
│   │   ├── single-scenario
│   │   └── testing-comparison
│   ├── ensemble.bak
│   │   ├── 2024-01-19_testing-comparison.zip
│   │   ├── constant-thresholds
│   │   ├── optimal-thresholds
│   │   ├── single-scenario
│   │   └── testing-comparison
│   ├── ensemble.zip
│   ├── optimal-thresholds_accuracy-plot_dynamical.svg
│   ├── optimal-thresholds_accuracy-plot_poisson.svg
│   ├── poisson-noise-schematic.svg
│   ├── schematic-plot_incidence-threshold-6.svg
│   ├── schematic-plot_incidence-threshold.svg
│   ├── schematic-plot_incidence.svg
│   ├── schematic-plot_no-threshold.svg
│   ├── schematic-plot_threshold-3.svg
│   ├── schematic-plot_threshold-5.svg
│   ├── schematic-plot_threshold-6.svg
│   ├── schematic-plot_threshold-8.svg
│   ├── schematic-plot.svg
│   ├── schematic-simulation_no-shade.png
│   ├── schematic-simulation_with-shade.png
│   ├── schematic-simulation.png
│   └── singlesim
│       ├── single-sim_beta.png
│       ├── single-sim_SI-state-space.png
│       └── single-sim_timeseries.png
├── profile.pb.gz
├── Project.toml
├── README.md
├── renv/
├── renv.lock
├── scripts
│   ├── calculate-dynamical-noise-vaccintion-rates.jl
│   ├── debugging.jl
│   ├── ensemble-diag-testing_constant-thresholds.jl
│   ├── ensemble-diag-testing_optimal-thresholds_single-timeseries.jl
│   ├── ensemble-diag-testing_optimal-thresholds.jl
│   ├── ensemble-diag-testing_scenarios_plots.jl
│   ├── ensemble-sim_noise-visualizations.jl
│   ├── ensemble-sim_optimal-accuracy-isocline.jl
│   ├── ensemble-sim_optimal-accuracy-lineplot.jl
│   ├── ensemble-sim_single-scenario.jl
│   ├── ensemble-sim.jl
│   ├── line-plots.jl
│   ├── outbreak-detection-schematic.jl
│   ├── schematic-plots.jl
│   ├── scratch.jl
│   ├── scratch.R
│   ├── single-sim_plots.jl
│   ├── single-sim_setup.jl
│   ├── single-sim.jl
│   └── threshold-optimization.jl
├── src
│   ├── ensemble-inspection_plots.jl
│   ├── ensemble-parameters.jl
│   ├── isocline_plots.jl
│   ├── line_plots.jl
│   ├── makie-plotting-setup.jl
│   ├── noise_plots.jl
│   ├── optimal-thresholds_plots.jl
│   ├── outbreak-threshold-chars_plots.jl
│   ├── OutbreakDetection.jl
│   ├── plotting-functions.jl
│   ├── plotting-helpers.jl
│   ├── quantile_plots.jl
│   ├── R
│   │   └── app.R
│   ├── single-scenario_plots.jl
│   ├── single-sim_plots.jl
│   └── threshold_comparison_plots.jl
├── test
│   └── runtests.jl
└── workflows
    └── CI.yml
```

- `_research`
- `data/` contains input and output data files
    - `CFR_2022.csv` contains CFR rates for representative countries in 2022
    - `input-populations.csv` contains population sizes for representative
    - `optimal-threshold-results` contains output excel tables of the optimal threshold results, separated into subdirectories by `R_0` of the simulation and the noise type
    - `seasonal-infectivity-import` contains the output data files of the outbreak detection characteristics for the ensemble simulations, separated into subdirectories by model specification. Files are saved in the Julia's HDF5-compliant `.jld2` format
    - `singlesim` contains data file for a single simulation (setup files and the output arrays)
- `notebooks` contains short notebooks to perform temporary analyses using Quarto and Rmarkdown documents
- `plots` contains all output plots
    - `ensemble` contains all plots related to the ensemble simulation
        - `optimal-thresholds` contains plots related to the optimal alert thresholds for each simulation type, separated by simulation `R_0` and noise type
            - `clinic-tested` contains the optimal threshold plots where each plot refers to a different level of the % of clinic visits that are tested, and the rows refer to the test type
            - `tests` contains the optimal threshold plots where each plot refers to a different test type and the rows refer to a different % of clinic visits that are tested
        - `single-scenario` contains the plots for a single scenario of the ensemble simulations, with subdirectories for noise type where appropriate (i.e., for alert-related metrics)
        - `testing-comparison` contains plots for alert metrics compared across test type and testing rate, separated into subdirectories by noise type. These figures are computationally expensive to compute so only produced for `R_0` = 16
    - `singlesim` contains plots for the single simulation
- `renv` contains the R package versions used for the prototype app examining the trade-off between test sensitivity and specificity and the detections ability when a number of true positives are tested
- `scripts` contains the Julia scripts used to examine single and ensemble simulations, using plotting and other functions defined in `src/*.jl` files
- `src` contains all Julia source files and functions used in the analysis pipeline and exploration scripts. These files are separated by purpose e.g., `cleaning-functions.jl` contains functions for cleaning the simulation arrays into dataframes for simpler plotting and manipulation, and `ensemble-functions.jl` contains all functions related to running the ensemble simulations. There is also an `R/` directory that runs the prototype R shiny app
- `test` contains all test scripts
- `workflows` contains the CI workflow using GitHub Actions. Currently it only contains a file that can run tests on push to the `main` branch, but it is not active
