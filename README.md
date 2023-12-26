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

To run the project, you can use the Makefile.
First, you should make sure you have Make installed on your computer.
If you have macOS, it should be pre-installed, but if you'd like to use the most up-to-date version, you can  install it using homebrew.

```bash
brew install make
```

> If you use the homebrew version of Make then you will need to use `gmake` instead of `make` to run the Makefile.

If you are on Linux, you can install Make with your package manager.
If you are on Windows, you can use the [Chocolatey](https://chocolatey.org/) package manager, download it [directly](https://gnuwin32.sourceforge.net/packages/make.html), or use the ubuntu package manager on [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) (recommended as it will give you access to a linux terminal).

Once you have a version of Make installed you can use the following terminal command to run the Makefile:

```bash
make
```

If you would like to run a specific target, you can specify it after the `make` command e.g. `make ensemble-sim`

Make will track when each file was last run and save a temporary file to the `tmp/` directory.
If you want to force a re-run then you can delete the associated temporary file, either manually, or by running the associated Make clean command e.g. `make clean-ensemble-sims` to delete all ensemble simulation files.
