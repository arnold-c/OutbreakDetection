# OutbreakDetection Agent Guidelines

## Build/Test Commands
- Run all tests: `just tests` or `julia --project=./OutbreakDetectionUtils -e "using Pkg; Pkg.test()"`
- Run single test file: `julia --startup-file=no --project=./OutbreakDetectionUtils -e 'using Test; @testset "filename" include("OutbreakDetectionUtils/test/filename.jl")'`
- Build manuscript: `just manuscript`
- Format code: `julia -e 'using JuliaFormatter; format(".")'` (uses `.JuliaFormatter.toml` config)

## Code Style
- **Formatting**: Use Runic (https://github.com/fredrikekre/Runic.jl) style (80 char margin). Format all files before committing.
- **Imports**: Use qualified imports (e.g., `using Distributions: Distributions`) to avoid namespace pollution.
- **Types**: Use parametric types for structs. Prefer concrete types in function signatures for performance.
- **Naming**: snake_case for functions/variables, PascalCase for types/modules. Use descriptive names (e.g., `seir_mod!` for in-place).
- **Functions**: In-place functions end with `!`. Use keyword arguments with defaults for clarity.
- **Error Handling**: Use `Try.jl` for operations with potentially invalid returns. Ensure type stability.
- **Testing**: Include Aqua.jl and JET.jl tests for code quality and type stability. All tests must pass before committing.
- **Documentation**: Document public functions with docstrings. Use triple-quoted strings with parameter descriptions.
- **Performance**: Use `@inbounds` for performance-critical loops. Prefer StaticArrays for small fixed-size arrays.
- **Project Structure**: Main package is `OutbreakDetectionUtils/`. Scripts in `scripts/`, plotting in `src/`, tests in `test/`.

## DrWatson Integration
- Scripts should start with `using DrWatson; @quickactivate "OutbreakDetection"` for path management.
- Use DrWatson helpers for data/output paths (see `OutbreakDetectionUtils/src/DrWatson-helpers.jl`).
