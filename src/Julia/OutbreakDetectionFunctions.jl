"""
Note that each file is a module that defines the exported functions.
They are just listed here for convenience of sourcing one file.
"""
module OutbreakDetectionFunctions

using Reexport

include("transmission-functions.jl")
@reexport using .TransmissionFunctions

include("structs.jl")
@reexport using .ODStructs

include("SEIR-model.jl")
@reexport using .SEIRModel

include("cleaning-functions.jl")
@reexport using .CleaningFunctions

include("detection-thresholds.jl")
@reexport using .DetectionThresholds

include("diag-testing-functions.jl")
@reexport using .DiagTestingFunctions

include("ensemble-functions.jl")
@reexport using .EnsembleFunctions

include("noise-functions.jl")
@reexport using .NoiseFunctions

include("plotting-functions.jl")
@reexport using .PlottingFunctions

end
