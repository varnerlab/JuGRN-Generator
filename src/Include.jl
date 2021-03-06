# julia modules -
using Dates
using JSON
using ProgressMeter
using Logging
using DataFrames
using FASTX
using BioSequences

# constants -
const path_to_package = dirname(pathof(@__MODULE__))

# my code -
include("Types.jl")
include("Extensions.jl")
include("Factory.jl")
include("Macros.jl")
include("Parser.jl")
include("Problem.jl")
include("Common.jl")
include("Main.jl")
include("./strategy/JuliaStrategy.jl")
