# Model Generation

The primary function in JuGRN is [`make_julia_model`](@ref), which generates a complete set of Julia source files from a network specification.

## From a File

```julia
using JuGRN

# From a .net file
make_julia_model("MyNetwork.net", "output_dir"; host_type=:bacteria)

# From a .json file (host_type inferred from JSON)
make_julia_model("MyNetwork.json", "output_dir")
```

## From a ProblemObject

If you have a saved model (see [Model Persistence (JLD2)](persistence.md)), you can regenerate code from the loaded object:

```julia
(model, metadata) = load_model("my_model.jld2")
make_julia_model(model, "regenerated_output"; host_type=:cell_free)
```

## Host Types

The `host_type` keyword controls biophysical parameters and gene abundance representations:

| Host Type | Description | Gene Abundance |
|-----------|-------------|----------------|
| `:bacteria` | Bacterial cells (default) | Copy number per cell (default: 2.0) |
| `:mammalian` | Mammalian cells | Copy number per cell (default: 2.0) |
| `:cell_free` | Cell-free expression (e.g., PURE system) | Concentration in nM (default: 5.0) |

## Control Function Generation

By default, JuGRN generates transcription control functions from the network topology using Hill-function transfer functions. To generate blank control stubs instead (for manual implementation):

```julia
make_julia_model("MyNetwork.net", "output_dir";
    host_type=:bacteria,
    control_function_generation=false)
```

## Pipeline

Internally, `make_julia_model` follows this pipeline:

1. **Parse** the input file ([`parse_grn_file`](@ref))
2. **Build** a `ProblemObject` containing species, connections, and configuration
3. **Generate** program components: Data.jl, Kinetics.jl, Control.jl, Inputs.jl, Network.dat
4. **Write** generated files to the output directory
5. **Transfer** distribution files (Balances.jl, SolveBalances.jl, utility functions)
6. **Generate** PTM.jl if post-translational modification reactions are present
