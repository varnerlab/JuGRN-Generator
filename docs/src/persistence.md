# Model Persistence (JLD2)

JuGRN can save and load model objects using the [JLD2](https://github.com/JuliaIO/JLD2.jl) format. This allows you to:

- Share model specifications as portable binary files
- Regenerate model code from a saved object without the original input file
- Attach metadata (author, description, version) to models

## Saving a Model

### From a File

```julia
using JuGRN

save_model("my_model.jld2", "MyNetwork.net";
    host_type=:bacteria,
    metadata=Dict{String,Any}(
        "author" => "J. Varner",
        "description" => "Three-gene memory network"
    ))
```

### From a ProblemObject

```julia
# Parse and build the problem object manually
parsed = parse_grn_file("MyNetwork.net")
# ... customize the problem object ...

save_model("my_model.jld2", problem_object;
    metadata=Dict{String,Any}("version" => "2.0"))
```

## Loading a Model

```julia
(model, metadata) = load_model("my_model.jld2")

# Inspect metadata
println(metadata["author"])         # "J. Varner"
println(metadata["save_timestamp"]) # "2024-01-15T10:30:00"
println(metadata["source_file"])    # "MyNetwork.net"

# Regenerate model code
make_julia_model(model, "regenerated_output"; host_type=:bacteria)
```

## What Gets Saved

The JLD2 file contains:

| Key | Description |
|-----|-------------|
| `problem_object` | Complete `ProblemObject` with species, connections, and configuration |
| `metadata` | User-provided metadata dictionary |
| `jugrn_version` | JuGRN version string |
| `save_timestamp` | ISO 8601 timestamp of when the model was saved |

The `ProblemObject` preserves all information needed to regenerate the model, including:
- Species list (genes, mRNAs, proteins, PTM species)
- Connection list (activation, inhibition, PTM reactions)
- Configuration dictionary (parameters, sequence lengths, host type)

See the [API Reference](@ref) for the full function signatures.
