# JSON Template Generation

JuGRN can generate a structured `.json` template from a `.net` file, bridging the gap between the quick sentence-based format and the richer JSON specification.

## Usage

```julia
using JuGRN

generate_json_template("MyNetwork.net", "MyNetwork.json"; host_type=:cell_free)
```

## What Gets Generated

The template contains:

1. **System section** with host type and default parameters
2. **Species list** with all gene/mRNA/protein triples inferred from the network
3. **Transcription models** with activators and repressors from the `.net` topology
4. **Translation models** for each mRNA/protein pair

Sequence paths are left as empty strings (`""`) for you to fill in with real FASTA file paths.

## Round-Trip Workflow

A common workflow is:

```julia
# 1. Start with a quick .net file
# 2. Generate a JSON template
generate_json_template("MyNetwork.net", "MyNetwork.json"; host_type=:cell_free)

# 3. Edit the JSON: add sequence paths, adjust parameters, change compartments

# 4. Generate the model from the enriched JSON
make_julia_model("MyNetwork.json", "my_model")
```

This lets you prototype quickly with `.net` and then refine with real data in `.json`.

See the [API Reference](@ref) for the full function signature.
