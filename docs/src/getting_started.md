# Getting Started

## Installation

JuGRN.jl can be installed from the Julia package manager. From the Julia REPL, enter package mode by pressing `]`, then:

```julia
pkg> add https://github.com/varnerlab/JuGRN-Generator.git
```

Or from a script:

```julia
using Pkg
Pkg.add(url="https://github.com/varnerlab/JuGRN-Generator.git")
```

## Requirements

- Julia 1.6 or higher
- Dependencies are installed automatically: JSON, DataFrames, FASTX, BioSequences, JLD2, ProgressMeter

## Quick Start

### 1. Define a network

Create a file called `MyNetwork.net`:

```
// A simple two-gene network
sigma70 activates gene_deGFP
sigma70 activates gene_sigma28
protein_sigma28 activates gene_deGFP
protein_deGFP inhibits gene_sigma28
```

### 2. Generate model code

```julia
using JuGRN

make_julia_model("MyNetwork.net", "my_model"; host_type=:cell_free)
```

### 3. Run the model

```julia
cd("my_model")
include("Include.jl")
include("Driver.jl")
```

The generated `Driver.jl` script solves the ODE system and produces time-course simulation data for all gene, mRNA, and protein species.

## Next Steps

- Learn about the [`.net` format](@ref net_format) for network specification
- Explore the [`.json` format](@ref json_format) for richer models with sequence data
- See [Post-Translational Modifications](ptm.md) for phosphorylation and binding events
- Use [Model Persistence (JLD2)](persistence.md) to save and share model objects
