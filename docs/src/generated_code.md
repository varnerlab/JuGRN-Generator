# Generated Code

JuGRN generates a complete, self-contained Julia project in the output directory. Here is a description of each generated file.

## Directory Structure

```
output_dir/
├── Include.jl          # Includes all source files
├── Driver.jl           # Example script to solve the ODE system
├── Project.toml        # Julia project file with dependencies
├── README.md           # Auto-generated README
└── src/
    ├── Balances.jl     # ODE balance equations
    ├── Control.jl      # Transcription control transfer functions
    ├── Data.jl         # Data dictionary with all model parameters
    ├── Error.jl        # Error handling utilities
    ├── Inputs.jl       # External input function (user-editable)
    ├── Kinetics.jl     # Transcription/translation rate calculations
    ├── PTM.jl          # Post-translational modification rates (if applicable)
    ├── SolveBalances.jl # ODE solver wrapper
    ├── Types.jl        # Type definitions for solver
    ├── Utility.jl      # Utility functions (Jacobian, etc.)
    ├── database/
    │   └── Database.db # Biophysical constants database
    └── network/
        └── Network.dat # Stoichiometric matrix
```

## Key Files

### Data.jl

Contains the `build_data_dictionary` function that returns a `Dict{String,Any}` with:

- **Stoichiometric matrix** loaded from `Network.dat`
- **Species arrays**: gene coding lengths, mRNA lengths, protein lengths, gene abundances
- **Initial conditions** for all species
- **Binding parameters**: Hill coefficient ``n`` and dissociation constant ``K`` for each regulatory connection
- **Control parameters**: weights ``W`` for the transcription control transfer function
- **Degradation and time constant modifiers**
- **PTM parameters** (if post-translational modifications are present)

### Control.jl

The `calculate_transcription_control_array` function computes a control value ``u_j \in [0, 1]`` for each gene ``j`` using a transfer function:

```math
u_j = \frac{W_{j,RNAP} + \sum_k f_k}{1 + W_{j,RNAP} + \sum_k f_k}
```

where each activator contributes:

```math
f_k = W_{j,k} \cdot \frac{p_k^{n_{j,k}}}{K_{j,k}^{n_{j,k}} + p_k^{n_{j,k}}}
```

and each repressor contributes:

```math
f_k = W_{j,k} \cdot \left(1 - \frac{p_k^{n_{j,k}}}{K_{j,k}^{n_{j,k}} + p_k^{n_{j,k}}}\right)
```

### Kinetics.jl

Contains rate functions for transcription, translation, and degradation. Transcription rates use saturation kinetics:

```math
r_{TX,i} = k_{cat} \cdot [RNAP] \cdot \frac{G_i}{K_{SAT} + G_i}
```

where the effective ``k_{cat}`` is adjusted by gene length relative to the average transcript length.

### Network.dat

The stoichiometric matrix ``S`` where rows are species and columns are reactions. For a model with ``n_g`` genes, ``n_m`` mRNAs, ``n_p`` proteins, and ``n_{PTM}`` post-translational reactions:

- Columns 1 to ``n_m``: transcription reactions
- Columns ``n_m + 1`` to ``n_m + n_p``: translation reactions
- Columns ``n_m + n_p + 1`` to end: PTM reactions (if present)

### PTM.jl

Generated only when post-translational modification reactions are present. Contains the `calculate_ptm_rates` function with:

- **Michaelis-Menten** kinetics for phosphorylation/dephosphorylation
- **Mass-action** kinetics for binding/unbinding
