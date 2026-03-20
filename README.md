## Gene Regulatory Network Generator in Julia (JuGRN)

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://varnerlab.github.io/JuGRN-Generator/stable/)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://varnerlab.github.io/JuGRN-Generator/dev/)

### Introduction
JuGRN is a code generation system that transforms simple descriptions of the connectivity of gene regulatory networks into model code written in the [Julia](https://julialang.org) programming language. JuGRN has been used in the publications:

1. [Tasseff R, Jensen H, Congleton J, Dai W, Rogers K, Sagar A, Yen A and J. Varner (2017) An Effective Model of the Retinoic Acid Induced Differentiation Program, Sci Reports, 7:14327 doi:10.1038/s41598-017-14523-5](https://www.nature.com/articles/s41598-017-14523-5)
2. [Gould R, Bassen DM, Chakrabarti A, Varner JD and Butcher J (2016) Population Heterogeneity in the Epithelial to Mesenchymal Transition Is Controlled by NFAT and Phosphorylated Sp1. PLoS Comput Biol 12(12): e1005251. doi:10.1371/journal.pcbi.1005251](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005251)

### Installation

`JuGRN.jl` requires Julia 1.6 or higher. Install from the Julia package manager:

```julia
pkg> add https://github.com/varnerlab/JuGRN-Generator.git
```

### Quick Start

Define a network file (`MyNetwork.net`):

```
sigma70 activates gene_deGFP
sigma70 activates gene_sigma28
protein_sigma28 activates gene_deGFP
protein_deGFP inhibits gene_sigma28
```

Generate model code:

```julia
using JuGRN
make_julia_model("MyNetwork.net", "my_model"; host_type=:cell_free)
```

### Features

- **Two input formats**: sentence-based `.net` for quick prototyping, structured `.json` with real sequence data
- **Compound regulatory syntax**: `(x, y) activates z` and `x inhibits (a, b)`
- **Post-translational modifications**: phosphorylation, dephosphorylation, binding, unbinding
- **JSON template generation**: convert `.net` to editable `.json` templates
- **Model persistence**: save/load via JLD2 for portable model objects
- **Multiple host types**: `:bacteria`, `:mammalian`, `:cell_free`

### Documentation

Full documentation is available at **[varnerlab.github.io/JuGRN-Generator](https://varnerlab.github.io/JuGRN-Generator/stable/)**.

### Format for the GRN model input file

JuGRN-Generator takes flat files of the form:

```
// three gene memory network
gene_1 induces (gene_2, gene_3)
gene_2 activates gene_3
gene_3 activates gene_2

// post-translational modifications
kinase_A phosphorylates protein_Y at S621 gives protein_Y_p
(protein_A, protein_B) bind complex_AB
complex_AB activates gene_target
```

See the [documentation](https://varnerlab.github.io/JuGRN-Generator/stable/) for full syntax details.
