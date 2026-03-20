# JuGRN.jl

*A code generation system for gene regulatory network models in Julia.*

## Overview

JuGRN transforms simple descriptions of gene regulatory network (GRN) connectivity into executable model code written in [Julia](https://julialang.org). Given a network specification file, JuGRN generates:

- Stoichiometric matrices for transcription, translation, and post-translational reactions
- Kinetic rate expressions (transcription, translation, degradation, and post-translational modification)
- Transcription control transfer functions (Hill-function style)
- Data dictionaries with default biophysical parameters
- ODE balance equations and solver scripts

JuGRN supports two input formats: a simple sentence-based `.net` format for quick prototyping, and a structured `.json` format for richer models with real sequence data.

## Publications

JuGRN has been used in the following publications:

1. [Tasseff R, Jensen H, Congleton J, Dai W, Rogers K, Sagar A, Yen A and J. Varner (2017) An Effective Model of the Retinoic Acid Induced Differentiation Program, *Sci Reports*, 7:14327](https://www.nature.com/articles/s41598-017-14523-5)
2. [Gould R, Bassen DM, Chakrabarti A, Varner JD and Butcher J (2016) Population Heterogeneity in the Epithelial to Mesenchymal Transition Is Controlled by NFAT and Phosphorylated Sp1. *PLoS Comput Biol* 12(12): e1005251](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005251)

## Features

- **Two input formats**: `.net` (sentence-based) and `.json` (structured with sequences)
- **Automatic species generation**: genes, mRNA, and proteins are created from network topology
- **Post-translational modifications**: phosphorylation, dephosphorylation, binding, and unbinding
- **Compound syntax**: many-to-one and one-to-many regulatory relationships
- **JSON template generation**: convert `.net` files to editable `.json` templates
- **Model persistence**: save/load models via JLD2 for portable model objects
- **Multiple host types**: bacteria, mammalian, and cell-free expression systems
