# Effective Biophysical Model

The `:effective` model type generates code based on the effective biophysical model of Adhikari et al. (2020), originally developed for cell-free transcription-translation (TX-TL) systems such as the PURE system.

> **Reference:** Adhikari A, Vilkhovoy M, Vadhin S, Lim HE and Varner JD (2020) A Genome-Scale Metabolic Reconstruction of *Escherichia coli* Producer Strains for the Cell-Free Expression System. *Front. Bioeng. Biotechnol.* 8:539081. doi: [10.3389/fbioe.2020.539081](https://doi.org/10.3389/fbioe.2020.539081)

## Usage

```julia
using JuGRN

make_julia_model("MyNetwork.net", "output_dir";
    host_type=:cell_free, model_type=:effective)
```

The `:effective` flag replaces the standard saturation-kinetics `Kinetics.jl`, `Control.jl`, and `Data.jl` with versions implementing the equations below.

## Transcription Kinetics (Eqs. 3-4)

The transcription rate for gene ``j`` is:

```math
r_{X,j} = V^{max}_{X,j} \frac{G_j}{\tau_{X,j} K_X + (1 + \tau_{X,j}) G_j + \mathcal{O}_{X,j}}
```

where:
- ``V^{max}_{X,j} = R_{X,T} \cdot v_X / \ell_{G,j}`` is the maximum transcription rate
- ``R_{X,T}`` is total RNAP concentration (default: 0.07 ``\mu``M)
- ``v_X`` is the elongation rate (default: 25 nt/s)
- ``\ell_{G,j}`` is the gene coding length (nt)
- ``\tau_{X,j}`` is a dimensionless time constant
- ``K_X`` is the saturation constant (default: 0.036 ``\mu``M)
- ``G_j`` is the gene concentration

The competition term ``\mathcal{O}_{X,j}`` captures resource sharing among genes:

```math
\mathcal{O}_{X,j} = \sum_{i \neq j} \frac{\tau_{X,j}}{\tau_{X,i}} (1 + \tau_{X,i}) G_i
```

## Translation Kinetics (Eqs. 5-6)

The translation rate for mRNA ``j`` is:

```math
r_{L,j} = V^{max}_{L,j} \frac{K_P \cdot m_j}{\tau_{L,j} K_L + (1 + \tau_{L,j}) K_P \cdot m_j + \mathcal{O}_{L,j}}
```

where:
- ``V^{max}_{L,j} = R_{L,T} \cdot v_L / \ell_{P,j}`` is the maximum translation rate
- ``R_{L,T}`` is total ribosome concentration (default: 2.3 ``\mu``M)
- ``v_L`` is the elongation rate (default: 1.5 aa/s)
- ``\ell_{P,j}`` is the protein coding length (aa)
- ``K_P`` is a polysome amplification factor (default: 10.0)
- ``K_L`` is the saturation constant (default: 450 ``\mu``M)
- ``m_j`` is the mRNA concentration

The translation competition term is:

```math
\mathcal{O}_{L,j} = \sum_{i \neq j} \frac{\tau_{L,j}}{\tau_{L,i}} (1 + \tau_{L,i}) K_P \cdot m_i
```

## Transcription Control (Eq. 9)

Transcription control uses a thermodynamic (Boltzmann-weight) formulation instead of the standard Hill-function transfer function:

```math
u_j = \frac{W_{j,RNAP} + \sum_{k \in \text{activators}} W_{j,k} f_{j,k}}{1 + W_{j,RNAP} + \sum_{k \in \text{activators}} W_{j,k} f_{j,k} + \sum_{k \in \text{repressors}} W_{j,k} f_{j,k}}
```

where the Boltzmann weights are:

```math
W = \exp\left(-\frac{\Delta G}{RT}\right)
```

and ``f_{j,k}`` is the fractional occupancy (Hill function) for regulator ``k`` at gene ``j``:

```math
f_{j,k} = \frac{p_k^{n_{j,k}}}{K_{j,k}^{n_{j,k}} + p_k^{n_{j,k}}}
```

Activators appear in both numerator and denominator; repressors appear only in the denominator.

## Translation Control (Eq. 10)

Translation capacity decays exponentially over time, modeling resource depletion in cell-free systems:

```math
w_j(t) = \exp\left(-\frac{0.693 \cdot t}{\tau_{L,1/2}}\right)
```

where ``\tau_{L,1/2}`` is the translation capacity half-life (default: 4.0 hr).

## Generated Parameters

The effective model `Data.jl` includes these additional parameters (defaults from Table 1 of Adhikari et al., 2020):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `R_XT` | 0.07 ``\mu``M | Total RNAP concentration |
| `R_LT` | 2.3 ``\mu``M | Total ribosome concentration |
| `v_X` | 25 nt/s | Transcription elongation rate |
| `v_L` | 1.5 aa/s | Translation elongation rate |
| `K_X` | 0.036 ``\mu``M | Transcription saturation constant |
| `K_L` | 450 ``\mu``M | Translation saturation constant |
| `K_P` | 10.0 | Polysome amplification factor |
| `temperature` | 310.15 K | System temperature (37 C) |
| `tau_L_half` | 4.0 hr | Translation capacity half-life |
| `tau_X_array` | 0.5 (per gene) | Transcription time constants |
| `tau_L_array` | 0.5 (per gene) | Translation time constants |
| `dG_dictionary` | -3.0 kJ/mol (per interaction) | Free energies for Boltzmann weights |

## Differences from Standard Model

| Aspect | `:standard` | `:effective` |
|--------|------------|-------------|
| Transcription rate | Simple saturation (``k_{cat} \cdot RNAP \cdot G / (K + G)``) | Competitive allocation with ``\mathcal{O}_{X,j}`` |
| Translation rate | Simple saturation | Competitive allocation with ``\mathcal{O}_{L,j}`` and polysome factor |
| Transcription control | Hill-function transfer function | Thermodynamic Boltzmann weights |
| Translation control | None (constant = 1) | Exponential capacity decay ``w(t)`` |
| Primary use case | General GRN modeling | Cell-free TX-TL systems |
