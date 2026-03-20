# Post-Translational Modifications

JuGRN supports post-translational modification (PTM) events in the `.net` format, extending the language beyond transcriptional regulation to include phosphorylation, dephosphorylation, and protein-protein binding.

## Phosphorylation

A kinase modifies a substrate protein at a specific site, producing a modified form:

```
kinase_A phosphorylates protein_Y at S621 gives protein_Y_p_S621
```

This generates:
- A Michaelis-Menten rate expression: ``v = k_{cat} \cdot E \cdot S / (K_m + S)``
- Parameters `kcat_kinase_A` and `Km_kinase_A` in the `ptm_parameter_dictionary`
- Stoichiometric coefficients: substrate consumed (``-1``), product formed (``+1``), enzyme unchanged

The `at S621` clause is optional and serves as annotation. The `gives` clause names the product species.

## Dephosphorylation

A phosphatase removes the modification:

```
phosphatase_B dephosphorylates protein_Y_p_S621 gives protein_Y
```

Uses the same Michaelis-Menten kinetics as phosphorylation.

## Binding

Two or more proteins form a complex:

```
(protein_A, protein_B) bind complex_AB
```

This generates:
- A mass-action rate expression: ``v = k_f \cdot [A] \cdot [B]``
- Parameter `kf_protein_A_protein_B` in the `ptm_parameter_dictionary`
- Stoichiometric coefficients: each reactant consumed (``-1``), complex formed (``+1``)

## Unbinding (Dissociation)

A complex dissociates into components:

```
complex_AB unbinds (protein_A, protein_B)
```

This generates:
- A first-order rate expression: ``v = k_r \cdot [complex]``
- Parameter `kr_complex_AB` in the `ptm_parameter_dictionary`
- Stoichiometric coefficients: complex consumed (``-1``), components formed (``+1``)

## Combining PTM with Transcriptional Regulation

Modified proteins and complexes can participate in transcriptional regulation:

```
// Phosphorylation creates an active form
kinase phosphorylates protein_TF at Y705 gives protein_TF_p

// The phosphorylated form activates a gene
protein_TF_p activates gene_target

// Binding creates a heterodimer
(protein_A, protein_B) bind complex_AB

// The complex activates transcription
complex_AB activates gene_X
```

## Species Generation

PTM sentences generate species differently from transcriptional regulation:

- **Transcriptional sentences** generate gene/mRNA/protein triples (e.g., `gene_X`, `mRNA_gene_X`, `protein_gene_X`)
- **PTM sentences** generate individual protein-level species for actors, targets, and products that aren't already covered by triples

## Generated Files

When PTM reactions are present, JuGRN generates an additional file:

| File | Description |
|------|-------------|
| `PTM.jl` | Rate expressions for all post-translational modification reactions |

The `Data.jl` file will also contain a `ptm_parameter_dictionary` with default kinetic parameters.

## Default Parameters

| Parameter | Default Value | Units | Description |
|-----------|--------------|-------|-------------|
| `kcat_*` | 10.0 | hr⁻¹ | Catalytic rate constant (phosphorylation/dephosphorylation) |
| `Km_*` | 0.1 | μM | Michaelis constant |
| `kf_*` | 1.0 | μM⁻¹ hr⁻¹ | Forward binding rate constant |
| `kr_*` | 0.1 | hr⁻¹ | Reverse (unbinding) rate constant |
