# [Network Files (.net)](@id net_format)

The `.net` format is a sentence-based specification for gene regulatory networks. Each line describes a regulatory relationship using natural language-like syntax.

## Basic Syntax

```
actor action target
```

where:
- **actor** is the regulatory species (e.g., a transcription factor protein)
- **action** is the type of regulation (e.g., `activates`, `inhibits`)
- **target** is the regulated gene

### Example

```
// Comments start with // or #
sigma70 activates gene_deGFP
GntR inhibits gene_Venus
sigma70 activates gene_Venus
```

## Action Verbs

### Transcriptional Regulation

| Verb | Type | Synonyms |
|------|------|----------|
| `activates` | Activation | `activate`, `activated`, `induces`, `induce`, `induced` |
| `inhibits` | Repression | `inhibit`, `inhibited`, `represses`, `repress`, `repressed` |

### Post-Translational Modifications

| Verb | Type | Synonyms |
|------|------|----------|
| `phosphorylates` | Phosphorylation | `phosphorylate`, `phosphorylated` |
| `dephosphorylates` | Dephosphorylation | `dephosphorylate`, `dephosphorylated` |
| `bind` | Binding | `binds`, `bound` |
| `unbinds` | Unbinding | `unbind`, `dissociate`, `dissociates` |

## Species Generation

From transcriptional regulation sentences, JuGRN automatically generates species triples:

| Symbol in `.net` | Generated Species |
|------------------|-------------------|
| `gene_X` | `gene_X` (gene), `mRNA_gene_X` (mRNA), `protein_gene_X` (protein) |

For example, the line `sigma70 activates gene_deGFP` generates **9 species** (3 triples for `sigma70`, `gene_deGFP`, and any other genes mentioned).

## Compound Syntax

### Many-to-One Regulation

Multiple actors can regulate a single target using parenthesized, comma-separated lists:

```
(GntR, sigma70) activates gene_Venus
```

This creates a compound regulatory connection where both `GntR` and `sigma70` jointly activate `gene_Venus`. The generated control function uses a combined Hill function with parameter names like `K_gene_Venus_GntR_sigma70`.

### One-to-Many Regulation

A single actor can regulate multiple targets:

```
LacI inhibits (gene_Venus, gene_GFP)
```

This expands into separate connections, each with its own parameters (`K_gene_Venus_LacI`, `K_gene_GFP_LacI`).

## Post-Translational Modification Syntax

See [Post-Translational Modifications](ptm.md) for the extended syntax supporting phosphorylation and binding events.

## Comments and Whitespace

- Lines starting with `//` or `#` are treated as comments
- Blank lines are ignored
- Whitespace is flexible (any amount of space between tokens)
