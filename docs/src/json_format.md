# [JSON Files (.json)](@id json_format)

The `.json` format provides a richer model specification with explicit species definitions, sequence data, and system parameters. Use this format when you need precise control over gene/protein lengths, compartmentalization, or system-level parameters.

## Structure

A JSON model file has four top-level sections:

```json
{
    "system": { ... },
    "list_of_species": [ ... ],
    "list_of_transcription_models": [ ... ],
    "list_of_translation_models": [ ... ]
}
```

## System Section

Defines the host organism and global parameters:

```json
"system": {
    "host": "CF_PURE",
    "parameters": {
        "volume": 15e-6,
        "temperature": 310.15
    }
}
```

Supported host values: `"CF_PURE"` (cell-free), `"bacteria"`, `"mammalian"`.

## Species List

Each species has a symbol, type, compartment, and optional path to a FASTA sequence file:

```json
"list_of_species": [
    {
        "symbol": "gene_deGFP",
        "type": "gene",
        "compartment": "cytoplasm",
        "sequence": "./sequences/gene_deGFP.fasta"
    },
    {
        "symbol": "mRNA_gene_deGFP",
        "type": "mRNA",
        "compartment": "cytoplasm",
        "sequence": ""
    },
    {
        "symbol": "protein_deGFP",
        "type": "protein",
        "compartment": "cytoplasm",
        "sequence": "./sequences/protein_deGFP.fasta"
    }
]
```

When a FASTA sequence path is provided, JuGRN computes the actual nucleotide/amino acid length and uses it in the generated model (instead of the default 1000 nt for genes).

## Transcription Models

Define which genes are transcribed and their regulatory inputs:

```json
"list_of_transcription_models": [
    {
        "input": "gene_deGFP",
        "output": "mRNA_gene_deGFP",
        "list_of_activators": [
            {"symbol": "protein_sigma70", "type": "protein"}
        ],
        "list_of_repressors": [
            {"symbol": "protein_GntR", "type": "protein"}
        ]
    }
]
```

## Translation Models

Define which mRNAs are translated:

```json
"list_of_translation_models": [
    {
        "input": "mRNA_gene_deGFP",
        "output": "protein_deGFP"
    }
]
```

## Generating a JSON Template from `.net`

You can generate a JSON template from an existing `.net` file using [`generate_json_template`](@ref):

```julia
generate_json_template("MyNetwork.net", "MyNetwork.json"; host_type=:cell_free)
```

This creates a fully populated JSON file with default parameters that you can customize with real sequence data and system parameters. See [JSON Template Generation](@ref) for details.
