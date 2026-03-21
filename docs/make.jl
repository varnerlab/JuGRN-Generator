using Documenter
using JuGRN

makedocs(
    sitename = "JuGRN.jl",
    modules = [JuGRN],
    authors = "Jeffrey Varner and contributors",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://varnerlab.github.io/JuGRN-Generator",
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Input Formats" => [
            "Network Files (.net)" => "net_format.md",
            "JSON Files (.json)" => "json_format.md",
        ],
        "Features" => [
            "Model Generation" => "model_generation.md",
            "Effective Biophysical Model" => "effective_model.md",
            "Post-Translational Modifications" => "ptm.md",
            "JSON Template Generation" => "json_template.md",
            "Model Persistence (JLD2)" => "persistence.md",
        ],
        "Generated Code" => "generated_code.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/varnerlab/JuGRN-Generator.git",
    devbranch = "master",
    push_preview = true,
)
