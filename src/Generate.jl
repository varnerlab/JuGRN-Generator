"""
    generate_json_template(path_to_net_file::String, path_to_json_output::String; host_type::Symbol=:bacteria)

Generate a JSON model template from a `.net` file. The template contains all species (gene/mRNA/protein triples),
transcription and translation models inferred from the network topology, and default system parameters.
Sequence paths are left as empty strings for the user to fill in.

This provides a bridge from the quick `.net` format to the richer `.json` format, giving users a starting
point they can customize with real sequence data, compartment information, and system parameters.
"""
function generate_json_template(path_to_net_file::String, path_to_json_output::String;
    host_type::Symbol=:bacteria)

    # Parse the .net file -
    statement_vector = _parse_net_file(path_to_net_file)

    # Build species and connections using the existing pipeline -
    problem_object = generate_problem_object(statement_vector)
    list_of_species = problem_object.list_of_species
    list_of_connections = problem_object.list_of_connections

    # Map host_type to string -
    host_type_map = Dict(:bacteria => "bacteria", :mammalian => "mammalian", :cell_free => "CF_PURE")
    host_string = get(host_type_map, host_type, "bacteria")

    # Build system parameters with defaults -
    system_parameters = Dict{String,Any}(
        "RNAP" => 0.07,
        "RIBOSOME" => 0.07
    )

    # Build species list -
    json_species = []
    for species_object in list_of_species

        species_type_str = species_object.species_type == :mrna ? "mRNA" : string(species_object.species_type)
        species_entry = Dict{String,Any}(
            "symbol" => species_object.species_symbol,
            "type" => species_type_str,
            "compartment" => "system",
            "sequence" => ""
        )
        push!(json_species, species_entry)
    end

    # Build transcription models from connections -
    # Group connections by target gene -
    genes = extract_species_of_type(list_of_species, :gene)
    json_transcription_models = []

    for gene_object in genes

        gene_symbol = gene_object.species_symbol

        # find the corresponding mRNA symbol -
        # The .net pipeline creates mRNA as "mRNA_<gene_symbol>" directly
        mRNA_symbol = "mRNA_$(gene_symbol)"
        # verify it exists in the species list -
        mRNA_exists = any(sp -> sp.species_type == :mrna && sp.species_symbol == mRNA_symbol, list_of_species)
        if !mRNA_exists
            mRNA_symbol = ""
            @warn "Could not find matching mRNA for gene $(gene_symbol)"
        end

        # find activating and inhibiting connections -
        activating = is_species_a_target_in_connection_list(list_of_connections, gene_object, :activate)
        inhibiting = is_species_a_target_in_connection_list(list_of_connections, gene_object, :inhibit)

        list_of_activators = []
        for conn in activating
            push!(list_of_activators, Dict{String,Any}(
                "symbol" => conn.connection_symbol,
                "type" => "positive"
            ))
        end

        list_of_repressors = []
        for conn in inhibiting
            push!(list_of_repressors, Dict{String,Any}(
                "symbol" => conn.connection_symbol,
                "type" => "negative"
            ))
        end

        tx_model = Dict{String,Any}(
            "title" => "$(gene_symbol)_promoter",
            "RNAP_symbol" => "RNAP",
            "input" => gene_symbol,
            "output" => mRNA_symbol,
            "compartment" => "system",
            "list_of_activators" => list_of_activators,
            "list_of_repressors" => list_of_repressors
        )
        push!(json_transcription_models, tx_model)
    end

    # Build translation models -
    mRNAs = extract_species_of_type(list_of_species, :mrna)
    proteins = extract_species_of_type(list_of_species, :protein)
    json_translation_models = []

    for (i, mRNA_object) in enumerate(mRNAs)

        protein_symbol = i <= length(proteins) ? proteins[i].species_symbol : "protein_unknown"
        tl_model = Dict{String,Any}(
            "ribosome_symbol" => "RIBOSOME",
            "title" => "translation_$(mRNA_object.species_symbol)",
            "input" => mRNA_object.species_symbol,
            "output" => protein_symbol,
            "compartment" => "system"
        )
        push!(json_translation_models, tl_model)
    end

    # Assemble the full JSON structure -
    json_model = Dict{String,Any}(
        "system" => Dict{String,Any}(
            "parameters" => system_parameters,
            "host" => host_string
        ),
        "list_of_species" => json_species,
        "list_of_transcription_models" => json_transcription_models,
        "list_of_translation_models" => json_translation_models
    )

    # Write -
    open(path_to_json_output, "w") do f
        JSON.print(f, json_model, 4)
    end

    @info "JSON template written to $(path_to_json_output). Edit to add sequence paths and tune parameters."
    return path_to_json_output
end

function _gene_name_from_symbol(symbol::String)::String
    # Extract the base gene name from species symbols like:
    # gene_Venus -> Venus, mRNA_Venus -> Venus, protein_Venus -> Venus
    for prefix in ["gene_", "mRNA_", "protein_", "P_"]
        if startswith(symbol, prefix)
            return symbol[length(prefix)+1:end]
        end
    end
    return symbol
end
