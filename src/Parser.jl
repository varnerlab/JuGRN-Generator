"""
    parse_grn_file(path_to_model_file::String)

Parse a gene regulatory network (GRN) specification file. Supports two formats:

- `.net` files: line-oriented sentence format (e.g., `GntR inhibits gene_Venus`).
  Returns `Array{VGRNSentence,1}` for the VGRNSentence-based pipeline.

- `.json` files: structured JSON format with species, transcription, and translation models.
  Returns `Dict{String,Any}` containing an intermediate representation with a species DataFrame.
"""
function parse_grn_file(path_to_model_file::String)

    # determine file type -
    file_extension = lowercase(splitext(path_to_model_file)[2])

    if file_extension == ".net"
        return _parse_net_file(path_to_model_file)
    elseif file_extension == ".json"
        return _parse_json_file(path_to_model_file)
    else
        throw(ArgumentError("Unsupported GRN file format: $(file_extension). Expected .net or .json"))
    end
end

function _parse_net_file(path_to_model_file::String)::Array{VGRNSentence,1}

    # known action verbs -
    action_verbs = Set(["activate", "activates", "activated", "induce", "induces", "induced",
                        "inhibit", "inhibits", "inhibited", "repress", "represses", "repressed",
                        "phosphorylate", "phosphorylates", "phosphorylated",
                        "dephosphorylate", "dephosphorylates", "dephosphorylated",
                        "bind", "binds", "bound",
                        "unbind", "unbinds", "dissociate", "dissociates"])

    statement_vector = VGRNSentence[]

    open(path_to_model_file, "r") do file
        for line in eachline(file)

            # skip comments and blank lines -
            stripped_line = strip(line)
            if isempty(stripped_line) || startswith(stripped_line, "//") || startswith(stripped_line, "#")
                continue
            end

            # find the action verb in the line -
            # this handles clauses like "(x, y, z) activates t" and "y inhibits (p, q, r)"
            tokens = split(stripped_line)
            action_index = findfirst(t -> lowercase(t) in action_verbs, tokens)
            if isnothing(action_index) || action_index == 1 || action_index == length(tokens)
                @warn "Skipping malformed line (no action verb found): $(stripped_line)"
                continue
            end

            # build the sentence object -
            sentence = VGRNSentence()
            sentence.original_sentence = stripped_line
            sentence.sentence_delimiter = ' '
            sentence.sentence_modifier_clause = ""
            sentence.sentence_product_clause = ""

            # everything before the action verb is the actor clause -
            sentence.sentence_actor_clause = String(join(tokens[1:action_index-1], " "))
            sentence.sentence_action_clause = String(tokens[action_index])

            # parse the remainder after the action verb for "at" and "gives"/"produces" keywords -
            remainder_tokens = tokens[action_index+1:end]
            _parse_post_action_clauses!(sentence, remainder_tokens)

            push!(statement_vector, sentence)
        end
    end

    return statement_vector
end

function _parse_post_action_clauses!(sentence::VGRNSentence, tokens::AbstractVector)

    # find keyword positions -
    at_index = findfirst(t -> lowercase(t) == "at", tokens)
    gives_index = findfirst(t -> lowercase(t) in ("gives", "produces", "forming"), tokens)

    if !isnothing(at_index) && !isnothing(gives_index)
        # pattern: target at site gives product
        sentence.sentence_target_clause = String(join(tokens[1:at_index-1], " "))
        sentence.sentence_modifier_clause = String(join(tokens[at_index+1:gives_index-1], " "))
        sentence.sentence_product_clause = String(join(tokens[gives_index+1:end], " "))
    elseif !isnothing(gives_index)
        # pattern: target gives product (no site)
        sentence.sentence_target_clause = String(join(tokens[1:gives_index-1], " "))
        sentence.sentence_product_clause = String(join(tokens[gives_index+1:end], " "))
    elseif !isnothing(at_index)
        # pattern: target at site (no explicit product)
        sentence.sentence_target_clause = String(join(tokens[1:at_index-1], " "))
        sentence.sentence_modifier_clause = String(join(tokens[at_index+1:end], " "))
    else
        # simple pattern: just target
        sentence.sentence_target_clause = String(join(tokens, " "))
    end
end

function _parse_json_file(path_to_model_file::String)::Dict{String,Any}

    intermediate_representation_dictionary = Dict{String,Any}()
    species_table = DataFrame(symbol=String[], type=Symbol[], compartment=String[], sequence=Union{Missing,String,FASTA.Record}[])

    # load the JSON GRN file -
    grn_model_dictionary = JSON.parsefile(path_to_model_file)

    # build the species table -
    list_of_species_dictionaries = grn_model_dictionary["list_of_species"]
    for species_dictionary in list_of_species_dictionaries

        local_data_row = Union{String,Symbol,Missing,FASTA.Record}[]
        push!(local_data_row, species_dictionary["symbol"])
        push!(local_data_row, Symbol(species_dictionary["type"]))
        push!(local_data_row, species_dictionary["compartment"])

        # load sequence if available -
        seq_path = species_dictionary["sequence"]
        if seq_path == ""
            push!(local_data_row, missing)
        else
            open(FASTA.Reader, seq_path) do reader
                for record in reader
                    push!(local_data_row, record)
                end
            end
        end

        push!(species_table, tuple(local_data_row...))
    end

    intermediate_representation_dictionary["model_species_table"] = species_table
    intermediate_representation_dictionary["raw_model"] = grn_model_dictionary

    return intermediate_representation_dictionary
end
