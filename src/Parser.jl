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

    statement_vector = VGRNSentence[]

    open(path_to_model_file, "r") do file
        for line in eachline(file)

            # skip comments and blank lines -
            stripped_line = strip(line)
            if isempty(stripped_line) || startswith(stripped_line, "//") || startswith(stripped_line, "#")
                continue
            end

            # split the sentence into tokens -
            tokens = split(stripped_line)
            if length(tokens) < 3
                @warn "Skipping malformed line: $(stripped_line)"
                continue
            end

            # build the sentence object -
            sentence = VGRNSentence()
            sentence.original_sentence = stripped_line
            sentence.sentence_delimiter = ' '

            # parse: actor action target
            sentence.sentence_actor_clause = String(tokens[1])
            sentence.sentence_action_clause = String(tokens[2])
            sentence.sentence_target_clause = String(join(tokens[3:end], " "))

            push!(statement_vector, sentence)
        end
    end

    return statement_vector
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
