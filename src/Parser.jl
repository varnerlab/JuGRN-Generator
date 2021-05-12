"""
    parse_grn_file(path_to_model_file::String)->Union{Array{VGRNSentence,1},Exception}

Tis comment will be awesome. The best every comment ever. 
"""
function parse_grn_file(path_to_model_file::String)::Dict{String,Any}

    # initialize -
    intermediate_representation_dictionary = Dict{String,Any}()
    species_table = DataFrame(symbol=String[], type=Symbol[], compartment=String[], sequence=Union{Missing,String,FASTA.Record}[])

    try

        # load the JSON GRN file -
        grn_model_dictionary = JSON.parsefile(path_to_model_file)

        # ok, so lets build an intermediate representation that is a bunch of DataFrames -
        list_of_species_dictionaries = grn_model_dictionary["list_of_species"]
        for species_dictionary in list_of_species_dictionaries
            
            # grab -
            local_data_row = Union{String,Symbol,Missing,FASTA.Record}[]
            push!(local_data_row, species_dictionary["symbol"])
            push!(local_data_row, Symbol(species_dictionary["type"]))
            push!(local_data_row, species_dictionary["compartment"])

            # Lets load the sequence -
            seq_path = species_dictionary["sequence"]
            if (seq_path == "")
                push!(local_data_row, missing)
            else 
                open(FASTA.Reader, seq_path) do reader
                    for record in reader
                        push!(local_data_row, record)
                    end
                end
            end
            
            # push into data frame -
            push!(species_table, tuple(local_data_row...))
        end

        # add to ir -
        intermediate_representation_dictionary["model_species_table"] = species_table

        # return -
        return intermediate_representation_dictionary
    catch error
        rethrow(error)
    end
end
