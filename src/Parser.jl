"""
    parse_grn_file(path_to_model_file::String)->Union{Array{VGRNSentence,1},Exception}

Tis comment will be awesome. The best every comment ever. 
"""
function parse_grn_file(path_to_model_file::String)::Union{Array{VGRNSentence,1},Exception}

    # We are going to load the sentences in the file into a vector
    # if not a valid model file, then throw an error -
    sentence_vector = VGRNSentence[]
    tmp_array::Array{String,1} = String[]

    try

        # Open the model file, and read each line into a vector -
        open(path_to_model_file, "r") do model_file
            for line in eachline(model_file)

                if (contains(line, "//") == false)
                    push!(tmp_array, chomp(line))
                end
            end
        end

        for sentence in tmp_array

            # make sentence a string -
            local_sentence = convert(String, sentence)

            if (length(local_sentence) > 1)

                # Ok, so now we have the array for sentences -
                grn_sentence = VGRNSentence()
                grn_sentence.original_sentence = local_sentence

                # split the sentence -
                split_array = split(local_sentence, " ")

                # sentence_actor_clause::String
                # sentence_action_clause::String
                # sentence_target_clause::String
                # sentence_delimiter::Char
                grn_sentence.sentence_actor_clause = split_array[1]
                grn_sentence.sentence_action_clause = split_array[2]
                grn_sentence.sentence_target_clause = split_array[3]
                grn_sentence.sentence_delimiter = ' '

                # add sentence to sentence_vector =
                push!(sentence_vector, grn_sentence)
            end
        end

        # return - (I know we don't need the return, but I *** hate *** the normal Julia convention)
        return sentence_vector
    catch err
        # showerror(stdout, err, backtrace());println()
        rethrow(error)
    end
end
