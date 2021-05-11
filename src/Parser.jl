"""
    parse_grn_file(path_to_model_file::String)->Union{Array{VGRNSentence,1},Exception}

Tis comment will be awesome. The best every comment ever. 
"""
function parse_grn_file(path_to_model_file::String)::Union{Dict{String,Any},Exception}

    try

        # load the JSON GRN file -
        grn_model_dictionary = JSON.parsefile(path_to_model_file)

        # return -
        return grn_model_dictionary
    catch error
        rethrow(error)
    end
end
