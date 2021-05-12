# abstract types
abstract type VLAbstractModelObject end

mutable struct VGRNSentence

    original_sentence::String
    sentence_actor_clause::String
    sentence_action_clause::String
    sentence_target_clause::String
    sentence_delimiter::Char

    function VGRNSentence()
        this = new()
    end
end

mutable struct ProgramComponent

    filename::String
    buffer::String

    function ProgramComponent()
        this = new()
    end

end

mutable struct SpeciesObject

    species_type::Symbol
    species_symbol::String

end

mutable struct ConnectionObject

    connection_symbol::String
    connection_actor_set::Array{SpeciesObject,1}
    connection_target_set::Array{SpeciesObject,1}
    connection_type::Symbol

    function ConnectionObject()
        this = new()
    end
end

mutable struct ProblemObject

    configuration_dictionary::Dict{String,Any}
    list_of_species::Array{SpeciesObject,1}
    list_of_connections::Array{ConnectionObject,1}

    function ProblemObject()
        this = new()
    end
end

struct VLJuliaModelObject <: VLAbstractModelObject

    # data -
    path_to_model_file::String
    path_to_output_dir::String
    defaults_file_name::String 

    # constructor -
    function VLJuliaModelObject(path_to_model_file, path_to_output_dir; 
        defaults_file_path::String="Defaults.toml")
        
        # build new model object -
        this = new(path_to_model_file, path_to_output_dir, defaults_file_path)
    end
end
