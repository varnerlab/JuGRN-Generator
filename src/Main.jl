"""
    make_julia_model(path_to_model_file::String, path_to_output_dir::String; host_type::Symbol=:bacteria, control_function_generation::Bool=true, model_type::Symbol=:standard)

Generate an executable Julia gene regulatory network (GRN) model from a network specification file.
Parses the network topology, generates stoichiometric matrices, kinetics, control logic, and data dictionary
functions, then writes the complete model to `path_to_output_dir`.

Supports `.net` (sentence-based) and `.json` (structured) input formats.

Input arguments:
- `path_to_model_file::String` - path to the GRN specification file (`.net` or `.json`).
- `path_to_output_dir::String` - path to the directory where generated model code will be written.
- `host_type::Symbol` - host organism type. Supported values: `:bacteria`, `:mammalian`, `:cell_free` (default: `:bacteria`).
- `control_function_generation::Bool` - if `true`, generates transcription control functions from the network topology; if `false`, generates blank control stubs (default: `true`).
- `model_type::Symbol` - model formulation type. `:standard` uses simple saturation kinetics; `:effective` uses the Adhikari et al. (2020) effective biophysical model with resource competition and thermodynamic control (default: `:standard`).
"""
function make_julia_model(path_to_model_file::String, path_to_output_dir::String;
    host_type::Symbol=:bacteria, control_function_generation::Bool=true, model_type::Symbol=:standard)

    # Load and parse the model file -
    parsed_result = parse_grn_file(path_to_model_file)

    # Build the problem object from the parsed result -
    problem_object = nothing
    if isa(parsed_result, Array{VGRNSentence,1})

        # .net file pipeline -
        problem_object = generate_problem_object(parsed_result)
    elseif isa(parsed_result, Dict{String,Any})

        # .json file pipeline -
        problem_object = generate_problem_object(parsed_result)

        # override host_type from JSON system.host if not explicitly set by caller -
        if haskey(problem_object.configuration_dictionary, "host_type")
            json_host = problem_object.configuration_dictionary["host_type"]
            host_type_map = Dict("CF_PURE" => :cell_free, "bacteria" => :bacteria, "mammalian" => :mammalian)
            host_type = get(host_type_map, json_host, host_type)
        end
    else
        throw(ArgumentError("Unsupported parse result type: $(typeof(parsed_result))"))
    end

    # generate from the problem object -
    _generate_model_from_problem_object(problem_object, path_to_output_dir;
        host_type=host_type, control_function_generation=control_function_generation, model_type=model_type)
end

"""
    make_julia_model(problem_object::ProblemObject, path_to_output_dir::String; host_type::Symbol=:bacteria, control_function_generation::Bool=true, model_type::Symbol=:standard)

Generate an executable Julia GRN model from a `ProblemObject`. Use this method to regenerate model code
from a previously saved (e.g., JLD2) model object.

- `model_type::Symbol` - model formulation type. `:standard` uses simple saturation kinetics; `:effective` uses the Adhikari et al. (2020) effective biophysical model (default: `:standard`).
"""
function make_julia_model(problem_object::ProblemObject, path_to_output_dir::String;
    host_type::Symbol=:bacteria, control_function_generation::Bool=true, model_type::Symbol=:standard)

    # check if host_type is stored in the problem object -
    if haskey(problem_object.configuration_dictionary, "host_type")
        json_host = problem_object.configuration_dictionary["host_type"]
        host_type_map = Dict("CF_PURE" => :cell_free, "bacteria" => :bacteria, "mammalian" => :mammalian)
        host_type = get(host_type_map, json_host, host_type)
    end

    _generate_model_from_problem_object(problem_object, path_to_output_dir;
        host_type=host_type, control_function_generation=control_function_generation, model_type=model_type)
end

function _generate_model_from_problem_object(problem_object::ProblemObject, path_to_output_dir::String;
    host_type::Symbol=:bacteria, control_function_generation::Bool=true, model_type::Symbol=:standard)

    # initialize -
    src_component_set = Set{ProgramComponent}()
    network_component_set = Set{ProgramComponent}()

    # create the output paths -
    _PATH_TO_OUTPUT_SRC_DIR = joinpath(path_to_output_dir, "src")
    _PATH_TO_ROOT_DIR = path_to_output_dir
    _PATH_TO_OUTPUT_DATABASE_DIR = joinpath(_PATH_TO_OUTPUT_SRC_DIR, "database")
    _PATH_TO_NETWORK_OUTPUT_DIR = joinpath(_PATH_TO_OUTPUT_SRC_DIR, "network")

    # Write the Inputs -
    program_component_inputs = build_inputs_buffer(problem_object)
    push!(src_component_set, program_component_inputs)

    if model_type == :effective

        # ---- Effective biophysical model (Adhikari et al., 2020) ---- #

        # Write the data_dictionary (effective version) -
        program_component_data_dictionary = build_data_dictionary_buffer_effective(problem_object, host_type)
        push!(src_component_set, program_component_data_dictionary)

        # Write the Kinetics (effective version) -
        program_component_kinetics = build_kinetics_buffer_effective(problem_object)
        push!(src_component_set, program_component_kinetics)

        # Write the Control functions (effective version with thermodynamic control) -
        program_component_control = build_control_buffer_effective(problem_object)
        push!(src_component_set, program_component_control)

    else

        # ---- Standard model ---- #

        # Write the data_dictionary -
        program_component_data_dictionary = build_data_dictionary_buffer(problem_object, host_type)
        push!(src_component_set, program_component_data_dictionary)

        # Write the Kinetics -
        program_component_kinetics = build_kinetics_buffer(problem_object)
        push!(src_component_set, program_component_kinetics)

        # Write the Control functions -
        if control_function_generation
            program_component_control = build_control_buffer(problem_object)
        else
            program_component_control = build_blank_control_buffer(problem_object)
        end
        push!(src_component_set, program_component_control)
    end

    # Write the stoichiometric_matrix -
    program_component_stoichiometric_matrix = generate_stoichiomteric_matrix_buffer(problem_object)
    push!(network_component_set, program_component_stoichiometric_matrix)

    # Dump the component_set to disk -
    write_program_components_to_disk(_PATH_TO_OUTPUT_SRC_DIR, src_component_set)
    write_program_components_to_disk(_PATH_TO_NETWORK_OUTPUT_DIR, network_component_set)

    # Transfer distribution files to the output -
    transfer_distribution_files("$(path_to_package)/distribution", _PATH_TO_OUTPUT_SRC_DIR, ".jl")

    # transfer continuous solver files -
    transfer_distribution_files("$(path_to_package)/distribution/continuous", _PATH_TO_OUTPUT_SRC_DIR, ".jl")
    transfer_distribution_files("$(path_to_package)/distribution/continuous", _PATH_TO_ROOT_DIR, ".toml")
    transfer_distribution_files("$(path_to_package)/distribution/continuous/root", _PATH_TO_ROOT_DIR, ".jl")
    transfer_distribution_files("$(path_to_package)/distribution/continuous/database", _PATH_TO_OUTPUT_DATABASE_DIR, ".db")

    # transfer the README -
    transfer_distribution_files("$(path_to_package)/distribution", _PATH_TO_ROOT_DIR, ".md")

    if model_type == :effective

        # For the effective model, overwrite the distribution Kinetics.jl and Types.jl
        # with the effective versions (distribution transfer overwrites our generated files) -
        effective_overwrite_set = Set{ProgramComponent}()
        push!(effective_overwrite_set, build_kinetics_buffer_effective(problem_object))
        push!(effective_overwrite_set, build_types_buffer_effective())
        write_program_components_to_disk(_PATH_TO_OUTPUT_SRC_DIR, effective_overwrite_set)
    end

    # Write PTM kinetics file (after distribution transfer to avoid overwrite) -
    program_component_ptm = build_ptm_kinetics_buffer(problem_object)
    if !isnothing(program_component_ptm)
        ptm_set = Set{ProgramComponent}()
        push!(ptm_set, program_component_ptm)
        write_program_components_to_disk(_PATH_TO_OUTPUT_SRC_DIR, ptm_set)
    end
end

"""
    save_model(path_to_file::String, problem_object::ProblemObject; metadata::Dict{String,Any}=Dict{String,Any}())

Save a `ProblemObject` to a JLD2 file. The saved file contains the complete model specification
(species, connections, configuration, sequence lengths) and can be loaded later to regenerate model code.

Optional `metadata` dictionary can store additional information (e.g., author, description, version).
"""
function save_model(path_to_file::String, problem_object::ProblemObject;
    metadata::Dict{String,Any}=Dict{String,Any}())

    # ensure .jld2 extension -
    if !endswith(path_to_file, ".jld2")
        path_to_file *= ".jld2"
    end

    jldsave(path_to_file;
        problem_object=problem_object,
        metadata=metadata,
        jugrn_version=string(@__MODULE__) * " v0.1.0",
        save_timestamp=string(Dates.now())
    )

    @info "Model saved to $(path_to_file)"
    return path_to_file
end

"""
    save_model(path_to_jld2::String, path_to_model_file::String; host_type::Symbol=:bacteria, metadata::Dict{String,Any}=Dict{String,Any}())

Parse a GRN specification file (`.net` or `.json`) and save the resulting `ProblemObject` to JLD2.
"""
function save_model(path_to_jld2::String, path_to_model_file::String;
    host_type::Symbol=:bacteria, metadata::Dict{String,Any}=Dict{String,Any}())

    # parse and build -
    parsed_result = parse_grn_file(path_to_model_file)
    problem_object = nothing
    if isa(parsed_result, Array{VGRNSentence,1})
        problem_object = generate_problem_object(parsed_result)
    elseif isa(parsed_result, Dict{String,Any})
        problem_object = generate_problem_object(parsed_result)
    else
        throw(ArgumentError("Unsupported parse result type: $(typeof(parsed_result))"))
    end

    # store source info in metadata -
    metadata["source_file"] = basename(path_to_model_file)
    metadata["host_type"] = string(host_type)

    return save_model(path_to_jld2, problem_object; metadata=metadata)
end

"""
    load_model(path_to_file::String) -> (problem_object::ProblemObject, metadata::Dict{String,Any})

Load a `ProblemObject` from a JLD2 file. Returns the problem object and any associated metadata.

The loaded problem object can be passed directly to `make_julia_model` to regenerate model code.
"""
function load_model(path_to_file::String)

    data = JLD2.load(path_to_file)
    problem_object = data["problem_object"]
    metadata = get(data, "metadata", Dict{String,Any}())

    # add load-time info -
    metadata["jugrn_version"] = get(data, "jugrn_version", "unknown")
    metadata["save_timestamp"] = get(data, "save_timestamp", "unknown")

    @info "Model loaded from $(path_to_file) (saved: $(metadata["save_timestamp"]))"
    return (problem_object, metadata)
end
