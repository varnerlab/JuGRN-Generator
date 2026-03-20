"""
    make_julia_model(path_to_model_file::String, path_to_output_dir::String; host_type::Symbol=:bacteria, control_function_generation::Bool=true)

Generate an executable Julia gene regulatory network (GRN) model from a network specification file.
Parses the network topology, generates stoichiometric matrices, kinetics, control logic, and data dictionary
functions, then writes the complete model to `path_to_output_dir`.

Supports `.net` (sentence-based) and `.json` (structured) input formats.

Input arguments:
- `path_to_model_file::String` - path to the GRN specification file (`.net` or `.json`).
- `path_to_output_dir::String` - path to the directory where generated model code will be written.
- `host_type::Symbol` - host organism type. Supported values: `:bacteria`, `:mammalian`, `:cell_free` (default: `:bacteria`).
- `control_function_generation::Bool` - if `true`, generates transcription control functions from the network topology; if `false`, generates blank control stubs (default: `true`).
"""
function make_julia_model(path_to_model_file::String, path_to_output_dir::String;
    host_type::Symbol=:bacteria, control_function_generation::Bool=true)

    # initialize -
    src_component_set = Set{ProgramComponent}()
    network_component_set = Set{ProgramComponent}()

    # create the output paths -
    _PATH_TO_OUTPUT_SRC_DIR = joinpath(path_to_output_dir, "src")
    _PATH_TO_ROOT_DIR = path_to_output_dir
    _PATH_TO_OUTPUT_DATABASE_DIR = joinpath(_PATH_TO_OUTPUT_SRC_DIR, "database")
    _PATH_TO_NETWORK_OUTPUT_DIR = joinpath(_PATH_TO_OUTPUT_SRC_DIR, "network")

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

    # Write the Inputs -
    program_component_inputs = build_inputs_buffer(problem_object)
    push!(src_component_set, program_component_inputs)

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
end
