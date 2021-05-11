"""
    make_julia_model(path_to_model_file::String,path_to_output_dir::String; host_type::Symbol = :bacteria, model_type::Symbol = :continuous)

This model generation will be awesome. super awesome. It will be great. The best ever
"""
function make_julia_model(path_to_model_file::String, path_to_output_dir::String; 
    host_type::Symbol=:bacteria, control_function_generation::Bool=true)

    # initialize -
    src_component_set = Set{ProgramComponent}()
    root_component_set = Set{ProgramComponent}()
    network_component_set = Set{ProgramComponent}()

    try

        # ok, need to create the paths -
        _PATH_TO_OUTPUT_SRC_DIR = joinpath(path_to_output_dir, "src")
        _PATH_TO_ROOT_DIR = path_to_output_dir
        _PATH_TO_OUTPUT_DATABASE_DIR = joinpath(_PATH_TO_OUTPUT_SRC_DIR, "database")
        _PATH_TO_NETWORK_OUTPUT_DIR = joinpath(_PATH_TO_OUTPUT_SRC_DIR, "network")

        # Load the statement_vector -
        statement_vector::Array{VGRNSentence,1} = parse_grn_file(path_to_model_file)

        # Generate the problem object -
        problem_object = generate_problem_object(statement_vector)

        # Write the Inputs -
        program_component_inputs = build_inputs_buffer(problem_object)
        push!(src_component_set, program_component_inputs)

        # Write the data_dictionary -
        program_component_data_dictionary = build_data_dictionary_buffer(problem_object, host_type)
        push!(src_component_set, program_component_data_dictionary)

        # Write the Kinetics ->
        program_component_kinetics = build_kinetics_buffer(problem_object)
        push!(src_component_set, program_component_kinetics)

        # Write the Control functions -
        program_component_control = nothing
        if (control_function_generation == true)
            program_component_control = build_control_buffer(problem_object)
        else
            program_component_control = build_blank_control_buffer(problem_object)
        end
        push!(src_component_set, program_component_control)

        # v1.0: We are generating this matrix in code, as part of the data_dictionary build process -
        # Write the dilution_matrix --
        # program_component_dilution_matrix = generate_dilution_matrix_buffer(problem_object)
        # push!(component_set,program_component_dilution_matrix)

        # Write the degradation_matrix --
        # program_component_degradation_matrix = generate_degradation_matrix_buffer(problem_object)
        # push!(component_set,program_component_degradation_matrix)

        # Write the stoichiometric_matrix --
        program_component_stoichiometric_matrix = generate_stoichiomteric_matrix_buffer(problem_object)
        push!(network_component_set, program_component_stoichiometric_matrix)

        # Dump the component_set to disk -
        write_program_components_to_disk(_PATH_TO_OUTPUT_SRC_DIR, src_component_set)
        write_program_components_to_disk(_PATH_TO_NETWORK_OUTPUT_DIR, network_component_set)

        # Transfer distrubtion jl files to the output -> these files are shared between model types
        transfer_distribution_files("$(path_to_package)/distribution", _PATH_TO_OUTPUT_SRC_DIR, ".jl")

        # transfer continuous files -
        transfer_distribution_files("$(path_to_package)/distribution/continuous", _PATH_TO_OUTPUT_SRC_DIR, ".jl")
        transfer_distribution_files("$(path_to_package)/distribution/continuous", _PATH_TO_ROOT_DIR, ".toml")
        transfer_distribution_files("$(path_to_package)/distribution/continuous/root", _PATH_TO_ROOT_DIR, ".jl")
        transfer_distribution_files("$(path_to_package)/distribution/continuous/database", _PATH_TO_OUTPUT_DATABASE_DIR, ".db")

        # transfer the README files -
        transfer_distribution_files("$(path_to_package)/distribution", _PATH_TO_ROOT_DIR, ".md")

    catch error
        rethrow(error)
    end
end

