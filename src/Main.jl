"""
    make_julia_model(path_to_model_file::String,path_to_output_dir::String; host_type::Symbol = :bacteria, model_type::Symbol = :continuous)

This model generation will be awesome. super awesome. It will be great. The best ever
"""
function make_julia_model(path_to_model_file::String, path_to_output_dir::String; 
    host_type::Symbol=:bacteria, model_type::Symbol=:continuous, control_function_generation::Bool=true)

    try

        # Load the statement_vector -
        statement_vector::Array{VGRNSentence,1} = parse_grn_file(path_to_model_file)

        # Generate the problem object -
        problem_object = generate_problem_object(statement_vector)

        # initialize the program set -
        component_set = Set{ProgramComponent}()

        # Write the Inputs -
        program_component_inputs = build_inputs_buffer(problem_object)
        push!(component_set, program_component_inputs)

        # Write the data_dictionary -
        program_component_data_dictionary = build_data_dictionary_buffer(problem_object, host_type)
        push!(component_set, program_component_data_dictionary)

        # v1.0: We transfer a pre-baked Kinetics.jl -
        # Write the Kinetics ->
        # program_component_kinetics = build_kinetics_buffer(problem_object)
        # push!(component_set,program_component_kinetics)

        # Write the Control functions -
        program_component_control = nothing
        if (control_function_generation == true)
            program_component_control = build_control_buffer(problem_object)
        else
            program_component_control = build_blank_control_buffer(problem_object)
        end
        push!(component_set, program_component_control)

        # v1.0: We are generating this matrix in code, as part of the data_dictionary build process -
        # Write the dilution_matrix --
        # program_component_dilution_matrix = generate_dilution_matrix_buffer(problem_object)
        # push!(component_set,program_component_dilution_matrix)

        # Write the degradation_matrix --
        # program_component_degradation_matrix = generate_degradation_matrix_buffer(problem_object)
        # push!(component_set,program_component_degradation_matrix)

        # Write the stoichiometric_matrix --
        program_component_stoichiometric_matrix = generate_stoichiomteric_matrix_buffer(problem_object)
        push!(component_set, program_component_stoichiometric_matrix)

        # Dump the component_set to disk -
        write_program_components_to_disk(path_to_output_dir, component_set)

        # Transfer distrubtion jl files to the output -> these files are shared between model types
        transfer_distribution_files("$(path_to_package)/distribution", path_to_output_dir, ".jl")

        # Transfer model type specific files -
        if model_type == :continuous

            # transfer continuous files -
            transfer_distribution_files("$(path_to_package)/distribution/continuous", path_to_output_dir, ".jl")
            transfer_distribution_files("$(path_to_package)/distribution/continuous", path_to_output_dir, ".toml")

        elseif model_type == :discrete

            # transfer discrete files -
            transfer_distribution_files("$(path_to_package)/distribution/discrete", path_to_output_dir, ".jl")
        else
            throw(error("unsupported model_type. Expected {continuous,discrete} got $(model_type)"))
        end


        # Transfer distibution JSON files to output -
        transfer_distribution_files("$(path_to_package)/distribution", path_to_output_dir, ".json")

        # transfer the README files -
        path_to_readme_file = splitdir(path_to_model_file)[1]
        transfer_distribution_files("$(path_to_package)/distribution", path_to_readme_file, ".md")

    catch error
        throw(error)
    end
end

