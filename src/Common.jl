function contains(string, token)
    return occursin(token, string)
end

function write_program_components_to_disk(file_path::String, set_of_program_components::Set{ProgramComponent})

    # check - do we have the file path?
    if (isdir(file_path) == false)
        mkpath(file_path)
    end

    # go through each component, and dump the buffer to disk -
    for program_component in set_of_program_components

        # Get the data -
        filename = program_component.filename
        program_buffer = program_component.buffer

        # build the path -
        path_to_program_file = file_path * "/" * filename

        # Write the file -
        outfile = open(path_to_program_file, "w")
        write(outfile, program_buffer);
        close(outfile);
    end

end

function transfer_distribution_file(path_to_distribution_files::String,
                                      input_file_name_with_extension::String,
                                      path_to_output_files::String,
                                      output_file_name_with_extension::String)

    # Load the specific file -
    # create src_buffer -
    src_buffer::Array{String} = String[]

    # check - do we have the file path?
    if (isdir(path_to_output_files) == false)
        mkpath(path_to_output_files)
    end

    # path to distrubtion -
    path_to_src_file = path_to_distribution_files * "/" * input_file_name_with_extension
    open(path_to_src_file, "r") do src_file
        for line in eachline(src_file)

        # need to add a new line for some reason in Julia 0.6
            new_line_with_line_ending = line * "\n"
            push!(src_buffer, new_line_with_line_ending)
        end
    end

  # Write the file to the output -
    path_to_program_file = path_to_output_files * "/" * output_file_name_with_extension
    outfile = open(path_to_program_file, "w")
    write(outfile, src_buffer);
    close(outfile);

end

function transfer_distribution_files(path_to_distribution_files::String,
                                      path_to_output_files::String,
                                      file_extension::String)


    # check - do we have the file path?
    if (isdir(path_to_output_files) == false)
        mkpath(path_to_output_files)
    end

    # Search the directory for src files -
    # load the files -
    searchdir(path, key) = filter(x -> contains(x, key), readdir(path))

    # build src file list -
    list_of_src_files = searchdir(path_to_distribution_files, file_extension)

    # go thru the src file list, and copy the files to the output path -
    for src_file in list_of_src_files

        # create src_buffer -
        src_buffer::Array{String,1} = String[]

        # path to distrubtion -
        path_to_src_file = path_to_distribution_files * "/" * src_file
        open(path_to_src_file, "r") do src_file
            for line in eachline(src_file)

                # need to add a new line for some reason in Julia 0.6
                new_line_with_line_ending = line * "\n"
                push!(src_buffer, new_line_with_line_ending)
            end
        end

        # Write the file to the output -
        path_to_program_file = path_to_output_files * "/" * src_file
        open(path_to_program_file, "w") do f
            for line in src_buffer
                write(f, line)
            end
        end
    end
end

function number_of_species_of_type(list_of_species::Array{SpeciesObject}, species_type::Symbol)

    number_of_species = 0
    for species_object in list_of_species
        local_species_type = species_object.species_type
        if (local_species_type == species_type)
            number_of_species = number_of_species + 1
        end
    end

    return number_of_species
end

function generate_degradation_matrix_buffer(problem_object::ProblemObject)

  # initialize the buffer -
    buffer = ""

  # list of species -
    list_of_species = problem_object.list_of_species

  # get dimension -
    number_of_genes = number_of_species_of_type(list_of_species, :gene)
    number_of_mRNA = number_of_species_of_type(list_of_species, :mrna)
    number_of_proteins = number_of_species_of_type(list_of_species, :protein)
    number_of_degrdation_reactions = number_of_mRNA + number_of_proteins

  # build the gene block -
    for gene_index = 1:number_of_genes
        for reaction_index = 1:number_of_degrdation_reactions
            buffer *= " 0.0 "
        end

        buffer *= "\n"
    end

    for outer_reaction_index = 1:number_of_degrdation_reactions

        for inner_reaction_index = 1:number_of_degrdation_reactions

            if (outer_reaction_index == inner_reaction_index)
                buffer *= " -1.0 "
            else
                buffer *= " 0.0 "
            end
        end

        buffer *= "\n"
    end

  # build the component -
    filename = "Degradation.dat"
    program_component::ProgramComponent = ProgramComponent()
    program_component.filename = filename
    program_component.buffer = buffer

  # return -
    return (program_component)
end

function generate_stoichiomteric_matrix_buffer(problem_object::ProblemObject)

  # list of species -
    list_of_species = problem_object.list_of_species

  # how many mRNA and protein species do we have?
    number_of_genes = number_of_species_of_type(list_of_species, :gene)
    number_of_mRNA = number_of_species_of_type(list_of_species, :mrna)
    number_of_proteins = number_of_species_of_type(list_of_species, :protein)

  # initialize the buffer -
    buffer = ""

  # how many species?
    number_of_species = length(list_of_species)
    for row_species_index = 1:number_of_species

    # what is the species type?
        species_object_row = list_of_species[row_species_index]
        species_type_row = species_object_row.species_type

        if (species_type_row == :gene)

            for txtl_index = 1:(number_of_mRNA + number_of_proteins)
                buffer *= " 0.0 "
            end

      # add a new line -
            buffer *= "\n"
        end
    end

    for row_species_index = 1:(number_of_species - number_of_genes)

        for txtl_index = 1:(number_of_mRNA + number_of_proteins)
            if (row_species_index == txtl_index)
                buffer *= " 1.0 "
            else
                buffer *= " 0.0 "
            end
        end

        buffer *= "\n"
    end


  # build the component -
    filename = "Network.dat"
    program_component::ProgramComponent = ProgramComponent()
    program_component.filename = filename
    program_component.buffer = buffer

  # return -
    return (program_component)
end

function extract_species_of_type(list_of_species::Array{SpeciesObject}, species_type::Symbol)

  # initialize -
    local_list_of_species = SpeciesObject[]
    for species_object in list_of_species

        local_species_type = species_object.species_type
        if (species_type == local_species_type)
            push!(local_list_of_species, species_object)
        end
    end

    return local_list_of_species
end

function is_species_a_target_in_connection_list(list_of_connections::Array{ConnectionObject}, target_species::SpeciesObject, connection_type::Symbol)

    target_connection_list = ConnectionObject[]
    for connection_object in list_of_connections

        @show connection_object

    # get targets -
        local_target_set = connection_object.connection_target_set
        local_connection_type = connection_object.connection_type
        if (local_connection_type == connection_type && does_set_contain_species(local_target_set, target_species) == true)

            push!(target_connection_list, connection_object)
        end
    end

    return target_connection_list
end


function does_set_contain_species(species_set::Array{SpeciesObject}, species_object::SpeciesObject)

  # initilize -
    does_set_contains_species = false
    for local_species_object in species_set

    # get type and name -
        local_species_type = local_species_object.species_type
        local_species_symbol = local_species_object.species_symbol

    # get the test symbol and type -
        test_species_type = species_object.species_type
        test_species_symbol = species_object.species_symbol

        if (local_species_type == test_species_type && test_species_symbol == local_species_symbol)

            @show (local_species_type, local_species_symbol, test_species_type, test_species_symbol)

            does_set_contains_species = true
            return does_set_contains_species
        end
    end

    return does_set_contains_species
end

function generate_dilution_matrix_buffer(problem_object::ProblemObject)

  # list of species -
    list_of_species = problem_object.list_of_species

  # initialize the buffer -
    buffer = ""

  # how many species?
    number_of_species = length(list_of_species)
    for row_species_index = 1:number_of_species

    # what is the species type?
        species_object = list_of_species[row_species_index]
        species_type = species_object.species_type

        for col_species_index = 1:number_of_species

            if (row_species_index == col_species_index)
                if (species_type == :gene)
                    buffer *= " 0.0 "
                else
                    buffer *= " -1.0 "
                end
            else
                buffer *= " 0.0 "
            end
        end # inner for

    # add a new line -
        buffer *= "\n"
    end # outer for

  # build the component -
    filename = "Dilution.dat"
    program_component::ProgramComponent = ProgramComponent()
    program_component.filename = filename
    program_component.buffer = buffer

  # return -
    return (program_component)
end

function partition!(list_of_species::Array{SpeciesObject})

  # ok, frist, we need to split into balanced and unbalanced lists -
    list_of_gene_indexes::Array{Int} = Int[]
    list_of_mRNA_indexes::Array{Int} = Int[]
    list_of_protein_indexes::Array{Int} = Int[]

    for (index, species_object) in enumerate(list_of_species)

        species_type::Symbol = species_object.species_type
        if (species_type == :gene)
            push!(list_of_gene_indexes, index)
        elseif (species_type == :mrna)
      push!(list_of_mRNA_indexes, index)
    elseif (species_type == :protein)
      push!(list_of_protein_indexes, index)
        end
    end

  # combine -
    permutation_index_array = vcat(list_of_gene_indexes, list_of_mRNA_indexes, list_of_protein_indexes)

  # permute the array -
    permute!(list_of_species, permutation_index_array)
end
