function partition!(list_of_species::Array{SpeciesObject})::Array{SpeciesObject}

    genes = filter(s -> s.species_type == :gene, list_of_species)
    mRNAs = filter(s -> s.species_type == :mrna, list_of_species)
    proteins = filter(s -> s.species_type == :protein, list_of_species)

    empty!(list_of_species)
    append!(list_of_species, genes)
    append!(list_of_species, mRNAs)
    append!(list_of_species, proteins)

    return list_of_species
end

function number_of_species_of_type(list_of_species::Array{SpeciesObject}, species_type::Symbol)::Int
    return count(s -> s.species_type == species_type, list_of_species)
end

function extract_species_of_type(list_of_species::Array{SpeciesObject}, species_type::Symbol)::Array{SpeciesObject}
    return filter(s -> s.species_type == species_type, list_of_species)
end

function is_species_a_target_in_connection_list(list_of_connections::Array{ConnectionObject}, target_species::SpeciesObject, connection_type::Symbol)::Array{ConnectionObject}

    matching_connections = ConnectionObject[]
    target_symbol = target_species.species_symbol
    for connection_object in list_of_connections
        if connection_object.connection_type == connection_type
            for target_object in connection_object.connection_target_set
                if target_object.species_symbol == target_symbol
                    push!(matching_connections, connection_object)
                    break
                end
            end
        end
    end

    return matching_connections
end

function is_file_path_ok(path_to_file::String)
    if isfile(path_to_file) == false
        throw(ArgumentError("File path $(path_to_file) does not exist or is not accessible."))
    end
    return nothing
end

function generate_default_project_file(path_to_defaults_file::String)

    default_content = """
    [time]
    start = 0.0
    stop = 10.0
    step = 0.1

    [solver]
    abstol = 1e-9
    reltol = 1e-6
    """

    open(path_to_defaults_file, "w") do f
        write(f, default_content)
    end

    return nothing
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

function move_existing_project_at_path(path_to_existing_project::String)::Bool

    # we are getting called *if* we already know there is a dir conflict -
    # if this is getting called, we have an existing dir where the user wants to write code.
    # we need then create a new dir called *.0, and mv the offending dir to this location?
    # return true if this worked, otherwise false -

    # parent and child dir that we are generating into -
    parent_dir = dirname(path_to_existing_project)
    child_dir = basename(path_to_existing_project)
    destination_path = ""

    # current backup index  -
    current_backup_index = 0

    # do we already have the destination?
    loop_flag = true
    while loop_flag

         # make a destination path -
        destination_path = joinpath(parent_dir, "$(child_dir).$(current_backup_index)")

        # we don't have this dir, we are done -
        if (isdir(destination_path) == false)
            loop_flag = false
        end

        # ok, looks like we already have this dir, update the counter -
        current_backup_index = current_backup_index + 1
    end

    # mv -
    mv(path_to_existing_project, destination_path)

    # check -
    if (isdir(destination_path) == false)
        return false
    end

    return true
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
    list_of_connections = problem_object.list_of_connections

    # how many mRNA and protein species do we have?
    number_of_genes = number_of_species_of_type(list_of_species, :gene)
    number_of_mRNA = number_of_species_of_type(list_of_species, :mrna)
    number_of_proteins = number_of_species_of_type(list_of_species, :protein)

    # count PTM reactions -
    ptm_types = Set([:phosphorylate, :dephosphorylate, :bind, :unbind])
    ptm_connections = filter(c -> c.connection_type in ptm_types, list_of_connections)
    number_of_ptm_reactions = length(ptm_connections)

    # total columns = TXTL rates + PTM rates -
    number_of_txtl_rates = number_of_mRNA + number_of_proteins
    number_of_total_rates = number_of_txtl_rates + number_of_ptm_reactions

    # build species index map -
    species_index_map = Dict{String,Int}()
    for (i, sp) in enumerate(list_of_species)
        species_index_map[sp.species_symbol] = i
    end

    # initialize stoichiometric matrix as zeros -
    number_of_species = length(list_of_species)
    S = zeros(Float64, number_of_species, number_of_total_rates)

    # fill TXTL block (existing logic) -
    # mRNA and protein species get +1 on the diagonal in the TXTL columns
    for row_species_index = 1:(number_of_species - number_of_genes)
        if row_species_index <= number_of_txtl_rates
            S[number_of_genes + row_species_index, row_species_index] = 1.0
        end
    end

    # fill PTM block -
    for (ptm_index, conn) in enumerate(ptm_connections)
        col = number_of_txtl_rates + ptm_index

        if conn.connection_type == :phosphorylate || conn.connection_type == :dephosphorylate
            # substrate is consumed (-1), product is formed (+1), enzyme unchanged -
            # target is the substrate -
            for target_obj in conn.connection_target_set
                substrate_symbol = target_obj.species_symbol
                # look for protein version first, then raw symbol -
                substrate_key = _find_species_key(species_index_map, substrate_symbol)
                if !isnothing(substrate_key)
                    S[species_index_map[substrate_key], col] = -1.0
                end
            end

            # product: need to check the sentence product clause -
            # product symbols are stored as species - find them -
            product_symbols = _get_product_symbols_for_ptm(conn, problem_object)
            for prod_sym in product_symbols
                prod_key = _find_species_key(species_index_map, prod_sym)
                if !isnothing(prod_key)
                    S[species_index_map[prod_key], col] = 1.0
                end
            end

        elseif conn.connection_type == :bind
            # reactants consumed (-1), complex formed (+1) -
            for actor_obj in conn.connection_actor_set
                actor_key = _find_species_key(species_index_map, actor_obj.species_symbol)
                if !isnothing(actor_key)
                    S[species_index_map[actor_key], col] = -1.0
                end
            end
            for target_obj in conn.connection_target_set
                target_key = _find_species_key(species_index_map, target_obj.species_symbol)
                if !isnothing(target_key)
                    S[species_index_map[target_key], col] = 1.0
                end
            end

        elseif conn.connection_type == :unbind
            # complex consumed (-1), components formed (+1) -
            for actor_obj in conn.connection_actor_set
                actor_key = _find_species_key(species_index_map, actor_obj.species_symbol)
                if !isnothing(actor_key)
                    S[species_index_map[actor_key], col] = -1.0
                end
            end
            for target_obj in conn.connection_target_set
                target_key = _find_species_key(species_index_map, target_obj.species_symbol)
                if !isnothing(target_key)
                    S[species_index_map[target_key], col] = 1.0
                end
            end
        end
    end

    # convert matrix to buffer string -
    buffer = ""
    for row = 1:number_of_species
        for col = 1:number_of_total_rates
            buffer *= " $(S[row, col]) "
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

function _find_species_key(species_index_map::Dict{String,Int}, symbol::String)::Union{String,Nothing}
    # try exact match first -
    if haskey(species_index_map, symbol)
        return symbol
    end
    # try protein_ prefix -
    protein_symbol = "protein_" * symbol
    if haskey(species_index_map, protein_symbol)
        return protein_symbol
    end
    return nothing
end

function _get_product_symbols_for_ptm(conn::ConnectionObject, problem_object::ProblemObject)::Vector{String}
    # for PTM connections, the product species are stored in config -
    ptm_products = get(problem_object.configuration_dictionary, "ptm_product_map", Dict{String,Vector{String}}())
    conn_key = conn.connection_symbol * "_" * string(conn.connection_type)
    return get(ptm_products, conn_key, String[])
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

