function build_function_header_buffer(comment_dictionary)

  # initialize -
    buffer = ""

  # get some data from the comment_dictionary -
    function_name = comment_dictionary["function_name"]
    function_description = comment_dictionary["function_description"]
    input_arg_array = comment_dictionary["input_args"]
    output_arg_array = comment_dictionary["output_args"]

    buffer *= "# ----------------------------------------------------------------------------------- #\n"
    buffer *= "# Function: $(function_name)\n"
    buffer *= "# Description: $(function_description)\n"
    buffer *= "# Generated on: $(now())\n"
    buffer *= "#\n"
    buffer *= "# Input arguments:\n"

    for argument_dictionary in input_arg_array

        arg_symbol = argument_dictionary["symbol"]
        arg_description = argument_dictionary["description"]

    # write the buffer -
        buffer *= "# $(arg_symbol) => $(arg_description) \n"
    end

    buffer *= "#\n"
    buffer *= "# Output arguments:\n"
    for argument_dictionary in output_arg_array

        arg_symbol = argument_dictionary["symbol"]
        arg_description = argument_dictionary["description"]

    # write the buffer -
        buffer *= "# $(arg_symbol) => $(arg_description) \n"
    end
    buffer *= "# ----------------------------------------------------------------------------------- #\n"

  # return the buffer -
    return buffer
end

function build_copyright_header_buffer()
    
    # What is the current year?
    current_year = string(Dates.year(now()))

    # Build header -
    buffer = ""
    buffer *= "# ----------------------------------------------------------------------------------- #\n"
    buffer *= "# Copyright (c) $(current_year) Varnerlab\n"
    buffer *= "# Robert Frederick Smith School of Chemical and Biomolecular Engineering\n"
    buffer *= "# Cornell University, Ithaca NY 14850\n"
    buffer *= "#\n"
    buffer *= "# Permission is hereby granted, free of charge, to any person obtaining a copy\n"
    buffer *= "# of this software and associated documentation files (the \"Software\"), to deal\n"
    buffer *= "# in the Software without restriction, including without limitation the rights\n"
    buffer *= "# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
    buffer *= "# copies of the Software, and to permit persons to whom the Software is\n"
    buffer *= "# furnished to do so, subject to the following conditions:\n"
    buffer *= "#\n"
    buffer *= "# The above copyright notice and this permission notice shall be included in\n"
    buffer *= "# all copies or substantial portions of the Software.\n"
    buffer *= "#\n"
    buffer *= "# THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
    buffer *= "# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
    buffer *= "# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
    buffer *= "# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
    buffer *= "# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
    buffer *= "# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n"
    buffer *= "# THE SOFTWARE.\n"
    buffer *= "# ----------------------------------------------------------------------------------- #\n"

    # return -
    return buffer
end

function build_copyright_header_buffer(problem_object::ProblemObject)

  # What is the current year?
    current_year = string(Dates.year(now()))

  # Get comment data from
    buffer = ""
    buffer *= "# ----------------------------------------------------------------------------------- #\n"
    buffer *= "# Copyright (c) $(current_year) Varnerlab\n"
    buffer *= "# Robert Frederick Smith School of Chemical and Biomolecular Engineering\n"
    buffer *= "# Cornell University, Ithaca NY 14850\n"
    buffer *= "#\n"
    buffer *= "# Permission is hereby granted, free of charge, to any person obtaining a copy\n"
    buffer *= "# of this software and associated documentation files (the \"Software\"), to deal\n"
    buffer *= "# in the Software without restriction, including without limitation the rights\n"
    buffer *= "# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
    buffer *= "# copies of the Software, and to permit persons to whom the Software is\n"
    buffer *= "# furnished to do so, subject to the following conditions:\n"
    buffer *= "#\n"
    buffer *= "# The above copyright notice and this permission notice shall be included in\n"
    buffer *= "# all copies or substantial portions of the Software.\n"
    buffer *= "#\n"
    buffer *= "# THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
    buffer *= "# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
    buffer *= "# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
    buffer *= "# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
    buffer *= "# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
    buffer *= "# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n"
    buffer *= "# THE SOFTWARE.\n"
    buffer *= "# ----------------------------------------------------------------------------------- #\n"

  # return -
    return buffer
end

function iterate_binding_control_connection(gene_object::SpeciesObject, list_of_connections::Array{ConnectionObject})

  # What is my gene_symbol -
    gene_symbol = gene_object.species_symbol


    buffer = ""
    for connection_object in list_of_connections

    # actor -
        connection_symbol = connection_object.connection_symbol

        buffer *= "\tn_$(gene_symbol)_$(connection_symbol) = binding_parameter_dictionary[\"n_$(gene_symbol)_$(connection_symbol)\"]\n"
        buffer *= "\tK_$(gene_symbol)_$(connection_symbol) = binding_parameter_dictionary[\"K_$(gene_symbol)_$(connection_symbol)\"]\n"

    end

    return buffer
end

function iterate_control_control_connection(gene_object::SpeciesObject, list_of_connections::Array{ConnectionObject})

  # What is my gene_symbol -
    gene_symbol = gene_object.species_symbol


    buffer = ""
    for connection_object in list_of_connections

    # actor -
        connection_symbol = connection_object.connection_symbol
        buffer *= "\tW_$(gene_symbol)_$(connection_symbol) = control_parameter_dictionary[\"W_$(gene_symbol)_$(connection_symbol)\"]\n"
    end

    return buffer
end

function build_blank_control_buffer(problem_object::ProblemObject)::Union{ProgramComponent,Exception}

  # initialize -
    filename = "Control.jl"

    try

        # build the header -
        header_buffer = build_copyright_header_buffer(problem_object)

        # get the comment buffer -
        comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["transcription_control_function"]
        function_comment_buffer = build_function_header_buffer(comment_header_dictionary)

        # extract list of genes -
        list_of_species = problem_object.list_of_species
        list_of_genes = extract_species_of_type(list_of_species, :gene)

        # initialize the buffer -
        buffer = ""
        buffer *= header_buffer
        buffer *= "#\n"
        buffer *= function_comment_buffer
        buffer *= "function calculate_transcription_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})\n"
        buffer *= "\n"
        buffer *= "\t# initialize the control - \n"
        buffer *= "\tcontrol_array = zeros($(length(list_of_genes)))\n"
        buffer *= "\n"

        buffer *= "\t# Alias the species - \n"
        for (index, species_object) in enumerate(list_of_species)

          # Grab the symbol -
            species_symbol = species_object.species_symbol

            # write the record -
            buffer *= "\t$(species_symbol) = x[$(index)]\n"
        end

        buffer *= "\n"
        buffer *= "\t#TODO: Replave this line with your transcription control logic\n"
        buffer *= "\n"

        buffer *= "\t# return - \n"
        buffer *= "\treturn control_array\n"
        buffer *= "end\n"
        buffer *= "\n"

        # we need to add the translation control -
        # get the comment buffer -
        translation_comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["translation_control_function"]
        translation_function_comment_buffer = build_function_header_buffer(translation_comment_header_dictionary)

        buffer *= "#\n"
        buffer *= translation_function_comment_buffer
        buffer *= "function calculate_translation_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})\n"
        buffer *= "\n"
        buffer *= "\t# initialize the control - \n"
        buffer *= "\tcontrol_array = ones($(length(list_of_genes)))\n"
        buffer *= "\n"
        buffer *= "\t# return - \n"
        buffer *= "\treturn control_array\n"
        buffer *= "end\n"

        # build the component -
        program_component::ProgramComponent = ProgramComponent()
        program_component.filename = filename
        program_component.buffer = buffer

        # return -
        return program_component
    catch error
        rethrow(error)
    end
end

function build_control_buffer(problem_object::ProblemObject)

    # initialize -
    filename = "Control.jl"

    # build the header -
    header_buffer = build_copyright_header_buffer(problem_object)

    # get the comment buffer -
    comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["transcription_control_function"]
    function_comment_buffer = build_function_header_buffer(comment_header_dictionary)

    # extract list of genes -
    list_of_species = problem_object.list_of_species
    list_of_genes = extract_species_of_type(list_of_species, :gene)

    # get list of connections -
    list_of_connections::Array{ConnectionObject} = problem_object.list_of_connections

    # initialize the buffer -
    buffer = ""
    buffer *= header_buffer
    buffer *= "#\n"
    buffer *= function_comment_buffer
    buffer *= "function calculate_transcription_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})\n"
    buffer *= "\n"
    buffer *= "\t# initialize the control - \n"
    buffer *= "\tcontrol_array = zeros($(length(list_of_genes)))\n"
    buffer *= "\n"

    buffer *= "\t# Alias the species - \n"
    for (index, species_object) in enumerate(list_of_species)
        species_symbol = species_object.species_symbol
        buffer *= "\t$(species_symbol) = x[$(index)]\n"
    end
    buffer *= "\n"

    buffer *= "\t# Get the binding and control parameter dictionaries - \n"
    buffer *= "\tbinding_parameter_dictionary = data_dictionary[\"binding_parameter_dictionary\"]\n"
    buffer *= "\tcontrol_parameter_dictionary = data_dictionary[\"control_parameter_dictionary\"]\n"
    buffer *= "\n"

    # generate control logic for each gene -
    for (gene_index, gene_object) in enumerate(list_of_genes)

        gene_symbol = gene_object.species_symbol
        activating_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :activate)
        inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :inhibit)

        # binding parameters -
        all_connections = ConnectionObject[]
        append!(all_connections, activating_connections)
        append!(all_connections, inhibiting_connections)

        if !isempty(all_connections)
            buffer *= iterate_binding_control_connection(gene_object, all_connections)
            buffer *= iterate_control_control_connection(gene_object, all_connections)
            buffer *= "\n"
        end

        # build the control transfer function -
        buffer *= "\t# Transfer function target: $(gene_symbol)\n"
        buffer *= "\tactor_set = [\n"

        for connection_object in activating_connections
            connection_symbol = connection_object.connection_symbol
            buffer *= "\t\tW_$(gene_symbol)_$(connection_symbol)*(protein_$(connection_symbol)^n_$(gene_symbol)_$(connection_symbol))/(K_$(gene_symbol)_$(connection_symbol)^n_$(gene_symbol)_$(connection_symbol)+protein_$(connection_symbol)^n_$(gene_symbol)_$(connection_symbol))\n"
        end

        for connection_object in inhibiting_connections
            connection_symbol = connection_object.connection_symbol
            buffer *= "\t\tW_$(gene_symbol)_$(connection_symbol)*(1 - protein_$(connection_symbol)^n_$(gene_symbol)_$(connection_symbol)/(K_$(gene_symbol)_$(connection_symbol)^n_$(gene_symbol)_$(connection_symbol)+protein_$(connection_symbol)^n_$(gene_symbol)_$(connection_symbol)))\n"
        end

        buffer *= "\t]\n"
        buffer *= "\tW_$(gene_symbol)_RNAP = control_parameter_dictionary[\"W_$(gene_symbol)_RNAP\"]\n"
        buffer *= "\tcontrol_array[$(gene_index)] = (W_$(gene_symbol)_RNAP + sum(actor_set))/(1 + W_$(gene_symbol)_RNAP + sum(actor_set))\n"
        buffer *= "\n"
    end

    buffer *= "\t# return - \n"
    buffer *= "\treturn control_array\n"
    buffer *= "end\n"
    buffer *= "\n"

    # translation control -
    translation_comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["translation_control_function"]
    translation_function_comment_buffer = build_function_header_buffer(translation_comment_header_dictionary)

    buffer *= "#\n"
    buffer *= translation_function_comment_buffer
    buffer *= "function calculate_translation_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})\n"
    buffer *= "\n"
    buffer *= "\t# initialize the control - \n"
    buffer *= "\tcontrol_array = ones($(length(list_of_genes)))\n"
    buffer *= "\n"
    buffer *= "\t# return - \n"
    buffer *= "\treturn control_array\n"
    buffer *= "end\n"

    # build the component -
    program_component::ProgramComponent = ProgramComponent()
    program_component.filename = filename
    program_component.buffer = buffer

    # return -
    return (program_component)
end

# this function needs to swicth based upon the host type -
function build_data_dictionary_buffer(intermediate_representation::Dict{String,Any})::ProgramComponent

    # initialize -
    buffer = Array{String,1}()
    filename = "Data.jl"

    try

        # build the header -
        header_buffer = build_copyright_header_buffer()
        +(buffer,header_buffer; suffix="\n")
        +(buffer,"function build_data_dictionary()::Dict{String,Any}";suffix="\n")
        +(buffer,"\n")
        +(buffer,"# initialize -";prefix="\t",suffix="\n")
        +(buffer,"data_dictionary = Dict{String,Any}()"; prefix="\t",suffix="\n")
        +(buffer,"\n")
        +(buffer,"try"; prefix="\t",suffix="\n")
        +(buffer,"\n")
        +(buffer,"# Load the stoichiometric_matrix -"; prefix="\t\t",suffix="\n")
        +(buffer,"stoichiometric_matrix = readdlm(\"./network/Network.dat\")";prefix="\t\t",suffix="\n")
        +(buffer,"data_dictionary[\"stoichiometric_matrix\"] = stoichiometric_matrix"; suffix="\n",prefix="\t\t")
        +(buffer,"return data_dictionary";prefix="\t\t",suffix="\n")
        +(buffer,"\n")
        +(buffer,"catch error"; prefix="\t",suffix="\n")
        +(buffer,"rethrow(error)"; prefix="\t\t",suffix="\n")
        +(buffer,"end"; prefix="\t",suffix="\n")
        +(buffer,"end";suffix="\n")

        # collapse -
        flat_buffer = ""
        [flat_buffer *= line for line in buffer]

        # build the component -
        program_component::ProgramComponent = ProgramComponent()
        program_component.filename = filename
        program_component.buffer = flat_buffer

        # return -
        return (program_component)
    catch error
        rethrow(error)
    end

end

function build_data_dictionary_buffer(problem_object::ProblemObject, host_flag::Symbol)

    filename = "Data.jl"

    # build the header -
    header_buffer = build_copyright_header_buffer(problem_object)

  # get the comment buffer -
    comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["data_dictionary_function"]
    function_comment_buffer = build_function_header_buffer(comment_header_dictionary)

  # default parameter default dictionary -
    parameter_value_default_dictionary = problem_object.configuration_dictionary["default_parameter_dictionary"]

  # get list of species from the po -
    list_of_species::Array{SpeciesObject} = problem_object.list_of_species

  # initialize the buffer -
    buffer = ""
    buffer *= header_buffer
    buffer *= "#\n"
    buffer *= function_comment_buffer
    buffer *= "function build_data_dictionary(time_span::Tuple{Float64,Float64,Float64}, path_to_biophysical_constants_file::String = \"./Default.json\", host_type::Symbol = :bacteria)::Dict{String,Any}\n"
    buffer *= "\n"
    buffer *= "\t# load the biophysical_constants dictionary \n"
    buffer *= "\tbiophysical_constants_dictionary = build_biophysical_dictionary(path_to_biophysical_constants_file, host_type)\n"
    buffer *= "\n"
    buffer *= "\t# stoichiometric_matrix and dilution_matrix - \n"
    buffer *= "\tstoichiometric_matrix = readdlm(\"./Network.dat\")\n"

  # v1.0: We combine these arrays into a single dilution_degradation_matrix, see below -
  # buffer *= "\tdilution_matrix = readdlm(\"./Dilution.dat\")\n"
  # buffer *= "\tdegradation_matrix = readdlm(\"./Degradation.dat\")\n"
    buffer *= "\n"
    buffer *= "\t# number of states, and rates - \n"
    buffer *= "\t(number_of_states,number_of_rates) = size(stoichiometric_matrix)\n"
    buffer *= "\n"

    buffer *= "\t# array of species types - \n"
    buffer *= "\tspecies_symbol_type_array = [\n"

  # write out the length of genes -
    for (index, species_object) in enumerate(list_of_species)

    # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

    # write the type -
        buffer *= "\t\t:$(species_type)\t;\t# $(index)\t$(species_symbol)\n"
    end
    buffer *= "\t]\n"
    buffer *= "\n"

  # for some whacky reason, we need to add the species_symbol_type_array to the biophysical dictionary -
    buffer *= "\t# we need to store the species symbol array for later - \n"
    buffer *= "\tbiophysical_constants_dictionary[\"species_symbol_type_array\"] = species_symbol_type_array\n"
    buffer *= "\n"

    buffer *= "\t# array of gene lengths - \n"
    buffer *= "\tgene_coding_length_array = [\n"

  # write out the length of genes -
    for (index, species_object) in enumerate(list_of_species)

    # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

        if (species_type == :gene)
            gene_seq_lengths = get(problem_object.configuration_dictionary, "gene_sequence_lengths", Dict{String,Float64}())
            gene_length = get(gene_seq_lengths, species_symbol, 1000.0)
            buffer *= "\t\t$(gene_length)\t;\t# $(index)\t$(species_symbol)\n"
        end
    end
    buffer *= "\t]\n"
    buffer *= "\n"

    buffer *= "\t# array of mRNA coding lengths - \n"
    buffer *= "\tmRNA_coding_length_array = [\n"

  # write out the length of genes -
    counter = 1
    for (index, species_object) in enumerate(list_of_species)

    # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

        if (species_type == :mrna)
            buffer *= "\t\tgene_coding_length_array[$(counter)]\t;\t# $(index)\t$(counter)\t$(species_symbol)\n"
            counter = counter + 1
        end
    end

    buffer *= "\t]\n"
    buffer *= "\n"

    buffer *= "\t# array of mRNA coding lengths - \n"
    buffer *= "\tprotein_coding_length_array = [\n"

  # write out the length of genes -
    counter = 1
    for (index, species_object) in enumerate(list_of_species)

    # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

        if (species_type == :protein)
            protein_seq_lengths = get(problem_object.configuration_dictionary, "protein_sequence_lengths", Dict{String,Float64}())
            if haskey(protein_seq_lengths, species_symbol)
                protein_length = protein_seq_lengths[species_symbol]
                buffer *= "\t\t$(protein_length)\t;\t# $(index)\t$(counter)\t$(species_symbol)\n"
            else
                buffer *= "\t\tround((0.33)*mRNA_coding_length_array[$(counter)])\t;\t# $(index)\t$(counter)\t$(species_symbol)\n"
            end
            counter = counter + 1
        end
    end

    buffer *= "\t]\n"
    buffer *= "\n"

  # ---- THIS LOGIC SWITCHES BASED ON HOST_TYPE FLAG (BELOW) ------------------- #
    if host_flag == :bacteria || host_flag == :mammalian

        buffer *= "\t# array of intracellular gene copy numbers - \n"
        buffer *= "\tgene_abundance_array = [\n"

      # write out the length of genes -
        for (index, species_object) in enumerate(list_of_species)

        # grab the species -
            species_symbol = species_object.species_symbol
            species_type = species_object.species_type

            if (species_type == :gene)
                buffer *= "\t\t2.0\t;\t# (number/cell) $(index)\t$(species_symbol)\n"
            end
        end
        buffer *= "\t]\n"
        buffer *= "\n"
    elseif host_flag == :cell_free

      buffer *= "\t# array of gene concentrations - \n"
      buffer *= "\tgene_abundance_array = [\n"

      # write out the length of genes -
      for (index, species_object) in enumerate(list_of_species)

        # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

        if (species_type == :gene)
          buffer *= "\t\t5.0\t;\t# (nM) $(index)\t$(species_symbol)\n"
        end
      end
      buffer *= "\t]\n"
      buffer *= "\n"
  else
      throw(error("unsupported host_type value. Expected {:bacteria,:mammalian,:cell_free}. Got $(host_flag)"))
    end
  # ---- THIS LOGIC SWITCHES BASED ON HOST_TYPE FLAG (ABOVE) ------------------- #


    buffer *= "\t# initial condition array - \n"
    buffer *= "\tinitial_condition_array = [\n"

  # write out species -
    for (index, species_object) in enumerate(list_of_species)

    # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

        if (species_type == :gene)
            buffer *= "\t\tgene_abundance_array[$(index)]\t;\t# $(index)\t$(species_symbol)\n"
        elseif (species_type == :mrna || species_type == :protein)
      buffer *= "\t\t0.0\t;\t# $(index)\t$(species_symbol)\n"
        end
    end

    buffer *= "\t]\n"
    buffer *= "\n"

    # Setup parameters dictionaries -
    list_of_genes = extract_species_of_type(list_of_species, :gene)

    # get list of connections -
    list_of_connections::Array{ConnectionObject} = problem_object.list_of_connections
    buffer *= "\tbinding_parameter_dictionary = Dict{String,Float64}()\n"
    default_TF_KD_parameter = parameter_value_default_dictionary["default_TF_KD_parameter"]
    for (index, gene_object) in enumerate(list_of_genes)

        # get gene symbol -
        gene_symbol = gene_object.species_symbol

        # connections -
        activating_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :activate)
        inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :inhibit)

        # Iterate through all connections -
        list_of_all_connections = ConnectionObject[]
        append!(list_of_all_connections, activating_connections)
        append!(list_of_all_connections, inhibiting_connections)
        for connection_object in list_of_all_connections

            # connection -
            connection_symbol = connection_object.connection_symbol

            # write the line -
            buffer *= "\tbinding_parameter_dictionary[\"n_$(gene_symbol)_$(connection_symbol)\"] = 1.0\n"
            buffer *= "\tbinding_parameter_dictionary[\"K_$(gene_symbol)_$(connection_symbol)\"] = $(default_TF_KD_parameter)\n"

        end
    end

    buffer *= "\n"
    buffer *= "\t# Alias the control function parameters - \n"
    buffer *= "\tcontrol_parameter_dictionary = Dict{String,Float64}()\n"
    default_W_gene_symbol_RNAP = parameter_value_default_dictionary["default_background_expression_parameter"]
    for (index, gene_object) in enumerate(list_of_genes)

        # get gene symbol -
        gene_symbol = gene_object.species_symbol

        # connections -
        activating_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :activate)
        inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :inhibit)

        # generate an RNAP term -
        buffer *= "\tcontrol_parameter_dictionary[\"W_$(gene_symbol)_RNAP\"] = $(default_W_gene_symbol_RNAP)\n"

        # process all connections -
        list_of_all_connections = ConnectionObject[]
        append!(list_of_all_connections, activating_connections)
        append!(list_of_all_connections, inhibiting_connections)
        for connection_object in list_of_all_connections

            # actor -
            connection_symbol = connection_object.connection_symbol
            buffer *= "\tcontrol_parameter_dictionary[\"W_$(gene_symbol)_$(connection_symbol)\"] = 1.0\n"
        end
    end

    buffer *= "\n"
    buffer *= "\t# degradation modifiers - \n"
    buffer *= "\tdegradation_modifier_array = [\n"

    # write out modifiers -
    for (index, species_object) in enumerate(list_of_species)

        # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

        if (species_type == :gene)
            buffer *= "\t\t0.0\t;\t# $(index)\t$(species_symbol)\n"
        elseif (species_type == :mrna || species_type == :protein)
            buffer *= "\t\t1.0\t;\t# $(index)\t$(species_symbol)\n"
        end
    end

    buffer *= "\t]\n"
    buffer *= "\n"

    buffer *= "\n"
    buffer *= "\t# time constant modifiers - \n"
    buffer *= "\ttime_constant_modifier_array = [\n"

    # write out modifiers -
    for (index, species_object) in enumerate(list_of_species)

        # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

        if (species_type == :gene)
            buffer *= "\t\t0.0\t;\t# $(index)\t$(species_symbol)\n"
        elseif (species_type == :mrna || species_type == :protein)
            buffer *= "\t\t1.0\t;\t# $(index)\t$(species_symbol)\n"
        end
    end

    buffer *= "\t]\n"
    buffer *= "\n"

    # build the degradation matrix -
    buffer *= "\t# Dilution degrdation matrix - \n"
    buffer *= "\tdilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary,species_symbol_type_array,degradation_modifier_array)\n"

    # buffer *= "\n"
    # buffer *= "\t# Precompute the translation parameters - \n"
    # buffer *= "\ttranslation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array, host_type)\n"

    # buffer *= "\n"
    # buffer *= "\t# Precompute the kinetic limit of transcription - \n"
    # buffer *= "\ttranscription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)\n"

    # put the misc dictionary -
    # buffer *= "\n"
    # buffer *= include_function("misc_parameter_dictionary","\t")
    # buffer *= "\n"

    # parameter name mapping -
    # PTM parameters -
    ptm_types = Set([:phosphorylate, :dephosphorylate, :bind, :unbind])
    ptm_connections = filter(c -> c.connection_type in ptm_types, list_of_connections)
    if !isempty(ptm_connections)
        buffer *= "\n"
        buffer *= "\t# Post-translational modification parameters - \n"
        buffer *= "\tptm_parameter_dictionary = Dict{String,Float64}()\n"

        for conn in ptm_connections
            conn_symbol = conn.connection_symbol
            if conn.connection_type == :phosphorylate || conn.connection_type == :dephosphorylate
                buffer *= "\tptm_parameter_dictionary[\"kcat_$(conn_symbol)\"] = 10.0\t# hr^-1\n"
                buffer *= "\tptm_parameter_dictionary[\"Km_$(conn_symbol)\"] = 0.1\t# muM\n"
            elseif conn.connection_type == :bind
                buffer *= "\tptm_parameter_dictionary[\"kf_$(conn_symbol)\"] = 1.0\t# muM^-1*hr^-1\n"
            elseif conn.connection_type == :unbind
                buffer *= "\tptm_parameter_dictionary[\"kr_$(conn_symbol)\"] = 0.1\t# hr^-1\n"
            end
        end
    end

    buffer *= "\n"
    buffer *= "\t# Parameter name index array - \n"
    name_parameter_mapping_buffer = generate_parameter_name_mapping(list_of_genes, list_of_connections::Array{ConnectionObject})
    buffer *= name_parameter_mapping_buffer

    # return block -
    buffer *= "\n"
    buffer *= "\t# =============================== DO NOT EDIT BELOW THIS LINE ============================== #\n"
    buffer *= "\tdata_dictionary = Dict{String,Any}()\n"
    buffer *= "\tdata_dictionary[\"number_of_states\"] = number_of_states\n"
    buffer *= "\tdata_dictionary[\"species_symbol_type_array\"] = species_symbol_type_array\n"
    buffer *= "\tdata_dictionary[\"initial_condition_array\"] = initial_condition_array\n"
    buffer *= "\tdata_dictionary[\"gene_coding_length_array\"] = gene_coding_length_array\n"
    buffer *= "\tdata_dictionary[\"mRNA_coding_length_array\"] = mRNA_coding_length_array\n"
    buffer *= "\tdata_dictionary[\"protein_coding_length_array\"] = protein_coding_length_array\n"

  # v1.0: Moved these constants into biophysical_constants_dictionary -
  # buffer *= "\tdata_dictionary[\"average_transcript_length\"] = average_transcript_length\n"
  # buffer *= "\tdata_dictionary[\"average_protein_length\"] = average_protein_length\n"
  # buffer *= "\tdata_dictionary[\"volume_of_single_cell\"] = V\n"
  # buffer *= "\tdata_dictionary[\"mass_of_single_cell\"] = mass_of_single_cell\n"
  # buffer *= "\tdata_dictionary[\"rnapII_concentration\"] = rnapII_concentration  # muM \n"
  # buffer *= "\tdata_dictionary[\"ribosome_concentration\"] = ribosome_concentration # muM \n"
  # buffer *= "\tdata_dictionary[\"degradation_constant_mRNA\"] = degradation_constant_mRNA  # hr^-1 \n"
  # buffer *= "\tdata_dictionary[\"degradation_constant_protein\"] = degradation_constant_protein  # hr^-1 \n"
  # buffer *= "\tdata_dictionary[\"kcat_transcription\"] = kcat_transcription  # hr^-1 \n"
  # buffer *= "\tdata_dictionary[\"kcat_translation\"] = kcat_translation  # hr^-1 \n"
  # buffer *= "\tdata_dictionary[\"kcat_transcription_initiation\"] = kcat_transcription_initiation  # hr^-1 \n"
  # buffer *= "\tdata_dictionary[\"kcat_translation_initiation\"] = kcat_translation_initiation  # hr^-1 \n"
  # buffer *= "\tdata_dictionary[\"maximum_specific_growth_rate\"] = maximum_specific_growth_rate  # hr^-1 \n"
  # buffer *= "\tdata_dictionary[\"death_rate_constant\"] = death_rate_constant \n"
  # buffer *= "\tdata_dictionary[\"avg_gene_concentration\"] = avg_gene_concentration \n"
  # buffer *= "\tdata_dictionary[\"saturation_constant_transcription\"] = saturation_transcription \n"
  # buffer *= "\tdata_dictionary[\"saturation_constant_translation\"] = saturation_translation \n"
  # buffer *= "\n"

    buffer *= "\tdata_dictionary[\"stoichiometric_matrix\"] = stoichiometric_matrix\n"
    buffer *= "\tdata_dictionary[\"dilution_degradation_matrix\"] = dilution_degradation_matrix\n"
  # buffer *= "\tdata_dictionary[\"degradation_matrix\"] = degradation_matrix\n"
    buffer *= "\tdata_dictionary[\"binding_parameter_dictionary\"] = binding_parameter_dictionary\n"
    buffer *= "\tdata_dictionary[\"control_parameter_dictionary\"] = control_parameter_dictionary\n"
    buffer *= "\tdata_dictionary[\"parameter_name_mapping_array\"] = parameter_name_mapping_array\n"
    buffer *= "\tdata_dictionary[\"transcription_kinetic_limit_array\"] = transcription_kinetic_limit_array\n"
    buffer *= "\tdata_dictionary[\"translation_parameter_array\"] = translation_parameter_array\n"
    buffer *= "\tdata_dictionary[\"degradation_modifier_array\"] = degradation_modifier_array\n"
    buffer *= "\tdata_dictionary[\"time_constant_modifier_array\"] = time_constant_modifier_array\n"
    buffer *= "\tdata_dictionary[\"biophysical_constants_dictionary\"] = biophysical_constants_dictionary\n"

    # include PTM parameters if present -
    if !isempty(ptm_connections)
        buffer *= "\tdata_dictionary[\"ptm_parameter_dictionary\"] = ptm_parameter_dictionary\n"
    end

    buffer *= "\t# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #\n"
    buffer *= "\treturn data_dictionary\n"
    buffer *= "end\n"

  # build the component -
    program_component::ProgramComponent = ProgramComponent()
    program_component.filename = filename
    program_component.buffer = buffer

  # return -
    return (program_component)
end

function build_kinetics_buffer(problem_object::ProblemObject)

    filename = "Kinetics.jl"

  # get list of species from the po -
    list_of_species::Array{SpeciesObject} = problem_object.list_of_species

  # build the header -
    header_buffer = build_copyright_header_buffer(problem_object)

  # initialize the buffer -
    buffer = ""
    buffer *= header_buffer

    comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["calculate_transcription_rates"]
    function_comment_buffer = build_function_header_buffer(comment_header_dictionary)
    buffer *= "#\n"
    buffer *= function_comment_buffer
    buffer *= "function calculate_transcription_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})\n"
    buffer *= "\n"

    buffer *= "\t# Alias the species - \n"
    number_of_genes = 0
    for (index, species_object) in enumerate(list_of_species)

    # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

    # grab -
        if (species_type == :gene)
            buffer *= "\t$(species_symbol) = x[$(index)]\n"
            number_of_genes = number_of_genes + 1
        end
    end
    buffer *= "\n"
    buffer *= "\t# Initialize the transcription rate - \n"
    buffer *= "\ttranscription_rate_array = zeros($(number_of_genes))\n"
    buffer *= "\tKSAT = data_dictionary[\"saturation_constant_transcription\"]\n"
    buffer *= "\tkcat_transcription = data_dictionary[\"kcat_transcription\"]\n"
    buffer *= "\tkcat_transcription_initiation = data_dictionary[\"kcat_transcription_initiation\"]\n"
    buffer *= "\tmugmax = data_dictionary[\"maximum_specific_growth_rate\"]\n"
    buffer *= "\trnapII_concentration = data_dictionary[\"rnapII_concentration\"]\n"
    buffer *= "\taverage_transcript_length = data_dictionary[\"average_transcript_length\"]\n"
    buffer *= "\tgene_coding_length_array = data_dictionary[\"gene_coding_length_array\"]\n"

    buffer *= "\n"
    buffer *= "\t# Populate the transcription rate array - \n"
    counter = 1
    for (index, species_object) in enumerate(list_of_species)

    # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

    # grab -
        if (species_type == :gene)
            buffer *= "\t# Gene: $(species_symbol)\n"
            buffer *= "\tgene_length = gene_coding_length_array[$(index)]\n"
            buffer *= "\tlength_factor = (average_transcript_length/gene_length)\n"
            buffer *= "\tkcat = (kcat_transcription*length_factor*kcat_transcription_initiation)/(kcat_transcription*length_factor+mugmax)\n"
            buffer *= "\ttranscription_rate_array[$(counter)] = kcat*(rnapII_concentration)*(($(species_symbol))/(KSAT+$(species_symbol)))\n"
            buffer *= "\n"
            counter = counter + 1
        end
    end

    buffer *= "\n"
    buffer *= "\t# return transcription_rate_array - \n"
    buffer *= "\treturn transcription_rate_array\n"

    buffer *= "end\n"
    buffer *= "\n"


    # calculate_background_transcription_rates -
    comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["calculate_background_transcription_rates"]
    function_comment_buffer = build_function_header_buffer(comment_header_dictionary)
    buffer *= function_comment_buffer
    buffer *= "function calculate_background_transcription_rates(t::Float64,x::Array{Float64,1},transcription_rate_array::Array{Float64,1},data_dictionary::Dict{String,Any})\n"
    buffer *= "\treturn zeros(length(x))\n"
    buffer *= "end\n"
    buffer *= "\n"
    buffer *= "\n"


    comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["calculate_translation_rates"]
    function_comment_buffer = build_function_header_buffer(comment_header_dictionary)
    buffer *= function_comment_buffer
    buffer *= "function calculate_translation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})\n"
    buffer *= "\n"

    buffer *= "\t# Alias the species - \n"
    number_of_mRNA = 0
    for (index, species_object) in enumerate(list_of_species)

    # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

    # grab -
        if (species_type == :mrna)
            buffer *= "\t$(species_symbol) = x[$(index)]\n"
            number_of_mRNA = number_of_mRNA + 1
        end
    end
    buffer *= "\n"
    buffer *= "\t# Initialize the translation rate - \n"
    buffer *= "\ttranslation_rate_array = zeros($(number_of_mRNA))\n"
    buffer *= "\tKSAT = data_dictionary[\"saturation_constant_translation\"]\n"
    buffer *= "\tkcat_translation = data_dictionary[\"kcat_translation\"]\n"
    buffer *= "\tribosome_concentration = data_dictionary[\"ribosome_concentration\"]\n"
    buffer *= "\taverage_protein_length = data_dictionary[\"average_protein_length\"]\n"
    buffer *= "\tprotein_coding_length_array = data_dictionary[\"protein_coding_length_array\"]\n"
    buffer *= "\n"
    buffer *= "\t# Populate the translation rate array - \n"
    counter = 1
    for (index, species_object) in enumerate(list_of_species)

    # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

    # grab -
        if (species_type == :mrna)
            buffer *= "\t# Transcript: $(species_symbol)\n"
            buffer *= "\tprotein_length = protein_coding_length_array[$(counter)]\n"
            buffer *= "\tscale_factor = (average_protein_length/protein_length)\n"
            buffer *= "\ttranslation_rate_array[$(counter)] = scale_factor*kcat_translation*(ribosome_concentration)*(($(species_symbol))/(KSAT+$(species_symbol)))\n"
            buffer *= "\n"
            counter = counter + 1
        end
    end

    buffer *= "\n"
    buffer *= "\t# return translation array - \n"
    buffer *= "\treturn translation_rate_array\n"
    buffer *= "end\n"
    buffer *= "\n"


  # calculate_mRNA_degradation_rates -
    comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["calculate_mRNA_degradation_rates"]
    function_comment_buffer = build_function_header_buffer(comment_header_dictionary)
    buffer *= function_comment_buffer
    buffer *= "function calculate_mRNA_degradation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})\n"
    buffer *= "\n"

    buffer *= "\t# Alias the species - \n"
    number_of_mRNA = 0
    for (index, species_object) in enumerate(list_of_species)

    # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

    # grab -
        if (species_type == :mrna)
            buffer *= "\t$(species_symbol) = x[$(index)]\n"
            number_of_mRNA = number_of_mRNA + 1
        end
    end
    buffer *= "\n"
    buffer *= "\t# Initialize the degrdation array - \n"
    buffer *= "\tdegradation_rate_array = zeros($(number_of_mRNA))\n"
    buffer *= "\tmRNA_degrdation_constant = data_dictionary[\"degradation_constant_mRNA\"]\n"
    buffer *= "\n"
    buffer *= "\t# Calculate the degradation_rate_array - \n"
    counter = 1
    for (index, species_object) in enumerate(list_of_species)

    # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

    # grab -
        if (species_type == :mrna)
            buffer *= "\tdegradation_rate_array[$(counter)] = (mRNA_degrdation_constant)*$(species_symbol)\n"
            counter = counter + 1
        end
    end
    buffer *= "\n"
    buffer *= "\t# return the degrdation rate array - \n"
    buffer *= "\treturn degradation_rate_array\n"
    buffer *= "end\n"
    buffer *= "\n"

  # calculate_protein_degradation_rates
    comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["calculate_protein_degradation_rates"]
    function_comment_buffer = build_function_header_buffer(comment_header_dictionary)
    buffer *= function_comment_buffer
    buffer *= "function calculate_protein_degradation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})\n"
    buffer *= "\n"
    buffer *= "\t# Alias the species - \n"
    number_of_proteins = 0
    for (index, species_object) in enumerate(list_of_species)

    # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

    # grab -
        if (species_type == :protein)
            buffer *= "\t$(species_symbol) = x[$(index)]\n"
            number_of_proteins = number_of_proteins + 1
        end
    end
    buffer *= "\n"
    buffer *= "\t# Initialize the degrdation array - \n"
    buffer *= "\tdegradation_rate_array = zeros($(number_of_proteins))\n"
    buffer *= "\tprotein_degrdation_constant = data_dictionary[\"degradation_constant_protein\"]\n"
    buffer *= "\n"
    buffer *= "\t# Calculate the degradation_rate_array - \n"
    counter = 1
    for (index, species_object) in enumerate(list_of_species)

    # grab the species -
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type

    # grab -
        if (species_type == :protein)
            buffer *= "\tdegradation_rate_array[$(counter)] = (protein_degrdation_constant)*$(species_symbol)\n"
            counter = counter + 1
        end
    end
    buffer *= "\n"
    buffer *= "\t# return the degrdation rate array - \n"
    buffer *= "\treturn degradation_rate_array\n"
    buffer *= "end\n"
    buffer *= "\n"
    buffer *= "\n"

  # build the component -
    program_component::ProgramComponent = ProgramComponent()
    program_component.filename = filename
    program_component.buffer = buffer

  # return -
    return (program_component)
end

function build_ptm_kinetics_buffer(problem_object::ProblemObject)::Union{ProgramComponent,Nothing}

    ptm_buffer = build_ptm_kinetics_section(problem_object)
    if isempty(ptm_buffer)
        return nothing
    end

    filename = "PTM.jl"

    header_buffer = build_copyright_header_buffer(problem_object)
    buffer = ""
    buffer *= header_buffer
    buffer *= ptm_buffer

    program_component = ProgramComponent()
    program_component.filename = filename
    program_component.buffer = buffer

    return program_component
end

function build_ptm_kinetics_section(problem_object::ProblemObject)::String

    list_of_connections = problem_object.list_of_connections
    list_of_species = problem_object.list_of_species
    ptm_types = Set([:phosphorylate, :dephosphorylate, :bind, :unbind])
    ptm_connections = filter(c -> c.connection_type in ptm_types, list_of_connections)

    if isempty(ptm_connections)
        return ""
    end

    buffer = ""
    buffer *= "# ----------------------------------------------------------------------------------- #\n"
    buffer *= "# Function: calculate_ptm_rates\n"
    buffer *= "# Description: Calculate post-translational modification rates at time t\n"
    buffer *= "# ----------------------------------------------------------------------------------- #\n"
    buffer *= "function calculate_ptm_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})\n"
    buffer *= "\n"

    # alias all species -
    buffer *= "\t# Alias the species - \n"
    for (index, species_object) in enumerate(list_of_species)
        species_symbol = species_object.species_symbol
        buffer *= "\t$(species_symbol) = x[$(index)]\n"
    end
    buffer *= "\n"

    buffer *= "\t# Get PTM kinetic parameters - \n"
    buffer *= "\tptm_parameter_dictionary = data_dictionary[\"ptm_parameter_dictionary\"]\n"
    buffer *= "\n"

    buffer *= "\t# Initialize the PTM rate array - \n"
    buffer *= "\tptm_rate_array = zeros($(length(ptm_connections)))\n"
    buffer *= "\n"

    for (ptm_index, conn) in enumerate(ptm_connections)

        conn_symbol = conn.connection_symbol

        if conn.connection_type == :phosphorylate || conn.connection_type == :dephosphorylate
            # Michaelis-Menten: v = kcat * E * S / (Km + S)
            enzyme_symbol = _resolve_species_symbol(conn.connection_actor_set[1], list_of_species)
            substrate_symbol = _resolve_species_symbol(conn.connection_target_set[1], list_of_species)
            label = conn.connection_type == :phosphorylate ? "Phosphorylation" : "Dephosphorylation"

            buffer *= "\t# $(label): $(conn_symbol) -> $(conn.connection_type)\n"
            buffer *= "\tkcat_$(conn_symbol) = ptm_parameter_dictionary[\"kcat_$(conn_symbol)\"]\n"
            buffer *= "\tKm_$(conn_symbol) = ptm_parameter_dictionary[\"Km_$(conn_symbol)\"]\n"
            buffer *= "\tptm_rate_array[$(ptm_index)] = kcat_$(conn_symbol)*($(enzyme_symbol))*(($(substrate_symbol))/(Km_$(conn_symbol)+$(substrate_symbol)))\n"
            buffer *= "\n"

        elseif conn.connection_type == :bind
            # mass action: v = kf * A * B * ...
            buffer *= "\t# Binding: $(conn_symbol)\n"
            buffer *= "\tkf_$(conn_symbol) = ptm_parameter_dictionary[\"kf_$(conn_symbol)\"]\n"
            rate_expr = "kf_$(conn_symbol)"
            for actor_obj in conn.connection_actor_set
                actor_sym = _resolve_species_symbol(actor_obj, list_of_species)
                rate_expr *= "*($(actor_sym))"
            end
            buffer *= "\tptm_rate_array[$(ptm_index)] = $(rate_expr)\n"
            buffer *= "\n"

        elseif conn.connection_type == :unbind
            # mass action: v = kr * Complex
            buffer *= "\t# Unbinding: $(conn_symbol)\n"
            buffer *= "\tkr_$(conn_symbol) = ptm_parameter_dictionary[\"kr_$(conn_symbol)\"]\n"
            complex_sym = _resolve_species_symbol(conn.connection_actor_set[1], list_of_species)
            buffer *= "\tptm_rate_array[$(ptm_index)] = kr_$(conn_symbol)*($(complex_sym))\n"
            buffer *= "\n"
        end
    end

    buffer *= "\t# return - \n"
    buffer *= "\treturn ptm_rate_array\n"
    buffer *= "end\n"
    buffer *= "\n"

    return buffer
end

function _resolve_species_symbol(species_obj::SpeciesObject, list_of_species::Array{SpeciesObject})::String
    # try to find the exact symbol in the species list -
    symbol = species_obj.species_symbol
    for sp in list_of_species
        if sp.species_symbol == symbol
            return symbol
        end
    end
    # try with protein_ prefix -
    protein_symbol = "protein_" * symbol
    for sp in list_of_species
        if sp.species_symbol == protein_symbol
            return protein_symbol
        end
    end
    return symbol
end

function build_inputs_buffer(problem_object::ProblemObject)

    filename = "Inputs.jl"

  # build the header -
    license_header_buffer = build_copyright_header_buffer(problem_object)

  # get the comment buffer -
    comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["input_function"]
    function_comment_buffer = build_function_header_buffer(comment_header_dictionary)

  # initialize the buffer -
    buffer = ""
    buffer *= license_header_buffer
    buffer *= "#\n"
    buffer *= function_comment_buffer
    buffer *= "function calculate_input_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})\n"
    buffer *= "\n"
    buffer *= "\t# Initialize default - \n"
    buffer *= "\tu_array = zeros(length(x))\n"
    buffer *= "\n"
    buffer *= "\t# return - \n"
    buffer *= "\treturn u_array\n"
    buffer *= "end\n"

  # build the component -
    program_component::ProgramComponent = ProgramComponent()
    program_component.filename = filename
    program_component.buffer = buffer

  # return -
    return (program_component)
end

function generate_parameter_name_mapping(list_of_genes::Array{SpeciesObject}, list_of_connections::Array{ConnectionObject})

  # iterate through the list of genes - add names to a string -
    list_of_names = String[]

    for (index, gene_object) in enumerate(list_of_genes)

    # get gene symbol -
        gene_symbol = gene_object.species_symbol

    # connections -
        activating_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :activate)
        inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :inhibit)

    # grab the list of connections -
        list_of_all_connections = ConnectionObject[]
        append!(list_of_all_connections, activating_connections)
        append!(list_of_all_connections, inhibiting_connections)
        for connection_object in list_of_all_connections

      # actor -
            connection_symbol = connection_object.connection_symbol

            tmp_name = "n_$(gene_symbol)_$(connection_symbol)"
            push!(list_of_names, tmp_name)

            tmp_name = "K_$(gene_symbol)_$(connection_symbol)"
            push!(list_of_names, tmp_name)
        end
    end

    for (index, gene_object) in enumerate(list_of_genes)

    # get gene symbol -
        gene_symbol = gene_object.species_symbol

    # connections -
        activating_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :activate)
        inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :inhibit)

    # generate an RNAP term -
        tmp_name = "W_$(gene_symbol)_RNAP"
        push!(list_of_names, tmp_name)

    # grab the list of connections -
        list_of_all_connections = ConnectionObject[]
        append!(list_of_all_connections, activating_connections)
        append!(list_of_all_connections, inhibiting_connections)
        for connection_object in list_of_all_connections

      # actor -
            connection_symbol = connection_object.connection_symbol
            tmp_name = "W_$(gene_symbol)_$(connection_symbol)"
            push!(list_of_names, tmp_name)
        end
    end

  # push some "hard coded" constants onto the end -
  # we should encode this list in a config file ...
  # "rnapII_concentration"	; 				# 1
  # "ribosome_concentration"	;				# 2
  # "degradation_constant_mRNA"	;			# 3
  # "degradation_constant_protein"	;	# 4
  # "kcat_transcription"	;						# 5
  # "kcat_translation"	;							# 6
  # "maximum_specific_growth_rate"	;	# 7
  # "saturation_constant_transcription"	;	# 8
  # "saturation_constant_translation"	;	# 9
    push!(list_of_names, "rnapII_concentration")
    push!(list_of_names, "ribosome_concentration")
    push!(list_of_names, "degradation_constant_mRNA")
    push!(list_of_names, "degradation_constant_protein")
    push!(list_of_names, "kcat_transcription")
    push!(list_of_names, "kcat_translation")
    push!(list_of_names, "maximum_specific_growth_rate")
    push!(list_of_names, "saturation_constant_transcription")
    push!(list_of_names, "saturation_constant_translation")

    buffer = ""
    buffer *= "\tparameter_name_mapping_array = [\n"
    for (index, symbol) in enumerate(list_of_names)

        buffer *= "\t\t\"$(symbol)\"\t;\t# $(index)\n"

    end
    buffer *= "\t]\n"

  # return -
    return buffer
end

# ====================================================================================== #
# Effective biophysical model builders (Adhikari et al., 2020)
# ====================================================================================== #

function build_kinetics_buffer_effective(problem_object::ProblemObject)

    filename = "Kinetics.jl"

    # get list of species from the po -
    list_of_species::Array{SpeciesObject} = problem_object.list_of_species

    # build the header -
    header_buffer = build_copyright_header_buffer(problem_object)

    # count genes and mRNAs -
    list_of_genes = extract_species_of_type(list_of_species, :gene)
    list_of_mRNAs = extract_species_of_type(list_of_species, :mrna)
    list_of_proteins = extract_species_of_type(list_of_species, :protein)
    number_of_genes = length(list_of_genes)

    # build the index arrays -
    mRNA_indices = Int[]
    protein_indices = Int[]
    gene_indices = Int[]
    for (index, species_object) in enumerate(list_of_species)
        if species_object.species_type == :mrna
            push!(mRNA_indices, index)
        elseif species_object.species_type == :protein
            push!(protein_indices, index)
        elseif species_object.species_type == :gene
            push!(gene_indices, index)
        end
    end

    # initialize the buffer -
    buffer = ""
    buffer *= header_buffer
    buffer *= "# ====================================================================================== #\n"
    buffer *= "# Effective biophysical TXTL kinetics (Adhikari et al., 2020)\n"
    buffer *= "# ====================================================================================== #\n"
    buffer *= "\n"

    # Main function: calculate_txtl_kinetics_array -
    buffer *= "function calculate_txtl_kinetics_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})::Array{Float64,1}\n"
    buffer *= "\n"
    buffer *= "\t# initialize - \n"
    buffer *= "\tcalculated_txtl_kinetics_array = Float64[]\n"
    buffer *= "\n"
    buffer *= "\t# Get effective model parameters from data dictionary - \n"
    buffer *= "\tR_XT = data_dictionary[\"R_XT\"]\t# RNAP concentration (muM)\n"
    buffer *= "\tR_LT = data_dictionary[\"R_LT\"]\t# ribosome concentration (muM)\n"
    buffer *= "\tv_X = data_dictionary[\"v_X\"]\t# transcription elongation rate (nt/s)\n"
    buffer *= "\tv_L = data_dictionary[\"v_L\"]\t# translation elongation rate (aa/s)\n"
    buffer *= "\tK_X = data_dictionary[\"K_X\"]\t# transcription saturation constant (muM)\n"
    buffer *= "\tK_L = data_dictionary[\"K_L\"]\t# translation saturation constant (muM)\n"
    buffer *= "\tK_P = data_dictionary[\"K_P\"]\t# polysome amplification factor\n"
    buffer *= "\ttau_X_array = data_dictionary[\"tau_X_array\"]\t# TX dimensionless time constants\n"
    buffer *= "\ttau_L_array = data_dictionary[\"tau_L_array\"]\t# TL dimensionless time constants\n"
    buffer *= "\tgene_coding_length_array = data_dictionary[\"gene_coding_length_array\"]\n"
    buffer *= "\tprotein_coding_length_array = data_dictionary[\"protein_coding_length_array\"]\n"
    buffer *= "\tspecies_symbol_type_array = data_dictionary[\"species_symbol_type_array\"]\n"
    buffer *= "\n"

    # Gene concentrations -
    buffer *= "\t# Gene concentrations - \n"
    buffer *= "\tgene_abundance_array = Float64[]\n"
    for (i, gi) in enumerate(gene_indices)
        buffer *= "\tpush!(gene_abundance_array, x[$(gi)])\t# $(list_of_genes[i].species_symbol)\n"
    end
    buffer *= "\n"

    # mRNA concentrations -
    buffer *= "\t# mRNA concentrations - \n"
    buffer *= "\tmRNA_concentration_array = Float64[]\n"
    for (i, mi) in enumerate(mRNA_indices)
        buffer *= "\tpush!(mRNA_concentration_array, x[$(mi)])\t# $(list_of_mRNAs[i].species_symbol)\n"
    end
    buffer *= "\n"

    # Number of genes -
    buffer *= "\t# number of genes - \n"
    buffer *= "\tnumber_of_genes = $(number_of_genes)\n"
    buffer *= "\n"

    # Convert elongation rates to per-hour -
    buffer *= "\t# Convert elongation rates to per-hour - \n"
    buffer *= "\tv_X_hr = v_X * 3600.0\t# nt/hr\n"
    buffer *= "\tv_L_hr = v_L * 3600.0\t# aa/hr\n"
    buffer *= "\n"

    # Transcription rates (Eq. 3 with competition Eq. 4) -
    buffer *= "\t# ---- Transcription kinetic limit (Eq. 3) with competition (Eq. 4) ---- #\n"
    buffer *= "\ttranscription_rate_array = zeros(number_of_genes)\n"
    buffer *= "\tfor j in 1:number_of_genes\n"
    buffer *= "\n"
    buffer *= "\t\t# gene concentration - \n"
    buffer *= "\t\tG_j = gene_abundance_array[j]\n"
    buffer *= "\t\tl_G_j = gene_coding_length_array[j]\n"
    buffer *= "\t\ttau_X_j = tau_X_array[j]\n"
    buffer *= "\n"
    buffer *= "\t\t# Vmax for transcription (Eq. 7) - \n"
    buffer *= "\t\tV_X_max_j = R_XT * (v_X_hr / l_G_j)\n"
    buffer *= "\n"
    buffer *= "\t\t# Competition term (Eq. 4) - \n"
    buffer *= "\t\tO_X_j = 0.0\n"
    buffer *= "\t\tfor i in 1:number_of_genes\n"
    buffer *= "\t\t\tif i != j\n"
    buffer *= "\t\t\t\ttau_X_i = tau_X_array[i]\n"
    buffer *= "\t\t\t\tG_i = gene_abundance_array[i]\n"
    buffer *= "\t\t\t\tO_X_j += (tau_X_j / tau_X_i) * (1.0 + tau_X_i) * G_i\n"
    buffer *= "\t\t\tend\n"
    buffer *= "\t\tend\n"
    buffer *= "\n"
    buffer *= "\t\t# Transcription rate (Eq. 3) - \n"
    buffer *= "\t\ttranscription_rate_array[j] = V_X_max_j * G_j / (tau_X_j * K_X + (1.0 + tau_X_j) * G_j + O_X_j)\n"
    buffer *= "\tend\n"
    buffer *= "\n"

    # Translation rates (Eq. 5 with competition Eq. 6) -
    buffer *= "\t# ---- Translation kinetic limit (Eq. 5) with competition (Eq. 6) ---- #\n"
    buffer *= "\ttranslation_rate_array = zeros(number_of_genes)\n"
    buffer *= "\tfor j in 1:number_of_genes\n"
    buffer *= "\n"
    buffer *= "\t\t# mRNA concentration - \n"
    buffer *= "\t\tm_j = mRNA_concentration_array[j]\n"
    buffer *= "\t\tl_P_j = protein_coding_length_array[j]\n"
    buffer *= "\t\ttau_L_j = tau_L_array[j]\n"
    buffer *= "\n"
    buffer *= "\t\t# Vmax for translation (Eq. 8) - \n"
    buffer *= "\t\tV_L_max_j = K_P * R_LT * (v_L_hr / l_P_j)\n"
    buffer *= "\n"
    buffer *= "\t\t# Competition term (Eq. 6) - \n"
    buffer *= "\t\tO_L_j = 0.0\n"
    buffer *= "\t\tfor i in 1:number_of_genes\n"
    buffer *= "\t\t\tif i != j\n"
    buffer *= "\t\t\t\ttau_L_i = tau_L_array[i]\n"
    buffer *= "\t\t\t\tm_i = mRNA_concentration_array[i]\n"
    buffer *= "\t\t\t\tO_L_j += (tau_L_j / tau_L_i) * (1.0 + tau_L_i) * m_i\n"
    buffer *= "\t\t\tend\n"
    buffer *= "\t\tend\n"
    buffer *= "\n"
    buffer *= "\t\t# Translation rate (Eq. 5) - \n"
    buffer *= "\t\ttranslation_rate_array[j] = V_L_max_j * m_j / (tau_L_j * K_L + (1.0 + tau_L_j) * m_j + O_L_j)\n"
    buffer *= "\tend\n"
    buffer *= "\n"

    # Package result -
    buffer *= "\t# Package: [TX_rates; TL_rates] - \n"
    buffer *= "\tfor value in transcription_rate_array\n"
    buffer *= "\t\tpush!(calculated_txtl_kinetics_array, value)\n"
    buffer *= "\tend\n"
    buffer *= "\tfor value in translation_rate_array\n"
    buffer *= "\t\tpush!(calculated_txtl_kinetics_array, value)\n"
    buffer *= "\tend\n"
    buffer *= "\n"
    buffer *= "\t# return - TX rates, then TL rates - \n"
    buffer *= "\treturn calculated_txtl_kinetics_array\n"
    buffer *= "end\n"

    # build the component -
    program_component::ProgramComponent = ProgramComponent()
    program_component.filename = filename
    program_component.buffer = buffer

    # return -
    return (program_component)
end

function build_control_buffer_effective(problem_object::ProblemObject)

    filename = "Control.jl"

    # build the header -
    header_buffer = build_copyright_header_buffer(problem_object)

    # extract lists -
    list_of_species = problem_object.list_of_species
    list_of_genes = extract_species_of_type(list_of_species, :gene)
    list_of_connections::Array{ConnectionObject} = problem_object.list_of_connections

    number_of_genes = length(list_of_genes)

    # initialize the buffer -
    buffer = ""
    buffer *= header_buffer
    buffer *= "# ====================================================================================== #\n"
    buffer *= "# Effective biophysical control functions (Adhikari et al., 2020)\n"
    buffer *= "# ====================================================================================== #\n"
    buffer *= "\n"

    # --- Transcription control (Eq. 9): thermodynamic/Boltzmann formulation --- #
    buffer *= "function calculate_transcription_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})\n"
    buffer *= "\n"
    buffer *= "\t# initialize the control - \n"
    buffer *= "\tcontrol_array = zeros($(number_of_genes))\n"
    buffer *= "\n"

    buffer *= "\t# Alias the species - \n"
    for (index, species_object) in enumerate(list_of_species)
        species_symbol = species_object.species_symbol
        buffer *= "\t$(species_symbol) = x[$(index)]\n"
    end
    buffer *= "\n"

    buffer *= "\t# Get the thermodynamic parameter dictionaries - \n"
    buffer *= "\tdG_dictionary = data_dictionary[\"dG_dictionary\"]\n"
    buffer *= "\tbinding_parameter_dictionary = data_dictionary[\"binding_parameter_dictionary\"]\n"
    buffer *= "\ttemperature = data_dictionary[\"temperature\"]\t# K\n"
    buffer *= "\tRT = 8.314e-3 * temperature\t# kJ/mol\n"
    buffer *= "\n"

    # generate control logic for each gene -
    for (gene_index, gene_object) in enumerate(list_of_genes)

        gene_symbol = gene_object.species_symbol
        activating_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :activate)
        inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :inhibit)

        buffer *= "\t# ---- Control for gene $(gene_index): $(gene_symbol) ---- #\n"

        # RNAP Boltzmann weight -
        buffer *= "\tdG_$(gene_symbol)_RNAP = dG_dictionary[\"dG_$(gene_symbol)_RNAP\"]\n"
        buffer *= "\tW_$(gene_symbol)_RNAP = exp(-dG_$(gene_symbol)_RNAP / RT)\n"

        # numerator starts with RNAP alone -
        buffer *= "\tnumerator_$(gene_index) = W_$(gene_symbol)_RNAP\n"
        buffer *= "\tdenominator_$(gene_index) = 1.0 + W_$(gene_symbol)_RNAP\n"

        # activators: contribute to both numerator and denominator -
        for connection_object in activating_connections
            connection_symbol = connection_object.connection_symbol

            buffer *= "\tn_$(gene_symbol)_$(connection_symbol) = binding_parameter_dictionary[\"n_$(gene_symbol)_$(connection_symbol)\"]\n"
            buffer *= "\tK_$(gene_symbol)_$(connection_symbol) = binding_parameter_dictionary[\"K_$(gene_symbol)_$(connection_symbol)\"]\n"
            buffer *= "\tdG_$(gene_symbol)_$(connection_symbol) = dG_dictionary[\"dG_$(gene_symbol)_$(connection_symbol)\"]\n"
            buffer *= "\tW_$(gene_symbol)_$(connection_symbol) = exp(-dG_$(gene_symbol)_$(connection_symbol) / RT)\n"
            buffer *= "\tf_$(gene_symbol)_$(connection_symbol) = (protein_$(connection_symbol)^n_$(gene_symbol)_$(connection_symbol))/(K_$(gene_symbol)_$(connection_symbol)^n_$(gene_symbol)_$(connection_symbol)+protein_$(connection_symbol)^n_$(gene_symbol)_$(connection_symbol))\n"
            buffer *= "\tnumerator_$(gene_index) += W_$(gene_symbol)_$(connection_symbol)*f_$(gene_symbol)_$(connection_symbol)\n"
            buffer *= "\tdenominator_$(gene_index) += W_$(gene_symbol)_$(connection_symbol)*f_$(gene_symbol)_$(connection_symbol)\n"
        end

        # repressors: contribute only to denominator -
        for connection_object in inhibiting_connections
            connection_symbol = connection_object.connection_symbol

            buffer *= "\tn_$(gene_symbol)_$(connection_symbol) = binding_parameter_dictionary[\"n_$(gene_symbol)_$(connection_symbol)\"]\n"
            buffer *= "\tK_$(gene_symbol)_$(connection_symbol) = binding_parameter_dictionary[\"K_$(gene_symbol)_$(connection_symbol)\"]\n"
            buffer *= "\tdG_$(gene_symbol)_$(connection_symbol) = dG_dictionary[\"dG_$(gene_symbol)_$(connection_symbol)\"]\n"
            buffer *= "\tW_$(gene_symbol)_$(connection_symbol) = exp(-dG_$(gene_symbol)_$(connection_symbol) / RT)\n"
            buffer *= "\tf_$(gene_symbol)_$(connection_symbol) = (protein_$(connection_symbol)^n_$(gene_symbol)_$(connection_symbol))/(K_$(gene_symbol)_$(connection_symbol)^n_$(gene_symbol)_$(connection_symbol)+protein_$(connection_symbol)^n_$(gene_symbol)_$(connection_symbol))\n"
            buffer *= "\tdenominator_$(gene_index) += W_$(gene_symbol)_$(connection_symbol)*f_$(gene_symbol)_$(connection_symbol)\n"
        end

        buffer *= "\tcontrol_array[$(gene_index)] = numerator_$(gene_index) / denominator_$(gene_index)\n"
        buffer *= "\n"
    end

    buffer *= "\t# return - \n"
    buffer *= "\treturn control_array\n"
    buffer *= "end\n"
    buffer *= "\n"

    # --- Translation control (Eq. 10): exponential decay of translation capacity --- #
    buffer *= "# ---- Translation control (Eq. 10): exponential decay of translation capacity ---- #\n"
    buffer *= "function calculate_translation_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})\n"
    buffer *= "\n"
    buffer *= "\t# initialize the control - \n"
    buffer *= "\tcontrol_array = ones($(number_of_genes))\n"
    buffer *= "\n"
    buffer *= "\t# Get the translation capacity half-life - \n"
    buffer *= "\ttau_L_half = data_dictionary[\"tau_L_half\"]\t# hr\n"
    buffer *= "\n"
    buffer *= "\t# Compute the translation capacity decay (Eq. 10) - \n"
    buffer *= "\tepsilon = exp(-0.693 * t / tau_L_half)\n"
    buffer *= "\n"
    buffer *= "\t# Apply to all genes - \n"
    buffer *= "\tfor j in 1:$(number_of_genes)\n"
    buffer *= "\t\tcontrol_array[j] = epsilon\n"
    buffer *= "\tend\n"
    buffer *= "\n"
    buffer *= "\t# return - \n"
    buffer *= "\treturn control_array\n"
    buffer *= "end\n"

    # build the component -
    program_component::ProgramComponent = ProgramComponent()
    program_component.filename = filename
    program_component.buffer = buffer

    # return -
    return (program_component)
end

function build_data_dictionary_buffer_effective(problem_object::ProblemObject, host_flag::Symbol)

    filename = "Data.jl"

    # build the header -
    header_buffer = build_copyright_header_buffer(problem_object)

    # get the comment buffer -
    comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["data_dictionary_function"]
    function_comment_buffer = build_function_header_buffer(comment_header_dictionary)

    # get list of species from the po -
    list_of_species::Array{SpeciesObject} = problem_object.list_of_species

    # initialize the buffer -
    buffer = ""
    buffer *= header_buffer
    buffer *= "#\n"
    buffer *= function_comment_buffer
    buffer *= "function build_data_dictionary(time_span::Tuple{Float64,Float64,Float64}, path_to_biophysical_constants_file::String = \"./Default.json\", host_type::Symbol = :cell_free)::Dict{String,Any}\n"
    buffer *= "\n"
    buffer *= "\t# load the biophysical_constants dictionary \n"
    buffer *= "\tbiophysical_constants_dictionary = build_biophysical_dictionary(path_to_biophysical_constants_file, host_type)\n"
    buffer *= "\n"
    buffer *= "\t# stoichiometric_matrix and dilution_matrix - \n"
    buffer *= "\tstoichiometric_matrix = readdlm(\"./Network.dat\")\n"
    buffer *= "\n"
    buffer *= "\t# number of states, and rates - \n"
    buffer *= "\t(number_of_states,number_of_rates) = size(stoichiometric_matrix)\n"
    buffer *= "\n"

    # species type array -
    buffer *= "\t# array of species types - \n"
    buffer *= "\tspecies_symbol_type_array = [\n"
    for (index, species_object) in enumerate(list_of_species)
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type
        buffer *= "\t\t:$(species_type)\t;\t# $(index)\t$(species_symbol)\n"
    end
    buffer *= "\t]\n"
    buffer *= "\n"

    buffer *= "\t# we need to store the species symbol array for later - \n"
    buffer *= "\tbiophysical_constants_dictionary[\"species_symbol_type_array\"] = species_symbol_type_array\n"
    buffer *= "\n"

    # gene coding length array -
    buffer *= "\t# array of gene lengths - \n"
    buffer *= "\tgene_coding_length_array = [\n"
    for (index, species_object) in enumerate(list_of_species)
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type
        if (species_type == :gene)
            gene_seq_lengths = get(problem_object.configuration_dictionary, "gene_sequence_lengths", Dict{String,Float64}())
            gene_length = get(gene_seq_lengths, species_symbol, 1000.0)
            buffer *= "\t\t$(gene_length)\t;\t# $(index)\t$(species_symbol)\n"
        end
    end
    buffer *= "\t]\n"
    buffer *= "\n"

    # mRNA coding length -
    buffer *= "\t# array of mRNA coding lengths - \n"
    buffer *= "\tmRNA_coding_length_array = [\n"
    counter = 1
    for (index, species_object) in enumerate(list_of_species)
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type
        if (species_type == :mrna)
            buffer *= "\t\tgene_coding_length_array[$(counter)]\t;\t# $(index)\t$(counter)\t$(species_symbol)\n"
            counter = counter + 1
        end
    end
    buffer *= "\t]\n"
    buffer *= "\n"

    # protein coding length -
    buffer *= "\t# array of protein coding lengths - \n"
    buffer *= "\tprotein_coding_length_array = [\n"
    counter = 1
    for (index, species_object) in enumerate(list_of_species)
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type
        if (species_type == :protein)
            protein_seq_lengths = get(problem_object.configuration_dictionary, "protein_sequence_lengths", Dict{String,Float64}())
            if haskey(protein_seq_lengths, species_symbol)
                protein_length = protein_seq_lengths[species_symbol]
                buffer *= "\t\t$(protein_length)\t;\t# $(index)\t$(counter)\t$(species_symbol)\n"
            else
                buffer *= "\t\tround((0.33)*mRNA_coding_length_array[$(counter)])\t;\t# $(index)\t$(counter)\t$(species_symbol)\n"
            end
            counter = counter + 1
        end
    end
    buffer *= "\t]\n"
    buffer *= "\n"

    # gene abundance array -
    buffer *= "\t# array of gene concentrations (muM) - \n"
    buffer *= "\tgene_abundance_array = [\n"
    for (index, species_object) in enumerate(list_of_species)
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type
        if (species_type == :gene)
            buffer *= "\t\t0.005\t;\t# (muM) $(index)\t$(species_symbol)\n"
        end
    end
    buffer *= "\t]\n"
    buffer *= "\n"

    # initial condition -
    buffer *= "\t# initial condition array - \n"
    buffer *= "\tinitial_condition_array = [\n"
    for (index, species_object) in enumerate(list_of_species)
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type
        if (species_type == :gene)
            buffer *= "\t\tgene_abundance_array[$(index)]\t;\t# $(index)\t$(species_symbol)\n"
        elseif (species_type == :mrna || species_type == :protein)
            buffer *= "\t\t0.0\t;\t# $(index)\t$(species_symbol)\n"
        end
    end
    buffer *= "\t]\n"
    buffer *= "\n"

    # --- Effective biophysical model parameters (Table 1, Adhikari et al., 2020) --- #
    buffer *= "\t# ====================================================================================== #\n"
    buffer *= "\t# Effective biophysical model parameters (Adhikari et al., 2020)\n"
    buffer *= "\t# ====================================================================================== #\n"
    buffer *= "\tR_XT = 0.07\t# RNAP concentration (muM)\n"
    buffer *= "\tR_LT = 2.3\t# ribosome concentration (muM)\n"
    buffer *= "\tv_X = 25.0\t# transcription elongation rate (nt/s)\n"
    buffer *= "\tv_L = 1.5\t# translation elongation rate (aa/s)\n"
    buffer *= "\tK_X = 0.036\t# transcription saturation constant (muM)\n"
    buffer *= "\tK_L = 450.0\t# translation saturation constant (muM)\n"
    buffer *= "\tK_P = 10.0\t# polysome amplification factor\n"
    buffer *= "\ttemperature = 310.15\t# temperature (K) - 37C for cell-free\n"
    buffer *= "\ttau_L_half = 4.0\t# translation capacity half-life (hr)\n"
    buffer *= "\n"

    # per-gene time constants -
    list_of_genes = extract_species_of_type(list_of_species, :gene)
    list_of_connections::Array{ConnectionObject} = problem_object.list_of_connections
    number_of_genes = length(list_of_genes)

    buffer *= "\t# Per-gene transcription time constants (dimensionless) - \n"
    buffer *= "\ttau_X_array = [\n"
    for (index, gene_object) in enumerate(list_of_genes)
        gene_symbol = gene_object.species_symbol
        buffer *= "\t\t0.5\t;\t# $(index)\t$(gene_symbol)\n"
    end
    buffer *= "\t]\n"
    buffer *= "\n"

    buffer *= "\t# Per-gene translation time constants (dimensionless) - \n"
    buffer *= "\ttau_L_array = [\n"
    for (index, gene_object) in enumerate(list_of_genes)
        gene_symbol = gene_object.species_symbol
        buffer *= "\t\t0.5\t;\t# $(index)\t$(gene_symbol)\n"
    end
    buffer *= "\t]\n"
    buffer *= "\n"

    # per-gene degradation constants -
    buffer *= "\t# Per-gene mRNA degradation rate constants (hr^-1) - \n"
    buffer *= "\ttheta_m_array = [\n"
    for (index, gene_object) in enumerate(list_of_genes)
        gene_symbol = gene_object.species_symbol
        buffer *= "\t\t$(round(log(2)/(15.0/60.0); digits=4))\t;\t# $(index)\t$(gene_symbol) (half-life = 15 min)\n"
    end
    buffer *= "\t]\n"
    buffer *= "\n"

    buffer *= "\t# Per-gene protein degradation rate constants (hr^-1) - \n"
    buffer *= "\ttheta_p_array = [\n"
    for (index, gene_object) in enumerate(list_of_genes)
        gene_symbol = gene_object.species_symbol
        buffer *= "\t\t$(round(log(2)/(10.0*24.0); digits=6))\t;\t# $(index)\t$(gene_symbol) (half-life = 10 days)\n"
    end
    buffer *= "\t]\n"
    buffer *= "\n"

    # free energy dictionary -
    buffer *= "\t# Free energy dictionary for thermodynamic control (kJ/mol) - \n"
    buffer *= "\tdG_dictionary = Dict{String,Float64}()\n"
    for (index, gene_object) in enumerate(list_of_genes)

        gene_symbol = gene_object.species_symbol

        # RNAP binding energy -
        buffer *= "\tdG_dictionary[\"dG_$(gene_symbol)_RNAP\"] = -30.0\t# RNAP binding energy\n"

        # connections -
        activating_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :activate)
        inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :inhibit)

        for connection_object in activating_connections
            connection_symbol = connection_object.connection_symbol
            buffer *= "\tdG_dictionary[\"dG_$(gene_symbol)_$(connection_symbol)\"] = -20.0\t# activator (RNAP + $(connection_symbol))\n"
        end

        for connection_object in inhibiting_connections
            connection_symbol = connection_object.connection_symbol
            buffer *= "\tdG_dictionary[\"dG_$(gene_symbol)_$(connection_symbol)\"] = 10.0\t# repressor ($(connection_symbol) alone)\n"
        end
    end
    buffer *= "\n"

    # binding parameter dictionary -
    parameter_value_default_dictionary = problem_object.configuration_dictionary["default_parameter_dictionary"]
    buffer *= "\t# Binding parameter dictionary (Hill coefficients and dissociation constants) - \n"
    buffer *= "\tbinding_parameter_dictionary = Dict{String,Float64}()\n"
    for (index, gene_object) in enumerate(list_of_genes)

        gene_symbol = gene_object.species_symbol

        activating_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :activate)
        inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :inhibit)

        list_of_all_connections = ConnectionObject[]
        append!(list_of_all_connections, activating_connections)
        append!(list_of_all_connections, inhibiting_connections)
        for connection_object in list_of_all_connections
            connection_symbol = connection_object.connection_symbol
            buffer *= "\tbinding_parameter_dictionary[\"n_$(gene_symbol)_$(connection_symbol)\"] = 1.5\n"
            buffer *= "\tbinding_parameter_dictionary[\"K_$(gene_symbol)_$(connection_symbol)\"] = 50.0\t# muM\n"
        end
    end
    buffer *= "\n"

    # control parameter dictionary (needed for compatibility) -
    buffer *= "\t# Control parameter dictionary - \n"
    buffer *= "\tcontrol_parameter_dictionary = Dict{String,Float64}()\n"
    buffer *= "\n"

    # degradation modifiers -
    buffer *= "\t# degradation modifiers - \n"
    buffer *= "\tdegradation_modifier_array = [\n"
    for (index, species_object) in enumerate(list_of_species)
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type
        if (species_type == :gene)
            buffer *= "\t\t0.0\t;\t# $(index)\t$(species_symbol)\n"
        elseif (species_type == :mrna || species_type == :protein)
            buffer *= "\t\t1.0\t;\t# $(index)\t$(species_symbol)\n"
        end
    end
    buffer *= "\t]\n"
    buffer *= "\n"

    # time constant modifiers -
    buffer *= "\t# time constant modifiers - \n"
    buffer *= "\ttime_constant_modifier_array = [\n"
    for (index, species_object) in enumerate(list_of_species)
        species_symbol = species_object.species_symbol
        species_type = species_object.species_type
        if (species_type == :gene)
            buffer *= "\t\t0.0\t;\t# $(index)\t$(species_symbol)\n"
        elseif (species_type == :mrna || species_type == :protein)
            buffer *= "\t\t1.0\t;\t# $(index)\t$(species_symbol)\n"
        end
    end
    buffer *= "\t]\n"
    buffer *= "\n"

    # build the degradation matrix -
    buffer *= "\t# Dilution degrdation matrix - \n"
    buffer *= "\tdilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary,species_symbol_type_array,degradation_modifier_array)\n"

    # PTM parameters (if any) -
    ptm_types = Set([:phosphorylate, :dephosphorylate, :bind, :unbind])
    ptm_connections = filter(c -> c.connection_type in ptm_types, list_of_connections)
    if !isempty(ptm_connections)
        buffer *= "\n"
        buffer *= "\t# Post-translational modification parameters - \n"
        buffer *= "\tptm_parameter_dictionary = Dict{String,Float64}()\n"

        for conn in ptm_connections
            conn_symbol = conn.connection_symbol
            if conn.connection_type == :phosphorylate || conn.connection_type == :dephosphorylate
                buffer *= "\tptm_parameter_dictionary[\"kcat_$(conn_symbol)\"] = 10.0\t# hr^-1\n"
                buffer *= "\tptm_parameter_dictionary[\"Km_$(conn_symbol)\"] = 0.1\t# muM\n"
            elseif conn.connection_type == :bind
                buffer *= "\tptm_parameter_dictionary[\"kf_$(conn_symbol)\"] = 1.0\t# muM^-1*hr^-1\n"
            elseif conn.connection_type == :unbind
                buffer *= "\tptm_parameter_dictionary[\"kr_$(conn_symbol)\"] = 0.1\t# hr^-1\n"
            end
        end
    end

    # parameter name mapping -
    buffer *= "\n"
    buffer *= "\t# Parameter name index array - \n"
    name_parameter_mapping_buffer = generate_parameter_name_mapping_effective(list_of_genes, list_of_connections)
    buffer *= name_parameter_mapping_buffer

    # return block -
    buffer *= "\n"
    buffer *= "\t# =============================== DO NOT EDIT BELOW THIS LINE ============================== #\n"
    buffer *= "\tdata_dictionary = Dict{String,Any}()\n"
    buffer *= "\tdata_dictionary[\"number_of_states\"] = number_of_states\n"
    buffer *= "\tdata_dictionary[\"species_symbol_type_array\"] = species_symbol_type_array\n"
    buffer *= "\tdata_dictionary[\"initial_condition_array\"] = initial_condition_array\n"
    buffer *= "\tdata_dictionary[\"gene_coding_length_array\"] = gene_coding_length_array\n"
    buffer *= "\tdata_dictionary[\"mRNA_coding_length_array\"] = mRNA_coding_length_array\n"
    buffer *= "\tdata_dictionary[\"protein_coding_length_array\"] = protein_coding_length_array\n"
    buffer *= "\tdata_dictionary[\"stoichiometric_matrix\"] = stoichiometric_matrix\n"
    buffer *= "\tdata_dictionary[\"dilution_degradation_matrix\"] = dilution_degradation_matrix\n"
    buffer *= "\tdata_dictionary[\"binding_parameter_dictionary\"] = binding_parameter_dictionary\n"
    buffer *= "\tdata_dictionary[\"control_parameter_dictionary\"] = control_parameter_dictionary\n"
    buffer *= "\tdata_dictionary[\"parameter_name_mapping_array\"] = parameter_name_mapping_array\n"
    buffer *= "\tdata_dictionary[\"degradation_modifier_array\"] = degradation_modifier_array\n"
    buffer *= "\tdata_dictionary[\"time_constant_modifier_array\"] = time_constant_modifier_array\n"
    buffer *= "\tdata_dictionary[\"biophysical_constants_dictionary\"] = biophysical_constants_dictionary\n"
    buffer *= "\n"
    buffer *= "\t# Effective biophysical model parameters - \n"
    buffer *= "\tdata_dictionary[\"R_XT\"] = R_XT\n"
    buffer *= "\tdata_dictionary[\"R_LT\"] = R_LT\n"
    buffer *= "\tdata_dictionary[\"v_X\"] = v_X\n"
    buffer *= "\tdata_dictionary[\"v_L\"] = v_L\n"
    buffer *= "\tdata_dictionary[\"K_X\"] = K_X\n"
    buffer *= "\tdata_dictionary[\"K_L\"] = K_L\n"
    buffer *= "\tdata_dictionary[\"K_P\"] = K_P\n"
    buffer *= "\tdata_dictionary[\"temperature\"] = temperature\n"
    buffer *= "\tdata_dictionary[\"tau_L_half\"] = tau_L_half\n"
    buffer *= "\tdata_dictionary[\"tau_X_array\"] = tau_X_array\n"
    buffer *= "\tdata_dictionary[\"tau_L_array\"] = tau_L_array\n"
    buffer *= "\tdata_dictionary[\"theta_m_array\"] = theta_m_array\n"
    buffer *= "\tdata_dictionary[\"theta_p_array\"] = theta_p_array\n"
    buffer *= "\tdata_dictionary[\"dG_dictionary\"] = dG_dictionary\n"

    if !isempty(ptm_connections)
        buffer *= "\tdata_dictionary[\"ptm_parameter_dictionary\"] = ptm_parameter_dictionary\n"
    end

    buffer *= "\t# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #\n"
    buffer *= "\treturn data_dictionary\n"
    buffer *= "end\n"

    # build the component -
    program_component::ProgramComponent = ProgramComponent()
    program_component.filename = filename
    program_component.buffer = buffer

    # return -
    return (program_component)
end

function generate_parameter_name_mapping_effective(list_of_genes::Array{SpeciesObject}, list_of_connections::Array{ConnectionObject})

    # iterate through the list of genes - add names to a string -
    list_of_names = String[]

    for (index, gene_object) in enumerate(list_of_genes)

        gene_symbol = gene_object.species_symbol

        activating_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :activate)
        inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :inhibit)

        list_of_all_connections = ConnectionObject[]
        append!(list_of_all_connections, activating_connections)
        append!(list_of_all_connections, inhibiting_connections)
        for connection_object in list_of_all_connections
            connection_symbol = connection_object.connection_symbol
            push!(list_of_names, "n_$(gene_symbol)_$(connection_symbol)")
            push!(list_of_names, "K_$(gene_symbol)_$(connection_symbol)")
        end
    end

    # free energy parameters -
    for (index, gene_object) in enumerate(list_of_genes)
        gene_symbol = gene_object.species_symbol

        push!(list_of_names, "dG_$(gene_symbol)_RNAP")

        activating_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :activate)
        inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections, gene_object, :inhibit)

        list_of_all_connections = ConnectionObject[]
        append!(list_of_all_connections, activating_connections)
        append!(list_of_all_connections, inhibiting_connections)
        for connection_object in list_of_all_connections
            connection_symbol = connection_object.connection_symbol
            push!(list_of_names, "dG_$(gene_symbol)_$(connection_symbol)")
        end
    end

    # effective model global parameters -
    push!(list_of_names, "R_XT")
    push!(list_of_names, "R_LT")
    push!(list_of_names, "v_X")
    push!(list_of_names, "v_L")
    push!(list_of_names, "K_X")
    push!(list_of_names, "K_L")
    push!(list_of_names, "K_P")
    push!(list_of_names, "tau_L_half")
    push!(list_of_names, "temperature")

    # per-gene parameters -
    for (index, gene_object) in enumerate(list_of_genes)
        gene_symbol = gene_object.species_symbol
        push!(list_of_names, "tau_X_$(gene_symbol)")
        push!(list_of_names, "tau_L_$(gene_symbol)")
        push!(list_of_names, "theta_m_$(gene_symbol)")
        push!(list_of_names, "theta_p_$(gene_symbol)")
    end

    buffer = ""
    buffer *= "\tparameter_name_mapping_array = [\n"
    for (index, symbol) in enumerate(list_of_names)
        buffer *= "\t\t\"$(symbol)\"\t;\t# $(index)\n"
    end
    buffer *= "\t]\n"

    return buffer
end

function build_types_buffer_effective()

    filename = "Types.jl"

    buffer = ""
    buffer *= "# Types for the effective biophysical model\n"
    buffer *= "# No custom types required - all parameters are stored in data dictionary\n"

    program_component::ProgramComponent = ProgramComponent()
    program_component.filename = filename
    program_component.buffer = buffer

    return program_component
end
