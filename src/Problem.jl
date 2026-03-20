function generate_problem_object(intermediate_representation::Dict{String,Any})

    # Initialize -
    problem_object::ProblemObject = ProblemObject()

    # Load the base configuration -
    path_to_configuration_file = "$(path_to_package)/configuration/Configuration.json"
    config_dict = Dict{String,Any}(JSON.parsefile(path_to_configuration_file))

    # Get the raw model from the JSON -
    raw_model = intermediate_representation["raw_model"]
    species_table = intermediate_representation["model_species_table"]

    # Store system parameters and host type in config -
    if haskey(raw_model, "system")
        system_dict = raw_model["system"]
        config_dict["system_parameters"] = Dict{String,Any}(system_dict["parameters"])
        config_dict["host_type"] = string(system_dict["host"])
    end

    # Build species list from JSON species table (already typed) -
    list_of_species = SpeciesObject[]
    for row in eachrow(species_table)
        species_type = row.type == :mRNA ? :mrna : row.type  # normalize mRNA -> mrna
        push!(list_of_species, SpeciesObject(species_type, string(row.symbol)))
    end
    partition!(list_of_species)

    # Compute sequence lengths from FASTA records -
    gene_sequence_lengths = Dict{String,Float64}()
    protein_sequence_lengths = Dict{String,Float64}()
    for row in eachrow(species_table)
        if !ismissing(row.sequence) && isa(row.sequence, FASTA.Record)
            seq_length = Float64(length(FASTA.sequence(row.sequence)))
            if row.type == :gene
                gene_sequence_lengths[string(row.symbol)] = seq_length
            elseif row.type == :protein
                protein_sequence_lengths[string(row.symbol)] = seq_length
            end
        end
    end
    config_dict["gene_sequence_lengths"] = gene_sequence_lengths
    config_dict["protein_sequence_lengths"] = protein_sequence_lengths

    # Build connections from transcription models -
    list_of_connections = ConnectionObject[]
    if haskey(raw_model, "list_of_transcription_models")
        for tx_model in raw_model["list_of_transcription_models"]

            # target gene is the input of this transcription model -
            target_gene_symbol = string(tx_model["input"])
            target_species = SpeciesObject(:gene, target_gene_symbol)

            # process activators -
            if haskey(tx_model, "list_of_activators")
                for activator in tx_model["list_of_activators"]
                    connection = ConnectionObject()
                    actor_symbol = string(activator["symbol"])
                    connection.connection_symbol = actor_symbol
                    connection.connection_actor_set = [SpeciesObject(:protein, actor_symbol)]
                    connection.connection_target_set = [target_species]
                    connection.connection_type = :activate
                    push!(list_of_connections, connection)
                end
            end

            # process repressors -
            if haskey(tx_model, "list_of_repressors")
                for repressor in tx_model["list_of_repressors"]
                    connection = ConnectionObject()
                    actor_symbol = string(repressor["symbol"])
                    connection.connection_symbol = actor_symbol
                    connection.connection_actor_set = [SpeciesObject(:protein, actor_symbol)]
                    connection.connection_target_set = [target_species]
                    connection.connection_type = :inhibit
                    push!(list_of_connections, connection)
                end
            end
        end
    end

    # Store transcription and translation models for richer code generation -
    if haskey(raw_model, "list_of_transcription_models")
        config_dict["list_of_transcription_models"] = raw_model["list_of_transcription_models"]
    end
    if haskey(raw_model, "list_of_translation_models")
        config_dict["list_of_translation_models"] = raw_model["list_of_translation_models"]
    end

    # Set on problem_object -
    problem_object.configuration_dictionary = config_dict
    problem_object.list_of_species = list_of_species
    problem_object.list_of_connections = list_of_connections

    return problem_object
end

function generate_problem_object(statement_vector::Array{VGRNSentence})

    # Initilize an empty problem object -
    problem_object::ProblemObject = ProblemObject()

    # Load the JSON configuration file -
    path_to_configuration_file = "$(path_to_package)/configuration/Configuration.json"
    config_dict = JSON.parsefile(path_to_configuration_file)

    # construct the array of species -
    species_array::Array{SpeciesObject} = build_species_list(statement_vector)

    # construct the array of reactions -
    connection_array::Array{ConnectionObject} = build_connection_list(statement_vector,config_dict)

    # build PTM product map from sentence product clauses -
    ptm_product_map = Dict{String,Vector{String}}()
    ptm_verbs = Set(["phosphorylate", "phosphorylates", "phosphorylated",
                      "dephosphorylate", "dephosphorylates", "dephosphorylated",
                      "bind", "binds", "bound",
                      "unbind", "unbinds", "dissociate", "dissociates"])

    for (i, sentence) in enumerate(statement_vector)
        action = lowercase(sentence.sentence_action_clause)
        if action in ptm_verbs && isdefined(sentence, :sentence_product_clause) && !isempty(sentence.sentence_product_clause)
            # find the matching connection -
            if i <= length(connection_array)
                conn = connection_array[i]
                conn_key = conn.connection_symbol * "_" * string(conn.connection_type)
                product_symbols = String[]
                recursive_species_parser!(reverse(collect(sentence.sentence_product_clause)), product_symbols)
                ptm_product_map[conn_key] = product_symbols
            end
        end
    end
    config_dict["ptm_product_map"] = ptm_product_map

    # set data on problem_object -
    problem_object.configuration_dictionary = config_dict
    problem_object.list_of_species = species_array
    problem_object.list_of_connections = connection_array

    # return the problem_object -
    return problem_object
end

function build_connection_list(statement_vector::Array{VGRNSentence,1}, configuration_dictionary::AbstractDict{String,Any})

  # initialize -
  list_of_connections = ConnectionObject[]

  # build synonym sets from configuration -
  induction_synonyms = _build_synonym_set(configuration_dictionary, "list_of_induction_synonyms")
  repression_synonyms = _build_synonym_set(configuration_dictionary, "list_of_repression_synonyms")
  phosphorylation_synonyms = _build_synonym_set(configuration_dictionary, "list_of_phosphorylation_synonyms")
  dephosphorylation_synonyms = _build_synonym_set(configuration_dictionary, "list_of_dephosphorylation_synonyms")
  binding_synonyms = _build_synonym_set(configuration_dictionary, "list_of_binding_synonyms")
  unbinding_synonyms = _build_synonym_set(configuration_dictionary, "list_of_unbinding_synonyms")

  # iterate through the statement vector -
  for vgrn_sentence in statement_vector

    # Create conenction object -
    connection_object = ConnectionObject()

    # Who are my actors?
    list_of_actor_symbols = String[]
    actor_string = vgrn_sentence.sentence_actor_clause
    recursive_species_parser!(reverse(collect(actor_string)),list_of_actor_symbols)
    actor_set::Array{SpeciesObject} = species_object_factory(list_of_actor_symbols)

    # Who are my targets?
    list_of_target_symbols = String[]
    target_string = vgrn_sentence.sentence_target_clause
    recursive_species_parser!(reverse(collect(target_string)),list_of_target_symbols)
    target_set::Array{SpeciesObject} = species_object_factory(list_of_target_symbols)

    # set them on the connection object -
    # ok - do we have a compound set of actors -
    if (length(actor_set) == 1)

      actor_object = actor_set[1]
      connection_object.connection_symbol = actor_object.species_symbol
    else

      local_buffer = ""
      number_of_actors = length(actor_set)
      for (index,actor_object::SpeciesObject) in enumerate(actor_set)

        local_buffer *= "$(actor_object.species_symbol)"
        if (index <= number_of_actors - 1)
          local_buffer *= "_"
        end
      end

      connection_object.connection_symbol = local_buffer
    end


    connection_object.connection_actor_set = actor_set
    connection_object.connection_target_set = target_set

    # What is my type?
    connection_action_string = vgrn_sentence.sentence_action_clause
    if (in(connection_action_string, induction_synonyms))
      connection_object.connection_type = :activate
    elseif (in(connection_action_string, repression_synonyms))
      connection_object.connection_type = :inhibit
    elseif (in(connection_action_string, phosphorylation_synonyms))
      connection_object.connection_type = :phosphorylate
    elseif (in(connection_action_string, dephosphorylation_synonyms))
      connection_object.connection_type = :dephosphorylate
    elseif (in(connection_action_string, binding_synonyms))
      connection_object.connection_type = :bind
    elseif (in(connection_action_string, unbinding_synonyms))
      connection_object.connection_type = :unbind
    end

    # cache this connection -
    push!(list_of_connections,connection_object)
  end

  # return -
  return list_of_connections
end

function _build_synonym_set(config::AbstractDict{String,Any}, key::String)::Set{String}
    synonyms = Set{String}()
    if haskey(config, key)
        for entry in config[key]
            push!(synonyms, entry["symbol"])
        end
    end
    return synonyms
end

function species_object_factory(list_of_symbols::Array{String})

  species_set = SpeciesObject[]
  for symbol_text in list_of_symbols
    push!(species_set,SpeciesObject(:gene,symbol_text))
  end

  return species_set
end

function build_species_list(statement_vector::Array{VGRNSentence})

  # separate transcriptional sentences from post-translational sentences -
  ptm_verbs = Set(["phosphorylate", "phosphorylates", "phosphorylated",
                    "dephosphorylate", "dephosphorylates", "dephosphorylated",
                    "bind", "binds", "bound",
                    "unbind", "unbinds", "dissociate", "dissociates"])

  # collect gene symbols from transcriptional regulation sentences -
  list_of_gene_symbols = String[]
  ptm_species_set = Set{String}()  # additional species from PTM events

  for vgrn_sentence in statement_vector

    action = lowercase(vgrn_sentence.sentence_action_clause)

    if action in ptm_verbs
      # PTM sentence: collect actor and target symbols as protein-level species -
      actor_symbols = String[]
      recursive_species_parser!(reverse(collect(vgrn_sentence.sentence_actor_clause)), actor_symbols)
      for s in actor_symbols
        push!(ptm_species_set, s)
      end

      target_symbols = String[]
      recursive_species_parser!(reverse(collect(vgrn_sentence.sentence_target_clause)), target_symbols)
      for s in target_symbols
        push!(ptm_species_set, s)
      end

      # product clause creates new species -
      if isdefined(vgrn_sentence, :sentence_product_clause) && !isempty(vgrn_sentence.sentence_product_clause)
        product_symbols = String[]
        recursive_species_parser!(reverse(collect(vgrn_sentence.sentence_product_clause)), product_symbols)
        for s in product_symbols
          push!(ptm_species_set, s)
        end
      end
    else
      # transcriptional regulation: actors and targets are gene-level symbols -
      species_string = vgrn_sentence.sentence_actor_clause * " " * vgrn_sentence.sentence_target_clause
      recursive_species_parser!(reverse(collect(species_string)), list_of_gene_symbols)
    end
  end

  # unique gene symbols -
  gene_symbol_set = Set{String}(list_of_gene_symbols)
  unique_gene_symbols = sort!(collect(gene_symbol_set))

  # build gene/mRNA/protein triples from gene symbols -
  list_of_species_objects = SpeciesObject[]
  all_protein_symbols = Set{String}()

  for species_symbol in unique_gene_symbols

    gene_object = SpeciesObject(:gene, species_symbol)
    mRNA_symbol = "mRNA_" * species_symbol
    mRNA_object = SpeciesObject(:mrna, mRNA_symbol)
    protein_symbol = "protein_" * species_symbol
    protein_object = SpeciesObject(:protein, protein_symbol)
    push!(all_protein_symbols, protein_symbol)

    push!(list_of_species_objects, gene_object)
    push!(list_of_species_objects, mRNA_object)
    push!(list_of_species_objects, protein_object)
  end

  # add PTM species that aren't already covered by the gene/mRNA/protein triples -
  for ptm_symbol in sort!(collect(ptm_species_set))
    # check if this is already a protein from the triple generation -
    if !("protein_" * ptm_symbol in all_protein_symbols) && !(ptm_symbol in all_protein_symbols)
      # determine type: binding products are complexes, modified proteins are proteins -
      push!(list_of_species_objects, SpeciesObject(:protein, ptm_symbol))
    end
  end

  # partition by type -
  return partition!(list_of_species_objects)
end

function recursive_species_parser!(sentence_char_array::Array{Char},list_of_symbols::Array{String})

  if (isempty(sentence_char_array) == true)
    return
  end

  # process each char -
  test_char = pop!(sentence_char_array)
  if (test_char == '(' || test_char == ' ')
    recursive_species_parser!(sentence_char_array,list_of_symbols)

  elseif (isletter(test_char) == true)

    # cache this letter -
    biological_symbol_cache = Char[]
    push!(biological_symbol_cache,test_char)

    # When should we stop?
    stop_set = Set{Char}()
    push!(stop_set,' ')
    push!(stop_set,',')
    push!(stop_set,')')
    push!(stop_set,'\n')
    push!(stop_set,'\r')

    # ok, so lets read until we hit a space or a comma -
    if (isempty(sentence_char_array))
      stop_flag = true
    else
      stop_flag = false
    end

    while (stop_flag == false)

      if (isempty(sentence_char_array) == false)
        next_test_char = pop!(sentence_char_array)

        if (in(next_test_char,stop_set) == false)
          push!(biological_symbol_cache,next_test_char)
        else
          stop_flag = true
        end
      else
        stop_flag = true
      end
    end

    # Store the symbol -
    if (isempty(biological_symbol_cache) == false)
      push!(list_of_symbols,String(biological_symbol_cache))
    end

    # ok, so we should be ready for another dive -
    recursive_species_parser!(sentence_char_array,list_of_symbols)
  end
end
