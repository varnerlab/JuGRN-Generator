# Alias the txtl parameters -
txtl_parameter_dictionary = Dict{String,Float64}()
txtl_parameter_dictionary["rnapII_concentration"] = rnapII_concentration  # muM
txtl_parameter_dictionary["ribosome_concentration"] = ribosome_concentration # muM
txtl_parameter_dictionary["degradation_constant_mRNA"] = degradation_constant_mRNA  # hr^-1
txtl_parameter_dictionary["degradation_constant_protein"] = degradation_constant_protein  # hr^-1
txtl_parameter_dictionary["kcat_transcription"] = kcat_transcription  # hr^-1
txtl_parameter_dictionary["kcat_translation"] = kcat_translation  # hr^-1
txtl_parameter_dictionary["maximum_specific_growth_rate"] = maximum_specific_growth_rate  # hr^-1
txtl_parameter_dictionary["avg_gene_concentration"] = avg_gene_concentration
txtl_parameter_dictionary["saturation_constant_transcription"] = saturation_transcription
txtl_parameter_dictionary["saturation_constant_translation"] = saturation_translation
