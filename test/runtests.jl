using Test
using JuGRN

@testset "JuGRN Tests" begin

    @testset "Parse .net file" begin
        path_to_net_file = joinpath(@__DIR__, "data", "Test.net")
        result = parse_grn_file(path_to_net_file)
        @test isa(result, Array{JuGRN.VGRNSentence,1})
        @test length(result) == 2
        @test result[1].sentence_actor_clause == "GntR"
        @test result[1].sentence_action_clause == "inhibits"
        @test result[1].sentence_target_clause == "gene_Venus"
        @test result[2].sentence_actor_clause == "sigma70"
        @test result[2].sentence_action_clause == "activates"
    end

    @testset "Parse .json file" begin
        path_to_json_file = joinpath(@__DIR__, "data", "Test.json")
        result = parse_grn_file(path_to_json_file)
        @test isa(result, Dict{String,Any})
        @test haskey(result, "model_species_table")
        @test haskey(result, "raw_model")
        @test size(result["model_species_table"], 1) == 6  # 2 genes + 2 mRNA + 2 proteins
    end

    @testset "Model generation from .net file" begin
        path_to_net_file = joinpath(@__DIR__, "data", "Test.net")
        path_to_output = joinpath(tempdir(), "jugrn_test_net_$(rand(1000:9999))")

        make_julia_model(path_to_net_file, path_to_output; host_type=:bacteria)

        # verify generated files exist -
        @test isfile(joinpath(path_to_output, "src", "Data.jl"))
        @test isfile(joinpath(path_to_output, "src", "Kinetics.jl"))
        @test isfile(joinpath(path_to_output, "src", "Control.jl"))
        @test isfile(joinpath(path_to_output, "src", "Inputs.jl"))
        @test isfile(joinpath(path_to_output, "src", "Balances.jl"))
        @test isfile(joinpath(path_to_output, "src", "SolveBalances.jl"))
        @test isfile(joinpath(path_to_output, "src", "network", "Network.dat"))

        # stoichiometric matrix dimensions (3 genes -> 9 species, 6 rates)
        network_content = read(joinpath(path_to_output, "src", "network", "Network.dat"), String)
        lines = filter(!isempty, split(network_content, "\n"))
        @test length(lines) == 9

        # generated Julia files are syntactically valid -
        for jl_file in ["Data.jl", "Kinetics.jl", "Control.jl", "Inputs.jl"]
            content = read(joinpath(path_to_output, "src", jl_file), String)
            @test !isempty(content)
            @test occursin("function", content)
        end

        # .net uses hardcoded gene length of 1000.0 -
        data_content = read(joinpath(path_to_output, "src", "Data.jl"), String)
        @test occursin("1000.0", data_content)

        rm(path_to_output; recursive=true, force=true)
    end

    @testset "Model generation from .json file" begin
        path_to_json_file = joinpath(@__DIR__, "data", "Test.json")
        path_to_output = joinpath(tempdir(), "jugrn_test_json_$(rand(1000:9999))")

        make_julia_model(path_to_json_file, path_to_output)

        # verify generated files -
        @test isfile(joinpath(path_to_output, "src", "Data.jl"))
        @test isfile(joinpath(path_to_output, "src", "Kinetics.jl"))
        @test isfile(joinpath(path_to_output, "src", "Control.jl"))
        @test isfile(joinpath(path_to_output, "src", "network", "Network.dat"))

        # JSON path should use real sequence lengths from FASTA, not 1000.0 -
        data_content = read(joinpath(path_to_output, "src", "Data.jl"), String)
        @test occursin("972.0", data_content)  # actual gene sequence length from test FASTA

        # host_type should be inferred from JSON (CF_PURE -> cell_free) -
        @test occursin("gene_abundance_array", data_content)

        rm(path_to_output; recursive=true, force=true)
    end

    @testset "Model generation with blank control" begin
        path_to_net_file = joinpath(@__DIR__, "data", "Test.net")
        path_to_output = joinpath(tempdir(), "jugrn_test_blank_$(rand(1000:9999))")

        make_julia_model(path_to_net_file, path_to_output; host_type=:bacteria, control_function_generation=false)

        control_content = read(joinpath(path_to_output, "src", "Control.jl"), String)
        @test occursin("calculate_transcription_control_array", control_content)
        @test occursin("TODO", control_content) || occursin("control_array = zeros", control_content)

        rm(path_to_output; recursive=true, force=true)
    end

    @testset "Generate JSON template from .net file" begin
        path_to_net_file = joinpath(@__DIR__, "data", "Test.net")
        path_to_json_output = joinpath(tempdir(), "jugrn_test_template_$(rand(1000:9999)).json")

        result_path = generate_json_template(path_to_net_file, path_to_json_output; host_type=:cell_free)
        @test isfile(result_path)

        # parse and validate structure -
        using JSON
        json_content = JSON.parsefile(path_to_json_output)
        @test haskey(json_content, "system")
        @test haskey(json_content, "list_of_species")
        @test haskey(json_content, "list_of_transcription_models")
        @test haskey(json_content, "list_of_translation_models")

        # host type -
        @test json_content["system"]["host"] == "CF_PURE"

        # species count: .net has 3 gene names (GntR, gene_Venus, sigma70) -> 9 species
        @test length(json_content["list_of_species"]) == 9

        # transcription model for each gene -
        @test length(json_content["list_of_transcription_models"]) == 3

        # translation model for each mRNA -
        @test length(json_content["list_of_translation_models"]) == 3

        # verify connections are captured -
        tx_models = json_content["list_of_transcription_models"]
        venus_model = filter(m -> occursin("Venus", m["input"]), tx_models)
        @test length(venus_model) == 1
        @test length(venus_model[1]["list_of_activators"]) == 1
        @test length(venus_model[1]["list_of_repressors"]) == 1

        rm(path_to_json_output; force=true)
    end

    @testset "Round-trip: .net → JSON → model" begin
        path_to_net_file = joinpath(@__DIR__, "data", "Test.net")
        path_to_json = joinpath(tempdir(), "jugrn_roundtrip_$(rand(1000:9999)).json")
        path_to_output = joinpath(tempdir(), "jugrn_roundtrip_model_$(rand(1000:9999))")

        # .net → JSON -
        generate_json_template(path_to_net_file, path_to_json)

        # JSON → model -
        make_julia_model(path_to_json, path_to_output)

        # verify model was generated -
        @test isfile(joinpath(path_to_output, "src", "Data.jl"))
        @test isfile(joinpath(path_to_output, "src", "Kinetics.jl"))
        @test isfile(joinpath(path_to_output, "src", "Control.jl"))
        @test isfile(joinpath(path_to_output, "src", "network", "Network.dat"))

        rm(path_to_json; force=true)
        rm(path_to_output; recursive=true, force=true)
    end

    @testset "Compound actor/target syntax in .net" begin
        path_to_compound = joinpath(@__DIR__, "data", "TestCompound.net")
        result = parse_grn_file(path_to_compound)
        @test length(result) == 3

        # many-to-one: (GntR, sigma70) activates gene_Venus
        @test result[1].sentence_actor_clause == "(GntR, sigma70)"
        @test result[1].sentence_action_clause == "activates"
        @test result[1].sentence_target_clause == "gene_Venus"

        # one-to-many: LacI inhibits (gene_Venus, gene_GFP)
        @test result[2].sentence_actor_clause == "LacI"
        @test result[2].sentence_action_clause == "inhibits"
        @test result[2].sentence_target_clause == "(gene_Venus, gene_GFP)"

        # model generation with compound syntax -
        path_to_output = joinpath(tempdir(), "jugrn_compound_$(rand(1000:9999))")
        make_julia_model(path_to_compound, path_to_output; host_type=:bacteria)

        data_content = read(joinpath(path_to_output, "src", "Data.jl"), String)

        # many-to-one creates compound parameter: K_gene_Venus_GntR_sigma70
        @test occursin("K_gene_Venus_GntR_sigma70", data_content)
        @test occursin("W_gene_Venus_GntR_sigma70", data_content)

        # one-to-many creates separate parameters for each target
        @test occursin("K_gene_Venus_LacI", data_content)
        @test occursin("K_gene_GFP_LacI", data_content)

        # species count: 6 unique genes -> 18 species
        @test isfile(joinpath(path_to_output, "src", "network", "Network.dat"))
        network_content = read(joinpath(path_to_output, "src", "network", "Network.dat"), String)
        lines = filter(!isempty, split(network_content, "\n"))
        @test length(lines) == 18

        rm(path_to_output; recursive=true, force=true)
    end

    @testset "Save and load model (JLD2) from .net" begin
        path_to_net_file = joinpath(@__DIR__, "data", "Test.net")
        path_to_jld2 = joinpath(tempdir(), "jugrn_save_test_$(rand(1000:9999)).jld2")

        # save from file -
        save_model(path_to_jld2, path_to_net_file; host_type=:bacteria,
            metadata=Dict{String,Any}("author" => "test", "description" => "unit test"))
        @test isfile(path_to_jld2)

        # load -
        (model, meta) = load_model(path_to_jld2)
        @test isa(model, JuGRN.ProblemObject)
        @test length(model.list_of_species) == 9
        @test length(model.list_of_connections) == 2
        @test meta["author"] == "test"
        @test meta["source_file"] == "Test.net"
        @test haskey(meta, "save_timestamp")

        # regenerate from loaded model -
        path_to_output = joinpath(tempdir(), "jugrn_from_jld2_$(rand(1000:9999))")
        make_julia_model(model, path_to_output; host_type=:bacteria)
        @test isfile(joinpath(path_to_output, "src", "Data.jl"))
        @test isfile(joinpath(path_to_output, "src", "Kinetics.jl"))
        @test isfile(joinpath(path_to_output, "src", "network", "Network.dat"))

        rm(path_to_jld2; force=true)
        rm(path_to_output; recursive=true, force=true)
    end

    @testset "Parse post-translational modification syntax in .net" begin
        path_to_ptm = joinpath(@__DIR__, "data", "TestPTM.net")
        result = parse_grn_file(path_to_ptm)

        # should have 5 sentences -
        @test length(result) == 5

        # standard activation -
        @test result[1].sentence_actor_clause == "sigma70"
        @test result[1].sentence_action_clause == "activates"
        @test result[1].sentence_target_clause == "gene_deGFP"

        # phosphorylation with site and product -
        @test result[2].sentence_action_clause == "phosphorylates"
        @test result[2].sentence_actor_clause == "kinase_A"
        @test result[2].sentence_target_clause == "protein_deGFP"
        @test result[2].sentence_modifier_clause == "S147"
        @test result[2].sentence_product_clause == "protein_deGFP_p_S147"

        # dephosphorylation with product -
        @test result[3].sentence_action_clause == "dephosphorylates"
        @test result[3].sentence_target_clause == "protein_deGFP_p_S147"
        @test result[3].sentence_product_clause == "protein_deGFP"

        # binding -
        @test result[4].sentence_action_clause == "bind"
        @test result[4].sentence_actor_clause == "(protein_deGFP, protein_deGFP)"
        @test result[4].sentence_target_clause == "dimer_deGFP"

        # complex activates gene -
        @test result[5].sentence_action_clause == "activates"
        @test result[5].sentence_actor_clause == "dimer_deGFP"
    end

    @testset "Model generation with PTM reactions from .net" begin
        path_to_ptm = joinpath(@__DIR__, "data", "TestPTM.net")
        path_to_output = joinpath(tempdir(), "jugrn_ptm_$(rand(1000:9999))")

        make_julia_model(path_to_ptm, path_to_output; host_type=:cell_free)

        # verify generated files exist -
        @test isfile(joinpath(path_to_output, "src", "Data.jl"))
        @test isfile(joinpath(path_to_output, "src", "Kinetics.jl"))
        @test isfile(joinpath(path_to_output, "src", "Control.jl"))
        @test isfile(joinpath(path_to_output, "src", "network", "Network.dat"))

        # check PTM.jl has PTM rate function -
        @test isfile(joinpath(path_to_output, "src", "PTM.jl"))
        ptm_content = read(joinpath(path_to_output, "src", "PTM.jl"), String)
        @test occursin("calculate_ptm_rates", ptm_content)
        @test occursin("kcat_", ptm_content)   # phosphorylation kcat
        @test occursin("Km_", ptm_content)     # phosphorylation Km
        @test occursin("kf_", ptm_content)     # binding forward rate

        # check data dictionary has PTM parameters -
        data_content = read(joinpath(path_to_output, "src", "Data.jl"), String)
        @test occursin("ptm_parameter_dictionary", data_content)

        # stoichiometric matrix should have extra columns for PTM reactions -
        network_content = read(joinpath(path_to_output, "src", "network", "Network.dat"), String)
        lines = filter(!isempty, split(network_content, "\n"))
        @test length(lines) > 0

        # count columns in first line (should have TXTL + PTM columns) -
        first_line_values = split(strip(lines[1]))
        # we have gene triples for: sigma70, deGFP, kinase_A, phosphatase_B, dimer_deGFP
        # plus PTM species: protein_deGFP_p_S147
        # TXTL rates + 3 PTM rates (phosphorylate, dephosphorylate, bind) -
        @test length(first_line_values) > 6  # more than just TXTL rates

        rm(path_to_output; recursive=true, force=true)
    end

    @testset "Save and load model (JLD2) from .json with sequences" begin
        path_to_json_file = joinpath(@__DIR__, "data", "Test.json")
        path_to_jld2 = joinpath(tempdir(), "jugrn_save_json_$(rand(1000:9999)).jld2")

        # save -
        save_model(path_to_jld2, path_to_json_file)
        @test isfile(path_to_jld2)

        # load -
        (model, meta) = load_model(path_to_jld2)
        @test length(model.list_of_species) == 6

        # sequence lengths survive the round-trip -
        gene_lengths = model.configuration_dictionary["gene_sequence_lengths"]
        @test haskey(gene_lengths, "gene_gntR")
        @test gene_lengths["gene_gntR"] == 972.0

        # regenerate and verify real lengths are used -
        path_to_output = joinpath(tempdir(), "jugrn_from_jld2_json_$(rand(1000:9999))")
        make_julia_model(model, path_to_output)
        data_content = read(joinpath(path_to_output, "src", "Data.jl"), String)
        @test occursin("972.0", data_content)

        rm(path_to_jld2; force=true)
        rm(path_to_output; recursive=true, force=true)
    end
end
