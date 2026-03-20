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
end
