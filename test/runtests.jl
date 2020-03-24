using JuGRN

# -- Model creation tests ------------------------------------------------------ #
function build_grn_model_test() 
    return true
end
#------------------------------------------------------------------------------- #


@testset "model_creation_test_set" begin
    @test build_grn_model_test() == true
end