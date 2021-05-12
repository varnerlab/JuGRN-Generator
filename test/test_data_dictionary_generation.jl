# load the generator package -
using JuGRN

# init -
test_set = Set{ProgramComponent}()

# path to test file -
path_to_net_file = "./test/data/Test.json"
model_dictionary = parse_grn_file(path_to_net_file)

# build the dd -
dd = build_data_dictionary_buffer(model_dictionary)
push!(test_set, dd)

# write -
write_program_components_to_disk("./test/generated_model",test_set)
