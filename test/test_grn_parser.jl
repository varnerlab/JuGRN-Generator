# load the generator package -
using JuGRN

# path to test file -
path_to_net_file = "./test/data/Test.json"

# parse -
model_dictionary = parse_grn_file(path_to_net_file)
