# load the generator package -
using JuGRN

# path to test file -
path_to_net_file = "./test/data/test.net"
path_to_output_dir = "./test/generated_code"

# go -
make_julia_model(path_to_net_file, path_to_output_dir; host_type=:cell_free, control_function_generation=true)

