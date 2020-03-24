function include_function(filename::String,pad_string)

    # create src_buffer -
    src_buffer::Array{String,1} = String[]

    # path to distrubtion -
    path_to_src_file = "$(path_to_package)/include/$(filename).jl"
    open(path_to_src_file,"r") do src_file
        for line in eachline(src_file)

            new_line_with_line_ending = line*"\n"
            push!(src_buffer,new_line_with_line_ending)
        end
    end

    string_value = ""
    for line in src_buffer
        string_value *= pad_string*line
    end

    return string_value
end
