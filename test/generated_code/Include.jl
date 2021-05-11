# ----------------------------------------------------------------------------------- #
# Copyright (c) 2019 Varnerlab
# Robert Frederick School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
function _include_my_codes(path_to_src::String, fileNameArray::Array{String,1})

    try

        for file_name in fileNameArray
            
            path_to_load = joinpath(path_to_src, file_name)
            include(path_to_load)
        end

    catch error
        rethrow(error)
    end

end

# get files from dir -
searchdir(path,key) = filter(x -> contains(x, key), readdir(path))

# what is the root directory?
_PATH_TO_ROOT = pwd()
_PATH_TO_SRC = joinpath(_PATH_TO_ROOT, "src")
_PATH_TO_DATABASE = joinpath(_PATH_TO_SRC, "database")
_PATH_TO_NETWORK = joinpath(_PATH_TO_SRC, "network")

# system packages - these are required to be installed
# if not then install them
import Pkg
Pkg.activate(_PATH_TO_ROOT)
Pkg.instantiate()
Pkg.update()

using LinearAlgebra # pre-installed w/Julia
using Statistics    # pre-installed w/Julia
using DifferentialEquations
using DelimitedFiles
using JSON
using ProgressMeter
using SQLite

# includes -
my_src_codes_array = searchdir(_PATH_TO_SRC, ".jl")
_include_my_codes(_PATH_TO_SRC, my_src_codes_array)
