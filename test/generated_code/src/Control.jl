# ----------------------------------------------------------------------------------- #
# Copyright (c) 2021 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
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
#
# ----------------------------------------------------------------------------------- #
# Function: calculate_transcription_control_array
# Description: Calculate the transcriptional control array at time t
# Generated on: 2021-05-11T16:19:52.053
#
# Input arguments:
# t::Float64 => Current time value (scalar) 
# x::Array{Float64,1} => State array (number_of_species x 1) 
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters 
#
# Output arguments:
# control_array::Array{Float64,1} => Transcriptional control array (number_of_genes x 1) at time t 
# ----------------------------------------------------------------------------------- #
function calculate_transcription_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})

	# initialize the control - 
	control_array = zeros(3)

	# Alias the species - 
	GntR = x[1]
	gene_Venus = x[2]
	sigma70 = x[3]
	mRNA_GntR = x[4]
	mRNA_gene_Venus = x[5]
	mRNA_sigma70 = x[6]
	protein_GntR = x[7]
	protein_gene_Venus = x[8]
	protein_sigma70 = x[9]

	# Alias the binding parameters - 
	binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
	n_gene_Venus_sigma70 = binding_parameter_dictionary["n_gene_Venus_sigma70"]
	K_gene_Venus_sigma70 = binding_parameter_dictionary["K_gene_Venus_sigma70"]
	n_gene_Venus_GntR = binding_parameter_dictionary["n_gene_Venus_GntR"]
	K_gene_Venus_GntR = binding_parameter_dictionary["K_gene_Venus_GntR"]

	# Alias the control function parameters - 
	control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
	W_GntR_RNAP = control_parameter_dictionary["W_GntR_RNAP"]
	W_gene_Venus_RNAP = control_parameter_dictionary["W_gene_Venus_RNAP"]
	W_gene_Venus_sigma70 = control_parameter_dictionary["W_gene_Venus_sigma70"]
	W_gene_Venus_GntR = control_parameter_dictionary["W_gene_Venus_GntR"]
	W_sigma70_RNAP = control_parameter_dictionary["W_sigma70_RNAP"]

	# Control function for GntR - 
	control_array[1] = (W_GntR_RNAP)/(1+W_GntR_RNAP)

	# Transfer function target:gene_Venus actor:sigma70
	actor_set_gene_Venus_sigma70 = [
		protein_sigma70
	]
	actor = prod(actor_set_gene_Venus_sigma70)
	b_gene_Venus_sigma70 = (actor^(n_gene_Venus_sigma70))/(K_gene_Venus_sigma70^(n_gene_Venus_sigma70)+actor^(n_gene_Venus_sigma70))

	# Transfer function target:gene_Venus actor:GntR
	actor_set_gene_Venus_GntR = [
		protein_GntR
	]
	actor = prod(actor_set_gene_Venus_GntR)
	b_gene_Venus_GntR = (actor^(n_gene_Venus_GntR))/(K_gene_Venus_GntR^(n_gene_Venus_GntR)+actor^(n_gene_Venus_GntR))

	# Control function for gene_Venus - 
	control_array[2] = (W_gene_Venus_RNAP+W_gene_Venus_sigma70*b_gene_Venus_sigma70)/(1+W_gene_Venus_RNAP+W_gene_Venus_sigma70*b_gene_Venus_sigma70+W_gene_Venus_GntR*b_gene_Venus_GntR)

	# Control function for sigma70 - 
	control_array[3] = (W_sigma70_RNAP)/(1+W_sigma70_RNAP)

	# return - 
	return control_array
end

#
# ----------------------------------------------------------------------------------- #
# Function: calculate_translation_control_array
# Description: Calculate the translation control array at time t
# Generated on: 2021-05-11T16:19:52.062
#
# Input arguments:
# t::Float64 => Current time value (scalar) 
# x::Array{Float64,1} => State array (number_of_species x 1) 
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters 
#
# Output arguments:
# control_array::Array{Float64,1} => Translation control array (number_of_genes x 1) at time t 
# ----------------------------------------------------------------------------------- #
function calculate_translation_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})

	# initialize the control - 
	control_array = ones(3)

	# return - 
	return control_array
end
