# ----------------------------------------------------------------------------------- #
# Copyright (c) 2019 Varnerlab
# Robert Frederick School of Chemical and Biomolecular Engineering
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
function generate_discrete_metabolic_system(data_dictionary,time_step_size)

    # Get stuff from data dictionary -
    mugmax = data_dictionary["maximum_specific_growth_rate"]
    S = data_dictionary["stoichiometric_matrix"]

    # what is the size of the intracellular state array?
    number_of_intracellular_states = 0      # Need to fig this out at gen time

    # -- Compute AHAT ----------------------------------------------------------- %
    # Dilution and degrdation terms (we need for A matrix)
    delta_term = (mugmax)
    AM = -1*eye(number_of_intracellular_states,number_of_intracellular_states)

    # mRNA -
    for index = 1:number_of_intracellular_states
        AM[index,index] = delta_term*AM[index,index]
    end

    # Compute AHAT -
    AHAT = exp(AM*time_step_size)
    # --------------------------------------------------------------------------- %

    # -- Compute BHAT ----------------------------------------------------------- %
    # Dilution and degrdation terms (we need for A matrix)
    SM = S[1:number_of_intracellular_states,:]
    IM = eye(number_of_intracellular_states,number_of_intracellular_states)
    BHAT = (AHAT - IM)*inv(AM)*SM
    # --------------------------------------------------------------------------- %

    # -- Compute DHAT ----------------------------------------------------------- %
    DHAT = inv((IM - AHAT))
    # --------------------------------------------------------------------------- %

    # return my discrete system arrays -
    return (AHAT,BHAT,DHAT)
end

function generate_discrete_grn_system(data_dictionary,time_step_size)

    # How many genes do we have?
    gene_coding_length_array = data_dictionary["gene_coding_length_array"]
    number_of_genes = length(gene_coding_length_array)

    # Get stuff from data dictionary -
    kDM = data_dictionary["degradation_constant_mRNA"]
    kDP = data_dictionary["degradation_constant_protein"]
    mugmax = data_dictionary["maximum_specific_growth_rate"]
    S = data_dictionary["stoichiometric_matrix"]


    # -- Compute AHAT ----------------------------------------------------------- %
    # Dilution and degrdation terms (we need for A matrix)
    mRNA_delta = (kDM+mugmax)
    P_delta = (kDP+mugmax)
    AM = -1*eye(2*number_of_genes,2*number_of_genes)

    # mRNA -
    for index = 1:number_of_genes
        AM[index,index] = mRNA_delta*AM[index,index]
    end

    # protein -
    for index = (number_of_genes+1):(2*number_of_genes)
        AM[index,index] = P_delta*AM[index,index]
    end

    # Compute AHAT -
    AHAT = exp(AM*time_step_size)
    # --------------------------------------------------------------------------- %

    # -- Compute BHAT ----------------------------------------------------------- %
    # Dilution and degrdation terms (we need for A matrix)
    SM = S[(number_of_genes+1):end,:]
    IM = eye(2*number_of_genes,2*number_of_genes)
    BHAT = (AHAT - IM)*inv(AM)*SM
    # --------------------------------------------------------------------------- %

    # -- Compute DHAT ----------------------------------------------------------- %
    DHAT = inv((IM - AHAT))
    # --------------------------------------------------------------------------- %

    return (AHAT,BHAT,DHAT)
end

function evaluate_discrete_system(t,state_array,AM,BM,data_dictionary)

    # we need to pre-prepend the genes -
    # How many genes do we have?
    gene_coding_length_array = data_dictionary["gene_coding_length_array"]
    number_of_genes = length(gene_coding_length_array)

    # For now, assume all genes are at avg_gene_concentration -
    avg_gene_concentration = data_dictionary["avg_gene_concentration"]
    gene_array = avg_gene_concentration*ones(number_of_genes)
    full_state_array = [gene_array ; state_array]

    # ok, evaluate the kinetics at this state -
    rT = calculate_transcription_rates(t,full_state_array,data_dictionary)
    rX = calculate_translation_rates(t,full_state_array,data_dictionary)

    # Call my control function -
    control_array = Control(t,full_state_array,data_dictionary)

    # build the overall rate vector -
    rV = [rT.*control_array ; rX]

    # calculate the new state -
    xnew = AM*state_array+BM*rV

    # return -
    return xnew
end

function archive_solution_array(state_archive,state_array)

    # add an addional row to the archive w/the new solution -
    (nrows,ncols) = size(state_archive)

    # initialize new archive -
    new_archive_array = zeros(nrows+1,ncols)

    # copy old data into new archive -
    for row_index = 1:nrows
        for col_index = 1:ncols
            new_archive_array[row_index,col_index] = state_archive[row_index,col_index]
        end
    end

    # add new row -
    for col_index = 1:ncols
        new_archive_array[end,col_index] = state_array[col_index]
    end

    # return -
    return new_archive_array
end
