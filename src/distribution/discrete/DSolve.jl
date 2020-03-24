# include -
include("Include.jl")

function main(time_start,time_stop,time_step_size,data_dictionary)

    # Get the system dimension -
    gene_coding_length_array = data_dictionary["gene_coding_length_array"]
    number_of_genes = length(gene_coding_length_array)

    # Generate the discrete system arrays -
    (AM,BM,DM) = generate_discrete_grn_system(data_dictionary,time_step_size)

    # initilize the state/time archives -
    (number_of_species,number_of_species) = size(AM)
    state_archive = zeros(1,number_of_species)
    time_archive = Float64[]

    # what is my global IC -
    IC = data_dictionary["initial_condition_array"][(number_of_genes+1):end]
    for state_index = 1:number_of_species
        state_archive[1,state_index] = IC[state_index]
    end

    # push the initial time -
    push!(time_archive,time_start)

    # ======================================================================== #
    # Run the model at the specified time points -
    time_array = collect(time_start:time_step_size:time_stop)

    # Get the initial state array -
    xold = state_archive[end,:]

    # main loop -
    for time_value in time_array

        # calculate the new state -
        xnew = evaluate_discrete_system(time_value,xold,AM,BM,data_dictionary)


        # cache the new state and time -
        state_archive = archive_solution_array(state_archive,xnew)
        push!(time_archive,time_value)

        # update the old state -
        xold = xnew
    end
    # ======================================================================== #

    # return the archives -
    return (time_archive,state_archive)
end

# Initailize -
time_start = 0.0
time_stop = 12.0
time_step_size = 0.01

# load the data dictionary -
data_dictionary = DataDictionary(0.0,0.0,0.0)

# call main -
(time_archive,state_archive) = main(time_start,time_stop,time_step_size,data_dictionary)

# dump results to the tmp folder -
data_archive = [time_archive state_archive]
writedlm("TODO:PATH_TO_FILE",data_archive)
