"""
    build_julia_model_object(path_to_model_file::String, path_to_output_dir::String;
        defaults_file_name::String="Defaults.toml")->VLJuliaModelObject

Build a JUGRN data object in Julia language. If a TOML file with default parameter values is not provided by the user, a `Defaults.toml` file is generated automatically. 
Additionally, the user can edit the generated `Defaults.toml` file with their own values. After generating code, if a directory already exists at the user specified location, 
it can be deleted or backed-up before new code is written based on user input.

Input arguments:
`path_to_model_file::String` - user-specified path to read the Network.vff file from as well as the Defaults.toml file (if it already exists)
`path_to_output_dir::String` - path to where generated code will be written
`defaults_file_name::String` - name of the TOML file with default parameter values (optional).

Output arguments:
`model_object::VLResult{VLJuliaModelObject}` - abstract Julia data object holding information for generating model code.
"""
function build_julia_model_object(path_to_model_file::String, path_to_output_dir::String;
    defaults_file_name::String="Defaults.toml")::VLJuliaModelObject

    # Check: is path_to_model_file legit?
    check_result = is_file_path_ok(path_to_model_file)
    if (isnothing(check_result.value) == false)
        return check_result
    end

    # Check: do we have a trailing slash on the path_to_output_dir?
    last_char = path_to_output_dir[end]
    if (last_char == '/')
        path_to_output_dir = path_to_output_dir[1:(end - 1)]
    end

    # Check: is path_to_output_dir legit?
    # before we get too far along, we need to check if the user already has code in the location that we want to generate at -
    # if they do, then move it -
    if (isdir(path_to_output_dir) == true)

        # let the user know that we found code -
        @info "Ooops! We found some code where you wanted to generate your model code."

        # ask the user: should we move or delete?
        user_input_result = _request_user_input("Overwrite the existing project files [Y]/N ? ")
        if (lowercase(user_input_result) == "no" || lowercase(user_input_result) == "n")

            @info "Got it. Backing up the existing files ..."

            # ok, looks like we may have a conflict - mv the offending code
            if (move_existing_project_at_path(path_to_output_dir) == false)

                # Something happend ... the world is ending ...
                error = ArgumentError("automatic directory conflict resolution failed. Unable to move existing directory $(path_to_output_dir)")
                
                # throw -
                throw(error)
            end
        else
            @info "Ok! Overwriting the existing files ..."
        end
    end

    # Check: do we have a defaults file?
    path_to_defaults_file = joinpath(dirname(path_to_model_file), defaults_file_name)
    if (isfile(path_to_defaults_file) == false)

        # ok, we do not have a defaults file, so lets build one, and let the user know -
        generate_defaults_result = generate_default_project_file(path_to_defaults_file)
        if (isnothing(generate_defaults_result.value) == false)
            return generate_defaults_result
        end

        # let the user know that we have generated a new Defaults.toml file -
        @info "We generated Defaults.toml using the generate_default_project_file function, see the help system in the REPL"
    end

    # build the model wrapper -
    model_object = VLJuliaModelObject(path_to_model_file, path_to_output_dir;
        defaults_file_path=path_to_defaults_file)

    # return -
    return model_object
end