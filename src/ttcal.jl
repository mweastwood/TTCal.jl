#!/usr/bin/env julia

println("""
     TTCal
    =======
    A calibration routine developed for the OVRO LWA.
    Written by Michael Eastwood (mweastwood@astro.caltech.edu).
    """)

import TTCal

immutable Command
    name::UTF8String
    help::UTF8String
end

immutable Option
    flag::UTF8String
    help::UTF8String
    T::Type
    required::Bool
    min_number::Number # The minimum number of arguments passed to this option
    max_number::Number # The maximum number of arguments passed to this option
end

function Base.print(io::IO,command::Command)
    print(io,"  ")
    print(io,rpad(command.name,13," "))
    print(io,replace(command.help,"\n","\n"*" "^15))
end

function Base.print(io::IO,option::Option)
    print(io,"  ")
    print(io,rpad(option.flag,13," "))
    print(io,replace(option.help,"\n","\n"*" "^15))
    option.required && print(io," (required)")
end

const commands = [Command("flagdata","Flag antennas within a measurement set."),
                  Command("bandpass","Solve for a bandpass calibration."),
                  Command("applycal","Apply a calibration.")]
const command_names = [command.name for command in commands]

const options = Dict{UTF8String,Vector{Option}}()
options["flagdata"] = [
    Option("--input","""
        The list of measurement sets to flag.""",
        UTF8String,true,1,Inf)
    Option("--antennas","""
        The list of antennas that will be flagged within the supllied
        list of measurement sets. The first antenna is number one.""",
        Int,true,1,Inf)]
options["bandpass"] = [
    Option("--input","""
        The measurement set used to solve for a bandpass calibration.""",
        UTF8String,true,1,1),
    Option("--output","""
        The output file to which the bandpass calibration will be written.
        This output file will be overwritten if it already exists.""",
        UTF8String,true,1,1),
    Option("--sources","""
        A JSON file describing the sources to be used for the sky model.""",
        UTF8String,true,1,1),
    Option("--maxiter","""
        Set the maximum number of Stefcal iterations to take on each
        frequency channel.""",
        Int,false,1,1),
    Option("--tolerance","""
        Set the relative tolerance used to determine convergence.""",
        Float64,false,1,1)]
options["applycal"] = [
    Option("--input","""
        The list of measurement sets that the calibration will be applied to.""",
        UTF8String,true,1,Inf),
    Option("--calibration","""
        The file containing the calibration solution.""",
        UTF8String,true,1,1)]

function print_command_help()
    print("""
        usage: ttcal.jl command options...

        commands:
        """)
    for command in commands
        println(command)
    end
    println("")
end

function print_option_help(command)
    command_options = options[command]
    println("$command options:")
    for option in command_options
        println(option)
    end
    println("")
end

function parse_option(opt::Option,args)
    if length(args) > opt.max_number
        error("Too many arguments passed to flag $(opt.flag) ($(opt.max_number) maximum).")
    end
    if length(args) < opt.min_number
        error("Too few arguments passed to flag $(opt.flag) ($(opt.min_number) minimum).")
    end
    if opt.max_number > 1
        return opt.T[parse_option(opt.T,arg) for arg in args]
    end
    parse_option(opt.T,args[1])
end

function parse_option(T::Type,arg)
    try
        return parse_option_helper(T,arg)
    catch
        error("Unable to convert $arg to type $T.")
    end
end

parse_option_helper(::Type{UTF8String},arg) = utf8(arg)
parse_option_helper(::Type{Int},arg) = int(arg)
parse_option_helper(::Type{Float64},arg) = float64(arg)

looks_like_flag(str) = startswith(str,"--")

function parse_args()
    if length(ARGS) < 1 || !(ARGS[1] in command_names)
        print_command_help()
        error("Please provide one of the listed commands.")
    end

    command = ARGS[1]
    option_flags = [option.flag for option in options[command]]

    if length(ARGS) < 2
        print_option_help(command)
        error("Please select from the list options.")
    end

    args = Dict{UTF8String,Any}()
    idx = 2
    while idx <= length(ARGS)
        if !(ARGS[idx] in option_flags)
            print_option_help(command)
            error("\"$(ARGS[idx])\" is not a recognized flag.")
        end

        option_flag = ARGS[idx]
        option = options[command][findfirst(option_flags,option_flag)]

        next_idx = findnext(looks_like_flag,ARGS,idx+1)
        if next_idx == 0 # This is the last option
            args[option_flag] = parse_option(option,ARGS[idx+1:end])
            break
        else
            args[option_flag] = parse_option(option,ARGS[idx+1:next_idx-1])
            idx = next_idx
        end
    end

    # Verify that the required options are provided
    for option in options[command]
        if option.required && !(option.flag in keys(args))
            error("Required flag $(option.flag) is missing.")
        end
    end

    command,args
end

# Catch exceptions to hide the verbose Julia output that is
# not especially helpful to users.
try
    command,args = parse_args()
    if     command == "flagdata"
        TTCal.run_flagdata(args)
    elseif command == "bandpass"
        TTCal.run_bandpass(args)
    elseif command == "applycal"
        TTCal.run_applycal(args)
    end
catch err
    if isa(err, ErrorException)
        println(err.msg)
    else
        throw(err)
    end
end

