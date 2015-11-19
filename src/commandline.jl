# Copyright (c) 2015 Michael Eastwood
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

const s = ArgParseSettings(description = "A calibration routine developed for the OVRO LWA.")

@add_arg_table s begin
    "gaincal"
        help = "Solve for a gain calibration."
        action = :command
    "polcal"
        help = "Solve for a polarization calibration."
        action = :command
    "peel"
        help = "Peel sources from the dataset."
        action = :command
    "applycal"
        help = "Apply a calibration."
        action = :command
end

@add_arg_table s["gaincal"] begin
    "--input"
        help = "The measurement set used to solve for the calibration."
        arg_type = ASCIIString
        required = true
    "--output"
        help = "The output file to which the calibration will be written. This output file will be overwritten if it already exists."
        arg_type = ASCIIString
        required = true
    "--sources"
        help = "A JSON file describing the sources to be used for the sky model."
        arg_type = ASCIIString
        required = true
    "--beam"
        help = "The name of the beam model to use ($(join(keys(beam_dictionary),",")))."
        arg_type = ASCIIString
        default = "sine"
    "--maxiter"
        help = "Set the maximum number of (Mitch|Stef)Cal iterations to take on each frequency channel."
        arg_type = Int
        default  = 20
    "--tolerance"
        help = "Set the relative tolerance used to determine convergence."
        arg_type = Float64
        default  = 1e-3
    "--force-imaging"
        help = "Create and use the MODEL_DATA and CORRECTED_DATA columns even if they do not already exist in the measurement set."
        action = :store_true
end

@add_arg_table s["polcal"] begin
    "--input"
        help = "The measurement set used to solve for the calibration."
        arg_type = ASCIIString
        required = true
    "--output"
        help = "The output file to which the calibration will be written. This output file will be overwritten if it already exists."
        arg_type = ASCIIString
        required = true
    "--sources"
        help = "A JSON file describing the sources to be used for the sky model."
        arg_type = ASCIIString
        required = true
    "--beam"
        help = "The name of the beam model to use ($(join(keys(beam_dictionary),",")))."
        arg_type = ASCIIString
        default = "sine"
    "--maxiter"
        help = "Set the maximum number of (Mitch|Stef)Cal iterations to take on each frequency channel."
        arg_type = Int
        default  = 20
    "--tolerance"
        help = "Set the relative tolerance used to determine convergence."
        arg_type = Float64
        default  = 1e-3
    "--force-imaging"
        help = "Create and use the MODEL_DATA and CORRECTED_DATA columns even if they do not already exist in the measurement set."
        action = :store_true
end

@add_arg_table s["peel"] begin
    "--input"
        help = "The measurement set that will have sources peeled."
        arg_type = ASCIIString
        required = true
    "--sources"
        help = "A JSON file describing the sources to be peeled from the given measurement set."
        arg_type = ASCIIString
        required = true
    "--output"
        help = "If provided, write the direction dependent calibration solutions to disk as NumPy arrays."
        arg_type = ASCIIString
        default = ""
    "--beam"
        help = "The name of the beam model to use ($(join(keys(beam_dictionary),",")))."
        arg_type = ASCIIString
        default = "sine"
    "--peeliter"
        help = "The number of iterations to take while peeling."
        arg_type = Int
        default = 3
    "--maxiter"
        help = "Set the maximum number of (Mitch|Stef)Cal iterations to take on each frequency channel."
        arg_type = Int
        default  = 20
    "--tolerance"
        help = "Set the relative tolerance used to determine convergence."
        arg_type = Float64
        default  = 1e-3
    "--minuvw"
        help = "The minimum baseline length (measured in wavelengths) to use while peeling sources. This parameter can be used to mitigate sensitivity to unmodeled diffuse emission."
        arg_type = Float64
        default  = 15.0
end

@add_arg_table s["applycal"] begin
    "--input"
        help = "The measurement set that the calibration will be applied to."
        arg_type = ASCIIString
        required = true
    "--calibration"
        help = "The file containing the calibration solution."
        arg_type = ASCIIString
        required = true
    "--force-imaging"
        help = "Write the calibrated visibilities to the CORRECTED_DATA column regardless of whether or not the measurement set already has this column (it will be created if it doesn't exist)."
        action = :store_true
    "--corrected"
        help = "Apply the calibration to the CORRECTED_DATA column instead of the DATA column."
        action = :store_true
end

function main(args)
    parsed_args = parse_args(args,s)
    command = parsed_args["%COMMAND%"]
    if     command == "gaincal"
        run_gaincal(parsed_args["gaincal"])
    elseif command == "polcal"
        run_polcal(parsed_args["polcal"])
    elseif command == "applycal"
        run_applycal(parsed_args["applycal"])
    elseif command == "peel"
        run_peel(parsed_args["peel"])
    end
end

function run_gaincal(args)
    println("Running `gaincal` on $(args["input"])")
    ms = MeasurementSet(ascii(args["input"]))
    sources = readsources(args["sources"])
    beam = beam_dictionary[args["beam"]]()
    cal = gaincal(ms,sources,beam,
                  maxiter=args["maxiter"],
                  tolerance=args["tolerance"],
                  force_imaging_columns=args["force-imaging"])
    write(args["output"],cal)
    cal
end

function run_polcal(args)
    println("Running `polcal` on $(args["input"])")
    ms = MeasurementSet(ascii(args["input"]))
    sources = readsources(args["sources"])
    beam = beam_dictionary[args["beam"]]()
    cal = polcal(ms,sources,beam,
                 maxiter=args["maxiter"],
                 tolerance=args["tolerance"],
                 force_imaging_columns=args["force-imaging"])
    write(args["output"],cal)
    cal
end

function run_peel(args)
    println("Running `peel` on $(args["input"])")
    ms = MeasurementSet(ascii(args["input"]))
    sources = readsources(args["sources"])
    beam = beam_dictionary[args["beam"]]()
    calibrations = peel!(GainCalibration,ms,sources,beam,
                         peeliter=args["peeliter"],
                         maxiter=args["maxiter"],
                         tolerance=args["tolerance"],
                         minuvw=args["minuvw"])
    if !isempty(args["output"])
        for i = 1:length(calibrations)
            filename = args["output"]*"-$i.npz"
            write_for_python(filename, calibrations[i])
        end
    end
end

function run_applycal(args)
    println("Running `applycal` on $(args["input"])")
    ms = MeasurementSet(ascii(args["input"]))
    cal = read(args["calibration"])
    force_imaging_columns = haskey(args,"force-imaging")
    apply_to_corrected = haskey(args,"corrected")
    applycal!(ms,cal,
              force_imaging_columns=args["force-imaging"],
              apply_to_corrected=args["corrected"])
    cal
end

precompile(main,(Vector{UTF8String},))
precompile(run_gaincal,(Dict{AbstractString,Any},))
precompile(run_polcal,(Dict{AbstractString,Any},))
precompile(run_peel,(Dict{AbstractString,Any},))
precompile(run_applycal,(Dict{AbstractString,Any},))

