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

using DocOpt

doc = """
TTCal is a calibration routine developed by Michael Eastwood for the Long Wavelength Array at the
Owens Valley Radio Observatory.

Online documentation is available at http://mweastwood.info/TTCal.jl/
Please file bug reports at https://github.com/mweastwood/TTCal.jl/issues

Usage:
    ttcal.jl applycal <ms> <cal>
        [--force-imaging] [--corrected]
    ttcal.jl gaincal <ms> <cal> <sources> [--beam BEAM]
        [--maxiter ITERS] [--tolerance TOL] [--force-imaging] [--minuvw MINUVW]
    ttcal.jl polcal <ms> <cal> <sources> [--beam BEAM]
        [--maxiter ITERS] [--tolerance TOL] [--force-imaging] [--minuvw MINUVW]
    ttcal.jl peel <ms> <sources> [--beam BEAM]
        [--peeliter PITERS] [--maxiter ITERS] [--tolerance TOL] [--minuvw MINUVW]
    ttcal.jl shave <ms> <sources> [--beam BEAM]
        [--peeliter PITERS] [--maxiter ITERS] [--tolerance TOL] [--minuvw MINUVW]
    ttcal.jl zest <ms> <sources> [--beam BEAM]
        [--peeliter PITERS] [--maxiter ITERS] [--tolerance TOL] [--minuvw MINUVW]
    ttcal.jl prune <ms> <sources> [--beam BEAM]
        [--peeliter PITERS] [--maxiter ITERS] [--tolerance TOL] [--minuvw MINUVW]

Options:
    --force-imaging
        Create and use the MODEL_DATA and CORRECTED_DATA columns even if they do not already exist
        in the measurement set.
    --corrected
        Apply the calibration to the CORRECTED_DATA column instead of the DATA column.
    --beam BEAM
        The name of the beam model to use. Documentation for the available beam models is available
        at http://mweastwood.info/TTCal.jl/beammodel/
    --maxiter ITERS
        The maximum number of (Mitch|Stef)cal iterations to take on solving for the calibration on a
        single frequency channel before moving on to the next frequency channel. Channels that fail
        to reach the convergence criterion before the maximum number of iterations are flagged.
        This parameter defaults to 20.
    --tolerance TOL
        Set the relative tolerance used to determine convergence. Setting this to a larger number
        will generate a rougher calibration (but in less time). If this is set too low, you may find
        that all of your data is flagged because nothing converged within the maximum number of
        iterations. This parameter defaults to 1e-3.
    --peeliter PITERS
        This number defines the number of passes that will be made through the source while peeling.
        If there are too few iterations there may be residual errors associated with the imperfect
        removal of bright sources on the first pass through the source list. This parameter defaults
        to 3.
    --minuvw MINUVW
        The minimum baseline length (measured in wavelengths) to use while peeling sources. This
        parameter can be used to mitigate sensitivity to unmodeled diffuse emission. This parameter
        defaults to 0.
"""

function main(args)
    parsed = docopt(doc, args, version=v"0.3.0")
    if     parsed["applycal"]
        run_applycal(parsed)
    elseif parsed["gaincal"]
        run_gaincal(parsed)
    elseif parsed["polcal"]
        run_polcal(parsed)
    elseif parsed["peel"]
        run_peel(parsed)
    elseif parsed["shave"]
        run_shave(parsed)
    elseif parsed["zest"]
        run_zest(parsed)
    elseif parsed["prune"]
        run_prune(parsed)
    else
        error("How did you get here?")
    end
end

macro cli_intro(name)
    quote
        println("Running `", $(string(name)), "` on ", args["<ms>"])
    end |> esc
end

macro cli_load_ms()
    quote
        ms = Table(ascii(args["<ms>"]))
        meta = Metadata(ms)
    end |> esc
end

macro cli_load_beam()
    quote
        if args["--beam"] === nothing
            beam = SineBeam()
        else
            beam = select_beam(args["--beam"])
        end
        dump(beam)
    end |> esc
end

macro cli_load_sources()
    quote
        sources = readsources(args["<sources>"])
    end |> esc
end

macro cli_convergence_criteria()
    quote
        maxiter = args["--maxiter"] === nothing? 20 : parse(Int, args["--maxiter"])
        tolerance = args["--tolerance"] === nothing? 1e-3 : parse(Float64, args["--tolerance"])
        peeliter = args["--peeliter"] === nothing? 3 : parse(Int, args["--peeliter"])
        minuvw = args["--minuvw"] === nothing? 0.0 : parse(Float64, args["--minuvw"])
    end |> esc
end

macro cli_cleanup()
    quote
        unlock(ms)
        nothing
    end
end

function run_applycal(args)
    @cli_intro applycal
    @cli_load_ms
    cal = read(args["<cal>"])
    write_to_corrected = args["--force-imaging"] || Tables.exists(ms, "CORRECTED_DATA")
    apply_to_corrected = args["--corrected"] && Tables.exists(ms, "CORRECTED_DATA")
    data = apply_to_corrected? read(ms, "CORRECTED_DATA") : read(ms, "DATA")
    applycal!(data, meta, cal)
    write_to_corrected? write(ms, "CORRECTED_DATA", data) : write(ms, "DATA", data)
    @cli_cleanup
end

function run_gaincal(args)
    @cli_intro gaincal
    @cli_load_ms
    @cli_load_sources
    @cli_load_beam
    @cli_convergence_criteria
    data = read(ms, "DATA")
    flag_short_baselines!(data, meta, minuvw)
    cal = gaincal(data, meta, beam, sources, maxiter=maxiter, tolerance=tolerance)
    write(args["<cal>"], cal)
    @cli_cleanup
end

function run_polcal(args)
    @cli_intro polcal
    @cli_load_ms
    @cli_load_sources
    @cli_load_beam
    @cli_convergence_criteria
    data = Tables.exists(ms, "CORRECTED_DATA")? read(ms, "CORRECTED_DATA") : read(ms, "DATA")
    flag_short_baselines!(data, meta, minuvw)
    cal = polcal(data, meta, beam, sources, maxiter=maxiter, tolerance=tolerance)
    write(args["<cal>"], cal)
    @cli_cleanup
end

macro peel_input()
    quote
        ms = Table(ascii(args["input"]))
        sources = readsources(args["sources"])
        beam = beam_dictionary[args["beam"]]()
        meta = Metadata(ms)
        data = Tables.exists(ms, "CORRECTED_DATA")? read(ms, "CORRECTED_DATA") : read(ms, "DATA")
        flag_short_baselines!(data, meta, args["minuvw"])
        peeliter = args["peeliter"]
        maxiter = args["maxiter"]
        tolerance = args["tolerance"]
    end |> esc
end

for (routine, T) in ((:peel, PeelingSource), (:shave, ShavingSource),
                     (:zest, ZestingSource), (:prune, PruningSource))
    func = Symbol("run_", routine)
    @eval function $func(args)
        @cli_intro $routine
        @cli_load_ms
        @cli_load_sources
        @cli_load_beam
        @cli_convergence_criteria
        data = Tables.exists(ms, "CORRECTED_DATA")? read(ms, "CORRECTED_DATA") : read(ms, "DATA")
        flag_short_baselines!(data, meta, minuvw)
        peelingsources = $T[$T(source) for source in sources]
        calibrations = peel!(data, meta, beam, peelingsources,
                             peeliter=peeliter, maxiter=maxiter, tolerance=tolerance)
        Tables.exists(ms, "CORRECTED_DATA")? write(ms, "CORRECTED_DATA", data) : write(ms, "DATA", data)
        @cli_cleanup
    end
end

function select_beam(str)
    dictionary = Dict("constant" => ConstantBeam,
                      "sin"      => SineBeam,
                      "sine"     => SineBeam,
                      "memo178"  => Memo178Beam)
    if haskey(dictionary, str)
        return dictionary[str]()
    else
        m = match(r"sin\^([-+]?\d*\.?\d+)", str)
        if m !== nothing
            return SineBeam(parse(Float64, m.captures[1]))
        end
        error("Unknown beam model.")
    end
end

precompile(main, (Dict{AbstractString,Any},))
precompile(run_applycal, (Dict{AbstractString,Any},))
precompile(run_gaincal, (Dict{AbstractString,Any},))
precompile(run_polcal, (Dict{AbstractString,Any},))
precompile(run_peel, (Dict{AbstractString,Any},))
precompile(run_shave, (Dict{AbstractString,Any},))
#precompile(run_zest, (Dict{AbstractString,Any},))
#precompile(run_prune, (Dict{AbstractString,Any},))

