#!/usr/bin/env julia

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

using CLI

CLI.set_name("ttcal.jl")
CLI.set_banner("""
     TTCal
    =======
    A calibration routine developed for the OVRO LWA.
    Written by Michael Eastwood (mweastwood@astro.caltech.edu).
    """)

CLI.print_banner()

push!(CLI.commands,Command("bandpass","Solve for a bandpass calibration."))
push!(CLI.commands,Command("polcal","Solve for a polarization calibration."))
push!(CLI.commands,Command("applycal","Apply a calibration."))

CLI.options["bandpass"] = [
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
        Float64,false,1,1),
    Option("--force-imaging","""
        Create and use the MODEL_DATA and CORRECTED_DATA columns even if
        they do not already exist in the measurement set.""",
        Nothing,false,0,0),
    Option("--model-present","""
        Use this flag to indicate that there is already a set of model
        visibilities present in the MODEL_DATA column of the measurement
        set. TTCal will use these model visibilities in place of the
        source model specified by the --sources flag (for now --sources
        is still mandatory even though it is effectively ignored).""",
        Nothing,false,0,0)]
CLI.options["polcal"] = [
    Option("--input","""
        The measurement set used to solve for a polarization calibration.""",
        UTF8String,true,1,1),
    Option("--output","""
        The output file to which the polarization calibration will be written.
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
CLI.options["applycal"] = [
    Option("--input","""
        The list of measurement sets that the calibration will be applied to.""",
        UTF8String,true,1,Inf),
    Option("--calibration","""
        The file containing the calibration solution.""",
        UTF8String,true,1,1),
    Option("--force-imaging","""
        Write the calibrated visibilities to the CORRECTED_DATA column
        regardless of whether or not the measurement set already has
        this column (it will be created if it doesn't exist).""",
        Nothing,false,0,0),
    Option("--corrected","""
        Apply the calibration to the CORRECTED_DATA column instead
        of the DATA column.""",
        Nothing,false,0,0)]

import TTCal

# Catch exceptions to hide the verbose Julia output that is
# not especially helpful to users.
# (comment out the try/catch block for debugging)
try
    command,args = CLI.parse_args(ARGS)
    if     command == "bandpass"
        TTCal.run_bandpass(args)
    elseif command == "polcal"
        TTCal.run_polcal(args)
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

