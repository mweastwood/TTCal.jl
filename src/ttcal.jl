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

push!(CLI.commands,Command("flagdata","Flag antennas within a measurement set."))
push!(CLI.commands,Command("bandpass","Solve for a bandpass calibration."))
push!(CLI.commands,Command("applycal","Apply a calibration."))

CLI.options["flagdata"] = [
    Option("--input","""
        The list of measurement sets to flag.""",
        UTF8String,true,1,Inf)
    Option("--antennas","""
        The list of antennas that will be flagged within the supllied
        list of measurement sets. The first antenna is number one.""",
        Int,true,1,Inf)]
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
        Float64,false,1,1)]
CLI.options["applycal"] = [
    Option("--input","""
        The list of measurement sets that the calibration will be applied to.""",
        UTF8String,true,1,Inf),
    Option("--calibration","""
        The file containing the calibration solution.""",
        UTF8String,true,1,1)]

import TTCal

# Catch exceptions to hide the verbose Julia output that is
# not especially helpful to users.
try
    command,args = CLI.parse_args(ARGS)
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

