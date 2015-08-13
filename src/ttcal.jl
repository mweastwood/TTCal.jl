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

push!(CLI.commands,Command("gaincal","Solve for a gain calibration."))
push!(CLI.commands,Command("polcal","Solve for a polarization calibration."))
push!(CLI.commands,Command("peel","Peel sources from the dataset."))
push!(CLI.commands,Command("applycal","Apply a calibration."))

CLI.options["gaincal"] = [
    Option("--input","""
        The measurement set used to solve for the calibration.""",
        T=UTF8String,
        required=true,
        min=1,
        max=1),
    Option("--output","""
        The output file to which the calibration will be written.
        This output file will be overwritten if it already exists.""",
        T=UTF8String,
        required=true,
        min=1,
        max=1),
    Option("--sources","""
        A JSON file describing the sources to be used for the sky model.""",
        T=UTF8String,
        required=true,
        min=1,
        max=1),
    Option("--maxiter","""
        Set the maximum number of (Stef|Mitch)Cal iterations to take on each
        frequency channel.""",
        T=Int,
        min=1,
        max=1),
    Option("--tolerance","""
        Set the relative tolerance used to determine convergence.""",
        T=Float64,
        min=1,
        max=1),
    Option("--force-imaging","""
        Create and use the MODEL_DATA and CORRECTED_DATA columns even if
        they do not already exist in the measurement set.""")]
CLI.options["polcal"] = CLI.options["gaincal"]
CLI.options["peel"] = [
    Option("--input","""
        The measurement set that will have sources peeled.""",
        T=UTF8String,
        required=true,
        min=1,
        max=1),
    Option("--sources","""
        A JSON file describing the sources to be peeled from the given
        measurement set.""",
        T=UTF8String,
        required=true,
        min=1,
        max=1),
    Option("--minuvw","""
        The minimum baseline length (measured in wavelengths) to use while
        peeling sources. This parameter can be used to mitigate sensitivity
        to unmodeled diffuse emission.""",
        T=Float64,
        min=1,
        max=1)]
CLI.options["applycal"] = [
    Option("--input","""
        The list of measurement sets that the calibration will be applied to.""",
        T=UTF8String,
        required=true,
        min=1,
        max=Inf),
    Option("--calibration","""
        The file containing the calibration solution.""",
        T=UTF8String,
        required=true,
        min=1,
        max=1),
    Option("--force-imaging","""
        Write the calibrated visibilities to the CORRECTED_DATA column
        regardless of whether or not the measurement set already has
        this column (it will be created if it doesn't exist)."""),
    Option("--corrected","""
        Apply the calibration to the CORRECTED_DATA column instead
        of the DATA column.""")]

import TTCal

# Catch exceptions to hide the verbose Julia output that is
# not especially helpful to users.
# (comment out the try/catch block for debugging)
try
    command,args = CLI.parse_args(ARGS)
    if     command == "gaincal"
        TTCal.run_gaincal(args)
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

