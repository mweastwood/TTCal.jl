using ArgParse

import TTCal

function parse()
    s = ArgParseSettings()
    s.description = "A calibration routine developed for the OVRO LWA."
    @add_arg_table s begin
        "gaintable"
            help = "Gains will be written to/read from this file."
            required = true
        "measurementsets"
            help = "List of measurement sets to use in the calibration"
            required = true
            nargs = '+'
        "--applycal", "-a"
            help = "Apply the given calibration without solving for a new calibration"
            action = :store_true
        "--doubleprecision"
            help = "Solve for 64+64 bit complex gains (instead of 32+32 bit complex gains)."
            action = :store_true
        "--flags"
            help = "List of antennas to flag. Antennas are numbered starting from 1."
            arg_type = Int
            nargs = '+'
        "--niter", "-n"
            help = "The maximum number of iterations used to achieve convergence on a given frequency channel."
            arg_type = Int
            default = 30
        "--refant", "-r"
            help = "The reference antenna."
            arg_type = Int
            default = 1
        "--RK"
            help = "The order of the Runge-Kutta method to use in the inner calibration loop."
            arg_type = Int
            default = 4
        "--tol", "-t"
            help = "The desired relative tolerance to be used while checking for convergence."
            arg_type = Float64
            default = 1e-5
    end
    parse_args(s)
end

function main()
    args = parse()
    TTCal.run(args)
    nothing
end
main()

