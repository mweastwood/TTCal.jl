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

"""
    GainCalibration <: Calibration

This type stores the information for calibrating the
electronic gains of the interferometer. That is, it stores
complex gains and flags for each antenna, frequency channel,
and polarization.
"""
immutable GainCalibration <: Calibration
    gains::Array{Complex64,3}
    flags::Array{Bool,3}
end

"""
    GainCalibration(Nant,Nfreq)

Create a calibration table for `Nant` antennas with
`Nfreq` frequency channels where all the gains are
initially set to unity.
"""
function GainCalibration(Nant,Nfreq)
    gains = ones(Complex64,Nant,Nfreq,2)
    flags = zeros(Bool,Nant,Nfreq,2)
    GainCalibration(gains,flags)
end

Nant(cal::GainCalibration) = size(cal.gains,1)
Nfreq(cal::GainCalibration) = size(cal.gains,2)

doc"""
    invert(cal::GainCalibration)

Returns the inverse of the given calibration.
The complex gain $g$ of each antenna is set to $1/g$.
"""
function invert(cal::GainCalibration)
    output = GainCalibration(Nant(cal),Nfreq(cal))
    for i in eachindex(cal.gains)
        output.gains[i] = inv(cal.gains[i])
    end
    output
end

"""
    fixphase!(cal::GainCalibration,reference_antenna)

Set the phase of the reference antenna to zero.
"""
function fixphase!(cal::GainCalibration,
                   reference_antenna)
    for pol = 1:2, β = 1:Nfreq(cal)
        factor = (conj(cal.gains[reference_antenna,β,pol])
                    / abs(cal.gains[reference_antenna,β,pol]))
        for ant = 1:Nant(cal)
            cal.gains[ant,β,pol] = cal.gains[ant,β,pol]*factor
        end
    end
end

"""
    gaincal(ms::Table,sources::Vector{PointSource};
            maxiter = 20, tolerance = 1e-5,
            force_imaging_columns = false,
            reference_antenna = 1)

Solve for the interferometer's electronic gains.
"""
function gaincal(ms::Table,
                 sources::Vector{PointSource};
                 maxiter::Int = 20,
                 tolerance::Float64 = 1e-5,
                 force_imaging_columns::Bool = false,
                 reference_antenna::Int = 1)
    phase_dir = MeasurementSets.phase_direction(ms)
    u,v,w = MeasurementSets.uvw(ms)
    ν = MeasurementSets.frequency(ms)
    ant1,ant2 = MeasurementSets.antennas(ms)

    frame = ReferenceFrame()
    set!(frame,MeasurementSets.position(ms))
    set!(frame,MeasurementSets.time(ms))
    sources = filter(source -> isabovehorizon(frame,source),sources)

    Nant  = numrows(Table(ms[kw"ANTENNA"]))
    Nfreq = length(ν)
    calibration = GainCalibration(Nant,Nfreq)

    data  = ms["DATA"]
    model = genvis(frame,phase_dir,sources,u,v,w,ν)
    flags = MeasurementSets.flags(ms)

    if force_imaging_columns || Tables.checkColumnExists(ms,"MODEL_DATA")
        ms["MODEL_DATA"] = model
    end

    gaincal!(calibration,data,model,flags,
             ant1,ant2,maxiter,tolerance,
             reference_antenna)
    calibration
end

function gaincal!(calibration,data,model,flags,
                  ant1,ant2,maxiter,tolerance,
                  reference_antenna)
    for β = 1:Nfreq(calibration)
        # Calibrate the X polarization
        gaincal_onechannel!(slice(calibration.gains,:,β,1),
                            slice(calibration.flags,:,β,1),
                            slice(data, 1,β,:),
                            slice(model,1,β,:),
                            slice(flags,1,β,:),
                            ant1, ant2,
                            maxiter, tolerance)
        # Calibrate the Y polarization
        gaincal_onechannel!(slice(calibration.gains,:,β,2),
                            slice(calibration.flags,:,β,2),
                            slice(data, 4,β,:),
                            slice(model,4,β,:),
                            slice(flags,4,β,:),
                            ant1, ant2,
                            maxiter, tolerance)
    end
    fixphase!(calibration,reference_antenna)
end

function gaincal_onechannel!(gains, gain_flags,
                             data, model, data_flags,
                             ant1, ant2,
                             maxiter, tolerance)
    # If the entire channel is flagged, don't bother calibrating.
    if all(data_flags)
        gain_flags[:] = true
        return
    end

    square_data  = gaincal_makesquare( data,data_flags,ant1,ant2)
    square_model = gaincal_makesquare(model,data_flags,ant1,ant2)
    best_gains   = gaincal_firstguess(square_data,square_model)
    converged = @iterate(GainCalStep(),RK4,maxiter,tolerance,
                         best_gains,square_data,square_model)

    # Flag the entire channel if the solution did not converge.
    if !converged
        gain_flags[:] = true
        return
    end

    # Output
    for ant = 1:length(gains)
        gains[ant] = best_gains[ant]
    end
end

function gaincal_makesquare(data,flags,ant1,ant2)
    Nbase = length(data)
    Nant  = round(Integer,div(sqrt(1+8Nbase)-1,2))
    output = zeros(Complex64,Nant,Nant)
    for α = 1:Nbase
        flags[α] && continue
        ant1[α] == ant2[α] && continue
        output[ant1[α],ant2[α]] = data[α]
        output[ant2[α],ant1[α]] = conj(data[α])
    end
    output
end

doc"""
    gaincal_firstguess(data,model)

Take a rough first guess at the visibilities.

This function works by looking for the principle eigenvector
of the matrix $G[i,j] = data[i,j]/model[i,j]$. In the
absence of noise, and if the model is complete, then $G = g g^{*}$,
where $g$ is the vector of complex gains.

The main confounding factor for this method is the presence of
an additive noise term on the diagonal (from the autocorrelations).
Most calibration routines neglect the autocorrelations for this
reason, but there is no obvious way to exclude them in this case.
Hence we use this as a first guess to an iterative method that
refines the calibration.
"""
function gaincal_firstguess(data,model)
    G = similar(data)
    for i in eachindex(G)
        G[i] = model[i] == 0? 0 : data[i]/model[i]
    end
    λ,v = eigs(G,nev=1,which=:LM)
    w = sqrt(λ[1])
    squeeze(v,2)*w
end

doc"""
    gaincal_step(gains,data,model) -> step

Given the `data` and `model` visibilities, and the current
guess for the electronic `gains`, solve for `step` such
that the new value of the gains is `gains+step`.

The update step is defined such that the new value of the
gains minimizes
$\sum_{i,j}\|data[i,j] - gains[i,j] conj(newgains[i,j]) model[i,j]\|^2$,
where $i$ and $j$ label the antennas.

*References*

* Michell, D. et al. 2008, JSTSP, 2, 5
* Salvini, S. & Wijnholds, S. 2014, A&A, 571, 97
"""
function gaincal_step(input,data,model)
    Nant = length(input)
    step = zeros(Complex64,Nant)
    @inbounds for j = 1:Nant
        numerator   = zero(Complex64)
        denominator = zero(Complex64)
        for i = 1:Nant
            GM = input[i]*model[i,j]
            numerator = numerator + GM'*data[i,j]
            denominator = denominator + GM'*GM
        end
        ok = abs(denominator) > eps(Float32)
        step[j] = ifelse(ok,conj(numerator/denominator) - input[j],0)
    end
    step
end

immutable GainCalStep <: StepFunction end
call(::GainCalStep,g,V,M) = gaincal_step(g,V,M)

function corrupt!(data::Array{Complex64,3},
                  flags::Array{Bool,3},
                  cal::GainCalibration,
                  ant1,ant2)
    Nbase = length(ant1)
    for α = 1:Nbase, β = 1:Nfreq(cal)
        if (cal.flags[ant1[α],β,1] || cal.flags[ant1[α],β,2]
                                   || cal.flags[ant2[α],β,1]
                                   || cal.flags[ant2[α],β,2])
            flags[:,β,α] = true
        end
        data[1,β,α] = cal.gains[ant1[α],β,1]*conj(cal.gains[ant2[α],β,1])*data[1,β,α]
        data[2,β,α] = cal.gains[ant1[α],β,1]*conj(cal.gains[ant2[α],β,2])*data[2,β,α]
        data[3,β,α] = cal.gains[ant1[α],β,2]*conj(cal.gains[ant2[α],β,1])*data[3,β,α]
        data[4,β,α] = cal.gains[ant1[α],β,2]*conj(cal.gains[ant2[α],β,2])*data[4,β,α]
    end
end

