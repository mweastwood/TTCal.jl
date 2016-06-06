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

abstract BeamModel

function call(beam::BeamModel, ν, direction::Direction)
    direction.sys == dir"AZEL" || error("Direction must be in the AZEL coordinate system.")
    azimuth   = longitude(direction)
    elevation =  latitude(direction)
    beam(ν, azimuth, elevation)
end

"""
    ConstantBeam <: BeamModel

In this beam model the Jones matrix is assumed to be the identity
in every direction. This is the simplest possible beam model
and should be used if you wish to avoid the application of
a beam model.

For example, if you wish to calculate model visibilities but
you have already manually attenuating the flux by a beam model,
you should use `ConstantBeam` to avoid applying a beam model again.

    beam = ConstantBeam()
    sources = readsources("mysources.json")
    # ...
    # attenuate the sources with your own custom beam model here
    # ...
    metadata = collect_metadata(measurement_set, beam)
    model_visibilities = genvis(metadata, sources)
"""
immutable ConstantBeam <: BeamModel end
call(::ConstantBeam, ν, az, el) = one(JonesMatrix)

doc"""
    SineBeam <: BeamModel

This beam is azimuthally symmetric and independent of frequency.
The gain of an antenna scales as $\sin(\rm el)^\alpha$.
For an LWA dipole a reasonable approximation to the the antenna gain
is $\sin(\rm el)^{1.6}$, so $1.6$ is the default value of $\alpha$.

Note that because Jones matrices operate on the electric field
which must be squared to get power, the diagonal elements of the Jones matrix
are $\sin(\rm el)^{\alpha/2}$.

    beam = SineBeam() # equivalent to SineBeam(1.6)
    beam = SineBeam(2.0)
"""
immutable SineBeam <: BeamModel
    α :: Float64
end

SineBeam() = SineBeam(1.6)

function call(beam::SineBeam, ν, az, el)
    # note that the factor of 1/2 in the exponent
    # comes from the fact that the Jones matrix
    # operates on electric fields, which must be
    # squared to get the power
    s = sin(el)^(beam.α/2)
    JonesMatrix(s, 0, 0, s)
end

doc"""
    Memo178Beam <: BeamModel

This beam is based on the parametric fit to EM simulations
presented in [LWA memo 178](http://www.faculty.ece.vt.edu/swe/lwa/memo/lwa0178a.pdf)
by Jayce Dowell.

The dipole gain is expressed as

\\\\[
    p(\theta) = \left[1-\left(\frac{\theta}{\pi/2}\right)^\alpha\right]\cos^\beta\theta
              + \gamma\left(\frac{\theta}{\pi/2}\right)\cos^\delta\theta,
\\\\]

where $\theta$ is the zenith angle, and $\alpha$, $\beta$, $\gamma$, and $\delta$ are
parameters that were fit for in the E- and H-planes of the dipole.
"""
immutable Memo178Beam <: BeamModel end

function call(::Memo178Beam, ν, az, el)
    # NOTE: I might have the x and y dipoles swapped here
    x = P178(ν, az, el) |> sqrt
    y = P178(ν, az+π/2, el) |> sqrt
    JonesMatrix(x, 0, 0, y)
end

const _E178 = [-4.529931167425190e+01 -3.066691727279789e+01 +7.111192148086860e+01 +1.131338637814271e+01
               +1.723596273204143e+02 +1.372536555724785e+02 -2.664504470520252e+02 -3.493942140370373e+01
               -2.722311453669980e+02 -2.368121452949910e+02 +4.343953141004854e+02 +5.159758406047006e+01
               +2.408402047219155e+02 +2.302040670131884e+02 -3.993175413350526e+02 -4.152532958513045e+01
               -1.334589043702679e+02 -1.414814115813401e+02 +2.312015990805256e+02 +2.016753495756661e+01
               +4.917278320096442e+01 +5.820262041648866e+01 -8.952191930147437e+01 -6.170044495806915e+00
               -1.244683117802972e+01 -1.652458005311876e+01 +2.396407975769994e+01 +1.185153343629695e+00
               +2.195017889732819e+00 +3.281357788749097e+00 -4.500463893128127e+00 -1.311234364955677e-01
               -2.690339381925372e-01 -4.547844648614167e-01 +5.919912952153034e-01 +4.850030701166225e-03
               +2.243604641077215e-02 +4.311447324076417e-02 -5.345256010140458e-02 +7.056536178294128e-04
               -1.211643367267544e-03 -2.665185970775750e-03 +3.157455237232649e-03 -1.102532965104807e-04
               +3.810179416998095e-05 +9.681257383030349e-05 -1.099243217364324e-04 +6.250327364174989e-06
               -5.277176362640874e-07 -1.567746596003443e-06 +1.710529173897111e-06 -1.350506188926788e-07]

const _H178 = [+4.062920357822495e+02 +3.038713068453467e+01
               -1.706845366994521e+03 -1.337217221068207e+02
               +3.095438045596764e+03 +2.567017051330929e+02
               -3.164514198869798e+03 -2.757538828204631e+02
               +2.035008485840167e+03 +1.850998851042284e+02
               -8.713893456840954e+02 -8.227855478515090e+01
               +2.564477521133275e+02 +2.502479051323083e+01
               -5.262251671400085e+01 -5.287316195660030e+00
               +7.519512104996896e+00 +7.754610944702455e-01
               -7.337746902310935e-01 -7.744604652809532e-02
               +4.663414310713179e-02 +5.024254273102379e-03
               -1.740005497709271e-03 -1.908989166200470e-04
               +2.892116885178882e-05 +3.223985512686652e-06]

@eval function E178(ν, el)
    x = ν / 10e6
    θ = π/2 - el
    α = @evalpoly(x, $(_E178[:,1]...))
    β = @evalpoly(x, $(_E178[:,2]...))
    γ = @evalpoly(x, $(_E178[:,3]...))
    δ = @evalpoly(x, $(_E178[:,4]...))
    (1-(θ/(π/2))^α)*cos(θ)^β + γ*(θ/(π/2))*cos(θ)^δ
end

@eval function H178(ν, el)
    x = ν / 10e6
    θ = π/2 - el
    α = @evalpoly(x, $(_H178[:,1]...))
    β = @evalpoly(x, $(_H178[:,2]...))
    # γ and δ are neglected on the H plane
    (1-(θ/(π/2))^α)*cos(θ)^β
end

function P178(ν, az, el)
    sqrt((E178(ν,el)*cos(az))^2 + (H178(ν,el)*sin(az))^2)
end

doc"""
    ZernikeBeam <: BeamModel

This beam is composed of Zernike polynomials.
"""
immutable ZernikeBeam <: BeamModel end

function call(::ZernikeBeam, ν, az, el)
    ρ = cos(el)
    θ = az
    value = (_Zcoeff[1]*zernike(0, 0, ρ, θ)
            + _Zcoeff[2]*zernike(2, 0, ρ, θ)
            + _Zcoeff[3]*zernike(4, 0, ρ, θ)
            + _Zcoeff[4]*zernike(4, 4, ρ, θ)
            + _Zcoeff[5]*zernike(6, 0, ρ, θ)
            + _Zcoeff[6]*zernike(6, 4, ρ, θ)
            + _Zcoeff[7]*zernike(8, 0, ρ, θ)
            + _Zcoeff[8]*zernike(8, 4, ρ, θ)
            + _Zcoeff[9]*zernike(8, 8, ρ, θ))
    JonesMatrix(sqrt(value), 0, 0, sqrt(value))
end

const _Zcoeff = [ 0.5925713994750834,
                 -0.4622486219893028,
                 -0.054924184973998307,
                 -0.0028805328944419696,
                 -0.02407776673368796,
                 -0.006155457593922782,
                 -0.023973603224075223,
                 -0.003090132046373044,
                  0.00497413312773207]

function zernike(n, m, ρ, θ)
    zernike_radial_part(n, abs(m), ρ) * zernike_azimuthal_part(m, θ)
end

function zernike_radial_part(n, m, ρ)
    R0 = ρ^m
    n == m && return R0
    R2 = ((m+2)*ρ^2 - (m+1))*R0
    for n′ = m+4:2:n
        recurrence_relation = (2*(n-1)*(2n*(n-2)*ρ^2-m^2-n*(n-2))*R2 - n*(n+m-2)*(n-m-2)*R0) / ((n+m)*(n-m)*(n-2))
        R0 = R2
        R2 = recurrence_relation
    end
    R2
end

function zernike_azimuthal_part(m, θ)
    if m == 0
        return 1.0
    elseif m > 0
        return cos(m*θ)
    else
        return sin(m*θ)
    end
end

const beam_dictionary = Dict("constant" => ConstantBeam,
                             "sine"     => SineBeam,
                             "memo178"  => Memo178Beam,
                             "zernike"  => ZernikeBeam)

