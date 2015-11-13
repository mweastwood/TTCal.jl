let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant

    name,ms = createms(Nant,Nfreq)

    # Put a source right at the phase center such that all the visibilities
    # should be unity

    source = Source("Cristiano",
                    Point("Ronaldo",
                          ms.phase_direction,
                          Spectrum(1,0,0,0,10e6,[0.0])))
    visibilities = genvis(ms, source, ConstantBeam())

    model = zeros(Complex64,4,Nfreq,Nbase)
    model[1,:,:] = 1
    model[4,:,:] = 1
    @test visibilities ≈ model
end

