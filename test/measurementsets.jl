@testset "measurementsets.jl" begin
    Nant  = 5
    Nfreq = 10
    name, ms = createms(Nant, Nfreq)
    meta = collect_metadata(ms, ConstantBeam())
    @test TTCal.Nant(meta)  == Nant
    @test TTCal.Nfreq(meta) == Nfreq
    @test TTCal.Nbase(meta) == (Nant*(Nant+1))÷2

    data_flags = zeros(Bool,4,5,6)
    row_flags = ones(Bool,6)
    @test TTCal.resolve_flags(data_flags, row_flags) == ones(Bool,6,5)
    data_flags = ones(Bool,4,5,6)
    row_flags = zeros(Bool,6)
    @test TTCal.resolve_flags(data_flags, row_flags) == ones(Bool,6,5)

    data = Visibilities(zeros(JonesMatrix, TTCal.Nbase(meta), Nfreq),
                        zeros(       Bool, TTCal.Nbase(meta), Nfreq))

    # no flags
    data.flags[:] = false
    flag_short_baselines!(data, meta, 0)
    @test sum(data.flags) == 0

    # flag just the autocorrelations
    data.flags[:] = false
    flag_short_baselines!(data, meta, 1e-5)
    @test sum(data.flags) == Nant*Nfreq

    # flag everything
    data.flags[:] = false
    flag_short_baselines!(data, meta, 100)
    @test sum(data.flags) == TTCal.Nbase(meta)*Nfreq
end

