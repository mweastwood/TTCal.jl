@testset "subsrc.jl" begin
    Nant = 10
    Nfreq = 2

    name,ms = createms(Nant,Nfreq)
    meta = Metadata(ms)
    beam = ConstantBeam()
    sources = readsources("sources.json")
    data = genvis(meta, beam, sources)
    data′ = deepcopy(data)

    subsrc!(data′, meta, beam, sources)
    @test all(data′.data .== zero(JonesMatrix))
    @test all(data′.flags .== false)

    putsrc!(data′, meta, beam, sources)
    @test all(data′.data .== data.data)
    @test all(data′.flags .== false)
end

