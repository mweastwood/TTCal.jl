@testset "subsrc.jl" begin
    Nant = 10
    Nfreq = 2
    name,ms = createms(Nant,Nfreq)
    meta = TTCal.collect_metadata(ms, ConstantBeam())
    sources = readsources("sources.json")
    data = genvis(meta, sources)
    subsrc!(data, meta, sources)
    @test data.data == zeros(JonesMatrix, TTCal.Nbase(meta), Nfreq)
end

