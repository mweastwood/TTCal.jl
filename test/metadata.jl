@testset "metadata.jl" begin
    Nant  = 5
    Nfreq = 10
    name, ms = createms(Nant, Nfreq)
    meta = Metadata(ms)
    @test TTCal.Nant(meta)  == Nant
    @test TTCal.Nfreq(meta) == Nfreq
    @test TTCal.Nbase(meta) == (Nant*(Nant+1))÷2
end

