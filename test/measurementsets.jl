@testset "measurementsets.jl" begin
    Nant  = 5
    Nfreq = 10
    name, ms = createms(Nant, Nfreq)
    meta = TTCal.collect_metadata(ms, ConstantBeam())
    @test TTCal.Nant(meta)  == Nant
    @test TTCal.Nfreq(meta) == Nfreq
    @test TTCal.Nbase(meta) == (Nant*(Nant+1))÷2

    data_flags = zeros(Bool,4,5,6)
    row_flags = ones(Bool,6)
    @test TTCal.resolve_flags(data_flags, row_flags) == ones(Bool,6,5)
    data_flags = ones(Bool,4,5,6)
    row_flags = zeros(Bool,6)
    @test TTCal.resolve_flags(data_flags, row_flags) == ones(Bool,6,5)
end

