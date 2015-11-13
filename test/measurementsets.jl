let
    u = [1.0; 10.0; 20.0]
    v = [0.0;  0.0;  0.0]
    w = [0.0;  0.0;  0.0]
    ν = [TTCal.c/6.0; TTCal.c/4.0]

    flags = zeros(4,length(ν),length(u))
    TTCal.flag_short_baselines!(flags,2.0,u,v,w,ν)

    manual_flags = zeros(4,length(ν),length(u))
    manual_flags[:,:,1] = true
    manual_flags[:,1,2] = true
    @test flags == manual_flags
end

