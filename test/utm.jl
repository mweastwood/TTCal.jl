@testset "utm.jl" begin
    # Coordinates for the CN tower
    zone = 17
    latitude = 43 + 38/60 + 33.24/3600
    longitude = -(79 + 23/60 + 13.7/3600)
    easting = 630084
    northing = 4833439

    E,N = TTCal.latlong_to_utm(zone,latitude,longitude)
    @test round(Int,E) == easting
    @test round(Int,N) == northing

    lat,long = TTCal.utm_to_latlong(zone,easting,northing)
    @test abs( lat -  latitude) < 0.1/3600
    @test abs(long - longitude) < 0.1/3600


    frame  = ReferenceFrame()
    position = observatory("VLA")
    time = Epoch(epoch"UTC", (2015-1858)*365days)
    set!(frame, position)
    set!(frame, time)
    zenith = Direction(dir"AZEL",   0degrees, 90degrees)
    north  = Direction(dir"AZEL",   0degrees,  0degrees)
    east   = Direction(dir"AZEL",  90degrees,  0degrees)
    south  = Direction(dir"AZEL", 180degrees,  0degrees)
    west   = Direction(dir"AZEL", 270degrees,  0degrees)
    @test measure(frame, zenith, dir"ITRF") ≈ TTCal.azel_to_itrf(zenith, position)
    @test measure(frame,  north, dir"ITRF") ≈ TTCal.azel_to_itrf(north, position)
    @test measure(frame,   east, dir"ITRF") ≈ TTCal.azel_to_itrf(east, position)
    @test measure(frame,  south, dir"ITRF") ≈ TTCal.azel_to_itrf(south, position)
    @test measure(frame,   west, dir"ITRF") ≈ TTCal.azel_to_itrf(west, position)
end

