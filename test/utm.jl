let
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
end

