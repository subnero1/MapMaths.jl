@doc """
    WebMercator(wmx, wmy)

Web Mercator coordinates, shifted and scaled such that we have the following
identities.

```
WebMercator(LonLat( 0, 0 ))[] == (0, 0)
WMX(Lon(180))[] == 1
WMX(Lon(-180))[] == -1
wMY(Lat(90))[] == Inf
wMY(Lat(-90))[] == -Inf
```
""" WebMercator

@doc """
    WMX(wmx)

East-west component of [`WebMercator`](@ref) coordinates.
""" WMX

@doc """
    WMY(wmy)

North-south component of [`WebMercator`](@ref) coordinates.
""" WMY

###############################################################################

@doc """
    LatLon(lat, lon)

Latitude and longitude in degrees.
""" LatLon

@doc """
    LonLat(lon, lat)

Longitude and latitude in degrees.
""" LonLat

@doc """
    Lat(lat)

Latitude in degrees.
""" Lat

@doc """
    Lon(lon)

Longitude in degrees.
""" Lon

###############################################################################

@doc """
    EastNorth(east, north)

Easting and northing in meters.

Easting denotes the distance one has to travel westward (positive easting) or
eastward (negative easting) to reach zero degree longitude.

Northing denotes the distance one has to travel southward (positive northing) or
northward (negative northing) to reach the equator.
""" EastNorth

@doc """
    East(east)

Easting in meters.

Easting denotes the distance one has to travel westward (positive easting) or
eastward (negative easting) to reach zero degree longitude.
""" East

@doc """
    North(north)

Northing in meters.

Northing denotes the distance one has to travel southward (positive northing) or
northward (negative northing) to reach the equator.
""" North

