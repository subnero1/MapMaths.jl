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

###############################################################################

@doc """
    ECEF(x,y,z)

Earth-centred earth-fixed coordinates in meters.

Cartesian coordinate system such that `{lon == 0} ⊂ {y == 0}` and `{lat == 0} ⊂
{z == 0}`.

# Example
```
julia> ECEF(LonLat(0,0))
ECEF{Float64}(6.378137e6, 0.0, 0.0)

julia> ECEF(LonLat(90,0))
ECEF{Float64}(0.0, 6.378137e6, 0.0)

julia> ECEF(LonLat(0,90))
ECEF{Float64}(0.0, 0.0, 6.356752314245179e6)
```
""" ECEF

@doc """
    tdist(c1, c2) -> Number

Compute the "tunnel distance" between points `c1` and `c2`.

Tunnel distance is the straight-line distance in Euclidean space, i.e.
`norm(ECEF(c1) - ECEF(c2))`.

Each of `c1` and `c2` may be any of the following:
- A pair of north-south and east-west coordinates.
- A two-dimensional coordinate.
- Either of the above with an additional altitude.

# Example
```
julia> tdist(EastNorth(0,0), EastNorth(3,4))
5.0000000000001075

julia> tdist(EastNorth(0,0), (EastNorth(3,0), Alt(4)))
5.0000005648476336

julia> tdist(LatLon(0,0), LatLon(0,180)) # Two times the earth radius
1.2756274e7
```
""" tdist