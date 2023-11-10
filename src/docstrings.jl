@doc """
    abstract type Coordinate{N, T <: Number}

`N`-dimensional coordinate with eltype `T`.
""" Coordinate

@doc """
    abstract type EastWestCoordinate{T} <: Coordinate{1, T}

One-dimensional coordinate in east-west direction.
""" EastWestCoordinate

@doc """
    abstract type NorthSouthCoordinate{T} <: Coordinate{1, T}

One-dimensional coordinate in north-south direction.
""" NorthSouthCoordinate

###############################################################################

@doc """
    WebMercator{T} <: Coordinate{2,T}

WebMercator coordinates, shifted and scaled such that we have the following
identities.

```jldoctest; output = false
julia> @assert WebMercator(LonLat( 0, 0 ))[] == (0, 0)

julia> @assert WMX(Lon(180))[] == +1

julia> @assert WMY(Lat(90))[] == +Inf
```
""" WebMercator

@doc """
    WMX{T} <: EastWestCoordinate{T}

x-component of [`WebMercator`](@ref) coordinates.
""" WMX

@doc """
    WMY{T} <: NorthSouthCoordinate{T}

y-component of [`WebMercator`](@ref) coordinates.
""" WMY

###############################################################################

@doc """
    LatLon{T} <: Coordinate{2,T}

Geodetic latitude and longitude in degrees.
""" LatLon

@doc """
    LonLat{T} <: Coordinate{2,T}

Geodetic longitude and latitude in degrees.
""" LonLat

@doc """
    Lat{T} <: NorthSouthCoordinate{T}

Geodetic latitude in degrees.
""" Lat

@doc """
    Lon{T} <: EastWestCoordinate{T}

Geodetic longitude in degrees.
""" Lon

###############################################################################

@doc """
    EastNorth{T} <: Coordinate{2,T}

Easting and northing in meters.

 - Easting denotes the signed distance (in meters) to the prime meridian.
 - Northing denotes the signed distance (in meters) to the equator.
""" EastNorth

@doc """
    East{T} <: EastWestCoordinate{T}

Easting, i.e. signed distance (in meters) to the prime meridian.
""" East

@doc """
    North{T} <: NorthSouthCoordinate{T}

Northing, i.e. signed distance (in meters) to the equator.
""" North

###############################################################################

@doc """
    ECEF{T} <: Coordinate{3}

Earth-centred earth-fixed coordinates in meters.

# Example

```jldoctest
julia> ECEF(LonLat(0,0))
ECEF{Float64}(6.378137e6, 0.0, 0.0)

julia> ECEF(LonLat(90,0))
ECEF{Float64}(0.0, 6.378137e6, 0.0)

julia> ECEF(LonLat(0,90))
ECEF{Float64}(0.0, 0.0, 6.356752314245179e6)
```
""" ECEF

@doc """
    tdist(a,b) -> Number

Compute the "tunnel distance" between points `a` and `b`.

Tunnel distance is the straight-line distance in Euclidean space, i.e.
`norm(ECEF(a) - ECEF(b))`.

Each of `a` and `b` may be any of the following:
- A pair of north-south and east-west coordinates.
- A two-dimensional coordinate.
- Either of the above with an additional altitude.

# Example

```jldoctest
julia> tdist(EastNorth(0,0), EastNorth(3,4))
5.0000000000001075

julia> tdist(EastNorth(0,0), (EastNorth(3,0), Alt(4)))
5.0000005648476336

julia> tdist(EastNorth(0,0), (East(3), North(0), Alt(4)))
5.0000005648476336
```
""" tdist