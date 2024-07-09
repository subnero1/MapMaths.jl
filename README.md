# MapMaths.jl

Utility tools for working with various planetary coordinate systems.

## Examples

- Convert latitude and longitude to WebMercator coordinates.

  ```jldoctest MapMaths
  julia> WebMercator(LatLon(1.3, 103.8))
  WebMercator{Float64}(0.5766666666666667, 0.007222841972793043)
  ```

- Move a latitude / longitude coordinate  1 km east and 0.5 km north.

  ```jldoctest MapMaths
  julia> LatLon(1.3, 103.8) + EastNorth(1000, 500)
  LatLon{Float64}(1.3045218239325143, 103.80898545013565)
  ```

- How non-conformal are the WebMercator coordinates at various latitudes?

  ```jldoctest MapMaths
  julia> lats = [0, 30, 60]
         ratios = [
             WMX(East(1), Lat(lat))[] / WMY(North(1), Lat(lat))[]
             for lat in lats
         ]
         pretty_table([ lats ratios ]; header = ["Latitude", "Scale ratio"])
  ┌──────────┬─────────────┐
  │ Latitude │ Scale ratio │
  ├──────────┼─────────────┤
  │      0.0 │    0.993306 │
  │     30.0 │    0.994971 │
  │     60.0 │    0.998318 │
  └──────────┴─────────────┘
  ```

## Alternatives

There are two other Julia package which address similar needs as this one, namely [CoordRefSystems](https://github.com/JuliaEarth/CoordRefSystems.jl) and [Geodesy](https://github.com/JuliaGeo/Geodesy.jl). The below table attempts to summarise the key differences between them.

| | MapMaths | [CoordRefSystems](https://github.com/JuliaEarth/CoordRefSystems.jl) | [Geodesy](https://github.com/JuliaGeo/Geodesy.jl) |
|:-:|:-:|:-:|:-:|
| **Ease of use** | ✅ | ✅ | Awkward syntax for coordinate conversions. <br> Doesn't automatically promote ints to floats. |
| **Global coordinate types** | LatLon, WebMercator, ECEF | Too many to list them all. | LatLon, ECEF, UTM |
| **Local coordinates** | Curved (tangential is WIP) | :x: | Tangential (ENU) |
| **Datum support** | WGS84 only | Many | Many |
| **Unitful** | :x: | :white_check_mark: | :x: |
| **Reference package** | :x: | [PROJ](https://github.com/OSGeo/PROJ) | [GeographicLib](http://geographiclib.sourceforge.net/) |

We hope that eventually these three packages can evolve into a single one that is best at everything. Until then, we recommend [CoordRefSystems](https://github.com/JuliaEarth/CoordRefSystems.jl) as the default choice due to its professional maintenance and long feature list, but we also warmly invite you to keep an eye on this package as a testbed for highly powerful and flexible coordinate arithmetic.


## API

### Types

##### Abstract coordinate types

- `Coordinate{N,T}`: `N`-dimensional coordinate with eltype `T`.
- `EastWestCoordinate{T} <: Coordinate{1,T}`: One-dimensional coordinate in east-west direction.
- `NorthSouthCoordinate{T} <: Coordinate{1,T}`: One-dimensional coordinate in north-south direction.

##### Latitude / longitude coordinate types

- `LatLon{T} <: Coordinate{2,T}`: [Geodetic](https://en.wikipedia.org/wiki/Geodetic_coordinates) latitude and longitude in degrees.
- `LonLat{T} <: Coordinate{2,T}`: As above, but with arguments reversed.
- `Lat{T} <: NorthSouthCoordinate{T}`: Latitude.
- `Lon{T} <: EastWestCoordinate{T}`: Longitude.

##### WebMercator coordinate types

- `WebMercator{T} <: Coordinate{2,T}`: [WebMercator](https://en.wikipedia.org/wiki/Web_Mercator_projection) coordinates.
- `WMX{T} <: EastWestCoordinate{T}`: x-component of WebMercator coordinates.
- `WMY{T} <: NorthSouthCoordinate{T}`: y-component of WebMercator coordinates.

WebMercator coordinates are shifted and scaled such that we have the following identities.
```jldoctest MapMaths
julia> @assert WebMercator(LonLat( 0, 0 ))[] == (0, 0)

julia> @assert WMX(Lon(180))[] == +1

julia> @assert WMY(Lat(90))[] == +Inf
```

##### Easting / northing coordinate types

- `EastNorth{T} <: Coordinate{2,T}`: Easting and northing in meters (see below).
- `East{T} <: EastWestCoordinate{T}`: Easting, i.e. signed distance (in meters) to the prime meridian.
- `North{T} <: NorthSouthCoordinate{T}`: Northing, i.e. signed distance (in meters) to the equator.

Example:

```jldoctest MapMaths
julia> EastNorth(LatLon(0,90))
EastNorth{Float64}(1.0018754171394622e7, 0.0)

julia> EastNorth(LatLon(90, 0))
EastNorth{Float64}(0.0, 1.0001965729312724e7)
```

(Aside: the above shows that the earth is not a sphere, and that the meter has diverged from its [historical definition](https://en.wikipedia.org/wiki/Metre#Meridional_definition).)

##### Further coordinate types

- `ECEF{T} <: Coordinate{3,T}`: [Earth-centered, earth-fixed](https://en.wikipedia.org/wiki/Earth-centered,_Earth-fixed_coordinate_system) coordinates in meters.
- `Alt{T} <: Coordinate{1,T}`: Altitude (height above the WGS84 ellipsoid) in meters.


### Functions

##### Coordinate conversion

- Any reasonable coordinate conversion can be performed using the constructor of the target coordinate type.

  ```jldoctest MapMaths
  julia> WebMercator(LatLon(0, 45)) # Convert surface coordinates
  WebMercator{Float64}(0.25, 0.0)

  julia> WMX(LatLon(0, 45)) # Extract only east-west coordinate
  WMX{Float64}(0.25)

  julia> WMX(Lon(45)) # Convert one-dimensional coordinates
  WMX{Float64}(0.25)
  ```

  (Actually, some coordinate conversions may be missing. PRs welcome!)

- Converting east-west coordinates sometimes requires latitude information.

  ```jldoctest MapMaths
  julia> East(Lon(45))
  ERROR: Cannot convert from Lon to East without knowing the latitude. Call Lon(::East, ::NorthSouthCoordinate) instead.

  julia> East(Lon(45), Lat(0))
  East{Float64}(5.009377085697311e6)
  ```

- Latitude information can also be provided to north-south coordinate conversions to compute
  ```julia
  X(y,z) = X(typeof(y)(z) + y) - X(z)
  ```

  In general, `X(y,z)` expresses "convert `y` to `X` at latitude `z`".

  Example: One degree latitude is about 110.5 km at the equator, but 111.7 km at the poles.

  ```jldoctest MapMaths
  julia> North(Lat(1), Lat(0))
  North{Float64}(110574.38855779877)

  julia> North(Lat(1), Lat(89))
  North{Float64}(111693.86491418444)
  ```

- Conversion to `ECEF` takes an optional altitude argument.

  ```jldoctest MapMaths
  julia> ECEF(LatLon(0, 90))  # Altitude defaults to 0
  ECEF{Float64}(0.0, 6.378137e6, 0.0)

  julia> ECEF(LatLon(0, 90), Alt(1e7))  # Altitude can also be provided explicitly
  ECEF{Float64}(0.0, 1.6378137e7, 0.0)
  ```

##### Computing distances

- `tdist(x,y) = norm(SVector(ECEF(x) - ECEF(y)))` computes the straight-line distance between two points. The `t` in `tdist()` stands for "tunnel" since traveling along the path measured by `tdist()` typically requires tunneling through the earth. Like `ECEF()`, `tdist()` accepts an optional altitude which defaults to 0.

  ```jldoctest MapMaths
  julia> tdist(EastNorth(0,0), EastNorth(3,0))  # Both altitudes default to 0
  2.9999999999999725

  julia> tdist(
             EastNorth(0,0),
             (EastNorth(3,0), Alt(4))  # Second altitude is set to 4
         )
  5.0000005648476336
  ```

- We may add a function `sdist(x,y)` at some point to compute surface distance (also known as great-circle or geodesic distance).

##### Coordinate arithmetic

- Coordinates of the same coordinate type can be added and subtracted, and these operations are performed element-wise. Eltypes will be promoted to a common type as usual.

  ```jldoctest MapMaths
  julia> LatLon{Int}(1,2) - LatLon(2,1)
  LatLon{Float64}(-1.0, 1.0)
  ```

- Coordinates of different coordinate types generally cannot be added and subtracted.

  ```jldoctest MapMaths
  julia> Lon(1) - WMX(0.5)
  ERROR: MethodError: no method matching -(::Lon{Float64}, ::WMX{Float64})
  ```

- The one exception to this rule is that you can add and subtract `EastNorth` to / from any `Coordinate{2}`. The `EastNorth` coordinate is then interpreted as a translation of the other coordinate.

  ```jldoctest MapMaths
  julia> LatLon(1,2) + EastNorth(3,4)
  LatLon{Float64}(1.0000361746684363, 2.0000269535362025)
  ```

##### Utility functions

- Construct `X <: Coordinate{N}` from `Vararg{Number, N}` or `NTuple{N, Number}`.

  ```jldoctest MapMaths
  julia> LatLon(1,2)
  LatLon{Float64}(1.0, 2.0)

  julia> LatLon((1,2))
  LatLon{Float64}(1.0, 2.0)
  ```

- Convert to number / tuple of numbers using `getindex()`.

  ```jldoctest MapMaths
  julia> Lat(1)[]
  1.0

  julia> LatLon(1,2)[]
  (1.0, 2.0)
  ```

- Unpack using `iterate()`.

  ```jldoctest MapMaths
  julia> (lat,) = Lat(1); @show lat;
  lat = 1.0

  julia> (lat,lon) = LatLon(1,2); @show lat; @show lon;
  lat = 1.0
  lon = 2.0
  ```

  A common use-case of this function is to strip the coordinate type from function arguments.

  ```julia
  print_degrees((x,)::Union{Lat, Lon}) = print(x, "°")
  ```

- Convert `x::X{T1}` to `X{T2}` using `SX{T2}(x)` for any supertype `SX` of `X`.

  ```jldoctest MapMaths
  julia> x = WMX(1)
  WMX{Float64}(1.0)

  julia> WMX{Float32}(x)
  WMX{Float32}(1.0f0)

  julia> EastWestCoordinate{Float32}(x)
  WMX{Float32}(1.0f0)

  julia> Coordinate{1, Float32}(x)
  WMX{Float32}(1.0f0)
  ```
