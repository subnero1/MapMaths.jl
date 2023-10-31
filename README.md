# MapMaths.jl

Utility tools for working with various coordinate systems.

## Examples

Convert latitude and longitude to WebMercator coordinates.

```julia
julia> WebMercator(LatLon(1.3, 103.8))
WebMercator{Float64}(0.5766666666666667, 0.007222841972793043)
```

How far from the equator is Singapore?

```julia
julia> North(LatLon(1.3, 103.8))
North{Float64}(143746.80623936345)
```

Move 1 km east and 0.5 km north.

```julia
julia> LatLon(1.3, 103.8) + EastNorth(1000, 500)
LatLon{Float64}(1.3045218239325143, 103.80898545013565)
```

## MapMaths vs Geodesy

MapMaths serves a purpose which is very similar to that of [Geodesy.jl](https://github.com/JuliaGeo/Geodesy.jl). We hope that at some point in the future the two packages could merge into a single "coordinate maths for Julia" package, but in the meantime the two packages differ in various aspects, and you may prefer to use one or the other depending on your application.

- MapMaths currently supports only the WGS84 datum.

- Geodesy is much older, more stable and has a much larger number of contributors.

- MapMaths has first-class support for WebMercator coordinates, i.e. working with WebMercator coordinates is as convenient as working with any other coordinate type when using MapMaths.
  ```julia
  julia> using Geodesy

  julia> ECEF(LLA(0,0), wgs84) # <- Convenient
  ECEF(6.378137e6, 0.0, 0.0)

  julia> ECEF(LLAfromWebMercator(wgs84)([0,0,0]), wgs84) # <- Clumsy
  ECEF(6.378137e6, 0.0, 0.0)
  ```
  ```julia
    julia> using MapMaths

    julia> ECEF(LatLon(0,0)) # <- Convenient
    ECEF{Float64}(6.378137e6, 0.0, 0.0)

    julia> ECEF(WebMercator(0,0)) # <- Also convenient
    ECEF{Float64}(6.378137e6, 0.0, 0.0)
  ```

- MapMaths defaults to floating-point numbers for better type-stability, and eltype conversions "just work".
  ```julia
  julia> using MapMaths

  julia> LatLon(0,0)
  LatLon{Float64}(0.0, 0.0)

  julia> LatLon{Float64}[ LatLon{Int}(0,0) ]
  1-element Vector{LatLon{Float64}}:
   LatLon{Float64}(0.0, 0.0)
  ```
  ```julia
  julia> using Geodesy

  julia> LatLon{Float64}[ LatLon(0,0) ]
  ERROR: MethodError: Cannot `convert` an object of type
    LatLon{Int64} to an object of type
    LatLon{Float64}
  ```

- MapMaths allows you do to a large number of coordinate conversion by basically just throwing all the required information at it.
  ```julia
  julia> using MapMaths

  julia> ECEF(LatLon(0, 90))  # Altitude defaults to 0...
  ECEF{Float64}(0.0, 6.378137e6, 0.0)

  julia> ECEF(LatLon(0, 90), Alt(1))  # ... but can be provided explicitly when desired.
  ECEF{Float64}(0.0, 6.378138e6, 0.0)

  julia> ECEF(LonLat(90, 0))  # We can also swap latitude and longitude...
  ECEF{Float64}(0.0, 6.378137e6, 0.0)

  julia> ECEF(Lon(90), Lat(0))  # ... or provide them separately.
  ECEF{Float64}(0.0, 6.378137e6, 0.0)
  ```

  Doing similar conversions in Geodesy can be quite cumbersome.

  ```julia
  julia> using Geodesy

  julia> c = LatLon(0, 90)  # Location provided from somewhere else
         ECEF(LLA(c.lat, c.lon, 0), wgs84)  # Cumbersome, and only works for c::LatLon
  ECEF(0.0, 6.378137e6, 0.0)
  ```
