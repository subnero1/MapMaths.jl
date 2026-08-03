# MapMaths.jl

Fast and flexible coordinate conversions for Julia.

**WARNING:** This package is currently in an experimental phase. See [below](#history-outlook--alternatives) for details. 


## How Does It Work?

1. Assemble a coordinate system. You may do so using the `&` operator...
    ```julia
    julia> Lat() & Lon() & Alt()
    GeodeticLatLonAlt()
    ```
    ... but more commonly you would use one of the built-in convenience constructors. 
    ```julia
    julia> LatLonAlt()
    GeodeticLatLonAlt()
    ```

2. Create a coordinate.
    ```julia
    julia> lla_coords = Coordinate(LatLonAlt(), 1.286770, 103.854307, 2.0)
    Coordinate(GeodeticLatLonAlt(), 1.28677, 103.854307, 2.0)
    ```

3. Convert.
    ```julia
    julia> cartesian_coords = cmap(Cartesian(), lla_coords, Wgs84())
    Coordinate(Cartesian(), -1.5268872205890939e6, 6.19103341816399e6, 142271.98546019828)
    ```

4. Extract.
    ```julia
    julia> Tuple(cartesian_coords)
    (-1.5268872205890939e6, 6.19103341816399e6, 142271.98546019828)
    ```
    Extraction is based on iterating the coordinate object and can also be used in other ways. 
    ```julia
    julia> x,y,z = cartesian_coords; (x,y,z)
    (-1.5268872205890939e6, 6.19103341816399e6, 142271.98546019828)

    julia> SVector(cartesian_coords...)
    3-element SVector{3, Float64} with indices SOneTo(3):
        -1.5268872205890939e6
          6.19103341816399e6
    142271.98546019828
    ```

## Coordinate Systems

### Earth-Centered, Earth-Fixed Coordinate Systems

- Cartesian: `Cartesian() = CartesianX() & CartesianY() & CartesianZ()`
- Cylindrical: `Cylindrical() = Longitude() & CylindricalRadius() & CartesianZ()`
- Spherical: `Spherical() = Longitude() & GeocentricLatitude() & Radius()`

- Geocentric LatLonAlt: `GeocentricLonLatAlt() = Longitude() & GeocentricLatitude() & GeocentricAltitude()`
  <br>
  Analogously: `GeocentricLatLonAlt()`, `GeocentricLonLat()`, `GeocentricLatLon()`, `GeocentricLatAlt()`

  `GeocentricAltitude()` is defined as the distance to the reference ellipsoid along the line to the centre of the earth. This is conceptually different from the `GeodeticAltitude()`, but in practice the difference is tiny. 

  ```julia
  julia> cmap(GeodeticAltitude(), Coordinate(GeocentricLatAlt(), 45, 100), Wgs84())
  Coordinate(GeodeticAltitude(), 99.99943604320288)
  ```

- Geodetic LatLonAlt: `GeodeticLonLatAlt() = Longitude() & GeodeticLatitude() & GeodeticAltitude()`
  <br>
  Analogously: `GeodeticLatLonAlt()`, `GeodeticLonLat()`, `GeodeticLatLon()`, `GeodeticLatAlt()`

  Abbreviations: `Lon() = Longitude()`, `Lat() = GeodeticLatitude()`, `Alt() = GeodeticAltitude()`, and all the above combinations without the `Geodetic` prefix (e.g. `LonLat()`). We recommend to omit the `Geodetic` prefix in codes that use only the geodetic coordinate systems but to add it for extra clarity in codes which use both geodetic and geocentric coordinate systems. 

- WebMercator: `WebMercatorAlt() = WmX() & WmY() & GeodeticAltitude()`
  <br>
  Analogously: `WebMercator() = WmX() & WmY()`. 

  The WebMercator coordinates are scaled so that `Coordinate(Lon(), 180)` is mapped to `Coordinate(WmX(), 1.0)`. 

### Local Coordinate Systems

- East, North, Up: `EastNorthUp() = East() & North() & Up()`
- Azimuth, Elevation, Range: `AzimuthElevationRange() = Azimuth() & Elevation() & Range()`. 

Local coordinate systems must be combined with an `Origin()` to be convertible to and from global coordinate systems. 

```julia
julia> # Local to local conversion - no Origin needed
       cmap(EastNorthUp(), Coordinate(AzimuthElevationRange(), 0, 0, 1))
Coordinate(EastNorthUp(), 1.0, 0.0, 0.0)

julia> # Local to global conversion - must provide an Origin
       coord = Origin(Wgs84(), LonLatAlt(), 0, 0, 0) + Coordinate(EastNorthUp(), 1, 0, 0)
       cmap(Cartesian(), coord)
Coordinate(Cartesian(), 6.378137e6, 1.0, 0.0)
```


## Data

Some coordinate transformations require a geodetic datum (a mathematical model of the Earth). MapMaths currently implements only one datum, namely `Wgs84()`. See [`src/data.jl`](src/data.jl) for how to add your own data. 

MapMaths currently does not support datum conversions. 


## Notes

Just like coordinate systems, you can also combine coordinates using the `&` operator. 

```julia
julia> Coordinate(LonLat(), 1, 2) & Coordinate(Alt(), 3)
Coordinate(GeodeticLonLatAlt(), 1.0, 2.0, 3.0)
```

Conversion to and from "standard" combinations work just fine, even if their order is not standard. 

```julia
julia> cmap(WmY() & WmX(), Coordinate(CartesianX() & CartesianZ() & CartesianY(), 0, 0, 1), Wgs84())
Coordinate(WmY() & WmX(), 0.0, 0.5)
```

However, conversions to and from "unusual" combinations may not always work. 

```julia
julia> cmap(Lat() & WmX(), Coordinate(CartesianX() & CartesianZ() & CartesianY(), 0, 0, 1), Wgs84())
ERROR: Cannot convert from CartesianX() & CartesianZ() & CartesianY() to GeodeticLatitude() & WmX(). Please raise an issue if you think this is a bug.
```


## History, Outlook & Alternatives

The main goal of MapMaths v0.1 was to sand away some of the rough edges of [Geodesy.jl](https://github.com/JuliaGeo/Geodesy.jl), but this goal has since been achieved with much better feature-completeness in [CoordRefSystems.jl](https://github.com/JuliaEarth/CoordRefSystems.jl). Correspondingly, the current version v0.2 instead focuses on prototyping a framework that makes it economical to support a much larger share of partial coordinate conversions (e.g. `GeocentricLatAlt()` to `GeodeticLat()`) while guaranteeing transitive closure (if you can convert `A` to `B` and `B` to `C`, then you can also directly convert `A` to `C`). This framework is currently still experimental. We hope that we can join forces with the Geodesy and CoordRefSystems communities if and when it stabilises. 


## Under The Hood

MapMaths' conversion mechanism is based on two primitives:

- `cmap_impl()` provides explicit conversion rules between concrete start- and endpoints (e.g. `CartesianXY()` to `Longitude()`). 

- `lossless_parent()` and `lossy_parent()` provide fallback coordinate systems via which to route the conversion in case there is no direct conversion rule (e.g. `lossless_parent(::WmX) = Longitude()` or `lossy_parent(::Longitude()) = CartesianXY()`). 

Conversion from `A` to `B` then proceeds as follows.

1. Check if there is a `cmap_impl()` method mapping any part of `A` to any part of `B`. If so, apply this method and then recursively convert the result to `B`. 

    If there are multiple applicable `cmap_impl()` methods, then `cmap()` will favor methods producing more parts of `B`. 

2. Otherwise, try converting `B` to `lossy_parent(A)` and then convert the result to `A`. Skip this step if `B` is already an ancestor of `A` according to the `lossless_parent()` relationship. 

3. If both of the above steps failed, then try converting `B` to `lossless_parent(B)` and then convert the result to `A`. 

4. If all of these steps fail, throw an error. 