"""
    East()

Meters east of the [`Origin`](@ref).

# Example
```jldoctest
julia> cmap(
           Origin(Wgs84(), LonLatAlt(), 0, 0, 0) + EastNorthUp(), 
           Coordinate(
               Origin(Wgs84(), LonLatAlt(), 90, 0, 0) + EastNorthUp(), 
               -Wgs84()[SemiMajorAxis], 0, 0,
           ),
       )
Coordinate(Origin(Wgs84(), GeodeticLonLatAlt(), 0.0, 0.0, 0.0) + EastNorthUp(), 6.378137e6, 0.0, 0.0)
"""
struct East <: ScalarCoordinateSystem end

"""
    North()

Meters north of the [`Origin`](@ref).

# Example
```jldoctest
julia> cmap(
           Origin(Wgs84(), LonLatAlt(), 0, 0, 0) + EastNorthUp(), 
           Coordinate(
               Origin(Wgs84(), LonLatAlt(), 0, 90, 0) + EastNorthUp(), 
               0, -Wgs84()[SemiMajorAxis], 0,
           ),
       )
Coordinate(Origin(Wgs84(), GeodeticLonLatAlt(), 0.0, 0.0, 0.0) + EastNorthUp(), 0.0, 6.356752314245179e6, 0.0)
```
"""
struct North <: ScalarCoordinateSystem end

"""
    Up()

Meters above the [`Origin`](@ref).

# Example
```jldoctest
julia> cmap(
           Alt(),
           Coordinate(
               Origin(Wgs84(), LonLatAlt(), 0, 0, 1) + EastNorthUp(), 
               0, 0, 2,
           ),
           Wgs84()
       )
Coordinate(GeodeticAltitude(), 3.0)
```
"""
struct Up <: ScalarCoordinateSystem end

"""
    Azimuth()

Angle in degrees relative to north. 

# Example
```jldoctest
julia> cmap(Azimuth(), Coordinate(EastNorthUp(), 1, 0, 0))
Coordinate(Azimuth(), 90.0)
```
"""
struct Azimuth <: ScalarCoordinateSystem end

"""
    Elevation()

Angle in degrees relative to the horizontal plane. 

# Example
```jldoctest
julia> cmap(Elevation(), Coordinate(EastNorthUp(), 0, 0, 1))
Coordinate(Elevation(), 90.0)
```
"""
struct Elevation <: ScalarCoordinateSystem end

"""
    HorizontalRange()

Distance in meters from [`Up`](@ref) axis. 

# Example
```jldoctest
julia> cmap(HorizontalRange(), Coordinate(EastNorthUp(), 1, 0, 1))
Coordinate(HorizontalRange(), 1.0)
```
"""
struct HorizontalRange <: ScalarCoordinateSystem end

"""
    Range()

Distance in meters from the [`Origin`](@ref). 

# Example
```jldoctest
julia> cmap(Range(), Coordinate(EastNorthUp(), 1, 0, 1))
Coordinate(Range(), 1.4142135623730951)
```
"""
struct Range <: ScalarCoordinateSystem end

"""
    EastNorth() = East() & North()
"""
@simple_vector_coordinate_system(EastNorth, East, North)

"""
    EastNorthUp() = East() & North() & Up()
"""
@simple_vector_coordinate_system(EastNorthUp, East, North, Up)

"""
    AzimuthElevation() = Azimuth() & Elevation()
"""
@simple_vector_coordinate_system(AzimuthElevation, Azimuth, Elevation)

"""
    AzimuthElevationRange() = Azimuth() & Elevation() & Range()
"""
@simple_vector_coordinate_system(AzimuthElevationRange, Azimuth, Elevation, Range)

###

lossy_parent(::Azimuth) = EastNorth()
lossy_parent(::Elevation) = HorizontalRange() & Up()
lossy_parent(::HorizontalRange) = EastNorth()
lossy_parent(::Range) = HorizontalRange() & Up()

@symmetric lossless_parent(::Azimuth, ::HorizontalRange) = EastNorth()
@symmetric lossless_parent(::Elevation, ::Range) = HorizontalRange() & Up()
@symmetric lossy_parent(::HorizontalRange, ::Up) = EastNorthUp()

@symmetric lossless_parent(::Azimuth, ::Elevation, ::Range) = Azimuth() & HorizontalRange() & Up()
@symmetric lossless_parent(::Azimuth, ::HorizontalRange, ::Up) = EastNorthUp()

###

cmap_impl(::Azimuth, (east, north)::Coordinate{EastNorth}) = atand(east, north)
cmap_impl(::HorizontalRange, (east, north)::Coordinate{EastNorth}) = hypot(east, north)
cmap_impl(::HorizontalRange, (el, r)::Coordinate(Elevation() & Range())) = r * cosd(el)
cmap_impl(::Up, (el, r)::Coordinate(Elevation() & Range())) = r * sind(el)
cmap_impl(::Range, (east, north, up)::Coordinate{EastNorthUp}) = hypot(east, north, up)
cmap_impl(::Range, (r, up)::Coordinate(HorizontalRange() & Up())) = hypot(r, up)
cmap_impl(::Elevation, (r, up)::Coordinate(HorizontalRange() & Up())) = atand(up, r)

cmap_impl(::EastNorth, (az, r)::Coordinate(Azimuth() & HorizontalRange())) = r .* sincosd(az)

function cmap_impl(::EastNorthUp, (az, el, r)::Coordinate(Azimuth() & Elevation() & Range()))
    sin_lon, cos_lon = sincosd(az)
    sin_lat, cos_lat = sincosd(el)
    return (
        r * cos_lat * sin_lon,
        r * cos_lat * cos_lon,
        r * sin_lat,
    )
end

###

"""
    Origin(datum, coord)
    Origin(datum, system, coords...)

Construct an origin for local-to-global coordinate conversions. 

Origins can be combined with a local coordinate or coordinate system using `+`.

# Example
```jldoctest
julia> cmap(Cartesian(), Coordinate(Origin(Wgs84(), LonLatAlt(), 0, 0, 0) + EastNorthUp(), 1, 0, 0))
Coordinate(Cartesian(), 6.378137e6, 1.0, 0.0)
```
"""
struct Origin{T, D}
    origin::SVector{3, T}
    ecef_from_enu::SMatrix{3, 3, T, 9}
    datum::D
end

function Origin(datum::Datum, coord::Coordinate)
    lon, lat, alt = cmap(LonLatAlt(), coord, datum)
    origin = SVector(cmap(Cartesian(), coord, datum)...)
    sin_lon, cos_lon = sincosd(lon)
    sin_lat, cos_lat = sincosd(lat)
    ecef_from_enu = @SMatrix [
        -sin_lon    -cos_lon * sin_lat    cos_lon * cos_lat
        +cos_lon    -sin_lon * sin_lat    sin_lon * cos_lat
        0                      cos_lat              sin_lat
    ]
    return Origin(origin, ecef_from_enu, datum)
end

Origin(datum::Datum, sys::CoordinateSystem, coords...) = Origin(datum, Coordinate(sys, coords))

function Base.show(io::IO, origin::Origin)
    print(io, Origin, "(", origin.datum, ", ", LonLatAlt(), ", ")
    join(io, cmap(LonLatAlt(), Coordinate(Cartesian(), origin.origin...), origin.datum), ", ")
    print(io, ")")
    return nothing
end

struct CoordinateSystemWithOrigin{S, T} <: CoordinateSystem
    system::S
    origin::Origin{T}
end

Base.iterate(sys::CoordinateSystemWithOrigin, state...) = iterate((sys.system, sys.origin), state...)
Base.show(io::IO, (sys, origin)::CoordinateSystemWithOrigin) = print(io, origin, " + ", sys)

Coordinate((system, origin)::CoordinateSystemWithOrigin, values...) = origin + Coordinate(system, values...)
Coordinate((system, origin)::CoordinateSystemWithOrigin, T::Type{<:Number}, values...) = origin + Coordinate(system, T, values...)

@symmetric Base.:+(origin::Origin, sys::CoordinateSystem) = CoordinateSystemWithOrigin(sys, origin)
@symmetric Base.:+(origin::Origin, coord::Coordinate) = Coordinate{:inner}(origin + coordinate_system(coord), coord.value)

lossy_parent((sys, origin)::CoordinateSystemWithOrigin) = origin + @something(lossy_parent(sys), return nothing)
lossless_parent((sys, origin)::CoordinateSystemWithOrigin) = origin + @something(lossless_parent(sys), return nothing)

lossless_parent(::CoordinateSystemWithOrigin{EastNorthUp}) = Cartesian()

for (maybe_datum, maybe_datum_arg) in (
        ((), ()),
        ((:datum,), (:(datum::Datum),)),
    )
    @eval begin
        function cmap(sys::CoordinateSystemWithOrigin, c::Coordinate{<:CoordinateSystemWithOrigin}, $(maybe_datum_arg...))
            sys_out, origin_out = sys
            sys_in, origin_in = coordinate_system(c)
            if origin_out == origin_in
                return origin_in + cmap(sys_out, Coordinate(sys_in, c...))
            else
                return cmap(sys, cmap(Cartesian(), c), $(maybe_datum...))
            end
        end
    end
end

cmap_impl((_, origin)::CoordinateSystemWithOrigin{EastNorthUp}, c::Coordinate{Cartesian}) = Tuple(origin.ecef_from_enu' * (SVector(c...) - origin.origin))
function cmap_impl(::Cartesian, c::Coordinate{<:CoordinateSystemWithOrigin{EastNorthUp}})
    _, origin = coordinate_system(c)
    return Tuple(origin.ecef_from_enu * SVector(c...) + origin.origin)
end
