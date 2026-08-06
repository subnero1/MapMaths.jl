module MapMaths

export
    CoordinateSystem,
    Coordinate,
    coordinate_system,
    Datum,
    cmap,
    SemiMajorAxis,
    SemiMinorAxis,
    Wgs84,
    CartesianX,
    CartesianY,
    CartesianZ,
    CartesianXY,
    Cartesian,
    Longitude,
    CylindricalRadius,
    Radius,
    GeocentricLatitude,
    GeocentricAltitude,
    Cylindrical,
    Spherical,
    GeocentricLonLat,
    GeocentricLatLon,
    GeocentricLonLatAlt,
    GeocentricLatLonAlt,
    Lon,
    GeocentricLat,
    GeocentricAlt,
    GeodeticLatitude,
    GeodeticAltitude,
    GeodeticLonLat,
    GeodeticLatLon,
    GeodeticLatAlt,
    GeodeticLonLatAlt,
    GeodeticLatLonAlt,
    GeodeticLat,
    GeodeticAlt,
    Lat,
    Alt,
    LonLat,
    LatLon,
    LatAlt,
    LonLatAlt,
    LatLonAlt,
    WmX,
    WmY,
    WebMercator,
    WebMercatorAlt,
    East,
    North,
    Up,
    Azimuth,
    Elevation,
    HorizontalRange,
    Range,
    EastNorth,
    EastNorthUp,
    AzimuthElevation,
    AzimuthElevationRange,
    Origin

using .Meta: isexpr, unescape
using Base.ScopedValues
using Combinatorics
using LinearAlgebra
using StaticArrays

include("utils.jl")
include("ellipses.jl")

# The following functions are faster but less accurate than their Base
# counterparts. Uncomment them for fair comparisons against packages which use
# these fast definitions (e.g. CoordRefSystems.jl).
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# sind(x) = sin(deg2rad(x))
# cosd(x) = cos(deg2rad(x))
# sincosd(x) = sincos(deg2rad(x))
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


#####################################
# Coordinates and coordinate systems

"""
    CoordinateSystem

Supertype for coordinate systems (e.g. `Cartesian()` or `LonLat()`).

Coordinate systems can be combined using `&`.

# Example
```jldoctest
julia> Lon() & Lat()
GeodeticLonLat()
```
"""
abstract type CoordinateSystem end

"""
    Coordinate{S <: CoordinateSystem}

Generic coordinate type. 

Notable operations: 

- `Coordinate` constructors are generally of the form `Coordinate(system,
  args...)`. Check the docstring of the coordinate system for the precise
  signature. 
- Like `CoordinateSystem`s, `Coordinate`s can be combined using `&`. 
- Iterating a `Coordinate` yields its components. 
- `Coordinate(s::CoordinateSystem) = Coordinate{typeof(s))}` is a convenience
  method for constructing `Coordinate` type signatures for dispatch. 

# Example
```jldoctest
julia> Coordinate(Lon(), 1) & Coordinate(Lat(), 2)
Coordinate(GeodeticLonLat(), 1.0, 2.0)

julia> (Coordinate(LonLat(), 1,2)...,)
(1.0, 2.0)
```
"""
struct Coordinate{S, V}
    system::S
    value::V

    # Dummy inner constructor to avoid ambiguities with the default outer constructor
    Coordinate{:inner}(system::CoordinateSystem, value::Tuple) = new{typeof(system), typeof(value)}(system, value)
end

Coordinate(system::CoordinateSystem) = Coordinate{typeof(system)}
Base.iterate(coord::Coordinate, state...) = iterate(coord.value, state...)
Base.length(coord::Coordinate) = length(coord.value)
Base.eltype(coord::Coordinate) = eltype(coord.value)
function Base.show(io::IO, coord::Coordinate)
    print(io, Coordinate, "(", coord.system, ", ")
    join(io, coord.value, ", ")
    print(io, ")")
    return
end

"""
    coordinate_system(coord)

Coordinate system of the given coordinate. Also works for types, in which case
it returns the type of the coordinate system.
"""
function coordinate_system end
coordinate_system(coord::Coordinate) = coord.system
coordinate_system(::Type{<:Coordinate{S}}) where {S} = S


###################
# Conversion graph

"""
    lossy_parent(system)

Coordinate system to use when converting to the given coordinate system.

Defining `lossy_parent(::Child) = parent` is similar to adding a method
```
cmap(system::Child, coord) = cmap(system, cmap(parent, coord))
```
but it avoids the ambiguities that the above method would create.

See also [`lossless_parent`](@ref).
"""
function lossy_parent end

"""
    lossless_parent(system)

Coordinate system to use when converting to or from the given coordinate system.

Defining `lossless_parent(::Child) = parent` is similar to adding the methods
```
cmap(system::Child, coord) = cmap(system, cmap(parent, coord))
cmap(system, coord::Coordinate{Child}) = cmap(system, cmap(parent, coord))
```
but it avoids the ambiguities that the above methods would create.

This function can also be called with a `Coordinate`, in which case it is
equivalent to `lossless_parent(coordinate_system(coord))`.

See also [`lossy_parent`](@ref).
"""
function lossless_parent end

lossy_parent(system::CoordinateSystem...) = lossless_parent(system...)
lossless_parent(system::CoordinateSystem...) = nothing
lossless_parent(coord::Coordinate) = lossless_parent(coordinate_system(coord))


####################
# Tuple coordinates

"""
    TupleCoordinateSystem <: CoordinateSystem

Combination of multiple coordinate systems (e.g. `LonLat`). See [`&`](@ref).
"""
struct TupleCoordinateSystem{P <: NTuple{<:Any, CoordinateSystem}} <: CoordinateSystem
    parts::P
end
TupleCoordinateSystem(parts::CoordinateSystem...) = TupleCoordinateSystem(parts)

"""
    parts(coord_or_system)

Decompose a coordinate or coordinate system into its parts. Also works for types. 

# Example
```jldoctest
julia> parts(LonLat)
(Longitude, GeodeticLatitude)
```
"""
function parts end
parts(system::CoordinateSystem) = (system,)
parts(System::Type{<:CoordinateSystem}) = (System,)
parts(coord::Coordinate) = (coord,)
parts(Coord::Type{<:Coordinate}) = (Coord,)
parts(tuple_system::TupleCoordinateSystem) = tuple_system.parts
parts(TupleSystem::Type{<:TupleCoordinateSystem}) = fieldtypes(fieldtype(TupleSystem, :parts))
parts(tuple_coord::Coordinate{<:TupleCoordinateSystem}) = map(
    (s, v) -> Coordinate{:inner}(s, (v,)),
    parts(coordinate_system(tuple_coord)),
    tuple_coord.value,
)
function parts(TupleCoord::Type{<:Coordinate{<:TupleCoordinateSystem}})
    return map(
        (S, V) -> Coordinate{S, Tuple{V}},
        parts(fieldtype(TupleCoord, :system)),
        fieldtypes(fieldtype(TupleCoord, :value)),
    )
end

"""
    n_parts(system_or_coord)

Number of parts in the given coordinate or coordinate system.
"""
n_parts(coord_or_system) = length(parts(coord_or_system))

"""
    (&)(coords...)
    (&)(systems...)

Combine coordinates or coordinate systems.  

See also [`join_coords`](@ref) and [`join_systems`](@ref) 
"""
Base.:&(systems::CoordinateSystem...) = TupleCoordinateSystem(flatten(parts.(systems)))
Base.:&(system::CoordinateSystem) = system

function Base.:&(coords::Coordinate...)
    return Coordinate{:inner}(
        (&)(coordinate_system.(coords)...),
        flatten(map(coord -> coord.value, coords))
    )
end
Base.:&(coord::Coordinate) = coord

lossy_parent(system::TupleCoordinateSystem) = lossy_parent(parts(system)...)
lossless_parent(system::TupleCoordinateSystem) = lossless_parent(parts(system)...)

function Base.show(io::IO, system::TupleCoordinateSystem)
    if length(parts(system)) == 0
        return print(io, join_systems, "()")
    end
    return join(io, parts(system), " & ")
end


########
# Datum

"""
    Datum

Supertype for coordinate datums (e.g. `Wgs84()`).

Notable operations:

- `Ellipse(datum)` returns the reference ellipse for the given datum.
- `datum[EllipseParameter]` is a convenience method for
  `Ellipse(datum)[EllipseParameter]`.
"""
abstract type Datum end
Base.getindex(datum::Datum, P::Type{<:EllipseParameter}) = Ellipse(datum)[P]

"""
    NoDatum <: Datum

Placeholder indicating that no datum has been provided to a conversion function.
"""
struct NoDatum <: Datum end
Ellipse(::NoDatum) = error("Conversion requires a datum.")

#####################
# Scalar coordinates

"""
    ScalarCoordinateSystem <: CoordinateSystem

Supertype for coordinate systems with a single scalar value (e.g. `Longitude` or
`CartesianX`).

Coordinates in a scalar coordinate system `system` can be constructed using
`Coordinate(system, [T], x::Number)` where `T` is an optional `Number` type. 
Similarly, Coordinates in a combinations of scalar coordinate systems can be
constructed using `Coordinate(system, [T], x::Number...)`.
"""
abstract type ScalarCoordinateSystem <: CoordinateSystem end

Coordinate(system::ScalarCoordinateSystem, value::Number) = Coordinate(system, default_numtype(value), value)
Coordinate(system::ScalarCoordinateSystem, T::Type{<:Number}, value::Number) = Coordinate{:inner}(system, (convert(T, value),))

numtype(Coord::Type{<:Coordinate{<:ScalarCoordinateSystem}}) = numtype(fieldtype(Coord, :value))
with_numtype(c::Coordinate{<:ScalarCoordinateSystem}, T::Type{<:Number}) = Coordinate(c.system, T, only(c.value))


#####################
# Vector coordinates

const VectorCoordinateSystem{N} = TupleCoordinateSystem{<:NTuple{N, ScalarCoordinateSystem}}

Coordinate(system::VectorCoordinateSystem{N}, values::Vararg{Number, N}) where {N} = Coordinate(system, values)
Coordinate(system::VectorCoordinateSystem{N}, values::NTuple{N, Number}) where {N} = Coordinate(system, default_numtype(values), values)
Coordinate(system::VectorCoordinateSystem{N}, T::Type{<:Number}, values::Vararg{Number, N}) where {N} = Coordinate(system, T, values)
Coordinate(system::VectorCoordinateSystem{N}, T::Type{<:Number}, values::NTuple{N, Number}) where {N} = Coordinate{:inner}(system, convert.(T, values))

numtype(Coord::Type{<:Coordinate{<:VectorCoordinateSystem}}) = numtype(fieldtype(Coord, :value))
with_numtype(c::Coordinate{<:VectorCoordinateSystem}, T::Type{<:Number}) = Coordinate(c.system, T, c.value)


########################
# Coordinate conversion

"""
    cmap(system, coord, [datum])

Map a coordinate to the given coordinate system, using `datum` if provided (and
needed).

The behaviour of `cmap()` can be extended by adding methods to
[`lossy_parent()`](@ref), [`lossless_parent()`](@ref), and
[`cmap_impl()`](@ref). See the "Under The Hood" section in the
[README](README.md#under-the-hood) for details.

# Example
```jldoctest
julia> cmap(Longitude(), Coordinate(CartesianXY(), 0, 1))
Coordinate(Longitude(), 90.0)
```
"""
function cmap(system::CoordinateSystem, coord::Coordinate, datum = NoDatum())
    return @something(
        permute_coords(system, coord),
        apply_cmap_impl(system, coord, datum),
        bump_system(system, coord, datum),
        bump_coord(system, coord, datum),
        error("Cannot convert from $(coordinate_system(coord)) to $system. Please raise an issue if you think this is a bug.")
    )
end

"""
    cmap_impl(system, coord, datum)

Coordinate conversion rules to be used by [`cmap()`](@ref).
"""
function cmap_impl end

"""
    permute_coords(system, coord)

Check if `system` is a permutation of `coordinate_system(coord)`, and if so,
apply that permutation to `coord`. Returns `nothing` otherwise.
"""
@generated function permute_coords(system::CoordinateSystem, coord::Coordinate)
    result = Expr(:block)
    for idx in permutations(1:n_parts(coord), n_parts(system))
        push!(
            result.args,
            quote
                let permuted_coord = (&)($(map(i -> :(parts(coord)[$i]), idx)...))
                    if system == coordinate_system(permuted_coord)
                        return permuted_coord
                    end
                end
            end
        )
    end
    push!(
        result.args,
        :(return nothing)
    )
    return result
end

"""
    apply_cmap_impl(system, coord, datum)

Check whether there exists a `cmap_impl()` method for any subset of
`system` and `coord` components. If so, apply it and return the
result, otherwise return `nothing`.
"""
@generated function apply_cmap_impl(system::CoordinateSystem, coord::Coordinate, datum::Datum)
    result = Expr(:block)
    for n_systems in reverse(1:n_parts(system)),
            n_coords in 1:n_parts(coord),
            system_idx in permutations(1:n_parts(system), n_systems),
            coord_idx in permutations(1:n_parts(coord), n_coords)
        remainder = if n_systems < n_parts(system)
            (
                :(
                    cmap(
                        (&)($(map(i -> :(parts(system)[$i]), setdiff(1:n_parts(system), system_idx))...)),
                        coord,
                        datum,
                    )
                ),
            )
        else
            ()
        end
        push!(
            result.args,
            quote
                let permuted_system = (&)($(map(i -> :(parts(system)[$i]), system_idx)...)),
                        permuted_coord = (&)($(map(i -> :(parts(coord)[$i]), coord_idx)...))
                    if Core._hasmethod( # https://discourse.julialang.org/t/method-lookup-in-generated-functions/138047/6
                            Tuple{
                                typeof(cmap_impl),
                                typeof(permuted_system),
                                typeof(permuted_coord),
                                typeof(datum),
                            }
                        )
                        return cmap(
                            system,
                            (&)(
                                Coordinate(
                                    permuted_system,
                                    cmap_impl(permuted_system, permuted_coord, datum),
                                ),
                                $(remainder...),
                            ),
                        )
                    end
                end
            end
        )
    end
    push!(
        result.args,
        :(return nothing)
    )
    return result
end

"""
    bump_system(system, coord, datum)

Route the conversion via `lossy_parent(system)`, if possible. Returns `nothing`
otherwise. 
"""
function bump_system(system::CoordinateSystem, coord::Coordinate, datum)
    if (
            !isnothing(lossy_parent(system))
                && !is_ancestor(lossless_parent, system, coordinate_system(coord))
                && !is_set_equal(parts(lossy_parent(system)), parts(coordinate_system(coord)))
        )
        return cmap(
            system,
            cmap(
                lossy_parent(system),
                coord,
                datum,
            ),
            datum,
        )
    end
    return nothing
end

"""
    bump_coord(system, coord, datum)

Route the conversion via `lossless_parent(coord)`, if possible. Returns `nothing`
otherwise. 
"""
function bump_coord(system::CoordinateSystem, coord::Coordinate, datum)
    if (
            !isnothing(lossless_parent(coord))
                && !is_ancestor(lossy_parent, coordinate_system(coord), system)
                && !is_set_equal(parts(system), parts(lossless_parent(coord)))
        )
        return cmap(
            system,
            cmap(
                lossless_parent(coord),
                coord,
                datum
            ),
            datum
        )
    end
    return nothing
end

"""
    is_ancestor(parent, sys1, sys2)

Check whether `sys1` is an ancestor of `sys2` according to the `parent` relation. 

See also [`lossy_parent`](@ref) and [`lossless_parent`](@ref).
"""
function is_ancestor(parent, sys1::CoordinateSystem, sys2::CoordinateSystem)
    if issubset(parts(sys2), parts(sys1))
        return true
    end
    if isnothing(parent(sys2))
        return false
    end
    return is_ancestor(parent, sys1, parent(sys2))
end
is_ancestor(system::CoordinateSystem, coord::Coordinate) = is_ancestor(system, coordinate_system(coord))

################
# Miscellaneous

macro simple_vector_coordinate_system(Vec, Parts...)
    Vec = esc(Vec)
    Parts = esc.(Parts)
    return quote
        Core.@__doc__ const $Vec = typeof((&)($(map(Part -> :($Part()), Parts)...)))
        $Vec() = (&)($(map(Part -> :($Part()), Parts)...))
        function Base.show(io::IO, ::$Vec)
            @print_var(io, $Vec)
            return print(io, "()")
        end
    end
end


###########################
# Concrete implementations

include("data.jl")

struct CartesianX <: ScalarCoordinateSystem end
struct CartesianY <: ScalarCoordinateSystem end
struct CartesianZ <: ScalarCoordinateSystem end

@simple_vector_coordinate_system(CartesianXY, CartesianX, CartesianY)
@simple_vector_coordinate_system(Cartesian, CartesianX, CartesianY, CartesianZ)

include("geocentric.jl")
include("geodetic.jl")
include("webmercator.jl")
include("local.jl")

end
