module MapMaths

using Elliptic
using StaticArrays
using LinearAlgebra

export
    Coordinate, EastWestCoordinate, NorthSouthCoordinate,
    WebMercator, WMX, WMY,
    LatLon, LonLat, Lat, Lon,
    EastNorth, East, North,
    ECEF, Alt, tdist

abstract type Coordinate{N, T <: Number} end
(::Type{C})(v::NTuple{N, Number}) where {N, C <: Coordinate{N}} = C(v...)
(::Type{C})(v::AbstractVector{<:Number}) where {C <: Coordinate} = C(v...)

abstract type EastWestCoordinate{T <: Number} <: Coordinate{1, T} end
abstract type NorthSouthCoordinate{T <: Number} <: Coordinate{1, T} end

for (C,D) in (
    (:WMX, :EastWestCoordinate), (:WMY, :NorthSouthCoordinate),
    (:Lon, :EastWestCoordinate), (:Lat, :NorthSouthCoordinate),
    (:East, :EastWestCoordinate), (:North, :NorthSouthCoordinate),
)
    @eval begin
        struct $C{T <: Number} <: $D{T}
            v::T
            $C{T}(v::Number) where {T <: Number} = new(v) # Remove default constructors
        end

        # Eltype conversion
        Coordinate{1,T}(c::$C) where {T <: Number} = $C{T}(c.v)
        $D{T}(c::$C) where {T <: Number} = $C{T}(c.v)
        $C{T}(c::$C) where {T <: Number} = $C{T}(c.v)

        # Construct from number
        $C(v::Number) = $C{float(typeof(v))}(v) # Convert ints to floats by default

        # Construct from single coordinate
        $C(c::$D{T}) where {T <: Number} = $C{T}(c) # Preserve eltype when converting
        $C{T}(c::$D) where {T <: Number} = $C{T}($D{T}(c)) # Match eltypes
        $C{T}(c::$D{T}) where {T <: Number} = $C{T}(_convert($C, c)) # Dispatch to _convert

        # Construct from coordinate pair
        $C(c::$D{T}, o::Coordinate{1}) where {T <: Number} = $C{T}(c, o) # Preserve eltype when converting
        $C{T}(c::$D, o::Coordinate{1}) where {T <: Number} = $C{T}($D{T}(c), Coordinate{1,T}(o)) # Match eltypes
        $C{T}(c::$D{T}, o::Coordinate{1,T}) where {T <: Number} = $C{T}(_convert($C, c, o)) # Dispatch to _convert
    end
end

for (C,(C1,C2)) in (
    (:WebMercator, (:WMX, :WMY)),
    (:LatLon, (:Lat, :Lon)),
    (:LonLat, (:Lon, :Lat)),
    (:EastNorth, (:East, :North)),
)
    @eval begin
        struct $C{T <: Number} <: Coordinate{2, T}
            c1::$C1{T}
            c2::$C2{T}
            $C{T}(c1::$C1{T}, c2::$C2{T}) where {T <: Number} = new(c1, c2) # Remove default constructors
        end

        # Eltype conversion
        Coordinate{2,T}(c::$C) where {T <: Number} = $C{T}(c...)
        $C{T}(c::$C) where {T <: Number} = $C{T}(c...)

        # Turn numbers into coordinates
        $C(v1::Number, v2::Number) = $C{float(promote_type(typeof(v1),typeof(v2)))}(v1, v2) # Convert ints to floats by default
        $C{T}(v1::Number, v2::Number) where {T <: Number} = $C{T}($C1{T}(v1), $C2{T}(v2)) # Turn numbers into coordinates
        $C(x1::Union{Number, Coordinate{1}}, x2::Union{Number, Coordinate{1}}) = $C{promote_type(eltype(x1),eltype(x2))}(x1, x2) # Infer eltype if at least one argument is a coordinate

        # Split two-dimensional coordinates
        $C(c::Coordinate{2,T}) where {T <: Number} = $C{T}(c) # Preserve eltype when converting
        $C{T}(c::Coordinate{2}) where {T <: Number} = $C{T}(c.c1, c.c2) # Dispatch to single-argument conversion

        # Standardize coordinate pairs
        $C{T}(c2::supertype($C2), c1::supertype($C1)) where {T <: Number} = $C{T}(c1, c2) # Swap arguments if needed
        $C{T}(c1::supertype($C1), c2::supertype($C2)) where {T <: Number} = $C{T}(supertype($C1){T}(c1), supertype($C2{T})(c2)) # Convert to coordinate type

        # Dispatch to convert
        function $C{T}(c1::supertype($C1){T}, c2::supertype($C2){T}) where {T <: Number}
            return $C{T}(
                _convert($C1, c1, c2),
                _convert($C2, c2, c1),
            )
        end

        # Extract one-dimensional coordinates
        (::Type{D1})(c::$C) where {D1 <: supertype($C1)} = D1(c.c1, c.c2)
        (::Type{D2})(c::$C) where {D2 <: supertype($C2)} = D2(c.c2, c.c1)
    end
end

struct Alt{T <: Number} <: Coordinate{1,T}
    v::T
end
Alt{T}(alt::Alt) where {T <: Number} = Alt{T}(alt.v)

struct ECEF{T <: Number} <: Coordinate{3,T}
    v1::T
    v2::T
    v3::T

    ECEF{T}(v1::Number, v2::Number, v3::Number) where {T <: Number} = new(v1, v2, v3)
end

# Eltype conversion
Coordinate{3,T}(c::ECEF) where {T <: Number} = ECEF{T}(c...)
ECEF{T}(c::ECEF) where {T <: Number} = ECEF{T}(c...)

# Construct from number
ECEF(v1::Number, v2::Number, v3::Number) = ECEF{float(promote_type(typeof.((v1,v2,v3))...))}(v1,v2,v3) # Convert ints to floats by default

# Construct from coordinates
ECEF(c1::Coordinate{1}, c2::Coordinate{1}) = ECEF{promote_type(eltype(c1), eltype(c2))}(c1, c2)
ECEF(c1::Coordinate{1}, c2::Coordinate{1}, alt::Alt) = ECEF{promote_type(eltype(c1), eltype(c2), eltype(alt))}(c1, c2, alt)
ECEF(c::Coordinate{2}) = ECEF{eltype(c)}(c)
ECEF(c::Coordinate{2}, alt::Alt) = ECEF{promote_type(eltype(c), eltype(alt))}(c, alt)

# Allow arguments to be passed as tuples (required e.g. for tdsit)
ECEF(c::Tuple{Coordinate{1}, Coordinate{1}}) = ECEF(c...)
ECEF(c::Tuple{Coordinate{1}, Coordinate{1}, Alt}) = ECEF(c...)
ECEF(c::Tuple{Coordinate{2}, Alt}) = ECEF(c...)

# Promote eltypes and convert to LonLat
ECEF{T}(c1::Coordinate{1}, c2::Coordinate{1}) where {T <: Number} = ECEF{T}(c1, c2, Alt{T}(0))
ECEF{T}(c1::Coordinate{1}, c2::Coordinate{1}, alt::Alt) where {T <: Number} = ECEF{T}(LonLat{T}(c1, c2), Alt{T}(alt))
ECEF{T}(c::Coordinate{2}) where {T <: Number} = ECEF{T}(LonLat{T}(c), Alt{T}(0))
ECEF{T}(c::Coordinate{2}, alt::Alt) where {T <: Number} = ECEF{T}(LonLat{T}(c), Alt{T}(alt))

# Construct from LonLat
function ECEF{T}((lon, lat)::LonLat{T}, (alt,)::Alt{T}) where {T <: Number}
    sind_lon, cosd_lon = sincosd(lon)
    sind_lat, cosd_lat = sincosd(lat)
    ECEF{T}(
        (N(lat) + alt) * cosd_lat * cosd_lon,
        (N(lat) + alt) * cosd_lat * sind_lon,
        ((1-f)^2 * N(lat) + alt) * sind_lat,
    )
end

const EcefArguments = Union{
    Coordinate{2},
    Tuple{Coordinate{2}, Alt},
    Tuple{Coordinate{1}, Coordinate{1}},
    Tuple{Coordinate{1}, Coordinate{1}, Alt},
}

tdist(c1::EcefArguments, c2::EcefArguments) = tdist(ECEF(c1), ECEF(c2))
function tdist(c1::ECEF, c2::ECEF)
    T = promote_type(eltype(c1), eltype(c2))
    tdist(ECEF{T}(c1), ECEF{T}(c2))
end
tdist(c1::ECEF{T}, c2::ECEF{T}) where {T <: Number} = norm(SVector(c1) - SVector(c2))

include("docstrings.jl")

"""
    _convert(To::Type{<:Coordinate}, from::Coordinate, other...)-> to::Number

Convert coordinate `from` to type `To`. This function is called after `To` and
`from` have been promoted to the same number types, so methods do not have to
worry about number type conversions. Instead, all they have to do is provide the
correct conversion formula.
"""
_convert(::Type{To}, ::From) where {To, From} = error("No known conversion from $From to $To")
_convert(::Type{To}, c::To) where {To} = c # Identity conversion
_convert(::Type{To}, c::EastWestCoordinate, ::NorthSouthCoordinate) where {To <: EastWestCoordinate} = _convert(To, c) # Drop north south coordinate if not needed
_convert(::Type{To}, c::NorthSouthCoordinate, ::EastWestCoordinate) where {To <: NorthSouthCoordinate} = _convert(To, c) # Drop east west coordinate if not needed
_convert(::Type{To}, c::From, origin::NorthSouthCoordinate) where {To <: NorthSouthCoordinate, From <: NorthSouthCoordinate} =
    _convert(To, From(origin) + c) - _convert(To, origin)

# All matching coordinate types are implicitly convertible
Base.convert(::Type{C}, c::Coordinate{N}) where {N, C <: Coordinate{N}} = C(c)

EastWestCoordinate(c::EastWestCoordinate) = c
EastWestCoordinate(c::EastWestCoordinate, ::NorthSouthCoordinate) = c
EastWestCoordinate(::NorthSouthCoordinate, c::EastWestCoordinate) = c
NorthSouthCoordinate(c::NorthSouthCoordinate) = c
NorthSouthCoordinate(c::NorthSouthCoordinate, ::EastWestCoordinate) = c
NorthSouthCoordinate(::EastWestCoordinate, c::NorthSouthCoordinate) = c

Base.eltype(::Type{<:Coordinate{N, T}}) where {N, T <: Number} = T

Base.length(::Coordinate{N}) where N = N

Base.iterate(c::Coordinate{1}, state...) = iterate((c.v,), state...)
Base.iterate(c::Coordinate{2}, state...) = iterate((c.c1[], c.c2[]), state...)
Base.iterate(c::ECEF, state...) = iterate((c.v1, c.v2, c.v3), state...)

Base.getindex(c::Coordinate{1}) = c.v
function Base.getindex(c::Coordinate{1}, i::Int)
    if i == 1; return c.v; end
    throw(BoundsError(c, i))
end
Base.getindex(c::Coordinate{2}) = (c.c1[], c.c2[])
function Base.getindex(c::Coordinate{2}, i::Int)
    if i == 1; return c.c1[]; end
    if i == 2; return c.c2[]; end
    throw(BoundsError(c, i))
end
Base.getindex(c::ECEF) = (c.v1, c.v2, c.v3)
function Base.getindex(c::ECEF, i::Int)
    if i == 1; return c.v1; end
    if i == 2; return c.v2; end
    if i == 3; return c.v3; end
    throw(BoundsError(c, i))
end

Base.show(io::IO, c::Coordinate{2}) = print(io, typeof(c), "(", c[1], ", ", c[2], ")")

# Arithmetic
Base.promote_rule(::Type{<:Coordinate{N,T1}}, ::Type{<:Coordinate{N,T2}}) where {N, T1 <: Number, T2 <: Number} = Coordinate{N,promote_type(T1, T2)}
for op in (:+, :-)
    @eval begin
        Base.$op(c1::Coordinate{N}, c2::Coordinate{N}) where {N} = $op(promote(c1, c2)...)
        Base.$op(c1::Coordinate{N,T}, c2::Coordinate{N,T}) where {N,T <: Number} = throw(MethodError($op, (c1, c2)))
        Base.$op(c1::C, c2::C) where {C <: Coordinate} = C(map($op, c1[], c2[]))
        Base.$op(c::C) where {C <: Coordinate} = C(map($op, c[]))
    end
end

function needs_latitude(A,B)
    for (To,From) in ((A,B), (B,A))
        @eval begin
            _convert(::Type{$To}, ::$From) = error("Cannot convert from $($From) to $($To) without knowing the latitude. Call $($From)(::$($To), ::NorthSouthCoordinate) instead.")
            _convert(::Type{$To}, from::$From, ns::NorthSouthCoordinate) = _convert($To, from, Lat(ns))
        end
    end
end


#########################
# WMX, WMY <--> Lat, Lon

_convert(::Type{WMX}, (lon,)::Lon) = lon / 180
_convert(::Type{Lon}, (wmx,)::WMX) = 180 * wmx

function _convert(::Type{WMY}, (lat,)::Lat)
    @assert abs(lat) <= 90
    return log(tand((lat+90)/2)) / π
end
_convert(::Type{Lat}, (wmy,)::WMY) = 2*atand(exp(π * wmy)) - 90


############################
# Lat, Lon <--> East, North

# WGS84 ellipsoid
const a = 6.378137e6
const inv_f = 298.257_223_563
const f = inv(inv_f)
const b = a * (1-f)
const e2 = (2 - f) / inv_f

function N(lat)
    s,c = sincosd(lat)
    return a^2/sqrt((a*c)^2 + (b*s)^2)
end
ror_from_lat(lat) = N(lat) * cosd(lat) # Radius of rotation, i.e. radius orthogonal to z-axis.


needs_latitude(Lon, East)
_convert(::Type{Lon}, (east,)::East, (lat,)::Lat) = 180*east / (π*ror_from_lat(lat))
_convert(::Type{East}, (lon,)::Lon, (lat,)::Lat) = π*ror_from_lat(lat)*lon / 180

function _convert(::Type{North}, (lat,)::Lat)
    @assert abs(lat) <= 90
    return a * (1-f)^2 * Elliptic.Pi(e2, lat*π/180, e2)
end

function _convert(::Type{Lat}, (north,)::North)
    @assert abs(north) <= North(Lat(90))[]
    return bisect(lat -> North(Lat(lat))[] - north, -90.0, 90.0)
end

function bisect(f, a, b)
    fa = f(a)
    fb = f(b)
    @assert sign(fa) != sign(fb)
    while true
        m = (a+b)/2
        fm = f(m)
        if a == m || b == m
            break
        end
        if sign(fa) == sign(fm)
            a,fa = m,fm
        else
            b,fb = m,fm
        end
    end
    return a
end

for op in (:+, :-)
    @eval begin
        function Base.$op(c::Coordinate{2,T}, d::EastNorth{T}) where {T <: Number}
            typeof(c)(
                EastWestCoordinate(c) + typeof(EastWestCoordinate(c))(East(d), Lat(c)),
                typeof(NorthSouthCoordinate(c))(North(c) + North(d)),
            )
        end
    end
end


############################
# East, North <--> WMX, WMY

needs_latitude(East, WMX)
_convert(::Type{East}, wmx::WMX, lat::Lat) = _convert(East, Lon(wmx, lat), lat)
_convert(::Type{WMX}, east::East, lat::Lat) = _convert(WMX, Lon(east, lat), lat)

_convert(::Type{North}, wmy::WMY) = _convert(North, Lat(wmy))
_convert(::Type{WMY}, north::North) = _convert(WMY, Lat(north))

end # module MapMaths
