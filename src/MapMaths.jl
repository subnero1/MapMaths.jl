module MapMaths

using Elliptic

export
    WebMercator, WMX, WMY,
    LatLon, LonLat, Lat, Lon,
    EastNorth, East, North

abstract type Coordinate{N, T <: Number} end

abstract type EastWestCoordinate{T <: Number} <: Coordinate{1, T} end
abstract type NorthSouthCoordinate{T <: Number} <: Coordinate{1, T} end

for (C,D) in (
    (:WMX, :EastWestCoordinate), (:WMY, :NorthSouthCoordinate),
    (:Lon, :EastWestCoordinate), (:Lat, :NorthSouthCoordinate),
    (:East, :EastWestCoordinate), (:North, :NorthSouthCoordinate),
)
    O = if D == :EastWestCoordinate; :NorthSouthCoordinate else :EastWestCoordinate end
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
        $C(c::$D{T}, o::$O) where {T <: Number} = $C{T}(c, o) # Preserve eltype when converting
        $C{T}(c::$D, o::$O) where {T <: Number} = $C{T}($D{T}(c), $O{T}(o)) # Match eltypes
        $C{T}(c::$D{T}, o::$O{T}) where {T <: Number} = $C{T}(_convert($C, c, o)) # Dispatch to _convert
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
    end
end

include("docstrings.jl")

"""
    _convert(To::Type{<:Coordinate}, from::Coordinate, other...)-> to::Number

Convert coordinate `from` to type `To`. This function is called after `To` and
`from` have been promoted to the same number types, so methods do not have to
worry about number type conversions. Instead, all they have to do is provide the
correct conversion formula.
"""
_convert(::Type{To}, ::From) where {To, From} = error("No known conversion from $From to $To")
_convert(::Type{To}, c::Coordinate{1}, ::Coordinate{1}) where {To <: Coordinate{1}} = _convert(To, c) # Drop second coordinate if not needed

# All matching coordinate types are implicitly convertible
Base.convert(::Type{C}, c::EastWestCoordinate) where {C <: EastWestCoordinate} = C(c)
Base.convert(::Type{C}, c::NorthSouthCoordinate) where {C <: NorthSouthCoordinate} = C(c)
Base.convert(::Type{C}, c::Coordinate{2}) where {C <: Coordinate{2}} = C(c)

Base.eltype(::Type{<:Coordinate{N, T}}) where {N, T <: Number} = T

Base.length(::Coordinate{N}) where N = N

Base.iterate(c::Coordinate{1}, state...) = iterate((c.v,), state...)
Base.iterate(c::Coordinate{2}, state...) = iterate((c.c1[], c.c2[]), state...)

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

Base.show(io::IO, c::Coordinate{2}) = print(io, typeof(c), "(", c[1], ", ", c[2], ")")

function needs_latitude(A,B)
    for (To,From) in ((A,B), (B,A))
        @eval begin
            _convert(::Type{$To}, ::$From) = error("Cannot convert from $($From) to $($To) without knowing the latitude. Call `$($From)(::$($To), ::NorthSouthCoordinate)` instead.")
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


############################
# East, North <--> WMX, WMY

needs_latitude(East, WMX)
_convert(::Type{East}, wmx::WMX, lat::Lat) = _convert(East, _convert(Lon, wmx), lat)
_convert(::Type{WMX}, east::East, lat::Lat) = _convert(WMX, _convert(Lon, east), lat)

_convert(::Type{North}, wmy::WMY) = _convert(North, _convert(Lat, wmy))
_convert(::Type{WMY}, north::North) = _convert(WMY, _convert(Lat, north))

end # module MapMaths
