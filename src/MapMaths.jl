module MapMaths

export GeoReference, Geoid, Ellipsoid, Spheroid, Georef, Georeffed, Datum, Wgs84Spheroid, Wgs84Georef, Wgs84

export SpaceCoordinate,
    SurfaceCoordinate,
    LongitudeCoordinate,
    LatitudeCoordinate,
    AltitudeCoordinate,
    Ecef,
    Alt,
    Lon,
    Lat,
    LonLat,
    LonLatAlt,
    LatLon,
    LatLonAlt,
    WmX,
    WmY,
    WebMercator,
    WebMercatorAlt

using .Meta: isexpr
using StaticArrays

include("utils.jl")
# include("ellipses.jl")

@conversion_constructible abstract type Geoid end
@conversion_constructible abstract type Ellipsoid <: Geoid end
@conversion_constructible abstract type Spheroid <: Ellipsoid end

@conversion_constructible abstract type Georef end

@ampersand @combo struct Datum
    (Geoid, Georef)
end
Base.iterate(datum::Datum, state...) = iterate((Geoid(datum), Georef(datum)), state...)

struct Wgs84Spheroid <: Spheroid end
struct Wgs84Georef <: Georef end
Wgs84() = Wgs84Spheroid() & Wgs84Georef()

macro static_length(Type, N)
    Type = esc(Type)
    N = esc(N)
    return :(Base.length(::Union{$Type, Type{<:$Type}}) = $N)
end

"""
    @coordinate_combo struct Type [<: Supertype]
        (PartType1, PartType2, ...)
    end

Extension of [`@combo`](@ref) which defines additional methods specific to
coordinate types.
"""
macro coordinate_combo(struct_def::Expr)
    @assert isexpr(struct_def, :struct)
    @assert struct_def.args[1] == false # mutability
    Type_maybe_sub_Supertype = struct_def.args[2]
    body = struct_def.args[3]
    @assert isexpr(body, :block)
    @assert isexpr(body.args[2], :tuple)
    Parts = esc.(body.args[2].args)
    PartTypeVars = ntuple(i -> Symbol("P$i"), length(Parts))
    PartTypeVarDefs = map((V, T) -> :($V <: $T), PartTypeVars, Parts)
    if Type_maybe_sub_Supertype isa Symbol
        Type = esc(Type_maybe_sub_Supertype)
        Supertype = nothing
        type_signature = :($Type{$(PartTypeVarDefs...)})
    elseif isexpr(Type_maybe_sub_Supertype, :(<:), 2)
        (Type, Supertype) = esc.(Type_maybe_sub_Supertype.args)
        type_signature = :($Type{$(PartTypeVarDefs...)} <: $Supertype)
    else
        error("Unexpected type signature: $Type_maybe_sub_Supertype")
    end
    return quote
        $(
            if !isnothing(Supertype)
                quote
                    if hasmethod(length, Tuple{Type{$Supertype}})
                        @assert sum(length.(tuple($(Parts...)))) == length($Supertype)
                    end
                end
            end
        )
        $(combo(struct_def))
        @static_length($Type, sum(length.(tuple($(Parts...)))))
        Base.Tuple(c::$Type) = flatten(Tuple.(parts(c)))
        function MapMaths.numtype(::Type{$Type{$(PartTypeVars...)}}) where {$(PartTypeVarDefs...)}
            return promote_numtype($(PartTypeVars...))
        end
    end
end

"""
    @named_coordinate_combo Name = ComboType{PartType1, PartType2, ...}

Add a type alias for a combination of coordinate types, along with several
method definitions to make the alias behave like a [`@number_tuple`](@ref).
"""
macro named_coordinate_combo(expr::Expr)
    @assert isexpr(expr, :(=), 2)
    @assert expr.args[1] isa Symbol
    Name = esc(expr.args[1])
    ComboType_and_PartTypes = expr.args[2]
    @assert isexpr(ComboType_and_PartTypes, :curly)
    ComboType = esc(ComboType_and_PartTypes.args[1])
    PartTypes = esc.(ComboType_and_PartTypes.args[2:end])
    part_vars = ntuple(i -> Symbol("part_$i"), length(PartTypes))
    part_var_defs = map((part, Part) -> :($part::supertype($Part)), part_vars, PartTypes)
    return quote
        const $Name{T <: Number} = $ComboType{$(map(PartType -> :($PartType{T}), PartTypes)...)}

        function $Name{T}(parts::Vararg{Number, length($ComboType)}) where {T <: Number}
            tail_1 = parts
            $(
                map(
                    ((i, PartType),) -> quote
                        $(Symbol("part_$i")) = $PartType{T}($(Symbol("tail_$i"))[1:length($PartType)])
                        $(Symbol("tail_$(i+1)")) = $(Symbol("tail_$i"))[length($PartType)+1:end]
                    end,
                    enumerate(PartTypes),
                )...
            )
            return $ComboType($(ntuple(i -> Symbol("part_$i"), length(PartTypes))...))
        end
        $Name(numbers::Vararg{Number, length($ComboType)}) = $Name{default_numtype(numbers)}(numbers)
        (C::Type{<:$Name})(numbers::NTuple{length($ComboType), Number}) = C(numbers...)
        (C::Type{<:$Name})(numbers::StaticVector{length($ComboType), <:Number}) = C(numbers...)
        (C::Type{<:$Name})(c::$Name) = C(parts(c)...)
        $Name($(part_var_defs...)) = $Name{promote_numtype($(part_vars...))}($(part_vars...))

        function Base.show(io::IO, c::$Name)
            print(io, $Name, "{", numtype(c), "}(")
            join(io, Tuple(c), ", ")
            return print(io, ")")
        end
    end
end

@conversion_constructible abstract type Coordinate end
@conversion_constructible abstract type GeoidlessCoordinate <: Coordinate end

@convertible abstract type SpaceCoordinate <: Coordinate end
@static_length(SpaceCoordinate, 3)

@conversion_constructible abstract type GeoidlessSpaceCoordinate <: GeoidlessCoordinate end

@convertible abstract type SurfaceCoordinate <: GeoidlessCoordinate end
@convertible abstract type LongitudeCoordinate <: GeoidlessCoordinate end
@convertible abstract type LatitudeCoordinate <: GeoidlessCoordinate end
@convertible abstract type AltitudeCoordinate <: GeoidlessCoordinate end

@static_length(SurfaceCoordinate, 2)
@static_length(LongitudeCoordinate, 1)
@static_length(LatitudeCoordinate, 1)
@static_length(AltitudeCoordinate, 1)

@ampersand @coordinate_combo struct LongitudeLatitudeCoordinate <: SurfaceCoordinate
    (LongitudeCoordinate, LatitudeCoordinate)
end
@coordinate_combo struct LatitudeLongitudeCoordinate <: SurfaceCoordinate
    (LatitudeCoordinate, LongitudeCoordinate)
end

@ampersand @coordinate_combo struct LongitudeAltitudeCoordinate <: GeoidlessCoordinate
    (LongitudeCoordinate, AltitudeCoordinate)
end
@ampersand @coordinate_combo struct LatitudeAltitudeCoordinate <: GeoidlessCoordinate
    (LatitudeCoordinate, AltitudeCoordinate)
end

@ampersand @coordinate_combo struct SurfaceAltitudeCoordinate <: GeoidlessSpaceCoordinate
    (SurfaceCoordinate, AltitudeCoordinate)
end
@symmetric Base.:&(c::LongitudeAltitudeCoordinate, lat::LatitudeCoordinate) =
    (LongitudeCoordinate(c) & lat) & AltitudeCoordinate(c)
@symmetric Base.:&(c::LatitudeAltitudeCoordinate, lon::LongitudeCoordinate) =
    (lon & LatitudeCoordinate(c)) & AltitudeCoordinate(c)

@ampersand @combo struct GeoidCoordinate <: Coordinate
    (GeoidlessCoordinate, Geoid)
end
Base.iterate(c::GeoidCoordinate, state...) = iterate((GeoidlessCoordinate(c), Geoid(c)), state...)
@symmetric Base.:&((c1, geoid)::GeoidCoordinate, c2::GeoidlessCoordinate) = (c1 & c2) & geoid

@ampersand @combo struct GeoidSpaceCoordinate <: SpaceCoordinate
    (GeoidlessSpaceCoordinate, Geoid)
end
(C::Type{<:GeoidlessCoordinate})(c::GeoidSpaceCoordinate) = C(GeoidlessSpaceCoordinate(c))
Base.iterate(c::GeoidSpaceCoordinate, state...) = iterate((GeoidlessSpaceCoordinate(c), Geoid(c)), state...)

@ampersand @combo struct Georeffed
    (Coordinate, Georef)
end
Base.iterate(c::Georeffed, state...) = iterate((Coordinate(c), Georef(c)), state...)
@symmetric Base.:&(coord::Coordinate, (geoid, georef)::Datum) = coord & georef
@symmetric Base.:&(coord::GeoidlessCoordinate, (geoid, georef)::Datum) = (coord & geoid) & georef
@symmetric Base.:&((c1, georef)::Georeffed, c2::Coordinate) = (c1 & c2) & georef
@symmetric Base.:&((coord, georef)::Georeffed{<:GeoidlessCoordinate}, geoid::Geoid) = (coord & geoid) & georef

(C::Type{<:LongitudeCoordinate})(c::GeoidlessSpaceCoordinate) = C(SurfaceCoordinate(c))
(C::Type{<:LatitudeCoordinate})(c::GeoidlessSpaceCoordinate) = C(SurfaceCoordinate(c))

(C::Type{<:LongitudeLatitudeCoordinate})(c::SurfaceCoordinate) = C(LongitudeCoordinate(c), LatitudeCoordinate(c))
(C::Type{<:LatitudeLongitudeCoordinate})(c::SurfaceCoordinate) = C(LatitudeCoordinate(c), LongitudeCoordinate(c))

@number_tuple(Ecef <: SpaceCoordinate, 3)
@number_tuple(Alt <: AltitudeCoordinate, 1)

@number_tuple(Lon <: LongitudeCoordinate, 1)
@number_tuple(Lat <: LatitudeCoordinate, 1)

@named_coordinate_combo LonLat = LongitudeLatitudeCoordinate{Lon, Lat}
@named_coordinate_combo LonLatAlt = SurfaceAltitudeCoordinate{LonLat, Alt}
@named_coordinate_combo LatLon = LatitudeLongitudeCoordinate{Lat, Lon}
@named_coordinate_combo LatLonAlt = SurfaceAltitudeCoordinate{LatLon, Alt}

@number_tuple(WmX <: LongitudeCoordinate, 1)
@number_tuple(WmY <: LatitudeCoordinate, 1)
@named_coordinate_combo WebMercator = LongitudeLatitudeCoordinate{WmX, WmY}
@named_coordinate_combo WebMercatorAlt = SurfaceAltitudeCoordinate{WebMercator, Alt}

(C::Type{<:Lon})((wmx,)::WmX) = C(wmx * 180)
(C::Type{<:WmX})((lon,)::Lon) = C(lon / 180)
(C::Type{<:WmX})(c::LongitudeCoordinate) = C(Lon(c))

(C::Type{<:Lat})((wmy,)::WmY) = C(2 * atand(exp(π * wmy)) - 90)
function (C::Type{<:WmY})((lat,)::Lat)
    @assert abs(lat) <= 90
    return C(log(tand((lat + 90) / 2)) / π)
end
(C::Type{<:WmY})(c::LatitudeCoordinate) = C(Lat(c))

end