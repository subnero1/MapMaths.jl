module MapMaths

export numtype,
    promote_numtype,
    default_numtype,
    SemiMajorAxis,
    SemiMinorAxis,
    Flattening,
    InverseFlattening,
    SquaredEccentricity,
    Eccentricity,
    Datum,
    Coordinate,
    Wgs84,
    SpaceCoordinate,
    SurfaceCoordinate,
    EastCoordinate,
    NorthCoordinate,
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

using StaticArrays

include("utils.jl")
include("ellipses.jl")

abstract type Geoid <: Combineable end
abstract type Georeference <: Combineable end

abstract type Coordinate <: Combineable end
Base.length(c::Coordinate) = length(typeof(c))
Base.iterate(c::Coordinate, state...) = iterate(Tuple(c), state...)

Base.getindex(c::Coordinate, C::Type{<:Coordinate}) = maybe_only(C(c)...)
Base.getindex(cd::Combo(Coordinate, Datum), (C, d)::Combo(Type{<:Coordinate}, Datum)) = maybe_only(C((C | d)(cd))...)

Base.:|(d::Datum, c::Coordinate) = c | d
Base.:|((c1, d)::Combo(Coordinate, Datum), c2::Coordinate) = (c1 | c2) | d

(C::Type{<:Coordinate})((c, d)::Combo(Coordinate, Datum)) = C(c)
((C, d_out)::Combo(Type{<:Coordinate}, Datum))((c, d_in)::Combo(Coordinate, Datum)) =
    if d_out == d_in
        return C(c) | d_out
    else
        error("TODO")
        # return C((Spatial | d_out)(c))
    end

macro coordinate_dimension(Dim::Symbol, N::Integer)
    return esc(quote
        abstract type $Dim <: Coordinate end
        (::Type{C})(c::C) where {C <: $Dim} = c
        # Base.convert(C::Type{<:$Dim}, c::$Dim) = C(c)
        # Base.convert(C::Type{<:Combo($Dim, Datum)}, c::Combo($Dim, Datum)) = C(c)
        Base.length(::Type{<:$Dim}) = $N
    end)
end

macro number_constructors(Type::Symbol, N)
    return esc(quote
        $Type(numbers::Vararg{Number, $N}) = $Type{default_numtype(numbers)}(numbers)
        (V::Type{<:$Type})(numbers::NTuple{$N, Number}) = V(numbers...)
        (V::Type{<:$Type})(numbers::StaticVector{$N, <:Number}) = V(numbers...)
    end)
end

macro base_coordinate(Type_sub_Dim::Expr)
    @assert Type_sub_Dim.head == :(<:)
    (Type::Symbol, Dim::Symbol) = Type_sub_Dim.args
    N = :(length($Dim))
    return esc(quote
        struct $Type{T <: Number} <: $Dim
            coords::NTuple{$N, T}
            $Type{T}(coords::Vararg{Number, $N}) where {T <: Number} = new{T}(coords)
        end
        @number_constructors($Type, $N)
        (C::Type{<:$Type})(c::$Type) = C(c.coords)
        Base.Tuple(c::$Type) = c.coords
        numtype(::Type{$Type{T}}) where {T <: Number} = T
        function Base.show(io::IO, c::$Type)
            print(io, typeof(c), "(")
            join(io, Tuple(c), ", ")
            return print(io, ")")
        end
    end)
end

macro combo_coordinate(Type_sub_Dim::Expr, SubDims::Symbol...)
    @assert Type_sub_Dim.head == :(<:)
    (Type::Symbol, Dim::Symbol) = Type_sub_Dim.args
    FieldTypeVars = ntuple(i -> Symbol("C$i"), length(SubDims))
    FieldTypeVarDefs = map((V, T) -> :($V <: $T), FieldTypeVars, SubDims)

    return esc(
        quote
            struct $Type{$(FieldTypeVarDefs...)} <: $Dim
                coords::Tuple{$(FieldTypeVars...)}
                function $Type{$(FieldTypeVars...)}(coords::Tuple{$(SubDims...)}) where {$(FieldTypeVarDefs...)}
                    return new{$(FieldTypeVars...)}(coords)
                end
            end

            $Type(coords::Tuple{$(SubDims...)}) = $Type{typeof.(coords)...}(coords)
            (C::Type{<:$Type})(c::$Type) = C(c.coords)
            $(map(((i, SubDim),) -> quote
                (C::Type{<:$SubDim})(c::$Type) = C(c.coords[$i])
            end, enumerate(SubDims))...)

            function numtype(::Type{$Type{$(FieldTypeVars...)}}) where {$(FieldTypeVarDefs...)}
                return promote_numtype($(FieldTypeVars...))
            end
            Base.Tuple(c::$Type) = flatten(Tuple.(c.coords))

            function Base.show(io::IO, c::$Type)
                print(io, $Type, "(")
                join(io, c.coords, ", ")
                return print(io, ")")
            end
        end,
    )
end

macro named_combo_coordinate(Name::Symbol, ComboType::Symbol, FieldTypes::Symbol...)
    return esc(
        quote
            const $Name{T <: Number} = $ComboType{$(map(FieldType -> :($FieldType{T}), FieldTypes)...)}

            function $Name(coords::Tuple{$(map(FieldType -> :(supertype($FieldType)), FieldTypes)...)})
                return $Name{promote_numtype(coords...)}(coords)
            end
            function $Name{T}(coords::Vararg{Number, length($Name)}) where {T <: Number}
                tail_1 = coords
                $(
                    map(
                        ((i, FieldType),) -> quote
                            $(Symbol("coord_$i")) = $FieldType{T}($(Symbol("tail_$i"))[1:length($FieldType)])
                            $(Symbol("tail_$(i+1)")) = $(Symbol("tail_$i"))[length($FieldType)+1:end]
                        end,
                        enumerate(FieldTypes),
                    )...
                )
                return $ComboType(tuple($(ntuple(i -> Symbol("coord_$i"), length(FieldTypes))...)))
            end
            @number_constructors($Name, length($Name))

            function Base.show(io::IO, c::$Name)
                print(io, $Name, "{", numtype(c), "}(")
                join(io, Tuple(c), ", ")
                return print(io, ")")
            end
        end,
    )
end

struct Wgs84 <: Datum end
AxisAlignedEllipse(::Wgs84) = SemiMajorAxis(6378137.0) | InverseFlattening(298.257223563)

@coordinate_dimension(SpaceCoordinate, 3)
@coordinate_dimension(SurfaceCoordinate, 2)
@coordinate_dimension(EastCoordinate, 1)
@coordinate_dimension(NorthCoordinate, 1)
@coordinate_dimension(AltitudeCoordinate, 1)

@combo_coordinate(EastNorthCoordinate <: SurfaceCoordinate, EastCoordinate, NorthCoordinate)
@combo_coordinate(NorthEastCoordinate <: SurfaceCoordinate, NorthCoordinate, EastCoordinate)
@combo_coordinate(SurfaceAltitudeCoordinate <: SpaceCoordinate, SurfaceCoordinate, AltitudeCoordinate)

@symmetric Base.:|(east::EastCoordinate, north::NorthCoordinate) = EastNorthCoordinate((east, north))
@symmetric Base.:|(surf::SurfaceCoordinate, alt::AltitudeCoordinate) = SurfaceAltitudeCoordinate((surf, alt))
@symmetric Base.:|(east::EastCoordinate, north::NorthCoordinate, alt::AltitudeCoordinate) =
    SurfaceAltitudeCoordinate((EastNorthCoordinate((east, north)), alt))

(C::Type{<:EastCoordinate})(c::SpaceCoordinate) = C(SurfaceCoordinate(c))
(C::Type{<:NorthCoordinate})(c::SpaceCoordinate) = C(SurfaceCoordinate(c))

(C::Type{<:NorthEastCoordinate})(c::EastNorthCoordinate) = C((c.coords[2], c.coords[1]))
(C::Type{<:EastNorthCoordinate})(c::NorthEastCoordinate) = C((c.coords[2], c.coords[1]))
(C::Type{<:NorthEastCoordinate})(c::SurfaceCoordinate) = C(EastNorthCoordinate(c))

@base_coordinate(Ecef <: SpaceCoordinate)
@base_coordinate(Alt <: AltitudeCoordinate)

@base_coordinate(Lon <: EastCoordinate)
@base_coordinate(Lat <: NorthCoordinate)
@named_combo_coordinate(LonLat, EastNorthCoordinate, Lon, Lat)
@named_combo_coordinate(LonLatAlt, SurfaceAltitudeCoordinate, LonLat, Alt)
@named_combo_coordinate(LatLon, NorthEastCoordinate, Lat, Lon)
@named_combo_coordinate(LatLonAlt, SurfaceAltitudeCoordinate, LatLon, Alt)

@base_coordinate(WmX <: EastCoordinate)
@base_coordinate(WmY <: NorthCoordinate)
@named_combo_coordinate(WebMercator, EastNorthCoordinate, WmX, WmY)
@named_combo_coordinate(WebMercatorAlt, SurfaceAltitudeCoordinate, WebMercator, Alt)

(C::Type{<:Lon})((wmx,)::WmX) = C(wmx * 180)
(C::Type{<:WmX})((lon,)::Lon) = C(lon / 180)
(C::Type{<:WmX})(c::EastCoordinate) = C(Lon(c))

(C::Type{<:Lat})((wmy,)::WmY) = C(2 * atand(exp(π * wmy)) - 90)
function (C::Type{<:WmY})((lat,)::Lat)
    @assert abs(lat) <= 90
    return C(log(tand((lat + 90) / 2)) / π)
end
(C::Type{<:WmY})(c::NorthCoordinate) = C(Lat(c))

end