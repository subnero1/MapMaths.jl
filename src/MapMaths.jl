module MapMaths

export numtype,
    promote_numtype,
    default_numtype,
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

flatten() = ()
flatten(x, y...) = (x, flatten(y...)...)
flatten(x::Tuple, y...) = (flatten(x...)..., flatten(y...)...)

permutations(a) = ((a,))
permutations(a, b) = ((a, b), (b, a))
permutations(a, b, c) = ((a, b, c), (a, c, b), (b, a, c), (b, c, a), (c, a, b), (c, b, a))

numtype(x) = numtype(typeof(x))
numtype(T::Type) = throw(MethodError(numtype, (T,)))
numtype(T::Type{<:Number}) = T
numtype(T::Type{<:NTuple{<:Any, Number}}) = promote_numtype(fieldtypes(T)...)

promote_numtype(args...) = promote_type(numtype.(args)...)
default_numtype(args...) = float(promote_numtype(args...))

abstract type Combineable end
struct TempCombo{T <: Tuple} <: Combineable
    parts::T
end
Base.:|(parts::Vararg{Combineable}) = TempCombo(parts)
Base.:|(combo::TempCombo, part::Combineable) = |(combo.parts..., part)
Base.show(io::IO, combo::TempCombo) = join(io, combo.parts, " | ")

macro symmetric(fun_def)
    @assert fun_def.head in (:function, :(=))
    @assert fun_def.args[1].head == :call
    name = fun_def.args[1].args[1]
    args = fun_def.args[1].args[2:end]
    return esc(
        quote
            $fun_def
            $(
                map(
                    permuted_args -> begin
                        if all(permuted_args .== args)
                            return nothing
                        end
                        return :($name($(permuted_args...)) = $name($(args...)))
                    end,
                    permutations(args...),
                )...
            )
        end,
    )
    return nothing
end

abstract type Coordinate <: Combineable end
Base.length(c::Coordinate) = length(typeof(c))
Base.iterate(c::Coordinate, state...) = iterate(Tuple(c), state...)

macro coordinate_dimension(Dim::Symbol, N::Integer)
    return esc(quote
        abstract type $Dim <: Coordinate end
        (::Type{C})(c::C) where {C <: $Dim} = c
        Base.convert(::Type{C}, c::$Dim) where {C <: $Dim} = C(c)
        Base.length(::Type{<:$Dim}) = $N
    end)
end

macro number_constructors(Type::Symbol, N)
    return esc(quote
        $Type(numbers::Vararg{Number, $N}) = $Type{default_numtype(numbers)}(numbers)
        (::Type{V})(numbers::NTuple{$N, Number}) where {V <: $Type} = V(numbers...)
        (::Type{V})(numbers::StaticVector{$N, <:Number}) where {V <: $Type} = V(numbers...)
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
        (::Type{C})(c::$Type) where {C <: $Type} = C(c.coords)
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
            (::Type{C})(c::$Type) where {C <: $Type} = C(c.coords)
            $(
                map(
                    ((i, SubDim),) -> quote
                        (::Type{C})(c::$Type) where {C <: $SubDim} = C(c.coords[$i])
                    end,
                    enumerate(SubDims),
                )...
            )

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

            function Base.show(io::IO, c::$Name{T}) where {T <: Number}
                print(io, $Name, "{", T, "}(")
                join(io, Tuple(c), ", ")
                return print(io, ")")
            end
        end,
    )
end

abstract type Datum end
struct WithDatum{T, D <: Datum}
    value::T
    datum::D
end
const MaybeWithDatum{T} = Union{T, WithDatum{T}}

Base.iterate((; value, datum)::WithDatum, state...) = iterate((value, datum), state...)
Base.:(|)(value, datum::Datum) = WithDatum(value, datum)
Base.show(io::IO, (value, datum)::WithDatum) = print(io, "WithDatum(", value, ", ", datum, ")")

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

(::Type{C})(c::SpaceCoordinate) where {C <: EastCoordinate} = C(SurfaceCoordinate(c))
(::Type{C})(c::SpaceCoordinate) where {C <: NorthCoordinate} = C(SurfaceCoordinate(c))

(::Type{C})(c::EastNorthCoordinate) where {C <: NorthEastCoordinate} = C((c.coords[2], c.coords[1]))
(::Type{C})(c::NorthEastCoordinate) where {C <: EastNorthCoordinate} = C((c.coords[2], c.coords[1]))
(::Type{C})(c::SurfaceCoordinate) where {C <: NorthEastCoordinate} = C(EastNorthCoordinate(c))

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

(::Type{C})((wmx,)::WmX) where {C <: Lon} = C(wmx * 180)
(::Type{C})((lon,)::Lon) where {C <: WmX} = C(lon / 180)
(::Type{C})(c::EastCoordinate) where {C <: WmX} = C(Lon(c))

(::Type{C})((wmy,)::WmY) where {C <: Lat} = C(2 * atand(exp(π * wmy)) - 90)
function (::Type{C})((lat,)::Lat) where {C <: WmY}
    @assert abs(lat) <= 90
    return C(log(tand((lat + 90) / 2)) / π)
end
(::Type{C})(c::NorthCoordinate) where {C <: WmY} = C(Lat(c))

end