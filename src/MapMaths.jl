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
using Base.ScopedValues
using StaticArrays

include("utils.jl")
include("coordinate_utils.jl")
# include("ellipses.jl")

##################
# Geoid & Georef

@conversion_constructible abstract type Geoid end
@conversion_constructible abstract type Ellipsoid <: Geoid end
@conversion_constructible abstract type Spheroid <: Ellipsoid end

@conversion_constructible abstract type Georef end

@iterable @ampersand @combo struct Datum
    (Geoid, Georef)
end

struct Wgs84Spheroid <: Spheroid end
struct Wgs84Georef <: Georef end
Wgs84() = Wgs84Spheroid() & Wgs84Georef()

##############################################
# Abstract coordinate types and combinations

@conversion_constructible abstract type Coordinate end
@conversion_constructible abstract type GeoidlessCoordinate <: Coordinate end

@convertible abstract type SpaceCoordinate <: Coordinate end
@static_length(SpaceCoordinate, 3)

@conversion_constructible abstract type GeoidlessSpaceCoordinate <: GeoidlessCoordinate end

@convertible abstract type SurfaceCoordinate <: GeoidlessCoordinate end
@static_length(SurfaceCoordinate, 2)

@convertible abstract type LongitudeCoordinate <: GeoidlessCoordinate end
@static_length(LongitudeCoordinate, 1)

@convertible abstract type LatitudeCoordinate <: GeoidlessCoordinate end
@static_length(LatitudeCoordinate, 1)

@convertible abstract type AltitudeCoordinate <: GeoidlessCoordinate end
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

##############################################
# Coordinate and geoid / georef combinations

@iterable @ampersand @combo struct GeoidCoordinate <: Coordinate
    (GeoidlessCoordinate, Geoid)
end
@symmetric Base.:&((c1, geoid)::GeoidCoordinate, c2::GeoidlessCoordinate) = (c1 & c2) & geoid
@symmetric Base.:&(c::Coordinate, ::Geoid) = c

@iterable @ampersand @combo struct GeoidSpaceCoordinate <: SpaceCoordinate
    (GeoidlessSpaceCoordinate, Geoid)
end
(C::Type{<:GeoidlessCoordinate})(c::GeoidSpaceCoordinate) = C(GeoidlessSpaceCoordinate(c))

@iterable @ampersand @combo struct Georeffed
    (Coordinate, Georef)
end
@symmetric Base.:&(coord::Coordinate, (geoid, georef)::Datum) = coord & georef
@symmetric Base.:&(coord::GeoidlessCoordinate, (geoid, georef)::Datum) = (coord & geoid) & georef
@symmetric Base.:&((c1, georef)::Georeffed, c2::Coordinate) = (c1 & c2) & georef
@symmetric Base.:&((coord, georef)::Georeffed{<:GeoidlessCoordinate}, geoid::Geoid) = (coord & geoid) & georef

struct CoordinateType{C <: Coordinate} <: Coordinate end
struct GeoidlessCoordinateType{C <: GeoidlessCoordinate} <: GeoidlessCoordinate end
CoordinateType(C::Type{<:Coordinate}) = CoordinateType{C}()
CoordinateType(C::Type{<:GeoidlessCoordinate}) = GeoidlessCoordinateType{C}()

@symmetric Base.:&(C::Type{<:Coordinate}, geoid::Geoid) = CoordinateType(C) & geoid
@symmetric Base.:&(C::Type{<:Coordinate}, georef::Georef) = CoordinateType(C) & georef
@symmetric Base.:&(C::Type{<:Coordinate}, (geoid, georef)::Datum) = (C & geoid) & georef

(::CoordinateType{C})(c::Coordinate) where {C} = C(c)
(::GeoidlessCoordinateType{C})(c::GeoidlessCoordinate) where {C} = C(c)

Base.show(io::IO, ::CoordinateType{C}) where {C} = print(io, C)
Base.show(io::IO, ::GeoidlessCoordinateType{C}) where {C} = print(io, C)

###########
# Routing

function convert_georef(georef1::Georef, (c, georef2)::Georeffed{SpaceCoordinate})
    if georef1 != georef2
        error("No conversion from $georef2 to $georef1")
    end
    return c & georef1
end

((C, georef)::Georeffed)(c::Georeffed{SpaceCoordinate}) = C(convert_georef(georef, c)) & georef
function ((C, georef1)::Georeffed)((c, georef2)::Georeffed)
    if georef1 != georef2
        error("Georef conversion can only be done on space coordinates")
    end
    return C(c) & georef1
end
(C::Coordinate)((c, georef2)::Georeffed) = C(c)

function ((C, geoid1)::GeoidCoordinate)((c, geoid2)::Union{GeoidCoordinate, GeoidSpaceCoordinate})
    if geoid1 == geoid2
        return C(c) & geoid1
    else
        return (C & geoid1)(Ecef(c))
    end
end

(C::Type{<:LongitudeCoordinate})(c::GeoidlessSpaceCoordinate) = C(SurfaceCoordinate(c))
(C::Type{<:LatitudeCoordinate})(c::GeoidlessSpaceCoordinate) = C(SurfaceCoordinate(c))

(C::Type{<:LongitudeLatitudeCoordinate})(c::SurfaceCoordinate) = C(LongitudeCoordinate(c), LatitudeCoordinate(c))
(C::Type{<:LatitudeLongitudeCoordinate})(c::SurfaceCoordinate) = C(LatitudeCoordinate(c), LongitudeCoordinate(c))

#############################
# Concrete coordinate types

@coordinate_number_tuple(Ecef <: SpaceCoordinate, 3)
@coordinate_number_tuple(Alt <: AltitudeCoordinate, 1)

@coordinate_number_tuple(Lon <: LongitudeCoordinate, 1)
@coordinate_number_tuple(Lat <: LatitudeCoordinate, 1)

@named_coordinate_combo LonLat = LongitudeLatitudeCoordinate{Lon, Lat}
@named_coordinate_combo LonLatAlt = SurfaceAltitudeCoordinate{LonLat, Alt}
@named_coordinate_combo LatLon = LatitudeLongitudeCoordinate{Lat, Lon}
@named_coordinate_combo LatLonAlt = SurfaceAltitudeCoordinate{LatLon, Alt}

@coordinate_number_tuple(WmX <: LongitudeCoordinate, 1)
@coordinate_number_tuple(WmY <: LatitudeCoordinate, 1)
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