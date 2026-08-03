# struct Longitude <: ScalarCoordinateSystem end   <---   Reused from geocentric.jl
"""
    GeodeticLatitude()

Angle in degrees between the equatorial plane and the normal vector to the
reference ellipsoid. North is positive, south is negative.
"""
struct GeodeticLatitude <: ScalarCoordinateSystem end

"""
    GeodeticAltitude()

Distance in meters from the surface of the reference ellipsoid along the normal
vector to the ellipsoid.
"""
struct GeodeticAltitude <: ScalarCoordinateSystem end

"""
    GeodeticLonLat() = Longitude() & GeodeticLatitude()
"""
@simple_vector_coordinate_system(GeodeticLonLat, Longitude, GeodeticLatitude)

"""
    GeodeticLatLon() = GeodeticLatitude() & Longitude()
"""
@simple_vector_coordinate_system(GeodeticLatLon, GeodeticLatitude, Longitude)

"""
    GeodeticLatAlt() = GeodeticLatitude() & GeodeticAltitude()
"""
@simple_vector_coordinate_system(GeodeticLatAlt, GeodeticLatitude, GeodeticAltitude)

"""
    GeodeticLonLatAlt() = Longitude() & GeodeticLatitude() & GeodeticAltitude()
"""
@simple_vector_coordinate_system(GeodeticLonLatAlt, Longitude, GeodeticLatitude, GeodeticAltitude)

"""
    GeodeticLatLonAlt() = GeodeticLatitude() & Longitude() & GeodeticAltitude()
"""
@simple_vector_coordinate_system(GeodeticLatLonAlt, GeodeticLatitude, Longitude, GeodeticAltitude)

const GeodeticLat = GeodeticLatitude
const GeodeticAlt = GeodeticAltitude
const Lat = GeodeticLatitude
const Alt = GeodeticAltitude
const LonLat = GeodeticLonLat
const LatLon = GeodeticLatLon
const LatAlt = GeodeticLatAlt
const LonLatAlt = GeodeticLonLatAlt
const LatLonAlt = GeodeticLatLonAlt

###

lossy_parent(::Lat) = CylindricalRadius() & CartesianZ()
lossy_parent(::Alt) = CylindricalRadius() & CartesianZ()

@symmetric lossy_parent(::Lon, ::Lat) = Cylindrical()
@symmetric lossy_parent(::Lon, ::Alt) = Cylindrical()
@symmetric lossy_parent(::Lat, ::Alt) = CylindricalRadius() & CartesianZ()

@symmetric lossless_parent(::Lon, ::Lat, ::Alt) = Cylindrical()

########################
# Cartesian from LonLat

function cmap_impl(
        ::CylindricalRadius,
        (lat, alt)::Coordinate{LatAlt},
        datum,
    )
    sin_lat, cos_lat = sincosd(lat)
    prime_vertical = datum[SemiMajorAxis] / sqrt(1 - datum[SquaredEccentricity] * sin_lat^2)
    return (prime_vertical + alt) * cos_lat
end

function cmap_impl(
        ::CartesianZ,
        (lat, alt)::Coordinate{LatAlt},
        datum,
    )
    sin_lat = sind(lat)
    prime_vertical = datum[SemiMajorAxis] / sqrt(1 - datum[SquaredEccentricity] * sin_lat^2)
    return ((1 - datum[SquaredEccentricity]) * prime_vertical + alt) * sin_lat
end

function cmap_impl(
        ::typeof(CylindricalRadius() & CartesianZ()),
        (lat, alt)::Coordinate{LatAlt},
        datum,
    )
    sin_lat, cos_lat = sincosd(lat)
    prime_vertical = datum[SemiMajorAxis] / sqrt(1 - datum[SquaredEccentricity] * sin_lat^2)
    return (
        (prime_vertical + alt) * cos_lat,
        ((1 - datum[SquaredEccentricity]) * prime_vertical + alt) * sin_lat,
    )
end

########################
# LatLon from Cartesian

function cmap_impl(::typeof(GeodeticLatitude()), (p, z)::Coordinate(CylindricalRadius() & CartesianZ()), datum)
    return atand(
        reverse(
            unnormalized_normal(
                Ellipse(datum),
                (p, z)
            )
        )...
    )
end

function cmap_impl(::typeof(GeodeticAltitude()), (p, z)::Coordinate(CylindricalRadius() & CartesianZ()), datum)
    cos_lat, sin_lat = normal(Ellipse(datum), (p, z))
    prime_vertical = datum[SemiMajorAxis] / sqrt(1 - datum[SquaredEccentricity] * sin_lat^2)
    return if abs(cos_lat) > abs(sin_lat)
        p / cos_lat - prime_vertical
    else
        z / sin_lat - prime_vertical * (1 - datum[SquaredEccentricity])
    end
    # return (p + abs(z) - prime_vertical * (cos_lat + (1 - datum[SquaredEccentricity]) * abs(sin_lat))) / (cos_lat + abs(sin_lat))
end

function cmap_impl(::typeof(GeodeticLatitude() & GeodeticAltitude()), (p, z)::Coordinate(CylindricalRadius() & CartesianZ()), datum)
    cos_lat, sin_lat = normal(Ellipse(datum), (p, z))
    prime_vertical = datum[SemiMajorAxis] / sqrt(1 - datum[SquaredEccentricity] * sin_lat^2)
    return (
        atand(sin_lat, cos_lat),
        if abs(cos_lat) > abs(sin_lat)
            p / cos_lat - prime_vertical
        else
            z / sin_lat - prime_vertical * (1 - datum[SquaredEccentricity])
        end,
    )
end
