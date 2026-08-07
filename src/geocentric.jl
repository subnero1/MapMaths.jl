"""
    Longitude()

Angle in degrees between the prime meridian and the meridian plane of the point.
East is positive, west is negative.
"""
struct Longitude <: ScalarCoordinateSystem end

"""
    CylindricalRadius()

Distance in meters from the Z-axis.
"""
struct CylindricalRadius <: ScalarCoordinateSystem end

"""
    Radius()

Distance in meters from the origin.
"""
struct Radius <: ScalarCoordinateSystem end

"""
    GeocentricLatitude()

Angle in degrees between the equatorial plane and the line from the origin.
North is positive, south is negative.
"""
struct GeocentricLatitude <: ScalarCoordinateSystem end

"""
    GeocentricAltitude()

Distance in meters from the surface of the reference ellipsoid along the line to the origin.
"""
struct GeocentricAltitude <: ScalarCoordinateSystem end

"""
    Cylindrical() = Longitude() & CylindricalRadius() & CartesianZ()
"""
@simple_vector_coordinate_system(Cylindrical, Longitude, CylindricalRadius, CartesianZ)

"""
    Spherical() = Longitude() & GeocentricLatitude() & Radius()
"""
@simple_vector_coordinate_system(Spherical, Longitude, GeocentricLatitude, Radius)

"""
    GeocentricLonLat() = Longitude() & GeocentricLatitude()
"""
@simple_vector_coordinate_system(GeocentricLonLat, Longitude, GeocentricLatitude)

"""
    GeocentricLatLon() = GeocentricLatitude() & Longitude()
"""
@simple_vector_coordinate_system(GeocentricLatLon, GeocentricLatitude, Longitude)

"""
    GeocentricLatAlt() = GeocentricLatitude() & GeocentricAltitude()
"""
@simple_vector_coordinate_system(GeocentricLatAlt, GeocentricLatitude, GeocentricAltitude)

"""
    GeocentricLonLatAlt() = Longitude() & GeocentricLatitude() & GeocentricAltitude()
"""
@simple_vector_coordinate_system(GeocentricLonLatAlt, Longitude, GeocentricLatitude, GeocentricAltitude)

"""
    GeocentricLatLonAlt() = GeocentricLatitude() & Longitude() & GeocentricAltitude()
"""
@simple_vector_coordinate_system(GeocentricLatLonAlt, GeocentricLatitude, Longitude, GeocentricAltitude)

const Lon = Longitude
const GeocentricLat = GeocentricLatitude
const GeocentricAlt = GeocentricAltitude

###

lossy_parent(::Longitude) = CartesianXY()
lossy_parent(::CylindricalRadius) = CartesianXY()
lossy_parent(::Radius) = CylindricalRadius() & CartesianZ()
lossy_parent(::GeocentricLatitude) = CylindricalRadius() & CartesianZ()
lossy_parent(::GeocentricAltitude) = GeocentricLatitude() & Radius()

@symmetric lossless_parent(::Longitude, ::CylindricalRadius) = CartesianXY()
@symmetric lossless_parent(::GeocentricLatitude, ::Radius) = CylindricalRadius() & CartesianZ()
@symmetric lossless_parent(::GeocentricLatitude, ::GeocentricAltitude) = GeocentricLatitude() & Radius()

@symmetric lossless_parent(::Longitude, ::GeocentricLatitude, ::GeocentricAltitude) = Longitude() & GeocentricLatitude() & Radius()
@symmetric lossless_parent(::Longitude, ::GeocentricLatitude, ::Radius) = Longitude() & CylindricalRadius() & CartesianZ()
@symmetric lossless_parent(::Longitude, ::CylindricalRadius, ::CartesianZ) = CartesianX() & CartesianY() & CartesianZ()

###

cmap_impl(::Longitude, (x, y)::Coordinate{CartesianXY}, datum) = atand(y, x)
cmap_impl(::CylindricalRadius, (x, y)::Coordinate{CartesianXY}, datum) = sqrt(x^2 + y^2)
cmap_impl(::CylindricalRadius, (lat, r)::Coordinate(GeocentricLatitude() & Radius()), datum) = r * cosd(lat)
cmap_impl(::CartesianZ, (lat, r)::Coordinate(GeocentricLatitude() & Radius()), datum) = r * sind(lat)
cmap_impl(::Radius, (x, y, z)::Coordinate{Cartesian}, datum) = sqrt(x^2 + y^2 + z^2)
cmap_impl(::Radius, (r, z)::Coordinate(CylindricalRadius() & CartesianZ()), datum) = sqrt(r^2 + z^2)

function cmap_impl(::Radius, (lat, alt)::Coordinate(GeocentricLatitude() & GeocentricAltitude()), datum)
    (s, c) = sincosd(lat)
    (a, b) = datum[SemiMajorAxis], datum[SemiMinorAxis]
    return alt + sqrt((a * b)^2 / ((b * c)^2 + (a * s)^2))
end

cmap_impl(::GeocentricLatitude, (r, z)::Coordinate(CylindricalRadius() & CartesianZ()), datum) = atand(z, r)

function cmap_impl(::GeocentricAltitude, (lat, r)::Coordinate(GeocentricLatitude() & Radius()), datum)
    (s, c) = sincosd(lat)
    (a, b) = datum[SemiMajorAxis], datum[SemiMinorAxis]
    return r - sqrt((a * b)^2 / ((b * c)^2 + (a * s)^2))
end

cmap_impl(::CartesianX, (lon, r)::Coordinate(Longitude() & CylindricalRadius()), datum) = r * cosd(lon)
cmap_impl(::CartesianY, (lon, r)::Coordinate(Longitude() & CylindricalRadius()), datum) = r * sind(lon)
cmap_impl(::CartesianXY, (lon, r)::Coordinate(Longitude() & CylindricalRadius()), datum) = r .* reverse(sincosd(lon))

function cmap_impl(::Cartesian, (lon, lat, r)::Coordinate(Longitude() & GeocentricLatitude() & Radius()), datum)
    sin_lon, cos_lon = sincosd(lon)
    sin_lat, cos_lat = sincosd(lat)
    return (
        r * cos_lat * cos_lon,
        r * cos_lat * sin_lon,
        r * sin_lat,
    )
end
