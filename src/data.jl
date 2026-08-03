"""
    Wgs84()

The World Geodetic System 1984 datum (EPSG code 4326).
"""
struct Wgs84 <: Datum end
Ellipse(::Wgs84) = Ellipse(SemiMajorAxis(6378137.0), InverseFlattening(298.257223563))
