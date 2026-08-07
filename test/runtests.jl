using MapMaths
using Test
using Documenter
using Combinatorics
using StaticArrays
using Logging

function test_oneway_cmap(input, output, maybe_datum...; atol = 0, rtol = 0)
    for n in 1:MapMaths.n_parts(output)
        for parts in combinations(MapMaths.parts(output), n)
            expected = (&)(parts...)
            @test isapprox(SVector(cmap(coordinate_system(expected), input, maybe_datum...)...), SVector(expected...); atol, rtol)
        end
    end
    return
end

function test_cmap(c1, c2, maybe_datum...)
    test_oneway_cmap(c1, c2, maybe_datum...)
    test_oneway_cmap(c2, c1, maybe_datum...)
    return
end

DocMeta.setdocmeta!(
    MapMaths,
    :DocTestSetup,
    quote
        using MapMaths
        using MapMaths: numtype, parts, flatten, @symmetric
    end;
    recursive = true,
)
with_logger(ConsoleLogger(stderr, Logging.Error)) do
    doctest(MapMaths; manual = false, testset = "doctests")
end

@testset "geocentric" begin
    test_cmap(Coordinate(Cylindrical(), 0, 1, 2), Coordinate(Cartesian(), 1, 0, 2))
    test_cmap(Coordinate(Cylindrical(), 90, 1, 2), Coordinate(Cartesian(), 0, 1, 2))

    test_cmap(Coordinate(Spherical(), 0, 0, 1), Coordinate(Cartesian(), 1, 0, 0))
    test_cmap(Coordinate(Spherical(), 90, 0, 1), Coordinate(Cartesian(), 0, 1, 0))
    test_cmap(Coordinate(Spherical(), 0, 90, 1), Coordinate(Cartesian(), 0, 0, 1))

    test_cmap(Coordinate(Spherical(), 0, 0, 1), Coordinate(Cylindrical(), 0, 1, 0))
    test_cmap(Coordinate(Spherical(), 0, 90, 1), Coordinate(Cylindrical(), 0, 0, 1))

    test_cmap(Coordinate(GeocentricLonLatAlt(), 0, 0, 0), Coordinate(Cartesian(), Wgs84()[SemiMajorAxis], 0, 0), Wgs84())
    test_cmap(Coordinate(GeocentricLonLatAlt(), 0, 0, 1), Coordinate(Cartesian(), Wgs84()[SemiMajorAxis] + 1, 0, 0), Wgs84())
end

@testset "geodetic" begin
    test_cmap(Coordinate(LonLatAlt(), 0, 0, 0), Coordinate(Cartesian(), Wgs84()[SemiMajorAxis], 0, 0), Wgs84())
    test_cmap(Coordinate(LonLatAlt(), 0, 0, 1), Coordinate(Cartesian(), Wgs84()[SemiMajorAxis] + 1, 0, 0), Wgs84())
    test_cmap(Coordinate(LonLatAlt(), 0, 90, 0), Coordinate(Cartesian(), 0, 0, Wgs84()[SemiMinorAxis]), Wgs84())
end

@testset "webmercator" begin
    test_cmap(Coordinate(WebMercator(), 0, 0), Coordinate(LonLat(), 0, 0))
    test_cmap(Coordinate(WebMercator(), 0.5, 0), Coordinate(LonLat(), 90, 0))
    test_cmap(Coordinate(WebMercator(), 0, Inf), Coordinate(LonLat(), 0, 90))
    test_cmap(Coordinate(WebMercatorAlt(), 0, 0, 0), Coordinate(LonLatAlt(), 0, 0, 0))
    test_cmap(Coordinate(WebMercatorAlt(), 0.5, 0, 0), Coordinate(LonLatAlt(), 90, 0, 0))
    test_cmap(Coordinate(WebMercatorAlt(), 0, Inf, 0), Coordinate(LonLatAlt(), 0, 90, 0))
end

@testset "local" begin
    test_cmap(Coordinate(AzimuthElevationRange(), 0, 0, 1), Coordinate(EastNorthUp(), 0, 1, 0))
    test_cmap(Coordinate(AzimuthElevationRange(), 90, 0, 1), Coordinate(EastNorthUp(), 1, 0, 0))
    test_cmap(Coordinate(AzimuthElevationRange(), 0, 90, 1), Coordinate(EastNorthUp(), 0, 0, 1))

    test_cmap(Coordinate(Origin(Wgs84(), LonLatAlt(), 0, 0, 0) + EastNorthUp(), 1, 0, 0), Coordinate(Cartesian(), Wgs84()[SemiMajorAxis], 1, 0), Wgs84())
    test_cmap(Coordinate(Origin(Wgs84(), LonLatAlt(), 0, 0, 0) + EastNorthUp(), 0, 1, 0), Coordinate(Cartesian(), Wgs84()[SemiMajorAxis], 0, 1), Wgs84())
    test_cmap(Coordinate(Origin(Wgs84(), LonLatAlt(), 0, 0, 0) + EastNorthUp(), 0, 0, 1), Coordinate(Cartesian(), Wgs84()[SemiMajorAxis] + 1, 0, 0), Wgs84())
    test_cmap(Coordinate(Origin(Wgs84(), LonLatAlt(), 90, 0, 0) + EastNorthUp(), 1, 0, 0), Coordinate(Cartesian(), -1, Wgs84()[SemiMajorAxis], 0), Wgs84())
    test_cmap(Coordinate(Origin(Wgs84(), LonLatAlt(), 90, 0, 0) + EastNorthUp(), 0, 1, 0), Coordinate(Cartesian(), 0, Wgs84()[SemiMajorAxis], 1), Wgs84())
    test_cmap(Coordinate(Origin(Wgs84(), LonLatAlt(), 90, 0, 0) + EastNorthUp(), 0, 0, 1), Coordinate(Cartesian(), 0, Wgs84()[SemiMajorAxis] + 1, 0), Wgs84())

    test_cmap(
        Coordinate(Origin(Wgs84(), LonLatAlt(), 0, 0, 0) + EastNorthUp(), Wgs84()[SemiMajorAxis], 0, 0),
        Coordinate(Origin(Wgs84(), LonLatAlt(), 90, 0, 0) + AzimuthElevationRange(), -90, 0, Wgs84()[SemiMajorAxis]),
    )
end
