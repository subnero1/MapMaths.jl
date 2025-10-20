using MapMaths
using Test

@testset "Number constructors" begin
    Lon(1)
    LonLat(1, 2)
    LonLatAlt(1, 2, 3)

    @test MapMaths.numtype(Lon(1)) == Float64
    @test MapMaths.numtype(Lon{Int}(1)) == Int
    @test MapMaths.numtype(Lon(1.0f0)) == Float32

    @test MapMaths.numtype(LonLat(1, 2)) == Float64
    @test MapMaths.numtype(LonLat(1, 2.0f0)) == Float32
    @test MapMaths.numtype(LonLat(1.0f0, 2.0f0)) == Float32
    @test MapMaths.numtype(LonLat(1.0f0, 2.0)) == Float64
    @test MapMaths.numtype(LonLat{Int}(1, 2)) == Int
    @test MapMaths.numtype(LonLat{Int}(1.0, 2.0)) == Int
end

@testset "Combinations" begin
    @test Lon(1) & Lat(2) == Lat(2) & Lon(1) == LonLat(1, 2)
    @test Lon(1) & Lat(2) & Alt(3) ==
          Lon(1) & Alt(3) & Lat(2) ==
          Lat(2) & Lon(1) & Alt(3) ==
          Lat(2) & Alt(3) & Lon(1) ==
          Alt(3) & Lon(1) & Lat(2) ==
          Alt(3) & Lat(2) & Lon(1) ==
          LonLatAlt(1, 2, 3)

    @test Wgs84() & Lon(1) == Lon(1) & Wgs84()
    @test Lon(1) & Wgs84() & Lat(2) == Wgs84() & Lon(1) & Lat(2) == LonLat(1, 2) & Wgs84()
    @test Lon(1) & Lat(2) & Wgs84() & Alt(3) ==
          Lon(1) & Wgs84() & Lat(2) & Alt(3) ==
          Wgs84() & Lon(1) & Lat(2) & Alt(3) ==
          LonLatAlt(1, 2, 3) & Wgs84()
    @test Ecef(1, 2, 3) & Wgs84() ==
          Ecef(1, 2, 3) & Wgs84Georef() ==
          Wgs84() & Ecef(1, 2, 3) ==
          Wgs84Georef() & Ecef(1, 2, 3)
end

@testset "Conversion constructors" begin
    @test Lon(LonLat(1, 2)) == Lon(1)
    @test Lon(LatLon(1, 2)) == Lon(2)
    @test Lon(LonLatAlt(1, 2, 3)) == Lon(1)
    @test Lon(LatLonAlt(1, 2, 3)) == Lon(2)
    @test Alt(LonLatAlt(1, 2, 3)) == Alt(3)
    @test LonLat(LonLatAlt(1, 2, 3)) == LonLat(1, 2)
    @test LonLat(LatLonAlt(1, 2, 3)) == LonLat(2, 1)

    @test Lon(LonLat(1, 2) & Wgs84()) == Lon(1)
    @test Lon(LatLon(1, 2) & Wgs84()) == Lon(2)
    @test Lon(LonLatAlt(1, 2, 3) & Wgs84()) == Lon(1)
    @test Lon(LatLonAlt(1, 2, 3) & Wgs84()) == Lon(2)
    @test Alt(LonLatAlt(1, 2, 3) & Wgs84()) == Alt(3)
    @test LonLat(LonLatAlt(1, 2, 3) & Wgs84()) == LonLat(1, 2)
    @test LonLat(LatLonAlt(1, 2, 3) & Wgs84()) == LonLat(2, 1)

    @test (Lon & Wgs84())(LonLat(1, 2) & Wgs84()) == Lon(1) & Wgs84()
    @test (Lon & Wgs84())(LatLon(1, 2) & Wgs84()) == Lon(2) & Wgs84()
    @test (Lon & Wgs84())(LonLatAlt(1, 2, 3) & Wgs84()) == Lon(1) & Wgs84()
    @test (Lon & Wgs84())(LatLonAlt(1, 2, 3) & Wgs84()) == Lon(2) & Wgs84()
    @test (Alt & Wgs84())(LonLatAlt(1, 2, 3) & Wgs84()) == Alt(3) & Wgs84()
    @test (LonLat & Wgs84())(LonLatAlt(1, 2, 3) & Wgs84()) == LonLat(1, 2) & Wgs84()
    @test (LonLat & Wgs84())(LatLonAlt(1, 2, 3) & Wgs84()) == LonLat(2, 1) & Wgs84()
end

# @testset "getindex" begin
#     @test Lon(1)[Lon] == 1
#     @test LonLat(1, 2)[Lon] == 1
#     @test LonLat(1, 2)[Lat] == 2
#     @test LonLatAlt(1, 2, 3)[Lon] == 1
#     @test LonLatAlt(1, 2, 3)[Lat] == 2
#     @test LonLatAlt(1, 2, 3)[Alt] == 3

#     @test (Lon(1)&Wgs84())[Lon&Wgs84()] == 1
#     @test (LonLat(1, 2)&Wgs84())[Lon&Wgs84()] == 1
#     @test (LonLat(1, 2)&Wgs84())[Lat&Wgs84()] == 2
#     @test (LonLatAlt(1, 2, 3)&Wgs84())[Lon&Wgs84()] == 1
#     @test (LonLatAlt(1, 2, 3)&Wgs84())[Lat&Wgs84()] == 2
#     @test (LonLatAlt(1, 2, 3)&Wgs84())[Alt&Wgs84()] == 3
# end

@testset "show" begin
    @test repr(Lon(1)) == "Lon{Float64}(1.0)"
    @test repr(LonLat(1, 2)) == "LonLat{Float64}(1.0, 2.0)"
    @test repr(LonLatAlt(1, 2, 3)) == "LonLatAlt{Float64}(1.0, 2.0, 3.0)"
    @test repr(Lon(1) & Wgs84()) == "Lon{Float64}(1.0) & Wgs84Spheroid() & Wgs84Georef()"
    @test repr(Lon & Wgs84()) == "Lon & Wgs84Spheroid() & Wgs84Georef()"
end
