using MapMaths
using Test

@testset "Number constructors" begin
    Lon(1)
    LonLat(1, 2)
    LonLatAlt(1, 2, 3)

    @test numtype(Lon(1)) == Float64
    @test numtype(Lon{Int}(1)) == Int
    @test numtype(Lon(1.0f0)) == Float32

    @test numtype(LonLat(1, 2)) == Float64
    @test numtype(LonLat(1, 2.0f0)) == Float32
    @test numtype(LonLat(1.0f0, 2.0f0)) == Float32
    @test numtype(LonLat(1.0f0, 2.0)) == Float64
    @test numtype(LonLat{Int}(1, 2)) == Int
    @test numtype(LonLat{Int}(1.0, 2.0)) == Int
end
@testset "Conversion constructors" begin
    @test Lon(LonLat(1, 2)) == Lon(1)
    @test Lon(LatLon(1, 2)) == Lon(2)
    @test Lon(LonLatAlt(1, 2, 3)) == Lon(1)
    @test Lon(LatLonAlt(1, 2, 3)) == Lon(2)
    @test Alt(LonLatAlt(1, 2, 3)) == Alt(3)
    @test LonLat(LonLatAlt(1, 2, 3)) == LonLat(1, 2)
    @test LonLat(LatLonAlt(1, 2, 3)) == LonLat(2, 1)
end

@testset "Combinations" begin
    @test Lon(1) | Lat(2) == LonLat(1, 2)
    @test Lat(2) | Lon(1) == LonLat(1, 2)
    @test Lon(1) | Lat(2) | Alt(3) == LonLatAlt(1, 2, 3)
    @test Lon(1) | Alt(3) | Lat(2) == LonLatAlt(1, 2, 3)
    @test Lat(2) | Lon(1) | Alt(3) == LonLatAlt(1, 2, 3)
    @test Lat(2) | Alt(3) | Lon(1) == LonLatAlt(1, 2, 3)
    @test Alt(3) | Lon(1) | Lat(2) == LonLatAlt(1, 2, 3)
    @test Alt(3) | Lat(2) | Lon(1) == LonLatAlt(1, 2, 3)
end

@testset "Coordinate show" begin
    repr(Lon(1)) == "Lon{Float64}(1.0)"
    repr(LonLat(1, 2)) == "LonLat{Float64}(1.0, 2.0)"
    repr(LonLatAlt(1, 2, 3)) == "LonLatAlt{Float64}(1.0, 2.0, 3.0)"
end