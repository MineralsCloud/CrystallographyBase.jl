using LinearAlgebra: dot, norm, diagm
using SymPy: symbols

@testset "Test consistency of constructors of `MetricTensor`" begin
    lattice = Lattice(
        4.59983732u"angstrom", 4.59983732u"angstrom", 2.95921356u"angstrom", 90, 90, 90
    )
    𝐚, 𝐛, 𝐜 = basisvectors(lattice)
    @test MetricTensor(lattice) ==
        MetricTensor(𝐚, 𝐛, 𝐜) ==
        MetricTensor(
            4.59983732u"angstrom", 4.59983732u"angstrom", 2.95921356u"angstrom", 90, 90, 90
        )
    @test Lattice(MetricTensor(lattice)) == lattice
    @testset "Test `zero`, `one`, and `oneunit`" begin
        @test one(MetricTensor(lattice)) == MetricTensor(one(lattice))
        @test oneunit(MetricTensor(lattice)) == MetricTensor(oneunit(lattice))
        @test zero(MetricTensor(lattice)) == MetricTensor(zero(lattice))
    end
end

@testset "Test lengths in a hexagonal lattice" begin
    g = MetricTensor(1, 1, 2, 90, 90, 120)  # Primitive hexagonal
    @test g ≈ MetricTensor([1, 0, 0], [-1//2, sqrt(3) / 2, 0], [0, 0, 2])
    @test g ≈ MetricTensor(Lattice([
        1.0 -1/2 0.0
        0.0 sqrt(3)/2 0.0
        0.0 0.0 2.0
    ]))
    𝐚 = ReducedCoordinates(1, 2, 1)
    @test lengthof(𝐚, g)^2 ≈ 7
end

@testset "Test distance between atoms in a hexagonal lattice" begin
    g = MetricTensor(1, 1, 2, 90, 90, 120)  # Primitive hexagonal
    a = ReducedCoordinates(1, 1, 1)
    b = ReducedCoordinates(1//3, 1//3, 1//2)
    @test distance(a, g, b)^2 == 13 / 9
end

@testset "Test symbolic calculation" begin
    a, c = symbols("a, c"; positive=true)
    @test MetricTensor(a, a, c, 90, 90, 120) ==
        MetricTensor([a^2 -0.5*a^2 0; -0.5*a^2 a^2 0; 0 0 c^2])  # Primitive hexagonal
end
