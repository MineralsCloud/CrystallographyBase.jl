using Mendeleev: elements
using Unitful: uconvert

# Data from https://materialsproject.org/materials/mp-2657/
@testset "Test `crystaldensity`" begin
    lattice = Lattice(
        4.59983732u"angstrom", 4.59983732u"angstrom", 2.95921356u"angstrom", 90, 90, 90
    )
    positions = [
        [0, 0, 0.5],
        [0.5, 0.5, 0.0],
        [0.19567869, 0.80432131, 0],
        [0.80432131, 0.19567869, 0],
        [0.30432131, 0.30432131, 0.5],
        [0.69567869, 0.69567869, 0.5],
    ]
    @testset "Test `Element`s" begin
        atoms = elements[[22, 22, 8, 8, 8, 8]]
        cell = Cell(lattice, positions, atoms)
        @test cellvolume(cell) ≈ 62.612530083185085u"angstrom^3"
        @test uconvert(u"g/cm^3", crystaldensity(lattice, atoms)) ≈
            4.236179319948116u"g/cm^3"
        @test uconvert(u"g/cm^3", crystaldensity(cell)) ≈ 4.236179319948116u"g/cm^3"
    end
    @testset "Test `AbstractString`s" begin
        atoms = elements[["Titanium", "Titanium", "Oxygen", "Oxygen", "Oxygen", "Oxygen"]]
        cell = Cell(lattice, positions, atoms)
        @test uconvert(u"g/cm^3", crystaldensity(cell)) ≈ 4.236179319948116u"g/cm^3"
    end
    @testset "Test `Integer`s" begin
        atoms = [22, 22, 8, 8, 8, 8]
        cell = Cell(lattice, positions, atoms)
        @test uconvert(u"g/cm^3", crystaldensity(cell)) ≈ 4.236179319948116u"g/cm^3"
    end
    @testset "Test `Symbol`s" begin
        atoms = [:Ti, :Ti, :O, :O, :O, :O]
        cell = Cell(lattice, positions, atoms)
        @test uconvert(u"g/cm^3", crystaldensity(cell)) ≈ 4.236179319948116u"g/cm^3"
    end
end
