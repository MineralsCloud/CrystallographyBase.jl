using Brillouin

@testset "Test `WignerSeitzCell` against `Brillouin`" begin
    # Data from https://thchr.github.io/Brillouin.jl/dev/wignerseitz/
    lattice = Lattice([0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0])
    ws = WignerSeitzCell(lattice)
    cᴿ = wignerseitz([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])
    @test ws.verts == [
        [0.25, -0.75, 0.25],
        [-0.5, -0.5, 0.5],
        [-0.75, 0.25, 0.25],
        [0.75, -0.25, -0.25],
        [-0.25, -0.25, -0.25],
        [0.25, 0.25, -0.75],
        [-0.5, 0.5, -0.5],
        [0.5, -0.5, -0.5],
        [-0.25, -0.25, 0.75],
        [0.5, -0.5, 0.5],
        [-0.25, 0.75, -0.25],
        [0.25, 0.25, 0.25],
        [0.5, 0.5, -0.5],
        [-0.5, 0.5, 0.5],
    ]
    @test ws.faces == [
        [3, 7, 5, 2],
        [6, 8, 5, 7],
        [7, 3, 14, 11],
        [11, 13, 6, 7],
        [14, 3, 2, 9],
        [1, 2, 5, 8],
        [1, 10, 9, 2],
        [13, 11, 14, 12],
        [14, 9, 10, 12],
        [4, 8, 6, 13],
        [4, 10, 1, 8],
        [13, 12, 10, 4],
    ]
    @test cartesianize(cᴿ).verts == WignerSeitzCell(lattice, true).verts
    @testset "Construct Brillouin zones from `WignerSeitzCell`" begin
        # Data from https://thchr.github.io/Brillouin.jl/dev/wignerseitz/
        reci_lattice = inv(lattice)
        ws = WignerSeitzCell(reci_lattice)
        cᴳ = wignerseitz([[-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0]])
        @test ws.verts == [
            [0.25, -0.5, -0.25],
            [0.5, -0.25, 0.25],
            [0.25, -0.25, 0.5],
            [-0.25, 0.25, -0.5],
            [0.25, -0.25, -0.5],
            [-0.25, -0.5, 0.25],
            [0.25, 0.75, 0.5],
            [0.75, 0.25, 0.5],
            [0.25, 0.5, 0.75],
            [0.5, 0.25, 0.75],
            [-0.25, -0.75, -0.5],
            [-0.5, -0.75, -0.25],
            [-0.5, -0.25, -0.75],
            [-0.25, -0.5, -0.75],
            [-0.75, -0.25, -0.5],
            [-0.75, -0.5, -0.25],
            [-0.25, 0.25, 0.5],
            [-0.5, -0.25, 0.25],
            [-0.5, 0.25, -0.25],
            [-0.25, 0.5, 0.25],
            [0.5, 0.25, -0.25],
            [0.25, 0.5, -0.25],
            [0.5, 0.75, 0.25],
            [0.75, 0.5, 0.25],
        ]
        @test ws.faces == [
            [9, 7, 20, 17],
            [24, 21, 22, 23],
            [8, 24, 23, 7, 9, 10],
            [10, 3, 2, 8],
            [5, 21, 24, 8, 2, 1],
            [1, 11, 14, 5],
            [2, 3, 6, 12, 11, 1],
            [5, 14, 13, 4, 22, 21],
            [19, 20, 7, 23, 22, 4],
            [4, 13, 15, 19],
            [14, 11, 12, 16, 15, 13],
            [9, 17, 18, 6, 3, 10],
            [17, 20, 19, 15, 16, 18],
            [18, 16, 12, 6],
        ]
        @test cartesianize(cᴳ).verts == WignerSeitzCell(reci_lattice, true).verts
    end
end
