@testset "Test `shift`" begin
    lattice = [
        -3.0179389205999998 -3.0179389205999998 0.0000000000000000
        -5.2272235447000002 5.2272235447000002 0.0000000000000000
        0.0000000000000000 0.0000000000000000 -9.7736219469000005
    ]
    positions = [[2 / 3, 1 / 3, 1 / 4], [1 / 3, 2 / 3, 3 / 4]]
    atoms = [1, 1]
    cell = Cell(lattice, positions, atoms)
    shifted = shift(cell, 1, 1, 1)
    @test shift(shifted, 2, 2, 2) == shift(cell, 3, 3, 3)
end
