@testset "Test `ismagnetic`" begin
    # Scalar (collinear) moments
    @test !ismagnetic(MagneticAtom(:Fe, 0.0))
    @test !ismagnetic(MagneticAtom(:Fe, 1e-9))
    @test ismagnetic(MagneticAtom(:Fe, 1e-7))
    # Vector (non-collinear) moments
    @test !ismagnetic(MagneticAtom(:Fe, [0.0, 0.0, 0.0]))
    @test !ismagnetic(MagneticAtom(:Fe, [1e-9, 0.0, 0.0]))
    @test ismagnetic(MagneticAtom(:Fe, [1e-7, 0.0, 0.0]))
    # Custom tolerance
    @test !ismagnetic(MagneticAtom(:Fe, 1e-5), 1e-3)
    @test ismagnetic(MagneticAtom(:Fe, 1e-2), 1e-3)
    # MagneticCell behaviour (sum of moments)
    lattice = Matrix{Float64}(I, 3, 3)
    positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    cell_nonmag = Cell(lattice, positions, [MagneticAtom(:Fe, 0.0), MagneticAtom(:Fe, 0.0)])
    cell_mag = Cell(lattice, positions, [MagneticAtom(:Fe, 1e-6), MagneticAtom(:Fe, 1e-6)])
    cell_cancel = Cell(
        lattice, positions, [MagneticAtom(:Fe, 1e-6), MagneticAtom(:Fe, -1e-6)]
    )
    @test !ismagnetic(cell_nonmag)
    @test ismagnetic(cell_mag)
    @test !ismagnetic(cell_cancel)
end
