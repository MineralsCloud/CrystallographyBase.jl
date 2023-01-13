@testset "Test the simplest `supercell` for a `Lattice`" begin
    lattice = Lattice(2, 1, 5, 90, 90, 90; axis=:c)
    a = supercell(lattice, 3)
    b = Lattice(6, 3, 15, 90, 90, 90; axis=:c)
    @test a == b
    @test reciprocal(a) == reciprocal(b)
end

@testset "Test `supercell` from primitive to body-centered for a `Lattice`" begin
    lattice = Lattice(2, 1, 5, 90, 90, 90; axis=:c)
    P = [
        0 1 1
        1 0 1
        1 1 0
    ]
    a = supercell(lattice, P)
    b = StandardizedFromPrimitive(P)(a)
    @test a == b
    @test reciprocal(a) == reciprocal(b)
end

# Example from http://staff.ustc.edu.cn/~zqj/posts/Transforming-between-supercell-and-primitive-cell/
@testset "Test `supercell` from hexagonal to orthorhombic for a `Lattice`" begin
    lattice = Lattice(2, 2, 3, 90, 90, 120; axis=:c)  # Primitive is hexagonal
    P = [
        2 1 0  # 𝐨₁ = 2𝐡₁ + 𝐡₂
        0 1 0  # 𝐨₂ = 𝐡₂
        0 0 1  # 𝐨₃ = 𝐡₃
    ]'  # Note the transpose!
    a = supercell(lattice, P)  # Hexagonal to orthorhombic
    b = Lattice([
        2√3 0 0
        0 2 0
        0 0 3
    ])  # The orthorhombic supercell
    @test a ≈ b
    @test reciprocal(a) ≈ reciprocal(b)
    transform = StandardizedFromPrimitive(P)
    c = transform(lattice)
    @test a ≈ c
    @test reciprocal(a) ≈ reciprocal(c)
    @test inv(transform)(c) ≈ lattice
end
