using EnumX: @enumx

export CrystalSystem, LatticeSystem, Bravais

"Represent the 7 lattice systems."
@enumx LatticeSystem begin
    Triclinic = 1
    Monoclinic = 2
    Orthorhombic = 3
    Tetragonal = 4
    Rhombohedral = 5
    Hexagonal = 6
    Cubic = 7
end

"Represent the 7 crystal systems."
@enumx CrystalSystem begin
    Triclinic = 1
    Monoclinic = 2
    Orthorhombic = 3
    Tetragonal = 4
    Trigonal = 5
    Hexagonal = 6
    Cubic = 7
end

"Represent the 14 Bravais lattices."
@enumx BravaisArithmeticClass begin
    PrimitiveTriclinic = 1
    PrimitiveMonoclinic = 2
    BaseCenteredMonoclinic = 3
    PrimitiveOrthorhombic = 4
    BaseCenteredOrthorhombic = 5
    BodyCenteredOrthorhombic = 6
    FaceCenteredOrthorhombic = 7
    PrimitiveTetragonal = 8
    BodyCenteredTetragonal = 9
    PrimitiveHexagonal = 10
    PrimitiveRhombohedral = 11
    RCentredHexagonal = PrimitiveRhombohedral
    PrimitiveCubic = 12
    BodyCenteredCubic = 13
    FaceCenteredCubic = 14
end
const Bravais = BravaisArithmeticClass
