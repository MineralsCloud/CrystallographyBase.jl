using LinearAlgebra: Diagonal
using Spglib: Cell

export CrystalSystem,
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Cubic,
    Trigonal,
    Hexagonal,
    Centering,
    BaseCentering,
    Primitive,
    BodyCentering,
    FaceCentering,
    RhombohedralCentering,
    BaseCentering,
    Bravais,
    PrimitiveTriclinic,
    PrimitiveMonoclinic,
    ACenteredMonoclinic,
    BCenteredMonoclinic,
    CCenteredMonoclinic,
    PrimitiveOrthorhombic,
    ACenteredOrthorhombic,
    BCenteredOrthorhombic,
    CCenteredOrthorhombic,
    BodyCenteredOrthorhombic,
    FaceCenteredOrthorhombic,
    PrimitiveTetragonal,
    BodyCenteredTetragonal,
    PrimitiveCubic,
    BodyCenteredCubic,
    FaceCenteredCubic,
    PrimitiveHexagonal,
    RCenteredHexagonal,
    Cell,
    Lattice
export centering, crystalsystem, basis_vectors, cellparameters

"Represent one of the seven crystal systems."
abstract type CrystalSystem end
"""
    Triclinic()

Represent the triclinic system.
"""
struct Triclinic <: CrystalSystem end
"""
    Monoclinic()

Represent the monoclinic system.
"""
struct Monoclinic <: CrystalSystem end
"""
    Orthorhombic()

Represent the orthorhombic system.
"""
struct Orthorhombic <: CrystalSystem end
"""
    Tetragonal()

Represent the tetragonal system.
"""
struct Tetragonal <: CrystalSystem end
"""
    Cubic()

Represent the cubic system.
"""
struct Cubic <: CrystalSystem end
"""
    Trigonal()

Represent the trigonal system.
"""
struct Trigonal <: CrystalSystem end
"""
    Hexagonal()

Represent the hexagonal system.
"""
struct Hexagonal <: CrystalSystem end

"Represent the centering types."
abstract type Centering end
"""
    Primitive()

Represent no centering.
"""
struct Primitive <: Centering end
"""
    BodyCentering()

Represent the body-centering.
"""
struct BodyCentering <: Centering end
"""
    FaceCentering()

Represent the face-centering.
"""
struct FaceCentering <: Centering end
"""
    RhombohedralCentering()

Represent the rhombohedral-centering of the hexagonal system.
"""
struct RhombohedralCentering <: Centering end
"""
    BaseCentering{:A}()
    BaseCentering{:B}()
    BaseCentering{:C}()

Represent the base-centering.
"""
struct BaseCentering{T} <: Centering end
const ACentering = BaseCentering{:A}
const BCentering = BaseCentering{:B}
const CCentering = BaseCentering{:C}

"""
    Bravais(a::CrystalSystem, b::Centering, obverse::Bool=true)

Represent a Bravais lattice type.
"""
struct Bravais{A<:CrystalSystem,B<:Centering}
    obverse::Bool
end
Bravais(a::CrystalSystem, b::Centering, obverse::Bool = true) =
    Bravais{typeof(a),typeof(b)}(obverse)

const PrimitiveTriclinic = Bravais{Triclinic,Primitive}
const PrimitiveMonoclinic = Bravais{Monoclinic,Primitive}
const ACenteredMonoclinic = Bravais{Monoclinic,ACentering}
const BCenteredMonoclinic = Bravais{Monoclinic,BCentering}
const CCenteredMonoclinic = Bravais{Monoclinic,CCentering}
const PrimitiveOrthorhombic = Bravais{Orthorhombic,Primitive}
const ACenteredOrthorhombic = Bravais{Orthorhombic,ACentering}
const BCenteredOrthorhombic = Bravais{Orthorhombic,BCentering}
const CCenteredOrthorhombic = Bravais{Orthorhombic,CCentering}
const BodyCenteredOrthorhombic = Bravais{Orthorhombic,BodyCentering}
const FaceCenteredOrthorhombic = Bravais{Orthorhombic,FaceCentering}
const PrimitiveTetragonal = Bravais{Tetragonal,Primitive}
const BodyCenteredTetragonal = Bravais{Tetragonal,BodyCentering}
const PrimitiveCubic = Bravais{Cubic,Primitive}
const BodyCenteredCubic = Bravais{Cubic,BodyCentering}
const FaceCenteredCubic = Bravais{Cubic,FaceCentering}
const PrimitiveHexagonal = Bravais{Hexagonal,Primitive}
const RCenteredHexagonal = Bravais{Hexagonal,RhombohedralCentering}

"Represent the real lattices and the reciprocal lattices."
abstract type AbstractLattice{T} end

struct Lattice{T} <: AbstractLattice{T}
    data::SMatrix{3,3,T,9}
end
"""
    Lattice(mat::AbstractMatrix)

Construct a `Lattice` from a matrix.

!!! note
    The basis vectors of the matrix are stored as columns.
"""
Lattice(mat::AbstractMatrix) = Lattice(SMatrix{3,3}(mat))
"""
    Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector)

Construct a `Lattice` from three basis vectors.
"""
Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector) = Lattice(hcat(𝐚, 𝐛, 𝐜))
Lattice(lattice::Lattice) = lattice
"""
    Lattice(cell::Cell)

Get the lattice of a `Cell`.
"""
Lattice(cell::Cell) = Lattice(cell.lattice)
"""
    Lattice(a, b, c, α, β, γ)

Construct a `Lattice` from the six cell parameters.
"""
function Lattice(a, b, c, α, β, γ)
    # From https://github.com/LaurentRDC/crystals/blob/dbb3a92/crystals/lattice.py#L321-L354
    v = cellvolume(1, 1, 1, α, β, γ)
    # reciprocal lattice
    a_recip = sind(α) / (a * v)
    csg = (cosd(α) * cosd(β) - cosd(γ)) / (sind(α) * sind(β))
    sg = sqrt(1 - csg^2)
    a1 = [1 / a_recip, -csg / sg / a_recip, cosd(β) * a]
    a2 = [0, b * sind(α), b * cosd(α)]
    a3 = [0, 0, c]
    return Lattice(a1, a2, a3)
end
@functor Lattice

"""
    basis_vectors(lattice::Lattice)

Get the three basis vectors from a `lattice`.
"""
basis_vectors(lattice::Lattice) = lattice[:, 1], lattice[:, 2], lattice[:, 3]

"""
    centering(bravais::Bravais)

Get the centering type of a Bravais type.
"""
centering(::Bravais{A,B}) where {A,B} = B()

"""
    crystalsystem(bravais::Bravais)

Get the crystal system of a Bravais type.
"""
crystalsystem(::Bravais{A,B}) where {A,B} = A()
"""
    crystalsystem(a, b, c, α, β, γ)

Guess the crystal system from the six cell parameters.
"""
function crystalsystem(a, b, c, α, β, γ)
    if a == b == c
        if α == β == γ
            α == 90 ? Cubic() : Trigonal()
        else
            α == β == 90 && γ == 120 ? Hexagonal() : Triclinic()
        end
    else
        if α == β == γ == 90
            a == b || a == c || b == c ? Tetragonal() : Orthorhombic()
        else
            α == β == 90 || β == γ == 90 || α == γ == 90 ? Monoclinic() : Triclinic()
        end
    end
end
"""
    crystalsystem(lattice::Lattice)

Get the crystal system of a `lattice`.
"""
crystalsystem(lattice::Lattice) = crystalsystem(cellparameters(lattice))

"""
    cellparameters(lattice::Lattice)

Get the six cell parameters from a `lattice`.
"""
function cellparameters(lattice::Lattice)
    𝐚, 𝐛, 𝐜 = basis_vectors(lattice)
    a, b, c = norm(𝐚), norm(𝐛), norm(𝐜)
    γ, β, α =
        acosd(dot(𝐚, 𝐛) / (a * b)), acosd(dot(𝐚, 𝐜) / (a * c)), acosd(dot(𝐛, 𝐜) / (b * c))
    return a, b, c, α, β, γ
end

# See https://en.wikipedia.org/wiki/Supercell_(crystal)
"""
    supercell(cell::Lattice, expansion::AbstractMatrix{<:Integer})

Allow the supercell to be a tilted extension of `cell`.
"""
function supercell(lattice::Lattice, expansion::AbstractMatrix)
    if any(!isinteger(x) for x in expansion)
        throw(ArgumentError("`expansion` must be an integer matrix!"))
    end
    @assert det(expansion) >= 1
    return Lattice(lattice.data * expansion)
end
"""
    supercell(cell::Lattice, expansion::AbstractVector{<:Integer})

Return a supercell based on `cell` and expansion coefficients.
"""
supercell(lattice::Lattice, expansion::AbstractVector) =
    supercell(lattice, Diagonal(expansion))
function supercell(cell::Cell, expansion) end

Base.iterate(lattice::AbstractLattice) = iterate(lattice.data)
Base.iterate(lattice::AbstractLattice, state) = iterate(lattice.data, state)

Base.eltype(::AbstractLattice{T}) where {T} = T

Base.length(::AbstractLattice) = 9

Base.size(::AbstractLattice) = (3, 3)
Base.size(::AbstractLattice, dim::Integer) = dim <= 2 ? 3 : 1

Base.IteratorSize(::Type{<:AbstractLattice}) = Base.HasShape{2}()

Base.axes(lattice::AbstractLattice, dim::Integer) = axes(lattice.data, dim)

Base.getindex(lattice::AbstractLattice, i) = getindex(lattice.data, i)
Base.getindex(lattice::AbstractLattice, I::Vararg) = getindex(lattice.data, I...)

Base.firstindex(::AbstractLattice) = 1

Base.lastindex(::AbstractLattice) = 9

for op in (:+, :-)
    @eval Base.broadcast(::typeof($op), lattice::AbstractLattice, number::Number) =
        Lattice(broadcast($op, lattice.data, number))
    @eval Base.broadcast(::typeof($op), number::Number, lattice::AbstractLattice) =
        broadcast($op, lattice, number)
end
for op in (:*, :/, ://)
    @eval Base.$op(lattice::AbstractLattice, number::Number) =
        Lattice(($op)(lattice.data, number))
    @eval Base.$op(number::Number, lattice::AbstractLattice) = ($op)(lattice, number)
end

function Base.show(io::IO, x::AbstractLattice)
    if get(io, :compact, false)
        print(io, string(typeof(x)), '(')
        print(io, x.data, ')')
    else
        println(io, string(nameof(typeof(x))))
        for row in eachrow(x.data)
            print(io, " ")
            println(io, row)
        end
    end
end
