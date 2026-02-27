using CrystallographyCore: Inverted
using StructEqualHash: @struct_equal_hash

import CrystallographyCore: basisvectors, natoms, atomtypes

export ShiftedLattice, ShiftedCell, shift

struct ShiftedLattice{T} <: AbstractLattice{T}
    original::Lattice{T}
    by::SVector{3,T}
end
function ShiftedLattice(lattice::Lattice, 𝐱::AbstractVector)
    T = Base.promote_eltype(lattice, 𝐱)
    return ShiftedLattice(convert(Lattice{T}, lattice), SVector{3,T}(𝐱))
end
@struct_equal_hash ShiftedLattice

basisvectors(shifted::ShiftedLattice) = basisvectors(parent(shifted))

shift(lattice::Lattice, 𝐱::AbstractVector) = ShiftedLattice(lattice, 𝐱)
shift(shifted::ShiftedLattice, 𝐱::AbstractVector) = shift(parent(shifted), shifted.by .+ 𝐱)
function shift(lattice::Lattice, x::Integer, y::Integer, z::Integer)
    𝐚, 𝐛, 𝐜 = basisvectors(lattice)
    return shift(lattice, x * 𝐚 + y * 𝐛 + z * 𝐜)
end
function shift(shifted::ShiftedLattice, x::Integer, y::Integer, z::Integer)
    𝐚, 𝐛, 𝐜 = basisvectors(shifted)
    return shift(shifted, x * 𝐚 + y * 𝐛 + z * 𝐜)
end

Base.parent(shifted::ShiftedLattice) = shifted.original

# See https://github.com/MineralsCloud/CrystallographyCore.jl/blob/d9b808c/src/transform.jl#L10C80-L10C80
(shifted::ShiftedLattice)(reduced::AbstractVector) = parent(shifted)(reduced) .+ shifted.by

(inverted::Inverted{<:ShiftedLattice})(cartesian::AbstractVector) =
    inv(parent(inverted.lattice))(cartesian .- inverted.lattice.by)

mutable struct ShiftedCell{L,P,T} <: AbstractCell
    original::Cell{L,P,T}
    by::SVector{3,L}
end
@struct_equal_hash ShiftedCell

Base.parent(shifted::ShiftedCell) = shifted.original

natoms(shifted::ShiftedCell) = length(parent(shifted).atoms)

atomtypes(shifted::ShiftedCell) = unique(parent(shifted).atoms)

ShiftedLattice(shifted::ShiftedCell) = shift(Lattice(parent(shifted)), shifted.by)

function shift(cell::Cell{L}, 𝐱::AbstractVector{X}) where {L,X}
    T = Base.promote_type(L, X)
    return ShiftedCell(
        Cell(convert(Lattice{T}, Lattice(cell)), cell.positions, cell.atoms),
        SVector{3,T}(𝐱),
    )
end
shift(shifted::ShiftedCell, 𝐱::AbstractVector) = shift(parent(shifted), shifted.by .+ 𝐱)
function shift(cell::Cell, x::Integer, y::Integer, z::Integer)
    𝐚, 𝐛, 𝐜 = basisvectors(Lattice(cell))
    return shift(cell, x * 𝐚 + y * 𝐛 + z * 𝐜)
end
function shift(shifted::ShiftedCell, x::Integer, y::Integer, z::Integer)
    𝐚, 𝐛, 𝐜 = basisvectors(ShiftedLattice(shifted))
    return shift(shifted, x * 𝐚 + y * 𝐛 + z * 𝐜)
end

function Base.getproperty(shifted::ShiftedCell, name::Symbol)
    if name in (:positions, :atoms)
        return getproperty(parent(shifted), name)
    elseif name == :lattice
        return ShiftedLattice(shifted)
    else
        return getfield(shifted, name)
    end
end
