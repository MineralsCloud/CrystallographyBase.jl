using StructEquality: @struct_hash_equal_isequal_isapprox

import CrystallographyCore: basisvectors

export ShiftedLattice

@struct_hash_equal_isequal_isapprox struct ShiftedLattice{T} <: AbstractLattice{T}
    original::Lattice{T}
    by::SVector{3,T}
end

basisvectors(shifted::ShiftedLattice) = basisvectors(parent(shifted))

function shift(lattice::Lattice, 𝐱::AbstractVector)
    T = Base.promote_eltype(lattice, 𝐱)
    return ShiftedLattice(convert(Lattice{T}, lattice), SVector{3,T}(𝐱))
end
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
