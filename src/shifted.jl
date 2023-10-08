import CrystallographyCore: basisvectors

export ShiftedLattice

struct ShiftedLattice{T} <: AbstractLattice{T}
    original::Lattice{T}
    by::SVector{3,T}
end

# basisvectors(lattice::ShiftedLattice) = basisvectors(lattice.original) .+ Ref(lattice.by)

Base.parent(lattice::ShiftedLattice) = lattice.original

Base.size(::ShiftedLattice) = (3, 3)

Base.getindex(lattice::ShiftedLattice, i::Int) = getindex(parent(lattice), i)

Base.setindex!(lattice::ShiftedLattice, v, i::Int) = setindex!(parent(lattice), v, i)

Base.IndexStyle(::Type{<:ShiftedLattice}) = IndexLinear()

function Base.:*(lattice::ShiftedLattice, x::AbstractVector)
    return parent(lattice) * x + lattice.by
end
