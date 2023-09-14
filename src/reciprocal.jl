using LinearAlgebra: cross
using StaticArrays: Size

import CrystallographyCore: basisvectors

export ReciprocalLattice, MonkhorstPackGrid, reciprocal, basisvectors

"""
    ReciprocalLattice(data::MMatrix)

Construct a `ReciprocalLattice`.

!!! warning
    You should not use this function directly, always use `reciprocal` of a `Lattice`.
"""
struct ReciprocalLattice{T} <: AbstractLattice{T}
    data::MMatrix{3,3,T,9}
end
ReciprocalLattice(data::AbstractMatrix) = ReciprocalLattice(MMatrix{3,3}(data))

Base.BroadcastStyle(::Type{<:ReciprocalLattice}) = Broadcast.ArrayStyle{ReciprocalLattice}()
Base.similar(
    bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{ReciprocalLattice}}, ::Type{S}
) where {S} = similar(ReciprocalLattice{S}, axes(bc))
ReciprocalLattice{S}(::UndefInitializer, dims) where {S} =
    ReciprocalLattice(Array{S,length(dims)}(undef, dims))

"""
    reciprocal(lattice::Lattice)
    reciprocal(lattice::ReciprocalLattice)

Get the reciprocal of a `Lattice` or a `ReciprocalLattice`.
"""
function reciprocal(lattice::Lattice)
    Ω = _det(lattice.data)  # Cannot use `cellvolume`, it takes the absolute value!
    𝐚, 𝐛, 𝐜 = basisvectors(lattice)
    return inv(Ω) * ReciprocalLattice(hcat(cross(𝐛, 𝐜), cross(𝐜, 𝐚), cross(𝐚, 𝐛)))
end
function reciprocal(lattice::ReciprocalLattice)
    Ω⁻¹ = _det(lattice.data)  # Cannot use `cellvolume`, it takes the absolute value!
    𝐚⁻¹, 𝐛⁻¹, 𝐜⁻¹ = eachbasisvector(lattice)
    return inv(Ω⁻¹) * Lattice(hcat(cross(𝐛⁻¹, 𝐜⁻¹), cross(𝐜⁻¹, 𝐚⁻¹), cross(𝐚⁻¹, 𝐛⁻¹)))
end

"""
    basisvectors(lattice::ReciprocalLattice)

Get the three basis vectors from a `ReciprocalLattice`.
"""
basisvectors(lattice::ReciprocalLattice) = Tuple(eachbasisvector(lattice))

eachbasisvector(lattice::ReciprocalLattice) = eachcol(lattice)

"""
    MonkhorstPackGrid(mesh, is_shift)

Represent the Monkhorst--Pack grid.

# Arguments
- `mesh`: A length-three vector specifying the k-point grid (``nk_1 × nk_2 × nk_3``) as in Monkhorst--Pack grids.
- `is_shift`: A length-three vector specifying whether the grid is displaced by half a grid step in the corresponding directions.
"""
struct MonkhorstPackGrid
    mesh::SVector{3,UInt}
    is_shift::SVector{3,Bool}
    function MonkhorstPackGrid(mesh, is_shift)
        @assert all(mesh .>= 1)
        if eltype(is_shift) != Bool
            is_shift = Bool.(is_shift)
        end
        return new(mesh, is_shift)
    end
end

Base.parent(lattice::ReciprocalLattice) = lattice.data

Base.getindex(lattice::ReciprocalLattice, i::Int) = getindex(parent(lattice), i)

Base.setindex!(lattice::ReciprocalLattice, v, i::Int) = setindex!(parent(lattice), v, i)

# Customizing broadcasting
# See https://github.com/JuliaArrays/StaticArraysCore.jl/blob/v1.4.2/src/StaticArraysCore.jl#L397-L398
# and https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/structuredbroadcast.jl#L7-L14
struct ReciprocalLatticeStyle <: Broadcast.AbstractArrayStyle{2} end
ReciprocalLatticeStyle(::Val{2}) = ReciprocalLatticeStyle()
ReciprocalLatticeStyle(::Val{N}) where {N} = Broadcast.DefaultArrayStyle{N}()

# Base.BroadcastStyle(::Type{<:ReciprocalLattice}) = ReciprocalLatticeStyle()

Base.similar(::Broadcast.Broadcasted{ReciprocalLatticeStyle}, ::Type{T}) where {T} =
    similar(Lattice{T})
# Override https://github.com/JuliaArrays/StaticArrays.jl/blob/v1.6.2/src/abstractarray.jl#L129
function Base.similar(lattice::ReciprocalLattice, ::Type{T}, _size::Size) where {T}
    if _size == size(lattice)
        ReciprocalLattice{T}(undef)
    else
        return similar(Array(lattice), T, _size)
    end
end
# Override https://github.com/JuliaLang/julia/blob/v1.10.0-beta2/base/abstractarray.jl#L839
function Base.similar(lattice::ReciprocalLattice, ::Type{T}, dims::Dims) where {T}
    if dims == size(lattice)
        ReciprocalLattice{T}(undef)
    else
        return similar(Array(lattice), T, dims)
    end
end
function Base.similar(::Type{<:ReciprocalLattice}, ::Type{T}, s::Size) where {T}
    if s == (3, 3)
        ReciprocalLattice{T}(undef)
    else
        return Array{T}(undef, Tuple(s))
    end
end
function Base.similar(::Type{<:ReciprocalLattice}, ::Type{T}, dim, dims...) where {T}
    if (dim, dims...) == (3, 3)
        ReciprocalLattice{T}(undef)
    else
        return Array{T}(undef, dim, dims...)
    end
end
