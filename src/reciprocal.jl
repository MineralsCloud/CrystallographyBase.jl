using LinearAlgebra: cross
using StaticArrays: Size

import CrystallographyCore: basisvectors

export ReciprocalPoint,
    ReciprocalLattice, MonkhorstPackGrid, reciprocal, coordinates, weights, basisvectors

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
    return ReciprocalLattice(
        inv(Ω) * transpose(hcat(cross(𝐛, 𝐜), cross(𝐜, 𝐚), cross(𝐚, 𝐛)))
    )
end
function reciprocal(lattice::ReciprocalLattice)
    Ω⁻¹ = _det(lattice.data)  # Cannot use `cellvolume`, it takes the absolute value!
    𝐚⁻¹, 𝐛⁻¹, 𝐜⁻¹ = basisvectors(lattice)
    return Lattice(inv(Ω⁻¹) * hcat(cross(𝐛⁻¹, 𝐜⁻¹), cross(𝐜⁻¹, 𝐚⁻¹), cross(𝐚⁻¹, 𝐛⁻¹)))
end

"""
    basisvectors(lattice::ReciprocalLattice)

Get the three basis vectors from a `ReciprocalLattice`.
"""
basisvectors(lattice::ReciprocalLattice) = lattice[1, :], lattice[2, :], lattice[3, :]

"""
    ReciprocalPoint(x, y, z, w)

Represent a special point of the 3D Brillouin zone. Each of them has a weight `w`.
"""
struct ReciprocalPoint{T}
    coord::SVector{3,T}
    weight::Float64
end
ReciprocalPoint(coord::AbstractVector{T}, weight) where {T} =
    ReciprocalPoint{T}(SVector{3}(coord), weight)
ReciprocalPoint(x, y, z, w) = ReciprocalPoint(SVector(x, y, z), w)

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

"""
    coordinates(arr::AbstractArray{<:ReciprocalPoint})

Get the coordinates of an array of `ReciprocalPoint`s.
"""
coordinates(arr::AbstractArray{<:ReciprocalPoint}) = map(x -> x.coord, arr)

"""
    weights(arr::AbstractArray{<:ReciprocalPoint})

Get the weights of an array of `ReciprocalPoint`s.
"""
weights(arr::AbstractArray{<:ReciprocalPoint}) = map(x -> x.weight, arr)

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
