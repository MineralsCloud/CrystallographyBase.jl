using StaticArrays: SHermitianCompact, SDiagonal

export MetricTensor, distance

struct MetricTensor{T} <: AbstractMatrix{T}
    data::SHermitianCompact{3,T,6}
end
MetricTensor(data::AbstractMatrix) = MetricTensor(SHermitianCompact{3}(data))
"""
    MetricTensor(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector)

Generate a `MetricTensor` from the three basis vectors.
"""
function MetricTensor(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector)
    𝐚𝐛𝐜 = (𝐚, 𝐛, 𝐜)
    return MetricTensor([dot(𝐱, 𝐲) for 𝐱 in 𝐚𝐛𝐜, 𝐲 in 𝐚𝐛𝐜])
end
"""
    MetricTensor(lattice::Lattice)

Generate a `MetricTensor` from a `Lattice`.
"""
function MetricTensor(lattice::AbstractLattice)
    data = parent(lattice)
    return MetricTensor(transpose(data) * data)
end
"""
    MetricTensor(a, b, c, α, β, γ)

Generate a `MetricTensor` from the six cell parameters.
"""
function MetricTensor(a, b, c, α, β, γ)
    g₁₂ = a * b * cosd(γ)
    g₁₃ = a * c * cosd(β)
    g₂₃ = b * c * cosd(α)
    return MetricTensor(SHermitianCompact(SVector(a^2, g₁₂, g₁₃, b^2, g₂₃, c^2)))
end

(g::MetricTensor)(𝐚::AbstractVector) = sqrt(dot(𝐚, g, 𝐚))
(g::MetricTensor)(𝐚::AbstractVector, 𝐛::AbstractVector) = g(𝐚 - 𝐛)

"""
    distance(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector)

Get the distance between two coordinates using a `MetricTensor`.
"""
distance(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector) = g(𝐚, 𝐛)

"""
    Lattice(g::MetricTensor)

Construct a `Lattice` from a `MetricTensor`.
"""
Lattice(g::MetricTensor) = Lattice(latticeconstants(g))

"""
    latticeconstants(g::MetricTensor)

Get the six lattice constants from a `MetricTensor`.
"""
function latticeconstants(g::MetricTensor)
    a², b², c², ab, ac, bc = g[1, 1], g[2, 2], g[3, 3], g[1, 2], g[1, 3], g[2, 3]
    a, b, c = map(sqrt, (a², b², c²))
    γ, β, α = acosd(ab / (a * b)), acosd(ac / (a * c)), acosd(bc / (b * c))
    return a, b, c, α, β, γ
end

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L130-L131
Base.one(::Type{MetricTensor{T}}) where {T} =
    MetricTensor(SDiagonal(one(T), one(T), one(T)))
Base.one(g::MetricTensor) = one(typeof(g))

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L132-L133
Base.oneunit(::Type{MetricTensor{T}}) where {T} =
    MetricTensor(SDiagonal(oneunit(T), oneunit(T), oneunit(T)))
Base.oneunit(g::MetricTensor) = oneunit(typeof(g))

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L134-L135
Base.zero(::Type{MetricTensor{T}}) where {T} = MetricTensor(zeros(T, 3, 3))
Base.zero(lattice::MetricTensor) = zero(typeof(lattice))

Base.parent(g::MetricTensor) = g.data

Base.size(::MetricTensor) = (3, 3)

Base.getindex(g::MetricTensor, i::Int) = getindex(parent(g), i)

Base.IndexStyle(::Type{<:MetricTensor}) = IndexLinear()

Base.inv(g::MetricTensor) = MetricTensor(inv(parent(g)))
