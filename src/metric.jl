using StaticArrays: SHermitianCompact

export MetricTensor
export directioncosine, directionangle, distance, interplanar_spacing

struct MetricTensor{T} <: AbstractMatrix{T}
    data::SHermitianCompact{3,T,6}
end
MetricTensor(m::AbstractMatrix) = MetricTensor(SHermitianCompact{3}(m))
"""
    MetricTensor(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector)

Generate a `MetricTensor` from the three basis vectors.
"""
function MetricTensor(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector)
    vecs = (𝐚, 𝐛, 𝐜)
    return MetricTensor([dot(vᵢ, vⱼ) for vᵢ in vecs, vⱼ in vecs])
end
"""
    MetricTensor(lattice::Lattice)

Generate a `MetricTensor` from a `Lattice`.
"""
function MetricTensor(lattice::Lattice)
    data = lattice.data
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
@functor MetricTensor

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

"""
    directioncosine(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector)

Get the direction cosine of two vectors and a `MetricTensor`.
"""
directioncosine(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector) =
    dot(𝐚, g, 𝐛) / (norm(𝐚, g) * norm(𝐛, g))

"""
    directionangle(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector)

Get the direction angle of two vectors and a `MetricTensor`.
"""
directionangle(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector) =
    acosd(directioncosine(𝐚, g, 𝐛))

"""
    distance(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector)

Get the distance between two coordinates using a `MetricTensor`.
"""
distance(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector) = norm(𝐛 - 𝐚, g)

"""
    interplanar_spacing(𝐚::AbstractVector, g::MetricTensor)

Get the interplanar spacing from a `MetricTensor`.
"""
interplanar_spacing(𝐚::AbstractVector, g::MetricTensor) = inv(norm(𝐚, g))

Base.size(::MetricTensor) = (3, 3)

Base.IndexStyle(::Type{<:MetricTensor}) = IndexLinear()

Base.getindex(g::MetricTensor, I::Vararg) = getindex(g.data, I...)

Base.inv(g::MetricTensor) = MetricTensor(inv(g.data))

dot(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector) = 𝐚' * g.data * 𝐛
norm(𝐚::AbstractVector, g::MetricTensor) = sqrt(dot(𝐚, g, 𝐚))
