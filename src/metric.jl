using StaticArrays: SHermitianCompact

export MetricTensor
export directioncosine, directionangle, distance, interplanar_spacing

struct MetricTensor{T} <: AbstractMatrix{T}
    data::SHermitianCompact{3,T,6}
end
MetricTensor(m::AbstractMatrix) = MetricTensor(SHermitianCompact{3}(m))
"""
    MetricTensor(ð::AbstractVector, ð::AbstractVector, ð::AbstractVector)

Generate a `MetricTensor` from the three basis vectors.
"""
function MetricTensor(ð::AbstractVector, ð::AbstractVector, ð::AbstractVector)
    vecs = (ð, ð, ð)
    return MetricTensor([dot(váµ¢, vâ±¼) for váµ¢ in vecs, vâ±¼ in vecs])
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
    MetricTensor(a, b, c, Î±, Î², Î³)

Generate a `MetricTensor` from the six cell parameters.
"""
function MetricTensor(a, b, c, Î±, Î², Î³)
    gââ = a * b * cosd(Î³)
    gââ = a * c * cosd(Î²)
    gââ = b * c * cosd(Î±)
    return MetricTensor(SHermitianCompact(SVector(a^2, gââ, gââ, b^2, gââ, c^2)))
end
@functor MetricTensor

"""
    Lattice(g::MetricTensor)

Construct a `Lattice` from a `MetricTensor`.
"""
Lattice(g::MetricTensor) = Lattice(cellparameters(g))

"""
    cellparameters(g::MetricTensor)

Get the six cell parameters from a `MetricTensor`.
"""
function cellparameters(g::MetricTensor)
    aÂ², bÂ², cÂ², ab, ac, bc = g[1, 1], g[2, 2], g[3, 3], g[1, 2], g[1, 3], g[2, 3]
    a, b, c = map(sqrt, (aÂ², bÂ², cÂ²))
    Î³, Î², Î± = acosd(ab / (a * b)), acosd(ac / (a * c)), acosd(bc / (b * c))
    return a, b, c, Î±, Î², Î³
end

"""
    directioncosine(ð::AbstractVector, g::MetricTensor, ð::AbstractVector)

Get the direction cosine of two vectors and a `MetricTensor`.
"""
directioncosine(ð::AbstractVector, g::MetricTensor, ð::AbstractVector) =
    dot(ð, g, ð) / (norm(ð, g) * norm(ð, g))

"""
    directionangle(ð::AbstractVector, g::MetricTensor, ð::AbstractVector)

Get the direction angle of two vectors and a `MetricTensor`.
"""
directionangle(ð::AbstractVector, g::MetricTensor, ð::AbstractVector) =
    acosd(directioncosine(ð, g, ð))

"""
    distance(ð::AbstractVector, g::MetricTensor, ð::AbstractVector)

Get the distance between two coordinates using a `MetricTensor`.
"""
distance(ð::AbstractVector, g::MetricTensor, ð::AbstractVector) = norm(ð - ð, g)

"""
    interplanar_spacing(ð::AbstractVector, g::MetricTensor)

Get the interplanar spacing from a `MetricTensor`.
"""
interplanar_spacing(ð::AbstractVector, g::MetricTensor) = inv(norm(ð, g))

Base.size(::MetricTensor) = (3, 3)

Base.IndexStyle(::Type{<:MetricTensor}) = IndexLinear()

Base.getindex(g::MetricTensor, I::Vararg) = getindex(g.data, I...)

Base.inv(g::MetricTensor) = MetricTensor(inv(g.data))

dot(ð::AbstractVector, g::MetricTensor, ð::AbstractVector) = ð' * g.data * ð
norm(ð::AbstractVector, g::MetricTensor) = sqrt(dot(ð, g, ð))
