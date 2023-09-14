export StandardizedFromPrimitive,
    PrimitiveFromStandardized, PrimitiveToStandardized, StandardizedToPrimitive

abstract type ChangeOfBasis{T} <: AbstractMatrix{T} end

# Idea from https://spglib.github.io/spglib/definition.html#transformation-to-the-primitive-cell
struct StandardizedFromPrimitive{T} <: ChangeOfBasis{T}
    tf::MMatrix{3,3,T,9}
end
struct PrimitiveFromStandardized{T} <: ChangeOfBasis{T}
    tf::MMatrix{3,3,T,9}
end
"""
    PrimitiveFromStandardized(tf::AbstractMatrix)

Construct the transformation from a standardized cell to a primitive cell.
"""
PrimitiveFromStandardized(tf::AbstractMatrix) = PrimitiveFromStandardized{eltype(tf)}(tf)
"""
    StandardizedFromPrimitive(tf::AbstractMatrix)

Construct the transformation from a primitive cell to a standardized cell.
"""
StandardizedFromPrimitive(tf::AbstractMatrix) = StandardizedFromPrimitive{eltype(tf)}(tf)
const PrimitiveToStandardized = StandardizedFromPrimitive
const StandardizedToPrimitive = PrimitiveFromStandardized

(x::Union{StandardizedFromPrimitive,PrimitiveFromStandardized})(v) = inv(x.tf) * collect(v)
(x::Union{StandardizedFromPrimitive,PrimitiveFromStandardized})(lattice::Lattice) =
    Lattice(parent(lattice) * x.tf)

Base.inv(x::StandardizedFromPrimitive) = PrimitiveFromStandardized(inv(x.tf))
Base.inv(x::PrimitiveFromStandardized) = StandardizedFromPrimitive(inv(x.tf))
Base.:∘(x::PrimitiveFromStandardized, y::StandardizedFromPrimitive) = ∘(y, x)
Base.:∘(x::StandardizedFromPrimitive, y::PrimitiveFromStandardized) =
    x.tf * y.tf ≈ I ? identity : error("undefined!")

Base.size(::ChangeOfBasis) = (3, 3)

Base.getindex(x::ChangeOfBasis, i) = getindex(x.tf, i)

Base.IndexStyle(::Type{<:ChangeOfBasis}) = IndexLinear()
