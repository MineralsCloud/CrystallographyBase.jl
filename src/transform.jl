export ReducedFromCartesian,
    CartesianFromReduced,
    ReducedToCartesian,
    CartesianToReduced,
    StandardizedFromPrimitive,
    PrimitiveFromStandardized,
    PrimitiveToStandardized,
    StandardizedToPrimitive

abstract type ChangeOfBasis{T} <: AbstractMatrix{T} end
struct CartesianFromReduced{T} <: ChangeOfBasis{T}
    tf::MMatrix{3,3,T,9}
end
struct ReducedFromCartesian{T} <: ChangeOfBasis{T}
    tf::MMatrix{3,3,T,9}
end
# This requires the a-vector is parallel to the Cartesian x-axis.
# See https://en.wikipedia.org/wiki/Fractional_coordinates
"""
    CartesianFromReduced(lattice::Union{Lattice,ReciprocalLattice})
    CartesianFromReduced(a, b, c, α, β, γ)

Get the transformation from fractional coordinates to Cartesian coordinates.
"""
CartesianFromReduced(lattice::Lattice) = CartesianFromReduced(parent(lattice))
CartesianFromReduced(lattice::ReciprocalLattice) =
    CartesianFromReduced(transpose(parent(lattice)))
function CartesianFromReduced(a, b, c, α, β, γ)
    Ω = cellvolume(a, b, c, α, β, γ)
    b_sinγ, b_cosγ = b .* (sind(γ), cosd(γ))
    return CartesianFromReduced(
        [
            a b_cosγ c*cosd(β)
            0 b_sinγ c*_auxiliary(α, β, γ)
            0 0 Ω/(a * b_sinγ)
        ]
    )
end
"""
    ReducedFromCartesian(lattice::Union{Lattice,ReciprocalLattice})
    ReducedFromCartesian(a, b, c, α, β, γ)

Get the transformation from Cartesian coordinates to fractional coordinates.
"""
ReducedFromCartesian(lattice::Lattice) = ReducedFromCartesian(inv(parent(lattice)))
ReducedFromCartesian(lattice::ReciprocalLattice) =
    ReducedFromCartesian(transpose(inv(parent(lattice))))
function ReducedFromCartesian(a, b, c, α, β, γ)
    Ω = cellvolume(a, b, c, α, β, γ)
    b_sinγ = b * sind(γ)
    return ReducedFromCartesian(
        [
            1/a -cotd(γ)/a -b * c * _auxiliary(β, α, γ)/Ω
            0 1/b_sinγ -a * c * _auxiliary(α, β, γ)/Ω
            0 0 a * b_sinγ/Ω
        ],
    )
end
const ReducedToCartesian = CartesianFromReduced
const CartesianToReduced = ReducedFromCartesian

# This is a helper function and should not be exported!
_auxiliary(α, β, γ) = (cosd(α) - cosd(β) * cosd(γ)) / sind(γ)

(x::Union{CartesianFromReduced,ReducedFromCartesian})(v) = x.tf * collect(v)

Base.inv(x::ReducedFromCartesian) = CartesianFromReduced(inv(x.tf))
Base.inv(x::CartesianFromReduced) = ReducedFromCartesian(inv(x.tf))
Base.:∘(x::CartesianFromReduced, y::ReducedFromCartesian) = ∘(y, x)
Base.:∘(x::ReducedFromCartesian, y::CartesianFromReduced) =
    x.tf * y.tf ≈ I ? identity : error("undefined!")

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
