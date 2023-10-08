using CrystallographyCore: AbstractLattice
using LinearAlgebra: Diagonal, norm

import CrystallographyCore: Lattice

export isrighthanded,
    islefthanded, latticesystem, latticeconstants, periodicity, super, shift

"""
    Lattice(a, b, c, α, β, γ; axis = :a)

Construct a `Lattice` from the six cell parameters.

The default convention we used here is that edge vector 𝐚 in the positive x-axis direction,
edge vector 𝐛 in the x-y plane with a positive y-axis component,
and edge vector 𝐜 with a positive z-axis component in the Cartesian system.
See [Wikipedia](https://en.wikipedia.org/w/index.php?title=Fractional_coordinates&oldid=961675499#In_crystallography).
You can also choose `axis = :c`.
"""
function Lattice(a, b, c, α, β, γ; axis=:a)
    Ω = cellvolume(a, b, c, α, β, γ)
    if axis == :a  # See https://en.wikipedia.org/w/index.php?title=Fractional_coordinates&oldid=961675499#In_crystallography
        sinγ, cosγ, cosα, cosβ, 𝟎 = sind(γ), cosd(γ), cosd(α), cosd(β), zero(a)
        return Lattice(
            [a, 𝟎, 𝟎],
            [b * cosγ, b * sinγ, zero(b)],
            [c * cosβ, c * (cosα - cosβ * cosγ) / sinγ, Ω / (a * b * sinγ)],
        )
    elseif axis == :c  # See https://github.com/LaurentRDC/crystals/blob/2d3a570/crystals/lattice.py#L356-L391
        sinα, cosα, sinβ, cosβ, 𝟎 = sind(α), cosd(α), sind(β), cosd(β), zero(c)
        x = Ω / (b * c * sinα)
        cos′ = (cosα * cosβ - cosd(γ)) / (sinα * sinβ)
        sin′ = sqrt(1 - cos′^2)
        return Lattice(
            [x, -x * cos′ / sin′, a * cosβ], [zero(b), b * sinα, b * cosα], [𝟎, 𝟎, c]
        )
    else
        error("aligning `$axis` axis is not supported!")
    end
end

"""
    isrighthanded(lattice::Lattice)

Test whether the basis vectors are defined to be right-handed.

The basis vectors are right-handed if and only if

```math
\\mathbf{a} \\cdot (\\mathbf{b} \\times \\mathbf{c}) > 0.
```
"""
function isrighthanded(lattice::Lattice)
    Δ = _det(parent(lattice))
    return Δ > zero(Δ)
end

"""
    islefthanded(lattice::Lattice)

Test whether the basis vectors are defined to be left-handed.

The basis vectors are left-handed if and only if

```math
\\mathbf{a} \\cdot (\\mathbf{b} \\times \\mathbf{c}) < 0.
```
"""
function islefthanded(lattice::Lattice)
    Δ = _det(parent(lattice))
    return Δ < zero(Δ)
end

"""
    latticesystem(bravais::Bravais)

Get the lattice system of a Bravais lattice.
"""
function latticesystem(bravais::Bravais)
    index = Int(bravais)
    if index == 1
        return LatticeSystem(:Triclinic)
    elseif 2 <= index <= 3
        return LatticeSystem(:Monoclinic)
    elseif 4 <= index <= 7
        return LatticeSystem(:Orthorhombic)
    elseif 8 <= index <= 9
        return LatticeSystem(:Tetragonal)
    elseif index == 10
        return LatticeSystem(:Hexagonal)
    elseif index == 11
        return LatticeSystem(:Rhombohedral)
    else  # 12 <= index <= 14
        return LatticeSystem(:Cubic)
    end
end
"""
    latticesystem(a, b, c, α, β, γ; kwargs...)

Guess the lattice system from the six lattice constants.

# Arguments
- `angletol=1e-5`: the absolute tolerance of angles (`α`, `β`, `γ`).
- `lengthtol=1e-5`: the absolute tolerance of edges (`a`, `b`, `c`).
"""
# See https://github.com/LaurentRDC/crystals/blob/2d3a570/crystals/lattice.py#L396-L475
function latticesystem(a, b, c, α, β, γ; angletol=1e-5, lengthtol=1e-5)
    lengths, angles = Base.vect(a, b, c), Base.vect(α, β, γ)
    ≊(θ, φ) = isapprox(θ, φ; atol=angletol)
    ≅(x, y) = isapprox(x, y; atol=lengthtol)
    unilength = all(x ≅ a for x in lengths)
    uniangle = all(θ ≊ α for θ in angles)
    function bilengths(vec)  # If and only if two lengths are equal.
        for x in vec
            if sum(isapprox(x, y; atol=lengthtol) for y in vec) == 2
                return true
            end
        end
        return false
    end
    # Checking for monoclinic system is generalized
    # to the case where a, b, and c can be cycled,
    # i.e., a != c && β != 90 && α == γ == 90
    #   || b != c && α != 90 && β == γ == 90
    #   || a != b && γ != 90 && α == β != 90
    for (lengths′, angles′) in zip(cyclic_perm(lengths), cyclic_perm(angles))
        (a′, _, c′), (α′, β′, γ′) = lengths′, angles′
        if !(a′ ≅ c′) && α′ ≊ 90 && γ′ ≊ 90 && !(β′ ≊ 90)
            return LatticeSystem(:Monoclinic)
        end
    end
    if unilength && uniangle
        return α ≊ 90 ? LatticeSystem(:Cubic) : LatticeSystem(:Rhombohedral)
    end
    if unilength && !uniangle  # Technically, a hexagonal system could have all 3 lengths equal.
        if any(θ ≊ 120 for θ in angles) && sum(θ ≊ 90 for θ in angles) == 2
            return LatticeSystem(:Hexagonal)
        end
    end
    if bilengths(lengths)  # At this point, two lengths are equal at most.
        if uniangle && α ≊ 90
            return LatticeSystem(:Tetragonal)
        elseif any(θ ≊ 120 for θ in angles) && sum(θ ≊ 90 for θ in angles) == 2
            return LatticeSystem(:Hexagonal)
        end
    else  # At this point, all lengths are not equal.
        return uniangle && α ≊ 90 ? LatticeSystem(:Orthorhombic) : LatticeSystem(:Triclinic)
    end
end
"""
    latticesystem(lattice::Lattice; angletol=1e-5, lengthtol=1e-5)

Get the lattice system of a `Lattice`.
"""
latticesystem(lattice::Lattice; kwargs...) =
    latticesystem(latticeconstants(lattice)...; kwargs...)
# Auxiliary functions
cyclic_perm(vec) = (circshift(vec, i) for i in 1:length(vec))  # See https://stackoverflow.com/a/43035441

"""
    latticeconstants(lattice::Lattice)

Get the six lattice constants from a `lattice`.
"""
function latticeconstants(lattice::AbstractLattice)  # Works for `ReciprocalLattice`s, too
    𝐚, 𝐛, 𝐜 = basisvectors(lattice)
    a, b, c = norm(𝐚), norm(𝐛), norm(𝐜)
    γ, β, α = acosd(dot(𝐚, 𝐛) / (a * b)),
    acosd(dot(𝐚, 𝐜) / (a * c)),
    acosd(dot(𝐛, 𝐜) / (b * c))
    return a, b, c, α, β, γ
end

# See https://github.com/LaurentRDC/crystals/blob/2d3a570/crystals/lattice.py#L161-L176
# Add the absolute value of the component of every lattice vector
# along the three euclidian vectors, which is effectively the sum of
# absolutes of rows.
"""
    periodicity(lattice::Lattice)

Get crystal periodicity in ``x``, ``y``, and ``z`` direction from the `Lattice`.
"""
periodicity(lattice::Lattice) = Tuple(sum(abs, parent(lattice); dims=2))

# See https://en.wikipedia.org/wiki/Supercell_(crystal)
"""
    super(lattice::Lattice, factors::AbstractMatrix{<:Integer})
    super(lattice::Lattice, factors::AbstractVector{<:Integer})
    super(lattice::Lattice, factor::Integer)

Create a supercell from `lattice`.
"""
function super(lattice::Lattice, factors::AbstractMatrix{<:Integer})
    if size(factors) != (3, 3)
        throw(ArgumentError("`factors` must be a 3×3 matrix!"))
    end
    # Sometimes the matrix can have negative determinant, see https://gitlab.com/ase/ase/-/issues/938
    @assert abs(_det(factors)) >= 1
    # See https://github.com/JuliaLang/julia/issues/51354
    return Lattice(parent(lattice) * factors)
end
super(lattice_or_cell, factors::AbstractVector{<:Integer}) =
    super(lattice_or_cell, Diagonal(factors))
# See https://stackoverflow.com/a/57270841
super(lattice_or_cell, factor::Integer) = super(lattice_or_cell, fill(factor, 3))
