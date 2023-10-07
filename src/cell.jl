using CrystallographyCore: Cell, eachatom
using LinearAlgebra: I, isdiag, diag

"""
    super(cell::Cell, factors::AbstractMatrix{<:Integer})
    super(cell::Cell, factors::AbstractVector{<:Integer})
    super(cell::Cell, factor::Integer)

Create a supercell from `cell`.

!!! note
    Currently, only integral replications are supported.
"""
function super(cell::Cell, factors::AbstractMatrix{<:Integer})
    if size(factors) != (3, 3)
        throw(ArgumentError("`factors` must be a 3×3 matrix!"))
    end
    @assert isdiag(factors) "currently not supported!"
    @assert all(factor >= 1 for factor in diag(factors)) "all factors must be greater than or equal to one!"
    l, m, n = diag(factors)
    𝐚, 𝐛, 𝐜 = eachcol(Matrix(I, 3, 3))
    new_atoms = eltype(cell.atoms)[]
    new_positions = eltype(cell.positions)[]
    for (atom, position) in eachatom(cell)
        for (i, j, k) in Iterators.product(0:(l - 1), 0:(m - 1), 0:(n - 1))
            # See https://doi.org/10.1186/s13321-016-0129-3 and #111
            new_position = (position + i * 𝐚 + j * 𝐛 + k * 𝐜) ./ (l, m, n)  # Make them within the boundary of the cell
            push!(new_positions, new_position)
            push!(new_atoms, atom)
        end
    end
    new_lattice = super(cell.lattice, factors)
    return Cell(new_lattice, new_positions, new_atoms)
end

function shift(cell::Cell, 𝐱::AbstractVector)
    new_lattice = shift(Lattice(cell), 𝐱)
    new_positions = Ref(ShiftedLattice(Lattice(cell), 𝐱)) .* cell.positions
    new_positions = Ref(new_lattice) .\ new_positions
    return Cell(new_lattice, new_positions, cell.atoms)
end
function shift(cell::Cell, x::Integer, y::Integer, z::Integer)
    𝐚, 𝐛, 𝐜 = basisvectors(Lattice(cell))
    return shift(cell, x * 𝐚 + y * 𝐛 + z * 𝐜)
end
