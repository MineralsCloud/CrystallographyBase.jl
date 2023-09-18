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
    @assert _det(factors) >= 1
    new_atoms = eltype(cell.atoms)[]
    new_positions = eltype(cell.positions)[]
    l, m, n = diag(factors)
    𝐚, 𝐛, 𝐜 = eachcol(Matrix(I, 3, 3))
    for (atom, position) in eachatom(cell)
        for (i, j, k) in Iterators.product(0:(l - 1), 0:(m - 1), 0:(n - 1))
            push!(new_atoms, atom)
            # See https://doi.org/10.1186/s13321-016-0129-3
            new_position = position + i * 𝐚 + j * 𝐛 + k * 𝐜
            new_position ./= (l, m, n)  # Make them within the boundary of the cell
            push!(new_positions, new_position)
        end
    end
    new_lattice = super(cell.lattice, factors)
    return Cell(new_lattice, new_positions, new_atoms)
end

function shift(cell::Cell, 𝐱::AbstractVector)
    new_lattice = shift(Lattice(cell), 𝐱)
    return Cell(new_lattice, cell.positions, cell.atoms)
end
