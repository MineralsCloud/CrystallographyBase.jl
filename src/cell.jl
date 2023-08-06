using CrystallographyCore: Cell, natoms, atomtypes, eachatom
using LinearAlgebra: isdiag, diag

export Cell, natoms, atomtypes, eachatom

"""
    supercell(cell::Cell, repfactors::AbstractMatrix{<:Integer})
    supercell(cell::Cell, repfactors::AbstractVector{<:Integer})
    supercell(cell::Cell, repfactor::Integer)

Create a supercell from `cell`.

!!! note
    Currently, only integral replications are supported.
"""
function supercell(cell::Cell, repfactors::AbstractMatrix{<:Integer})
    if size(repfactors) != (3, 3)
        throw(ArgumentError("`repfactors` must be a 3×3 matrix!"))
    end
    @assert isdiag(repfactors) "currently not supported!"
    @assert _det(repfactors) >= 1
    new_atoms = eltype(cell.atoms)[]
    new_positions = eltype(cell.positions)[]
    l, m, n = diag(repfactors)
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
    new_lattice = supercell(cell.lattice, repfactors)
    return Cell(new_lattice, new_positions, new_atoms)
end
