using CrystallographyCore: Cell, eachatom
using LinearAlgebra: I, isdiag, diag

export MagneticAtom, MagneticCell, ismagnetic, magnetization

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
    new_positions = float(eltype(cell.positions))[]
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

"""
    MagneticAtom(label, magnetic_moment)

An atom label paired with its initial magnetic moment.

# Arguments
- `label`: The atomic label (e.g., a `Symbol`, `String`, or custom type).
- `magnetic_moment`: The magnetic moment of the atom (scalar for collinear, vector for non-collinear).
"""
struct MagneticAtom{L,M}
    label::L
    magnetic_moment::M
end

const MagneticCell = Cell{L,P,<:MagneticAtom} where {L,P}

"""
    ismagnetic(atom::MagneticAtom, tol=1e-8)

Return `true` if the magnetic moment of `atom` is larger than `tol`.

The function is tolerant of small numerical noise and treats moments with magnitude
≤ `tol` as non-magnetic.
"""
ismagnetic(atom::MagneticAtom, tol=1e-8) = norm(atom.magnetic_moment) > tol
"""
    ismagnetic(cell::MagneticCell, tol=1e-8)

Return `true` if the total magnetic moment in `cell` (sum over all atoms) is
larger than `tol` in magnitude.
"""
ismagnetic(cell::MagneticCell, tol=1e-8) =
    norm(sum(atom.magnetic_moment for atom in cell.atoms)) > tol

"""
    magnetization(cell::Cell)

Return the magnetization (magnetic moment density) of a cell containing `MagneticAtom`s.

The magnetization is defined as the sum of all atomic magnetic moments divided by the cell volume.
"""
magnetization(cell::MagneticCell) =
    sum(atom.magnetic_moment for atom in cell.atoms) / cellvolume(cell)
