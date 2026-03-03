module SpglibExt

using CrystallographyBase: MagneticAtom

import CrystallographyCore: Cell
import Spglib: SpglibCell, unwrap_convert

function unwrap_convert(cell::Cell{L,P,<:MagneticAtom}) where {L,P}
    lattice, positions = cell.lattice, cell.positions
    clattice = Base.cconvert(Matrix{Cdouble}, permutedims(lattice))
    cpositions = Base.cconvert(Matrix{Cdouble}, reduce(hcat, positions))
    atoms = [atom.label for atom in cell.atoms]
    magmoms = [atom.magnetic_moment for atom in cell.atoms]
    atomtypes = unique(atoms)
    catoms = collect(Cint, findfirst(==(atom), atomtypes) for atom in atoms)  # Mapping between unique atom types and atom indices
    cmagmoms = if eltype(magmoms) <: AbstractVector
        Base.cconvert(Matrix{Cdouble}, reduce(hcat, magmoms))
    else
        Base.cconvert(Vector{Cdouble}, magmoms)  # `magmoms` could be empty!
    end
    return clattice, cpositions, catoms, cmagmoms
end

function SpglibCell(cell::Cell{L,P,<:MagneticAtom}) where {L,P}
    atoms = [atom.label for atom in cell.atoms]
    magmoms = [atom.magnetic_moment for atom in cell.atoms]
    return SpglibCell(cell.lattice, cell.positions, atoms, magmoms)
end

function Cell(cell::SpglibCell)
    if isempty(cell.magmoms)
        return Cell(cell.lattice, cell.positions, cell.atoms)  # Non-magnetic case
    else
        atoms = MagneticAtom.(cell.atoms, cell.magmoms)  # Magnetic case
        return Cell(cell.lattice, cell.positions, atoms)
    end
end

end
