module SpglibExt

using CrystallographyBase: Cell, MagneticAtom

import Spglib: unwrap_convert

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

end
