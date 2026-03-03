module SpglibExt

using CrystallographyBase: MagneticAtom

import CrystallographyCore: Cell
import Spglib: SpglibCell, is_spin_collinear, unwrap_convert

is_spin_collinear(cell::Cell{L,P,<:MagneticAtom}) where {L,P} =
    is_spin_collinear(SpglibCell(cell))

unwrap_convert(cell::Cell{L,P,<:MagneticAtom}) where {L,P} =
    unwrap_convert(SpglibCell(cell))

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
