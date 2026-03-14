using CrystallographyCore: AbstractCell, EachAtom

import CrystallographyCore: eachatom

export eachatomgroup

"""
    AtomGroup(atom, indices)

Container for one atom type and the indices of all atoms of that type in a cell.

# Fields
- `atom`: The atom label/value shared by the group (for example `:Si`, `"O"`,
    or any custom atom type).
- `indices`: An `NTuple` of integer indices into `cell.atoms` or `cell.positions`
    identifying the atoms belonging to this group.

Use `cell.positions[group.indices]` to obtain the corresponding fractional coordinates.

See also [`eachatomgroup`](@ref), [`eachatom`](@ref).
"""
struct AtomGroup{N,T}
    atom::T
    indices::NTuple{N,Int64}
end

function eachatom(group::AtomGroup, cell::AbstractCell)
    return EachAtom(
        ntuple(Returns(group.atom), length(group.indices)), cell.positions[group.indices]
    )
end

function eachatomgroup(cell::AbstractCell)
    types = atomtypes(cell)
    return Iterators.map(types) do type
        indices = findall(==(type), cell.atoms)
        AtomGroup(type, Tuple(cell.positions[indices]))
    end
end
