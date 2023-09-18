using CrystallographyCore: EachAtom

import CrystallographyCore: eachatom

export eachatomgroup

struct AtomGroup{A,B}
    atom::A
    positions::Vector{B}
end

eachatom(group::AtomGroup) =
    EachAtom(fill(group.atom, length(group.positions)), group.positions)

function eachatomgroup(cell::Cell)
    types = atomtypes(cell)
    return Iterators.map(types) do type
        indices = findall(==(type), cell.atoms)
        AtomGroup(type, cell.positions[indices])
    end
end
