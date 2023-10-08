using CrystallographyCore: AbstractCell, EachAtom

import CrystallographyCore: eachatom

export eachatomgroup

struct AtomGroup{N,A,B}
    atom::A
    positions::NTuple{N,B}
end

eachatom(group::AtomGroup) =
    EachAtom(ntuple(Returns(group.atom), length(group.positions)), group.positions)

function eachatomgroup(cell::AbstractCell)
    types = atomtypes(cell)
    return Iterators.map(types) do type
        indices = findall(==(type), cell.atoms)
        AtomGroup(type, Tuple(cell.positions[indices]))
    end
end
