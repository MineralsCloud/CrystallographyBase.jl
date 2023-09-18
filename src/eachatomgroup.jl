export eachatomgroup

struct AtomGroup{A,B}
    atom::A
    positions::Vector{B}
end
function eachatomgroup(cell::Cell)
    types = atomtypes(cell)
    return Iterators.map(types) do type
        indices = findall(==(type), cell.atoms)
        AtomGroup(type, cell.positions[indices])
    end
end
