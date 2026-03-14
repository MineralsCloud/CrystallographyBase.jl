using CrystallographyCore: AbstractCell, EachAtom

import CrystallographyCore: eachatom

export eachatomtype

"""
    AtomTypeGroup(atom, indices)

Container for one atom type and the indices of all atoms of that type in a cell.

Use `cell.positions[group.indices]` to obtain the corresponding fractional coordinates.

See also [`eachatomtype`](@ref), [`eachatom`](@ref).
"""
struct AtomTypeGroup{N,T}
    atom::T
    indices::NTuple{N,Int64}
end

function eachatom(group::AtomTypeGroup, cell::AbstractCell)
    return EachAtom(
        ntuple(Returns(group.atom), length(group.indices)), cell.positions[group.indices]
    )
end

"""
    eachatomtype(cell::AbstractCell)

Iterate over atom types in `cell`, grouped by atom value/type.

For each distinct value in `cell.atoms`, this function returns an
[`AtomTypeGroup`](@ref) with:
- `atom`: that distinct atom value,
- `indices`: indices into `cell.atoms` / `cell.positions` of all atoms of that type.

The output order follows `atomtypes(cell)` (first appearance order in `cell.atoms`).

# Examples
```jldoctest
julia> cell = Cell(rand(3, 3), [[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]], [:Si, :Si, :O])

julia> group = collect(eachatomtype(cell))[1]
AtomTypeGroup{2, Symbol}(:Si, (1, 2))

julia> cell.positions[collect(group.indices)]
2-element Array{ReducedCoordinates{Float64}, 1}:
 [0.0, 0.0, 0.0]
 [0.5, 0.5, 0.5]
```
"""
function eachatomtype(cell::AbstractCell)
    types = atomtypes(cell)
    return Iterators.map(types) do type
        indices = Tuple(findall(==(type), cell.atoms))
        AtomTypeGroup(type, indices)
    end
end
