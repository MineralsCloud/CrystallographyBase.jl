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

"""
    eachatomgroup(cell::AbstractCell)

Iterate over atom groups in `cell`, grouped by atom value/type.

For each distinct value in `cell.atoms`, this function returns an
[`AtomGroup`](@ref) with:
- `atom`: that distinct atom value,
- `indices`: indices into `cell.atoms` / `cell.positions` of all atoms of that type.

The output order follows `atomtypes(cell)` (first appearance order in
`cell.atoms`).

# Examples
```jldoctest
julia> cell = Cell(rand(3, 3), [[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]], [:Si, :Si, :O])

julia> group = collect(eachatomgroup(cell))[1]
AtomGroup{2, Symbol}(:Si, (1, 2))

julia> cell.positions[collect(group.indices)]
2-element Array{ReducedCoordinates{Float64}, 1}:
 [0.0, 0.0, 0.0]
 [0.5, 0.5, 0.5]
```
"""
function eachatomgroup(cell::AbstractCell)
    types = atomtypes(cell)
    return Iterators.map(types) do type
        indices = Tuple(findall(==(type), cell.atoms))
        AtomGroup(type, indices)
    end
end
