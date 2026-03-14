using CrystallographyCore: AbstractCell, EachAtom, atomtypes

import CrystallographyCore: eachatom

export eachatomtype, each_equivalent_atom

"""
    AtomTypeGroup(atom, indices)

Container for one atom type and the indices of all atoms of that type in a cell.

Use `cell.positions[group.indices]` to obtain the corresponding fractional coordinates.

See also [`eachatomtype`](@ref), [`eachatom`](@ref).
"""
struct AtomTypeGroup{N,T}
    atom::T
    indices::SVector{N,Int64}
end
AtomTypeGroup(atom, indices::AbstractVector) =
    AtomTypeGroup(atom, SVector{length(indices)}(indices))

eachatom(group::AtomTypeGroup, cell::AbstractCell) =
    EachAtom(fill(group.atom, length(group.indices)), cell.positions[group.indices])

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

julia> cell.positions[group.indices]
...
 [0.0, 0.0, 0.0]
 [0.5, 0.5, 0.5]
```
"""
function eachatomtype(cell::AbstractCell)
    types = atomtypes(cell)
    return Iterators.map(types) do type
        indices = findall(==(type), cell.atoms)
        AtomTypeGroup(type, indices)
    end
end

"""
    EquivalentAtomGroup(atom, group, indices)

Represent a single equivalence class of atoms.

- `atom` is the atom value shared by the group.
- `group` is the equivalence label (group ID) taken directly from the input
  `equivalence` vector passed to [`each_equivalent_atom`](@ref).
  It is **not** renumbered or compacted; it is simply the original label.
- `indices` are the positions in `cell.atoms` / `cell.positions` that belong
  to this equivalence class.

!!! info
    Using `group` as a direct label allows callers to recover the original
    equivalence classification (e.g. symmetry operation index, order number,
    etc.) without remapping.
"""
struct EquivalentAtomGroup{T}
    atom::T
    group::Int64
    indices::Vector{Int64}
end

"""
    each_equivalent_atom(cell, equivalence)

Iterate over symmetry-equivalent atom groups.

`equivalence` must be a vector of integer labels of the same length as
`cell.atoms`. Atoms with the same label are considered equivalent and
are grouped together.

The returned iterator yields [`EquivalentAtomGroup`](@ref) objects.
"""
function each_equivalent_atom(cell::AbstractCell, equivalence::AbstractVector{<:Integer})
    if length(equivalence) != length(cell.atoms)
        throw(DimensionMismatch("`equivalence` must have the same length as `cell.atoms`"))
    end
    # Build groups of indices by equivalence label.
    # We maintain a deterministic output order by recording group labels in
    # the order they first appear in `equivalence`.
    groups = Vector{Pair{Int,Vector{Int}}}()
    group_index = Dict{Int,Int}()
    for (i, g) in enumerate(equivalence)
        if haskey(group_index, g)
            push!(groups[group_index[g]].second, i)
        else
            push!(groups, g => [i])
            group_index[g] = length(groups)
        end
    end
    return Iterators.map(groups) do (g, indices)
        # Ensure this equivalence group refers to a single atom label.
        # `only` throws if there is more than one distinct label.
        atom = only(unique(cell.atoms[indices]))
        EquivalentAtomGroup(atom, g, indices)
    end
end

"""
    eachatom(group::EquivalentAtomGroup, cell::AbstractCell)

Iterate over `(atom, position)` pairs for a given `EquivalentAtomGroup`.
"""
function eachatom(group::EquivalentAtomGroup, cell::AbstractCell)
    indices = group.indices
    positions = cell.positions[indices]
    return EachAtom(fill(group.atom, length(indices)), positions)
end
