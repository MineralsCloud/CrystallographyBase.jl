export ReferenceAtom, formation_energy

"""
    ReferenceAtom(label, energy)

An atom label paired with its reference energy per atom.

# Arguments
- `label`: The atomic label (e.g. a `Symbol`, `String`, or custom type).
- `energy`: The reference energy per atom.
"""
struct ReferenceAtom{L,E}
    label::L
    energy::E
end

"""
    formation_energy(cell::Cell, energy, component_energies)

Calculate the formation energy per atom for a `Cell`.

The formation energy per atom is:

```math
E_f = \\frac{E_{\\text{total}} - \\sum_i E_i}{N}
```

where ``E_{\\text{total}}`` is the total energy of the compound, ``E_i`` is the
reference energy per atom for the ``i``-th atom in the cell, and ``N`` is the
total number of atoms.

# Arguments
- `cell`: A `Cell` whose `atoms` field provides the atom labels.
- `energy`: Total energy of the compound.
- `component_energies`: A mapping from atom label to reference energy per atom.
"""
function formation_energy(cell::Cell, energy, component_energies)
    ref_sum = sum(get(component_energies, atom, missing) for atom in cell.atoms)
    return (energy - ref_sum) / natoms(cell)
end

"""
    formation_energy(cell::Cell{L,P,<:ReferenceAtom}, energy)

Calculate the formation energy per atom when each atom in `cell` is a
[`ReferenceAtom`](@ref), so that reference energies are embedded in the cell
itself and no separate `component_energies` mapping is required.

# Arguments
- `cell`: A `Cell` whose `atoms` are [`ReferenceAtom`](@ref) values.
- `energy`: Total energy of the compound.
"""
function formation_energy(cell::Cell{L,P,<:ReferenceAtom}, energy) where {L,P}
    ref_sum = sum(atom.energy for atom in cell.atoms)
    return (energy - ref_sum) / natoms(cell)
end
