export formation_energy

"""
    formation_energy(cell, energy, component_energies)

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

# Returns
Formation energy per atom.
"""
function formation_energy(cell::Cell, energy, component_energies)
    ref_sum = sum(component_energies[atom] for atom in cell.atoms)
    return (energy - ref_sum) / natoms(cell)
end
