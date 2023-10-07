export cellvolume, crystaldensity, atomicmass

"""
    cellvolume(a, b, c, α, β, γ)

Calculate the cell volume from 6 cell parameters.
"""
cellvolume(a, b, c, α, β, γ) =
    a * b * c * sqrt(sind(α)^2 - cosd(β)^2 - cosd(γ)^2 + 2 * cosd(α) * cosd(β) * cosd(γ))
"""
    cellvolume(l::Lattice)
    cellvolume(c::Cell)

Calculate the cell volume from a `Lattice` or a `Cell`.
"""
cellvolume(lattice::AbstractLattice) = abs(_det(parent(lattice)))
cellvolume(cell::Cell) = cellvolume(Lattice(cell))
"""
    cellvolume(g::MetricTensor)

Calculate the cell volume from a `MetricTensor`.
"""
cellvolume(g::MetricTensor) = sqrt(_det(parent(g)))  # `sqrt` is always positive!

"""
    crystaldensity(lattice::Lattice, atoms)
    crystaldensity(cell::Cell)

Calculate the density of a crystal structure.

Here, `atoms` is an iterable of atomic numbers, element names, symbols, or `Mendeleev.Element`s.
You can extend the `atomicmass` method to work with custom types.
"""
function crystaldensity(lattice::Lattice, atoms)
    mass = sum(atomicmass, atoms)
    volume = cellvolume(lattice)
    return mass / volume
end
crystaldensity(cell::Cell) = crystaldensity(Lattice(cell), cell.atoms)

"""
    atomicmass(element::Element)
    atomicmass(i::Union{AbstractString,Integer,Symbol})

Return the atomic mass of an element.

!!! warning
    This function is by default not implemented. You need to load either package
    [`Mendeleev.jl`](https://github.com/Eben60/Mendeleev.jl) or
    [`PeriodicTable.jl`](https://github.com/JuliaPhysics/PeriodicTable.jl) to use it.
"""
function atomicmass end
