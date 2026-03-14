# Public API

## Contents

```@contents
Pages = ["public.md"]
Depth = 3
```

## Index

```@index
Pages = ["public.md"]
```

## Public interface

### Lattice and Cell

```@docs
Lattice
Base.:*(::AbstractMatrix, ::Lattice)
Base.:*(::Lattice, ::AbstractMatrix)
isrighthanded
islefthanded
basisvectors(::Lattice)
latticesystem
latticeconstants
periodicity
Cell
ReferenceAtom
MagneticAtom
super
eachatom
natoms
atomtypes
eachatomtype
each_equivalent_atom
atomcounts
cellvolume
crystaldensity
atomicmass
```

#### Other properties

```@docs
formation_energy
ismagnetic
magnetization
```

### Reciprocal space

Note that we take ``2\pi`` as ``1``, not the solid-state physics convention.

```@docs
ReciprocalLattice
basisvectors(::ReciprocalLattice)
MonkhorstPackGrid
reciprocal
```

### Metric tensor

```@docs
MetricTensor
lengthof
distance
```

### Transformations

```@docs
PrimitiveFromStandardized
StandardizedFromPrimitive
```
