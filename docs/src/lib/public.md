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
isrighthanded
islefthanded
basisvectors
latticesystem
latticeconstants
periodicity
Cell
super
eachatom
```

### Reciprocal space

Note that we take ``2\pi`` as ``1``, not the solid-state physics convention.

```@docs
ReciprocalLattice
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

### Others

```@docs
cellvolume
crystaldensity
atomicmass
```
