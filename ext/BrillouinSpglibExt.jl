module BrillouinSpglibExt

# See https://github.com/thchr/Brillouin.jl/blob/v0.5.13/ext/BrillouinSpglibExt.jl
using Brillouin: cartesianize!, points, latticize
using CrystallographyBase: reciprocal, basisvectors
using CrystallographyCore: eachbasisvector
using Spglib: SpglibCell, get_dataset
using StaticArrays: SVector

import Brillouin: irrfbz_path

function irrfbz_path(cell::SpglibCell, symprec=1e-5)
    dataset = get_dataset(cell, symprec)
    basisvectors = map(SVector{3,Float64}, eachcol(dataset.std_lattice))
    kpaths = irrfbz_path(dataset.spacegroup_number, basisvectors)
    reciprocal_basis = basisvectors(reciprocal(cell.lattice))
    cartesianize!(kpaths)
    for (label, k) in kpaths.points
        kpaths.points[label] = dataset.std_rotation_matrix * k
    end
    return latticize(kpaths, reciprocal_basis)
end

end
