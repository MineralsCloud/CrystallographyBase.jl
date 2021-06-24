var documenterSearchIndex = {"docs":
[{"location":"api/","page":"API","title":"API","text":"CurrentModule = CrystallographyBase","category":"page"},{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Pages = [\"api.md\"]\nDepth = 3","category":"page"},{"location":"api/#Lattice","page":"API","title":"Lattice","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"CrystalSystem\nTriclinic\nMonoclinic\nOrthorhombic\nTetragonal\nCubic\nTrigonal\nHexagonal\nCentering\nBaseCentering\nPrimitive\nBodyCentering\nFaceCentering\nRhombohedralCentering\nBaseCentering\nBravais\nLattice\ncentering\ncrystalsystem\nbasis_vectors\ncellparameters\nsupercell","category":"page"},{"location":"api/#CrystallographyBase.CrystalSystem","page":"API","title":"CrystallographyBase.CrystalSystem","text":"Represent one of the seven crystal systems.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.Triclinic","page":"API","title":"CrystallographyBase.Triclinic","text":"Triclinic()\n\nRepresent the triclinic system.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.Monoclinic","page":"API","title":"CrystallographyBase.Monoclinic","text":"Monoclinic()\n\nRepresent the monoclinic system.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.Orthorhombic","page":"API","title":"CrystallographyBase.Orthorhombic","text":"Orthorhombic()\n\nRepresent the orthorhombic system.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.Tetragonal","page":"API","title":"CrystallographyBase.Tetragonal","text":"Tetragonal()\n\nRepresent the tetragonal system.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.Cubic","page":"API","title":"CrystallographyBase.Cubic","text":"Cubic()\n\nRepresent the cubic system.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.Trigonal","page":"API","title":"CrystallographyBase.Trigonal","text":"Trigonal()\n\nRepresent the trigonal system.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.Hexagonal","page":"API","title":"CrystallographyBase.Hexagonal","text":"Hexagonal()\n\nRepresent the hexagonal system.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.Centering","page":"API","title":"CrystallographyBase.Centering","text":"Represent the centering types.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.BaseCentering","page":"API","title":"CrystallographyBase.BaseCentering","text":"BaseCentering{:A}()\nBaseCentering{:B}()\nBaseCentering{:C}()\n\nRepresent the base-centering.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.Primitive","page":"API","title":"CrystallographyBase.Primitive","text":"Primitive()\n\nRepresent no centering.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.BodyCentering","page":"API","title":"CrystallographyBase.BodyCentering","text":"BodyCentering()\n\nRepresent the body-centering.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.FaceCentering","page":"API","title":"CrystallographyBase.FaceCentering","text":"FaceCentering()\n\nRepresent the face-centering.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.RhombohedralCentering","page":"API","title":"CrystallographyBase.RhombohedralCentering","text":"RhombohedralCentering()\n\nRepresent the rhombohedral-centering of the hexagonal system.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.Bravais","page":"API","title":"CrystallographyBase.Bravais","text":"Bravais(a::CrystalSystem, b::Centering, obverse::Bool=true)\n\nRepresent a Bravais lattice type.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.Lattice","page":"API","title":"CrystallographyBase.Lattice","text":"Lattice(mat::AbstractMatrix)\n\nConstruct a Lattice from a matrix.\n\nnote: Note\nThe basis vectors of the matrix are stored as columns.\n\n\n\n\n\nLattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector)\n\nConstruct a Lattice from three basis vectors.\n\n\n\n\n\nLattice(cell::Cell)\n\nGet the lattice of a Cell.\n\n\n\n\n\nLattice(a, b, c, α, β, γ)\n\nConstruct a Lattice from the six cell parameters.\n\nThe convention we used here is that edge vector 𝐚 in the positive x-axis direction, edge vector 𝐛 in the x-y plane with positive y-axis component, and edge vector 𝐜 with positive z-axis component in the Cartesian-system. See Wikipedia.\n\n\n\n\n\nLattice(g::MetricTensor)\n\nConstruct a Lattice from a MetricTensor.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.centering","page":"API","title":"CrystallographyBase.centering","text":"centering(bravais::Bravais)\n\nGet the centering type of a Bravais type.\n\n\n\n\n\n","category":"function"},{"location":"api/#CrystallographyBase.crystalsystem","page":"API","title":"CrystallographyBase.crystalsystem","text":"crystalsystem(bravais::Bravais)\n\nGet the crystal system of a Bravais type.\n\n\n\n\n\ncrystalsystem(a, b, c, α, β, γ)\n\nGuess the crystal system from the six cell parameters.\n\n\n\n\n\ncrystalsystem(lattice::Lattice)\n\nGet the crystal system of a lattice.\n\n\n\n\n\n","category":"function"},{"location":"api/#Spglib.basis_vectors","page":"API","title":"Spglib.basis_vectors","text":"basis_vectors(lattice::Lattice)\n\nGet the three basis vectors from a lattice.\n\n\n\n\n\nbasis_vectors(lattice::ReciprocalLattice)\n\nGet the three basis vectors from a ReciprocalLattice.\n\n\n\n\n\n","category":"function"},{"location":"api/#CrystallographyBase.cellparameters","page":"API","title":"CrystallographyBase.cellparameters","text":"cellparameters(lattice::Lattice)\n\nGet the six cell parameters from a lattice.\n\n\n\n\n\ncellparameters(g::MetricTensor)\n\nGet the six cell parameters from a MetricTensor.\n\n\n\n\n\n","category":"function"},{"location":"api/#CrystallographyBase.supercell","page":"API","title":"CrystallographyBase.supercell","text":"supercell(lattice::Lattice, expansion::AbstractMatrix{<:Integer})\n\nAllow the supercell to be a tilted extension of cell.\n\n\n\n\n\nsupercell(lattice::Lattice, expansion::AbstractVector{<:Integer})\n\nReturn a supercell based on cell and expansion coefficients.\n\n\n\n\n\n","category":"function"},{"location":"api/#Reciprocal-space","page":"API","title":"Reciprocal space","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Note that we take 2pi as 1, not the solid-state physics convention.","category":"page"},{"location":"api/","page":"API","title":"API","text":"ReciprocalPoint\nReciprocalLattice\ninv\nreciprocal_mesh\ncoordinates\nweights","category":"page"},{"location":"api/#CrystallographyBase.ReciprocalPoint","page":"API","title":"CrystallographyBase.ReciprocalPoint","text":"ReciprocalPoint(x, y, z, w)\n\nRepresent a special point of the 3D Brillouin zone. Each of them has a weight w.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.ReciprocalLattice","page":"API","title":"CrystallographyBase.ReciprocalLattice","text":"ReciprocalLattice(mat::SMatrix)\n\nConstruct a ReciprocalLattice.\n\nwarning: Warning\nYou should not use this function directly, always use inv of a Lattice.\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.inv","page":"API","title":"Base.inv","text":"inv(lattice::Lattice)\ninv(lattice::ReciprocalLattice)\n\nGet the reciprocal of a Lattice or a ReciprocalLattice.\n\n\n\n\n\n","category":"function"},{"location":"api/#CrystallographyBase.reciprocal_mesh","page":"API","title":"CrystallographyBase.reciprocal_mesh","text":"reciprocal_mesh(cell::Cell, mesh, is_shift; kwargs...)\n\nList the ReciprocalPoints from the mesh of the reciprocal space of a Cell.\n\nArguments\n\ncell::Cell: the cell.\nmesh: a vector of three integers which specify the mesh numbers along reciprocal primitive axis.\nis_shift=falses(3): a vector of three elements specifying whether the mesh is shifted along the axis in half of adjacent mesh points irrespective of the mesh numbers. The allowed values are 0, 1, true, and false.\nis_time_reversal=true: whether to impose the time reversal symmetry on the mesh.\nsymprec=1e-5: distance tolerance in Cartesian coordinates to find crystal symmetry.\ncartesian=false: whether to return the reciprocal points in Cartesian coordinates.\nir_only=true: whether to return the reciprocal points only in the irreducible Brillouin zone.\n\n\n\n\n\n","category":"function"},{"location":"api/#CrystallographyBase.coordinates","page":"API","title":"CrystallographyBase.coordinates","text":"coordinates(arr::AbstractArray{<:ReciprocalPoint})\n\nGet the coordinates of an array of ReciprocalPoints.\n\n\n\n\n\n","category":"function"},{"location":"api/#CrystallographyBase.weights","page":"API","title":"CrystallographyBase.weights","text":"weights(arr::AbstractArray{<:ReciprocalPoint})\n\nGet the weights of an array of ReciprocalPoints.\n\n\n\n\n\n","category":"function"},{"location":"api/#Miller-and-Miller–Bravais-indices","page":"API","title":"Miller and Miller–Bravais indices","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Miller\nMillerBravais\nReciprocalMiller\nReciprocalMillerBravais\nfamily\n@m_str","category":"page"},{"location":"api/#CrystallographyBase.Miller","page":"API","title":"CrystallographyBase.Miller","text":"Miller(i, j, k)\n\nRepresent the Miller indices in the real space (crystal directions).\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.MillerBravais","page":"API","title":"CrystallographyBase.MillerBravais","text":"MillerBravais(i, j, k, l)\n\nRepresent the Miller–Bravais indices in the real space (crystal directions).\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.ReciprocalMiller","page":"API","title":"CrystallographyBase.ReciprocalMiller","text":"ReciprocalMiller(i, j, k)\n\nRepresent the Miller indices in the reciprocal space (planes).\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.ReciprocalMillerBravais","page":"API","title":"CrystallographyBase.ReciprocalMillerBravais","text":"ReciprocalMillerBravais(i, j, k, l)\n\nRepresent the Miller–Bravais indices in the reciprocal space (planes).\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.family","page":"API","title":"CrystallographyBase.family","text":"family(x::Union{Miller,MillerBravais,ReciprocalMiller,ReciprocalMillerBravais})\n\nList the all the directions/planes that are equivalent to x by symmetry.\n\n\n\n\n\n","category":"function"},{"location":"api/#CrystallographyBase.@m_str","page":"API","title":"CrystallographyBase.@m_str","text":"m_str(s)\n\nGenerate the Miller indices or Miller–Bravais indices quickly.\n\nExamples\n\njulia> m\"[-1, 0, 1]\"\n3-element Miller:\n -1\n  0\n  1\n\njulia> m\"<2, -1, -1, 3>\"\n4-element MillerBravais:\n  2\n -1\n -1\n  3\n\njulia> m\"(-1, 0, 1)\"\n3-element ReciprocalMiller:\n -1\n  0\n  1\n\njulia> m\"(1, 0, -1, 0)\"\n4-element ReciprocalMillerBravais:\n  1\n  0\n -1\n  0\n\n\n\n\n\n","category":"macro"},{"location":"api/#Metric-tensor","page":"API","title":"Metric tensor","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetricTensor\ndirectioncosine\ndirectionangle\ndistance\ninterplanar_spacing","category":"page"},{"location":"api/#CrystallographyBase.MetricTensor","page":"API","title":"CrystallographyBase.MetricTensor","text":"MetricTensor(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector)\n\nGenerate a MetricTensor from the three basis vectors.\n\n\n\n\n\nMetricTensor(lattice::Lattice)\n\nGenerate a MetricTensor from a Lattice.\n\n\n\n\n\nMetricTensor(a, b, c, α, β, γ)\n\nGenerate a MetricTensor from the six cell parameters.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.directioncosine","page":"API","title":"CrystallographyBase.directioncosine","text":"directioncosine(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector)\n\nGet the direction cosine of two vectors and a MetricTensor.\n\n\n\n\n\n","category":"function"},{"location":"api/#CrystallographyBase.directionangle","page":"API","title":"CrystallographyBase.directionangle","text":"directionangle(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector)\n\nGet the direction angle of two vectors and a MetricTensor.\n\n\n\n\n\n","category":"function"},{"location":"api/#CrystallographyBase.distance","page":"API","title":"CrystallographyBase.distance","text":"distance(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector)\n\nGet the distance between two coordinates using a MetricTensor.\n\n\n\n\n\n","category":"function"},{"location":"api/#CrystallographyBase.interplanar_spacing","page":"API","title":"CrystallographyBase.interplanar_spacing","text":"interplanar_spacing(𝐚::AbstractVector, g::MetricTensor)\n\nGet the interplanar spacing from a MetricTensor.\n\n\n\n\n\n","category":"function"},{"location":"api/#Transformations","page":"API","title":"Transformations","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"CartesianFromFractional\nFractionalFromCartesian\nPrimitiveFromStandardized\nStandardizedFromPrimitive","category":"page"},{"location":"api/#CrystallographyBase.CartesianFromFractional","page":"API","title":"CrystallographyBase.CartesianFromFractional","text":"CartesianFromFractional(lattice::Union{Lattice,ReciprocalLattice})\nCartesianFromFractional(a, b, c, α, β, γ)\n\nGet the transformation from fractional coordinates to Cartesian coordinates.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.FractionalFromCartesian","page":"API","title":"CrystallographyBase.FractionalFromCartesian","text":"FractionalFromCartesian(lattice::Union{Lattice,ReciprocalLattice})\nFractionalFromCartesian(a, b, c, α, β, γ)\n\nGet the transformation from Cartesian coordinates to fractional coordinates.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.PrimitiveFromStandardized","page":"API","title":"CrystallographyBase.PrimitiveFromStandardized","text":"PrimitiveFromStandardized(tf::AbstractMatrix)\n\nConstruct the transformation from a standardized cell to a primitive cell.\n\n\n\n\n\n","category":"type"},{"location":"api/#CrystallographyBase.StandardizedFromPrimitive","page":"API","title":"CrystallographyBase.StandardizedFromPrimitive","text":"StandardizedFromPrimitive(tf::AbstractMatrix)\n\nConstruct the transformation from a primitive cell to a standardized cell.\n\n\n\n\n\n","category":"type"},{"location":"api/#Others","page":"API","title":"Others","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"cellvolume","category":"page"},{"location":"api/#CrystallographyBase.cellvolume","page":"API","title":"CrystallographyBase.cellvolume","text":"cellvolume(a, b, c, α, β, γ)\n\nCalculates the cell volume from 6 cell parameters.\n\n\n\n\n\ncellvolume(l::Lattice)\ncellvolume(c::Cell)\n\nCalculates the cell volume from a Lattice or a Cell.\n\n\n\n\n\ncellvolume(g::MetricTensor)\n\nCalculates the cell volume from a MetricTensor.\n\n\n\n\n\n","category":"function"},{"location":"develop/#How-to-contribute","page":"Development","title":"How to contribute","text":"","category":"section"},{"location":"develop/#Download-the-project","page":"Development","title":"Download the project","text":"","category":"section"},{"location":"develop/","page":"Development","title":"Development","text":"Similar to section \"Installation\", run","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"julia> using Pkg\n\njulia> pkg\"dev CrystallographyBase\"","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"in Julia REPL.","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"Then the package will be cloned to your local machine at a path. On macOS, by default is located at ~/.julia/dev/CrystallographyBase unless you modify the JULIA_DEPOT_PATH environment variable. (See Julia's official documentation on how to do this.) In the following text, we will call it PKGROOT.","category":"page"},{"location":"develop/#instantiating","page":"Development","title":"Instantiate the project","text":"","category":"section"},{"location":"develop/","page":"Development","title":"Development","text":"Go to PKGROOT, start a new Julia session and run","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"julia> using Pkg; Pkg.instantiate()","category":"page"},{"location":"develop/#How-to-build-docs","page":"Development","title":"How to build docs","text":"","category":"section"},{"location":"develop/","page":"Development","title":"Development","text":"Usually, the up-to-state doc is available in here, but there are cases where users need to build the doc themselves.","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"After instantiating the project, go to PKGROOT, run (without the $ prompt)","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"$ julia --color=yes docs/make.jl","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"in your terminal. In a while a folder PKGROOT/docs/build will appear. Open PKGROOT/docs/build/index.html with your favorite browser and have fun!","category":"page"},{"location":"develop/#How-to-run-tests","page":"Development","title":"How to run tests","text":"","category":"section"},{"location":"develop/","page":"Development","title":"Development","text":"After instantiating the project, go to PKGROOT, run (without the $ prompt)","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"$ julia --color=yes test/runtests.jl","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"in your terminal.","category":"page"},{"location":"installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Pages = [\"installation.md\"]\nDepth = 5","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"To install this package, first, you need to install a julia executable from its official website. Version v1.0.0 and above is required. This package may not work on v0.7 and below.","category":"page"},{"location":"installation/#Installing-Julia","page":"Installation","title":"Installing Julia","text":"","category":"section"},{"location":"installation/#on-macOS","page":"Installation","title":"on macOS","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"If you are using a Mac, and have Homebrew installed, open Terminal.app and type:","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"brew cask install julia","category":"page"},{"location":"installation/#on-Linux","page":"Installation","title":"on Linux","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"On Linux, the best way to install Julia is to use the Generic Linux Binaries. The JILL script does this for you. Just run","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"bash -ci \"$(curl -fsSL https://raw.githubusercontent.com/abelsiqueira/jill/master/jill.sh)\"","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"installs Julia into $HOME/.local/bin. This script also has a Python version, JILL.py. It can also be used on macOS.","category":"page"},{"location":"installation/#Installing-the-package","page":"Installation","title":"Installing the package","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Now I am using macOS as a standard platform to explain the following steps:","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Open Terminal.app, and type julia to start a Julia session.\nRun\njulia> using Pkg; Pkg.update()\n\njulia> Pkg.add(\"CrystallographyBase\")\nand wait for its finish.\nRun\njulia> using CrystallographyBase\nand have fun!\nWhile using, please keep this Julia session alive. Restarting may recompile the package and cost some time.","category":"page"},{"location":"installation/#Reinstall","page":"Installation","title":"Reinstall","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"In the same Julia session, run\njulia> Pkg.rm(\"CrystallographyBase\"); Pkg.gc()\nPress ctrl+d to quit the current session. Start a new Julia session and repeat the above steps.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = CrystallographyBase","category":"page"},{"location":"#CrystallographyBase","page":"Home","title":"CrystallographyBase","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for CrystallographyBase.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Provides some basic types and methods for crystallography calculations. For more features, see Crystallography.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the Index for the complete list of documented functions and types.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The code is hosted on GitHub, with some continuous integration services to test its validity.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This repository is created and maintained by singularitti. You are very welcome to contribute.","category":"page"},{"location":"#Compatibility","page":"Home","title":"Compatibility","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Julia version: v1.3.0 to v1.6.1\nDependencies:\nCombinatorics.jl v0.7.0 and above\nCompat.jl v2.2.0 and above\nCoordinateTransformations.jl v0.5.1 and above\nCounters.jl v0.3.0 and above\nFunctors.jl v0.1.0 and above\nSpglib.jl v0.2.0 and above\nStaticArrays.jl v0.8.3 and above\nOS: macOS, Linux, Windows, and FreeBSD\nArchitecture: x86, x64, ARM","category":"page"},{"location":"#Manual-Outline","page":"Home","title":"Manual Outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"installation.md\",\n    \"develop.md\",\n    \"api.md\",\n]\nDepth = 3","category":"page"},{"location":"#main-index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
