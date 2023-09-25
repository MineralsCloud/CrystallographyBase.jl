export MonkhorstPackGrid,
    BrillouinZone,
    ReciprocalPaths,
    DispersionRelation,
    BandStructure,
    PhononSpectrum,
    specialpoints,
    suggestedpath,
    interpolate,
    eachpoint,
    eachpath,
    eachchain

"""
    MonkhorstPackGrid(mesh, is_shift)

Represent the Monkhorst--Pack grid.

# Arguments
- `mesh`: A length-three vector specifying the k-point grid (``nk_1 × nk_2 × nk_3``) as in Monkhorst--Pack grids.
- `is_shift`: A length-three vector specifying whether the grid is displaced by half a grid step in the corresponding directions.
"""
struct MonkhorstPackGrid
    mesh::SVector{3,UInt}
    is_shift::SVector{3,Bool}
    function MonkhorstPackGrid(mesh, is_shift)
        @assert all(mesh .>= 1)
        if eltype(is_shift) != Bool
            is_shift = Bool.(is_shift)
        end
        return new(mesh, is_shift)
    end
end

@anonymousenum BrillouinZone::UInt8 begin
    # PrimitiveTriclinic = 1
    # PrimitiveMonoclinic = 2
    # BaseCenteredMonoclinic = 3
    # PrimitiveOrthorhombic = 4
    # BaseCenteredOrthorhombic = 5
    # BodyCenteredOrthorhombic = 6
    # FaceCenteredOrthorhombic = 7
    # PrimitiveTetragonal = 8
    # BodyCenteredTetragonal = 9
    PrimitiveHexagonal = 10
    # PrimitiveRhombohedral = 11
    PrimitiveCubic = 12
    BodyCenteredCubic = 13
    FaceCenteredCubic = 14
end

specialpoints(bz::BrillouinZone) = _specialpoints(Val(Int(bz)))
specialpoints(bz::Union{Integer,Symbol}) = specialpoints(BrillouinZone(bz))
_specialpoints(::Val{10}) = (
    Γ=ReducedCoordinates(zeros(Rational{Int64}, 3)),
    A=ReducedCoordinates(0, 0, 1//2),
    K=ReducedCoordinates(2//3, 1//3, 0),
    H=ReducedCoordinates(2//3, 1//3, 1//2),
    M=ReducedCoordinates(1//2, 0, 0),
    L=ReducedCoordinates(1//2, 0, 1//2),
)
_specialpoints(::Val{12}) = (
    Γ=ReducedCoordinates(zeros(Rational{Int64}, 3)),
    X=ReducedCoordinates(0, 1//2, 0),
    M=ReducedCoordinates(1//2, 1//2, 0),
    R=ReducedCoordinates(1//2, 1//2, 1//2),
)
_specialpoints(::Val{13}) = (
    Γ=ReducedCoordinates(zeros(Rational{Int64}, 3)),
    N=ReducedCoordinates(0, 1//2, 0),
    P=ReducedCoordinates(1//4, 1//4, 1//4),
    H=ReducedCoordinates(-1//2, 1//2, 1//2),
)
_specialpoints(::Val{14}) = (
    Γ=ReducedCoordinates(zeros(Rational{Int64}, 3)),
    X=ReducedCoordinates(0, 1//2, 1//2),
    L=ReducedCoordinates(1//2, 1//2, 1//2),
    W=ReducedCoordinates(1//4, 3//4, 1//2),
    U=ReducedCoordinates(1//4, 5//8, 5//8),
    K=ReducedCoordinates(3//8, 3//4, 3//8),
)

suggestedpath(bz::BrillouinZone) = _suggestedpath(Val(Int(bz)))
suggestedpath(bz::Union{Integer,Symbol}) = suggestedpath(BrillouinZone(bz))
_suggestedpath(::Val{10}) = (:Γ, :M, :K, :Γ, :A, :L, :H, :A), (:L, :M), (:K, :H)
_suggestedpath(::Val{12}) = (:Γ, :X, :M, :Γ, :R, :X), (:M, :R)
_suggestedpath(::Val{13}) = (:Γ, :H, :N, :Γ, :P, :H), (:P, :N)
_suggestedpath(::Val{14}) = (:Γ, :X, :W, :K, :Γ, :L, :U, :W, :L, :K), (:U, :X)

struct ReciprocalPath{T}
    start_node::ReducedCoordinates{T}
    end_node::ReducedCoordinates{T}
    density::Int64
end
function ReciprocalPath(bz::BrillouinZone, start_node::Symbol, end_node::Symbol, density)
    points = specialpoints(bz)
    return ReciprocalPath(points[start_node], points[end_node], density)
end

struct DispersionRelation{S,T}
    path::ReciprocalPath{S}
    bands::Matrix{T}
    function DispersionRelation{S,T}(path, bands) where {S,T}
        if length(interpolate(path)) != size(bands, 1)
            throw(
                DimensionMismatch(
                    "the number of interpolated reciprocal points and bands are different!"
                ),
            )
        end
        return new(path, bands)
    end
end
DispersionRelation(path::ReciprocalPath{S}, bands::AbstractMatrix{T}) where {S,T} =
    DispersionRelation{S,T}(path, bands)
const BandStructure = DispersionRelation
const PhononSpectrum = DispersionRelation

eachpoint(path::ReciprocalPath) = (point for point in interpolate(path))
eachpoint(dispersion::DispersionRelation) =
    zip(interpolate(dispersion.path), eachrow(dispersion.bands))

function interpolate(path::ReciprocalPath)
    iter = (
        range(aᵢ; stop=bᵢ, length=path.density) for
        (aᵢ, bᵢ) in zip(path.start_node, path.end_node)
    )
    return map(ReducedCoordinates, zip(iter...))
end
interpolate(dispersion::DispersionRelation) = collect(eachpoint(dispersion))
