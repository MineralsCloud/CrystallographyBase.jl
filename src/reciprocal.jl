export MonkhorstPackGrid

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
_specialpoints(::Val{10}) = (
    Γ=[0, 0, 0],
    A=[0, 0, 1//2],
    K=[2//3, 1//3, 0],
    H=[2//3, 1//3, 1//2],
    M=[1//2, 0, 0],
    L=[1//2, 0, 1//2],
)
_specialpoints(::Val{12}) =
    (Γ=[0, 0, 0], X=[0, 1//2, 0], M=[1//2, 1//2, 0], R=[1//2, 1//2, 1//2])
_specialpoints(::Val{13}) =
    (Γ=[0, 0, 0], N=[0, 1//2, 0], P=[1//4, 1//4, 1//4], H=[-1//2, 1//2, 1//2])
_specialpoints(::Val{14}) = (
    Γ=[0, 0, 0],
    X=[0, 1//2, 1//2],
    L=[1//2, 1//2, 1//2],
    W=[1//4, 3//4, 1//2],
    U=[1//4, 5//8, 5//8],
    K=[3//8, 3//4, 3//8],
)

suggestedpath(bz::BrillouinZone) = _suggestedpath(Val(Int(bz)))
_suggestedpath(::Val{10}) = (:Γ, :M, :K, :Γ, :A, :L, :H, :A), (:L, :M), (:K, :H)
_suggestedpath(::Val{12}) = (:Γ, :X, :M, :Γ, :R, :X), (:M, :R)
_suggestedpath(::Val{13}) = (:Γ, :H, :N, :Γ, :P, :H), (:P, :N)
_suggestedpath(::Val{14}) = (:Γ, :X, :W, :K, :Γ, :L, :U, :W, :L, :K), (:U, :X)
