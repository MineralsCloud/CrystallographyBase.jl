function Base.show(io::IO, ::MIME"text/plain", path::ReciprocalPath)
    println(io, summary(path))
    println(io, " from: ", path.start_node)
    println(io, " to: ", path.end_node)
    println(io, " number of points: ", path.density)
    return nothing
end
function Base.show(io::IO, ::MIME"text/plain", atom::MagneticAtom)
    println(io, summary(atom))
    println(io, " label=", atom.label)
    println(io, " magmom=", atom.magnetic_moment)
    return nothing
end
function Base.show(io::IO, atom::MagneticAtom)
    if get(io, :compact, false)
        # Inside a nested container or when compact printing is explicitly requested
        print(io, atom.label, "(", atom.magnetic_moment, ")")
    else
        # Standard one-line representation (should ideally be a valid Julia expression)
        print(io, atom.label, " → ", atom.magnetic_moment)
    end
end
