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
