function Base.show(io::IO, ::MIME"text/plain", path::ReciprocalPath)
    println(io, summary(path))
    println(io, " from: ", path.start_node)
    println(io, " to: ", path.end_node)
    println(io, " number of points: ", path.density)
    return nothing
end
