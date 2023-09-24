function Base.show(io::IO, ::MIME"text/plain", A::ChangeOfBasis)
    println(io, string(typeof(A)))
    for row in eachrow(A.tf)
        println(io, ' ', join(row, "  "))
    end
end
function Base.show(io::IO, ::MIME"text/plain", path::ReciprocalPath)
    println(io, summary(path))
    println(io, " Brillouinzone(:", path.bz, ')')
    print(
        io, " from ", path.start_node, " to ", path.end_node, " (", path.density, " points)"
    )
    return nothing
end
function Base.show(io::IO, ::MIME"text/plain", paths::ReciprocalPaths)
    println(io, summary(paths))
    println(io, " Brillouin zone: ", paths.bz)
    for path in eachpath(paths)
        println(
            io,
            " from ",
            path.start_node,
            " to ",
            path.end_node,
            " (",
            path.density,
            " points)",
        )
    end
    return nothing
end
