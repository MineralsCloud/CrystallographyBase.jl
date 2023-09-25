function Base.show(io::IO, ::MIME"text/plain", A::ChangeOfBasis)
    println(io, string(typeof(A)))
    for row in eachrow(A.tf)
        println(io, ' ', join(row, "  "))
    end
end
function Base.show(io::IO, ::MIME"text/plain", path::ReciprocalPath)
    println(io, summary(path))
    println(io, " from: ", path.start_node)
    println(io, " to: ", path.end_node)
    println(io, " number of points: ", path.density)
    return nothing
end
