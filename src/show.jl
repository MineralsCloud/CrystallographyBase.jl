function Base.show(io::IO, ::MIME"text/plain", point::ReciprocalPoint)
    println(io, string(typeof(point)))
    println(io, " coord=", point.coord)
    print(io, " weight=", point.weight)
    return nothing
end
function Base.show(io::IO, ::MIME"text/plain", A::ChangeOfBasis)
    println(io, string(typeof(A)))
    for row in eachrow(A.tf)
        println(io, ' ', join(row, "  "))
    end
end
