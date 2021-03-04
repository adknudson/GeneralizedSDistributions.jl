function Base.show(io::IO, ::MIME"text/plain", G::GSDist)
    α, g, k, γ = params(G)
    println(io, "Generalized S-Distribution ($(G.dist))")
    println(io, " α: $α")
    println(io, " g: $g")
    println(io, " k: $k")
    print(io, " γ: $γ")
end