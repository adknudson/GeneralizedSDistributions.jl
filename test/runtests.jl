using GeneralizedSDistributions
using Distributions
using Test

@testset "GeneralizedSDistributions.jl" begin
    # Test distributions
    A = LogNormal(3, 1)
    B = NegativeBinomial(20, 0.2)
    
    # GeneralizedSDistributions using the median as initial values
    F = GSDist(A)
    G = GSDist(B)

    # The median for test/gsdist must be the same
    @test ==(map(D -> quantile(D, 0.5), [A, F])...)
    @test ==(map(D -> quantile(D, 0.5), [B, G])...)

    # Quantiles must be approximately the same
    for p in 0.05:0.1:0.95
        @test ≈(map(D -> quantile(D, p), [A, F])..., rtol=0.05)
        @test ≈(map(D -> quantile(D, p), [B, G])..., rtol=0.05)
    end
end
