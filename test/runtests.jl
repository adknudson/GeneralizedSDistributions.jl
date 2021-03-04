using GeneralizedSDistributions
using Distributions
using Test

@testset "Implementation" begin

    D = Gamma(2, 2)
    G = GSDist(D)

    @testset "Parameter retrieval" begin
        @test_nowarn params(G)
    end

    @testset "Computation of statistics" begin
        @test_nowarn maximum(G)
        @test_nowarn minimum(G)
        @test_nowarn extrema(G)
        @test_nowarn mean(G)
        @test_nowarn var(G)
        @test_nowarn std(G)
        @test_nowarn median(G)
        @test_nowarn modes(G)
        @test_nowarn mode(G)

        @test_throws MethodError skewness(G)

        @test_throws MethodError kurtosis(G)
        @test_throws MethodError kurtosis(G, true)
        @test_throws MethodError kurtosis(G, false)

        @test_throws MethodError isplatykurtic(G)
        @test_throws MethodError isleptokurtic(G)
        @test_throws MethodError ismesokurtic(G)

        @test_throws MethodError entropy(G)

        @test_throws MethodError mgf(G, 3)

        @test_throws MethodError cf(G, 3)
    end

    @testset "Probability evaluation" begin
        x = mode(G)
        y = mode(G) + 1

        @test_nowarn insupport(G, x)

        @test_nowarn pdf(G, x)
        @test_nowarn logpdf(G, x)
        @test_nowarn loglikelihood(G, x)

        @test_nowarn cdf(G, x)
        @test_nowarn logcdf(G, x)
        @test_nowarn logdiffcdf(G, y, x)
        @test_nowarn ccdf(G, x)
        @test_nowarn logccdf(G, x)

        @test_nowarn quantile(G, 0.3)
        @test_nowarn cquantile(G, 0.3)

        @test_nowarn invlogcdf(G, log(0.3))
        @test_nowarn invlogccdf(G, log(0.3))
    end

    @testset "Sampling" begin
        @test_nowarn sampler(G)
        @test_nowarn rand(G, 10)
    end

end
