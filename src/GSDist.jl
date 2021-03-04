"""
    GSDist <: ContinuousUnivariateDistribution

Generalized S-Distribution for approximating univariate distributions.
"""
struct GSDist <: ContinuousUnivariateDistribution
    F₀::Real
    x₀::Real
    α::Real
    g::Real
    k::Real
    γ::Real
    dist::Union{UnivariateDistribution, Nothing} # underlying distribution
    F::Union{Function, Nothing}                  # CDF
    function GSDist(F₀, x₀, α, g, k, γ, dist, F)
        F₀ < 0 || F₀ > 1 && throw(DomainError(F₀, "F₀ must be between 0 and 1"))
        α ≤ 0 && throw(DomainError(α, "α must be a positive real number"))
        g ≤ 0 && throw(DomainError(g, "g must be a positive real number"))
        k ≤ 0 && throw(DomainError(k, "k must be a positive real number"))
        γ ≤ 0 && throw(DomainError(γ, "γ must be a positive real number"))

        new(F₀, x₀, α, g, k, γ, dist, F)
    end
end
GSDist(F₀, x₀, α, g, k, γ) = GSDist(F₀, x₀, α, g, k, γ, nothing, nothing)
GSDist(F₀, x₀, α, g, k, γ, dist) = GSDist(F₀, x₀, α, g, k, γ, dist, _piecewise_cdf(α, g, k, γ, dist))

function GSDist(dist::UnivariateDistribution, F₀::Real=0.5;
    n::Int=21, diff::Function=three_point_midpoint, h::Real=1.0)
    F₀ < 0 || F₀ > 1 && throw(DomainError(F₀, "F₀ must be between 0 and 1"))
    F₀ = Float64(F₀)

    l,u = extrema(dist)
    if isinf(l) l = quantile(dist, 1e-6) end
    if isinf(u) u = quantile(dist, 1 - 1e-6) end
    ql, qu = cdf.(dist, (l,u))

    q = range(ql, qu, length=n)
    X = typeof(dist) <: DiscreteUnivariateDistribution ? Set(quantile.(dist, q)) : quantile.(dist, q)

    FX = cdf.(dist, X)
    dF = y->diff(x->cdf(dist,x), y, h)
    fX = typeof(dist) <: DiscreteUnivariateDistribution ? dF.(X) : pdf.(dist, X)

    f = (t, p) -> p[1] * t.^p[2] .* (1.0 .- t.^p[3]).^p[4]
    
    p0 = [1, 0.5, 1, 0.5]
    fit = curve_fit(f, FX, fX, p0)
    p = coef(fit)

    F = _piecewise_cdf(p..., dist)

    GSDist(F₀, quantile(dist, F₀), p..., dist, F)
end

params(G::GSDist) = (G.α, G.g, G.k, G.γ)

function quantile(D::GSDist, p::Real)
    q = _gsd_quantile(p, D.F₀, D.x₀, D.α, D.g, D.k, D.γ)
    if isnan(q)
        @warn "Unable to calculate the quantile of the GSDist. Falling back to the quantile of the underlying distribution."
        quantile(D.dist, p)
    else
        q
    end
end

function mean(G::GSDist)
    m = _gsd_mean(G.F₀, G.x₀, G.α, G.g, G.k, G.γ)
    if isnan(m)
        @warn "Unable to calculate the mean of the GSDist. Falling back to the mean of the underlying distribution."
        mean(G.dist)
    else
        m
    end
end

function var(G::GSDist)
    m1 = mean(G)
    m2 = try
        _gsd_moment(G.F₀, G.x₀, G.α, G.g, G.k, G.γ; moment=2)
    catch y
        if isa(y, DomainError)
            @warn "Unable to calculate the variance of the GSDist. Falling back to the variance of the underlying distribution."
            return var(G.dist)
        end
    end
    m2 - m1^2
end

std(G::GSDist) = sqrt(var(G))

cdf(G::GSDist, x::Real) = G.F(x)

function pdf(G::GSDist, x::Real)
    F = cdf(G, x)
    α, g, k, γ = params(G)
    α * F^g * (1 - F^k)^γ
end

function logpdf(G::GSDist, x::Real)
    F = cdf(G, x)
    α, g, k, γ = params(G)
    log(α) + g*log(F) + γ * log(1 - F^k)
end

minimum(G::GSDist) = minimum(G.dist)
maximum(G::GSDist) = maximum(G.dist)

insupport(G::GSDist, x::Real) = minimum(G) ≤ x ≤ maximum(G)
insupport(G::GSDist, x::VecOrMat{<:Real}) = insupport.(G, x)

function mode(G::GSDist)
    α, g, k, γ = params(G)
    x = (g / (k*γ + g))^(1/k)
    quantile(G, x)
end

# skewness(G::GSDist) = _gsd_moment(G.F₀, G.x₀, G.α, G.g, G.k, G.γ; moment=3)
# kurtosis(G::GSDist) = _gsd_moment(G.F₀, G.x₀, G.α, G.g, G.k, G.γ; moment=4)
# entropy(G::GSDist, b::Real)
# mgf(G::GSDist, t::Any)
# cf(G::GSDist, t::Any)
