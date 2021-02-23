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
    F::Union{UnivariateDistribution, Nothing}
    function GSDist(F₀, x₀, α, g, k, γ, F)
        F₀ < 0 || F₀ > 1 && throw(DomainError(F₀, "F₀ must be between 0 and 1"))
        α ≤ 0 && throw(DomainError(α, "α must be a positive real number"))
        g ≤ 0 && throw(DomainError(g, "g must be a positive real number"))
        k ≤ 0 && throw(DomainError(k, "k must be a positive real number"))
        γ ≤ 0 && throw(DomainError(γ, "γ must be a positive real number"))

        new(F₀, x₀, α, g, k, γ, F)
    end
end

# given the parameters directly
GSDist(F₀::Real, x₀::Real, α::Real, g::Real, k::Real, γ::Real) = GSDist(F₀, x₀, α, g, k, γ, nothing)

# given a univariate distribution
function GSDist(F::UnivariateDistribution, F₀::Real=0.5; n::Int=21, diff::Function=three_point_midpoint, h::Real=1.0)
    F₀ < 0 || F₀ > 1 && throw(DomainError(F₀, "F₀ must be between 0 and 1"))
    F₀ = Float64(F₀)

    l,u = extrema(F)
    if isinf(l) l = quantile(F, 1e-6) end
    if isinf(u) u = quantile(F, 1 - 1e-6) end
    ql, qu = cdf.(F, (l,u))

    q = range(ql, qu, length=n)
    X = typeof(F) <: DiscreteUnivariateDistribution ? Set(quantile.(F, q)) : quantile.(F, q)

    FX = cdf.(F, X)
    dF = y->diff(x->cdf(F,x), y, h)
    fX = typeof(F) <: DiscreteUnivariateDistribution ? dF.(X) : pdf.(F, X)

    f = (t, p) -> p[1] * t.^p[2] .* (1.0 .- t.^p[3]).^p[4]
    
    p0 = [1, 0.5, 1, 0.5]
    fit = curve_fit(f, FX, fX, p0)
    p = coef(fit)

    GSDist(F₀, quantile(F, F₀), p..., F)
end

# given data
# function GSDist(x::Vector{<:Real}) end


# Quantile function
_gsd_quantile(p::Real, F₀::Real, x₀::Real, α::Real, g::Real, k::Real, γ::Real) = x₀ + _beta_inc(F₀^k, p^k, (1-g)/k, 1-γ) / (α*k)
function quantile(D::GSDist, p::Real)
    q = _gsd_quantile(p, D.F₀, D.x₀, D.α, D.g, D.k, D.γ)
    if isnan(q)
        @warn "Unable to calculate the quantile of the GSDist. Falling back to the quantile of the underlying distribution."
        quantile(D.F, p)
    else
        q
    end
end


function _gsd_mean(F₀::Real, x₀::Real, α::Real, g::Real, k::Real, γ::Real)
	z1 = F₀^k
	a1 = (1 - g) / k
	b1 = 1 - γ
	a2 = (2 - g) / k
	b2 = 1 - γ
	x₀ + (_beta_inc(z1, prevfloat(1.0), a1, b1) - _beta_inc(a2, b2)) * inv(α * k)
end
function mean(D::GSDist)
    m = _gsd_mean(D.F₀, D.x₀, D.α, D.g, D.k, D.γ)
    if isnan(m)
        @warn "Unable to calculate the mean of the GSDist. Falling back to the mean of the underlying distribution."
        mean(D.F)
    else
        m
    end
end

function _gsd_moment(F₀::Real, x₀::Real, α::Real, g::Real, k::Real, γ::Real; moment::Int=1)
    z1 = F₀^k
    a = (1 - g) / k
    b = 1 - γ
    c = inv(α * k)
    
    f = q -> (x₀ + c*_beta_inc(z1, q^k, a, b))^moment
    quadgk(f, 0, 1, atol=1e-4)[1]
end



function var(D::GSDist)
    m1 = mean(D)
    m2 = try
        _gsd_moment(D.F₀, D.x₀, D.α, D.g, D.k, D.γ; moment=2)
    catch y
        if isa(y, DomainError)
            @warn "Unable to calculate the variance of the GSDist. Falling back to the variance of the underlying distribution."
            return var(D.F)
        end
    end
    m2 - m1^2
end
std(D::GSDist) = sqrt(var(D))
