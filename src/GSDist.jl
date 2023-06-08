two_point_forward(f, x::Real, h::Real=1.0) = (f(x+h) - f(x)) / h
three_point_endpoint(f, x::Real, h::Real=1.0) = (-3f(x) + 4f(x+h)-f(x+2h)) / 2h
three_point_midpoint(f, x::Real, h::Real=1.0) = (f(x+h) - f(x-h)) / 2h
five_point_endpoint(f, x::Real, h::Real=1.0) = (-25f(x) + 48f(x+h) - 36f(x+2h) + 16f(x+3h) - 3f(x+4h)) / 12h
five_point_midpoint(f, x::Real, h::Real=1.0) = (f(x-2h) - 8f(x-h) + 8f(x+h)-f(x+2h)) / 12h

"""
Unregularized incomplete beta function.
"""
function _beta_inc end
_beta_inc(a::Real, b::Real, z1::Real, z2::Real) = _beta_inc(a, b, z2) - _beta_inc(a, b, z1)
_beta_inc(a::Real, b::Real, x::Real) = beta_inc(a, b, x)[1] * beta(a, b)
_beta_inc(a::Real, b::Real) = beta(a, b)


_gsd_quantile(p::Real, F₀::Real, x₀::Real, α::Real, g::Real, k::Real, γ::Real) = x₀ + _beta_inc((1-g)/k, 1-γ, F₀^k, p^k) / (α*k)

function _gsd_mean(F₀::Real, x₀::Real, α::Real, g::Real, k::Real, γ::Real)
	z1 = F₀^k
	a1 = (1 - g) / k
	b1 = 1 - γ
	a2 = (2 - g) / k
	b2 = 1 - γ
	x₀ + (_beta_inc(a1, b1, z1, 1) - _beta_inc(a2, b2)) / (α * k)
end

function _gsd_moment(F₀::Real, x₀::Real, α::Real, g::Real, k::Real, γ::Real; moment::Int=1)
    z1 = F₀^k
    a = (1 - g) / k
    b = 1 - γ
    c = 1 / (α * k)
    
    f = q -> (x₀ + c * _beta_inc(a, b, z1, q^k))^moment
    quadgk(f, 0, 1, atol=1e-4)[1]
end


function _piecewise_cdf(α, g, k, γ, dist; probs=[0.01, 0.15, 0.85, 0.99])
	p = [α, g, k, γ]
	cutpoint = Float64.(quantile(dist, probs))
	d = length(probs) - 1

    function gsd!(du, u, p, t)
        α, g, k, γ = p
        du[1] = α * u[1]^g * (1 - u[1]^k)^γ
    end

	P = [ODEProblem(gsd!, [probs[i]], (cutpoint[i], cutpoint[i+1]), p) for i in 1:d]
	S = [solve(x, AutoTsit5(Rosenbrock23())) for x in P]
	
	# Make a piecewise CDF
	x -> begin
		x < cutpoint[1]   && return max(S[1](x)[1],   0.0)
		x > cutpoint[end] && return min(S[end](x)[1], 1.0)
		idx = findfirst(cutpoint[1:end-1] .≤ x .≤ cutpoint[2:end])
		return S[idx](x)[1]
	end
end



"""
    GSDist <: ContinuousUnivariateDistribution

Generalized S-Distribution for approximating univariate distributions.

# Arguments
- `F₀`: the reference quantile.
- `x₀`: the value at the reference quantile.
- `α, g, k, γ`: parameters of the Generalized S-Distribution
- `dist`: the underlying probability distribution.
- `F`: the cdf of the distribution.
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

"""
    GSDist(dist::UnivariateDistribution, F₀::Real=0.5;
    n::Int=21, diff::Function=three_point_midpoint, h::Real=1.0)

# Arguments
- `dist::UnivariateDistribution`: a univariate distribution.
- `F₀::Real`: the reference quantile to center the distribution around.
- `n::Integer`: the number of points on the distribution to use for least squares fit.
- `diff::Function`: one of the finite difference methods for gradient approximation.
- `h::Real`: the denominator used in the diff function.
"""
function GSDist(dist::UnivariateDistribution, F₀::Real=0.5;
    n::Int=21, diff::Function=three_point_midpoint, h::Real=1.0)
    F₀ < 0 || F₀ > 1 && throw(DomainError(F₀, "F₀ must be between 0 and 1"))
    F₀ = Float64(F₀)

    l, u = extrema(dist)
    l = isinf(l) ? quantile(dist, 1e-6) : l
    u = isinf(u) ? quantile(dist, 1 - 1e-6) : u

    ql, qu = cdf.(dist, (l,u))

    q = range(ql, qu, length=n)
    X = typeof(dist) <: DiscreteUnivariateDistribution ? Set(quantile.(dist, q)) : quantile.(dist, q)

    FX = cdf.(dist, X)
    dF = y -> diff(x -> cdf(dist, x), y, h)
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



function Base.show(io::IO, ::MIME"text/plain", G::GSDist)
    α, g, k, γ = params(G)
    println(io, "Generalized S-Distribution ($(G.dist))")
    println(io, " α: $α")
    println(io, " g: $g")
    println(io, " k: $k")
    print(io, " γ: $γ")
end