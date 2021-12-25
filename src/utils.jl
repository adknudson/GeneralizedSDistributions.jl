two_point_forward(f, x::Real, h::Real=1.0) = (f(x+h) - f(x)) / h
three_point_endpoint(f, x::Real, h::Real=1.0) = (-3f(x) + 4f(x+h)-f(x+2h)) / 2h
three_point_midpoint(f, x::Real, h::Real=1.0) = (f(x+h) - f(x-h)) / 2h
five_point_endpoint(f, x::Real, h::Real=1.0) = (-25f(x) + 48f(x+h) - 36f(x+2h) + 16f(x+3h) - 3f(x+4h)) / 12h
five_point_midpoint(f, x::Real, h::Real=1.0) = (f(x-2h) - 8f(x-h) + 8f(x+h)-f(x+2h)) / 12h


_beta_inc(z1::Real, z2::Real, a::Real, b::Real) = inv(a) * (z2^a * _₂F₁(a,1-b,a+1,z2) - z1^a * _₂F₁(a,1-b,a+1,z1))
_beta_inc(x::Real, a::Real, b::Real) = _beta_inc(0.0, x, a, b)
_beta_inc(a::Real, b::Real) = _beta_inc(0.0, prevfloat(1.0), a, b)


_gsd_quantile(p::Real, F₀::Real, x₀::Real, α::Real, g::Real, k::Real, γ::Real) = x₀ + _beta_inc(F₀^k, p^k, (1-g)/k, 1-γ) / (α*k)

function _gsd_mean(F₀::Real, x₀::Real, α::Real, g::Real, k::Real, γ::Real)
	z1 = F₀^k
	a1 = (1 - g) / k
	b1 = 1 - γ
	a2 = (2 - g) / k
	b2 = 1 - γ
	x₀ + (_beta_inc(z1, prevfloat(1.0), a1, b1) - _beta_inc(a2, b2)) * inv(α * k)
end

function _gsd_moment(F₀::Real, x₀::Real, α::Real, g::Real, k::Real, γ::Real; moment::Int=1)
    z1 = F₀^k
    a = (1 - g) / k
    b = 1 - γ
    c = inv(α * k)
    
    f = q -> (x₀ + c*_beta_inc(z1, q^k, a, b))^moment
    quadgk(f, 0, 1, atol=1e-4)[1]
end


function _piecewise_cdf(α, g, k, γ, dist; probs=[0.01, 0.15, 0.85, 0.99])
	p = [α, g, k, γ]
	cutpoint = Float64.(quantile.(dist, probs))
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