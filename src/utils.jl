two_point_forward(f, x::Real, h::Real=1.0) = (f(x+h) - f(x)) / h
three_point_endpoint(f, x::Real, h::Real=1.0) = (-3f(x) + 4f(x+h)-f(x+2h)) / 2h
three_point_midpoint(f, x::Real, h::Real=1.0) = (f(x+h) - f(x-h)) / 2h
five_point_endpoint(f, x::Real, h::Real=1.0) = (-25f(x) + 48f(x+h) - 36f(x+2h) + 16f(x+3h) - 3f(x+4h)) / 12h
five_point_midpoint(f, x::Real, h::Real=1.0) = (f(x-2h) - 8f(x-h) + 8f(x+h)-f(x+2h)) / 12h


_beta_inc(z1::Real, z2::Real, a::Real, b::Real) = inv(a) * (z2^a * _₂F₁(a,1-b,a+1,z2) - z1^a * _₂F₁(a,1-b,a+1,z1))
_beta_inc(x::Real, a::Real, b::Real) = _beta_inc(0.0, x, a, b)
_beta_inc(a::Real, b::Real) = _beta_inc(0.0, prevfloat(1.0), a, b)