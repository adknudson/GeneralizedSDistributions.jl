module GeneralizedSDistributions

import Distributions: UnivariateDistribution, 
                      ContinuousUnivariateDistribution, 
                      DiscreteUnivariateDistribution
import Distributions: rand, sampler, logpdf, cdf, quantile, minimum, maximum, insupport,
                      mean, var, std, modes, mode, skewness, kurtosis, entropy, mgf, cf,
                      pdf, params

import LsqFit: curve_fit, coef
import HypergeometricFunctions: _₂F₁
import QuadGK: quadgk
import OrdinaryDiffEq: ODEProblem, solve, AutoTsit5, Rosenbrock23


export GSDist


include("GSDist.jl")
include("show.jl")
include("utils.jl")

end
