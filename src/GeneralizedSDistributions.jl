module GeneralizedSDistributions

using Distributions: UnivariateDistribution, ContinuousUnivariateDistribution, DiscreteUnivariateDistribution
using LsqFit: curve_fit, coef
using QuadGK: quadgk
using SpecialFunctions: beta, beta_inc
using OrdinaryDiffEq: ODEProblem, solve, AutoTsit5, Rosenbrock23

import Distributions: 
rand, sampler, logpdf, cdf, quantile, minimum, maximum, insupport, mean, var, std, modes, mode, skewness, kurtosis, 
entropy, mgf, cf, pdf, params


export GSDist


include("GSDist.jl")


end
