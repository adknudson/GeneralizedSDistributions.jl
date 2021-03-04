module GeneralizedSDistributions

import Distributions: UnivariateDistribution, 
                      ContinuousUnivariateDistribution, 
                      DiscreteUnivariateDistribution
import Distributions: rand, sampler, logpdf, cdf, quantile, minimum, maximum, insupport,
                      mean, var, modes, mode, skewness, kurtosis, entropy, mgf, cf,
                      pdf, params

import LsqFit: curve_fit, coef
import HypergeometricFunctions: _₂F₁
import QuadGK: quadgk
import DifferentialEquations: ODEProblem, solve


export
    GSDist,
    rand, 
    sampler, 
    logpdf, pdf,
    cdf,
    quantile, 
    minimum, 
    maximum, 
    insupport,
    mean, 
    var, 
    std,
    mode


include("GSDist.jl")
include("show.jl")
include("utils.jl")

end
