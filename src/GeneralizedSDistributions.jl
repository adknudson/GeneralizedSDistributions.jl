module GeneralizedSDistributions

import Distributions: UnivariateDistribution, ContinuousUnivariateDistribution, DiscreteUnivariateDistribution
import Distributions: rand, sampler, logpdf, cdf, quantile, minimum, maximum, insupport
import Distributions: mean, var, modes, mode, skewness, kurtosis, entropy, mgf, cf
import Distributions: pdf

import LsqFit: curve_fit, coef
import HypergeometricFunctions: _₂F₁
import QuadGK: quadgk


export
    GSDist,
    rand, sampler, quantile, mean, var, std


include("GSDist.jl")
include("utils.jl")

end
