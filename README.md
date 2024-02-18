# GeneralizedSDistributions

Generalized S-Distributions: a family of distributions that can serve as statistical models for unimodal distributions.

The inspiration for this package stems from the problem of matching Pearson correlation coefficients for use in Gaussian copulas, a problem discussed by Xiao and Zhou (2019). Matching coefficients for pairs of continuous distributions is made tractable by using a weighted sum of Hermite polynomials, but for discrete distributions (especially discrete distributions with large support sets or infinite support), the problem is not so simple or efficient. 

The solution is to approximate discrete distributions by a continous distribution. However the question is then what distribution to use. Generalized S-Distributions solve this problem by approximating a wide range of unimodal distributions, both continuous and discrete.

## Examples

```julia-repl
julia> using GeneralizedSDistributions, Distributions

julia> D = NegativeBinomial(20, 0.002)
NegativeBinomial{Float64}(r=20.0, p=0.002)

julia> G = GSDist(D)
Generalized S-Distribution (NegativeBinomial{Float64}(r=20.0, p=0.002))
 α: 0.010708744772447428
 g: 1.0034154213417692
 k: 0.026052263652808192
 γ: 0.8423059246978783

julia> mean(D), mean(G)
(9980.0, 9973.383253892538)

julia> var(D), var(G)
(4.99e6, 4.914780667295575e6)

julia> quantile(D, 0.8), quantile(G, 0.8)
(11795, 11801.44831411299)

julia> cdf(D, 10000), cdf(G, 10000)
(0.5333848577165967, 0.5341049437932197)
```

## Functionality and Limitations

Functionality implemented:

* `rand` / `sampler`
* `cdf` / `pdf`
* `quantile`
* `mean` / `median`
* `modes` / `mode`
* `var` / `std`
* `minimum` / `maximum`
* `insupport`

Not yet implemented:

* `skewness` / `kurtosis` / `entropy`

*Almost* implemented:

* `moment` (the $j^{th}$ moment of the distribution)

When the mean/variance/quantile of a GSDist cannot be computed, then the fallback is to return the measure from the underlying distribution. The measures often cannot be computed if the median of the GSDist is very far from zero (i.e. a median above `1e6`). Implementing using BigFloats may solve this in the future.

Additionally, there are plans to allow for the estimation of a GSDist using sample data (the primary motivation for studying Generalized S-Distributions).

## References

* Voit, E. O. (1992). The S‐Distribution A Tool for Approximation and Classification of Univariate, Unimodal Probability Distributions. Biometrical Journal, 34(7), 855-878.
* Voit, E. O., & Yu, S. (1994). The S‐distribution: approximation of discrete distributions. Biometrical Journal, 36(2), 205-219.
* Hernández–Bermejo, B., & Sorribas, A. (2001). Analytical Quantile Solution for the S‐distribution, Random Number Generation and Statistical Data Modeling. Biometrical Journal, 43(8), 1007-1025.
* Muino, J. M., Voit, E. O., & Sorribas, A. (2006). GS-distributions: A new family of distributions for continuous unimodal variables. Computational statistics & data analysis, 50(10), 2769-2798.
* Xiao, Q., & Zhou, S. (2019). Matching a correlation coefficient by a Gaussian copula. Communications in Statistics-Theory and Methods, 48(7), 1728-1747.
