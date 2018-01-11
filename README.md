## mixedscores
code to implement the algorithm described in "Joint regression analysis of mixed-type outcome data via efficient scores".

### Usage + requirements
A reasonably current version of [R](https://cran.r-project.org/), along with the package RcppEigen and its dependencies, which would include [Rtools](https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows) for Windows users.

Alternatively, an R-only implementation is provided in `mbTest.R`. It is significantly slower than the C++ version, but this would likely cause an issue only when running extensive simulations.

The primary exposed (to R) function is `mbTest` (that is, the so-called *multiplier bootstrap* test). Three required arguments are:
1. `uu` The `n*q` matrix of scores, where each column represents the scores from a different margin.
2. `B`  The number of bootstrap iterations.
3. `dd` The norms to calculate.  May be an integer vector. Currently the supremum norm (requested via `-1`) and 1, 2, or 3 are meaningful arguments. Naturally, a positive number `p` induces calculations with the corresponding L-p norm.

`mbTest` returns a list of length 3 with the following elements:

* `obs.norm` The norm(s) of the observed data (one for each element of `dd`)
* `p.value`  The p-values corresponding to each element of `dd`
* `norms`    The norms generated via bootstrap sampling, since they may contain useful information beyond p-values. One column per element of `dd`.