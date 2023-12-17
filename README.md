
<!-- README.md is generated from README.Rmd. Please edit that file -->

# extBatchMarking: Implementation of Extended Batch Marking Model

<!-- badges: start -->
<!-- badges: end -->

The primary objective of `extBatchMarking` is to facilitate the fitting
of models developed by Cowen et al., $2017$ for the ecologist. The
marked models can be seamlessly integrated with unmarked models to
estimate population size. The combined model harnesses the power of both
the N-mixture model and the Viterbi algorithm for hidden markov model to
provide accurate population size estimates.

The primary objective of `extBatchMarking` is to facilitate the fitting
of models developed by Cowen et al., $2017$ for the ecologist. The
marked models can be seamlessly integrated with unmarked models to
estimate population size. The combined model harnesses the power of both
the N-mixture model and the Viterbi algorithm for hidden markov model to
provide accurate population size estimates.

In ecological research, it’s often challenging to directly count every
individual of a species due to various factors such as their elusive
nature or inaccessible habitats. As a result, Cowen et al., $2017$
employ two distinct modeling approaches: marked and unmarked models.

1.  **Marked Models:** These models focus on individuals that have been
    uniquely identified or ‘marked’ in some way, such as through
    tagging, banding, or other identification methods. These marked
    individuals are tracked over time, and their data is used to
    estimate parameters related to the population, such as survival
    rates, population growth, or movement patterns.

2.  **Unmarked Models:** Unmarked models, on the other hand, are
    designed to estimate population parameters without relying on
    individual identification.

The beauty of combining these two modeling approaches lies in their
synergy. By leveraging both marked and unmarked data, ecologists can
achieve more accurate and robust estimates of population abundance.
Marked data provide insights into specific individuals, while unmarked
data give a broader perspective on the entire population.

In practice, this combination is achieved through sophisticated
statistical techniques, often utilizing concepts like the N-mixture
model and algorithms like the Viterbi algorithm. These methods allow
ecologists to integrate data from marked and unmarked individuals,
resulting in more comprehensive and reliable population abundance
estimates.

The models showcased in this example represent the foundational
instances drawn from a set of four complex models available within the
`extBatchMarking` package. These models serve as essential building
blocks for understanding the advanced functionalities and capabilities
offered by the package.

The `extBatchMarking` package is designed to empower researchers with a
powerful tool set for the analysis of batch-marked data in ecological
and population studies. It allows users to efficiently fit and assess
batch-marked models, aiding in the estimation of critical population
parameters such as survival and capture probabilities.

In particular, the models illustrated here provide a comprehensive
introduction to the core concepts and methodologies underpinning the
package’s functionality. They are intended to facilitate an initial
grasp of how to work with batch-marked data, offering insights into the
modeling techniques used in the field of population ecology.

It’s worth noting that the results obtained using the `extBatchMarking`
package align with the findings presented in Cowen et al. (2017). This
alignment demonstrates the package’s reliability and ability to
replicate established research outcomes. The Cowen et al. (2017) results
section serves as a benchmark against which the package’s performance
can be validated, providing users with confidence in the accuracy of
their analyses.

By starting with these basic examples, users can progressively delve
into more intricate and tailored analyses within the `extBatchMarking`
package, ultimately enabling them to make meaningful contributions to
the understanding of population dynamics and ecology. The package’s
versatility and fidelity to established research findings make it a
valuable resource for both novice and experienced researchers in the
field.

The example will guide you through the steps of how to employ this
approach effectively, demonstrating its relevance and importance in
ecological and wildlife studies. It showcases the power of merging
marked and unmarked models to gain a deeper understanding of species
populations and their dynamics within natural ecosystems.”

## Installation

You can install the released version of `extBatchMarking` from
[CRAN](https://CRAN.R-project.org) with:

``` r
devtools::load_all(".")
#> ℹ Loading extBatchMarking
devtools::document()
#> ℹ Updating extBatchMarking documentation
#> ℹ Loading extBatchMarking
#> Writing 'batchMarkHmmLL.Rd'
#> Writing 'batchUnmarkHmmLL.Rd'
#> Writing 'batchUnmarkViterbi.Rd'
#> Writing 'batchMarkUnmarkHmmLL.Rd'
#> Writing 'batchUnmark2Viterbi.Rd'
# usethis::use_package_doc()
devtools::load_all()
#> ℹ Loading extBatchMarking
```

## Example 1

This is a basic example which shows how to fit a Batch marking model
with constant `phi` and `p`. Example 1 can also be found in the Cowen et
al., 2017 results using the `WeatherLoach` data using in the same paper:

``` r

library(extBatchMarking)
```

Load the data `WeatherLoach` from the `extBatchMarking` package: Here is
the step-by-step guide on how to load data directly from
`extBatchMarking` package. The defult data discussed in
`Cowen et al. 2017`.

``` r

data("WeatherLoach", package = "extBatchMarking")
```

First, we show with an example how to fit the `batchMarkHmmLL` and
`batchMarkUnmarkHmmLL` functions. `batchMarkHmmLL` and
`batchMarkUnmarkHmmLL` functions output the unoptimized log-likelihood
values of marked only model and the combined models. These allow users
know if the likelihood functions can be computed at the specified
initial values. Otherwise, `NAN` or `Inf` will be returned. If so, the
arguments of the functions should be revisited.

``` r

# Initial parameter
theta <- c(0, -1)

cores <- detectCores()-1

res1 <- batchMarkHmmLL(par         = theta,
                       data        = WeatherLoach,
                       choiceModel = "model4",
                       cores)

res1
#> [1] 132.3349
```

## Example 2

``` r

thet <- c(0.1, 0.1, 7, -1.5)

res3 <- batchMarkUnmarkHmmLL(par         = thet,
                             data        = WeatherLoach,
                             choiceModel = "model4",
                             Umax        = 1800,
                             nBins       = 20,
                             cores       = cores)
```

## Example 3

    #> initial  value 132.334856 
    #> final  value 124.984186 
    #> converged

Model survival probability value:

``` r

res$phi
#> [1] 0.6
```

Model detection probability value:

``` r

res$p
#> [1] 0.2
```

Model log-likelihood value:

``` r

res$ll
#> [1] 124.98
```

Model AIC value

``` r

res$AIC
#> [1] 253.97
```

If `heesian = TRUE` the hessian matrix will be outputed

``` r

res$hessian
#>          [,1]     [,2]
#> [1,] 216.6291 141.3538
#> [2,] 141.3538 150.7772
```

## Example 4

This example serves as a fundamental illustration of the process of
combining both marked and unmarked models to estimate the population
abundance of a species. It demonstrates a key approach used in
ecological and wildlife studies to gain insights into the size of a
specific species population within a given habitat.

``` r

theta <- c(0.1, 0.1, 7, -1.5)

res2 <- batchMarkUnmarkOptim(par=theta,
                            data=WeatherLoach,
                            Umax=1800,
                            nBins=20,
                            choiceModel="model4",
                            popSize = "Horvitz_Thompson",
                            method="BFGS",
                            parallel=FALSE,
                            control=list(trace = 1))
#> initial  value 394.338558 
#> iter  10 value 202.900359
#> final  value 202.899020 
#> converged
```

results

``` r

res2$phi
#> [1] 0.59

res2$p
#> [1] 0.2

res2$lambda
#> [1] 1318.6

res2$gam
#> [1] 0.2

res2$ll
#> [1] 202.9

res2$AIC
#> [1] 413.8

res2$U
#>  [1] 1310 1010  790  630  490  370  290  250  230  210  170

res2$M
#>  [1]   0 160 250 300 110 100  50  80 100 130  60

res2$N
#>  [1] 1310 1170 1040  930  600  470  340  330  330  340  230
```
