# extBatchMarking: Implementation of Extended Batch Marking Model

The primary objective of `extBatchMarking` is to facilitate the fitting
of models developed by Cowen et al., 2017 for the ecologist. The marked
models can be seamlessly integrated with unmarked models to estimate
population size. The combined model harnesses the power of both the
N-mixture model and the Viterbi algorithm for hidden markov model to
provide accurate population size estimates.

The primary objective of `extBatchMarking` is to facilitate the fitting
of models developed by Cowen et al., 2017 for the ecologist. The marked
models can be seamlessly integrated with unmarked models to estimate
population size. The combined model harnesses the power of both the
N-mixture model and the Viterbi algorithm for hidden markov model to
provide accurate population size estimates.

In ecological research, it’s often challenging to directly count every
individual of a species due to various factors such as their elusive
nature or inaccessible habitats. As a result, Cowen et al., 2017 employ
two distinct modeling approaches: marked and unmarked models.

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

    devtools::load_all(".")
    #> ℹ Loading extBatchMarking
    devtools::document()
    #> ℹ Updating extBatchMarking documentation
    #> ℹ Loading extBatchMarking
    devtools::load_all()
    #> ℹ Loading extBatchMarking

## Example 1

This is a basic example which shows how to fit a Batch marking model
with constant `phi` and `p`. Example 1 can also be found in the Cowen et
al., 2017 results using the `WeatherLoach` data used in the same paper:

Load the data `WeatherLoach` from the `extBatchMarking` package: Here is
the step-by-step guide on how to load data directly from
`extBatchMarking` package. The default data discussed in
`Cowen et al. 2017`.


    data("WeatherLoach", package = "extBatchMarking")

First, we show with an example how to fit the `batchMarkHmmLL` and
`batchMarkUnmarkHmmLL` functions. `batchMarkHmmLL` and
`batchMarkUnmarkHmmLL` functions output the unoptimized log-likelihood
values of marked only model and the combined models. These allow users
know if the likelihood functions can be computed at the specified
initial values. Otherwise, `NAN` or `Inf` will be returned. If so, the
arguments of the functions should be revisited.


    # Initial parameter
    theta <- c(0, -1)

    res1 <- batchMarkHmmLL(par           = theta,
                           data          = WeatherLoach,
                           choiceModel   = "model4",
                           covariate_phi = NULL,
                           covariate_p   = NULL)

    res1
    #> [1] 132.3349

\[1\] 132.3349

## Example 2


    thet <- c(0.1, 0.1, 7, -1.5)

    res3 <- batchMarkUnmarkHmmLL(par           = thet,
                                 data          = WeatherLoach,
                                 choiceModel   = "model4",
                                 Umax          = 1800,
                                 nBins         = 600,
                                 covariate_phi = NULL,
                                 covariate_p   = NULL)

    res3
    #> [1] 870.0261

\[1\] 870.0261

## Example 3

    #> initial  value 132.334856 
    #> final  value 124.984186 
    #> converged

initial value 132.334856 final value 124.984186 converged

Print the results with the `print.batchMarkOptim()` function: This
prints to the console the log-likelihood, AIC, recapture probability and
survival probability with their corresponding standard errors.


    print(res)
    #> [[1]]
    #> 
    #> 
    #> | log.likelihood|      AIC|
    #> |--------------:|--------:|
    #> |       124.9842| 253.9684|
    #> 
    #> [[2]]
    #> 
    #> 
    #> |      p| p_S.Error|    phi| phi_S.Error|
    #> |------:|---------:|------:|-----------:|
    #> | 0.1964|    0.0206| 0.6042|      0.0261|

\[\[1\]\]

<table>
<thead>
<tr class="header">
<th style="text-align: right;">log.likelihood</th>
<th style="text-align: right;">AIC</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">124.9842</td>
<td style="text-align: right;">253.9684</td>
</tr>
</tbody>
</table>

\[\[2\]\]

<table>
<thead>
<tr class="header">
<th style="text-align: right;">p</th>
<th style="text-align: right;">p_S.Error</th>
<th style="text-align: right;">phi</th>
<th style="text-align: right;">phi_S.Error</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">0.1964</td>
<td style="text-align: right;">0.0206</td>
<td style="text-align: right;">0.6042</td>
<td style="text-align: right;">0.0261</td>
</tr>
</tbody>
</table>

Plot the results with the `plot.batchMarkOptim()` function: This plots
the goodness-of-fit plot to perform a goodness-of-fit test for the model
fit.


    plot(res)

![](README_files/figure-markdown_strict/unnamed-chunk-4-1.png)

## Example 4

This example serves as a fundamental illustration of the process of
combining both marked and unmarked models to estimate the population
abundance of a species. It demonstrates a key approach used in
ecological and wildlife studies to gain insights into the size of a
specific species population within a given habitat.


    theta <- c(0.1, 0.1, 7, -1.5)

    res2 <- batchMarkUnmarkOptim(par=theta,
                                data=WeatherLoach,
                                Umax=1800,
                                nBins=600,
                                choiceModel="model4",
                                popSize = "Horvitz_Thompson",
                                method="BFGS",
                                control=list(trace = 1),
                                covariate_phi = NULL,
                                covariate_p   = NULL)
    #> initial  value 870.026082 
    #> iter  10 value 205.041347
    #> final  value 205.031738 
    #> converged

initial value 870.026082 iter 10 value 205.041347 final value 205.031738
converged

Print the results with the `print.batchMarkUnmarkOptim()` function: This
prints to the console the log-likelihood, AIC, recapture probability and
survival probability with their corresponding standard errors. The
species abundances are also provided which include number of unmarked
individuals, number of marked individuals, and Total abundance.


    print(res2)
    #> [[1]]
    #> 
    #> 
    #> | log.likelihood|      AIC|
    #> |--------------:|--------:|
    #> |       205.0317| 418.0635|
    #> 
    #> [[2]]
    #> 
    #> 
    #> |      p| p_S.Error|   phi| phi_S.Error|
    #> |------:|---------:|-----:|-----------:|
    #> | 0.1866|    0.0048| 0.614|      0.0169|
    #> 
    #> [[3]]
    #> 
    #> 
    #> | No.of.Unmarked.U.| No.of.Marked.M.| Abundance.N.|
    #> |-----------------:|---------------:|------------:|
    #> |              1500|               0|         1500|
    #> |               900|             171|         1071|
    #> |               900|             268|         1168|
    #> |               900|             322|         1222|
    #> |               300|             118|          418|
    #> |               300|             107|          407|
    #> |               300|              54|          354|
    #> |               300|              86|          386|
    #> |               300|             107|          407|
    #> |               300|             139|          439|
    #> |               300|              64|          364|
    #> 
    #> [[4]]
    #> 
    #> 
    #> |   Lambda| Lambda_S.Error|    Gam| Gam_S.Error|
    #> |--------:|--------------:|------:|-----------:|
    #> | 1346.844|       1903.361| 0.0599|       0.022|

\[\[1\]\]

<table>
<thead>
<tr class="header">
<th style="text-align: right;">log.likelihood</th>
<th style="text-align: right;">AIC</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">205.0317</td>
<td style="text-align: right;">418.0635</td>
</tr>
</tbody>
</table>

\[\[2\]\]

<table>
<thead>
<tr class="header">
<th style="text-align: right;">p</th>
<th style="text-align: right;">p_S.Error</th>
<th style="text-align: right;">phi</th>
<th style="text-align: right;">phi_S.Error</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">0.1866</td>
<td style="text-align: right;">0.0048</td>
<td style="text-align: right;">0.614</td>
<td style="text-align: right;">0.0169</td>
</tr>
</tbody>
</table>

\[\[3\]\]

<table>
<thead>
<tr class="header">
<th style="text-align: right;">No.of.Unmarked.U.</th>
<th style="text-align: right;">No.of.Marked.M.</th>
<th style="text-align: right;">Abundance.N.</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">1500</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1500</td>
</tr>
<tr class="even">
<td style="text-align: right;">900</td>
<td style="text-align: right;">171</td>
<td style="text-align: right;">1071</td>
</tr>
<tr class="odd">
<td style="text-align: right;">900</td>
<td style="text-align: right;">268</td>
<td style="text-align: right;">1168</td>
</tr>
<tr class="even">
<td style="text-align: right;">900</td>
<td style="text-align: right;">322</td>
<td style="text-align: right;">1222</td>
</tr>
<tr class="odd">
<td style="text-align: right;">300</td>
<td style="text-align: right;">118</td>
<td style="text-align: right;">418</td>
</tr>
<tr class="even">
<td style="text-align: right;">300</td>
<td style="text-align: right;">107</td>
<td style="text-align: right;">407</td>
</tr>
<tr class="odd">
<td style="text-align: right;">300</td>
<td style="text-align: right;">54</td>
<td style="text-align: right;">354</td>
</tr>
<tr class="even">
<td style="text-align: right;">300</td>
<td style="text-align: right;">86</td>
<td style="text-align: right;">386</td>
</tr>
<tr class="odd">
<td style="text-align: right;">300</td>
<td style="text-align: right;">107</td>
<td style="text-align: right;">407</td>
</tr>
<tr class="even">
<td style="text-align: right;">300</td>
<td style="text-align: right;">139</td>
<td style="text-align: right;">439</td>
</tr>
<tr class="odd">
<td style="text-align: right;">300</td>
<td style="text-align: right;">64</td>
<td style="text-align: right;">364</td>
</tr>
</tbody>
</table>

\[\[4\]\]

<table>
<thead>
<tr class="header">
<th style="text-align: right;">Lambda</th>
<th style="text-align: right;">Lambda_S.Error</th>
<th style="text-align: right;">Gam</th>
<th style="text-align: right;">Gam_S.Error</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">1346.844</td>
<td style="text-align: right;">1903.361</td>
<td style="text-align: right;">0.0599</td>
<td style="text-align: right;">0.022</td>
</tr>
</tbody>
</table>

Plot the results with the `plot.batchMarkOptim()` function: This plots
the googdness-of-fit plot to perform a goodness-of-fit test for the
model fit.


    plot(res2)

![](README_files/figure-markdown_strict/unnamed-chunk-7-1.png)

### Model 2 for marked model with covariate for phi

This provides some level of flexibility to the marked model by adding a
placeholder for time-dependent covariate to better estimate parameter
phi. For example, we might want to understand how temperature affects
the survival rate of a species.


    #-------------------------------------------------
    # Model 2: 10 phis and 1 prob
    #-------------------------------------------------

    theta <- c(-1, rep(0, 10))

    cv <- matrix(seq(2, 10, length = 10), ncol = 1)
    batchMarkHmmLL(theta, WeatherLoach, "model2", covariate_phi = cv, covariate_p = NULL)
    #> [1] 132.3349
    res <- batchMarkOptim(par=theta, data=WeatherLoach, covariate_phi = cv, covariate_p = NULL, 
                          choiceModel = "model2", method="BFGS", control = list(trace = 1))
    #> initial  value 132.334856 
    #> iter  10 value 97.474508
    #> iter  20 value 95.823647
    #> iter  30 value 95.787362
    #> iter  40 value 95.784462
    #> iter  50 value 95.783391
    #> iter  60 value 95.783072
    #> final  value 95.783042 
    #> converged

initial value 132.334856 iter 10 value 97.474508 iter 20 value 95.823647
iter 30 value 95.787362 iter 40 value 95.784462 iter 50 value 95.783391
iter 60 value 95.783072 final value 95.783042 converged


    #-------------------------------------------------
    # Model 2: 10 phis and 1 prob
    #-------------------------------------------------

    print(res)
    #> [[1]]
    #> 
    #> 
    #> | log.likelihood|      AIC|
    #> |--------------:|--------:|
    #> |         95.783| 213.5661|
    #> 
    #> [[2]]
    #> 
    #> 
    #> |     p| p_S.Error|    phi| phi_S.Error|
    #> |-----:|---------:|------:|-----------:|
    #> | 0.181|    0.0316| 0.3268|      0.0000|
    #> | 0.181|    0.0316| 0.0000|      0.0000|
    #> | 0.181|    0.0316| 1.0000|      0.8321|
    #> | 0.181|    0.0316| 0.1082|      0.1306|
    #> | 0.181|    0.0316| 0.5680|      0.1388|
    #> | 0.181|    0.0316| 0.4722|      0.0009|
    #> | 0.181|    0.0316| 0.9973|      8.2224|
    #> | 0.181|    0.0316| 0.6458|      0.0985|
    #> | 0.181|    0.0316| 0.6619|      0.1316|
    #> | 0.181|    0.0316| 0.4397|      0.0285|

\[\[1\]\]

<table>
<thead>
<tr class="header">
<th style="text-align: right;">log.likelihood</th>
<th style="text-align: right;">AIC</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">95.783</td>
<td style="text-align: right;">213.5661</td>
</tr>
</tbody>
</table>

\[\[2\]\]

<table>
<thead>
<tr class="header">
<th style="text-align: right;">p</th>
<th style="text-align: right;">p_S.Error</th>
<th style="text-align: right;">phi</th>
<th style="text-align: right;">phi_S.Error</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">0.181</td>
<td style="text-align: right;">0.0316</td>
<td style="text-align: right;">0.3268</td>
<td style="text-align: right;">0.0000</td>
</tr>
<tr class="even">
<td style="text-align: right;">0.181</td>
<td style="text-align: right;">0.0316</td>
<td style="text-align: right;">0.0000</td>
<td style="text-align: right;">0.0000</td>
</tr>
<tr class="odd">
<td style="text-align: right;">0.181</td>
<td style="text-align: right;">0.0316</td>
<td style="text-align: right;">1.0000</td>
<td style="text-align: right;">0.8321</td>
</tr>
<tr class="even">
<td style="text-align: right;">0.181</td>
<td style="text-align: right;">0.0316</td>
<td style="text-align: right;">0.1082</td>
<td style="text-align: right;">0.1306</td>
</tr>
<tr class="odd">
<td style="text-align: right;">0.181</td>
<td style="text-align: right;">0.0316</td>
<td style="text-align: right;">0.5680</td>
<td style="text-align: right;">0.1388</td>
</tr>
<tr class="even">
<td style="text-align: right;">0.181</td>
<td style="text-align: right;">0.0316</td>
<td style="text-align: right;">0.4722</td>
<td style="text-align: right;">0.0009</td>
</tr>
<tr class="odd">
<td style="text-align: right;">0.181</td>
<td style="text-align: right;">0.0316</td>
<td style="text-align: right;">0.9973</td>
<td style="text-align: right;">8.2224</td>
</tr>
<tr class="even">
<td style="text-align: right;">0.181</td>
<td style="text-align: right;">0.0316</td>
<td style="text-align: right;">0.6458</td>
<td style="text-align: right;">0.0985</td>
</tr>
<tr class="odd">
<td style="text-align: right;">0.181</td>
<td style="text-align: right;">0.0316</td>
<td style="text-align: right;">0.6619</td>
<td style="text-align: right;">0.1316</td>
</tr>
<tr class="even">
<td style="text-align: right;">0.181</td>
<td style="text-align: right;">0.0316</td>
<td style="text-align: right;">0.4397</td>
<td style="text-align: right;">0.0285</td>
</tr>
</tbody>
</table>


    plot(res)

![](README_files/figure-markdown_strict/unnamed-chunk-10-1.png)

### Model 3 for marked model with covariate for p

This provides some level of flexibility to the marked model by adding a
placeholder for time-dependent covariate to better estimate parameter
p. For example, we might want to understand how temperature affects the
survival rate of a species.


    #-------------------------------------------------
    # Model 3: 1 phis and 10 prob
    #-------------------------------------------------

    theta <- c(-1, rep(0, 10))

    cv <- matrix(seq(2, 10, length = 10), ncol = 1)
    batchMarkHmmLL(theta, WeatherLoach, "model3", covariate_phi = NULL, covariate_p = cv)
    #> [1] 200.0016
    res <- batchMarkOptim(par=theta, data=WeatherLoach, covariate_phi = NULL, covariate_p = cv, 
                          choiceModel = "model3", method="BFGS", control = list(trace = 1))
    #> initial  value 200.001592 
    #> iter  10 value 114.293180
    #> final  value 97.288967 
    #> converged

initial value 200.001592 iter 10 value 114.293180 final value 97.288967
converged


    #-------------------------------------------------
    # Model 2: 1 phis and 10 prob
    #-------------------------------------------------

    print(res)
    #> [[1]]
    #> 
    #> 
    #> | log.likelihood|      AIC|
    #> |--------------:|--------:|
    #> |         97.289| 216.5779|
    #> 
    #> [[2]]
    #> 
    #> 
    #> |      p| p_S.Error|   phi| phi_S.Error|
    #> |------:|---------:|-----:|-----------:|
    #> | 0.9559|    0.0182| 0.627|      0.0297|
    #> | 0.9991|    0.0011| 0.627|      0.0297|
    #> | 0.0084|    0.0099| 0.627|      0.0297|
    #> | 0.0156|    0.0078| 0.627|      0.0297|
    #> | 0.0793|    0.0241| 0.627|      0.0297|
    #> | 0.0933|    0.0260| 0.627|      0.0297|
    #> | 0.2411|    0.0410| 0.627|      0.0297|
    #> | 0.3234|    0.0388| 0.627|      0.0297|
    #> | 0.3903|    0.0352| 0.627|      0.0297|
    #> | 0.3164|    0.0311| 0.627|      0.0297|
    plot(res)

![](README_files/figure-markdown_strict/unnamed-chunk-12-1.png)

\[\[1\]\]

<table>
<thead>
<tr class="header">
<th style="text-align: right;">log.likelihood</th>
<th style="text-align: right;">AIC</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">97.289</td>
<td style="text-align: right;">216.5779</td>
</tr>
</tbody>
</table>

\[\[2\]\]

<table>
<thead>
<tr class="header">
<th style="text-align: right;">p</th>
<th style="text-align: right;">p_S.Error</th>
<th style="text-align: right;">phi</th>
<th style="text-align: right;">phi_S.Error</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">0.9559</td>
<td style="text-align: right;">0.0182</td>
<td style="text-align: right;">0.627</td>
<td style="text-align: right;">0.0297</td>
</tr>
<tr class="even">
<td style="text-align: right;">0.9991</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.627</td>
<td style="text-align: right;">0.0297</td>
</tr>
<tr class="odd">
<td style="text-align: right;">0.0084</td>
<td style="text-align: right;">0.0099</td>
<td style="text-align: right;">0.627</td>
<td style="text-align: right;">0.0297</td>
</tr>
<tr class="even">
<td style="text-align: right;">0.0156</td>
<td style="text-align: right;">0.0078</td>
<td style="text-align: right;">0.627</td>
<td style="text-align: right;">0.0297</td>
</tr>
<tr class="odd">
<td style="text-align: right;">0.0793</td>
<td style="text-align: right;">0.0241</td>
<td style="text-align: right;">0.627</td>
<td style="text-align: right;">0.0297</td>
</tr>
<tr class="even">
<td style="text-align: right;">0.0933</td>
<td style="text-align: right;">0.0260</td>
<td style="text-align: right;">0.627</td>
<td style="text-align: right;">0.0297</td>
</tr>
<tr class="odd">
<td style="text-align: right;">0.2411</td>
<td style="text-align: right;">0.0410</td>
<td style="text-align: right;">0.627</td>
<td style="text-align: right;">0.0297</td>
</tr>
<tr class="even">
<td style="text-align: right;">0.3234</td>
<td style="text-align: right;">0.0388</td>
<td style="text-align: right;">0.627</td>
<td style="text-align: right;">0.0297</td>
</tr>
<tr class="odd">
<td style="text-align: right;">0.3903</td>
<td style="text-align: right;">0.0352</td>
<td style="text-align: right;">0.627</td>
<td style="text-align: right;">0.0297</td>
</tr>
<tr class="even">
<td style="text-align: right;">0.3164</td>
<td style="text-align: right;">0.0311</td>
<td style="text-align: right;">0.627</td>
<td style="text-align: right;">0.0297</td>
</tr>
</tbody>
</table>
