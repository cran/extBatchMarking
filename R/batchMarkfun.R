#' @title Marked model only.

#' @description batchMarkOptim function provides the batch marking function to be optimized.
#'
#' @param par Initial values for the parameters to be optimized over.
#' @param data A capture-recapture data matrix or data frame
#' @param choiceModel This chooses among different models and allow for model selection
#' @param cores The number of cores for parallelization
#' @param method The method to be used. See \code{\link{optim}} for details.
#' @param lowerBound Lower bounds on the variables for the "L-BFGS-B" \code{method}.
#' @param parallel Logical. Should the algorithm be run in parallel? This will be implemented in a future version.
#' @param control A list of control parameters. See optim for details.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' @param ... Further arguments to be passed by user which goes into the \code{\link{optim}} function.

#' @details
#' Note that arguments after ... must be matched exactly.
#' \code{\link{batchMarkOptim}} depends on \code{\link{optim}} function to optimize the parameters of the marked model only. By default optim performs minimization.
#'
#' @references Laura L. E. Cowen, Panagiotis Besbeas, Byron J. T. Morgan, 2017.: Hidden Markov Models for Extended Batch Data,
#' Biometrics, 73, 1321-1331. DOI: 10.1111/biom.12701.
#'
#' @rdname batchMarkOptim
#' @importFrom stats optim dpois ppois
#' @importFrom parallel makeCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom utils install.packages installed.packages
#'
#' @return For \code{batchMarkOptim}, a list with components:
#' \describe{
#'  \item{phi}{The survival probability and remaining in the population between occasion t and t+1.}
#'  \item{p}{The capture probability at occasion time t.}
#'  \item{ll}{The optimized log-likelihood value of marked model.}
#'  \item{hessian}{The hessian matrix.}
#'  \item{AIC}{The Akaike Information Criteria for model selection.}
#' }
#'
#' @export
#'
#' @examples
#' # Load the package
#' library(extBatchMarking)
#'
#' # Load the WeatherLoach data from Cowen et al., 2017.
#' data(WeatherLoach)
#'
#' # Initial parameter values
#' theta <- c(0, -1)
#'
#' \donttest{
#' mod1 <- batchMarkOptim(
#'            par         = theta,
#'            data        = WeatherLoach,
#'            choiceModel = "model4",
#'            method      = "BFGS",
#'            parallel    = FALSE,
#'            hessian     = TRUE,
#'            control     = list(trace = 1)
#'      )
#'
#'  # Survival probability
#'  mod1$phi
#'  # Capture probability
#'  mod1$p
#'  # Optimized log-likelihood
#'  mod1$ll
#'  # The Hessian matrix
#'  mod1$hessian
#'  # The Aikaike Information Criteria
#'  mod1$AIC
#'  }
#'
#'  \donttest{
#'  mod2 <- batchMarkOptim(
#'            par         = theta,
#'            data        = WeatherLoach,
#'            choiceModel = "model4",
#'            method      = "L-BFGS-B",
#'            parallel    = FALSE,
#'            hessian     = TRUE,
#'            control     = list(trace = 1))
#'
#'  # Survival probability
#'  mod2$phi
#'  # Capture probability
#'  mod2$p
#'  # Optimized log-likelihood
#'  mod2$ll
#'  # The Hessian matrix
#'  mod2$hessian
#'  # The Akaike Information Criteria
#'  mod2$AIC
#'  }


batchMarkOptim <- function(par=NULL, data, choiceModel=c("model1", "model2", "model3", "model4"),
                           method=c("Nelder-Mead","BFGS", "CG", "L-BFGS-B"),parallel=FALSE, lowerBound=-Inf,
                           cores=1, hessian=FALSE, control, ...){

  # my_packages <- c("Rcpp", "parallel", "optimParallel", "RcppArmadillo")                  # Specify your packages
  # not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]    # Extract not installed packages
  # if(length(not_installed)) install.packages(not_installed) # Stop if not installed

  if(parallel) {

    cores   <- parallel::detectCores()-1
    cluster <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cluster, cores = cores)
    parallel::clusterExport(cluster, "batchLogit",     envir=environment())
    parallel::clusterExport(cluster, "probs",          envir=environment())
    parallel::clusterExport(cluster, "delta_g",        envir=environment())
    parallel::clusterExport(cluster, "batchLL",        envir=environment())
    parallel::clusterExport(cluster, "batchMarkHmmLL", envir=environment())
    parallel::clusterExport(cluster, "gamma_gt",       envir=environment())
    parallel::clusterExport(cluster, "data",           envir=environment())
    parallel::clusterExport(cluster, "%dopar%",        envir=environment())
    parallel::clusterExport(cluster, "%do%",           envir=environment())
    parallel::clusterExport(cluster, "%:%",            envir=environment())
    parallel::clusterExport(cluster, "choiceModel",    envir=environment())

  }

  #choiceModel <- match.arg(choiceModel)

  # ChoiceModel must be one of the list
  if(!any(choiceModel == c("model1", "model2", "model3", "model4")))
    stop("choiceModel must be one of the following: 'model1', 'model2', 'model3', or 'model4'")

  # Check if length of par is appropriate for each model
  if(is.null(par) & choiceModel == "model1") {par <- c(rep(0, ncol(data)-1), rep(-1, ncol(data)-1))}
  else if(!is.null(par) & choiceModel == "model1"){if(length(par) != 2*(ncol(data)-1)) stop(paste0("ERROR: par must be length: ", 2*(ncol(data)-1), ", not: ", length(par)))}
  if(is.null(par) & choiceModel == "model2") {par <- c(0, rep(-1, ncol(data)))}
  else if(!is.null(par) & choiceModel == "model2"){if(length(par) != ncol(data)) stop(paste0("ERROR: par must be length: ", ncol(data), ", not: ", length(par)))}
  if(is.null(par) & choiceModel == "model3") {par <- c(0, rep(-1, ncol(data)-1))}
  else if(!is.null(par) & choiceModel == "model3"){if(length(par) != ncol(data)) stop(paste0("ERROR: par must be length: ", ncol(data), ", not: ", length(par)))}
  if(is.null(par) & choiceModel == "model4") {par <- c(0, -1)}
  else if(!is.null(par) & choiceModel == "model4"){if(length(par) != 2) stop(paste0("ERROR: par must be length: ", 2, ", not: ", length(par)))}


  if(length(par) != length(lowerBound) && method == "L-BFGS-B"){

    warning("'lowerBound' must be length: ",length(par), ", not: ", length(lowerBound))
    lowerBound <- rep(-Inf, length(par))

  }

  if(method == "Nelder-Mead") {
    maxit = maxit
  }

  opt_ma <- optim(par=par, fn=batchMarkHmmLL, data=data, choiceModel=choiceModel,
                 cores=cores, method=method,lower=lowerBound, control = control, hessian = hessian,...)

  res <- list()

  if(choiceModel == "model1") {

    res$phi <- round(batchLogit(opt_ma$par[1:(ncol(data)-1)]),2)
    res$p   <- round(batchLogit(opt_ma$par[ncol(data):(2*(ncol(data)-1))]),2)
    res$ll  <- round(opt_ma$value, 2)
    res$AIC <- round(2*opt_ma$value+2*length(opt_ma$par),2)

  }

  else if(choiceModel == "model2") {

    res$phi <- round(batchLogit(opt_ma$par[2:ncol(data)]),2)
    res$p   <- round(batchLogit(opt_ma$par[1]),2)
    res$ll  <- round(opt_ma$value, 2)
    res$AIC <- round(2*opt_ma$value+2*length(opt_ma$par),2)

  }

  else if(choiceModel == "model3") {

    res$phi <- round(batchLogit(opt_ma$par[1]),2)
    res$p   <- round(batchLogit(opt_ma$par[-1]),2)
    res$ll  <- round(opt_ma$value, 2)
    res$AIC <- round(2*opt_ma$value+2*length(opt_ma$par),2)

  }

  else {

    res$phi <- round(batchLogit(opt_ma$par[1]),2)
    res$p   <- round(batchLogit(opt_ma$par[2]),2)
    res$ll  <- round(opt_ma$value,2)
    res$AIC <- round(2*opt_ma$value+2*length(opt_ma$par),2)

  }

  if(hessian) res$hessian <- opt_ma$hessian
  if(parallel) {

    parallel::stopCluster(cl = cluster)
    doParallel::stopImplicitCluster()
  }

  return(res)

}


##################################################################################

#' @title Combined Marked and Unmarked models.

#' @description batchMarkUnmarkOptim function provides the batch marking and unmarked function to be optimized.
#'
#' @param par Initial values for the parameters to be optimized over.
#' @param data A capture-recapture data matrix or data frame
#' @param choiceModel This chooses among different models and allow for model selection
#' @param cores The number of cores for parallelization
#' @param method The method to be used. See optim for details.
#' @param popSize The Horvitz_Thompson method or Model-Based to compute population size.
#' @param lowerBound Lower bounds on the variables for the "L-BFGS-B" method.
#' @param parallel Logical. Should the algorithm be run in parallel? This will be implemented in a future version.
#' @param Umax The maximum number of the unmarked individuals in the population for capture on any occasion.
#' @param nBins The number of bin size into which the matrix will be divided.
#' @param control a list of control parameters. See \code{\link{optim}} for details.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' @param ... Further arguments to be passed by user which goes into the optim function.
#' @return A list of the following optimized parameters will be returned.
#' \describe{
#'  \item{phi}{The survival probability and remaining in the population between occasion t and t+1.}
#'  \item{p}{The capture probability at occasion time t.}
#'  \item{ll}{The optimized log-likelihood value of marked model.}
#'  \item{hessian}{The hessian matrix.}
#'  \item{AIC}{The Akaike Information Criteria for model selection.}
#'  \item{lambda}{Initial mean abundance at occasion t = 1.}
#'  \item{gam}{Recruitment rate of individual into the unmarked population.}
#'  \item{M}{Total number of marked individual in the population.}
#'  \item{U}{Total number of unmarked individuals in the population available for capture at occasion t = 1,..., T.}
#'  \item{N}{Total population size at time t = 1, ..., T.}
#' }

#' @details
#' Note that arguments after ... must be matched exactly.
#'
#' batchMarkUnmarkOptim depends on optim function to optimize the parameters of the combined model. By default optim performs minimization.
#'
#' Example on Umax and nBins: Umax = 1800 has a matrix of 1801 x 1801 and nBins = 20, reduces the matrix to 90 x 90. This is done in Cowen et al., 2017 to reduce the computing time when dealing with large matrix.
#'
#' @references Laura L. E. Cowen, Panagiotis Besbeas, Byron J. T. Morgan, 2017.: Hidden Markov Models for Extended Batch Data,
#' Biometrics, 73, 1321-1331. DOI: 10.1111/biom.12701.
#'
#' @rdname batchMarkUnmarkOptim
#' @importFrom stats optim dpois ppois
#' @importFrom parallel makeCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom utils install.packages installed.packages
#' @export

#' @examples
#' # Load the package
#' library(extBatchMarking)
#'
#' # Load the WeatherLoach data from Cowen et al., 2017.
#' data(WeatherLoach)
#'
#' # Initial parameter values
#' theta <- c(0.1, 0.1, 7, -1.5)
#'
#' \donttest{
#' mod1 <- batchMarkUnmarkOptim(
#'            par         = theta,
#'            data        = WeatherLoach,
#'            Umax        = 1800,
#'            nBins       = 20,
#'            choiceModel = "model4",
#'            popSize    = "Horvitz_Thompson",
#'            method      = "CG",
#'            parallel    = FALSE,
#'            control     = list(trace = 1))
#'
#'  # Survival probability
#'  mod1$phi
#'  # Capture probability
#'  mod1$p
#'  # Optimized log-likelihood
#'  mod1$ll
#'  # The Hessian matrix
#'  mod1$hessian
#'  # The Aikaike Information Criteria
#'  mod1$AIC
#'  # The initial mean abundance
#'  mod1$lambda
#'  # Recruitment rate into the population
#'  mod1$gam
#'  # The estimated abundance of unmarked animals
#'  mod1$U
#'  # The estimated abundance of marked animals
#'  mod1$M
#'  # The estimated total abundance of marked and unmarked animals
#'  mod1$N
#'  }
#'
#'  \donttest{
#' mod2 <- batchMarkUnmarkOptim(
#'            par         = theta,
#'            data        = WeatherLoach,
#'            Umax        = 1800,
#'            nBins       = 20,
#'            choiceModel = "model4",
#'            popSize    = "Model-Based",
#'            method      = "L-BFGS-B",
#'            parallel    = FALSE,
#'            control     = list(trace = 1))
#'
#'  # Survival probability
#'  mod2$phi
#'  # Capture probability
#'  mod2$p
#'  # Optimized log-likelihood
#'  mod2$ll
#'  # The Hessian matrix
#'  mod2$hessian
#'  # The Akaike Information Criteria
#'  mod2$AIC
#'  # The initial mean abundance
#'  mod2$lambda
#'  # Recruitment rate into the population
#'  mod2$gam
#'  # The estimated abundance of unmarked animals
#'  mod2$U
#'  # The estimated abundance of marked animals
#'  mod2$M
#'  # The estimated total abundance of marked and unmarked animals
#'  mod2$N
#'  }


batchMarkUnmarkOptim <- function(par=NULL, data, choiceModel=c("model1", "model2", "model3", "model4"),
                                 method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B"), Umax=1800, nBins=20,
                                 popSize=c("Horvitz_Thompson", "Model-Based"), parallel=FALSE, lowerBound=-Inf,
                                 cores=1, hessian=FALSE, control,...){

  # my_packages   <- c("Rcpp", "parallel", "optimParallel", "RcppArmadillo")                # Specify your packages
  # not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]    # Extract not installed packages
  # if(length(not_installed)) utils::install.packages(not_installed)

  if(parallel) {

    cores   <- parallel::detectCores()-1
    cluster <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cluster, cores = cores)
    parallel::clusterExport(cluster, "batchLogit",           envir=environment())
    parallel::clusterExport(cluster, "probs",                envir=environment())
    parallel::clusterExport(cluster, "delta_g",              envir=environment())
    parallel::clusterExport(cluster, "batchLL",              envir=environment())
    parallel::clusterExport(cluster, "batchMarkHmmLL",       envir=environment())
    parallel::clusterExport(cluster, "gamma_gt",             envir=environment())
    parallel::clusterExport(cluster, "data",                 envir=environment())
    parallel::clusterExport(cluster, "%dopar%",              envir=environment())
    parallel::clusterExport(cluster, "%do%",                 envir=environment())
    parallel::clusterExport(cluster, "%:%",                  envir=environment())
    parallel::clusterExport(cluster, "batchUnmarkHmmLL",     envir=environment())
    parallel::clusterExport(cluster, "batchMarkUnmarkHmmLL", envir=environment())
    parallel::clusterExport(cluster, "choiceModel",          envir = environment())
    parallel::clusterExport(cluster, "popSize",              envir = environment())
    parallel::clusterExport(cluster, "batchUnmark2Viterbi",  envir = environment())


  }

  # popSize must be one of the list
  if(!any(popSize == c("Horvitz_Thompson", "Model-Based")))
    stop("popSize must be one of the following: 'Horvitz_Thompson', 'Model-Based'")

  # ChoiceModel must be one of the list
  if(!any(choiceModel == c("model1", "model2", "model3", "model4")))
    stop("choiceModel must be one of the following: 'model1', 'model2', 'model3', or 'model4'")

  if(is.null(par) & choiceModel == "model1") {par <- c(rep(0.1, ncol(data)+ncol(data)-1), 7, -1.5)}
  else if(!is.null(par) & choiceModel == "model1"){if(length(par) != ncol(data)+ncol(data)+1) stop(paste0("ERROR: par must be length: ", ncol(data)+ncol(data)+1, ", not: ", length(par)))}
  if(is.null(par) & choiceModel == "model2") {par <- c(rep(0.1, ncol(data)), 7, -1.5)}
  else if(!is.null(par) & choiceModel == "model2"){if(length(par) != ncol(data)+2) stop(paste0("ERROR: par must be length: ", ncol(data)+2, ", not: ", length(par)))}
  if(is.null(par) & choiceModel == "model3") {par <- c(rep(0.1, ncol(data)+1), 7, -1.5)}
  else if(!is.null(par) & choiceModel == "model3"){if(length(par) != ncol(data)+3) stop(paste0("ERROR: par must be length: ", ncol(data)+3, ", not: ", length(par)))}
  if(is.null(par) & choiceModel == "model4") {par <- c(rep(0.1, 2), 7, -1.5)}
  else if(!is.null(par) & choiceModel == "model4"){if(length(par) != 2*2) stop(paste0("ERROR: par must be length: ", 2*2, ", not: ", length(par)))}

  if(length(par) != length(lowerBound) && method == "L-BFGS-B"){

    warning("'lowerBound' must be length: ",length(par), ", not: ", length(lowerBound))
    lowerBound <- c(rep(-Inf, (length(par)-2)), 1, -Inf)

  }

  if(method == "Nelder-Mead") {
    maxit = maxit
  }

  opt_ma <- optim(par=par, fn=batchMarkUnmarkHmmLL, data=data, choiceModel=choiceModel, Umax=Umax, nBins=nBins,
                  cores=cores, method=method, lower=lowerBound, control=control, hessian=hessian, ...)

  res <- list()

  if(choiceModel == "model1") {

    res$phi    <- round(batchLogit(opt_ma$par[1:(ncol(data)-1)]),2)
    res$p      <- round(batchLogit(opt_ma$par[ncol(data):(2*ncol(data)-1)]),2)
    res$lambda <- round(exp(opt_ma$par[2*ncol(data)]),2)
    res$gam    <- round(exp(opt_ma$par[2*ncol(data)+1]),2)
    res$ll     <- round(opt_ma$value,2)
    res$AIC    <- round(2*opt_ma$value+2*length(opt_ma$par),2)

  }

  else if(choiceModel == "model2") {

    res$phi    <- round(batchLogit(opt_ma$par[2:ncol(data)]),2)
    res$p      <- round(batchLogit(opt_ma$par[1]),2)
    res$lambda <- round(exp(opt_ma$par[ncol(data)+1]),2)
    res$gam    <- round(exp(opt_ma$par[ncol(data)+2]),2)
    res$ll     <- round(opt_ma$value, 2)
    res$AIC    <- round(2*opt_ma$value+2*length(opt_ma$par),2)

  }

  else if(choiceModel == "model3") {

    res$phi    <- round(batchLogit(opt_ma$par[1]),2)
    res$p      <- round(batchLogit(opt_ma$par[2:(ncol(data)+1)]),2)
    res$lambda <- round(exp(opt_ma$par[ncol(data)+2]),2)
    res$gam    <- round(exp(opt_ma$par[ncol(data)+3]),2)
    res$ll     <- round(opt_ma$value, 2)
    res$AIC    <- round(2*opt_ma$value+2*length(opt_ma$par),2)

  }

  else {

    res$phi    <- round(batchLogit(opt_ma$par[1]),2)
    res$p      <- round(batchLogit(opt_ma$par[2]),2)
    res$lambda <- round(exp(opt_ma$par[3]),2)
    res$gam    <- round(exp(opt_ma$par[4]),2)
    res$ll     <- round(opt_ma$value, 2)
    res$AIC    <- round(2*opt_ma$value+2*length(opt_ma$par),2)

  }


  res$U <- batchUnmark2Viterbi(par = opt_ma$par, data = data, Umax = Umax, nBins = nBins,
                               choiceModel = choiceModel)

  # Total abundance: Marked + Unmarked abundance
  if(popSize == "Horvitz_Thompson"){

    dat   <- cbind(rep(0, nrow(data[-1,-1])), data[-1,-1])
    res$M <- colSums(dat)/res$p

  }else{

    phi <- res$phi

    if(length(phi) == 1) phi <- rep(phi,(ncol(data)-1))

    R     <- data[-1,1]
    Times <- ncol(data)
    M     <- matrix(0, nrow = length(R), ncol = Times)

    for (times in 2:Times) {
      for (i in 1:(times-1)) M[i,times] <- R[i] * prod(phi[i:(times-1)])
    }

    res$M <- colSums(M)
  }

  res$N <- res$U + res$M

  if(hessian) res$hessian <- opt_ma$hessian
  if(parallel) {

    parallel::stopCluster(cl = cluster)
    doParallel::stopImplicitCluster()
  }

  return(res)

}
