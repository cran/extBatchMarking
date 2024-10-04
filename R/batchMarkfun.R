#' @title Marked model only.

#' @description batchMarkOptim function optimizes \code{\link{batchMarkHmmLL}} function.
#'
#' @param par Initial values for the parameters to be optimized over.
#' @param data A capture-recapture data matrix or data frame
#' @param choiceModel This chooses among different models and allow for model selection
#' @param covariate_phi This covariate placeholder for the parameter phi_t
#' @param covariate_p This covariate placeholder for the parameter p_t
#' @param method The method to be used. See \code{\link{optim}} for details.
#' @param lowerBound Lower bounds on the variables for the "L-BFGS-B" \code{method}.
#' @param control A list of control parameters. See optim for details.
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
#'  \item{SE}{The standard error for each parameter.}
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
#'            par           = theta,
#'            data          = WeatherLoach,
#'            choiceModel   = "model4",
#'            method        = "BFGS",
#'            control       = list(trace = 1),
#'            covariate_phi = NULL,
#'            covariate_p   = NULL)
#'
#'  # print(mod1)
#'
#'  # Survival probability
#'  mod1$phi
#'  # Capture probability
#'  mod1$p
#'  # Optimized log-likelihood
#'  mod1$ll
#'  # The Aikaike Information Criteria
#'  mod1$AIC
#'  }
#'
#'  \donttest{
#'  mod2 <- batchMarkOptim(
#'            par           = theta,
#'            data          = WeatherLoach,
#'            choiceModel   = "model4",
#'            method        = "L-BFGS-B",
#'            control       = list(trace = 1),
#'            covariate_phi = NULL,
#'            covariate_p   = NULL)
#'
#'  # print(mod2)
#'  # Survival probability
#'  mod2$phi
#'  # Capture probability
#'  mod2$p
#'  # Optimized log-likelihood
#'  mod2$ll
#'  # The Akaike Information Criteria
#'  mod2$AIC
#'  }


batchMarkOptim <- function(par=NULL, data, covariate_phi = NULL, covariate_p = NULL, choiceModel=c("model1", "model2", "model3", "model4"),
                           method=c("Nelder-Mead","BFGS", "CG", "L-BFGS-B"), lowerBound=-Inf, control, ...){

  # ChoiceModel must be one of the list
  if(!any(choiceModel == c("model1", "model2", "model3", "model4")))
    stop("choiceModel must be one of the following: 'model1', 'model2', 'model3', or 'model4'")

  # Check if length of par is appropriate for each model
  if(is.null(par) & choiceModel == "model1") {par <- c(rep(0, ncol(data)-1), rep(-1, ncol(data)-1))}
  else if(!is.null(par) & choiceModel == "model1"){if(length(par) != 2*(ncol(data)-1)) stop(paste0("par must be length: ", 2*(ncol(data)-1), ", not: ", length(par)))}
  if(is.null(par) & choiceModel == "model2") {par <- c(0, rep(-1, ncol(data)))}
  else if(!is.null(par) & choiceModel == "model2"){if(length(par) != ncol(data)) stop(paste0("par must be length: ", ncol(data), ", not: ", length(par)))}
  if(is.null(par) & choiceModel == "model3") {par <- c(0, rep(-1, ncol(data)-1))}
  else if(!is.null(par) & choiceModel == "model3"){if(length(par) != ncol(data)) stop(paste0("par must be length: ", ncol(data), ", not: ", length(par)))}
  if(is.null(par) & choiceModel == "model4") {par <- c(0, -1)}
  else if(!is.null(par) & choiceModel == "model4"){if(length(par) != 2) stop(paste0("par must be length: ", 2, ", not: ", length(par)))}


  if(length(par) != length(lowerBound) && method == "L-BFGS-B"){

    warning("'lowerBound' must be length: ",length(par), ", not: ", length(lowerBound))
    lowerBound <- rep(-Inf, length(par))

  }

  if(method == "Nelder-Mead") {
    maxit = maxit
  }

  opt_ma <- optim(par=par, fn=batchMarkHmmLL, data=data, covariate_phi = covariate_phi, covariate_p = covariate_p,
                  choiceModel=choiceModel, method=method,lower=lowerBound, control = control, hessian = TRUE,...)

  res <- list()

  if(choiceModel == "model1") {

    res$phi <- round(batchLogit(opt_ma$par[1:(ncol(data)-1)]),4)
    res$p   <- round(batchLogit(opt_ma$par[ncol(data):(2*(ncol(data)-1))]),4)
    res$ll  <- round(opt_ma$value, 4)
    res$AIC <- round(2*opt_ma$value+2*length(opt_ma$par),4)

  }

  else if(choiceModel == "model2") {

    res$phi <- round(batchLogit(opt_ma$par[2:ncol(data)]),4)
    res$p   <- round(batchLogit(opt_ma$par[1]),4)
    res$ll  <- round(opt_ma$value, 4)
    res$AIC <- round(2*opt_ma$value+2*length(opt_ma$par),4)

  }

  else if(choiceModel == "model3") {

    res$phi <- round(batchLogit(opt_ma$par[1]),4)
    res$p   <- round(batchLogit(opt_ma$par[-1]),4)
    res$ll  <- round(opt_ma$value, 4)
    res$AIC <- round(2*opt_ma$value+2*length(opt_ma$par),4)

  }

  else {

    res$phi <- round(batchLogit(opt_ma$par[1]),4)
    res$p   <- round(batchLogit(opt_ma$par[2]),4)
    res$ll  <- round(opt_ma$value,4)
    res$AIC <- round(2*opt_ma$value+2*length(opt_ma$par),4)

  }

  serr_trans <- sqrt(diag(solve(opt_ma$hessian)))
  serr       <- c(res$phi, res$p) * (1 - c(res$phi, res$p)) * serr_trans

  if(choiceModel == "model1"){

    res$SE_phi <- round(serr[1:(ncol(data)-1)], 4)
    res$SE_p   <- round(serr[ncol(data):(2*(ncol(data)-1))], 4)

  }else if(choiceModel == "model2"){

    res$SE_phi <- round(serr[2:ncol(data)], 4)
    res$SE_p   <- round(serr[1], 4)

  }else if(choiceModel == "model3"){

    res$SE_phi <- round(serr[1], 4)
    res$SE_p   <- round(serr[-1], 4)

  }else{

    res$SE_phi <- round(serr[1], 4)
    res$SE_p   <- round(serr[2], 4)

  }

  # Goodness of Fit ---------------------------------------------------------

  # Compute Q_gt

  phi <- res$phi
  p   <- res$p

  if(length(phi) == 1) phi <- rep(phi,(ncol(data)-1))
  if(length(p)   == 1) p   <- rep(p, ncol(data))

  R     <- data[-1,1]
  Times <- ncol(data)
  Q_gt     <- matrix(0, nrow = length(R), ncol = Times)

  for (times in 2:Times) {
    for (i in 1:(times-1)) Q_gt[i,times] <- p[times] * prod(phi[i:(times-1)])
  }

  exp_rgt   <- R * Q_gt[, 2:ncol(Q_gt)]
  r_gt      <- data[2:nrow(data), 2:ncol(data)]
  rgt_error <- (r_gt - exp_rgt) / sqrt(exp_rgt)

  res$exp_rgt   <- exp_rgt
  res$rgt_error <- rgt_error

  class(res) <- "batchMarkOptim"

  return(res)

}


##################################################################################

#' @title Combined Marked and Unmarked models.

#' @description batchMarkUnmarkOptim function optimizes \code{\link{batchMarkUnmarkHmmLL}} function.
#'
#' @param par Initial values for the parameters to be optimized over.
#' @param data A capture-recapture data matrix or data frame
#' @param choiceModel This chooses among different models and allow for model selection
#' @param covariate_phi This covariate placeholder for the parameter phi_t
#' @param covariate_p This covariate placeholder for the parameter p_t
#' @param method The method to be used. See optim for details.
#' @param popSize The Horvitz_Thompson method or Model-Based to compute population size.
#' @param lowerBound Lower bounds on the variables for the "L-BFGS-B" method.
#' @param Umax The maximum number of the unmarked individuals in the population for capture on any occasion.
#' @param nBins The number of bin size into which the matrix will be divided.
#' @param control a list of control parameters. See \code{\link{optim}} for details.
#' @param ... Further arguments to be passed by user which goes into the optim function.
#' @return A list of the following optimized parameters will be returned.
#' \describe{
#'  \item{phi}{The survival probability and remaining in the population between occasion t and t+1.}
#'  \item{p}{The capture probability at occasion time t.}
#'  \item{ll}{The optimized log-likelihood value of marked model.}
#'  \item{SE}{The standard error for each parameter.}
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
#' @importFrom stats optim dpois ppois pbinom qnorm qqline qqnorm
#' @importFrom graphics segments par
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
#'            nBins       = 600,
#'            covariate_phi = NULL,
#'            covariate_p   = NULL,
#'            choiceModel = "model4",
#'            popSize    = "Horvitz_Thompson",
#'            method      = "CG",
#'            control     = list(trace = 1))
#'
#'  # Survival probability
#'  mod1$phi
#'  # Capture probability
#'  mod1$p
#'  # Optimized log-likelihood
#'  mod1$ll
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
#'            nBins       = 600,
#'            choiceModel = "model4",
#'            covariate_phi = NULL,
#'            covariate_p   = NULL,
#'            popSize    = "Model-Based",
#'            method      = "L-BFGS-B",
#'            control     = list(trace = 1))
#'
#'  # print(mod2)
#'  # plot(mod2)
#'  # Survival probability
#'  mod2$phi
#'  # Capture probability
#'  mod2$p
#'  # Optimized log-likelihood
#'  mod2$ll
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
                                 covariate_phi = NULL, covariate_p = NULL,
                                 method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B"), Umax=1800, nBins=20,
                                 popSize=c("Horvitz_Thompson", "Model-Based"), lowerBound=-Inf, control,...){

  # popSize must be one of the list
  if(!any(popSize == c("Horvitz_Thompson", "Model-Based")))
    stop("popSize must be one of the following: 'Horvitz_Thompson', 'Model-Based'")

  # ChoiceModel must be one of the list
  if(!any(choiceModel == c("model1", "model2", "model3", "model4")))
    stop("choiceModel must be one of the following: 'model1', 'model2', 'model3', or 'model4'")

  if(is.null(par) & choiceModel == "model1") {par <- c(rep(0.1, ncol(data)+ncol(data)-1), 7, -1.5)}
  else if(!is.null(par) & choiceModel == "model1"){if(length(par) != ncol(data)+ncol(data)+1) stop(paste0("par must be length: ", ncol(data)+ncol(data)+1, ", not: ", length(par)))}
  if(is.null(par) & choiceModel == "model2") {par <- c(rep(0.1, ncol(data)), 7, -1.5)}
  else if(!is.null(par) & choiceModel == "model2"){if(length(par) != ncol(data)+2) stop(paste0("par must be length: ", ncol(data)+2, ", not: ", length(par)))}
  if(is.null(par) & choiceModel == "model3") {par <- c(rep(0.1, ncol(data)+1), 7, -1.5)}
  else if(!is.null(par) & choiceModel == "model3"){if(length(par) != ncol(data)+3) stop(paste0("par must be length: ", ncol(data)+3, ", not: ", length(par)))}
  if(is.null(par) & choiceModel == "model4") {par <- c(rep(0.1, 2), 7, -1.5)}
  else if(!is.null(par) & choiceModel == "model4"){if(length(par) != 2*2) stop(paste0("par must be length: ", 2*2, ", not: ", length(par)))}

  if(length(par) != length(lowerBound) && method == "L-BFGS-B"){

    warning("'lowerBound' must be length: ",length(par), ", not: ", length(lowerBound))
    lowerBound <- c(rep(-Inf, (length(par)-2)), 1, -Inf)

  }

  if(method == "Nelder-Mead") {
    maxit = maxit
  }

  opt_ma <- optim(par=par, fn=batchMarkUnmarkHmmLL, data=data, covariate_phi = covariate_phi, covariate_p = covariate_p,
                  choiceModel=choiceModel, Umax=Umax, nBins=nBins, method=method, lower=lowerBound, control=control, hessian=TRUE, ...)

  res <- list()

  if(choiceModel == "model1") {

    res$phi    <- round(batchLogit(opt_ma$par[1:(ncol(data)-1)]),4)
    res$p      <- round(batchLogit(opt_ma$par[ncol(data):(2*ncol(data)-1)]),4)
    res$lambda <- round(exp(opt_ma$par[2*ncol(data)]),4)
    res$gam    <- round(exp(opt_ma$par[2*ncol(data)+1]),4)
    res$ll     <- round(opt_ma$value,2)
    res$AIC    <- round(2*opt_ma$value+2*length(opt_ma$par),4)

  }

  else if(choiceModel == "model2") {

    res$phi    <- round(batchLogit(opt_ma$par[2:ncol(data)]),4)
    res$p      <- round(batchLogit(opt_ma$par[1]),4)
    res$lambda <- round(exp(opt_ma$par[ncol(data)+1]),4)
    res$gam    <- round(exp(opt_ma$par[ncol(data)+2]),4)
    res$ll     <- round(opt_ma$value, 4)
    res$AIC    <- round(2*opt_ma$value+2*length(opt_ma$par),4)

  }

  else if(choiceModel == "model3") {

    res$phi    <- round(batchLogit(opt_ma$par[1]),4)
    res$p      <- round(batchLogit(opt_ma$par[2:(ncol(data)+1)]),4)
    res$lambda <- round(exp(opt_ma$par[ncol(data)+2]),4)
    res$gam    <- round(exp(opt_ma$par[ncol(data)+3]),4)
    res$ll     <- round(opt_ma$value, 4)
    res$AIC    <- round(2*opt_ma$value+2*length(opt_ma$par),4)

  }

  else {

    res$phi    <- round(batchLogit(opt_ma$par[1]),4)
    res$p      <- round(batchLogit(opt_ma$par[2]),4)
    res$lambda <- round(exp(opt_ma$par[3]),4)
    res$gam    <- round(exp(opt_ma$par[4]),4)
    res$ll     <- round(opt_ma$value, 4)
    res$AIC    <- round(2*opt_ma$value+2*length(opt_ma$par),4)

  }


  res$U <- batchUnmark2Viterbi(par = opt_ma$par, data = data, Umax = Umax, nBins = nBins, choiceModel = choiceModel)


  # Total Abundance ------------------------------------------------------------


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

  serr_trans <- sqrt(diag(solve(opt_ma$hessian)))
  serr       <- c((c(res$phi, res$p) * (1 - c(res$phi, res$p))), res$lambda, res$gam) * serr_trans
  res$SE     <- serr

  if(choiceModel == "model1"){

    res$SE_phi    <- round(serr[1:(ncol(data)-1)], 4)
    res$SE_p      <- round(serr[ncol(data):(2*(ncol(data)-1))], 4)
    res$SE_lambda <- round(serr[2*ncol(data)], 4)
    res$SE_gam    <- round(serr[2*ncol(data)+1], 4)

  }else if(choiceModel == "model2"){

    res$SE_phi    <- round(serr[2:ncol(data)], 4)
    res$SE_p      <- round(serr[1], 4)
    res$SE_lambda <- round(serr[ncol(data)+1], 4)
    res$SE_gam    <- round(serr[ncol(data)+2], 4)

  }else if(choiceModel == "model3"){

    res$SE_phi    <- round(serr[1], 4)
    res$SE_p      <- round(serr[2:(ncol(data)+1)], 4)
    res$SE_lambda <- round(serr[ncol(data)+2], 4)
    res$SE_gam    <- round(serr[ncol(data)+3], 4)

  }else{

    res$SE_phi    <- round(serr[1], 4)
    res$SE_p      <- round(serr[2], 4)
    res$SE_lambda <- round(serr[3], 4)
    res$SE_gam    <- round(serr[4], 4)

  }

  # Goodness of Fit ---------------------------------------------------------

  # Compute Q_gt

  phi <- res$phi
  p   <- res$p

  if(length(phi) == 1) phi <- rep(phi,(ncol(data)-1))
  if(length(p)   == 1) p   <- rep(p, ncol(data))

  R     <- data[-1,1]
  Times <- ncol(data)
  Q_gt     <- matrix(0, nrow = length(R), ncol = Times)

  for (times in 2:Times) {
    for (i in 1:(times-1)) Q_gt[i,times] <- p[times] * prod(phi[i:(times-1)])
  }

  exp_rgt   <- R * Q_gt[, 2:ncol(Q_gt)]
  r_gt      <- data[2:nrow(data), 2:ncol(data)]
  rgt_error <- (r_gt - exp_rgt) / sqrt(exp_rgt)

  # For the Unmarked
  u <- data[1,] - c(0, colSums(r_gt))

  F_lower <- vector(length = length(u))
  F_upper <- vector(length = length(u))

  p <- res$p
  if(length(p) == 1) p <- rep(p, ncol(data))

  for (t in 1:length(u)) {

    F_lower[t] <- pbinom(round(u[t])-1, size = round(res$U[t]), prob = p[t])
    F_upper[t] <- 1 - pbinom(round(u[t]), size = round(res$U[t]), prob = p[t])

  }

  z_lower       <- qnorm(F_lower)
  z_upper       <- qnorm(F_upper)
  low_seg       <- qqnorm(z_lower, plot.it = FALSE)
  low_mat       <- cbind(X=low_seg$x, Y=low_seg$y)
  res$low_clean <- low_mat[order(low_mat[,1]),]
  up_seg        <- qqnorm(z_upper, plot.it = FALSE)
  up_mat        <- cbind(X=up_seg$x, Y=up_seg$y)
  res$up_clean  <- up_mat[order(up_mat[,1]),]
  res$z_ave     <- (z_lower[order(low_mat[,1])] + z_upper[order(up_mat[,1])]) / 2

  res$exp_rgt   <- exp_rgt
  res$rgt_error <- rgt_error


  class(res) <- "batchMarkUnmarkOptim"

  return(res)

}


#' Print Method for batchMarkOptim Objects
#'
#' This function defines how objects of class "Employee" are printed.
#'
#' @param x An object of class "Employee".
#' @param ... Additional arguments passed to the print method.
#' @rdname print.batchMarkOptim

#' @export
print.batchMarkOptim <- function(x, ...) {

  tt1 <- data.frame("log-likelihood" = x$ll, "AIC" =  x$AIC)

  if(length(x$p) == 1 & length(x$phi) != 1) {

    x$p    <- rep(x$p, length(x$phi))
    x$SE_p <- rep(x$SE_p, length(x$phi))
  }

  if(length(x$phi) == 1 & length(x$p) != 1){

    x$phi    <- rep(x$phi, length(x$p))
    x$SE_phi <- rep(x$SE_phi, length(x$p))
  }

  tt2 <- data.frame("p" = x$p, "p_S Error" = x$SE_p, "phi" = x$phi, "phi_S Error" = x$SE_phi)
  tt3 <- list(knitr::kable(tt1), knitr::kable(tt2))

  return(tt3)
}

#' Print Method for batchMarkUnmarkOptim Objects
#'
#' This function defines how objects of class "Employee" are printed.
#'
#' @param x An object of class "Employee".
#' @param ... Additional arguments passed to the print method.
#' @rdname print.batchMarkUnmarkOptim
#' @export

print.batchMarkUnmarkOptim <- function(x, ...) {

  tt1 <- data.frame("log-likelihood" = x$ll, "AIC" =  x$AIC)

  if(length(x$p) == 1 & length(x$phi) != 1) {

    x$p       <- rep(x$p, length(x$phi))
    x$SE_p    <- rep(x$SE_p, length(x$phi))

  }

  if(length(x$phi) == 1 & length(x$p) != 1){

    x$phi     <- rep(x$phi, length(x$p))
    x$SE_phi  <- rep(x$SE_phi, length(x$p))

  }

  tt2 <- data.frame("p" = x$p, "p_S Error" = x$SE_p, "phi" = x$phi, "phi_S Error" = x$SE_phi)
  tt3 <- data.frame("No of Unmarked(U)" = round(x$U), "No of Marked(M)" = round(x$M), "Abundance(N)" = round(x$N))
  tt4 <- data.frame("Lambda" = x$lambda, "Lambda_S.Error" = x$SE_lambda, "Gam" = x$gam, "Gam_S.Error" = x$SE_gam)
  tt5 <- list(knitr::kable(tt1), knitr::kable(tt2), knitr::kable(tt3), knitr::kable(tt4))

  return(tt5)

}


#' Plot Method for batchMarkUnmarkOptim Objects
#'
#' This function defines how objects of class "Employee" are printed.
#'
#' @param x An object of class "Employee".
#' @param ... Additional arguments passed to the print method.
#' @rdname plot.batchMarkOptim
#' @export

plot.batchMarkOptim <- function(x, ...) {

  #par(mfrow = c(1,2))

  plot(x=as.vector(log(x$exp_rgt)), y=as.vector(x$rgt_error), type="p",
       xlab="log(Expected)", ylab = "(Observed-Expected) / sqrt(Expected)",
       main="Goodness-of-fit plots", pch=1, cex=1)

}


#' Plot Method for batchMarkUnmarkOptim Objects
#'
#' This function defines how objects of class "Employee" are printed.
#'
#' @param x An object of class "Employee".
#' @param ... Additional arguments passed to the print method.
#' @rdname plot.batchMarkUnmarkOptim
#' @export

plot.batchMarkUnmarkOptim <- function(x, ...) {

  par(mfrow = c(1,2))

  plot(x=as.vector(log(x$exp_rgt)), y=as.vector(x$rgt_error), type="p",
       xlab="log(Expected)", ylab = "(Observed-Expected) / sqrt(Expected)",
       main="Goodness-of-fit plots", pch=1, cex=1)

  {plot(range(x$low_clean[,1], x$up_clean[,1]), range(x$low_clean[,2], x$up_clean[,2]), type = "n",
        xlab="Theoretical Quantities", ylab="Sample Quantities", main="QQ plot of the unmarked")
    segments(x0 = x$low_clean[,1],y0 = x$low_clean[,2], x1 = x$up_clean[,1], y1 = x$up_clean[,2])
    qqline(x$z_ave)}

}

