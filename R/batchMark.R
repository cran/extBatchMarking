# '%dopar%' <- foreach::`%do%`
# '%do%'    <- foreach::`%dopar%`
# '%:%'     <- foreach::`%:%`
# '%knitr%'   <- knitr::'%knitr%'
# '%installed.packages%' <- utils::'%installed.packages%'
# '%install.packages%'   <- utils::'%install.packages%'
# '%ppois%'   <- stats::'%ppois%'
# '%dpois%'   <- stats::'%dpois%'
# '%makeCluster%' <- parallel::'%makeCluster%'
# '%stopCluster%' <- parallel::'%stopCluster%'
# '%registerDoParallel%' <- doParallel::'%registerDoParallel%'
# '%stopImplicitCluster%' <- doParallel::'%stopImplicitCluster%'


#'@title batchLogit function

#' @description
#' 'batchLogit' provides the number between 0 and 1.
#' @param x This is an input numerical value i.e double.
#' @return Returns a number between 0 and 1.
#' @rdname batchLogit
#'
#'

# Logit function for the phi's and probs
batchLogit <- function(x){1 / (1+exp(-x))}; # H(x)


############################################

#------------------------------------------
# Transition Probability Matrices from Rcpp
#-----------------------------------------

# Check the bmarkedhmmv13.Rmd



############################################
#-------------------------------------------
# State-dependent probability matrices
#-------------------------------------------


#' State-dependent probability function
#'
#' 'probs' computes the state-dependent transition matrix
#'
#' @param r The number of individuals from batch group "g" recaptured at recapture occasion t; g = 1,2,...,G, t = g+1,...,T. This must be an integer.
#' @param p The probability of capture at occasion t. This must be a number between 0 and 1.
#' @param R The number of individuals marked and released at sampling occasion g from batch group g; g = 1,2,...,G. This must be an integer.
#' @return PR diagonal matrix of the state-dependent probability.
#' @rdname probs




probs <- function(r, p, R){

  kl <- rep(-Inf, R[1]+1)
  i <- r:R[1];
  nku1 <- lgamma(i+1) - lgamma(r+1) - lgamma(i-r+1)
  kl[i+1] <- nku1 + r*log(p) + (i-r)*log1p(-p)
  kl[which(is.na(kl))] <- 0
  PR <- diag(exp(kl))
  return(PR)
}



#############################################

#-------------------------------------------
# Initial Distribution at t = 0
#-------------------------------------------


#' initial probability function
#'
#'
#' @param R The number of individuals marked and released at sampling occasion g from batch group g; g = 1,2,...,G. This must be an integer.
#' @return A vector of initial value with 1 at the observed position
#' @rdname delta_g


delta_g <- function(R){
  de <- rep(0, R+1)
  de[R+1] <- 1
  return(de)
}

##################################################


#------------------------------------------
# Batch Marking Log-likelihood
#-----------------------------------------

#'  batchLL function provides the batch marking log-likelihood
#'
#' @param phi The probability of surviving and remaining in the population between occasions t and t +1, given an individual was alive and in the population at occasion t. This must be a number between 0 and 1.
#' @param p The probability of capture at occasion t. This must be a number between 0 and 1.
#' @param R The number of individuals marked and released at sampling occasion g from batch group g; g = 1,2,...,G. This must be an integer.
#' @param begin_g The beginning of the occasion.
#' @param end_g The end of the occasion.
#' @return fr returns the log sum of the Hidden Markov Model.
#' @rdname batchLL
#' @importFrom parallel mclapply


batchLL <- function(phi, p, R, begin_g, end_g){


  U <- R[1]; u <- R[-1]

  #--------------------------------------
  # Rcpp_Armadillo code embedded in R

  f1 <- function(x){gamma_gt(U, x)}
  gm <- parallel::mclapply(phi, f1)

  #---------------------------------------

  #--------------------------------------
  # Compute when for marked group 1

  Gamma1 <- gm[[1]]
  delta0 <- delta_g(U)
  pu1    <- probs(u[1], p[1], U)
  a1     <- delta0%*%Gamma1%*%pu1

  #--------------------------------------

  #--------------------------------------

  # compute for t > 1

  if(begin_g < end_g){

    for (t in 2:length(u)) {

      Gammat <- gm[[t]]
      put <- probs(u[t], p[t], U)
      a1 <- a1%*%Gammat%*%put
    }
  }

  #----------------------------------------

  # Log-likelihood
  fr <- log(sum(a1))
  return(fr)
}

######################################################################


#---------------------------------------------------
# Batch Marking log-likelihood for different models
#---------------------------------------------------

#' @title Log-likelihood function for marked model.
#' @description
#' This helps users check whether the function can be optimized at the given initial values before optimizing using \code{\link{batchMarkOptim}}. After a quick check, if \code{NAN} or \code{Inf} is returned, the initial values should be revisited.
#'
#' @param par Initial values for the parameters to be optimized over.
#' @param data A capture-recapture data matrix or data frame.
#' @param covariate_phi This covariate placeholder for the parameter phi_t
#' @param covariate_p This covariate placeholder for the parameter p_t
#' @param choiceModel This chooses among different models and allows for model selection
#' @return Negative Log-likelihood value of the likelihood function
#' @rdname batchMarkHmmLL
#' @export
#'
#' @importFrom parallel detectCores makePSOCKcluster
#' @importFrom doParallel registerDoParallel
#' @import foreach foreach
#'
#' @examples
#' library(extBatchMarking)
#' # Initial parameter
#' theta <- c(0, -1)
#' res1  <- batchMarkHmmLL(par          = theta,
#'                        data          = WeatherLoach,
#'                        choiceModel   = "model4",
#'                        covariate_phi = NULL,
#'                        covariate_p   = NULL)
#' res1
#'


batchMarkHmmLL <- function(par = NULL, data, covariate_phi = NULL, covariate_p = NULL, choiceModel = c("model1", "model2", "model3", "model4")){

  if(!is.matrix(data)) data <- as.matrix(unname(data))
  if(nrow(data) != ncol(data)) stop(paste0("nrow must be equal to the ncol of the data"))

  # Number of capture and recapture
  R <- data[2:nrow(data), 1:ncol(data)]
  Ringed <- R[,1]; Mobs <- R[,-1]

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


  if(choiceModel == "model1"){

    # (phi_t, p_t)

    if(is.null(covariate_phi) & is.null(covariate_p)){

      Xphi <- diag(nrow(R))
      Xp   <- diag(nrow(R))

      # Parameters
      betaphi <- par[1:ncol(Xphi)]
      betap   <- par[(ncol(Xphi)+1):(ncol(Xp)+ncol(Xphi))]
      phi     <- batchLogit(t(Xphi) %*% betaphi)
      p       <- batchLogit(t(Xp) %*% betap)

    }else if(is.null(covariate_phi) & !is.null(covariate_p)){

      Xphi      <- diag(nrow(R))
      betaphi   <- par[1:ncol(Xphi)]
      phi       <- batchLogit(t(Xphi) %*% betaphi)
      covariate <- scale(covariate_p)
      Xp        <- matrix(c(rep(1, nrow(R)), covariate), nrow = nrow(R))
      betap     <- par[(ncol(Xphi)+1):(nrow(Xp)+ncol(Xphi))]
      p         <- batchLogit(as.matrix(rowSums(Xp * betap)))

    }else if(!is.null(covariate_phi) & is.null(covariate_p)){

      Xp        <- diag(nrow(R))
      covariate <- scale(covariate_phi)
      Xphi      <- matrix(c(rep(1, nrow(R)), covariate), nrow = nrow(R))
      betap     <- par[(ncol(Xphi)+1):(ncol(Xp)+ncol(Xphi))]
      p         <- batchLogit(t(Xp) %*% betap)
      betaphi   <- par[(ncol(Xp)+1):(nrow(Xphi)+ncol(Xp))]
      phi       <- batchLogit(as.matrix(rowSums(Xphi * betaphi)))

    }else{

      covariate1 <- scale(covariate_p)
      covariate2 <- scale(covariate_phi)
      Xphi      <- matrix(c(rep(1, nrow(R)), covariate2), nrow = nrow(R))
      Xp        <- matrix(c(rep(1, nrow(R)), covariate1), nrow = nrow(R))
      betap     <- par[(ncol(Xphi)+1):(nrow(Xp)+ncol(Xphi))]
      p         <- batchLogit(as.matrix(rowSums(Xp * betap)))
      betaphi   <- par[(ncol(Xp)+1):(nrow(Xphi)+ncol(Xp))]
      phi       <- batchLogit(as.matrix(rowSums(Xphi * betaphi)))

    }

  }

  else if(choiceModel == "model2"){

    # (phi_t, p_1) is model 2
    # Change the dimension of the Parameters

    Xp     <- matrix(1, ncol = 1, nrow = nrow(R))
    betap  <- par[ncol(Xp)]
    p      <- batchLogit(Xp %*% t(betap))

    if(is.null(covariate_phi) & is.null(covariate_p)){

      Xphi    <- diag(nrow(R))
      betaphi <- par[(ncol(Xp)+1):(ncol(Xphi)+ncol(Xp))]
      phi     <- batchLogit(t(Xphi) %*% betaphi)

    }else if(is.null(covariate_phi) & !is.null(covariate_p)){
      stop(paste0("Only supply covariate_phi for model 2. Hint: Set covariate_p = NULL"))
    }else if(!is.null(covariate_phi) & !is.null(covariate_p)){
      stop(paste0("You cannot supply covariate_p for 1-dimensional p in model 2."))
    }else{

      covariate <- scale(covariate_phi)
      Xphi      <- matrix(c(rep(1, nrow(R)), covariate), nrow = nrow(R))
      betaphi   <- par[(ncol(Xp)+1):(nrow(Xphi)+ncol(Xp))]
      phi       <- batchLogit(as.matrix(rowSums(Xphi * betaphi)))

    }

  }

  else if(choiceModel == 'model3'){

    # (phi_1, p_t) is model 3
    # Parameters

    Xphi        <- matrix(1, ncol = 1, nrow = nrow(R))
    betaphi     <- par[ncol(Xphi)]
    phi         <- batchLogit(Xphi %*% t(betaphi))

    if(is.null(covariate_phi) & is.null(covariate_p)) {
      Xp        <- diag(nrow(R))
      betap     <- par[(ncol(Xphi)+1):(ncol(Xp)+ncol(Xphi))]
      p         <- batchLogit(t(Xp) %*% betap)
    }else if(!is.null(covariate_phi) & is.null(covariate_p)){
      stop(paste0("Only supply covariate_p for model 3. Hint: Set covariate_phi = NULL"))
    }else if(!is.null(covariate_phi) & !is.null(covariate_p)){
      stop(paste0("You cannot supply covariate_phi for 1-dimensional phi in model 3."))
    }else{
      covariate <- scale(covariate_p)
      Xp        <- matrix(c(rep(1, nrow(R)), covariate), nrow = nrow(R))
      betap     <- par[(ncol(Xphi)+1):(nrow(Xp)+ncol(Xphi))]
      p         <- batchLogit(as.matrix(rowSums(Xp * betap)))
    }

  }else{

    if(!is.null(covariate_phi) & is.null(covariate_p)){
      stop(paste0("You cannot supply covariate_phi for 1-dimwnsional phi in model 4. Hint: Set covariate_phi = NULL"))
    }else if(is.null(covariate_phi) & !is.null(covariate_p)){
      stop(paste0("You cannot supply covariate_p for 1-dimwnsional p in model 4. Hint: Set covariate_p = NULL"))
    }else if(!is.null(covariate_phi) & !is.null(covariate_p)){
      stop(paste0("You cannot supply covariate_phi and covariate_p for 1-dimensional phi and p repspectively in model 4.
                  Hint: Set covariate_phi = NULL and covariate_p = NULL"))
    }else{

      # (phi_1, p_1) is model 4
      Xphi <- matrix(1, ncol = 1, nrow = nrow(R))
      Xp   <- matrix(1, ncol = 1, nrow = nrow(R))

      # Parameters
      betaphi <- par[ncol(Xphi)]
      betap   <- par[ncol(Xphi)+ncol(Xp)]
      phi     <- batchLogit(Xphi %*% t(betaphi))
      p       <- batchLogit(Xp %*% t(betap))

    }

  }

  logL <- 0

  for (i in 1:nrow(R)) {

    qw <- batchLL(phi[i:nrow(R)], p[i:nrow(R)], c(Ringed[i], Mobs[i, i:nrow(R)]), i, nrow(R))
    logL <- logL + qw

  }

  return(-logL)
}



#------------------------------------------
# Convolution of Poisson and Binomial
#-----------------------------------------

#'  Convolution of Poisson and Binomial for Batch
#'
#' This is the convolution of Poisson and Binomial distributions
#'
#' The convolution of Poisson and Binomial distribution helps us to compute the number of individuals that have survived from t-1 to t in the combined model while simultaneously computing the number of individuals recruited into the population at occasion t.
#'
#' The survival is modeled as Binomial distribution and the recruitment as the Poisson distirubiton
#' @param z This is the vector of numerical values
#' @param n The `nrow` of capture-recapture data matrix or data frame
#' @param par This is the vector of parameter values: average from Poisson distribution and probability of success from Binomial distribution
#' @return f This is the output of the convolution from the Binomial and Poisson distributions
#' @rdname dbinpois



dbinpois <- function(z, n, par){

  lambda <- par[1]; p <- par[2]
  zmax <- ceiling(max(z))
  x <- 0:zmax
  y <- 0:n

  nky <- lgamma(n+1) - lgamma(y+1) -lgamma(n-y+1)
  lnpfy <- nky + y*log(p) + (n-y)*log1p(-p)
  pz <- stats::convolve(dpois(x, lambda), rev(exp(lnpfy)), type = 'open')
  pz = pz[1:(zmax+1)]

  # Selected probabilities
  zvals <- 0:zmax
  iz <- z %in% zvals
  f <- iz - 0
  f[iz] <- pz[z[iz] + 1]

  return(f)

}

#---------------------------------------------------
# Unmarked model log-likelihood for different models
#---------------------------------------------------


#'  batchUnmarkHmmLL function provides the unmarked function to be optimized
#'
#' @param phi The probability of surviving and remaining in the population between occasions t and t +1, given an individual was alive and in the population at occasion t. This must be a number between 0 and 1.
#' @param p The probability of capture at occasion t. This must be a number between 0 and 1.
#' @param lambda The initial mean abundance (at occasion 1) for the unmarked population.
#' @param gam The recruitment rate into the unmarked population.
#' @param Umax The maximum number of the unmarked individuals in the population for capture on any occasion.
#' @param nBins The number of bin size into which the matrix will be divided.
#' @param u The number of individuals captured at sampling occasion t that were not marked; t = 1,...,T.
#' @return Negative Log-likelihood value of the likelihood function
#' @rdname batchUnmarkHmmLL
#' @importFrom optimbase ones zeros
#' @importFrom parallel detectCores makePSOCKcluster
#' @importFrom doParallel registerDoParallel
#' @import foreach foreach
#'
#'

batchUnmarkHmmLL <- function(phi, p, lambda, gam, Umax, nBins, u){

  # Initialization
  Nmax <- max(u) + 0.5 * (max(u) - min(u)); Nmax <- max(Umax, Nmax)
  gri  <- seq(from = 0, to = Nmax, by = nBins)
  m    <- length(gri) - 1
  mid  <- (gri[1:m] + gri[2:(m+1)]) / 2
  T1   <- length(u)
  gam  <- gam %*% optimbase::ones(1,T1)

  # initial distribution (at t = 1)
  delta <- ppois(gri[2:(m+1)] - 1, lambda) - ppois(gri[1:m] - 1, lambda)

  # first forward probability vector
  lnpfy1 <- rep(-Inf, m)
  indx1  <- mid > u[1]
  nky1   <- lgamma(mid[indx1]+1)-lgamma(u[1]+1)-lgamma(mid[indx1]-u[1]+1)
  lnpfy1[indx1] <- nky1 + u[1] * log(p[1]) + (mid[indx1]-u[1]) *log1p(-p[1])

  Py1 <- diag(exp(lnpfy1))
  a1  <- delta %*% Py1

  Gammat1 <- optimbase::zeros(m)

  for(t in 2:T1){

    for(indx in 1:m){

      thetconv <- c(gam[t-1] %*% mid[indx], phi[t-1])
      probs <- dbinpois(0:(gri[m+1] - 1), mid[indx], thetconv)
      Gammat1[indx,] <- colSums(matrix(probs, nrow = nBins, byrow = FALSE))

    }

  }

  for (t in 1:(T1-1)){

    lnpfyt <- rep(-Inf, m)
    indxt  <- mid > u[t+1]
    nkyt   <- lgamma(mid[indxt] + 1) - lgamma(u[t+1] + 1) - lgamma(mid[indxt] - u[t+1] + 1)
    lnpfyt[indxt] <- nkyt + u[t+1] * log(p[t+1]) + (mid[indxt]-u[t+1]) * log1p(-p[t+1])
    pyt    <- diag(exp(lnpfyt))
    a1     <- a1 %*% Gammat1 %*% pyt

  }

  # Remove any negative numbers
  a1 <- a1[a1 > 0]

  #parallel::stopCluster(cl)
  return(log(sum(a1, na.rm = TRUE)))

}

#---------------------------------------------------
# Viterbi algorithm for the unmarked model
#---------------------------------------------------

#'  batchUnmarkViterbi function provides the implementation of the Viterbi alogrithm for the unmarked model
#'
#' @param phi The probability of surviving and remaining in the population between occasions t and t +1, given an individual was alive and in the population at occasion t. This must be a number between 0 and 1.
#' @param p The probability of capture at occasion t. This must be a number between 0 and 1.
#' @param lambda the initial mean abundance (at occasion 1) for the unmarked population.
#' @param gam The recruitment rate into the unmarked population
#' @param Umax The maximum number of the unmarked individuals in the population for capture on any occasion.
#' @param nBins The number of bin size into which the matrix will be divided.
#' @param u The number of individuals captured at sampling occasion t that were not marked; t = 1,...,T.
#' @return Negative Log-likelihood value of the likelihood function
#' @rdname batchUnmarkViterbi
#' @importFrom optimbase ones zeros
#' @importFrom parallel detectCores makePSOCKcluster
#' @importFrom doParallel registerDoParallel
#' @import foreach foreach



batchUnmarkViterbi <- function(phi, p, lambda, gam, Umax, nBins, u){

  Nmax <- max(u) + 0.5 * (max(u) - min(u)); Nmax <- max(Umax, Nmax)
  gri  <- seq(from = 0, to = Nmax, by = nBins)
  m    <- length(gri) - 1
  mid  <- (gri[1:m] + gri[2:(m+1)]) / 2
  T1   <- length(u)
  gam  <- gam %*% optimbase::ones(1,T1)

  # initial distribution (at t = 1)
  delta <- ppois(gri[2:(m+1)] - 1, lambda) - ppois(gri[1:m] - 1, lambda)


  # first forward probability vector
  lnpfy1 <- rep(-Inf, m)
  indx1  <- mid > u[1]
  nky1   <- lgamma(mid[indx1]+1)-lgamma(u[1]+1)-lgamma(mid[indx1]-u[1]+1)
  lnpfy1[indx1] <- nky1 + u[1] * log(p[1]) + (mid[indx1]-u[1]) *log1p(-p[1])

  Py1 <- diag(exp(lnpfy1)); zeta <- optimbase::zeros(T1, m)
  zeta[1,] <- delta %*% Py1

  Gammat1 <- optimbase::zeros(m)

  for(t in 2:T1){

    for(indx in 1:m){

      thetconv <- c(gam[t-1] %*% mid[indx], phi[t-1])
      probs <- dbinpois(0:(gri[m+1] - 1), mid[indx], thetconv)
      Gammat1[indx,] <- colSums(matrix(probs, nrow = nBins, byrow = FALSE))

    }

  }

  for (t in 1:(T1-1)){

    lnpfyt <- rep(-Inf, m)
    indxt  <- mid > u[t+1]
    nkyt   <- lgamma(mid[indxt] + 1) - lgamma(u[t+1] + 1) - lgamma(mid[indxt] - u[t+1] + 1)
    lnpfyt[indxt] <- nkyt + u[t+1] * log(p[t+1]) + (mid[indxt]-u[t+1]) * log1p(-p[t+1])
    pyt    <- diag(exp(lnpfyt))
    zeta[t+1,] <- apply((Gammat1 * zeta[t,]), 2, max) %*% pyt
  }

  # Compute maximizing sequence of states
  imax <- which.max(zeta[T1,]); midhat <- NULL; midhat[T1] <- mid[imax]

  for (t in (T1-1):1) {

    imax <- which.max(zeta[t,] * Gammat1[,imax])
    midhat[t] <- mid[imax]

  }

  return(midhat)

}

#---------------------------------------------------
# Marked+Unmarked models for four different models
#---------------------------------------------------

#' @title Log-likelihood function for combined model.
#' @description
#' This helps users check whether the function can be optimized at the given initial values before optimizing using \code{\link{batchMarkUnmarkOptim}}. After a quick check, if \code{NAN} or \code{Inf} is returned, the initial values should be revisited.
#'
#' @param par Initial values for the parameters to be optimized over.
#' @param data A capture-recapture data matrix or data frame.
#' @param Umax The maximum number of the unmarked individuals in the population for capture on any occasion.
#' @param covariate_phi This covariate placeholder for the parameter phi_t
#' @param covariate_p This covariate placeholder for the parameter p_t
#' @param nBins The number of bin size into which the matrix will be divided.
#' @param choiceModel This chooses among different models and allow for model selection.
#' @return Negative Log-likelihood value of the likelihood function.
#' @rdname batchMarkUnmarkHmmLL
#' @export
#'
#' @importFrom optimbase ones zeros
#' @importFrom parallel detectCores makePSOCKcluster
#' @importFrom doParallel registerDoParallel
#' @import foreach foreach

#'
#' @examples
#' library(extBatchMarking)
#' theta <- c(0.1, 0.1, 7, -1.5)
#' res3  <- batchMarkUnmarkHmmLL(par        = theta,
#'                              data        = WeatherLoach,
#'                              choiceModel = "model4",
#'                              Umax        = 1800,
#'                              nBins       = 600,
#'                        covariate_phi = NULL,
#'                        covariate_p   = NULL)
#' res3


batchMarkUnmarkHmmLL <- function(par, data, Umax, nBins, covariate_phi = NULL, covariate_p = NULL,
                                 choiceModel = c("model1", "model2", "model3", "model4")){

  if(!is.matrix(data)) data <- as.matrix(unname(data))
  if(nrow(data) != ncol(data)) stop(paste0("nrow must be equal to the ncol of the data"))

  # Number of capture and recapture
  R <- data[2:nrow(data), 1:ncol(data)]
  Ringed <- R[,1]; Mobs <- R[,-1]
  u <- data[1,] - c(0, colSums(Mobs))

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

  if(choiceModel == "model1"){

    # (phi_t, p_t) model 1

    if(is.null(covariate_phi) & is.null(covariate_p)){

      Xphi <- diag(nrow(R))
      Xp   <- diag(nrow(R) + 1)

      # Parameters
      betaphi <- par[1:ncol(Xphi)]
      betap   <- par[(ncol(Xphi)+1):(ncol(Xp)+ncol(Xphi))]
      phi     <- batchLogit(t(Xphi) %*% betaphi)
      p       <- batchLogit(t(Xp) %*% betap)

      loglambda <- par[ncol(Xphi)+ncol(Xp)+1]; lambda <- exp(loglambda)
      loggamma  <- par[ncol(Xphi)+ncol(Xp)+2]; gam    <- exp(loggamma)

    }else if(is.null(covariate_phi) & !is.null(covariate_p)){

      Xphi      <- diag(nrow(R))
      betaphi   <- par[1:ncol(Xphi)]
      phi       <- batchLogit(t(Xphi) %*% betaphi)
      covariate <- scale(covariate_p)
      if(nrow(covariate) == nrow(R)) stop("nrow of covariate for parameter 'p' must be : ", nrow(R)+1, " and not ", nrow(covariate))
      Xp        <- matrix(c(rep(1, (nrow(R)+1)), covariate), nrow = (nrow(R)+1))
      betap     <- par[(ncol(Xphi)+1):(nrow(Xp)+ncol(Xphi))]
      p         <- batchLogit(as.matrix(rowSums(Xp * betap)))

      loglambda <- par[ncol(Xphi)+nrow(Xp)+1]; lambda <- exp(loglambda)
      loggamma  <- par[ncol(Xphi)+nrow(Xp)+2]; gam    <- exp(loggamma)

    }else if(!is.null(covariate_phi) & is.null(covariate_p)){

      Xp        <- diag(nrow(R)+1)
      covariate <- scale(covariate_phi)
      if(nrow(covariate) != nrow(R)) stop("covariate for parameter 'phi' must be : ", nrow(R), " and not ", nrow(covariate))
      Xphi      <- matrix(c(rep(1, nrow(R)), covariate), nrow = nrow(R))
      betap     <- par[(nrow(Xphi)+1):(ncol(Xp)+nrow(Xphi))]
      p         <- batchLogit(t(Xp) %*% betap)
      betaphi   <- betaphi <- par[1:nrow(Xphi)]
      phi       <- batchLogit(as.matrix(rowSums(Xphi * betaphi)))

      loglambda <- par[nrow(Xphi)+ncol(Xp)+1]; lambda <- exp(loglambda)
      loggamma  <- par[nrow(Xphi)+ncol(Xp)+2]; gam    <- exp(loggamma)

    }else{

      covariate1 <- scale(covariate_p)
      covariate2 <- scale(covariate_phi)
      if(nrow(covariate2) != nrow(R) | nrow(covariate1) == nrow(R)) stop("nrow of covariate for parameter 'phi' must be : ", nrow(R), " and not ", nrow(covariate2),
                                                                         "; nrow of covariate for parameter 'p' must be : ", nrow(R)+1, " and not ", nrow(covariate1))
      Xphi      <- matrix(c(rep(1, nrow(R)), covariate2), nrow = nrow(R))
      Xp        <- matrix(c(rep(1, nrow(R)+1), covariate1), nrow = nrow(R)+1)
      betap     <- par[(nrow(Xphi)+1):(nrow(Xp)+nrow(Xphi))]
      p         <- batchLogit(as.matrix(rowSums(Xp * betap)))
      betaphi   <- par[1:nrow(Xphi)]
      phi       <- batchLogit(as.matrix(rowSums(Xphi * betaphi)))

      loglambda <- par[nrow(Xphi)+nrow(Xp)+1]; lambda <- exp(loglambda)
      loggamma  <- par[nrow(Xphi)+nrow(Xp)+2]; gam    <- exp(loggamma)

    }

  }

  else if(choiceModel == 'model2'){

    # (phi_t, p_1) is model 2

    if(is.null(covariate_phi) & is.null(covariate_p)){

      Xphi <- diag(nrow(R))
      Xp   <- matrix(1, ncol = 1, nrow = nrow(R) + 1)

      # Change the dimension of the Parameters
      betap   <- par[ncol(Xp)]
      betaphi <- par[(ncol(Xp)+1):(ncol(Xphi)+ncol(Xp))]
      p       <- batchLogit(Xp %*% t(betap))
      phi     <- batchLogit(t(Xphi) %*% betaphi)

      loglambda <- par[ncol(Xphi)+ncol(Xp)+1]; lambda <- exp(loglambda)
      loggamma  <- par[ncol(Xphi)+ncol(Xp)+2]; gam    <- exp(loggamma)

    }else if(is.null(covariate_phi) & !is.null(covariate_p)){
      stop(paste0("Only supply covariate_phi for model 2. Hint: Set covariate_p = NULL"))
    }else if(!is.null(covariate_phi) & !is.null(covariate_p)){
      stop(paste0("You cannot supply covariate_p for 1-dimensional p in model 2."))
    }else{

      covariate <- scale(covariate_phi)
      if(nrow(covariate) != nrow(R)) stop("covariate for parameter 'phi' must be : ", nrow(R), " and not ", nrow(covariate))
      Xp        <- matrix(1, ncol = 1, nrow = nrow(R) + 1)
      Xphi      <- matrix(c(rep(1, nrow(R)), covariate), nrow = nrow(R))
      betap     <- par[ncol(Xp)]
      betaphi   <- par[(ncol(Xp)+1):(nrow(Xphi)+ncol(Xp))]
      p         <- batchLogit(Xp %*% t(betap))
      phi       <- batchLogit(as.matrix(rowSums(Xphi * betaphi)))

      loglambda <- par[nrow(Xphi)+ncol(Xp)+1]; lambda <- exp(loglambda)
      loggamma  <- par[nrow(Xphi)+ncol(Xp)+2]; gam    <- exp(loggamma)

    }

  }

  else if(choiceModel == "model3"){

    # (phi_1, p_t) is model 3

    if(is.null(covariate_phi) & is.null(covariate_p)){

      Xphi <- matrix(1, ncol = 1, nrow = nrow(R))
      Xp   <- diag(nrow(R) + 1)

      # Parameters
      betaphi <- par[ncol(Xphi)]
      betap   <- par[(ncol(Xphi)+1):(ncol(Xp)+ncol(Xphi))]
      phi     <- batchLogit(Xphi %*% t(betaphi))
      p       <- batchLogit(t(Xp) %*% betap)

      loglambda <- par[ncol(Xphi)+ncol(Xp)+1]; lambda <- exp(loglambda)
      loggamma  <- par[ncol(Xphi)+ncol(Xp)+2]; gam    <- exp(loggamma)

    }else if(!is.null(covariate_phi) & is.null(covariate_p)){
      stop(paste0("Only supply covariate_p for model 3. Hint: Set covariate_phi = NULL"))
    }else if(!is.null(covariate_phi) & !is.null(covariate_p)){
      stop(paste0("You cannot supply covariate_phi for 1-dimensional phi in model 3."))
    }else{

      Xphi    <- matrix(1, ncol = 1, nrow = nrow(R))
      betaphi <- par[ncol(Xphi)]
      covariate <- scale(covariate_p)
      if(nrow(covariate) == nrow(R)) stop("nrow of covariate for parameter 'p' must be : ", nrow(R)+1, " and not ", nrow(covariate))
      Xp      <- matrix(c(rep(1, (nrow(R)+1)), covariate), nrow = (nrow(R)+1))
      betap   <- par[(ncol(Xphi)+1):(nrow(Xp)+ncol(Xphi))]
      phi     <- batchLogit(Xphi %*% t(betaphi))
      p       <- batchLogit(as.matrix(rowSums(Xp * betap)))

      loglambda <- par[ncol(Xphi)+nrow(Xp)+1]; lambda <- exp(loglambda)
      loggamma  <- par[ncol(Xphi)+nrow(Xp)+2]; gam    <- exp(loggamma)

    }


  }else{

    # (phi_1, p_1) is model 4

    if(!is.null(covariate_phi) & is.null(covariate_p)){
      stop(paste0("You cannot supply covariate_phi for 1-dimwnsional phi in model 4. Hint: Set covariate_phi = NULL"))
    }else if(is.null(covariate_phi) & !is.null(covariate_p)){
      stop(paste0("You cannot supply covariate_p for 1-dimwnsional p in model 4. Hint: Set covariate_p = NULL"))
    }else if(!is.null(covariate_phi) & !is.null(covariate_p)){
      stop(paste0("You cannot supply covariate_phi and covariate_p for 1-dimensional phi and p repspectively in model 4.
                  Hint: Set covariate_phi = NULL and covariate_p = NULL"))
    }else{

      Xphi <- matrix(1, ncol = 1, nrow = nrow(R))
      Xp   <- matrix(1, ncol = 1, nrow = nrow(R) + 1)

      # Parameters
      betaphi <- par[ncol(Xphi)]
      betap   <- par[ncol(Xphi)+ncol(Xp)]
      phi     <- batchLogit(Xphi %*% t(betaphi))
      p       <- batchLogit(Xp %*% t(betap))

      loglambda <- par[ncol(Xphi)+ncol(Xp)+1]; lambda <- exp(loglambda)
      loggamma  <- par[ncol(Xphi)+ncol(Xp)+2]; gam    <- exp(loggamma)
    }

  }


  pm <- p[2:length(p)]; pu <- p[1:length(p)]

  # log-likelihood for marked
  logLm <- 0

  # Log-likelihood for the Marked
  for (i in 1:nrow(R)) {

    qw    <- batchLL(phi[i:nrow(R)], pm[i:nrow(R)], c(Ringed[i], Mobs[i, i:nrow(R)]), i, nrow(R))
    logLm <- logLm + qw

  }

  # log-likelihood for unmarked
  logLu <- batchUnmarkHmmLL(phi, pu, lambda, gam, Umax, nBins, u)

  f1 <- -logLm - logLu

  return(f1)

}


#---------------------------------------------------
# Viterbi Algorithm for four different models
#---------------------------------------------------


#'  batchUnmark2Viterbi function provides a wrapper for the batchUnmarkViterbi to compute the popuation abundance
#'
#' @param par Initial values for the parameters to be optimized over.
#' @param data A capture-recapture data matrix or data frame.
#' @param Umax The maximum number of the unmarked individuals in the population for capture on any occasion.
#' @param nBins The number of bin size into which the matrix will be divided.
#' @param choiceModel This chooses among different models and allows for model selection
#' @return Negative Log-likelihood value of the likelihood function
#' @rdname batchUnmark2Viterbi



batchUnmark2Viterbi <- function(par, data, Umax, nBins, choiceModel = c("model1", "model2", "model3", "model4")){

  if(!is.matrix(data)) data <- as.matrix(unname(data))

  # Number of capture and recapture
  R <- data[2:nrow(data), 1:ncol(data)]
  Ringed <- R[,1]; Mobs <- R[,-1]
  u <- data[1,] - c(0, colSums(Mobs))

  if(choiceModel == "model1"){

    # (phi_t, p_t) model 1

    Xphi <- diag(nrow(R))
    Xp   <- diag(nrow(R) + 1)

    # Parameters
    betaphi <- par[1:ncol(Xphi)]
    betap   <- par[(ncol(Xphi)+1):(ncol(Xp)+ncol(Xphi))]
    phi     <- batchLogit(t(Xphi) %*% betaphi)
    p       <- batchLogit(t(Xp) %*% betap)

  }

  else if(choiceModel == 'model2'){

    # (phi_t, p_1) is model 2

    Xphi <- diag(nrow(R))
    Xp   <- matrix(1, ncol = 1, nrow = nrow(R) + 1)

    # Change the dimension of the Parameters
    betap   <- par[ncol(Xp)]
    betaphi <- par[(ncol(Xp)+1):(ncol(Xphi)+ncol(Xp))]
    p       <- batchLogit(Xp %*% t(betap))
    phi     <- batchLogit(t(Xphi) %*% betaphi)

  }

  else if(choiceModel == "model3"){

    # (phi_1, p_t) is model 3
    Xphi <- matrix(1, ncol = 1, nrow = nrow(R))
    Xp   <- diag(nrow(R) + 1)

    # Parameters
    betaphi <- par[ncol(Xphi)]
    betap   <- par[(ncol(Xphi)+1):(ncol(Xp)+ncol(Xphi))]
    phi     <- batchLogit(Xphi %*% t(betaphi))
    p       <- batchLogit(t(Xp) %*% betap)

  }else{

    # (phi_1, p_1) is model 4
    Xphi <- matrix(1, ncol = 1, nrow = nrow(R))
    Xp   <- matrix(1, ncol = 1, nrow = nrow(R) + 1)

    # Parameters
    betaphi <- par[ncol(Xphi)]
    betap   <- par[ncol(Xphi)+ncol(Xp)]
    phi     <- batchLogit(Xphi %*% t(betaphi))
    p       <- batchLogit(Xp %*% t(betap))

  }

  loglambda <- par[ncol(Xphi)+ncol(Xp)+1]; lambda <- exp(loglambda)
  loggamma  <- par[ncol(Xphi)+ncol(Xp)+2]; gam    <- exp(loggamma)

  pm <- p[2:length(p)]
  pu <- p[1:length(p)]
  U  <- batchUnmarkViterbi(phi, p, lambda, gam, Umax, nBins, u)

  return(U)

}

