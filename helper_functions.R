################################################################################
# helper functions
#
# 1) analytical bias-correction of WZ 2021
# 2) split-panel jackknife bias-correction of WZ 2021
# 3) drop perfectly predicted/separated observations
# 4) data generating process
#
# date: 2022/07/13
################################################################################

# required packages
library(alpaca)
library(data.table)
library(forcats)
library(lfe)
library(MASS)

# analytical bias-correction of WZ 2021
biasCorrPPML <- function(obj, adj_df = TRUE) {
  # quantities from uncorrected ppml estimator
  dt <- copy(obj$data)
  fam <- poisson()
  W <- obj$Hessian
  beta_uncorr <- coef(obj)
  dt[, eta := predict(obj)]
  dt <- dt[, lambda := fam$linkinv(eta)]
  dt[, xtilde := as.vector(residuals(
    felm(x ~ 0 | exp_time + imp_time + pair, dt, weights = dt$lambda)
  ))]

  # prepare for unbalanced panels
  dt[, c("exp", "time") := tstrsplit(as.character(exp_time), "\\.")]
  dt[, c("imp", "time") := tstrsplit(as.character(imp_time), "\\.")]
  dt[, exp := as.integer(factor(exp))]
  dt[, imp := as.integer(factor(imp))]
  dt[, time := as.integer(factor(time))]
  dt <- dt[
    CJ(exp = exp, imp = imp, time = time, unique = TRUE),
    on = .(exp, imp, time)
  ]
  setnafill(dt, fill = 0, cols = c("y", "eta", "xtilde", "lambda"))

  # pair-wise quantities
  setkey(dt, imp, exp, time)
  nexp <- dt[, max(exp)]
  nimp <- dt[, max(imp)]
  ntime <- dt[, max(time)]
  y <- array(dt$y, c(ntime, nexp, nimp))
  sum_y <- apply(y, 2:3, sum)
  xtilde <- array(dt$xtilde, c(ntime, nexp, nimp))
  lambda <- array(dt$lambda, c(ntime, nexp, nimp))
  sum_lambda <- apply(lambda, 2:3, sum)
  theta <- array(NA_real_, c(ntime, nexp, nimp))
  S <- array(NA_real_, c(ntime, nexp, nimp))
  for (t in 1:ntime) {
    theta[t, , ] <- lambda[t, , ] / sum_lambda
    S[t, , ] <- y[t, , ] - theta[t, , ] * sum_y
  }
  theta[is.na(theta)] <- 0
  S[is.na(S)] <- 0
  H <- array(NA_real_, c(ntime, ntime, nexp, nimp))
  Hbar <- array(NA_real_, c(ntime, ntime, nexp, nimp))
  for (j in 1:nimp) {
    for (i in 1:nexp) {
      theta_ij <- theta[, i, j]
      H[, , i, j] <- (diag(theta_ij) - tcrossprod(theta_ij)) * sum_y[i, j]
      Hbar[, , i, j] <- (diag(theta_ij) - tcrossprod(theta_ij)) * sum_lambda[i, j]
    }
  }
  Gbar <- array(NA_real_, c(ntime, ntime, ntime, nexp, nimp))
  for (j in 1:nimp) {
    for (i in 1:nexp) {
      theta_ij <- theta[, i, j]
      sum_lambda_ij <- sum_lambda[i, j]
      for (t in 1:ntime) {
        for (s in 1:ntime) {
          for (r in 1:ntime) {
            # case 1: t=s=r
            if (t == s && s == r) {
              Gbar[t, s, r, i, j] <- -theta_ij[t] * (1 - theta_ij[t]) * (1 - 2 * theta_ij[t]) * sum_lambda_ij
            }
            # case 2: s=r!=t
            if (s == r && r != t) {
              Gbar[t, s, r, i, j] <- theta_ij[s] * (1 - 2 * theta_ij[s]) * theta_ij[t] * sum_lambda_ij
            }
            # case 3: t=s!=r
            if (t == s && s != r) {
              Gbar[t, s, r, i, j] <- theta_ij[s] * (1 - 2 * theta_ij[s]) * theta_ij[r] * sum_lambda_ij
            }
            # case 4: r=t!=s
            if (r == t && t != s) {
              Gbar[t, s, r, i, j] <- theta_ij[t] * (1 - 2 * theta_ij[t]) * theta_ij[s] * sum_lambda_ij
            }
            # case 5: r!=s!=t!=r
            if (r != s && s != t && t != r) {
              Gbar[t, s, r, i, j] <- -2 * theta_ij[r] * theta_ij[s] * theta_ij[t] * sum_lambda_ij
            }
          }
        }
      }
    }
  }
  SS <- array(NA_real_, c(ntime, ntime, nexp, nimp))
  HXS <- array(NA_real_, c(ntime, ntime, nexp, nimp))
  GbarX <- array(0, c(ntime, ntime, nexp, nimp))
  for (j in 1:nimp) {
    for (i in 1:nexp) {
      S_ij <- S[, i, j]
      xtilde_ij <- xtilde[, i, j]
      SS[, , i, j] <- tcrossprod(S_ij)
      HXS[, , i, j] <- tcrossprod(H[, , i, j] %*% xtilde_ij, S_ij)
      for (r in 1:ntime) {
        GbarX[, , i, j] <- GbarX[, , i, j] + Gbar[, , r, i, j] * xtilde_ij[r]
      }
    }
  }
  rm(y, xtilde, sum_y, sum_lambda, S, H, Gbar)

  # bias as sum over exporters
  B1 <- numeric(nexp)
  B2 <- numeric(nexp)
  for (i in 1:nexp) {
    Hbar_i <- matrix(0, ntime, ntime)
    SS_i <- matrix(0, ntime, ntime)
    HXS_i <- matrix(0, ntime, ntime)
    GbarX_i <- matrix(0, ntime, ntime)
    for (j in 1:nimp) {
      Hbar_i <- Hbar_i + Hbar[, , i, j]
      SS_i <- SS_i + SS[, , i, j]
      HXS_i <- HXS_i + HXS[, , i, j]
      GbarX_i <- GbarX_i + GbarX[, , i, j]
    }
    Hbar_i_inv <- ginv(Hbar_i)
    B1[i] <- sum(diag(Hbar_i_inv %*% HXS_i))
    B2[i] <- sum(diag(GbarX_i %*% Hbar_i_inv %*% SS_i %*% Hbar_i_inv)) / 2
  }
  bias_exp <- solve(W, sum(-B1 + B2))

  # bias as sum over importers
  D1 <- numeric(nimp)
  D2 <- numeric(nimp)
  for (j in 1:nimp) {
    Hbar_j <- matrix(0, ntime, ntime)
    SS_j <- matrix(0, ntime, ntime)
    HXS_j <- matrix(0, ntime, ntime)
    GbarX_j <- matrix(0, ntime, ntime)
    for (i in 1:nexp) {
      Hbar_j <- Hbar_j + Hbar[, , i, j]
      SS_j <- SS_j + SS[, , i, j]
      HXS_j <- HXS_j + HXS[, , i, j]
      GbarX_j <- GbarX_j + GbarX[, , i, j]
    }
    Hbar_j_inv <- ginv(Hbar_j)
    D1[j] <- sum(diag(Hbar_j_inv %*% HXS_j))
    D2[j] <- sum(diag(GbarX_j %*% Hbar_j_inv %*% SS_j %*% Hbar_j_inv)) / 2
  }
  bias_imp <- solve(W, sum(-D1 + D2))

  # adjust degrees of freedom
  if (adj_df) {
    bias_exp <- bias_exp * nexp / (nexp - 1)
    bias_imp <- bias_imp * nimp / (nimp - 1)
  }

  # debias \beta
  beta_uncorr - bias_exp - bias_imp
}

# split-panel bias correction
biasCorrSPJ <- function(obj, R) {
  fm <- obj$formula
  dt <- copy(obj$data)
  dt[, c("exp", "time") := tstrsplit(as.character(exp_time), "\\.")]
  dt[, c("imp", "time") := tstrsplit(as.character(imp_time), "\\.")]
  dt[, exp := factor(exp)]
  dt[, imp := factor(imp)]
  nexp <- dt[, length(unique(exp))]
  nimp <- dt[, length(unique(imp))]
  bias <- sapply(1:R, function(r) {
    dt[, exp_r := as.integer(fct_shuffle(exp))]
    dt[, imp_r := as.integer(fct_shuffle(imp))]
    dt[, exp_group := "a"]
    dt[exp_r > floor(nexp / 2), exp_group := "b"]
    dt[, imp_group := "a"]
    dt[imp_r > floor(nimp / 2), imp_group := "b"]
    mod1 <- feglm(fm, dt[exp_group == "a" & imp_group == "a"], poisson())
    mod2 <- feglm(fm, dt[exp_group == "a" & imp_group == "b"], poisson())
    mod3 <- feglm(fm, dt[exp_group == "b" & imp_group == "a"], poisson())
    mod4 <- feglm(fm, dt[exp_group == "b" & imp_group == "b"], poisson())
    (coef(mod1) + coef(mod2) + coef(mod3) + coef(mod4)) / 4
  })
  bias
}

# drop uninformative observations
dropObservations <- function(
    dt,
    depvar   = "y",
    fevars   = c("exp_time", "imp_time", "pair"),
    tempvar1 = "ybar",
    tempvar2 = "nk"
    ) {
  setDT(dt)
  n <- nrow(dt)
  for (iter in 1:100) {
    n0 <- n
    for (k in fevars) {
      dt[, (tempvar1) := mean(get(depvar)), by = k]
      dt <- dt[get(tempvar1) > 0]
    }
    for (k in fevars) {
      dt[, (tempvar2) := .N, by = k]
      dt <- dt[get(tempvar2) > 1]
    }
    n <- nrow(dt)
    if (n == n0) break
  }
  dt[, (tempvar1) := NULL]
  dt[, (tempvar2) := NULL]
  dt
}

# data generating process
# note 1: conditional missing at random, given \eta_{ij}
# note 2: WZ 2021 define N(0, 1) as mean = 0 and sd = 1
simulateData <- function(n, t, share_missing, dgp = 1) {
  exp <- rep(rep(1:n, each = n), t + 1)
  imp <- rep(rep(1:n, n), t + 1)
  time <- rep(0:t, each = n^2)
  dt <- data.table(exp, imp, time)
  dt[, exp_time := interaction(exp, time)]
  dt[, imp_time := interaction(imp, time)]
  dt[, pair := interaction(exp, imp)]
  setkey(dt, exp_time, imp_time, pair)
  dt[, alpha := rnorm(1, 0, 1 / 16), by = exp_time]
  dt[, gamma := rnorm(1, 0, 1 / 16), by = imp_time]
  dt[, eta := rnorm(1, 0, 1 / 16), by = pair]
  dt[, pi := alpha + gamma + eta]
  setkey(dt, time, pair)
  X <- matrix(NA_real_, n^2, t + 1)
  X[, 1] <- dt[time == 0, eta] + rnorm(n^2, 0, 0.5)
  Z <- matrix(0, n^2, t + 1)
  Z[, 1] <- rnorm(n^2)
  for (s in 2:(t + 1)) {
    pi <- dt[time == s - 1, pi]
    X[, s] <- X[, s - 1] / 2 + pi + rnorm(n^2, 0, 0.5)
    Z[, s] <- 0.3 * Z[, s - 1] + rnorm(n^2, 0, sqrt(1 - 0.3^2))
  }
  x <- as.vector(X[, -1])
  z <- as.vector(Z[, -1])
  dt <- dt[time != 0]
  dt[, x := x]
  dt[, z := z]
  # Here we set all trade flows of country pairs with the smallest
  # value of \eta_ij to zero. The corresponding observations become
  # uninformative and can be dropped.
  eta_thresh <- dt[time == 1, quantile(eta, share_missing, type = 1)]
  dt <- dt[eta > eta_thresh]
  dt[, lambda := exp(x + pi)]
  if (dgp == 1) {
    dt[, sigma2 := 1 / lambda^2]
  } else if (dgp == 2) {
    dt[, sigma2 := 1 / lambda]
  } else if (dgp == 3) {
    dt[, sigma2 := 1]
  } else {
    dt[, sigma2 := 0.5 / lambda + 0.5 * exp(2 * x)]
  }
  dt[, omega := exp(-0.5 * log(1 + sigma2) + sqrt(log(1 + sigma2)) * z)]
  dt[, y := lambda * omega]
  dt[, c("alpha", "gamma", "eta", "pi", "lambda", "sigma2", "omega", "z") := NULL]
  dt
}
