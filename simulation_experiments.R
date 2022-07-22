################################################################################
# simulation experiments
#
# date: 2022/07/13
################################################################################

# load required packages and helper functions
library(parallel)
source("helper_functions.R")

# setting
dgp <- 1
R <- 10000
psi_vector <- seq(0.0, 0.9, by = 0.1)
ncores <- 15
date <- paste0(Sys.Date())
file <- paste0(date, "_simulation_experiments_results.rds")

# simulation experiments
reslist <- lapply(psi_vector, function(psi, dgp, R, ncores) {
  cat(
    "--- ", paste0(Sys.time()),
    " - Fixed N= 200, T= 5",
    " - share uninformative data (psi)= ", psi,
    " ---\n",
    sep = ""
  )
  reslist <- mclapply(1:R, function(r, psi, dgp) {
    fm <- y ~ x | exp_time + imp_time + pair
    dt_sim <- simulateData(200, 5, psi, dgp)
    mod <- feglm(fm, dt_sim, poisson())
    b <- coef(mod)
    b_abc <- biasCorrPPML(mod, adj_df = TRUE)
    bias_spj <- biasCorrSPJ(mod, 1)
    b_spj <- 2 * b - bias_spj[1]
    data.frame(
      psi      = psi,
      beta     = b,
      beta_abc = b_abc,
      beta_spj = b_spj
    )
  }, psi = psi, dgp = dgp, mc.cores = ncores)
  do.call(rbind, reslist)
}, dgp = dgp, R = R, ncores = ncores)

# store results
res <- do.call(rbind, reslist)
saveRDS(res, file = file, compress = "xz")
