################################################################################
# simulation results
#
# date: 2022/07/13
################################################################################

## load required packages
library(data.table)
library(forcats)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(RColorBrewer)
library(xtable)

## results paper
# prepare results
# note: zzz.rds has to be adjusted
res <- readRDS("zzz.rds")
setDT(res)
setkey(res, psi)
measure_vars <- c("beta", paste0("beta_", c("abc", "spj")))
res <- melt(
  res,
  id.vars       = "psi",
  measure.vars  = measure_vars,
  variable.name = "estimator",
  value.name    = "beta"
)
new_levels <- c(
  "FEPPML" = "beta",
  "ABC"    = "beta_abc",
  "SPJ"    = "beta_spj"
)
res[, estimator := fct_recode(estimator, !!!new_levels)]
res2 <- res[, .(
  bias     = abs(mean(beta) - 1),
  sd       = sd(beta)
), by = .(estimator, psi)]
res2[, ratio := bias / sd]
res2[, bias := bias * 100]
id_vars <- c("psi", "estimator")
measure_vars <- c("bias", "ratio")
res3 <- melt(
  res2,
  id.vars       = id_vars,
  measure.vars  = measure_vars,
  variable.name = "statistic",
  value.name    = "value"
)
new_levels <- c(
  "Bias (in %)" = "bias",
  "Bias / SD"   = "ratio"
)
res3[, statistic := fct_recode(statistic, !!!new_levels)]

# figure
ggplot(res3, aes(x = factor(psi), y = value, color = estimator, group = estimator)) +
  facet_wrap(. ~ statistic, scales = "free_y") +
  geom_line() +
  geom_point() +
  labs(
    x     = "Share of Uninformative Observations",
    y     = NULL,
    color = "Estimator"
  ) +
  scale_color_manual(
    values = c("FE-PPML" = "#D9D9D9", "ABC" = "#969696", "SPJ" = "#525252")
  ) +
  theme_few() +
  theme(legend.position = "bottom")
ggsave("results_simulation.pdf", width = 8.27, height = 4.13)

# table supplement
res2 <- dcast(res2, psi ~ estimator, value.var = c("bias", "ratio"))
res2 <- res2[, .(psi, bias_FEPPML, bias_ABC, bias_SPJ, ratio_FEPPML, ratio_ABC, ratio_SPJ)]
print(xtable(res2, digits = 3), include.rownames = FALSE, booktabs = TRUE)
