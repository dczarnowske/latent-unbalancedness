################################################################################
# empirical example
#
# date: 2022/09/09
################################################################################

# required packages
library(data.table)
library(forcats)
library(ggplot2)
library(ggthemes)
library(haven)
library(xtable)
source("helper_functions.R")

# prepare data set
# note: data from the WZ replication package required
# link: https://www.dropbox.com/s/8zchjamfg2bqipd/Code%20repository%20for%20Weidner%20and%20Zylkin%20JIE%202021.zip?dl=0
dt_full <- read_dta("baci_example.dta")
setDT(dt_full)
names_old <- c("trade", "FTA", "isoexp", "isoimp", "year")
names_new <- c("y", "x", "exp", "imp", "time")
setnames(dt_full, names_old, names_new)
dt_full[, exp := factor(exp)]
dt_full[, imp := factor(imp)]
dt_full[, time := factor(time)]
dt_full[, exp_time := interaction(exp, time)]
dt_full[, imp_time := interaction(imp, time)]
dt_full[, pair := interaction(exp, imp)]
dt_full <- dt_full[!(isic_id %in% c(5, 6, 31:35))]
dt_full[, isic_id := as.integer(factor(isic_id))]

# text -- footnote 10 (coal industry (10))
dt_coal <- dt_full[isic_id == 4]
dt_coal[, any(y > 0), by = .(exp, time)][, all(V1), by = exp][, summary(V1)]
dt_coal[, any(y > 0), by = .(imp, time)][, all(V1), by = imp][, summary(V1)]
dt_coal[, sum(y > 0), by =  .(imp, time)][, summary(V1)]
dt_coal <- dropObservations(dt_coal)
dt_coal[, exp_time := droplevels(exp_time)]
dt_coal[, imp_time := droplevels(imp_time)]
nrow(dt_coal) / dt_coal[, length(unique(imp_time))] # Ibar 
nrow(dt_coal) / dt_coal[, length(unique(exp_time))] # Jbar 

# for all industries plus aggregate trade
dt_list <- lapply(1:max(dt_full$isic_id), function(id, dt_full) {
  dt <- dt_full[isic_id == id]
  desc <- tolower(dt[, first(isicdescr)])
  id <- as.integer(dt[, first(isiccode)])
  exp <- dt[, length(unique(exp))]
  imp <- dt[, length(unique(imp))]
  time <- dt[, length(unique(time))]
  N <- max(exp, imp)
  n <- N^2 * time
  pa <- dt[, length(unique(exp_time))]
  pg <- dt[, length(unique(imp_time))]
  dt <- dropObservations(dt)
  dt[, exp_time := droplevels(exp_time)]
  dt[, imp_time := droplevels(imp_time)]
  pa_ast <- dt[, length(unique(exp_time))]
  pg_ast <- dt[, length(unique(imp_time))]
  n_ast <- nrow(dt)
  c_alpha <- pa_ast / pa
  c_gamma <- pg_ast / pg
  d <- n_ast / n
  data.frame(
    description = desc,
    isiccode    = id,
    n           = n,
    N           = N,
    time        = time,
    c_alpha     = c_alpha,
    c_gamma     = c_gamma,
    d           = d
  )
}, dt_full = dt_full)
dt <- do.call(rbind, dt_list)
setDT(dt)
dt[, isic := as.character(isiccode)]
dt[isic == "100", isic := "Agg"]
dt[, isic := factor(isic, levels = isic)]

# table -- supplement
table_vars <- c("description", "isic", "c_alpha", "c_gamma", "d")
dt_table <- dt[, table_vars, with = FALSE]
print(xtable(dt_table), include.rownames = FALSE, booktabs = TRUE)

# figure -- inverse bias
fig_vars <- c("isic", "c_alpha", "c_gamma", "d")
dt_fig <- dt[, fig_vars, with = FALSE]
dt_fig[, B := c_alpha / d]
dt_fig[, D := c_gamma / d]
dt_fig[, V := 1 / sqrt(d)]
dt_fig[, BD_max := pmax(B, D)]
dt_fig[isic != "Agg", summary(BD_max)]
dt_fig[isic != "Agg", summary(V)]
dt_fig[isic == "Agg", BD_max]
dt_fig[isic == "Agg", V]
dt_fig <- melt(
  dt_fig,
  id.vars       = "isic",
  measure.vars  = c("B", "D", "V"),
  variable.name = "component",
  value.name    = "scale"
)
new_levels <- c(
  "Exporter (B)"        = "B",
  "Importer (D)"        = "D",
  "Std. Deviation (V)" = "V"
)
dt_fig[, component := fct_recode(component, !!!new_levels)]
ggplot(dt_fig, aes(x = isic, y = scale, fill = component)) +
  geom_col(position = "dodge") +
  labs(
    x    = "2 Digit ISIC Industry Codes",
    y    = "Scaling Factor",
    fill = NULL
  ) +
  geom_hline(yintercept = 8, alpha = 0) +
  scale_fill_manual(
    values = c("Exporter (B)" = "#D9D9D9", "Importer (D)" = "#969696", "Std. Deviation (V)" = "#525252")
  ) +
  scale_y_continuous(breaks = seq(0, 9, 1)) +
  theme_few() +
  theme(legend.position = "bottom")
ggsave("scaling_bias.pdf", width = 8.27, height = 4.13)

# figure -- patterns
# the three industries with the largest share of missing trading pairs
dt_list <- lapply(c(10, 13, 16), function(id, dt_full) {
  dt <- dt_full[isiccode == id]
  dt <- dropObservations(dt)
  dt
}, dt_full = dt_full)
dt_fig <- rbindlist(dt_list)
dt_fig[, exp := droplevels(exp)]
dt_fig[, imp := droplevels(imp)]
dt_fig[, pair := droplevels(pair)]
dt_fig[, time := droplevels(time)]
dt_fig[, exp_id := as.integer(exp)]
dt_fig[, imp_id := as.integer(imp)]
dt_fig[, T_ij := .N, by = .(isiccode, pair)]
dt_fig[, isic := factor(isiccode)]
new_levels <- c(
  "Coal (10)"       = "10",
  "Metal Ores (13)" = "13",
  "Tobacco (16)"    = "16"
)
dt_fig[, isic := fct_recode(isic, !!!new_levels)]
dt_fig <- dt_fig[, .SD[1], by = .(isic, pair)]
ggplot(dt_fig, aes(exp_id, imp_id, fill = T_ij)) +
  facet_grid(~isic) +
  geom_point(shape = "circle filled") +
  labs(
    x    = "Exporter (i)",
    y    = "Importer (j)",
    fill = expression(T[ij])
  ) +
  scale_fill_distiller(palette = "Spectral", direction = 1) +
  theme_few() +
  theme(legend.position = "bottom")
ggsave("patterns.pdf", width = 8.27, height = 4.13)
