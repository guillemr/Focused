library(microbenchmark)
library(furrr)

options(future.globals.maxSize = 8000 * 1024^2)

CORES <- 30
plan(multicore, workers = CORES)

source("Section_4_3_2_OCDlike_Simulation/helper_functions.R")
source("MdFOCuS_R_implementation/MdFocus_MeanGaussian_md.R")


ocd_runtime <- function(y, p) {
    ocd_det <- ocd_known(c(Inf, Inf, Inf), rep(0, p), rep(1, p))
    r <- ocd_detecting(y, ocd_det)
    r$t
}

run_one_sim <- function(n, p, seq) {
    y <- generate_sequence(n = n, p = p, cp = 500, magnitude = 0, dens = 0, seed = seq)
    y_t <- t(y)

    if(p < 5) {
    out <- microbenchmark(
        "MDFOCuS MR"=FocusCH(y_t, get_opt_cost = \(...) get_partial_opt(..., cost=cost_lr_partial0), threshold = rep(Inf, 3)),
        "MDFOCuS"=FocusCH(y_t, get_opt_cost = \(...) get_glo_opt(..., cost=lr_Focus0), threshold = Inf),
        "ocd"=ocd_runtime(y, p),
        times = 10,
        unit = "ms"
        )
    } else {
    out <- microbenchmark(
        "MDFOCuS 2d MR"=FocusCH_HighDim(y_t, get_opt_cost = \(...) get_partial_opt(..., cost=cost_lr_partial0, which_par = c(5, 25, p)), threshold = rep(Inf, 5)),
        "MDFOCuS 2d"=FocusCH_HighDim(y_t, get_opt_cost = \(...) get_glo_opt(..., cost=lr_Focus0), threshold = Inf),
        "ocd"=ocd_runtime(y, p),
        times = 10,
        unit = "ms"
        )
    }
    return(tibble(algorithm = out$expr, time = out$time, p = p, n = n))
}


### running a whole bunch of sims ###
library(tidyr)

run_lengths <- c(500, 1000, 5000, 1e4, 2e4)

tot_sims <- expand_grid(
    p = c(3, 50, 100),
    n = run_lengths,
    seq = 1:5
)

tot_res <- future_pmap(tot_sims, run_one_sim, .progress = T)
tot_res <- tot_res |> reduce(rbind)


save(tot_res, file = "runtime_sim.RData")

library(ggplot2)
library(scales)
out <- ggplot(data = tot_res, aes(x = n, y = time, group = algorithm, color = algorithm)) +
  geom_line(stat = "summary", fun.data = "mean_se") +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", width = 0.1) +
  facet_wrap(~p) +
  scale_x_log10(breaks = run_lengths) +
  scale_y_log10(name = "Time (nanoseconds)") +
  theme_minimal()

ggsave(out, file="runtimes.pdf", width = 10, height = 4)
