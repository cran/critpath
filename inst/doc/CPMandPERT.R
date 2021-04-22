## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo = FALSE------------------------------------------------------------
from <- c(1, 1, 1, 2, 3, 4)
to <- c(2, 3, 4, 4, 4, 5)
label <- c("A", "B", "C", "D", "E", "F")
time <- c(3, 2, 5, 1, 1, 3)
knitr::kable(cbind(as.numeric(from), to, label, time), col.names = c("Start. node", 
                                                         "End. node",
                                                         "Name",
                                                         "Duration"),
             align = "cccc", caption = "Tab. 1. Data for the CPM model")

## ---- echo = FALSE, include = FALSE-------------------------------------------
library(critpath)
plot_graphAOA(cpmexample1)

## ---- fig.align = 'center'----------------------------------------------------
plot_graphAOA(cpmexample1)

## -----------------------------------------------------------------------------
x <- solve_pathAOA(cpmexample1, deterministic = TRUE)

## -----------------------------------------------------------------------------
# Schedule
x[2]
# Directive term
x[3]
# Critical activities
x[4]

## ---- fig.align = 'center'----------------------------------------------------
plot_crit_pathAOA(x)

## ---- fig.align = 'center', fig.width = 6, fig.cap = "Gantt chart"------------
plot_gantt(x)

## ---- echo = FALSE------------------------------------------------------------
from <- c(1, 2, 3, 3, 3, 4, 5, 6, 7)
to <- c(2, 3, 4, 5, 6, 7, 7, 7, 8)
label <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")
opt_time <- c(3, 5, 5, 1, 6, 2, 5, 3, 4)
likely_time <- c(5, 7, 5, 6, 8, 6, 6, 5, 6)
pes_time <- c(7, 9, 8, 8, 10, 7, 7, 7, 8)
knitr::kable(cbind(as.numeric(from), to, label, opt_time, likely_time, pes_time), 
             col.names = c("Start. node", "End. node", "Name", "Opt. time", 
                           "Most like. time", "Pess. time"),
             align = "cccccc", caption = "Tab. 2. Data for the PERT model")

## -----------------------------------------------------------------------------
y <- solve_pathAOA(pertexample1, deterministic = FALSE)

## -----------------------------------------------------------------------------
# Schedule
y[2]
# Expected directive term
y[3]
# Standard deviation of the directive term
y[4]
# Critical activities
y[5]

## ---- fig.align = 'center',  fig.width = 6------------------------------------
plot_crit_pathAOA(y)

## ---- fig.align = 'center', fig.width = 6-------------------------------------
plot_gantt(y)

## -----------------------------------------------------------------------------
# Risk-taker's schedule
qnorm(0.3, y[[3]], y[[4]])

# Belayer's schedule
qnorm(0.6, y[[3]], y[[4]])

## -----------------------------------------------------------------------------
pnorm(30, mean = y[[3]], sd = y[[4]])

## ---- fig.align = 'center', fig.width = 6-------------------------------------
plot_norm(y)

