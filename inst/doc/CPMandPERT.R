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

## ---- fig.align = 'center', fig.cap = "Fig. 1. Graph for the cpmexample1 dataset"----
plot_graphAOA(cpmexample1)

## -----------------------------------------------------------------------------
x <- solve_pathAOA(cpmexample1, deterministic = TRUE)

## -----------------------------------------------------------------------------
# Schedule
x[2]
# Directive deadline
x[3]
# Critical activities
x[4]
# Free float and conditional float values
x[5]

## ---- fig.align = 'center', fig.cap = "Fig. 2. Critical path for the cpmexample1 dataset"----
plot_graphAOA(solved = x)

## ---- fig.align = 'center', fig.width = 6, fig.cap = "Fig. 3. Gantt chart for the cpmexample1 dataset"----
plot_gantt(x)

## ---- fig.align = 'center', fig.width = 6, fig.cap = "Fig. 4. ASAP chart for the cpmexample1 dataset"----
plot_asap(x)

## ---- fig.align = 'center', fig.width = 6, fig.cap = "Fig. 5. ALAP chart for the cpmexample1 dataset"----
plot_alap(x)

## ---- echo = FALSE------------------------------------------------------------
from <- c(1, 2, 3, 3, 3, 4, 5, 6, 7)
to <- c(2, 3, 4, 5, 6, 7, 7, 7, 8)
label <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")
opt_time <- c(3, 5, 5, 1, 6, 2, 5, 3, 4)
likely_time <- c(5, 7, 5, 6, 8, 6, 6, 5, 6)
pes_time <- c(7, 9, 8, 8, 10, 7, 7, 7, 8)
knitr::kable(cbind(as.numeric(from), to, label, opt_time, likely_time, pes_time), 
             col.names = c("Start. node", "End. node", "Name", "Optimistic dur.", 
                           "Most likely dur.", "Pessimistic dur."),
             align = "cccccc", caption = "Tab. 2. Data for the PERT model")

## -----------------------------------------------------------------------------
y <- solve_pathAOA(pertexample1, deterministic = FALSE)

## -----------------------------------------------------------------------------
# Schedule
y[2]
# Expected completion time
y[3]
# Standard deviation of the completion time
y[4]
# Critical activities
y[5]
# Free float and conditional float values
y[6]

## ---- fig.align = 'center',  fig.width = 6, fig.cap = "Fig. 6. Critical path for the pertexample1 dataset"----
plot_graphAOA(solved = y)

## ---- fig.align = 'center', fig.width = 6, fig.cap = "Fig. 7. Gantt chart for the pertexample1 dataset"----
plot_gantt(y)

## ---- fig.align = 'center', fig.width = 6, fig.cap = "Fig. 8. ASAP chart for the pertexample1 dataset"----
plot_asap(y)

## ---- fig.align = 'center', fig.width = 6, fig.cap = "Fig. 9. ALAP chart for the pertexample1 dataset"----
plot_alap(y)

## -----------------------------------------------------------------------------
# Risk-taker's schedule
PERT_newtime(new_prob = 0.3, y)

# Belayer's schedule
PERT_newtime(new_prob = 0.6, y)

## -----------------------------------------------------------------------------
PERT_newprob(new_DT = 30, y)

## ---- fig.align = 'center', fig.width = 6, fig.cap = "Fig. 10. Normal distribution for the pertexample1 dataset"----
plot_norm(y)

