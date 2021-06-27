## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo = FALSE, include = FALSE-------------------------------------------
library(critpath)
plot_graphAOA(lessexample1)

## ---- fig.align = 'center'----------------------------------------------------
plot_graphAOA(lessexample1)

## ---- echo = FALSE------------------------------------------------------------
from <- c(1, 1, 2, 2, 3, 3, 4)
to <- c(2, 3, 3, 4, 4, 5, 5)
label <- c("A", "B", "C", "D", "E", "F", "G")
time <- c(4, 6, 2, 6, 3, 3, 5)
bound_time <- c(2, 3, 1, 4, 2, 3, 3)
norm_cost <- c(50, 40, 20, 100, 50, 25, 60)
bound_cost <- c(70, 55, 24, 130, 60, 25, 75)
knitr::kable(cbind(from, to, label, time, bound_time, norm_cost, bound_cost), col.names = c("Start. node", 
                                                         "End. node",
                                                         "Name",
                                                         "Normal duration",
                                                         "Boundary duration",
                                                         "Normal cost",
                                                         "Boundary cost"),
             align = "cccc", caption = "Tab. 1. Data for the LESS model")

## -----------------------------------------------------------------------------
z <- solve_lessAOA(lessexample1, 50, 15)

## -----------------------------------------------------------------------------
# Data frame with some data and partial results
z[2]
# List of critical activities
z[3]
# The total cost vector of all iterations
z[4]
# Minimum total cost
z[5]
# Directive term for normal time
z[6]
# The shortest directive term
z[7]

## ---- fig.align = 'center'----------------------------------------------------
plot_crit_pathAOA(z)

## ---- fig.align = 'center'----------------------------------------------------
plot_TC(z)

