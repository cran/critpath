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
             align = "cccc", caption = "Tab. 1. Data loaded based on the start and end nodes")

## ---- echo = FALSE, include = FALSE-------------------------------------------
library(critpath)
plot_graphAOA(cpmexample1)

## ---- fig.align = 'center', fig.cap = "Fig. 1. Graph for the cpmexample1 dataset"----
plot_graphAOA(cpmexample1)

## ---- echo = FALSE------------------------------------------------------------
label <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
pred <- c(NA, NA, "A", "A", "A", "C,D", "D", "B,E", "H", "F,G,I")
time <- c(6, 2, 4, 6, 3, 2, 5, 3, 2, 2)
knitr::kable(cbind(label, pred, time), col.names = c("Name", 
                                                         "Predecessors",
                                                         "Duration"),
             align = "ccc", caption = "Tab. 2. Data loaded from the list of activity predecessors")

## ---- echo = FALSE, include = FALSE-------------------------------------------
library(critpath)
plot_graphAOA(cpmexample2, predecessors = TRUE)

## ---- fig.align = 'center', fig.cap = "Fig. 2. Graph for the cpmexample2 dataset"----
plot_graphAOA(cpmexample2, predecessors = TRUE)

