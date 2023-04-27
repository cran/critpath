#' Dataset for the PERT method
#'
#' Example from  Miszczyńska D., Miszczyński M. "Wybrane metody badań operacyjnych" (2000, ISBN:83-907712-0-9).
#' 10 activities, 8 nodes.
#' In this dataset, the activities occur on the edges and a list of direct predecessors has been added.
#'
#' @format A data frame composed of predetermined columns:
#' \describe{
#'   \item{label}{activity label}
#'   \item{pred}{preceding activities}
#'   \item{opt_time}{optimistic duration of activity}
#'   \item{likely_time}{the most likely duration of the activity}
#'   \item{pes_time}{pesimistic duration of activity}
#' }
"pertexample2"

