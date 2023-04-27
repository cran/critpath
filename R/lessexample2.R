#' Dataset for the LESS method
#'
#' Example from  Miszczyńska D., Miszczyński M. "Wybrane metody badań operacyjnych" (2000, ISBN:83-907712-0-9).
#' In this dataset, the activities occur on the edges and a list of direct predecessors has been added.
#'
#' @format A data frame composed of predetermined columns:
#' \describe{
#'   \item{label}{activity label}
#'   \item{pred}{preceding activities}
#'   \item{time}{normal duration of the activity}
#'   \item{bound_time}{the shortest duration of the activity}
#'   \item{norm_cost}{normal cost of the activity}
#'   \item{bound_cost}{boundary cost of the activity}
#' }
"lessexample2"

