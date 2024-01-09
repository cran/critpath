#===============================================================================
# Function that reads data in the DiagrammeR package format.
# Changing labels to character mode to display them on the graph.
# Data frame for the CPM method. TF stands for total float.

read_cpmAOA <- function(yourdata){
  # Check if the DiagrammeR package is loaded. If not, load it.
  pckg_check("DiagrammeR")

  stopifnot("The data frame for CPM has the wrong number of columns" = ncol(yourdata) == 4)

  create_edge_df(from = yourdata[,1], to = yourdata[,2], label = yourdata[,3],
                   time = yourdata[,4], TF = rep(0, nrow(yourdata)))
}
#===============================================================================
# Function reads data in the DiagrammeR package format.
# Change labels to character mode to display them on the graph.
# Data frame for the PERT method. TF stands for total float.

read_pertAOA <- function(yourdata, pert_param){
  # Check if the DiagrammeR package is loaded. If not, load it.
  pckg_check("DiagrammeR")

  stopifnot("The data frame for PERT has the wrong number of columns" = ncol(yourdata) == 6)

  create_edge_df(from = as.integer(yourdata[,1]), to = as.integer(yourdata[,2]),
                   label = as.character(yourdata[,3]),
                   time = PERT_mu(pert_param, yourdata[,4], yourdata[,5], yourdata[,6]),
                   timevar = PERT_var(pert_param, yourdata[,4], yourdata[,5], yourdata[,6]),
                   TF = rep(0, nrow(yourdata)),)
}
#===============================================================================
# Creates relations from an input data frame.

make_relationsAOA <- function(yourdata, deterministic, predecessors, pert_param){
  if (predecessors == FALSE){
    # Sort the data frame by node numbers.
    yourdata <- yourdata[order(yourdata[,1], yourdata[,2]),]
    if (deterministic == TRUE){
      read_cpmAOA(yourdata)
    }else{
      read_pertAOA(yourdata, pert_param)
    }
  }else{
    # After switching from the activities immediately preceding.
    AOA_df <- input_predAOA(yourdata)
    if (deterministic == TRUE){
      # New dataframe after concatenating the predecessor and duration.
      mergedAOA_df <- merge_pred_detAOA(yourdata, AOA_df)
      read_cpmAOA(mergedAOA_df)
    }else{
      # New dataframe after concatenating the predecessor and duration.
      mergedAOA_df <- merge_pred_probAOA(yourdata, AOA_df)
      read_pertAOA(mergedAOA_df, pert_param)
    }
  }

}
#===============================================================================
# Creates a data frame to hold the vertices. Initially, they are all zero.

make_nodesAOA <- function(input_tab){
  # Determining the number of vertices.
  nodes_num <- max(input_tab$to)

  create_node_df(
    n = nodes_num,
    type = c(1:nodes_num),
    label = TRUE,
    ES = rep(0, nodes_num),
    LF = rep(0, nodes_num)
  )
}
#===============================================================================
# Prepares an empty frame to store the schedule. There are 2 variants of the data frame.

create_empty_schedule <- function(relations, deterministic){
  if (deterministic == TRUE){
    data.frame(
      Name = relations$label,
      Duration = relations$time,
      ESij = rep(0, nrow(relations)),
      LSij = rep(0, nrow(relations)),
      EFij = rep(0, nrow(relations)),
      LFij = rep(0, nrow(relations)),
      TFij = rep(0, nrow(relations)),
      Crit = rep(c(" "), nrow(relations))
    )
  }else{
    data.frame(
      Name = relations$label,
      Duration = relations$time,
      Var = relations$timevar,
      ESij = rep(0, nrow(relations)),
      LSij = rep(0, nrow(relations)),
      EFij = rep(0, nrow(relations)),
      LFij = rep(0, nrow(relations)),
      TFij = rep(0, nrow(relations)),
      Crit = rep(c(" "), nrow(relations))
    )
  }
}

#===============================================================================
# Prepares and fills a data frame that stores the schedule of the CPM method.

scheduleAOA <- function(relations, deterministic){
  ES <- LF <- TF <- NULL
  # Create a data frame to hold the vertices.
  vertices <- make_nodesAOA(relations)

  # Create a data frame graph from vertices and relations.
  yourgraph <- create_graph(nodes_df = vertices,
                           edges_df = relations, directed = TRUE)

  yourschedule <- create_empty_schedule(relations, deterministic)
  # Calculates ES values for vertices. The loop goes through SINGLE activities, from the first to the last one.
  # Compares the ES value at a given vertex with the sum of the ES value at the predecessor vertex
  # and the activity duration time.
  for(i in  1:c(nrow(relations))){
    yourgraph <- set_node_attrs(yourgraph,
                           node_attr = ES,
                           values = max(get_node_attrs(yourgraph,
                                                       node_attr = ES,
                                                       nodes = c(relations$to[i])),
                                        get_node_attrs(yourgraph, node_attr = ES,
                                                       nodes = c(relations$from[i])) + c(relations$time[i])),
                           nodes = c(relations$to[i]))
  }

  # In all vertices, we change the value of the LF attribute to ES from the last vertex.
  yourgraph <- set_node_attrs(yourgraph,
                         node_attr = LF,
                         values = get_node_attrs(yourgraph,
                                                 node_attr = ES,
                                                 nodes = c(nrow(vertices))))

  # Calculates LF. The loop goes through SINGLE activities, from the last to the first one.
  # Compares the LF value in the previous vertex with the difference between
  # the LF value of the current vertex and the duration of the activity.
  for(i in  c(nrow(relations)):1){
    yourgraph <- set_node_attrs(yourgraph,
                           node_attr = LF,
                           values = min(get_node_attrs(yourgraph,
                                                       node_attr = LF,
                                                       nodes = c(relations$from[i])),
                                        get_node_attrs(yourgraph,
                                                       node_attr = LF,
                                                       nodes = c(relations$to[i])) - c(relations$time[i])),
                           nodes = c(relations$from[i]))
  }
 # Completes ES and LF values in the schedule.
  for (i in 1:c(nrow(relations))){
    yourschedule$ESij[i] <- yourgraph %>% get_node_attrs(node_attr = ES,
                                                    nodes = c(relations$from[i]))
    yourschedule$LFij[i] <- yourgraph %>% get_node_attrs(node_attr = LF,
                                                    nodes = c(relations$to[i]))
  }
  # Completes the rest of the schedule.
  yourschedule$LSij <- yourschedule$LFij - yourschedule$Duration
  yourschedule$EFij <- yourschedule$ESij + yourschedule$Duration
  yourschedule$TFij <- yourschedule$LFij - yourschedule$ESij - yourschedule$Duration
  yourschedule$Crit[which(yourschedule$TFij == 0)] <- c("*")

  # Completes TF for the edges data frame.
  yourgraph <- set_edge_attrs(yourgraph, edge_attr = TF, values = yourschedule$TFij)

  # Extract values of the ES attributes for all nodes
  ESnodes <- get_node_attrs(yourgraph, node_attr = ES)

  # Extract values of the LF attributes for all nodes
  LFnodes <- get_node_attrs(yourgraph, node_attr = LF)

  # Calculate spare/free and conditional slack of time
  AddInfo <- data.frame(
    Name = relations$label,
    FST = ESnodes[relations$to] - ESnodes[relations$from] - yourschedule$Duration,
    CST = LFnodes[relations$to] - LFnodes[relations$from] - yourschedule$Duration
  )

  # Messages after the computation is completed.
  if (deterministic == TRUE){
    # CPM
    cat("Completion time: ", max(yourschedule$LFij), "\n")
  }else{
    # PERT
    cat("Expected compl. time distribution: N(",max(yourschedule$LFij),",", sqrt(sum(yourschedule$Var[which(yourschedule$TFij == 0)])),")\n")
  }

  # Creates a list keeping the graph and schedule.
  if (deterministic == TRUE){
    # CPM.
    list(graphAOA = yourgraph, schedule = yourschedule,
      ComplTi = max(yourschedule$LFij),
      CritAct = yourschedule$Name[which(yourschedule$TFij == 0)],
      AddSlacks = AddInfo
      )
  }else{
    # PERT.
    list(graphAOA = yourgraph, schedule = yourschedule,
      ComplTi = max(yourschedule$LFij),
      SDevTi = sqrt(sum(yourschedule$Var[which(yourschedule$TFij == 0)])),
      CritAct = yourschedule$Name[which(yourschedule$TFij == 0)],
      AddSlacks = AddInfo)
  }
}

#===============================================================================
# Function that checks if the necessary packages are installed.

pckg_check <- function(pckg_name){
  if (!requireNamespace(pckg_name, quietly = TRUE)){
    napis <- paste("Package",pckg_name,"is needed for this function to work. Please install it.")
    stop(napis, call. = FALSE)
  }
}

#===============================================================================

# Function solving CPM and PERT methods. The argument is a data frame with information about the problem

#' Finds a solution using CPM and PERT methods. Relationships between activities can be given as a list of predecessors or start and end node numbers.
#'
#' @param input_data Data frame containing the structure of the graph and the duration of the activity.
#'   For the CPM method and start/end nodes you need 4 columns (the order is important, not the name of the column):
#'   \enumerate{
#'   \item \code{from} The number of the node where the activity starts.
#'   \item \code{to} The number of the node where the activity ends.
#'   \item \code{label} Activity labels.
#'   \item \code{time} Activities durations.
#'   }
#'   For the CPM method and predecessors list you need 3 columns (the order is important, not the name of the column):
#'   \enumerate{
#'   \item \code{label} Activity labels.
#'   \item \code{pred} List of predecessors.
#'   \item \code{time} Activities durations.
#'   }
#'   For the PERT method and start/end nodes you need 6 columns (the order is important, not the name of the column):
#'   \enumerate{
#'   \item \code{from} The number of the node where the activity starts.
#'   \item \code{to} The number of the node where the activity ends.
#'   \item \code{label} Activity labels.
#'   \item \code{opt_time} Optimistic duration of activities.
#'   \item \code{likely_time} The most likely duration of the activity.
#'   \item \code{pes_time} Pessimistic duration of activities.
#'   }
#'   For the PERT method and predecessors list you need 5 columns (the order is important, not the name of the column):
#'   \enumerate{
#'   \item \code{label} Activity labels.
#'   \item \code{pred} List of predecessors.
#'   \item \code{opt_time} Optimistic duration of activities.
#'   \item \code{likely_time} The most likely duration of the activity.
#'   \item \code{pes_time} Pessimistic duration of activities.
#'   }
#' @param deterministic A logical parameter specifying the solution method.
#'   If set to \code{TRUE} (default), the CPM method is used. If is set to \code{FALSE}, the PERT method is used.
#' @param predecessors TRUE if the user data contains a list of immediately preceding activities
#'   If set to \code{FALSE} (default), start nad end nodes are used. If is set to \code{TRUE}, predecessors list is used.
#' @param pert_param A parameter that controls the method of calculating the expected value and variance in the PERT method.
#'   0 - classic formula (default), 1 - 1st and 99th percentile of the beta distribution, 2 - 5th and 95th percentile of the beta distribution,
#'   3 - 5th and 95th percentiles of the beta distribution with modification by (Perry and Greig, 1975), 4 - Extended Pearson's and Tukey's formula
#'   (Pearson and Tukey, 1965), 5 - Golenko-Ginzburg's full formula (Golenko-Ginzburg, 1988), 6 - Golenko-Ginzburg's reduced formula
#'   (Golenko-Ginzburg, 1988), 7 - Farnum's and Stanton's formula (Farnum and Stanton, 1987).
#' @return The list is made of a graph, schedule and selected partial results.
#' @examples
#' x <- solve_pathAOA(cpmexample1, deterministic = TRUE)
#' y <- solve_pathAOA(pertexample1, deterministic = FALSE)
#' x <- solve_pathAOA(cpmexample2, deterministic = TRUE, predecessors = TRUE)
#' y <- solve_pathAOA(pertexample2, deterministic = FALSE, predecessors = TRUE)
#' @import DiagrammeR
#' @import utils
#' @export
solve_pathAOA <- function(input_data, deterministic = TRUE, predecessors = FALSE, pert_param = 0){
  # Check if the DiagrammeR package is loaded. If not, load it.
  pckg_check("DiagrammeR")

  solved_pathAOA <- vector("list", length = 4)
  # Loads data and creates a relations data frame.
  relations <- make_relationsAOA(input_data, deterministic, predecessors, pert_param)
  # Creates graph and schedule and writes them into list.
  scheduleAOA(relations, deterministic)
}

#===============================================================================

# Presentation of the critical path on a graph.

plot_crit_pathAOA <- function(yourlist, fixed_seed = 23){
  TF <- color <- time <- style <- NULL
  # Check if DiagrammeR package is loaded. If not, load it.
  pckg_check("DiagrammeR")

  # Check if yourlist is a list
  stopifnot("The function requires a list" = is.list(yourlist))

  # Temporarily remembered list items.
  yourgraph <- yourlist[[1]]

  # Marking of dummy activities dotted line
  if (min(get_edge_attrs(yourgraph, edge_attr = time)) == 0){
    yourgraph <- yourgraph %>%
      select_edges(conditions = time == 0) %>%
      set_edge_attrs_ws(edge_attr = style, value = "dotted")
  }

  # Color marking of critical activities.
  yourgraph <- yourgraph %>%
    clear_selection() %>%
    select_edges(conditions = TF == 0) %>%
    set_edge_attrs_ws(edge_attr = color, value = "red")
  # Set random seed to fixed value
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(fixed_seed)
  render_graph(yourgraph, layout = "fr")
}
#===============================================================================

# A function that draws a Gantt chart.

#'A Gantt chart
#'
#' @param yourlist List of objects that make up the solution to the project management problem.
#' @param show_dummy Decides whether dummy activities should be included in the chart. If so, set it to TRUE (set to FALSE by default).
#' @param bar_size Thickness of the bar drawn for activity (set to 10 by default).
#' @return Draws a Gantt chart broken down into critical ("CR") and non-critical ("NC") activities.
#'   Marks total floats.
#' @examples
#' x <- solve_pathAOA(cpmexample1, deterministic = TRUE)
#' plot_gantt(x)
#' @import ggplot2
#' @import reshape2
#' @export
plot_gantt <- function(yourlist, show_dummy = FALSE, bar_size = 10){
  Name <- TFij <- value <- NULL
  # Check if ggplot2 package is loaded. If not, load it.
  pckg_check("ggplot2")

  # Check if reshape2 package is loaded. If not, load it.
  pckg_check("reshape2")

  # Check if yourlist is a list
  stopifnot("The function requires a list" = is.list(yourlist))

  # Identify dummy activities
  dummy_rows <- which(yourlist[[2]]$Duration == 0)

  if (length(dummy_rows) == 0){
    schedule <- yourlist[[2]]
  }else if (show_dummy == TRUE){
    schedule <- yourlist[[2]]
  }else{
    schedule <- yourlist[[2]][-dummy_rows,]
  }

  # Sorting the schedule according to total float.
  dftmp <- schedule[order(schedule$TFij),]

  # The TFij column gives the critical activities the symbol "CR" and the non-critical activities "NC".
  dftmp$TFij[dftmp$TFij > 0] <- c("NC")
  dftmp$TFij[dftmp$TFij == 0] <- c("CR")

  # Gantt chart
  melthar <- melt(dftmp, measure.vars = c("ESij", "EFij"))
  melthar2 <- melt(dftmp, measure.vars = c("EFij", "LFij"))

  ggplot(melthar, aes(value, Name, colour = TFij)) +
    geom_line(size = bar_size) +
    ylab(NULL) +
    xlab(NULL) +
    theme_bw() +
    theme(legend.title=element_blank())

  ggplot() +
    geom_line(melthar, mapping = aes(value, Name, colour = TFij), size = bar_size) +
    geom_line(melthar2, mapping = aes(value, Name), colour ="cyan", size = bar_size) +
    ylab(NULL) +
    xlab(NULL) +
    theme_bw() +
    theme(legend.title=element_blank())
}

#===============================================================================
# Normal distribution plot for the PERT Method

#' The cumulative distribution function of the normal distribution
#'
#' @param yourlist List of objects making up the solution to the project management problem
#' @return Draws a graph of the normal distribution with the expected directive term from
#'   the PERT method and the standard deviation for this term. The chart also includes lines indicating
#'   the schedules of the risk-taker and the belayer.
#' @examples
#' y <- solve_pathAOA(pertexample1, deterministic = FALSE)
#' plot_norm(y)
#' @import ggplot2
#' @export
plot_norm <- function(yourlist){
  z <- df <- x <- NULL

  # Check if ggplot2 package is loaded. If not, load it.
  pckg_check("ggplot2")

  # Check if yourlist is a list
  stopifnot("The function requires a list" = is.list(yourlist))

  # Risk-taker schedule
  hryz <- stats::qnorm(0.3, yourlist[[3]], yourlist[[4]])
  #Belayer schedule
  hase <- stats::qnorm(0.6, yourlist[[3]], yourlist[[4]])
  # Check if ggplot2 package is loaded. If not, load it.
  pckg_check("ggplot2")
  # Variable that holds the random values plus/minus 3 variations
  z <- seq(yourlist[[3]]-3*yourlist[[4]], yourlist[[3]]+3*yourlist[[4]], 0.01)
  df <- data.frame(x = z, y = stats::dnorm(z))
  df <- cbind(df, stats::pnorm(z))
  ggplot(df, aes(x, y = stats::pnorm(z, yourlist[[3]], yourlist[[4]]))) +
    geom_line() +
    labs(y = "P(DT<=x)") +
    geom_vline(xintercept = hryz, color='blue')+
    geom_vline(xintercept = hase, color='blue')+
    theme_bw()
}

#===============================================================================
# Function that draws a graph before solving a problem.

plot_dataAOA <- function(input_data, predecessors, fixed_seed){
  time <- style <- NULL

  if (predecessors == TRUE){
    result_df <- input_predAOA(input_data)
    input_data <- merge_pred_detAOA(input_data, result_df)
  }
  relations <- create_edge_df(from = input_data[,1], to = input_data[,2],
                              label = as.character(input_data[,3]), time = input_data[,4])
  vertices <- make_nodesAOA(relations)
  # Build a graph from vertex and relationship data frames.
  yourgraph <- create_graph(nodes_df = vertices,
                            edges_df = relations, directed = TRUE)
  # Mark the dummy activities with a dotted line
  if (min(input_data[,4])==0){
    yourgraph <- yourgraph %>%
      select_edges(conditions = time == 0) %>%
      set_edge_attrs_ws(edge_attr = style, value = "dotted")
  }
  # Set random seed to fixed value
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(fixed_seed)
  render_graph(yourgraph, layout = "fr")
}

#===============================================================================

# A function that draws an ASAP (As Soon As Possible) chart.

#'An ASAP chart
#'
#' @param yourlist List of objects that make up the solution to the project management problem.
#' @param show_dummy Decides whether dummy activities should be included in the chart. If so, set it to TRUE (set to FALSE by default).
#' @param bar_size Thickness of the bar drawn for activity (set to 10 by default).
#' @return Draws an ASAP (activities start and finish As Soon As Possible) chart broken down into critical ("CR") and non-critical ("NC") activities.
#'   Marks total floats.
#' @examples
#' x <- solve_pathAOA(cpmexample1, deterministic = TRUE)
#' plot_asap(x)
#' @import ggplot2
#' @import reshape2
#' @export
plot_asap <- function(yourlist, show_dummy = FALSE, bar_size = 10){
  Name <- FST <- TFij <- value <- NULL
  # Check if ggplot2 package is loaded. If not, load it.
  pckg_check("ggplot2")

  # Check if reshape2 package is loaded. If not, load it.
  pckg_check("reshape2")

  # Check if yourlist is a list
  stopifnot("The function requires a list" = is.list(yourlist))

  # Identify dummy activities
  dummy_rows <- which(yourlist[[2]]$Duration == 0)

  if (length(dummy_rows) == 0){
    schedule <- yourlist[[2]]
    addslacks <- yourlist[["AddSlacks"]]
  }else if (show_dummy == TRUE){
    schedule <- yourlist[[2]]
    addslacks <- yourlist[["AddSlacks"]]
  }else{
    schedule <- yourlist[[2]][-dummy_rows,]
    addslacks <- yourlist[["AddSlacks"]][-dummy_rows,]
  }

  # Create temporary schedule with additional column
  schedule <- data.frame(schedule, ResTimeij = schedule$EFij + addslacks$FST)

  # Sorting the schedule according to total float.
  dftmp <- schedule[order(schedule$TFij),]

  # The TFij column gives the critical activities the symbol "CR" and the non-critical activities "NC".
  dftmp$TFij[dftmp$TFij > 0] <- c("NC")
  dftmp$TFij[dftmp$TFij == 0] <- c("CR")

  # Gantt chart as ASAP
  melthar <- melt(dftmp, measure.vars = c("ESij", "EFij"))
  melthar2 <- melt(dftmp, measure.vars = c("EFij", "ResTimeij"))

  ggplot(melthar, aes(value, Name, colour = TFij)) +
    geom_line(size = bar_size) +
    ylab(NULL) +
    xlab(NULL) +
    theme_bw() +
    theme(legend.title=element_blank())

  ggplot() +
    geom_line(melthar, mapping = aes(value, Name, colour = TFij), size = bar_size) +
    geom_line(melthar2, mapping = aes(value, Name), colour ="cyan", size = bar_size) +
    ylab(NULL) +
    xlab(NULL) +
    theme_bw() +
    theme(legend.title=element_blank())

}

#===============================================================================

# A function that draws an ALAP (As Late As Possible) chart.

#'An ALAP chart
#'
#' @param yourlist List of objects that make up the solution to the project management problem.
#' @param show_dummy Decides whether dummy activities should be included in the chart. If so, set it to TRUE (set to FALSE by default).
#' @param bar_size Thickness of the bar drawn for activity (set to 10 by default).
#' @return Draws an ALAP (activities start and finish As Late As Possible) chart broken down into critical ("CR") and non-critical ("NC") activities.
#'   Marks total float.
#' @examples
#' x <- solve_pathAOA(cpmexample1, deterministic = TRUE)
#' plot_alap(x)
#' @import ggplot2
#' @import reshape2
#' @export
plot_alap <- function(yourlist, show_dummy = FALSE, bar_size = 10){
  Name <- CST <- TFij <- value <- NULL
  # Check if ggplot2 package is loaded. If not, load it.
  pckg_check("ggplot2")

  # Check if reshape2 package is loaded. If not, load it.
  pckg_check("reshape2")

  # Check if yourlist is a list
  stopifnot("The function requires a list" = is.list(yourlist))

  # Identify dummy activities
  dummy_rows <- which(yourlist[[2]]$Duration == 0)

  if (length(dummy_rows) == 0){
    schedule <- yourlist[[2]]
    addslacks <- yourlist[["AddSlacks"]]
  }else if (show_dummy == TRUE){
    schedule <- yourlist[[2]]
    addslacks <- yourlist[["AddSlacks"]]
  }else{
    schedule <- yourlist[[2]][-dummy_rows,]
    addslacks <- yourlist[["AddSlacks"]][-dummy_rows,]
  }

  # Create temporary schedule with additional column
  schedule <- data.frame(schedule, ResTimeij = schedule$LSij - addslacks$CST)

  # Sorting the schedule according to total float.
  dftmp <- schedule[order(schedule$TFij),]

  # The TFij column gives the critical activities the symbol "CR" and the non-critical activities "NC".
  dftmp$TFij[dftmp$TFij > 0] <- c("NC")
  dftmp$TFij[dftmp$TFij == 0] <- c("CR")

  # Gantt chart as ASAP
  melthar <- melt(dftmp, measure.vars = c("LFij", "LSij"))
  melthar2 <- melt(dftmp, measure.vars = c("ResTimeij", "LSij"))

  ggplot(melthar, aes(value, Name, colour = TFij)) +
    geom_line(size = bar_size) +
    ylab(NULL) +
    xlab(NULL) +
    theme_bw() +
    theme(legend.title=element_blank())

  ggplot() +
    geom_line(melthar, mapping = aes(value, Name, colour = TFij), size = bar_size) +
    geom_line(melthar2, mapping = aes(value, Name), colour ="cyan", size = bar_size) +
    ylab(NULL) +
    xlab(NULL) +
    theme_bw() +
    theme(legend.title=element_blank())

}
#===============================================================================
# Function that draws a graph before or after solving a problem.

#' A graph of connections between nodes
#'
#' @param input_data Data frame describing the problem.
#' @param solved List of objects that make up the solution to the project management problem.
#' @param predecessors TRUE if the user data contains a list of immediately preceding activities
#' @param fixed_seed Optional parameter setting random seed to user value to get similar looking plots each time the function is run (set to 23 by default).
#' @return The function draws a graph showing dependencies between nodes. The "solved" parameter determines whether there is a critical path in the graph.
#'   In that case, you must solve the problem first. In the examples below, the function first draws the graph only on the basis of the data frame and then
#'   after determining the critical path.
#' @examples
#' plot_graphAOA(cpmexample1)
#' x <- solve_pathAOA(cpmexample1, TRUE)
#' plot_graphAOA(solved = x)
#' @import DiagrammeR
#' @export
plot_graphAOA <- function(input_data, predecessors = FALSE, solved = NULL, fixed_seed = 23){
  # Check if DiagrammeR package is loaded. If not, load it.
  pckg_check("DiagrammeR")

  if(is.null(solved)){
    plot_dataAOA(input_data, predecessors, fixed_seed)
  }else{
    plot_crit_pathAOA(solved, fixed_seed)
  }
}

#===============================================================================

# A function that computes a new directive term for a given probability

#' A new directive term for any probability
#'
#' @param new_prob Probability of the project completion. Default set to 0.5.
#' @param yourlist List of objects that make up the solution to the project management problem.
#' @return This function computes a new directive term for a probability given by the user. A normal distribution was assumed.
#' @examples
#' y <- solve_pathAOA(pertexample1, deterministic = FALSE)
#' PERT_newtime(new_prob = 0.3, y)
#' @export

PERT_newtime <- function(new_prob = 0.5, yourlist){
  # Check the entered parameters
  stopifnot("Value not in range [0,1]" = new_prob >= 0,
            "Value not in range [0,1]" = new_prob <= 1,
            "The function requires a list" = is.list(yourlist))
  DT <- stats::qnorm(new_prob, yourlist[[3]], yourlist[[4]])
  cat("New expected compl. time: ", DT, "\n")
  list(probDT = new_prob, newDT = DT)
}
#===============================================================================

# A function that calculates the probability of completing a project based on a given deadline.

#' Probability for the given directive term
#'
#' @param new_DT The given project completion date. The parameter must be greater than zero.
#' @param yourlist List of objects that make up the solution to the project management problem.
#' @return This function calculates the probability of completing the project within the time specified by the user. A normal distribution was assumed.
#' @examples
#' y <- solve_pathAOA(pertexample1, deterministic = FALSE)
#' PERT_newprob(new_DT = 30, y)
#' @export

PERT_newprob <- function(new_DT, yourlist){
  # Check the entered parameters
  stopifnot("New term <= 0" = new_DT > 0, "The function requires a list" = is.list(yourlist))
  probDT <- stats::pnorm(new_DT, mean = yourlist[[3]], sd = yourlist[[4]])
  cat("Prob. of completion: ", probDT, "\n")
  list(newDT = new_DT, prob_compl = probDT)
}
#===============================================================================
# Creates a data frame from predecessors list.

#' @import stringr

# Prefer over specified package or packages
# pckg_check("conflicted")
# conflict_prefer("filter", "dplyr", "stats")
# conflict_prefer("lag", "dplyr", "stats")
#conflicts_prefer(
#  dplyr::filter(),
#  dplyr::lag(),
#)

input_predAOA <- function(yourdata){
  occur_num <- Step2Col3 <- Step5Col4 <- dup_id <- NULL

  # Check if the strigr package is loaded. If not, load it.
  pckg_check("stringr")

  # Check if the dplyr package is loaded. If not, load it.
  pckg_check("dplyr")

  # Check if the names in the first column are unique
  stopifnot("Not all names are unique" = length(unique(yourdata[,1])) == nrow(yourdata))

  # Remove spaces from all columns and create a new dataframe.
  temp_df  <- as.data.frame(sapply(yourdata[,1:2], function(x) gsub(" " , "",x)))
  # Rename the dataframe columns.
  colnames(temp_df)[1:2] <- c("Activity", "Step1Col1")

  # Duplicate the immediate predecessor column (column 1) into column 2.
  temp_df <- data.frame(temp_df, Step2Col2 = temp_df$Step1Col1)

  # Create a list of all predecessors without the NA values, which means no predecessor.
  preced_constr_list <- unique(temp_df$Step1Col1[!is.na(temp_df$Step1Col1)])

  # Which activities have no predecessors?
  nopredec <- which(is.na(temp_df$Step2Col2))
  temp_df$Step2Col2[nopredec[-1]] <- c("-")

  # Add a column for start node numbers.
  temp_df <- data.frame(temp_df, Step2Col3 = 0)

  # Activities with NA will get number 1.
  temp_df$Step2Col3[nopredec] <- 1

  # Actions that have predecessors will be replaced with a "-" symbol.
  for (i in seq_along(preced_constr_list)) {
    # Row numbers where the given predecessor occurred.
    predec <- which(temp_df$Step2Col2 == preced_constr_list[i])
    temp_df$Step2Col2[predec[-1]] <- c("-")
    # starting node numbers
    temp_df$Step2Col3[predec] <- i + 1
  }

  # Add a column for dummy activities.
  temp_df <- data.frame(temp_df, Step3Col4 = temp_df$Step2Col2)

  # Auxiliary variables.
  col4 <- NULL
  col5 <- NULL
  col6 <- NULL
  counter <- 1
  dummy_names <- NULL
  # Point out and add dummy activities.
  for (i in 1:c(nrow(yourdata)-1)){
    # The number of occurrences of the action name in each row of the column.
    occur_num <- str_count(temp_df$Step2Col2, temp_df$Activity[i])

    # each activity can only occur once in a row
    if(sum(occur_num, na.rm = TRUE) > 0){
      stopifnot("Validate the list of predecessors" = max(occur_num, na.rm = TRUE) == 1)
    }

    # Did the activity only appear once?
    if(sum(occur_num, na.rm = TRUE) > 1){
      # In which rows did the activity occur?
      # Check the length of the expression (after removing the commas) for the indicated activities.
      # Replace the number of occurrences with the length of the expression.
      occur_num[which(occur_num == 1)] <- str_length(str_remove(temp_df$Step2Col2[which(occur_num == 1)], ","))

      # Replace dummy activity name.
      temp_df$Step3Col4[which(occur_num > 1)] <- str_replace(temp_df$Step3Col4[which(occur_num > 1)], temp_df$Activity[i], str_c("Dummy",temp_df$Activity[i]))

      # Remember the number and name of the activity that became dummy
      dummy_names[counter] <- str_c("Dummy",temp_df$Activity[i])
      col4[counter] <- i
      col5[counter] <- temp_df$Activity[i]
      # Remember the node number of the successor of the dummy activity.
      col6[counter] <- temp_df$Step2Col3[which(occur_num > 1)]
      counter <- counter + 1
    }
  }

  # Add an empty column for the numbers of the nodes that complete the activity.
  temp_df <- data.frame(temp_df, Step5Col4 = 0)

  # Remove Dummy's name so it doesn't appear in search results.
  bufcol <- temp_df$Step3Col4
  bufcol <- str_replace(bufcol, "Dummy.", "-")

  # Determine successor node numbers.
  for (i in 1:c(nrow(yourdata)-1)){
    # Does the activity appear in the list at least once?
    occur_list <- str_which(bufcol, temp_df$Activity[i])
    if(length(occur_list) > 0){
      # Remember the successor node number.
      temp_df$Step5Col4[i] <- temp_df$Step2Col3[occur_list[1]]
    }else{
      # The absence of a node for an activity means that it is dummy or ending.
      if(length(str_which(col5, temp_df$Activity[i])) > 0){
        temp_df$Step5Col4[i] <- max(c(temp_df$Step2Col3,temp_df$Step5Col4), na.rm = TRUE)+1
      }
    }
  }

  # Add the node number that ends the last activity.
  last_node <- which(temp_df$Step5Col4 == 0)
  temp_df$Step5Col4[last_node] <- max(c(temp_df$Step2Col3,temp_df$Step5Col4), na.rm = TRUE)+1

  # Create a dummy name vector.
  # dummy_names <- rep(c("Dummy"), length(col4))

  # Add new rows to the dataframe.
  if(!is.null(dummy_names)){
    temp_df[c(nrow(temp_df)+1):c(nrow(temp_df) + length(dummy_names)),] <- NA

    # Complete with data on dummy activities.
    temp_df[c(nrow(yourdata)+1):nrow(temp_df), 1] <- dummy_names
    temp_df[c(nrow(yourdata)+1):nrow(temp_df), 4] <- col4
    temp_df[c(nrow(yourdata)+1):nrow(temp_df), 5] <- col5
    temp_df[c(nrow(yourdata)+1):nrow(temp_df), 6] <- col6
  }

  # Additional dummy activities.
  # Add columns to indicate duplicate rows with start and end node numbers.
  temp_df <- temp_df %>%
    dplyr::group_by(Step2Col3, Step5Col4) %>%
    dplyr::mutate(num_dups = dplyr::n(),
           dup_id = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(is_duplicated = dup_id > 1)

  # If there are duplicate rows, change the node numbers.
  if(length(which(temp_df$is_duplicated)) > 0){
    # How many rows contain duplicate values?
    dup_rows <- c(1:length(which(temp_df$is_duplicated)))

    # Save the numbers of dummy activity nodes before changing their values.
    nodes_before <- temp_df$Step5Col4[which(temp_df$is_duplicated)]

    # Substitute node numbers for dummy activities identified.
    temp_df$Step5Col4[which(temp_df$is_duplicated)] <- max(temp_df$Step5Col4)+dup_rows

    # Generate a vector of dummy activity names.
    # dummy_names <- rep(c("Dummy"), length(which(temp_df$is_duplicated)))
    dummy_names <- str_c("Dummy",temp_df$Activity[nodes_before])

    # Add new rows to the dataframe.
    # temp_df[nrow(temp_df) + length(dummy_names),] <- NA
    temp_df[c(nrow(temp_df)+1):c(nrow(temp_df) + length(dummy_names)),] <- NA

    startrow <- nrow(temp_df) - length(dummy_names) + 1

    # Complete with data on dummy activities.
    temp_df[startrow:nrow(temp_df), 1] <- dummy_names
    temp_df[startrow:nrow(temp_df), 4] <- temp_df$Step5Col4[which(temp_df$is_duplicated)]
    temp_df[startrow:nrow(temp_df), 5] <- temp_df$Activity[which(temp_df$is_duplicated)]
    temp_df[startrow:nrow(temp_df), 6] <- nodes_before
  }

  # Create a dataframe for AOA with start and end nodes.
  # result_df <- data.frame(from = temp_df[, 4], to = temp_df[, 6], name = temp_df[, 1])
  fromto_df <- data.frame(from = temp_df[, 4], to = temp_df[, 6], name = temp_df[, 1])
  # Pass to node numbering.
  set_numberAOA(fromto_df)
}
#===============================================================================
# Creates dataframes from the input_predAOA function and user's data.
# Duration is deterministic.

merge_pred_detAOA <- function(yourdata, predecdata){
  # If dummy actions have been added, the number of rows is different.
  dummy_rows <- nrow(predecdata) - nrow(yourdata)
  if (dummy_rows == 0){
    predecdata <- data.frame(predecdata, duration = yourdata[,3])
  }else{
    tmp <- c(yourdata[,3], rep(0, dummy_rows))
    predecdata <- data.frame(predecdata, duration = tmp)
  }
  # Sort the data frame by node numbers.
  predecdata <- predecdata[with(predecdata, order(from, to)),]
}
#===============================================================================
# Creates dataframes from the input_predAOA function and user's data.
# Duration is probabilistic.

merge_pred_probAOA <- function(yourdata, predecdata){
  # If dummy actions have been added, the number of rows is different.
  dummy_rows <- nrow(predecdata) - nrow(yourdata)
  if (dummy_rows == 0){
    predecdata <- data.frame(predecdata, yourdata[,3:5])
  }else{
    dummy_mat <- matrix(0, nrow = dummy_rows, ncol = 3)
    colnames(dummy_mat) <- colnames(yourdata[,3:5])
    tmp <- rbind(yourdata[,3:5], dummy_mat)
    predecdata <- data.frame(predecdata, tmp)
  }
  # Sort the data frame by node numbers.
  predecdata <- predecdata[with(predecdata, order(from, to)),]
}
#===============================================================================
# Sets the graph node numbers.

set_numberAOA <- function(input_dt){
  # Set column names.
  colnames(input_dt) <- c("from", "to", "name")

  # Add columns for new node numbers.
  input_dt <- data.frame(input_dt, from2 = 0, to2 = 0)

  # Create a temporary graph in the DiagrammeR package.
  myedges_df <- create_edge_df(from = input_dt$from, to = input_dt$to, label = input_dt$name)
  nodes_num <- max(myedges_df$to)
  mynodes_df <- create_node_df(n = nodes_num, type = c(1:nodes_num), label = TRUE)
  tempgraph <- create_graph(nodes_df = mynodes_df, edges_df = myedges_df, directed = TRUE)

  # Initial parameter values.
  i <- 1
  node_no <- 1
  nodes_num <- count_nodes(tempgraph)

  # A loop that changes node numbers so that the start number is lower than the end number.
  while (i <= nodes_num){
    # Which nodes have no predecessors?
    nodes_selected <- tempgraph %>%
      select_nodes_by_degree(
        expressions = "indeg == 0") %>%
      get_selection()

    for (j in 1:length(nodes_selected)){
      if (nrow(input_dt[input_dt$from == nodes_selected[j],]) > 0){
        input_dt[input_dt$from == nodes_selected[j],]$from2 <- node_no
      }
      if (nrow(input_dt[input_dt$to == nodes_selected[j],]) > 0){
        input_dt[input_dt$to == nodes_selected[j],]$to2 <- node_no
      }
      node_no <- node_no + 1
    }
    # Delete used nodes.
    tempgraph <- tempgraph %>%
      select_nodes_by_degree(
        expressions = "indeg == 0") %>%
      delete_nodes_ws()
    # How many nodes have been examined?
    i <- i + length(nodes_selected)
  }
  # Create a data frame for AOA with start and end nodes.
  result_df <- data.frame(from = input_dt$from2, to = input_dt$to2, name = input_dt$name)
}

#===============================================================================
# Different approximation formulas for expected value in PERT.

PERT_mu <- function(mtd, tija, tijm, tijb){
  # Checks if mtd has the right value.
  stopifnot("There is no approximation method with the given number." = mtd >=0 & mtd <= 7)

  switch(as.character(mtd),
         # Classic formula
         "0" = (tija + 4* tijm + tijb)/6,
         # 1st and 99th percentiles
         "1" = (betainv(0.01, tija, tijm, tijb) + 4*tijm + betainv(0.99, tija, tijm, tijb))/6,
         # 5th and 95th percentiles (Moder and Rodgers, 1968)
         "2" = (betainv(0.05, tija, tijm, tijb) + 4*tijm + betainv(0.95, tija, tijm, tijb))/6,
         # 5th and 95th percentiles by (Perry and Greig, 1975)
         "3" = (betainv(0.05, tija, tijm, tijb) + 0.95*tijm + betainv(0.95, tija, tijm, tijb))/2.95,
         # Extended Pearson-Tukey (Pearson and Tukey, 1965)
         "4" = 0.185*betainv(0.05, tija, tijm, tijb) + 0.63*betainv(0.5, tija, tijm, tijb) + 0.185*betainv(0.95, tija, tijm, tijb),
         # 1st version (Golenko-Ginzburg, 1988)
         "5" = (2*tija + 9*tijm + 2*tijb)/13,
         # 2nd version (Golenko-Ginzburg, 1988)
         "6" = (3*tija + 2*tijb)/5,
         # (Farnum and Stanton, 1987)
         "7" = farnum_mu(tija, tijm, tijb)
         )
}

#===============================================================================
# Different approximation formulas for variance in PERT.

PERT_var <- function(mtd, tija, tijm, tijb){
  # Checks if mtd has the right value.
  stopifnot("There is no approximation method with the given number." = mtd >=0 & mtd <= 7)

  switch(as.character(mtd),
         # Classic formula
         "0" = (tijb - tija)^2/36,
         # 1st and 99th percentiles
         "1" = (betainv(0.99, tija, tijm, tijb) - betainv(0.01, tija, tijm, tijb))^2/36,
         # 5th and 95th percentiles (Moder and Rodgers, 1968)
         "2" = (betainv(0.95, tija, tijm, tijb) - betainv(0.05, tija, tijm, tijb))^2/10.2,
         # 5th and 95th percentiles by (Perry and Greig, 1975)
         "3" = (betainv(0.95, tija, tijm, tijb) - betainv(0.05, tija, tijm, tijb))^2/3.25^2,
         # Extended Pearson-Tukey (Pearson and Tukey, 1965)
         "4" = 0.185*betainv(0.05, tija, tijm, tijb)^2 + 0.63*betainv(0.5, tija, tijm, tijb)^2 + 0.185*betainv(0.95, tija, tijm, tijb)^2 - (0.185*betainv(0.05, tija, tijm, tijb) + 0.63*betainv(0.5, tija, tijm, tijb) + 0.185*betainv(0.95, tija, tijm, tijb))^2,
         # 1st version (Golenko-Ginzburg, 1988)
         "5" = ((tijb - tija)^2/1268)*(22 + 81*(tijm-tija)/(tijb-tija) - 81*((tijm - tija)/(tijb - tija))^2),
         # 2nd version (Golenko-Ginzburg, 1988)
         "6" = (tijb - tija)^2/25,
         # (Farnum and Stanton, 1987)
         "7" = farnum_var(tija, tijm, tijb)
  )
}
#===============================================================================
# Percentiles of the beta distribution.

betainv <- function(p, a, b, c){
  alpha <- beta <- qbeta <- NULL

  alpha <- (4*(b-a)/(c-a)) + 1
  beta <- (4*(c-b)/(c-a)) + 1
  qbeta(p, alpha, beta)*(c - a) + a
}

#===============================================================================
# Expected value by Farnum and Stanton

farnum_mu <- function(a, b, c){

  lim_left <- a + 0.13*(c - a)
  lim_right <- a + 0.87*(c - a)

  ifelse(b < lim_left, a + 2*(b - a)*(c - a)/(c - 3*a + 2*b),
         ifelse(b > lim_right, a + (c - a)^2/(3*c - a - 2*b), (a + 4*b + c)/6))
}

#===============================================================================
# Expected value by Farnum and Stanton

farnum_var <- function(a, b, c){

  lim_left <- a + 0.13*(c - a)
  lim_right <- a + 0.87*(c - a)

  ifelse(b < lim_left, ((b - a)^2)*(c - b)/(c - 2*a + b),
         ifelse(b > lim_right, (b - a)*(c - b)^2/(2*c - a - b), (c - a)^2/36))
}
