#===============================================================================
# Function that reads data in the DiagrammeR package format.
# Changing labels to character mode to display them on the graph.
# Data frame for the CPM method. TS stands for slack of time.

read_cpmAOA <- function(yourdata){
  if (ncol(yourdata) == 4){
    create_edge_df(from = yourdata[,1], to = yourdata[,2], label = yourdata[,3],
                   time = yourdata[,4], TS = rep(0, nrow(yourdata)))
  }else{
    print("The data frame for the CPM method should have 4 columns. Check the documentation.")
  }
}
#===============================================================================
# Function reads data in the DiagrammeR package format.
# Change labels to character mode to display them on the graph.
# Data frame for the PERT method. TS stands for slack of time.

read_pertAOA <- function(yourdata){
  if (ncol(yourdata) == 6){
    create_edge_df(from = as.integer(yourdata[,1]), to = as.integer(yourdata[,2]),
                   label = as.character(yourdata[,3]),
                   time = (yourdata[,4]+4*yourdata[,5]+yourdata[,6])/6,
                   timevar = (yourdata[,6]-yourdata[,4])^2/36,
                   TS = rep(0, nrow(yourdata)),)
  }else{
    print("The data frame for the PERT method should have 6 columns. Check the documentation.")
  }
}
#===============================================================================
# Creates relations from an input data frame.

make_relationsAOA <- function(yourdata, deterministic){
  if (deterministic == TRUE){
    read_cpmAOA(yourdata)
  }else{
    read_pertAOA(yourdata)
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
      Time = relations$time,
      ESij = rep(0, nrow(relations)),
      LSij = rep(0, nrow(relations)),
      EFij = rep(0, nrow(relations)),
      LFij = rep(0, nrow(relations)),
      TSij = rep(0, nrow(relations)),
      Crit = rep(c(" "), nrow(relations))
    )
  }else{
    data.frame(
      Name = relations$label,
      Time = relations$time,
      Var = relations$timevar,
      ESij = rep(0, nrow(relations)),
      LSij = rep(0, nrow(relations)),
      EFij = rep(0, nrow(relations)),
      LFij = rep(0, nrow(relations)),
      TSij = rep(0, nrow(relations)),
      Crit = rep(c(" "), nrow(relations))
    )
  }
}

#===============================================================================
# Prepares and fills a data frame that stores the schedule of the CPM method.

scheduleAOA <- function(relations, deterministic){
  ES <- LF <- TS <- NULL
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
  yourschedule$LSij <- yourschedule$LFij - yourschedule$Time
  yourschedule$EFij <- yourschedule$ESij + yourschedule$Time
  yourschedule$TSij <- yourschedule$LFij - yourschedule$ESij - yourschedule$Time
  yourschedule$Crit[which(yourschedule$TSij == 0)] <- c("*")

  # Completes TS for the edges data frame.
  yourgraph <- set_edge_attrs(yourgraph, edge_attr = TS, values = yourschedule$TSij)

  # Messages after the computation is completed.
  if (deterministic == TRUE){
    # CPM
    cat("Completion time: ", max(yourschedule$LFij), "\n")
  }else{
    # PERT
    cat("Expected compl. time distribution: N(",max(yourschedule$LFij),",", sqrt(sum(yourschedule$Var[which(yourschedule$TSij == 0)])),")\n")
  }

  # Creates a list keeping the graph and schedule.
  if (deterministic == TRUE){
    # CPM.
    list(graphAOA = yourgraph, schedule = yourschedule,
      ComplTi = max(yourschedule$LFij),
      CritAct = yourschedule$Name[which(yourschedule$TSij == 0)])
  }else{
    # PERT.
    list(graphAOA = yourgraph, schedule = yourschedule,
      ComplTi = max(yourschedule$LFij),
      SDevTi = sqrt(sum(yourschedule$Var[which(yourschedule$TSij == 0)])),
      CritAct = yourschedule$Name[which(yourschedule$TSij == 0)])
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

#' Finds a solution using CPM and PERT methods
#'
#' @param input_data Data frame containing the structure of the graph and the duration of the activity.
#'   For the CPM method there will be 4 columns (the order is important, not the name of the column):
#'   \enumerate{
#'   \item \code{from} The number of the node where the activity starts.
#'   \item \code{to} The number of the node where the activity ends.
#'   \item \code{label} Activity labels.
#'   \item \code{time} Activities duration.
#'   }
#'   For the PERT method there will be 4 columns (the order is important, not the name of the column):
#'   \enumerate{
#'   \item \code{from} The number of the node where the activity starts.
#'   \item \code{to} The number of the node where the activity ends.
#'   \item \code{label} Activity labels.
#'   \item \code{opt_time} Optimistic duration of activities.
#'   \item \code{likely_time} The most likely duration of the activity.
#'   \item \code{pes_time} Pessimistic duration of activities.
#'   }
#' @param deterministic A logical parameter specifying the solution method.
#'   If set to \code{TRUE} (default), the CPM method is used. If is set to \code{FALSE}, the PERT method is used.
#' @return The list is made of a graph, schedule and selected partial results.
#' @examples
#' x <- solve_pathAOA(cpmexample1, deterministic = TRUE)
#' y <- solve_pathAOA(pertexample1, deterministic = FALSE)
#' @import DiagrammeR
#' @import utils
#' @export
solve_pathAOA <- function(input_data, deterministic = TRUE){
  # Check if the DiagrammeR package is loaded. If not, load it.
  pckg_check("DiagrammeR")
  solved_pathAOA <- vector("list", length = 4)
  # Loads data and creates a relations data frame.
  relations <- make_relationsAOA(input_data, deterministic)
  # Creates graph and schedule and writes them into list.
  scheduleAOA(relations, deterministic)
}

#===============================================================================

# Presentation of the critical path on a graph.

#' Graph with marked critical path
#'
#' @param yourlist Data frame describing the problem
#' @param fixed_seed Optional parameter setting random seed to user value to get similar looking plots each time the function is run (set to      23 by default).
#' @return The function draws the graph along with the critical path by means of the DiagrammeR package functions.
#' @examples
#' x <- solve_pathAOA(cpmexample1, TRUE)
#' plot_crit_pathAOA(x)
#' @import DiagrammeR
#' @export
plot_crit_pathAOA <- function(yourlist, fixed_seed = 23){
  TS <- color <- time <- style <- NULL
  # Check if DiagrammeR package is loaded. If not, load it.
  pckg_check("DiagrammeR")
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
    select_edges(conditions = TS == 0) %>%
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
#' @param bar_size Thickness of the bar drawn for activity (set to 10 by default).
#' @return Draws a Gantt chart broken down into critical ("CR") and non-critical ("NC") activities.
#'   Marks the slack of time.
#' @examples
#' x <- solve_pathAOA(cpmexample1, deterministic = TRUE)
#' plot_gantt(x)
#' @import ggplot2
#' @import reshape2
#' @export
plot_gantt <- function(yourlist, bar_size = 10){
  Name <- TSij <- value <- NULL
  # Check if ggplot2 package is loaded. If not, load it.
  pckg_check("ggplot2")

  # Check if reshape2 package is loaded. If not, load it.
  pckg_check("reshape2")

  schedule <- yourlist[[2]]

  # Sorting the schedule according to slack of time.
  dftmp <- schedule[order(schedule$TSij),]

  # The TSij column gives the critical activities the symbol "CR" and the non-critical activities "NC".
  dftmp$TSij[dftmp$TSij > 0] <- c("NC")
  dftmp$TSij[dftmp$TSij == 0] <- c("CR")

  # Gantt chart
  melthar <- melt(dftmp, measure.vars = c("ESij", "EFij"))
  melthar2 <- melt(dftmp, measure.vars = c("EFij", "LFij"))

  ggplot(melthar, aes(value, Name, colour = TSij)) +
    geom_line(size = bar_size) +
    ylab(NULL) +
    xlab(NULL) +
    theme_bw() +
    theme(legend.title=element_blank())

  ggplot() +
    geom_line(melthar, mapping = aes(value, Name, colour = TSij), size = bar_size) +
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
#' @import stats
#' @export
plot_norm <- function(yourlist){
  z <- df <- x <- NULL
  # Risk-taker schedule
  hryz <- qnorm(0.3, yourlist[[3]], yourlist[[4]])
  #Belayer schedule
  hase <- qnorm(0.6, yourlist[[3]], yourlist[[4]])
  # Check if ggplot2 package is loaded. If not, load it.
  pckg_check("ggplot2")
  # Variable that holds the random values plus/minus 3 variations
  z <- seq(yourlist[[3]]-3*yourlist[[4]], yourlist[[3]]+3*yourlist[[4]], 0.01)
  df <- data.frame(x = z, y = dnorm(z))
  df <- cbind(df, pnorm(z))
  ggplot(df, aes(x, y = pnorm(z, yourlist[[3]], yourlist[[4]]))) +
    geom_line() +
    labs(y = "P(DT<=x)") +
    geom_vline(xintercept = hryz, color='blue')+
    geom_vline(xintercept = hase, color='blue')+
    theme_bw()
}

#===============================================================================
# Function that draws a graph before solving a problem.

#' Graph without critical path
#'
#' @param input_data Data frame describing the problem.
#' @param fixed_seed Optional parameter setting random seed to user value to get similar looking plots each time the function is run (set to      23 by default).
#' @return The function draws a relationship graph between activities without solving the problem
#'   and thus without marking critical activities.
#' @examples
#' plot_graphAOA(cpmexample1)
#' @import DiagrammeR
#' @export
plot_graphAOA <- function(input_data, fixed_seed = 23){
  time <- style <- NULL
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
