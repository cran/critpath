# Function that reads data in the DiagrammeR package format.
# Changing labels to character mode to display them on the graph.
# Data frame for the CPM method.

read_lessAOA <- function(yourdata){
  if (ncol(yourdata) == 7){
    create_edge_df(from = as.integer(yourdata[,1]), to = as.integer(yourdata[,2]),
                   label = as.character(yourdata[,3]),
                   time = yourdata[,4], bound_time = yourdata[,5],
                   TS = rep(0, nrow(yourdata)),
                   accel_cost = (yourdata[,7]-yourdata[,6])/(yourdata[,4]-yourdata[,5]))
  }else{
    print("The data frame for the LESS method should have 7 columns. Check the documentation.")
  }
}
#===============================================================================
# Function that creates a data frame to hold the vertices. Initially, they are all zero.

make_nodeslessAOA <- function(input_tab){
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
# Function that completes the graph in the current iteration

fill_grapfAOA <- function(actgraph, frame){
  ES <- LF <- NULL

  # Calculates ES values for vertices. The loop goes through SINGLE activities, from the first to the last one.
  # Compares the ES value at a given vertex with the sum of the ES value at the predecessor vertex
  # and the activity duration time.

  for(i in  1:c(nrow(frame))){
    actgraph <- set_node_attrs(actgraph, node_attr = ES,
                              values = max(get_node_attrs(actgraph,
                                                          node_attr = ES,
                                                          nodes = c(frame$to[i])),
                                           get_node_attrs(actgraph,
                                                          node_attr = ES,
                                                          nodes = c(frame$from[i])) + c(frame$time[i])),
                              nodes = c(frame$to[i])
                              )
  }

  # In all vertices, we change the value of the LF attribute to ES from the last vertex.
  actgraph <- set_node_attrs(actgraph, node_attr = LF,
                            values = get_node_attrs(actgraph,
                                                    node_attr = ES,
                                                    nodes = max(frame$to))
                            )
  # Calcultaes LF. The loop goes through SINGLE activities, from the last to the first one.
  # Compares the LF value in the previous vertex with the difference between
  # the LF value of the current vertex and the duration of the activity. In addition, we complete the schedule.

  for(i in  c(nrow(frame)):1){
    actgraph <- set_node_attrs(actgraph, node_attr = LF,
                              values = min(get_node_attrs(actgraph,
                                                          node_attr = LF,
                                                          nodes = c(frame$from[i])),
                                           get_node_attrs(actgraph, node_attr = LF,
                                                          nodes = c(frame$to[i])) - c(frame$time[i])),
                              nodes = c(frame$from[i])
                              )
  }
  return(actgraph)
}
#===============================================================================
# Function completes the TS value for the activities from the graph data frame.

fill_slackAOA <- function(actgraph, frame){
  TS <- LF <- ES <- time <- NULL
  for(i in 1:c(nrow(frame))){
    new_value <- get_node_attrs(actgraph, node_attr = LF,
                              nodes = frame$to[i]) - get_node_attrs(actgraph,
                                                                    node_attr = ES,
                                                                    nodes = frame$from[i]) - get_edge_attrs(actgraph,
                                                                                                            edge_attr = time,
                                                                                                            from = frame$from[i], to = frame$to[i])

    actgraph <- set_edge_attrs(actgraph, edge_attr = TS, values = new_value,
                              from = frame$from[i], to = frame$to[i])
  }
  return(actgraph)
}
#===============================================================================
# Function that calculates the direct cost
calc_DC <- function(input_cost){
  DC <- DC + input_cost
}
#===============================================================================
# Funcyion that calculates the indirect cost
calc_IC <- function(ICconst, ICslope, tfinal){
  ICconst + ICslope * tfinal
}
#===============================================================================
# Solving the LESS method

#' Determines the solution using the LESS method
#'
#' @param input_data Data frame containing the graph structure and activity durations.
#'   The frame consists of 7 columns (the order matters):
#'   \enumerate{
#'   \item \code{from} The number of the node where the activity starts.
#'   \item \code{to} The number of the node where the activity ends.
#'   \item \code{label} Activity labels.
#'   \item \code{time} Normal duration of the activity.
#'   \item \code{bound_time} Boundary (the shortest possible) duration of activities.
#'   \item \code{norm_cost} Normal costs.
#'   \item \code{bound_cost} Boundary costs.
#'   }
#' @param ICconst Intercept of the indirect cost function.
#' @param ICslope Slope of the indirect cost function.
#' @return A list made of a graph and a result set.
#' @examples
#' z <-  solve_lessAOA(lessexample1, 50, 15)
#' @import DiagrammeR
#' @export
solve_lessAOA <- function(input_data, ICconst, ICslope){
  LF <- ES <- TS <- nodes_num <- accel_cost <- time <- bound_time <- label <- NULL
# Reading data and creating a graph
  relations <- read_lessAOA(input_data)
  vertices <- make_nodeslessAOA(relations)
# the number of vertices in the graph
  nodes_num <- max(relations$to)

  # Constructing a graph from vertices and relations.
  yourgraph <- create_graph(nodes_df = vertices, edges_df = relations,
                       directed = TRUE)

  # Graph for normal times
  yourgraph <- fill_grapfAOA(yourgraph, relations)
  yourgraph <- fill_slackAOA(yourgraph, relations)

  # total cost for tn
  DC <- sum(input_data[,6])
  TC <- DC + ICconst + ICslope * as.numeric(get_node_attrs(yourgraph,
                                                           LF,
                                                           nodes = nodes_num))

  # starting iteration number
  iter <- 1
  repeat{
    # Save the graph with the best results
    best_graph <- yourgraph

    # Create a subgraph of critical activities.
    crit_graph <- yourgraph %>%
      select_edges(conditions = TS == 0) %>%
      transform_to_subgraph_ws() %>%
      clear_selection()

    # Make a list of possible paths for the subgraph.
    path_list <- get_paths(crit_graph, from = 1,  to = nodes_num)
    accel_buf <- rep(0,length(path_list))

    for (j in 1:length(path_list)) {
      # select a graph created from one of the critical paths
      graphtmp <- crit_graph %>%
        select_nodes_by_id(nodes = path_list[[j]]) %>%
        transform_to_subgraph_ws() %>%
        clear_selection()

      # Modification of the subgraph so that it consists of only the edges with sij> 0 and time> bound_time
      graphtmp <- graphtmp %>%
        select_edges(conditions = accel_cost > 0) %>%
        transform_to_subgraph_ws() %>%
        clear_selection()

      graphtmp <- graphtmp %>%
        select_edges(conditions = time > bound_time) %>%
        transform_to_subgraph_ws() %>%
        clear_selection()

      # remember the id of the edge with the lowest acceleration cost for a given path

      accel_buf[j] <- graphtmp %>%
        get_edge_ids(conditions = accel_cost == min(get_edge_attrs(graphtmp, accel_cost)))
    }

    # Create a vector of unique edges to accelerate
    accel_act <- unique(accel_buf)

    # Change the durations for the relations data frame
    relations$time[accel_act] <- relations$time[accel_act] - 1

    # Creating the graph from scratch
    yourgraph <- create_graph(nodes_df = vertices, edges_df = relations,
                         directed = TRUE)

    # Complete the main graph and determine TS
    yourgraph <- fill_grapfAOA(yourgraph, relations)
    yourgraph <- fill_slackAOA(yourgraph, relations)

    # acceleration cost for selected edges
    buffer<- yourgraph %>%
      select_edges_by_edge_id(edges = accel_act) %>%
      get_edge_attrs_ws(edge_attr = accel_cost)
    buffer <- sum(as.numeric(buffer))

    iter <- iter + 1
    # New total cost
    DC <- DC + buffer
    TC[iter] <- DC + ICconst + ICslope * as.numeric(get_node_attrs(yourgraph, LF,
                                                                   nodes = nodes_num))

    # check the loop stop condition
    if(TC[iter] > TC[iter-1]){
      break
    }
  }
  summ_graph <- data.frame(get_edge_info(best_graph)[,-4],
                            name = get_edge_attrs(best_graph, edge_attr = label),
                            time = get_edge_attrs(best_graph, edge_attr = time),
                            bound.time = get_edge_attrs(best_graph, edge_attr = bound_time),
                            TS = get_edge_attrs(best_graph, edge_attr = TS),
                            accel.cost = get_edge_attrs(best_graph, edge_attr = accel_cost)
  )

  cat("Reduced completion time: ", get_node_attrs(best_graph, node_attr = LF, nodes = nodes_num), "\n")
  min_time <- get_node_attrs(best_graph, node_attr = LF, nodes = nodes_num)

  list(graphAOA = best_graph, summary_less = summ_graph,
       critical = summ_graph$name[summ_graph$TS == 0],
       TC_iter = TC,
       min_cost = TC[length(TC)-1],
       min_time = unname(min_time))
}
#===============================================================================
# Function that plots total cost values

#' Total cost change plot
#'
#' @param your_list List containing solved problem
#' @return Based on the results of the LESS method, a graph of the total cost value of all iterations is created
#' @examples
#' z <- solve_lessAOA(lessexample1, 50, 15)
#' plot_TC(z)
#' @import ggplot2
#' @export
plot_TC <- function(your_list){
  iter <- NULL
  TC <- your_list[[4]]
  iter_res <- data.frame(iter = c(1:length(TC)), TC)
  ggplot(iter_res, aes(x = iter, y = TC)) + geom_point() + geom_line() +
    geom_text(aes(label = TC), vjust = -1)
}

