edge2HPD_FIX <- function (edge_df = NULL, axis.cols = NULL, type = "2D", desc = NULL, 
          ...) 
{
  require(fastmatch)
  if (is.null(edge_df)) {
    stop("No edge data provided")
  }
  if (!is.data.frame(edge_df)) {
    stop("edge_df is not a data frame")
  }
  lab1 <- unlist(edge_df[, 1])
  lab1 <- as.character(lab1)
  lab2 <- unlist(edge_df[, 2])
  lab2 <- as.character(lab2)
  nn <- length(unique(c(lab1, lab2)))
  size <- rep(1, nn)
  id <- 1:nn
  axis <- rep(1, nn)
  color <- as.character(rep("black", nn))
  radius <- rep(1, nn)
  HPD <- list()
  HPD$nodes$id <- id
  HPD$nodes$lab <- unique(c(lab1, lab2))
  HPD$nodes$axis <- axis
  HPD$nodes$radius <- radius
  HPD$nodes$size <- size
  HPD$nodes$color <- color
  ne <- nrow(edge_df)
  edge_df[, 1] <- as.character(edge_df[, 1])
  edge_df[, 2] <- as.character(edge_df[, 2])
  HPD$edges$id1 <- rep(NA, ne)
  HPD$edges$id2 <- rep(NA, ne)
  
  HPD$edges$id1 <- fmatch(edge_df[, 1], HPD$nodes$lab)
  HPD$edges$id2 <- fmatch(edge_df[, 2], HPD$nodes$lab)
  
  # for (n in 1:ne) {
  #   pat1 <- paste("\\b", edge_df[n, 1], "\\b", sep = "")
  #   pat2 <- paste("\\b", edge_df[n, 2], "\\b", sep = "")
  #   HPD$edges$id1[n] <- grep(pat1, HPD$nodes$lab)
  #   HPD$edges$id2[n] <- grep(pat2, HPD$nodes$lab)
  # }
  if (ncol(edge_df) > 2) {
    if (is.numeric(edge_df[, 3]) | is.integer(edge_df[, 3])) {
      edge_weight <- edge_df[, 3]
    }
    else {
      warning("No edge weight column detected. Setting default edge weight to 1")
      edge_weight <- rep(1, ne)
    }
  }
  HPD$edges$weight <- edge_weight
  HPD$edges$color <- rep("gray", ne)
  HPD$nodes <- as.data.frame(HPD$nodes)
  HPD$edges <- as.data.frame(HPD$edges)
  if (is.null(desc)) {
    desc <- "No description provided"
  }
  HPD$desc <- desc
  if (is.null(axis.cols)) {
    axis.cols <- brewer.pal(length(unique(HPD$nodes$axis)), 
                            "Set1")
  }
  HPD$axis.cols <- axis.cols
  HPD$nodes$axis <- as.integer(HPD$nodes$axis)
  HPD$nodes$size <- as.numeric(HPD$nodes$size)
  HPD$nodes$color <- as.character(HPD$nodes$color)
  HPD$nodes$lab <- as.character(HPD$nodes$lab)
  HPD$nodes$id <- as.integer(HPD$nodes$id)
  HPD$edges$id1 <- as.integer(HPD$edges$id1)
  HPD$edges$id2 <- as.integer(HPD$edges$id2)
  HPD$edges$weight <- as.numeric(HPD$edges$weight)
  HPD$edges$color <- as.character(HPD$edges$color)
  HPD$type <- type
  class(HPD) <- "HivePlotData"
  chkHPD(HPD)
  return(HPD)
}
