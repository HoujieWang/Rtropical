#' Plot the Tropical Principal Components with Data Projections
#'
#' Visualize the second order tropical principle components in \code{troppca}
#' as a tropical triangle with projections on a two-dimensional plot via tropical isometry.
#'
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics segments
#' @importFrom graphics par
#' @importFrom graphics legend
#' @importFrom graphics mtext
#' @importFrom latex2exp TeX
#' @importFrom phangorn upgma
#' @importFrom phangorn RF.dist
#' @param x a fitted \code{troppca} object.
#' @param point_class a vector of labels of all points in the given data matrix.
#' Not needed for unlabeled data. (default: NULL)
#' @param fw a logical variable to determine if to add
#' Fermat-Weber point of the data projection. (default: FALSE)
#' @param plot.tree a logical value indicating if to plot tree topology (default: FALSE)
#' @param leaf.names the name of the leaves of the first tree in the data set (default: NULL)
#' @param ntopology a numeric value indicating  the number of most frequent tree topology to plot (default: 6)
#' @param \dots Not used. Other arguments to plot
#'
#' @return \code{plot.troppca} does not return anything other than the plot.
#' @method plot troppca
#' @export
#' @export plot.troppca
plot.troppca <- function(x, point_class = NULL, fw = FALSE, plot.tree = FALSE, leaf.names = NULL, ntopology = 6, ...) {
  # x = pca_fit
  if (x$type == "linear space") {
    stop("Only principal component by tropical polytope is plottable.")
  }
  D <- eachrow(x$pc, x$pc[1, ], "-")[-1, ]

  if (plot.tree){
    x$projection = t(apply(x$projection, 1, function(xx){2 - max(xx) + xx}))
    numvectors <- nrow(x$projection)
    proj_trees_matrix<-lapply(1: numvectors, function(i){
      # i = 1
      p = x$projection[i, ]
      make.matrix(p, length(leaf.names), leaf.names)
    })
    proj_trees <- lapply(proj_trees_matrix, upgma)

    type <- seq(numvectors)
    for(i in seq(numvectors)){
      if(i == type[i]) {
        for(j in seq(i, numvectors)) {
          if(RF.dist(proj_trees[[i]], proj_trees[[j]]) == 0) {
            type[j] <- i
          }
        }
      }
    }

    ntopology = min(ntopology, length(unique(type)))
    # plot_type = as.numeric(names(table(type)[order(table(type), decreasing = TRUE)[1: ntopology]]))
    plot_type = table(type)[order(table(type), decreasing = TRUE)[1: ntopology]]
    type_names = as.numeric(names(plot_type))
    x$projection = x$projection[type %in% type_names, ]
    # type = type[type %in% type_names]
    # point_class = type
    point_class = plot_type[as.character(type[type %in% type_names])]
  }

  if (fw){
    if (is.null(point_class)) point_class <- as.factor(c(rep(1, nrow(x$projection))))
    proj_points_plot <- t(apply(x$projection, 1, polytope_iso, D = x$pc))
    fw_point <- tropFW(proj_points_plot)
    proj_points_plot <- rbind(proj_points_plot, fw_point$fw)
    proj_2D_plot_m <- proj_points_plot - proj_points_plot[, 1]
  } else{
    if (is.null(point_class)) point_class <- as.factor(rep(1, nrow(x$projection)))
    proj_points_plot <- t(apply(x$projection, 1, polytope_iso, D = x$pc))
    proj_2D_plot_m <- proj_points_plot - proj_points_plot[, 1]
  }
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(xpd = TRUE, mar = c(5.1, 4.1, 4.1, 5.1))
  k <- ncol(D)
  plot(D[1, ], D[2, ], xlab = "x1", ylab = "x2")
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      tseg1 <- tropseg(D[, i], D[, j])
      tseg2 <- tropseg(D[, i], D[, j], flag = TRUE)
      if (tseg1[[2]] < tseg2[[2]]) {
        tseg <- tseg1
      } else {
        tseg <- tseg2
      }
      segments(tseg[[1]][1, 1], tseg[[1]][2, 1], tseg[[1]][1, 2], tseg[[1]][2, 2], col = "black")
      segments(tseg[[1]][1, 2], tseg[[1]][2, 2], tseg[[1]][1, 3], tseg[[1]][2, 3], col = "black")
    }
  }
  labs <- unique(point_class)
  if (plot.tree) labs = plot_type
  for (i in 1:length(labs)) {
    points(x = proj_2D_plot_m[point_class == labs[i], 2], y = proj_2D_plot_m[point_class == labs[i], 3], pch = 16, cex = .75, col = (i + 1))
  }
  coord <- par("usr")
  if (fw){
    labs <- c(labs, "FW")
    if (plot.tree) labs <- c(plot_type, "FW")
    points(x = proj_2D_plot_m[nrow(proj_2D_plot_m), 2], y = proj_2D_plot_m[nrow(proj_2D_plot_m), 3], pch = 18, cex = 2, col = "black")
    legend(x = coord[2] * 1.05, y = coord[4], legend = labs, pch = c(rep(16, length(labs)-1), 18),
           col = c(2:(length(labs)), "black"))
  }else{
    legend(x = coord[2] * 1.05, y = coord[4], legend = labs, pch = rep(16, length(labs)),
           col = c(2:(length(labs) + 1)))
  }
  if (plot.tree){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(ask = TRUE)
    par(mfrow = c(ceiling(ntopology/2), 2))
    for (i in 1: ntopology) {
      tree_idx = type_names[i]
      plot(proj_trees[[i]], "c", FALSE,
           # main=paste0("(",i, ")"),
           cex = 1.3)

      mtext(TeX("\\bullet", bold = TRUE), col = i+1)
      mtext(paste0("         ", plot_type[i]))
    }
  }
}
