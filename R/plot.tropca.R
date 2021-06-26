#' Plot the Tropical Principal Components with Data Projections
#'
#' Visialize the second order tropical prinrciple components in \code{tropca}
#' as a tropical triangle with projections on a two-dimensional plot via tropical isometry.
#'
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics segments
#' @importFrom graphics par
#' @importFrom graphics legend
#'
#' @param x a fitted \code{tropca} object.
#' @param plab a vector of labels of all points in the given data matrix.
#' Not needed for unlabeled data. (default: NULL)
#' @param \dots Not used. Other arguments to plot
#'
#' @author Robert Page and Houjie Wang
#'
#' @keywords Tropical Geometry, Supervised Learning, Non-Euclidean Data
#'
#' @method plot tropca
#' @export
#' @export plot.tropca
plot.tropca <- function(x, plab = NULL, ...) {
  if (x$type == "linear space") {
    stop("Only principal components by tropical polytopes are plottable.")
  }
  object <- x
  if (is.null(plab)) plab <- as.factor(rep(1, nrow(object$projection)))
  D <- eachrow(object$pc, object$pc[1, ], "-")[-1, ]
  proj_points_plot <- do.call("rbind", lapply(lapply(1:nrow(object$projection), function(i) {
    object$projection[i, ]
  }), polytope_iso, D = object$pc))
  proj_2D_plot_m <- proj_points_plot - proj_points_plot[, 1]
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
  for (i in 1:length(unique(plab))) {
    points(x = proj_2D_plot_m[plab == unique(plab)[i], 2], y = proj_2D_plot_m[plab == unique(plab)[i], 3], pch = 16, cex = .75, col = (i + 1))
  }
  coord <- par("usr")
  legend(x = coord[2] * 1.05, y = coord[4], legend = unique(plab), pch = rep(16, 2), col = 2:(length(unique(plab)) + 1))
}
