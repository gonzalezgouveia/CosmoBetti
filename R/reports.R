# graphical reports and aggregate tables

#' Plot Betti Curves
#'
#' Plot Betti Curves from persistence diagram of desired dimension
#' in a ggplot style with facets and free scales
#'
#' @param betticurve is a betti matrix from bettimatrix
#' @examples
#' plot_betti(betticurve)
#' @export
plot_betti_curve <- function(betticurve){
  ggplot2::ggplot(betticurve, ggplot2::aes(x = radius, y = value, color = dim)) +
    ggplot2::geom_line() +
    ggplot2::facet_grid(rows = ggplot2::vars(dim), scales = "free")
}
