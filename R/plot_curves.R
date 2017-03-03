#' Make a plot of the ROC and PR curves.
#'
#' @return The evaluation AU plots.
#' @export
#'
#' @importFrom cowplot plot_grid
#' @import ggplot2
#'
#' @seealso \code{\link{run_genie3}}
plot_curves <- function(eval_object) {
  metrics <- eval_object$metrics

  if (!"name" %in% colnames(metrics)) {
    aes_1 <- aes(1 - spec, tpr)
    aes_2 <- aes(tpr, spec)
  } else {
    aes_1 <- aes(1 - spec, tpr, colour = name)
    aes_2 <- aes(tpr, spec, colour = name)
  }
  g1 <- ggplot(metrics, aes_1) + geom_step() + labs(x = "False positive rate", y = "Precision")
  g2 <- ggplot(metrics, aes_2) + geom_path() + labs(x = "Precision", y = "Recall")
  cowplot::plot_grid(g1, g2, nrow = 1)
}
