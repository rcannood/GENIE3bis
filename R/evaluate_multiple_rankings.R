#' Evaluate and compare multiple rankings
#'
#' @param rankings A list of data frames as returned by \code{\link{run_genie3}}.
#'   Each data frame must consist of three columns:
#'   \code{regulator} (factor or character),
#'   \code{target} (factor or character), and
#'   \code{importance} (numeric).
#' @param true_interactions A data frame containing the true interactionss.
#'   This data frame must consist of two columns:
#'   \code{regulator} (factor or character),
#'   \code{target} (factor or character).
#' @param extend_by The number of steps with which to fill the ranking
#'   as if random, if only a part of the ranking is given.
#'
#' @return a list containing 2 items, the ranked evaluation and the area
#'   under the curve scores.
#' @import dplyr
#' @importFrom tibble add_column
#' @export
#'
#' @seealso \code{\link{run_genie3}}
evaluate_multiple_rankings <- function(rankings, true_interactions, regulators, targets, extend_by = 10000) {
  requireNamespace("dplyr")

  if (is.null(names(rankings))) {
    ranking_names <- seq_along(rankings)
  } else {
    ranking_names <- names(rankings)
  }

  evals <- lapply(
    X = rankings,
    FUN = GENIE3::evaluate_ranking,
    true_interactions = true_interactions,
    regulators = regulators,
    targets = targets,
    extend_by = extend_by)

  metrics <- dplyr::bind_rows(lapply(ranking_names, function(rn) {
    df <- evals[[rn]]$metrics
    df %>% tibble::add_column(name = rep(rn, nrow(df)), .before = 1)
  }))

  area_under <- dplyr::bind_rows(lapply(ranking_names, function(rn) {
    df <- evals[[rn]]$area_under
    df %>% tibble::add_column(name = rep(rn, nrow(df)), .before = 1)
  }))

  list(metrics = metrics, area_under = area_under)
}
