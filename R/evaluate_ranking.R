#' Evaluate a network ranking
#'
#' @param ranking The data frame as returned by \code{\link{run_genie3}}.
#'   This data frame must consist of three columns:
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
#'
#' @import dplyr
#' @export
#'
#' @seealso \code{\link{run_genie3}}
evaluate_ranking <- function(ranking, true_interactions, regulators, targets, extend_by = 10000) {
  requireNamespace("dplyr")

  # join ranking with gold
  true_interactions1 <- true_interactions %>% mutate(true = T)
  ranking_joined <-
    ranking %>%
    left_join(true_interactions1, by = c("regulator", "target")) %>%
    mutate(true = !is.na(true))

  num_positive_interactions <- nrow(true_interactions)
  num_possible_interactions <- length(regulators) * length(targets) - sum(regulators %in% targets)

  evaluate_ranking_direct(
    values = ranking_joined$importance,
    are_true = ranking_joined$true,
    num_positive_interactions = num_positive_interactions,
    num_possible_interactions = num_possible_interactions,
    extend_by = extend_by)
}
