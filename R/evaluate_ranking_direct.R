#' Evaluate a ranking
#'
#' @param values A vector of importance values of predicted interactions.
#' @param are_true A vector denoting whether the corresponding predicted interactions are true.
#' @param num_positive_interactions The total number of positives.
#' @param num_possible_interactions The total number ranked values.
#' @param extend_by The number of steps with which to fill the ranking as if random, if only a part of the ranking is given.
#'
#' @return A list containing two items, the ranked evaluation and the area under the curve scores
#'
#' @import dplyr
#' @importFrom pracma trapz
#'
#' @export
#'
#' @seealso \code{\link{run_genie3}}
evaluate_ranking_direct <- function(values, are_true, num_positive_interactions,
                                    num_possible_interactions, extend_by = 10000) {
  requireNamespace("dplyr")
  requireNamespace("pracma")

  ord <- order(rank(values, ties.method = "random"))
  values <- values[ord]
  are_true <- are_true[ord]

  # calculate base statistics
  num_selected <- seq_along(are_true)
  tp <- cumsum(are_true)
  fp <- num_selected - tp
  length_ranking <- length(tp)
  num_negative_interactions <- num_possible_interactions - num_positive_interactions

  # extend base statistics, if necessary
  if (extend_by > 0 && length_ranking != num_possible_interactions) {
    diff.predictions <- num_possible_interactions - length_ranking
    diff.trues <- num_positive_interactions - tail(tp, 1)
    diff.negs <- num_negative_interactions - tail(fp, 1)

    multiplier <- seq_len(extend_by) / extend_by

    extra_num_selected <- multiplier * diff.predictions + tail(num_selected, 1)
    extra_tp <- multiplier * diff.trues + tail(tp, 1)
    extra_fp <- multiplier * diff.negs + tail(fp, 1)

    num_selected <- c(num_selected, extra_num_selected)
    are_true <- c(are_true, rep(NA, extend_by))
    tp <- c(tp, extra_tp)
    fp <- c(fp, extra_fp)
  }

  # calculate extended statistics
  metrics <- dplyr::data_frame(
    num_selected = c(0, num_selected),
    are_true = c(NA, are_true),
    tp = c(0, tp),
    fp = c(0, fp),
    fn = num_positive_interactions - tp,
    tn = num_negative_interactions - fp,
    acc = (tp + tn) / (num_positive_interactions + num_negative_interactions),
    tpr = tp / num_positive_interactions,
    spec = tn / num_negative_interactions,
    prec = ifelse(num_selected == 0, 1, tp / (tp + fp)),
    npv = tn / (tn + fn),
    f1 = 2 * tp / (2 * tp + fp + fn),
    mcc = ifelse(num_selected == 0, 0, (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))),
    informedness = tpr + spec - 1,
    markedness = prec + npv - 1
  )

  # calculate area under the curves
  area_under <- dplyr::data_frame(
    auroc = pracma::trapz(1 - metrics$spec, metrics$tpr),
    aupr = abs(pracma::trapz(metrics$tpr, metrics$prec)),
    F1 = ifelse(auroc + aupr != 0, 2 * auroc * aupr / (auroc + aupr), 0)
  )

  # generate output
  list(metrics = metrics, area_under = area_under)
}
