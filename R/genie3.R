#' Compute weighted adjacency matrix of inferred network
#'
#' @param data A matrix of observations of the different features. The rows must contain the observations.
#' @param K The choice of number of input genes randomly selected as candidates at each node. Must be \code{"all"} for all input features, \code{"sqrt"} for the square root of all input features (default), or an integer.
#' @param nb.trees The number of trees in ensemble for each target gene (default 1000).
#' @param regulators A set of indices or column names of entities whose observed values regulate the observed values of the targets.
#' @param targets A set of indices or column names of entities whose observed values are regulated by the regulators.
#' @param importance.measure Type of variable importance measure. Must be either \code{"IncNodePurity"} for the importance measure based on decrease of residual sum of squares, or \code{"\%IncMSE"} for importance measure optained by permutation of OOB data.
#' @param seed A random number generator seed for replication of analyses. NULL means the seed is not set.
#' @param trace Output additional information.
#' @param mc.cores The number of threads to use for parallel execution.\code{\link{genie3}}
#'
#' @references Huynh-Thu, V. A. et al. (2010) Inferring Regulatory Networks from Expression Data Using Tree-Based Methods. PLoS ONE.
#'
#' @return The weighted adjacency matrix of inferred network.
#' @export
#'
#' @importFrom randomForest randomForest
#' @importFrom parallel mclapply
#'
#' @examples
#' library(GENIE3)
#' library(ggplot2)
#' library(cowplot)
#'
#' data <- matrix(runif(100*100), ncol=100)
#' true.matrix <- matrix(sample(0:1, 20*100, replace=TRUE, prob=c(.9, .1)), ncol=10)
#' diag(true.matrix) <- 0
#' weights <- genie3(data, regulators=1:20, targets=1:100, mc.cores=8)
#' ranking <- get.ranking(weights)
#' eval <- evaluate.ranking(ranking, true.matrix=true.matrix)
#' eval$au.score
#'
#' # draw a ROC and PR curves
#' g1 <- ggplot(eval$metrics, aes(1 - spec, tpr)) + geom_step()
#' g2 <- ggplot(eval$metrics, aes(tpr, prec)) + geom_path()
#' plot_grid(g1, g2, nrow = 1)
#'
#' # Evaluate multiple rankings at the same time
#' weights.cor <- as.matrix(abs(cor(weights[,1:20], weights[,1:100])))
#' ranking.cor <- get.ranking(weights.cor)
#' rankings <- list(GENIE3=ranking, Correlation=ranking.cor)
#' eval <- evaluate.multiple.rankings(rankings, true.matrix=true.matrix)
#' eval$au.score
#'
#' # draw a ROC and PR curves
#' library(cowplot)
#' g1 <- ggplot(eval$metrics, aes(1 - spec, tpr, colour=name)) + geom_step()
#' g2 <- ggplot(eval$metrics, aes(tpr, prec, colour=name)) + geom_path()
#' plot_grid(g1, g2, nrow = 1)
genie3 <- function(data, regulators=seq_len(ncol(data)), targets=seq_len(ncol(data)),
                   K="sqrt", nb.trees=1000, importance.measure="IncNodePurity",
                   seed=NULL, trace=TRUE, mc.cores=1) {
  requireNamespace("randomForest")

  if (!(is.matrix(data) || is.data.frame(data)) || any(!apply(data, 2, is.finite))) {
    stop("Parameter \"data\" must be a matrix or a data frame consisting of finite values.")
  }

  num.samples <- nrow(data)
  num.genes <- ncol(data)
  feature.names <- colnames(data)

  # check nb.trees parameter
  if (nb.trees < 1) {
    stop("Parameter \"nb.trees\" must be larger than or equal to 0")
  }

  # check regulators parameter
  if (is.null(regulators)) {
    regulators <- seq_len(num.genes)
  } else if (class(regulators) == "character") {
    regulators <- match(regulators, feature.names)
  }
  if ((class(regulators) != "numeric" && class(regulators) != "integer") ||
      !all(regulators %in% seq_len(num.genes)) ||
      length(unique(regulators)) != length(regulators)) {
    stop("Parameter \"regulators\" must either be NULL, or be a subset of either colnames(data) or seq_len(ncol(data)), and may not contain any repeated elements.")
  }
  num.regulators <- length(regulators)
  regulator.names <- feature.names[regulators]

  # check targets parameter
  if (is.null(targets)) {
    targets <- seq_len(num.genes)
  } else if (class(targets) == "character") {
    targets <- match(targets, feature.names)
  }
  if ((class(targets) != "numeric" && class(targets) != "integer") ||
      !all(targets %in% seq_len(num.genes)) ||
      length(unique(targets)) != length(targets)) {
    stop("Parameter \"targets\" must either be NULL, or be a subset of either colnames(data) or seq_len(ncol(data)), and may not contain any repeated elements.")
  }
  num.targets <- length(targets)
  target.names <- feature.names[targets]

  # check K parameter
  if (class(K) == "numeric") {
    mtry <- K
  } else if (K == "sqrt") {
    mtry <- round(sqrt(num.regulators))
  } else if (K == "all") {
    mtry <- num.regulators-1
  } else {
    stop("Parameter \"K\" must be \"sqrt\", or \"all\", or an integer")
  }

  # check importance.measure parameter
  if (importance.measure != "IncNodePurity" && importance.measure != "%IncMSE") {
    stop("Parameter \"importance.measure\" must be \"IncNodePurity\" or \"%IncMSE\"")
  }

  if (importance.measure == "IncNodePurity") {
    importance.parameter <- F
  } else if (importance.measure == "%IncMSE") {
    importance.parameter <- T
  }

  # check seed parameter, use the current seed if none is given.
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # check mc.cores parameter
  if (mc.cores == "qsub") {
    requireNamespace("PRISM")
    qsub.environment <- c("data", "feature.names", "regulators", "nb.trees", "importance.measure", "regulator.names", "num.regulators", "importance.parameter")
    lapplyfun <- function(X, FUN, ...) PRISM::qsub.lapply(X = X, FUN = FUN, qsub.environment = qsub.environment)
  } else if (mc.cores < 1) {
    stop("Parameter \"mc.cores\" must be larger than or equal to 0")
  } else if (mc.cores > 1) {
    requireNamespace("parallel")
    lapplyfun <- function(...) parallel::mclapply(..., mc.cores=mc.cores)
  } else {
    lapplyfun <- lapply
  }

  if (trace) {
    cat("GENIE3 parameter checks are OK!\n")
    cat("Starting GENIE3 computations.\n")
    flush.console()
  }

  # scale the data matrix
  data <- scale(data)

  # compute importances for every target
  output <- do.call("cbind", lapplyfun(targets, function(target.index) {
    if (trace) {
      cat("Computing for target ", ifelse(is.null(feature.names), target.index, feature.names[[target.index]]), "\n", sep="")
      flush.console()
    }
    x <- data[,setdiff(regulators, target.index),drop=F]
    y <- data[,target.index]
    rf <- randomForest::randomForest(x, y, mtry = mtry, ntree = nb.trees, importance = importance.parameter)
    im <- rf$importance[,importance.measure]

    out <- setNames(rep(0, num.regulators), regulator.names)
    out[target.index != regulators] <- im
    out
  }))

  rownames(output) <- regulator.names
  colnames(output) <- target.names

  output / num.samples
}

#' Make an ordered data frame of regulatory links
#'
#' @param weights A weighted adjacency matrix as returned by \code{\link{genie3}}
#' @param max.links The maximum number of links to output.
#'
#' @return A data frame of links.
#' @import reshape2 dplyr
#' @export
#'
#' @seealso \code{\link{genie3}}
get.ranking <- function(weights, max.links = 100000) {
  requireNamespace("reshape2")
  requireNamespace("dplyr")

  # This following line is added just to appease R check. When looking at this code, pretend this line doesn't exist.
  regulator <- target <- value <- NULL

  ranking <- reshape2::melt(weights, varnames=c("regulator", "target"))                # transform into an expanded data frame
  ranking <- dplyr::filter(ranking, as.character(regulator) != as.character(target))   # remove self loops
  ranking <- dplyr::arrange(ranking, dplyr::desc(value))                               # order by value

  if (is.numeric(max.links) && nrow(ranking) > max.links) {
    ranking <- dplyr::top_n(ranking, max.links, value)                                 # take top max.links edges
  }

  ranking
}

#' Evaluate a network ranking
#'
#' @param ranking the data frame as returned by \code{\link{get.ranking}}. This data frame must contain at least 2 columns, called \code{regulator} and \code{target}.
#' @param true.matrix a matrix with 0's and 1's, representing the golden standard. The rownames and colnames must be the same as the names used in the regulator and target columns in \code{ranking}.
#' @param num.extend.steps The number of steps with which to fill the ranking as if random, if only a part of the ranking is given
#'
#' @return a list containing 2 items, the ranked evaluation and the area under the curve scores
#' @import ROCR pracma dplyr
#' @export
#'
#' @seealso \code{\link{genie3}}
evaluate.ranking <- function(ranking, true.matrix, num.extend.steps = 10000) {
  requireNamespace("ROCR")
  requireNamespace("pracma")

  # Check whether columns are already numeric
  if (!is.numeric(ranking$regulator)) {
    regulator <- match(as.character(ranking$regulator), rownames(true.matrix))
  } else {
    regulator <- ranking$regulator
  }

  if (!is.numeric(ranking$target)) {
    target <- match(as.character(ranking$target), colnames(true.matrix))
  } else {
    target <- ranking$target
  }

  # Retrieve true positives
  tp <- ifelse(!is.na(regulator) & !is.na(target), true.matrix[cbind(regulator, target)], 0)

  # Retrieve value
  if ("value" %in% colnames(ranking)) {
    value <- ranking$value
  } else {
    value <- dplyr::percent_rank(-seq_along(tp))
  }

  num.positives <- sum(true.matrix, na.rm = T)
  num.predictions <- sum(!is.na(true.matrix))

  evaluate.ranking.direct(tp, num.positives, num.predictions, num.extend.steps)
}

#' Evaluate a ranking
#'
#' @param is.true A vector with 0's and 1's, representing the golden standard.
#' @param num.positives The total number of positives
#' @param num.predictions The total number ranked values
#' @param num.extend.steps The number of steps with which to fill the ranking as if random, if only a part of the ranking is given
#'
#' @return a list containing two items, the ranked evaluation and the area under the curve scores
#' @export
evaluate.ranking.direct <- function(is.true, num.positives, num.predictions, num.extend.steps = 10000) {
  requireNamespace("dplyr")
  requireNamespace("pracma")

  # calculate base statistics
  num.selected <- seq_along(is.true)
  tp <- cumsum(is.true)
  fp <- num.selected - tp
  num.predicted <- length(tp)
  num.negatives <- num.predictions - num.positives

  # extend base statistics, if necessary
  if (num.extend.steps > 0 && num.predicted != num.predictions) {
    diff.predictions <- num.predictions - num.predicted
    diff.trues <- num.positives - tail(tp, 1)
    diff.negs <- num.negatives - tail(fp, 1)

    multiplier <- seq_len(num.extend.steps) / num.extend.steps

    extra.num.selected <- multiplier * diff.predictions + tail(num.selected, 1)
    extra.tp <- multiplier * diff.trues + tail(tp, 1)
    extra.fp <- multiplier * diff.negs + tail(fp, 1)

    num.selected <- c(num.selected, extra.num.selected)
    is.true <- c(is.true, rep(NA, num.steps))
    tp <- c(tp, extra.tp)
    fp <- c(fp, extra.fp)
  }

  # calculate extended statistics
  metrics <- dplyr::data_frame(
    num.selected = c(0, num.selected),
    is.true = c(NA, is.true),
    tp = c(0, tp),
    fp = c(0, fp),
    fn = num.positives - tp,
    tn = num.negatives - fp,
    acc = (tp + tn) / (num.positives + num.negatives),
    tpr = tp / num.positives,
    spec = tn / num.negatives,
    prec = ifelse(num.selected == 0, 1, tp / (tp + fp)),
    npv = tn / (tn + fn),
    f1 = 2 * tp / (2 * tp + fp + fn),
    mcc = ifelse(num.selected == 0, 0, (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))),
    informedness = tpr + spec - 1,
    markedness = prec + npv - 1
  )

  # calculate area under the curves
  area.under <- dplyr::data_frame(
    auroc = pracma::trapz(1 - metrics$spec, metrics$tpr),
    aupr = abs(pracma::trapz(metrics$tpr, metrics$prec)),
    F1 = ifelse(auroc + aupr != 0, 2 * auroc * aupr / (auroc + aupr), 0)
  )

  # generate output
  list(metrics = metrics, area.under = area.under)
}

#' Evaluate and compare multiple rankings
#'
#' @param rankings a list of rankings as preduced by \code{\link{get.ranking}}. See \code{\link{evaluate.ranking}} for more information.
#' @param true.matrix a matrix with 0's and 1's, representing the golden standard. The rownames and colnames must me the same as the names used in the regulator and target columns in \code{ranking}.
#' @param num.extend.steps The number of steps with which to fill the ranking as if random, if only a part of the ranking is given
#'
#' @return a list containing 2 items, the ranked evaluation and the area under the curve scores
#' @import dplyr
#' @export
#'
#' @seealso \code{\link{genie3}}
evaluate.multiple.rankings <- function(rankings, true.matrix, num.extend.steps = 10000) {
  requireNamespace("dplyr")

  if (is.null(names(rankings))) {
    ranking.names <- seq_along(rankings)
  } else {
    ranking.names <- names(rankings)
  }

  evals <- lapply(rankings, evaluate.ranking, true.matrix = true.matrix, num.extend.steps = num.extend.steps)

  metrics <- dplyr::bind_rows(lapply(ranking.names, function(rn) {
    data.frame(name = rn, evals[[rn]]$metrics, check.names = F, stringsAsFactors = F)
  }))

  area.under <- dplyr::bind_rows(lapply(ranking.names, function(rn) {
    data.frame(name = rn, evals[[rn]]$area.under, check.names = F, stringsAsFactors = F)
  }))

  list(metrics = metrics, area.under = area.under)
}
