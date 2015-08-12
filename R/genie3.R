#' Compute weighted adjacency matrix of inferred network
#'
#' @param data A matrix of observations of the different features. The rows must contain the observations.
#' @param K The choice of number of input genes randomly selected as candidates at each node. Must be "all" for all input features, "sqrt" for the square root of all input features (default), or an integer.
#' @param nb.trees The number of trees in ensemble for each target gene (default 1000).
#' @param regulators A set of indices or column names of entities whose observed values regulate the observed values of the targets.
#' @param targets A set of indices or column names of entities whose observed values are regulated by the regulators.
#' @param importance.measure Type of variable importance measure. Must be either "IncNodePurity" for the importance measure based on decrease of residual sum of squares, or "%IncMSE" for importance measure optained by permutation of OOB data."
#' @param seed A random number generator seed for replication of analyses. NULL means the seed is not set.
#' @param trace Output additional information.
#' @param mc.cores The number of threads to use for parallel execution.
#'
#' @return The weighted adjacency matrix of inferred network.
#' @export
#'
#' @examples
#' data <- matrix(runif(100*100), ncol=100)
#' true.matrix <- matrix(sample(0:1, 20*100, replace=T, prob=c(.9, .1)), ncol=10)
#' diag(true.matrix) <- 0
#' weights <- genie3(data, regulators=1:20, targets=1:100, mc.cores=8)
#' ranking <- get.ranking(weights)
#' evaluation <- evaluate.ranking(ranking, true.matrix=true.matrix)
#' evaluation$au.score
#' 
#' # ROC curve
#' library(ggplot2)
#' ggplot(evaluation$metrics, aes(fpr, rec)) + geom_line() + coord_cartesian(xlim = c(0, 1), ylim=c(0, 1))
#' 
#' # PR curve
#' ggplot(evaluation$metrics, aes(rec, prec)) + geom_line() + coord_cartesian(xlim = c(0, 1), ylim=c(0, 1))
genie3 <- function(data, K="sqrt", nb.trees=1000, regulators=seq_len(ncol(data)), targets=seq_len(ncol(data)), importance.measure="IncNodePurity", seed=NULL, trace=TRUE, mc.cores=1) {
  require(randomForest)
  
  if (mc.cores > 1) {
    require(parallel)
    lapplyfun <- function(...) mclapply(..., mc.cores=mc.cores)
  } else {
    lapplyfun <- lapply
  }
  
  # set random number generator seed if seed is given
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (importance.measure != "IncNodePurity" && importance.measure != "%IncMSE") {
    stop("Parameter importance.measure must be \"IncNodePurity\" or \"%IncMSE\"")
  }
  
  num.samples <- nrow(data)
  num.genes <- ncol(data)
  feature.names <- colnames(data)
  
  # check regulators parameter
  if (is.null(regulators)) {
    regulators <- seq_len(num.genes)
  } else if (class(regulators) == "character") {
    regulators <- match(regulators, colnames(data))
  }
  num.regulators <- length(regulators)
  regulator.names <- feature.names[regulators]
  
  # check targets parameter
  if (is.null(targets)) {
    targets <- seq_len(num.genes)
  } else if (class(targets) == "character") {
    targets <- match(targets, colnames(data))
  }
  num.targets <- length(targets)
  target.names <- feature.names[targets]
  
  # check mtry parameter
  if (class(K) == "numeric") {
    mtry <- K
  } else if (K == "sqrt") {
    mtry <- round(sqrt(num.regulators))
  } else if (K == "all") {
    mtry <- num.regulators-1
  } else {
    stop("Parameter K must be \"sqrt\", or \"all\", or an integer")
  }
  
  if (trace) {
    cat("Starting RF computations with ", nb.trees, " trees/target, and ", mtry, " candidate input genes/tree node\n", sep="")
    flush.console()
  }
  
  # normalise expression matrix
  data <- scale(data)
  
  # compute importances for every target 
  execute <- function(target.index) {
    if (trace) {
      cat("Computing for target ", ifelse(is.null(feature.names), target.index, feature.names[[target.index]]), "\n", sep="")
      flush.console()
    }
    x <- data[,setdiff(regulators, target.index),drop=F]
    y <- data[,target.index]
    rf <- randomForest(x, y, mtry=mtry, ntree=nb.trees, importance=TRUE)
    im <- rf$importance[,importance.measure]
    
    out <- setNames(rep(0, num.regulators), regulator.names)
    out[target.index != regulators] <- im
    out
  }
  
  output <- do.call("cbind", lapplyfun(targets, execute))
  
  rownames(output) <- regulator.names
  colnames(output) <- target.names
  
  output / num.samples
}

#' Make an ordered data frame of regulatory links
#'
#' @param matrix A weighted adjacency matrix as returned by genie3.
#' @param max.links The maximum number of links to output.
#'
#' @return A data frame of links.
#' @export
get.ranking <- function(weights, max.links=100000) {
  require(reshape2)
  require(dplyr)
  ranking <- 
    weights %>% 
    melt(varnames=c("regulator", "target")) %>%                     # transform into a data frame
    filter(as.character(regulator) != as.character(target)) %>%     # remove self loops
    arrange(desc(value))                                            # order by value
  
  if (is.integer(max.links) && nrow(ranking) > max.links) {
    ranking <- ranking %>% top_n(value, max.links)                  # take top max.links edges
  }
    
  ranking
}

#' Title Evaluate a network ranking
#'
#' @param ranking the data frame as returned by 'get.ranking'. This data frame must contain at least 2 columns, called 'regulator' and 'target'.
#' @param true.matrix a matrix with 0's and 1's, representing the golden standard. The rownames and colnames must me the same as the names used in the regulator and target columns in 'ranking'.
#' @param method 
#' @param perf.measures the performance measured in ROCR. Must at least contain 'fpr', 'rec', 'spec', and 'prec'.
#'
#' @return a list containing 2 items, the ranked evaluation and the area under the curve scores
#' @export
evaluate.ranking <- function(ranking, true.matrix, perf.measures=c("acc", "rec", "prec", "fpr", "spec", "phi", "f")) {
  require(ROCR)
  require(pracma)
  require(dplyr)
  
  # get TPs along ranking
  reg.names <- if (class(ranking$regulator) == "integer") seq_len(nrow(true.matrix)) else rownames(true.matrix)
  tar.names <- if (class(ranking$target) == "integer") seq_len(ncol(true.matrix)) else colnames(true.matrix)
  ranking2 <- ranking %>% mutate(
    tp=mapply(regulator, target, FUN=function(r, t) {
      if (r %in% reg.names & t %in% tar.names) true.matrix[[r, t]] else 0
    }),
    value=if("value" %in% colnames(ranking)) value else percent_rank(-seq_len(nrow(ranking)))
  )
  
  # iterate over ranking and calculate the eval df using given metrics
  rocr.pred <- prediction(ranking2$value, ranking2$tp)
  metrics <- do.call("cbind", lapply(perf.measures, function(p) performance(rocr.pred, p, x.measure="cutoff")@y.values[[1]]))
  colnames(metrics) <- perf.measures
  
  eval.df <- data.frame(
    cutoff=performance(rocr.pred, "acc", x.measure="cutoff")@x.values[[1]],
    balanced.acc=rowMeans(metrics[,c("rec", "spec")]),
    metrics
  )
  
  # calculate AUs
  auroc <- trapz(eval.df$fpr, eval.df$rec)
  aupr <- abs(trapz(eval.df$rec[-1], eval.df$prec[-1]))
  P <- sum(true.matrix, na.rm=T)
  if (P > 1) {
    aupr <- aupr / (1.0 - 1.0 / P)
  }
  F1 <- ifelse(auroc + aupr != 0, 2 * auroc * aupr / auroc + aupr, 0)
  au.df <- data.frame(auroc=auroc, aupr=aupr, F1=F1)
  
  # generate output
  list(metrics=eval.df, au.score=au.df)
}
