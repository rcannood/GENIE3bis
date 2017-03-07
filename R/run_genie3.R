#' Compute weighted adjacency matrix of inferred network
#'
#' @param data A data_frame of observations of the different genes. The rows must contain the observations.
#' @param num_candidate_genes The choice of number of input genes randomly selected as candidates at each node.
#'   Must be \code{"all"} for all input features, \code{"sqrt"} for the square root of all input features (default), or an integer.
#' @param num_trees The number of trees in ensemble for each target gene (default 1000).
#' @param regulators A set of indices or column names of entities whose observed values regulate the observed values of the targets.
#' @param targets A set of indices or column names of entities whose observed values are regulated by the regulators.
#' @param max_interactions The maximum number of interactions to be returned by GENIE3.
#' @param importance_measure Variable importance mode, one of impurity' or 'permutation'.
#'   The 'impurity' measure calculates the variance of the responses in each tree node and
#'   the 'permutation' calculates the increase of MSE after permutations of the regulators.
#' @param seed A random number generator seed for replication of analyses. NULL means the seed is not set.
#' @param verbose Output additional information.
#' @param scale_genes Whether the genes should be scaled. This is recommended in order to make importance values comparable to one another.
#' @param parallel_type Either the number of threads to use for parallel execution, or a \code{\link[PRISM]{qsub_configuration}} object.
#'
#' @references Huynh-Thu, V. A. et al. (2010) Inferring Regulatory Networks from Expression Data Using Tree-Based Methods. PLoS ONE.
#'
#' @return The weighted adjacency matrix of inferred network.
#' @export
#'
#' @importFrom ranger ranger
#' @import dplyr
#'
#'
#' @examples
#'
#' library(GENIE3)
#' library(dplyr)
#'
#' # generate random data
#' data <- dplyr::as_data_frame(matrix(runif(100 * 100), ncol = 100))
#' genes <- colnames(data)
#' regulators <- genes[1:20]
#' targets <- genes
#'
#' true_interactions <-
#'   expand.grid(
#'     regulator = factor(regulators, levels = genes),
#'     target = factor(targets, levels = genes)) %>%
#'   sample_n(100)
#'
#' # run GENIE3
#' ranking <- run_genie3(data, regulators, targets)
#'
#' # evaluate performance
#' eval <- evaluate_ranking(ranking, true_interactions, regulators, targets)
#' eval$area_under
#' plot_curves(eval)
#'
#' # evaluate multiple rankings at the same time
#' ranking_cor <- cor(data[,regulators], data[,targets]) %>%
#'   reshape2::melt(varnames = c("regulator", "target"), value.name = "importance") %>%
#'   arrange(desc(importance)) %>%
#'   mutate(regulator = factor(as.character(regulator), levels = genes),
#'          target = factor(as.character(target), levels = genes)) %>%
#'   filter(regulator != target)
#' rankings <- list(GENIE3=ranking, Correlation=ranking_cor)
#' evals <- evaluate_multiple_rankings(rankings, true_interactions, regulators, targets)
#' evals$area_under
#' plot_curves(evals)
#'
#' # run GENIE3 in parallel
#' ranking <- run_genie3(data, regulators, targets, parallel_type = 8)
#'
#' # run GENIE3 with PRISM
#' ranking <- run_genie3(data, regulators, targets, parallel_type = PRISM::override_qsub_config())
#'
#' # run GENIE3 with PRISM without waiting
#' handle <- run_genie3(data, regulators, targets, parallel_type = PRISM::override_qsub_config(wait = F))
#' ranking <- retrieve_genie3_output(handle)
run_genie3 <- function(data, regulators, targets,
                       num_candidate_genes = "sqrt", num_trees = 1000,
                       max_interactions = 100000, importance_measure = "impurity",
                       seed = NULL, verbose = T, scale_genes = T,
                       parallel_type = 1) {
  requireNamespace("randomForest")

  # Check data ----------------------------------------------------------------
  if (is.matrix(data)) {
    data <- as_data_frame(data)
  }
  if (!is.data.frame(data) || any(!apply(data, 2, is.finite))) {
    stop("Parameter ", sQuote("data"), " must be a data frame consisting of finite numerical values.")
  }
  if(length(unique(colnames(data))) != ncol(data)) {
    stop("Colnames of ", sQuote("data"), " must be unique.")
  }

  num_samples <- nrow(data)
  num_genes <- ncol(data)
  gene_names <- colnames(data)

  # Check num_trees -----------------------------------------------------------
  if (!is.numeric(num_trees) || num_trees < 1) {
    stop("Parameter ", sQuote("num_trees"), " must be a numeric larger than or equal to 1")
  }

  # Check max_interactions ----------------------------------------------------
  if (!is.numeric(max_interactions) || max_interactions < 1) {
    stop("Parameter ", sQuote("max_interactions"), " must be a numeric larger than or equal to 1")
  }

  # Check regulators ----------------------------------------------------------
  if (is.null(regulators)) {
    regulators <- seq_len(num_genes)
  }
  if (is.character(regulators)) {
    regulators <- match(regulators, gene_names)
  }
  if (!is.numeric(regulators) || !all(regulators %in% seq_len(num_genes)) ||
      length(unique(regulators)) != length(regulators)) {
    stop("Parameter ", sQuote("regulators"), " must either be NULL, or ",
         "be a subset of either colnames(data) or seq_len(ncol(data)), ",
         "and may not contain any repeated elements.")
  }
  num_regulators <- length(regulators)
  regulator_names <- gene_names[regulators]

  # Check targets -------------------------------------------------------------
  if (is.null(targets)) {
    targets <- seq_len(num_genes)
  }
  if (is.character(targets)) {
    targets <- match(targets, gene_names)
  }
  if (!is.numeric(targets) || !all(targets %in% seq_len(num_genes)) ||
      length(unique(targets)) != length(targets)) {
    stop("Parameter ", sQuote("targets"), " must either be NULL, or ",
         "be a subset of either colnames(data) or seq_len(ncol(data)), ",
         "and may not contain any repeated elements.")
  }
  num_targets <- length(targets)
  target_names <- gene_names[targets]

  # Check num_candidate_genes -------------------------------------------------
  if (is.numeric(num_candidate_genes)) {
    mtry <- num_candidate_genes
  } else if (num_candidate_genes == "sqrt") {
    mtry <- round(sqrt(num_regulators))
  } else if (num_candidate_genes == "all") {
    mtry <- num_regulators-1
  } else {
    stop("Parameter ", sQuote("num_candidate_genes"), " must be \"sqrt\", or \"all\", or an integer.")
  }

  # Check importance_measure --------------------------------------------------
  if (!importance_measure %in% c("impurity", "permutation")) {
    stop("Parameter ", sQuote("importance_measure"), " must be \"impurity\" or \"permutation\".")
  }

  # Check seed ----------------------------------------------------------------
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Check parallel_type -------------------------------------------------------
  if (is.null(parallel_type)) {
    parallel_type <- 1
  }
  if (is.numeric(parallel_type)) {
    if (parallel_type < 1) {
      stop("Parameter ", sQuote("mc.cores"), " must be larger than or equal to 0.")
    } else if (parallel_type == 1) {
      lapply_function <- lapply
    } else {
      requireNamespace("parallel")
      lapply_function <- function(...) parallel::mclapply(..., mc.cores = parallel_type)
    }
  } else if (PRISM:::is_qsub_config(parallel_type)) {
    requireNamespace("PRISM")
    qsub_environment <- c("data", "gene_names", "regulators", "num_trees",
                          "importance_measure", "regulator_names", "num_regulators",
                          "verbose")

    lapply_function <- function(X, FUN) {
      qsub_config <- PRISM::override_qsub_config(parallel_type, name = "GENIE3")
      PRISM::qsub_lapply(
        X = X,
        FUN = FUN,
        qsub_environment = qsub_environment,
        qsub_config = qsub_config
      )
    }
  }

  if (verbose) {
    cat("GENIE3 parameter checks are OK!\n")
    cat("Starting GENIE3 computations.\n")
    flush.console()
  }

  # Scale the data
  if (scale_genes) {
    data <- scale(data)
  }

  # Compute importances for every target
  lapply_output <- lapply_function(targets, function(target_index) {
    target_name <- gene_names[[target_index]]
    regs <- setdiff(regulators, target_index)

    if (verbose) {
      cat("Computing for target ", ifelse(is.null(gene_names), target.index, target_name), "\n", sep="")
      flush.console()
    }

    if (requireNamespace("ranger", quietly = T)) {
      if (verbose) {
        cat("Using package 'ranger'\n")
        flush.console()
      }
      data_ix <- data[,c(regs, target_index)]
      target_formula <- as.formula(paste0("`", target_name, "` ~ ."))

      rf <- ranger::ranger(
        formula = target_formula,
        data = as.data.frame(data_ix),
        importance = importance_measure,
        num.trees = num_trees,
        mtry = mtry
      )

      importance_values <- rf$variable.importance
    } else if (requireNamespace("randomForest", quietly = T)) {
      if (verbose) {
        cat("Using package 'randomForest'\n")
        flush.console()
      }

      x <- data[, regs, drop = F]
      y <- data[, target_index]

      importance <- c("impurity"="IncNodePurity", "permutation"="%IncMSE")[importance_measure]
      calculate_IncMSE <- importance == "%IncMSE"
      rf <- randomForest::randomForest(
        x = x,
        y = y,
        ntree = num_trees,
        mtry = mtry,
        importance = calculate_IncMSE
      )

      importance_values <- rf$importance[,importance]
    } else {
      stop("Neither packages ", sQuote("ranger"), " nor ", sQuote("randomForest"), " are installed.")
    }

    tibble::data_frame(
      regulator = factor(gene_names[regs], levels = gene_names),
      target = factor(target_name, levels = gene_names),
      importance = importance_values
    )
  })

  handle <- list(
    lapply_output = lapply_output,
    parallel_type = parallel_type,
    gene_names = gene_names,
    regulator_names = regulator_names,
    target_names = target_names,
    max_interactions = max_interactions
  )

  if (requireNamespace("PRISM", quietly = T) &&
      PRISM:::is_qsub_config(parallel_type) &&
      !parallel_type$wait)  {
    handle
  } else {
    retrieve_genie3_output(handle)
  }
}

#' Retrieve GENIE3 output
#'
#' Necessary when retrieving the output in combination with the PRISM package.
#'
#' @param handle
#'
#' @import dplyr
#'
#' @export
#'
#' @seealso \code{\link{run_genie3}}
retrieve_genie3_output <- function(handle) {
  if (requireNamespace("PRISM", quietly = T) &&
      PRISM:::is_qsub_config(handle$parallel_type) &&
      !handle$parallel_type$wait) {
    outputs <- PRISM::qsub_retrieve(handle$lapply_output)
  } else {
    outputs <- handle$lapply_output
  }

  ranking <- bind_rows(outputs) %>% arrange(desc(importance))

  if (nrow(ranking) > handle$max_interactions) {
    ranking <- ranking %>% top_n(handle$max_interactions, importance)
  }

  ranking
}
