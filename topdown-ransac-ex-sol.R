##### ransaclm #####

ransaclm <- function(formula, data, error_threshold, inlier_threshold, iterations = 100, seed, early_exit = FALSE, early_exit_proportion = 0.7){
  
  
  # formula has to be formula
  checkmate::assert_class(x = formula, "formula")
  
  # all further checks build on formula
  attributes_formula <- attributes(terms.formula(formula, data = data))
  used_variables <- attributes_formula$term.labels
  response_variable <- dimnames(attributes_formula[["factors"]])[[1]][1]
  number_parameters <- length(used_variables)
  
  # only data used in model has to be checked.
  # Number of rows has to be greater then number of columns due to p < n
  # Number of cols have to be at least 2 for lm. otherwise intercept only model 
  checkmate::assert_data_frame(x = data[,c(response_variable,used_variables) , drop = FALSE],
                               types = "numeric", min.cols = 2, min.rows = (number_parameters + 1))
  sample_size <- nrow(data)
  # error_threshold has to be positiv 
  checkmate::assert_numeric(x = error_threshold, finite = TRUE, max.len = 1, lower = 0)
  
  # inlier_threshold must be count data and smaller then sample size
  checkmate::assert_integerish(x = inlier_threshold, lower = 1,
                               upper = sample_size, len = 1, any.missing = FALSE)
  
  # iterations also count data
  checkmate::assert_count(x = iterations, positive = TRUE)
  # seed also count
  checkmate::assert_count(x = seed, positive = TRUE)
  checkmate::check_logical(early_exit)
  checkmate::assert_numeric(x = early_exit_proportion, finite = TRUE, max.len = 1, lower = 0, upper = 1)
  
  
  set.seed(seed = seed)
  
  sample_size_inliers <- ncol(data)
  nrows <- seq_len(nrow(data))
  
  results <-  list()
  early_exit_index <- 0
  
  for (rep in seq_len(iterations)) {
    
    hypothetical_inliers_rows <- get_hypothetical_inliers_rows(sample_size_inliers, data)
    hypothetical_inliers <- data[hypothetical_inliers_rows, , drop =FALSE]
    maybe_model <- get_maybe_model(formula, hypothetical_inliers)
    consensus_set <- is_covered(maybe_model, data, hypothetical_inliers_rows, error_threshold)
    good_fit <- check_consensus_set(consensus_set, inlier_threshold)
    
    ifelse(good_fit == TRUE,
           yes = consensus_set_size <- sum(consensus_set == 1),
           no = consensus_set_size <- NA )
    
    results[[paste("rep", rep , sep = "_")]] <- list(matrix = hypothetical_inliers,
                                                     model = maybe_model, 
                                                     good_fit = good_fit, 
                                                     .consensus_set = consensus_set,
                                                     set_size = consensus_set_size)
    early_exit_index <- early_exit_index + 1
    
    if (early_exit == TRUE && good_fit == TRUE && (consensus_set_size / sample_size) >= early_exit_proportion) {
      break
    }
  }
  
  consensus_set_size_vector <- unlist(lapply(results, `[`, "set_size"))
  best <- which.max(consensus_set_size_vector)
  
  if (length(best) == 0) {
    return( warning( "No best model found. Adjust thresholds."))
  }
  
  best_fit <- results[[best]]
  consensus_set <- best_fit[[".consensus_set", drop = FALSE]]
  
  final_model <- lm(y~ . - inlier, data = subset(data, consensus_set))
  final_consensus_set <- is_covered(maybe_model = final_model, data = data,
                                    hypothetical_inliers_rows , error_threshold = error_threshold)
  
  final_results <- list( model = final_model,
                         data = cbind( data, .consensus_set = final_consensus_set ),
                         early_exit_index = early_exit_index)
  
  return(final_results)
}


# n Parameter (Spalten) => n Gleichungen (Zeilen)
get_hypothetical_inliers_rows <- function(sample_size_inliers, data){
  nrows <- seq_len(nrow(data))
  rows <- sample(nrows, replace = FALSE, size = sample_size_inliers)
}


get_maybe_model <- function(formula, hypothetical_inliers){
  lm(formula = formula, data = hypothetical_inliers)
}


is_covered <- function(maybe_model, data, hypothetical_inliers_rows, error_threshold){
  y_hat <- predict(maybe_model, newdata = data[-hypothetical_inliers_rows, , drop = FALSE])
  covered <- which(abs(data[-hypothetical_inliers_rows,"y"] - y_hat) < error_threshold)
  seq_len(nrow(data)) %in% covered
}


check_consensus_set <- function(consensus_set, inlier_threshold){
  sum(consensus_set == 1) >= inlier_threshold
}


# generate toy data with very clear inlier / outlier distinction to
# validate ransaclm. intercept is 0, coefs have value 1.
# inputs: as named -- number  of observations, number of coefs, fraction of inliers
# output: a data.frame with columns y, x.1, ..., x.<n_coef>, inlier
make_ransac_data <- function(n_obs, n_coef, inlier_fraction = 0.7) {
  coef <- rep(1, n_coef)
  design <- matrix(runif(n_obs * n_coef), n_obs, n_coef)
  inlier <- sample(
    c(TRUE, FALSE)[c(
      rep(1, n_obs * inlier_fraction),
      rep(2, n_obs - n_obs * inlier_fraction)
    )]
  )
  # E(y) = design %*% coef if inlier, else random draw from
  # uniform dist. on [-10 * n_coef, -2 * n_coef] and [2 * n_coef, 10 * n_coef].
  # error variance = n_coef * Var(U[0,1])
  inlier_expected_y <- design %*% coef
  outlier_expected_y <- sample(c(-1, 1), n_obs, replace = TRUE) *
    runif(n_obs, min = 2 * n_coef, max = 10 * n_coef)
  y <- ifelse(inlier, inlier_expected_y, outlier_expected_y) +
    rnorm(n_obs, sd = sqrt(1 / 12 * n_coef))
  data.frame(y, x = design, inlier)
}

#-------------------------------------------------------------------------------

# summarize & visualize ransaclm() results on data from "make_ransac_data". Only
# works if ransaclm's output list has a "data" entry with columns "inlier" as
# generated by make_ransac_data() as well as a column ".consensus_set" flagging
# the observations in the consensus set found by RANSAC.
validate_ransac <- function(ransacmodel, plot = TRUE) {
  checkmate::assert_class(ransacmodel[["model"]], "lm")
  checkmate::assert_data_frame(ransacmodel[["data"]])
  data <- ransacmodel[["data"]]
  checkmate::assert_logical(data[, "inlier"], any.missing = FALSE)
  checkmate::assert_logical(data[, ".consensus_set"], any.missing = FALSE)
  
  consensus_set <- data[, ".consensus_set"]
  true_inliers <- data[, "inlier"]
  cat("Inlier:\n")
  print(table(true = true_inliers, estimated = consensus_set))
  cat("\nCoefficients: (should be intercept ~0, effects ~1)\n")
  print(summary(ransacmodel[["model"]])$coefficients)
  if (plot &
      (!all(c("x", "y") %in% names(data)))) {
    warning("Can't plot this data, expecting columns 'x' and 'y'.")
    plot <- FALSE
  }
  if (plot) {
    plot_ransac(ransacmodel, data)
  }
  invisible(ransacmodel)
}

plot_ransac <- function(ransacmodel, data,
                        colors = c(
                          black = rgb(0, 0, 0, .5),
                          red   = rgb(1, 0, 0, .5),
                          blue  = rgb(0, 0, 1, .5)
                        )) {
  # scatterplot with true inliers in red, estimated inliers in blue
  # true regression line in black, RANSAC model in blue
  with(data,
       plot(x, y, col = colors[inlier + 1], pch = 19, bty = "n"))
  abline(c(0, 1), lwd = 2, col = colors[1])
  abline(ransacmodel$model, col = colors[3], lty = 2, lwd = 2)
  with(subset(data, .consensus_set),
       points(x, y, pch = 4, cex = 1.5, col = colors[3]))
  legend("top",
         lty = c(NA, NA, 1, 2, NA), lwd = c(NA, NA, 2, 2, NA),
         pch = c(19, 19, NA, NA, 4), col = colors[c(1, 2, 1, 3, 3)],
         legend = c(
           "'true' outliers", "'true' inliers", "true regression line",
           "RANSAC estimate", "RANSAC consensus set"
         ),
         cex = .7, bg = NA, inset = c(0, -.15), ncol = 3, xpd = NA
  )
}


# immer set.seed() um Ergebnisse reproduzierbar zu machen...
set.seed(1874374111)
data_simple <- make_ransac_data(n_obs = 100, n_coef = 1, inlier_fraction = 0.7)

# univariate example:
ransac_simple <- ransaclm(y ~ . - inlier,
                          data = data_simple, error_threshold = 2,
                          inlier_threshold = 50, seed = 20171111)

validate_ransac(ransac_simple)



# Tests
library(testthat)

test_that("basic implementation works", {
  data_simple <- make_ransac_data(n_obs = 100, n_coef = 1, inlier_fraction = 0.7)
  ransaclm_simple <-  ransaclm(y ~ . - inlier,
                              data = data_simple, error_threshold = 2,
                              inlier_threshold = 50, seed = 20171111)
  expect_equal(class(ransac_simple[["model"]]), "lm")
  expect_equal(class(ransac_simple[["data"]]), "data.frame")
  expect_vector(ransac_simple[["data"]][,".consensus_set"], size = nrow(data_simple), ptype = logical())
  })

test_that("early exit works", {
  data_simple <- make_ransac_data(n_obs = 100, n_coef = 1, inlier_fraction = 0.7)
  ransaclm_simple <-  ransaclm(y ~ . - inlier,
                               data = data_simple, error_threshold = 2,
                               inlier_threshold = 45, early_exit = TRUE, early_exit_proportion = 0.68,
                               seed = 20171111)
  ransaclm_early_false <-  ransaclm(y ~ . - inlier,
                               data = data_simple, error_threshold = 2,
                               inlier_threshold = 45, early_exit = FALSE,
                               seed = 20171111)
  expect_lt(ransaclm_simple[["early_exit_index"]], expected = ransaclm_early_false[["early_exit_index"]])
  
})

test_that("ransaclm works for more then one coef", {
  data_simple <- make_ransac_data(n_obs = 100, n_coef = 2, inlier_fraction = 0.7)
  two_coef_ransaclm <- ransaclm(y ~ . - inlier,
           data = data_simple, error_threshold = 2,
           inlier_threshold = 50, seed = 20171111)
  expect_equal(class(two_coef_ransaclm[["model", drop = FALSE]]), "lm")
  expect_equal(class(ransac_simple[["data"]]), "data.frame")
})

test_that("ransaclm gives warnings if no model found", {
  data_simple <- make_ransac_data(n_obs = 100, n_coef = 2, inlier_fraction = 0.7)
  expect_warning(ransaclm(y ~ . - inlier,
                          data = data_simple, error_threshold = 2,
                          inlier_threshold = 90, seed = 20171111))
})

test_that("ransaclm gives error if used data is not appropriate", {
  data_simple_with_factor <- cbind(make_ransac_data(n_obs = 100, n_coef = 1, inlier_fraction = 0.7), factor_var = as.factor(c(rep(0,50), rep(1,50))))
  expect_error(ransaclm(y ~ . - inlier,
                          data = data_simple_with_factor, error_threshold = 2,
                          inlier_threshold = 90, seed = 20171111))
})

test_that("ransaclm gives warnings if p > n", {
  data_simple <- data.frame(y = rnorm(2), x1 = rnorm(2), x2 = rnorm(2))
  expect_error(ransaclm(y ~ . ,
                          data = data_simple, error_threshold = 4,
                          inlier_threshold = 2, seed = 20171111))
})


# Parallelisation schwer zu implementieren wenn early exit. Also ohne
library(foreach)
ransaclm_parallel <- function(formula, data, error_threshold, inlier_threshold, iterations = 100, seed){
  
  
  # formula has to be formula
  checkmate::assert_class(x = formula, "formula")
  
  # all further checks build on formula
  attributes_formula <- attributes(terms.formula(formula, data = data))
  used_variables <- attributes_formula$term.labels
  response_variable <- dimnames(attributes_formula[["factors"]])[[1]][1]
  number_parameters <- length(used_variables)
  
  # only data used in model has to be checked.
  # Number of rows has to be greater then number of columns due to p < n
  # Number of cols have to be at least 2 for lm. otherwise intercept only model 
  checkmate::assert_data_frame(x = data[,c(response_variable,used_variables) , drop = FALSE],
                               types = "numeric", min.cols = 2, min.rows = (number_parameters + 1))
  sample_size <- nrow(data)
  # error_threshold has to be positiv 
  checkmate::assert_numeric(x = error_threshold, finite = TRUE, max.len = 1, lower = 0)
  
  # inlier_threshold must be count data and smaller then sample size
  checkmate::assert_integerish(x = inlier_threshold, lower = 1,
                               upper = sample_size, len = 1, any.missing = FALSE)
  
  # iterations also count data
  checkmate::assert_count(x = iterations, positive = TRUE)
  # seed also count
  checkmate::assert_count(x = seed, positive = TRUE)
  
  set.seed(seed = seed)
  
  sample_size_inliers <- ncol(data)
  nrows <- seq_len(nrow(data))
  
  results <-  list()
  doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  
  foreach(rep = seq_len(iterations)) %do% {
    
    hypothetical_inliers <- get_hypothetical_inliers(sample_size_inliers, data)
    maybe_model <- get_maybe_model(formula, hypothetical_inliers)
    consensus_set <- is_covered(maybe_model, data, error_threshold)
    good_fit <- check_consensus_set(consensus_set, inlier_threshold)
    
    ifelse(good_fit == TRUE,
           yes = consensus_set_size <- sum(consensus_set == 1),
           no = consensus_set_size <- NA )
    
    results[[paste("rep", rep , sep = "_")]] <- list(matrix = hypothetical_inliers,
                                                     model = maybe_model, 
                                                     good_fit = good_fit, 
                                                     .consensus_set = consensus_set,
                                                     set_size = consensus_set_size)
    
    }
  
  consensus_set_size_vector <- unlist(lapply(results, `[`, "set_size"))
  best <- which.max(consensus_set_size_vector)
  
  if (length(best) == 0) {
    return( warning( "No best model found. Adjust thresholds."))
  }
  
  best_fit <- results[[best]]
  consensus_set <- best_fit[[".consensus_set", drop = FALSE]]
  
  final_model <- lm(y~ . - inlier, data = subset(data, consensus_set))
  final_consensus_set <- is_covered(final_model, data, error_threshold)
  
  final_results <- list( model = final_model,
                         data = cbind( data, .consensus_set = final_consensus_set ))
  
  return(final_results)
}

ransac_para <- ransaclm_parallel(y ~ . - inlier,
                          data = data_simple, error_threshold = 2,
                          inlier_threshold = 50, seed = 20171111)

validate_ransac(ransac_para) # same results as above
data_complex <- make_ransac_data(n_obs = 100000, n_coef = 10, inlier_fraction = 0.7)

rbenchmark::benchmark(ransaclm_parallel(y ~ . - inlier,
                                           data = data_complex, error_threshold = 2,
                                           inlier_threshold = 50, seed = 20171111),
                      ransaclm(y ~ . - inlier,
                                        data = data_complex, error_threshold = 2,
                                        inlier_threshold = 50, seed = 20171111),
                      replications = 100, columns = c("test", "elapsed", "relative"),
                      order = "elapsed")

