##### ransaclm #####

ransaclm <- function(formula, data, error_threshold, inlier_threshold, iterations, seed){
  
  # Anfoderungen an Inputs
  # formula muss formula sein
  # data muss numerisch sein, keine Faktoren, wegen samplen
  # error_threshold muss positiv sein
  # inlier_threshold muss count data sein
  # Interartion auch count data
  # seed auch count
  sample_size <- ncol(data)
  nrows <- seq_len(nrow(data))
  
  results <-  rep(list(list(model = NA, data = NA)), iterations)
  for(rep in seq_len(iterations)){
  hypothetical_inliers <- get_hypothetical_inliers(sample_size, data) # sample from each column one value
  maybe_model <- get_maybe_model(formula, hypothetical_inliers) 
  results[[rep]]["model"] <- maybe_model
  consensus_set <- is_covered(maybe_model, data, error_threshold)
  good_fit <- check_consensus_set(consensus_set, inlier_threshold)
  if(good_fit = TRUE) {
  consensus_set_size <- nrow(consensus_set)
  } else next
  
  
  best_fit <- get_best_fit(consensus_set_size)
  
  return(results)
}

test <- data.frame(x = rnorm(10) + 5, y = rnorm(10,5))

# n Parameter (Spalten) => n Gleichungen (Zeilen)
get_hypothetical_inliers <- function(sample_size, data){
  nrows <- seq_len(nrow(data))
  rows <- sample(nrows, replace = FALSE, size = sample_size)
  data[rows, ,drop = FALSE]
}

test_hypo <- get_hypothetical_inliers(sample_size = 2, data = test)

get_maybe_model <- function(formula, hypothetical_inliers){
  lm(formula = formula, data = hypothetical_inliers)
}

maybe_model <- get_maybe_model(formula = y~x, hypothetical_inliers = test_hypo)

is_covered <- function(maybe_model, data, error_threshold){
  upr <- (coef(maybe_model)[1] + error_threshold) + coef(maybe_model)[2] * data[, "x" , drop = FALSE]
  lwr <- (coef(maybe_model)[1] - error_threshold) + coef(maybe_model)[2] * data[, "x" , drop = FALSE]
  covered <- which(lwr <= data[, "y", drop = FALSE] & data[, "y", drop = FALSE] <= upr)
  data[covered, , drop = FALSE]
}



consensus_set <- is_covered(maybe_model, test, error_threshold = 1)


check_consensus_set <- function(consensus_set, inlier_threshold){
  nrow(consensus_set) >= inlier_threshold
}

check_consensus_set(consensus_set = consensus_set, inlier_threshold = 5)

c(test[[2]], best_fit = TRUE)
