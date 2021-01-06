##### ransaclm #####

ransaclm <- function(formula, data, error_threshold, inlier_threshold, iterations = 100, seed){
  
  # Anfoderungen an Inputs
  # formula muss formula sein
  # data muss numerisch sein, keine Faktoren, wegen samplen
  # error_threshold muss positiv sein
  # inlier_threshold muss count data sein
  # Interartion auch count data
  # seed auch count
  
  sample_size <- ncol(data)
  nrows <- seq_len(nrow(data))
  
  results <-  list()
  
  for (rep in seq_len(iterations)) {
    
    hypothetical_inliers <- get_hypothetical_inliers(sample_size, data) # sample from each column one value
    maybe_model <- get_maybe_model(formula, hypothetical_inliers)
    consensus_set <- is_covered(maybe_model, data, error_threshold)
    good_fit <- check_consensus_set(consensus_set, inlier_threshold)
    
    if (good_fit == TRUE) {
      consensus_set_size <- nrow(consensus_set)
    } else {
      consensus_set_size <- NA
    }
    results[[paste("rep", rep , sep = "_")]] <- list(matrix = hypothetical_inliers,
                                                     model = maybe_model, 
                                                     good_fit = good_fit, 
                                                     set_size = consensus_set_size)
  }
  
  consensus_set_size_vector <- unlist(lapply(results, `[`, "set_size"))
  best_fit <- which.max(consensus_set_size_vector)  
  final_model <- lm(y~ . - inlier, data = subset(data, consensus_set))
  consensus_set <- best_fit[[".consensus_set", drop = FALSE]]
  final_results <- list( model = final_model,
                         data = cbind( data, .consensus_set = consensus_set ))
  
  return(final_results)
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
  y_hat <- predict(maybe_model, newdata = data)
  covered <- which(abs(data[,"y"] - y_hat) < error_threshold)
  seq_len(nrow(data)) %in% covered
}



consensus_set <- is_covered(maybe_model, test, error_threshold = 1)


check_consensus_set <- function(consensus_set, inlier_threshold){
  nrow(consensus_set) >= inlier_threshold
}

check_consensus_set(consensus_set = consensus_set, inlier_threshold = 5)

get_best_fit <- function(consensus_set_size_vector){
  which.max(consensus_set_size_vector)
}

c(test[[2]], best_fit = TRUE)
