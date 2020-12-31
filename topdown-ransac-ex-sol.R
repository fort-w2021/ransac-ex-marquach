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
  rows <- seq_len(nrow(data))
  
  results <-  rep(list(list(model = NULL, data = NULL)), iterations)
  
  hypothetical_inliers <- get_hypothetical_inliers(sample_size, rows, data) # sample from each column one value
  maybe_model <- get_maybe_model(formula, hypothetical_inliers) 
  consensus_set <- is_covered(maybe_model, data, error_threshold)
  good_fit <- check_consensus_set(consensus_set, inlier_threshold)
  
  best_fit <- get_best_fit(good_fit)
  
  return(results)
}

# Zwei Parameter Zwei Gleichungen
get_hypothetical_inliers <- function(sample_size, data){
  
  hypothetical_inliers <- matrix(data = NA, nrow = sample_size, ncol = sample_size)
  colnames(hypothetical_inliers) <- colnames(data)
  rows <- seq_len(nrow(data))
  
  for(s in seq_len(sample_size)){
    hypothetical_inliers[,s  ] <- data[sample(rows, size = sample_size), s]   
    # Warum funktioniert hier drop = FALSE nicht? (Fehler: falsche Anzahl von Indizes für Matrix). In test[1:2, 1, drop = FALSE] funktioniert es ganz normal
  }
  
  hypothetical_inliers
}

get_hypothetical_inliers(sample_size = 2, data = test)
test
get_maybe_model(formula, hypothetical_inliers) 
get_maybe_model <- function(formula, hypothetical_inliers){
  coefs <- coef(lm(formula = formula, data = hypothetical_inliers))
}
get_maybe_model(formula = y~x, hypothetical_inliers = as.data.frame(test_hypo))

hypothetical_inliers <- matrix(data = NA, nrow = 2, ncol = 2)
for(s in seq_len(ncol(test))){
  hypothetical_inliers[,s  ] <- test[sample(1:5, size = 2), s, drop = FALSE]
}
test[sample(1:5,2), 2, drop = FALSE]
