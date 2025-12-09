randi <- function(posi, num, n, times, true_cps) {
  # Construct true labels
  true_labels <- rep(0, n)
  segs <- c(0, true_cps, n)
  for (i in 1:(length(segs) - 1)) {
    true_labels[(segs[i] + 1):segs[i + 1]] <- i
  }
  
  rand_indices <- numeric(times)
  total_pairs <- n * (n - 1) / 2
  
  for (j in 1:times) {
    pred_labels <- rep(0, n)
    
    if (num[j] == 0) {
      pred_labels[1:n] <- 1
    } else if (num[j] == 1) {
      pred_labels[1:posi[j, 1]] <- 1
      if (posi[j, 1] < n) {
        pred_labels[(posi[j, 1] + 1):n] <- 2
      }
    } else {
      pred_labels[1:posi[j, 1]] <- 1
      if (posi[j, num[j]] < n) {
        pred_labels[(posi[j, num[j]] + 1):n] <- num[j] + 1
      }
      for (i in 1:(num[j] - 1)) {
        pred_labels[(posi[j, i] + 1):posi[j, i + 1]] <- i + 1
      }
    }
    
    tp <- 0
    tn <- 0
    for (a in 2:n) {
      for (b in 1:(a - 1)) {
        same_pred <- pred_labels[a] == pred_labels[b]
        same_true <- true_labels[a] == true_labels[b]
        if (same_pred && same_true) {
          tp <- tp + 1
        } else if (!same_pred && !same_true) {
          tn <- tn + 1
        }
      }
    }
    
    rand_indices[j] <- (tp + tn) / total_pairs
  }
  
  # Calculate mean and 95% confidence interval
  mean_ri <- mean(rand_indices)
  ci <- quantile(rand_indices, probs = c(0.025, 0.975))
  
  # Return a list with all three
  return(list(mean = mean_ri, CI_lower = ci[1], CI_upper = ci[2]))
}