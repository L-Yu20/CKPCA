# Parameters
n <- 800

# Y-axis data (updated)
y <- c(0.925, 0.975, 0.989, 0.989, 0.991, 0.985, 0.982, 0.316, 0.536)

# X-axis data (actual values)
x <- c(
  2,
  floor(n^(1/5)), floor(n^(1/3)), floor(n^(3/7)), floor(n^(1/2)),
  floor(n^(4/7)), floor(n^(3/5)), floor(n^(3/4)), floor(n / 3)
)

# Plot
plot(x, y,
     type = "o", lwd = 2,
     xaxt = "n", ylim = c(0.3, 1),  # adjust y-limit for new data
     xlab = expression(beta[n]), 
     ylab = "Rand Index",
     main = "Rand Index in Case 2"
)

# X-axis
axis(1, at = x, labels = x)

# Emphasize the point at beta_n = 28
mtext("28", side = 1, at = 28, line = 1, col = "black", font = 2)
abline(v = 28, col = "red", lty = 2, lwd = 2)

# Add grid
grid()
