# Parameters
n <- 800

# Y-axis data
y <- c(0.996, 0.996, 0.997, 0.996, 0.997, 0.996, 0.996, 0.986, 0.485)

# X-axis labels
x_labels <- c(
  2,
  floor(n^(1/5)), floor(n^(1/3)), floor(n^(3/7)), floor(n^(1/2)),
  floor(n^(4/7)), floor(n^(3/5)), floor(n^(3/4)), floor(n / 3)
)

# Equally spaced x positions
x_index <- 1:length(y)

# Focus on the region around the optimum
ylim_focus <- c(0.9957, 0.9972)

# Do not connect points outside the displayed y-range
y_plot <- y
y_plot[y < ylim_focus[1] | y > ylim_focus[2]] <- NA

par(mar = c(5, 5.2, 4, 2) + 0.1)

# Plot
plot(x_index, y_plot,
     type = "o",
     lwd = 2,
     pch = 1,
     cex = 1.1,
     xaxt = "n",
     yaxt = "n",
     ylim = ylim_focus,
     xlab = expression(beta[n]),
     ylab = "",
     main = "Rand Index in Case 1"
)

axis(1, at = x_index, labels = x_labels)
axis(2, at = c(0.9960, 0.9965, 0.9970), las = 1, cex.axis = 0.9)

mtext("Rand Index", side = 2, line = 3.8)

# Highlight beta_n = 28
abline(v = 5, col = "red", lty = 2, lwd = 2)
points(5, y[5], pch = 16, cex = 1.2)
mtext("28", side = 1, at = 5, line = 1, font = 2)

# Mark truncated points at the lower boundary
clipped_idx <- which(y < ylim_focus[1])
points(clipped_idx, rep(ylim_focus[1], length(clipped_idx)),
       pch = 25, bg = "black", cex = 1)

# Grid
grid(nx = NA, ny = NULL, lty = 3, col = "grey80")
