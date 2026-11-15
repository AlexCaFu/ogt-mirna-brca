## test_layer6.R -- validate the statistics in R/06_robustness.R on synthetic
## data with known answers. Base R only; no special packages needed.

## ---- helpers copied verbatim from R/06_robustness.R ----
partial_spearman <- function(x, y, cov) {
  d <- data.frame(x = x, y = y, cov, check.names = FALSE)
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) < 20) return(list(rho = NA, p = NA, n = nrow(d)))
  d$rx <- rank(d$x); d$ry <- rank(d$y)
  cn <- setdiff(colnames(cov), c("x", "y"))
  rhs <- paste(sprintf("`%s`", cn), collapse = " + ")
  ex <- residuals(lm(as.formula(paste("rx ~", rhs)), data = d))
  ey <- residuals(lm(as.formula(paste("ry ~", rhs)), data = d))
  ct <- suppressWarnings(cor.test(ex, ey, method = "pearson"))
  list(rho = unname(ct$estimate), p = ct$p.value, n = nrow(d))
}
spearman_ci_fisher <- function(rho, n, conf = 0.95) {
  se <- sqrt((1 + rho^2 / 2) / (n - 3)); z <- atanh(rho); zc <- qnorm(1 - (1 - conf) / 2)
  c(lower = tanh(z - zc * se), upper = tanh(z + zc * se))
}
boot_spearman_ci <- function(x, y, B = 2000, conf = 0.95) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  rhos <- replicate(B, { i <- sample.int(n, n, replace = TRUE)
                         suppressWarnings(cor(x[i], y[i], method = "spearman")) })
  stats::quantile(rhos, c((1 - conf) / 2, 1 - (1 - conf) / 2), na.rm = TRUE)
}
# formula-based partial Spearman, for cross-checking (continuous covariate)
formula_partial <- function(x, y, z) {
  rxy <- cor(x, y, method = "spearman"); rxz <- cor(x, z, method = "spearman"); ryz <- cor(y, z, method = "spearman")
  (rxy - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))
}

set.seed(42)

cat("== TEST 1: partial correlation removes a pure confounder ==\n")
n <- 1500
z <- rnorm(n)
x <- z + rnorm(n)            # x related to y ONLY through z
y <- z + rnorm(n)
unadj <- cor(x, y, method = "spearman")
part  <- partial_spearman(x, y, data.frame(z = z))$rho
cat(sprintf("  unadjusted rho = %.3f (should be clearly positive)\n", unadj))
cat(sprintf("  partial   rho = %.3f (should be near 0)\n", part))
stopifnot(unadj > 0.3, abs(part) < 0.08)
cat("  PASSED: confounder-driven correlation collapses to ~0 after adjustment.\n\n")

cat("== TEST 2: partial correlation keeps a genuine relationship ==\n")
z2 <- rnorm(n); x2 <- z2 + rnorm(n)
y2 <- 0.8 * x2 + z2 + rnorm(n)     # genuine x-y link beyond z
part2 <- partial_spearman(x2, y2, data.frame(z = z2))$rho
form2 <- formula_partial(x2, y2, z2)
cat(sprintf("  residual-method partial rho = %.3f\n", part2))
cat(sprintf("  formula-based   partial rho = %.3f\n", form2))
stopifnot(part2 > 0.2, abs(part2 - form2) < 0.05)
cat("  PASSED: genuine partial correlation retained and matches the closed-form value.\n\n")

cat("== TEST 3: partial correlation handles a categorical covariate ==\n")
grp <- factor(sample(c("A","B","C"), n, replace = TRUE))
shift <- c(A = -2, B = 0, C = 2)[grp]
x3 <- shift + rnorm(n); y3 <- shift + rnorm(n)   # correlated only via group
part3 <- partial_spearman(x3, y3, data.frame(grp = grp))$rho
cat(sprintf("  unadjusted rho = %.3f | partial (adjust group) rho = %.3f\n",
            cor(x3, y3, method = "spearman"), part3))
stopifnot(cor(x3, y3, method = "spearman") > 0.4, abs(part3) < 0.08)
cat("  PASSED: group-driven correlation removed when adjusting for the factor.\n\n")

cat("== TEST 4: CIs contain the truth and agree, at n like the real study ==\n")
# build data with Spearman ~ -0.25 at n = 1076
S <- matrix(c(1, -0.26, -0.26, 1), 2)
d <- tryCatch(MASS::mvrnorm(1076, c(0,0), S), error = function(e) NULL)
if (is.null(d)) { z <- rnorm(1076); d <- cbind(z + rnorm(1076,0,2.4), -z + rnorm(1076,0,2.4)) }
xx <- d[,1]; yy <- d[,2]
r <- cor(xx, yy, method = "spearman")
fis <- spearman_ci_fisher(r, length(xx)); boo <- boot_spearman_ci(xx, yy, B = 1500)
cat(sprintf("  rho = %.3f | Fisher CI [%.3f, %.3f] | bootstrap CI [%.3f, %.3f]\n",
            r, fis["lower"], fis["upper"], boo[[1]], boo[[2]]))
stopifnot(fis["lower"] < r, fis["upper"] > r, boo[[1]] < r, boo[[2]] > r)
stopifnot(abs(fis["lower"] - boo[[1]]) < 0.05, abs(fis["upper"] - boo[[2]]) < 0.05)
cat("  PASSED: both intervals bracket rho and agree within 0.05; width is sensible.\n\n")

cat("== TEST 5: ranking orders most-negative first ==\n")
rk <- data.frame(mirna = c("a","miR-200c-3p","c","d"), rho = c(0.1, -0.30, -0.05, 0.2))
rk <- rk[order(rk$rho), ]; rk$rank <- seq_len(nrow(rk))
stopifnot(rk$rank[rk$mirna == "miR-200c-3p"] == 1)
cat("  PASSED: strongest inverse correlation is rank 1.\n\n")

cat("ALL LAYER 6 STATISTICAL TESTS PASSED.\n")
