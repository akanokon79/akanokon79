# Monte Carlo Simulation Study: Evaluating QR, LQR, FLQR, and BFLQR
# Model Performance, Variable Selection, and Robustness to Outliers
# Required packages

library(quantreg)
library(GIGrvg)
library(MASS)

# Set seed for exact numerical reproducibility
set.seed(2026)

# 1. Simulation Parameters & Settings

N    <- c(50, 100, 200)       # Sample sizes
TAU  <- c(0.25, 0.50, 0.75)   # Quantiles
R    <- 500                   # Replications per scenario
ITER <- 6000                  # MCMC total iterations
BURN <- 2000                  # MCMC burn-in iterations
p    <- 10                    # Number of predictors

# True Parameter Configuration
b0 <- 1.5
bt <- c(2, 2, 1.5, 0, 0, -1.5, -1.5, 0, 1, 1)

true_active <- (bt != 0)
true_zero   <- (bt == 0)

# AR(1) Predictor Correlation Structure
rho <- 0.5
S   <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))


# 2. Data Generation & Contamination Mechanisms

# Asymmetric Laplace Distribution (ALD) Random Generation
rALD <- function(n, tau, sigma = 1) {
  phi1 <- (1 - 2 * tau) / (tau * (1 - tau))
  phi2 <- 2 / (tau * (1 - tau))
  v <- rexp(n, 1 / sigma)
  z <- rnorm(n)
  phi1 * v + sqrt(phi2 * sigma * v) * z
}

# Generate Clean Synthetic Dataset
gendata <- function(n, tau) {
  X <- matrix(rnorm(n * p), n, p) %*% chol(S)
  y <- as.numeric(b0 + X %*% bt + rALD(n, tau))
  list(X = X, y = y)
}

# Apply 10% High-Leverage/Large-Scale Contamination
contaminate <- function(y) {
  m  <- max(1, round(0.10 * length(y)))
  id <- sample(seq_along(y), m, replace = FALSE)
  y[id] <- y[id] + 8 * sd(y) * sample(c(-1, 1), m, replace = TRUE)
  y
}

# Pinball (Check) Loss Function
QL <- function(y, yp, tau) {
  u <- y - yp
  mean(ifelse(u >= 0, tau * u, (tau - 1) * u))
}

# 3. Model Estimation Algorithms


# Unpenalized Quantile Regression (QR)
fitQR <- function(X, y, tau) {
  coef(rq(y ~ X, tau = tau, method = "fn"))
}

# Lasso Quantile Regression with Cross-Validation (LQR)
fitLQR <- function(X, y, tau) {
  lam <- seq(0.01, 1, length.out = 30)
  cv  <- sapply(lam, function(L) {
    z <- rq.fit.lasso(cbind(1, X), y, tau = tau, lambda = L)
    QL(y, cbind(1, X) %*% z$coefficients, tau)
  })
  L <- lam[which.min(cv)]
  rq.fit.lasso(cbind(1, X), y, tau = tau, lambda = L)$coefficients
}

# Frequentist Fused Lasso Quantile Regression (FLQR)
fitFLQR <- function(X, y, tau, l1 = 0.1, l2 = 0.1) {
  pp <- ncol(X)
  D  <- matrix(0, pp - 1, pp)
  for (j in 1:(pp - 1)) {
    D[j, j]     <- -1
    D[j, j + 1] <-  1
  }
  obj <- function(z) {
    e <- y - z[1] - X %*% z[-1]
    sum(ifelse(e >= 0, tau * e, (tau - 1) * e)) +
      l1 * sum(abs(z[-1])) +
      l2 * sum(abs(D %*% z[-1]))
  }
  optim(
    fitQR(X, y, tau),
    obj,
    method  = "Nelder-Mead",
    control = list(maxit = 10000, reltol = 1e-8)
  )$par
}

# Bayesian Fused Lasso Quantile Regression Gibbs Sampler (BFLQR)
fitBFLQR <- function(X, y, tau, iter = ITER, burn = BURN) {
  n  <- nrow(X)
  pp <- ncol(X)
  
  phi1 <- (1 - 2 * tau) / (tau * (1 - tau))
  phi2 <- 2 / (tau * (1 - tau))
  Z    <- cbind(1, X)
  
  D <- matrix(0, pp - 1, pp)
  for (j in 1:(pp - 1)) {
    D[j, j]     <- -1
    D[j, j + 1] <-  1
  }
  
  beta0 <- median(y)
  beta  <- rep(0, pp)
  sigma <- 1
  
  v      <- rep(1, n)
  tau2   <- rep(1, pp)
  omega2 <- rep(1, pp - 1)
  
  lambda1sq <- 1
  lambda2sq <- 1
  
  K  <- iter - burn
  S0 <- numeric(K)
  S1 <- matrix(0, K, pp)
  
  for (m in 1:iter) {
    r   <- as.numeric(y - beta0 - X %*% beta)
    chi <- pmax(r^2 / (phi2 * sigma), 1e-12)
    psi <- (phi1^2 / phi2 + 2) / sigma
    
    for (i in 1:n) {
      v[i] <- rgig(1, 0.5, chi[i], psi)
    }
    
    Qbeta <- diag(1 / pmax(tau2, 1e-12)) +
      t(D) %*% diag(1 / pmax(omega2, 1e-12)) %*% D
    
    Qhat <- matrix(0, pp + 1, pp + 1)
    Qhat[2:(pp + 1), 2:(pp + 1)] <- Qbeta
    
    Prec     <- (1 / phi2) * crossprod(Z, Z / v) + Qhat + 1e-8 * diag(pp + 1)
    Sigma.xi <- sigma * solve(Prec)
    
    mu.xi <- solve(Prec, (1 / phi2) * crossprod(Z, (y - phi1 * v) / v))
    
    xi    <- as.numeric(mvrnorm(1, mu.xi, Sigma.xi))
    beta0 <- xi[1]
    beta  <- xi[-1]
    
    for (j in 1:pp) {
      tau2[j] <- 1 / rgig(1, -0.5, max(beta[j]^2 / sigma, 1e-12), lambda1sq)
    }
    
    d <- as.numeric(D %*% beta)
    
    for (j in 1:(pp - 1)) {
      omega2[j] <- 1 / rgig(1, -0.5, max(d[j]^2 / sigma, 1e-12), lambda2sq)
    }
    
    lambda1sq <- rgamma(1, 1 + pp, rate = 1 + 0.5 * sum(tau2))
    lambda2sq <- rgamma(1, 1 + pp - 1, rate = 1 + 0.5 * sum(omega2))
    
    delta  <- y - beta0 - X %*% beta - phi1 * v
    a.star <- 1 + (3 * n + pp) / 2
    b.star <- 1 + sum(delta^2 / (2 * phi2 * v)) + sum(v) +
      0.5 * as.numeric(t(beta) %*% Qbeta %*% beta)
    
    sigma <- 1 / rgamma(1, a.star, rate = b.star)
    
    if (m > burn) {
      k <- m - burn
      S0[k]    <- beta0
      S1[k, ]  <- beta
    }
  }
  
  list(
    coef  = c(mean(S0), colMeans(S1)),
    draws = S1
  )
}

# 4. Variable Selection Performance Metrics

# Bayesian 95% Posterior Credible Interval Selection Rule
selectionMetrics <- function(draws) {
  lower    <- apply(draws, 2, quantile, 0.025)
  upper    <- apply(draws, 2, quantile, 0.975)
  selected <- (lower > 0 | upper < 0)
  
  TP <- sum(selected & true_active)
  TN <- sum(!selected & true_zero)
  FP <- sum(selected & true_zero)
  FN <- sum(!selected & true_active)
  
  list(
    Selection_Accuracy    = 100 * (TP + TN) / p,
    True_Zero_Recovery    = 100 * TN / sum(true_zero),
    True_Nonzero_Recovery = 100 * TP / sum(true_active),
    False_Positive_Rate   = 100 * FP / sum(true_zero),
    False_Negative_Rate   = 100 * FN / sum(true_active)
  )
}

# Absolute Threshold Selection Rule for Frequentist Penalized Models
selectionFromCoef <- function(coef) {
  beta     <- coef[-1]
  selected <- abs(beta) > 1e-6
  
  TP <- sum(selected & true_active)
  TN <- sum(!selected & true_zero)
  FP <- sum(selected & true_zero)
  FN <- sum(!selected & true_active)
  
  list(
    Selection_Accuracy    = 100 * (TP + TN) / p,
    True_Zero_Recovery    = 100 * TN / sum(true_zero),
    True_Nonzero_Recovery = 100 * TP / sum(true_active),
    False_Positive_Rate   = 100 * FP / sum(true_zero),
    False_Negative_Rate   = 100 * FN / sum(true_active)
  )
}

# 5. Single Monte Carlo Replication Execution

oneRun <- function(n, tau, r) {
  # Clean Training Data & Test Sets
  tr <- gendata(n, tau)
  te <- gendata(n, tau)
  
  # Contaminated Training Set (Same Predictors X, Contaminated Response y)
  tro   <- tr
  tro$y <- contaminate(tr$y)
  
  # Clean Model Fits
  q <- fitQR(tr$X, tr$y, tau)
  l <- fitLQR(tr$X, tr$y, tau)
  f <- fitFLQR(tr$X, tr$y, tau)
  b <- fitBFLQR(tr$X, tr$y, tau)
  
  # Contaminated Model Fits
  qo <- fitQR(tro$X, tro$y, tau)
  lo <- fitLQR(tro$X, tro$y, tau)
  fo <- fitFLQR(tro$X, tro$y, tau)
  bo <- fitBFLQR(tro$X, tro$y, tau)
  
  # Parameter MSE
  MSE <- function(est) mean((est[-1] - bt)^2)
  mse <- c(MSE(q), MSE(l), MSE(f), MSE(b$coef))
  
  # Variable Selection Performance
  QRsel <- list(
    Selection_Accuracy    = NA,
    True_Zero_Recovery    = NA,
    True_Nonzero_Recovery = NA,
    False_Positive_Rate   = NA,
    False_Negative_Rate   = NA
  )
  
  LQRsel   <- selectionFromCoef(l)
  FLQRsel  <- selectionFromCoef(f)
  BFLQRsel <- selectionMetrics(b$draws)
  
  sel <- rbind(
    unlist(QRsel),
    unlist(LQRsel),
    unlist(FLQRsel),
    unlist(BFLQRsel)
  )
  
  # Clean Test Loss
  pred <- list(
    QR    = cbind(1, te$X) %*% q,
    LQR   = cbind(1, te$X) %*% l,
    FLQR  = cbind(1, te$X) %*% f,
    BFLQR = cbind(1, te$X) %*% b$coef
  )
  
  cleanQL <- c(
    QL(te$y, pred$QR, tau),
    QL(te$y, pred$LQR, tau),
    QL(te$y, pred$FLQR, tau),
    QL(te$y, pred$BFLQR, tau)
  )
  
  # Contaminated-Training Test Loss
  pred.o <- list(
    QR    = cbind(1, te$X) %*% qo,
    LQR   = cbind(1, te$X) %*% lo,
    FLQR  = cbind(1, te$X) %*% fo,
    BFLQR = cbind(1, te$X) %*% bo$coef
  )
  
  outQL <- c(
    QL(te$y, pred.o$QR, tau),
    QL(te$y, pred.o$LQR, tau),
    QL(te$y, pred.o$FLQR, tau),
    QL(te$y, pred.o$BFLQR, tau)
  )
  
  deterioration <- 100 * (outQL - cleanQL) / cleanQL
  
  data.frame(
    n                             = n,
    tau                           = tau,
    Replication                   = r,
    Model                         = c("QR", "LQR", "FLQR", "BFLQR"),
    Parameter_Recovery_MSE        = mse,
    Selection_Accuracy            = sel[, 1],
    True_Zero_Recovery            = sel[, 2],
    True_Nonzero_Recovery         = sel[, 3],
    False_Positive_Rate           = sel[, 4],
    False_Negative_Rate           = sel[, 5],
    Clean_Test_QL                 = cleanQL,
    Outlier_Deterioration_Percent = deterioration
  )
}

# 6. Execute Full Simulation Loop

ALL     <- list()
k       <- 1
TOTAL   <- length(N) * length(TAU) * R
counter <- 0

for (n in N) {
  for (tau in TAU) {
    for (r in 1:R) {
      counter <- counter + 1
      
      cat(sprintf(
        "Progress: %d/%d | n=%d | tau=%.2f | rep=%d\n",
        counter, TOTAL, n, tau, r
      ))
      
      ALL[[k]] <- tryCatch(
        oneRun(n, tau, r),
        error = function(e) {
          message(
            "ERROR: n=", n, " tau=", tau, " rep=", r, " : ", e$message
          )
          NULL
        }
      )
      k <- k + 1
    }
  }
}

ALL <- ALL[!sapply(ALL, is.null)]

if (length(ALL) == 0) {
  stop("No successful replications were completed.")
}

SIM <- do.call(rbind, ALL)

# 7. Validation & Completeness Diagnostics

EXPECTED <- length(N) * length(TAU) * R
ACTUAL   <- length(unique(paste(SIM$n, SIM$tau, SIM$Replication, sep = "_")))

cat("\n\n")
cat("SIMULATION COMPLETENESS CHECK\n")
cat("\n")
cat("Expected replications: ", EXPECTED, "\n", sep = "")
cat("Completed replications: ", ACTUAL, "\n", sep = "")

if (ACTUAL < EXPECTED) {
  warning("Some replications failed. DO NOT use these results as final results.")
} else {
  cat("All requested replications completed successfully.\n")
}


# 8. Aggregation & Export

FINAL <- aggregate(
  cbind(
    Parameter_Recovery_MSE,
    Selection_Accuracy,
    True_Zero_Recovery,
    True_Nonzero_Recovery,
    False_Positive_Rate,
    False_Negative_Rate,
    Clean_Test_QL,
    Outlier_Deterioration_Percent
  ) ~ n + tau + Model,
  data = SIM,
  FUN  = function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
)

FINAL$Model <- factor(FINAL$Model, levels = c("QR", "LQR", "FLQR", "BFLQR"))
FINAL       <- FINAL[order(FINAL$n, FINAL$tau, FINAL$Model), ]

numcols <- c(
  "Parameter_Recovery_MSE",
  "Selection_Accuracy",
  "True_Zero_Recovery",
  "True_Nonzero_Recovery",
  "False_Positive_Rate",
  "False_Negative_Rate",
  "Clean_Test_QL",
  "Outlier_Deterioration_Percent"
)

FINAL[numcols] <- lapply(FINAL[numcols], function(x) round(x, 4))

# Filtered Subset Tables
T50  <- FINAL[FINAL$n == 50, ]
T100 <- FINAL[FINAL$n == 100, ]
T200 <- FINAL[FINAL$n == 200, ]

cat("\n================ TABLE: n = 50 ================\n")
print(T50, row.names = FALSE)

cat("\n================ TABLE: n = 100 ================\n")
print(T100, row.names = FALSE)

cat("\n================ TABLE: n = 200 ================\n")
print(T200, row.names = FALSE)

cat("\n================ COMPLETE TABLE ================\n")
print(FINAL, row.names = FALSE)

# Export Datasets for Submission
write.csv(T50,   "Model_Performance_n50.csv",                   row.names = FALSE)
write.csv(T100,  "Model_Performance_n100.csv",                  row.names = FALSE)
write.csv(T200,  "Model_Performance_n200.csv",                  row.names = FALSE)
write.csv(FINAL, "Complete_Model_Performance_All_Models.csv", row.names = FALSE)


# 9. Simulation Design Summary Output

cat("\n\n")
cat("SIMULATION DESIGN SUMMARY\n")
cat("\n")
cat("Sample sizes: n = 50, 100, 200\n")
cat("Quantiles: tau = 0.25, 0.50, 0.75\n")
cat("Replications per scenario: R = ", R, "\n", sep = "")
cat("Total scenarios: ", length(N) * length(TAU), "\n", sep = "")
cat("Total replication runs: ", length(N) * length(TAU) * R, "\n", sep = "")
cat("Predictors: p = 10\n")
cat("True-zero variables: X4, X5, X8\n")
cat("True-nonzero variables: X1, X2, X3, X6, X7, X9, X10\n")
cat("Outlier contamination: 10%\n")
cat("Outlier magnitude: 8 standard deviations\n")
cat("MCMC iterations: ", ITER, "\n", sep = "")
cat("Burn-in: ", BURN, "\n", sep = "")
cat("\n")