# ---- Load Required Libraries ----
library(readxl)       
library(forecast)     
library(tseries)      
library(rugarch)     
library(fBasics)      
library(FinTS)       
library(PerformanceAnalytics)
library(MTS)
library(zoo)
library(ggplot2)
library(dplyr)

# ---- Load and Prepare Data ----
setwd("~/data")
da <- read_excel("dataset_bk.xlsx", col_names = TRUE)
dates <- as.Date(da$`Date`, format = "%Y-%m-%d")
st <- zoo(da$ST, order.by = dates)

# Remove any zero or negative values before log transformation
st_clean <- st[st > 0 & !is.na(st)]
logchg <- na.omit(diff(log(st_clean)))  # log returns
logchg_vec <- coredata(logchg)

# ---- Enhanced Summary Statistics ----
cat("=== DESCRIPTIVE STATISTICS ===\n")
print(basicStats(logchg_vec))
cat("\nMean test (H0: mean = 0):\n")
print(t.test(logchg_vec))

# ---- Time Series Analysis ----
cat("\n=== TIME SERIES DIAGNOSTICS ===\n")
windows()
par(mfrow=c(2,2))
plot(logchg, main="Log Returns Time Series", ylab="Log Returns")
Acf(logchg_vec, main="ACF of Log Returns")
Pacf(logchg_vec, main="PACF of Log Returns")
qqnorm(logchg_vec); qqline(logchg_vec)

# Serial correlation tests
ljung_result <- Box.test(logchg_vec, lag = 20, type = "Ljung")
cat("Ljung-Box test (serial correlation):", ljung_result$p.value, "\n")

# ---- ARIMA MODELING ----
cat("\n=== ARIMA MODEL SELECTION ===\n")
st_log <- coredata(logchg)

# Compare multiple models
m1 <- auto.arima(st_log, ic = "aic", stepwise = FALSE, approximation = FALSE)
m2 <- auto.arima(st_log, ic = "bic", stepwise = FALSE, approximation = FALSE)

cat("AIC Model:", m1$arma, " AIC:", m1$aic, "\n")
cat("BIC Model:", m2$arma, " BIC:", m2$bic, "\n")

# Select best model (lower AIC generally preferred for forecasting)
best_model <- m1
summary(best_model)

# ---- Enhanced Residual Diagnostics ----
cat("\n=== RESIDUAL DIAGNOSTICS ===\n")
residuals_std <- residuals(best_model)

# Comprehensive residual tests
ljung_resid <- Box.test(residuals_std, lag = 20, type = "Ljung")
ljung_abs <- Box.test(abs(residuals_std), lag = 20, type = "Ljung")
ljung_sq <- Box.test(residuals_std^2, lag = 20, type = "Ljung")
jb_test <- jarqueberaTest(residuals_std)
arch_test <- archTest(residuals_std, lag = 20)

cat("Ljung-Box (levels):", ljung_resid$p.value, "\n")
cat("Ljung-Box (absolute):", ljung_abs$p.value, "\n")  
cat("Ljung-Box (squared):", ljung_sq$p.value, "\n")
cat("Jarque-Bera normality:", jb_test@test$p.value, "\n")
cat("ARCH test:", arch_test$p.value, "\n")

# Check for GARCH effects - with error handling
use_garch <- FALSE
arch_p_value <- NA

# Safe extraction of ARCH test p-value
tryCatch({
  if("ArchTest" %in% class(arch_test)) {
    arch_p_value <- arch_test$p.value
  } else if(is.list(arch_test) && "p.value" %in% names(arch_test)) {
    arch_p_value <- arch_test$p.value
  } else {
    # Alternative ARCH test using Box-Ljung on squared residuals
    arch_alt <- Box.test(residuals_std^2, lag = 12, type = "Ljung")
    arch_p_value <- arch_alt$p.value
    cat("Using alternative ARCH test (Ljung-Box on squared residuals)\n")
  }
}, error = function(e) {
  cat("Error in ARCH test, using alternative method\n")
  arch_alt <- Box.test(residuals_std^2, lag = 12, type = "Ljung")
  arch_p_value <- arch_alt$p.value
})

cat("ARCH test p-value:", arch_p_value, "\n")

# Proceed with GARCH if ARCH effects detected
if(!is.na(arch_p_value) && arch_p_value < 0.05) {
  cat("ARCH effects detected. Fitting GARCH model...\n")
  
  tryCatch({
    # Fit GARCH model
    arima_order <- arimaorder(best_model)
    spec <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(arima_order[1], arima_order[3]), include.mean = TRUE)
    )
    garch_fit <- ugarchfit(spec, st_log, solver = "hybrid")
    
    # Check if GARCH fit was successful
    if(convergence(garch_fit) == 0) {
      residuals_std <- residuals(garch_fit, standardize = TRUE)
      sigma_forecast <- sigma(garch_fit)
      last_sigma <- tail(sigma_forecast, 1)
      use_garch <- TRUE
      cat("GARCH model fitted successfully.\n")
    } else {
      cat("GARCH model failed to converge. Using ARIMA only.\n")
      use_garch <- FALSE
      last_sigma <- sd(residuals_std)
    }
  }, error = function(e) {
    cat("Error fitting GARCH model:", e$message, "\nUsing ARIMA only.\n")
    use_garch <- FALSE
    last_sigma <- sd(residuals_std)
  })
} else {
  cat("No significant ARCH effects detected. Using ARIMA only.\n")
  use_garch <- FALSE
  last_sigma <- sd(residuals_std)
}

# ---- Improved Monte Carlo Simulation ----
cat("\n=== MONTE CARLO SIMULATION ===\n")

# Parameters
n_simulations <- 100000
n_quarters <- 12  # 3 years of quarterly data
last_value <- tail(coredata(st_clean), 1)  # Use actual ST value, not log return

# Extract model parameters safely
mu <- mean(st_log, na.rm = TRUE)
sigma <- sd(residuals_std, na.rm = TRUE)

if(use_garch) {
  # Extract GARCH parameters
  tryCatch({
    mu <- fitted(garch_fit)[length(fitted(garch_fit))]
    model_order <- arimaorder(best_model)[c(1,3)]  # AR and MA orders
  }, error = function(e) {
    cat("Error extracting GARCH parameters, using ARIMA only\n")
    use_garch <<- FALSE
  })
} 

# Get ARIMA order safely
arima_order <- tryCatch({
  arimaorder(best_model)[c(1,3)]  # c(AR, MA)
}, error = function(e) {
  c(0, 0)  # Default to white noise if extraction fails
})

# Enhanced simulation function with better error handling
generate_path <- function() {
  tryCatch({
    if(use_garch && exists("garch_fit")) {
      # Use GARCH for volatility forecasting
      garch_sim <- ugarchsim(garch_fit, n.sim = n_quarters, m.sim = 1)
      log_returns <- as.numeric(fitted(garch_sim))
    } else {
      # Use ARIMA simulation or simple approach
      if(arima_order[2] > 0 && !is.null(coef(best_model))) {
        # Try to simulate from ARIMA model
        sim_result <- tryCatch({
          simulate(best_model, nsim = n_quarters)
        }, error = function(e) {
          # Fallback to simple random generation
          rnorm(n_quarters, mu, sigma)
        })
        log_returns <- as.numeric(sim_result)
      } else {
        # Simple random walk
        log_returns <- rnorm(n_quarters, mu, sigma)
      }
    }
    
    # Ensure log_returns is valid
    if(any(is.na(log_returns)) || length(log_returns) != n_quarters) {
      log_returns <- rnorm(n_quarters, mu, sigma)
    }
    
    # Convert log returns to price path
    price_path <- last_value * cumprod(exp(log_returns))
    return(tail(price_path, 1))
    
  }, error = function(e) {
    # Ultimate fallback
    log_returns <- rnorm(n_quarters, mu, sigma)
    price_path <- last_value * cumprod(exp(log_returns))
    return(tail(price_path, 1))
  })
}

# Run simulation
cat("Running", n_simulations, "simulations...\n")
set.seed(123)  # For reproducibility
mc_results <- replicate(n_simulations, generate_path())

# ---- Results Analysis ----
cat("\n=== SIMULATION RESULTS ===\n")
cat("Starting Value:", last_value, "\n")
cat("Mean Final Value:", mean(mc_results), "\n")
cat("Median Final Value:", median(mc_results), "\n")
cat("Standard Deviation:", sd(mc_results), "\n")

# Risk metrics (adjusted for short-term borrowings where higher = worse)
percentiles <- quantile(mc_results, probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99))
cat("\nPercentiles:\n")
print(percentiles)

# Value at Risk (VaR) - Upper tail risk for short-term borrowings
var_95 <- quantile(mc_results, 0.95)  # 95th percentile (upper tail)
var_99 <- quantile(mc_results, 0.99)  # 99th percentile (upper tail)
cat("\nValue at Risk - 95% worst case (upper 5%):", var_95, "\n")
cat("Value at Risk - 99% worst case (upper 1%):", var_99, "\n")

# Expected Shortfall (Conditional VaR) - Mean of worst outcomes
es_95 <- mean(mc_results[mc_results >= var_95])
es_99 <- mean(mc_results[mc_results >= var_99])
cat("Expected Shortfall - 95% level:", es_95, "\n")
cat("Expected Shortfall - 99% level:", es_99, "\n")

# Risk-adjusted probabilities
prob_increase <- mean(mc_results > last_value)
prob_stress_increase <- mean(mc_results > last_value * 1.5)  # 50% increase threshold
prob_severe_stress <- mean(mc_results > last_value * 2.0)    # 100% increase threshold

cat("Probability of increase (deterioration):", prob_increase, "\n")
cat("Probability of >50% increase:", prob_stress_increase, "\n")
cat("Probability of >100% increase (crisis):", prob_severe_stress, "\n")

# Additional risk metrics specific to borrowings
funding_stress_threshold <- last_value * 1.25  # 25% increase as stress threshold
prob_funding_stress <- mean(mc_results > funding_stress_threshold)
cat("Probability of funding stress (>25% increase):", prob_funding_stress, "\n")

# Calculate funding stability metrics
volatility_metric <- sd(mc_results) / mean(mc_results)  # Coefficient of variation
cat("Relative volatility (CV):", volatility_metric, "\n")

# ---- Visualization ----
windows()
par(mfrow=c(2,2))

# Histogram with risk-focused annotations
hist(mc_results, breaks = 100, main = "Distribution of Short-Term Borrowings (12Q Forecast)", 
     xlab = "Short-term Borrowings Level", col = "lightcyan", border = "white")
abline(v = last_value, col = "blue", lwd = 2, lty = 2)
abline(v = mean(mc_results), col = "darkgreen", lwd = 2)
abline(v = var_95, col = "orange", lwd = 2, lty = 3)
abline(v = var_99, col = "red", lwd = 2, lty = 3)
abline(v = funding_stress_threshold, col = "purple", lwd = 2, lty = 4)
legend("topright", c("Current Level", "Mean Forecast", "95% VaR", "99% VaR", "Stress Threshold"), 
       col = c("blue", "darkgreen", "orange", "red", "purple"), 
       lty = c(2, 1, 3, 3, 4), cex = 0.8)

# Q-Q plot
qqnorm(mc_results, main = "Q-Q Plot of Simulation Results")
qqline(mc_results)

# Time series of sample paths (first 100)
sample_paths <- matrix(NA, nrow = n_quarters, ncol = 100)
set.seed(456)  # Different seed for visualization

for(i in 1:100) {
  path_result <- tryCatch({
    if(use_garch && exists("garch_fit")) {
      garch_sim <- ugarchsim(garch_fit, n.sim = n_quarters, m.sim = 1)
      log_returns <- as.numeric(fitted(garch_sim))
    } else {
      # Use same logic as main simulation
      if(arima_order[2] > 0 && !is.null(coef(best_model))) {
        sim_result <- tryCatch({
          as.numeric(simulate(best_model, nsim = n_quarters))
        }, error = function(e) {
          rnorm(n_quarters, mu, sigma)
        })
        log_returns <- sim_result
      } else {
        log_returns <- rnorm(n_quarters, mu, sigma)
      }
    }
    
    # Ensure valid log returns
    if(any(is.na(log_returns)) || length(log_returns) != n_quarters) {
      log_returns <- rnorm(n_quarters, mu, sigma)
    }
    
    last_value * cumprod(exp(log_returns))
  }, error = function(e) {
    # Fallback path
    log_returns <- rnorm(n_quarters, mu, sigma)
    last_value * cumprod(exp(log_returns))
  })
  
  sample_paths[,i] <- path_result
}

matplot(1:n_quarters, sample_paths, type = "l", col = alpha("gray", 0.3), 
        main = "Sample Simulation Paths", xlab = "Quarter", ylab = "ST Value")
lines(1:n_quarters, rowMeans(sample_paths), col = "red", lwd = 2)

# Convergence plot
running_mean <- cumsum(mc_results) / (1:n_simulations)
plot(running_mean[1:min(10000, n_simulations)], type = "l", 
     main = "Convergence of Mean", xlab = "Simulation", ylab = "Running Mean")

cat("\n=== SIMULATION COMPLETE ===\n")
cat("Model used:", ifelse(use_garch, "ARIMA-GARCH", "ARIMA"), "\n")
cat("Simulation horizon:", n_quarters, "quarters\n")
cat("Number of simulations:", n_simulations, "\n")
cat("\n=== RISK ASSESSMENT SUMMARY ===\n")
cat("Current ST Borrowings Level:", last_value, "\n")
cat("Expected Level (12Q):", round(mean(mc_results), 2), "\n")
cat("Funding Risk (99% VaR):", round(var_99, 2), "\n")
cat("Crisis Probability (>100% increase):", round(prob_severe_stress * 100, 2), "%\n")
cat("Stress Probability (>25% increase):", round(prob_funding_stress * 100, 2), "%\n")