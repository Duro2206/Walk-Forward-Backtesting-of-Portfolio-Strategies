
rm(list=ls(all.names = TRUE))
graphics.off()

library(pacman)
p_load(riskParityPortfolio)
p_load(portfolioBacktest)
p_load(quantmod)
p_load(PerformanceAnalytics)
p_load(CVXR)
p_load(DT)
p_load(xts)
p_load(ggplot2)
p_load(scales)
p_load(corrplot)

from <- as.Date("2010-01-01")
to   <- as.Date("2025-08-31")
freq <- "daily"

Tickers <- c("MSFT","GS","CVX","COST","CAT","SAP","UL","ASML","SHEL","TLT","GLD","LQD")

# Prices for the universe
prices <- xts()
for (tk in Tickers) {
  px <- Ad(getSymbols(Symbols = tk, from = from, to = to, periodicity = freq, auto.assign = FALSE))
  px <- na.approx(px, na.rm = FALSE)
  prices <- cbind(prices, px)
}
colnames(prices) <- Tickers
tclass(prices) <- "Date"

# S&P 500 index (for benchmark = "index")
spx <- Ad(getSymbols("^GSPC", from = from, to = to, periodicity = freq, auto.assign = FALSE))
colnames(spx) <- "SP500"

# Wrap to portfolioBacktest dataset format, now with index
universe <- list(adjusted = prices, index = spx)

# ----------------------------
# Randomized 2-year datasets
# ----------------------------
set.seed(123)
dataset_list <- financialDataResample(
  universe,
  N_sample     = 12, 
  T_sample     = 252*2,                 
  num_datasets = 100
)

# ----------------------------
# Train/Test split helper
# ----------------------------
split_dataset_list <- function(ds_list, prop = 0.70) {
  trn <- list()
  tst <- list()
  for (nm in names(ds_list)) {
    ds <- ds_list[[nm]]
    A  <- ds$adjusted
    I  <- ds$index
    n  <- nrow(A)
    n_trn <- max(1, min(n-1, floor(n * prop)))
    trn[[nm]] <- list(adjusted = A[1:n_trn, ], index = I[1:n_trn, ])
    tst[[nm]] <- list(adjusted = A[(n_trn+1):n, ], index = I[(n_trn+1):n, ])
  }
  list(train = trn, test = tst)
}
splits <- split_dataset_list(dataset_list, prop = 0.70)
dataset_list_trn <- splits$train
dataset_list_tst <- splits$test

# ----------------------------
# Portfolio definitions 
# ----------------------------
# 1. Equally Weighted Porfolio
EWP_portfolio_fun <- function(dataset, ...) {
  X <- CalculateReturns(dataset$adjusted, "log")[-1] # computed the log returns
  N <- ncol(X)
  weights <- rep(1 /N, N)  # equal weights
  return(weights)
}

# 2. Global Minimum-Variance Portfolio with no short-selling
GMVP_portfolio_fun <- function(dataset, ...) {
  X <- CalculateReturns(dataset$adjusted,"log")[-1] # computed the log returns
  mu <- colMeans(X) # mean return
  Sigma <- cov(X) # covariance matrix
  ## Defining the Optimization variables and problem for convex optimization
  w <- Variable(nrow(Sigma)) # defining the optimization variables i.e. the weight of each asset
  objective <- Minimize(quad_form(w, Sigma)) # setting objective function
  constraints <- list(w >= 0, sum(w) == 1) # setting the constraints
  prob <- Problem(objective,constraints) # defining the problem
  result <- CVXR::solve(prob) # solve the optimization problem
  weights <- as.vector(result$getValue(w))
  return(weights)
}

# Sum must be equal to 1 ensures that entire portfolio is fully invested
sum(GMVP_portfolio_fun(dataset_list$`dataset 1`))

# GMVP's weight and risk distribution
w <- GMVP_portfolio_fun(dataset_list$`dataset 1`)
names(w) <- Tickers
X <- CalculateReturns(dataset_list$`dataset 1`$adjusted,"log")[-1]
Sigma <- cov(X)
colnames(Sigma) <- rownames(Sigma) <- Tickers
barplotPortfolioRisk(w,Sigma) + scale_y_continuous(labels = scales::percent) + 
  ggtitle("GMVP Portfolio Weight and Risk Distribution")

# 3. Markowitz’s mean-variance portfolio (MVP) with no short-selling
Markowitz_portfolio_fun <- function(dataset, lambda = 0.5, ...) {
  X <- CalculateReturns(dataset$adjusted,"log")[-1] # computed the log returns
  mu <- colMeans(X) # mean log return
  Sigma <- cov(X) # covariance matrix
  # design mean-variance portfolio
  w <- Variable(nrow(Sigma))
  objective <- Maximize(t(mu) %*% w - lambda*quad_form(w, Sigma))
  constraints = list(w >= 0, sum(w) == 1)
  prob <- Problem(objective,constraints)
  result <- solve(prob)
  weights <- as.vector(result$getValue(w))
  return(weights)
}

# Sum must be equal to 1 ensures that entire portfolio is fully invested
sum(GMVP_portfolio_fun(dataset_list$`dataset 1`))

# Markowitz’s Mean Variance Portfolio's weight and risk distribution
w <- Markowitz_portfolio_fun(dataset_list$`dataset 1`)
names(w) <- Tickers
X <- CalculateReturns(dataset_list$`dataset 1`$adjusted,"log")[-1]
Sigma <- cov(X)
colnames(Sigma) <- rownames(Sigma) <- Tickers
barplot(w, main="MVP Weights Allocation", names.arg=colnames(X), las=2)

# 4. Risk Parity Portfolio
RPP_portfolio_fun <- function(dataset, ...) {
  X <- CalculateReturns(dataset$adjusted, "log")[-1] # computed the log returns
  Sigma <- cov(X) # covariance matrix
  result <- riskParityPortfolio(Sigma)
  weights <- as.vector(result$w)
  return(weights)
}

# RPP's weight and risk distribution
w <- RPP_portfolio_fun(dataset_list$`dataset 1`)
names(w) <- Tickers
X <- CalculateReturns(dataset_list$`dataset 1`$adjusted,"log")[-1]
Sigma <- cov(X)
colnames(Sigma) <- rownames(Sigma) <- Tickers
barplotPortfolioRisk(w,Sigma) + scale_y_continuous(labels = scales::percent) + 
  ggtitle("Portfolio weight and risk distribution")

# 5. Inverse Volatility Portfolio
IVP_portfolio_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma <- cov(X)  # compute SCM
  # design IVP
  weights <- riskParityPortfolio(Sigma, formulation='diag')$w
  return(weights)
}

# IVP's weight and risk distribution
w <- IVP_portfolio_fun(dataset_list$`dataset 1`)
names(w) <- Tickers
X <- CalculateReturns(dataset_list$`dataset 1`$adjusted,"log")[-1]
Sigma <- cov(X)
colnames(Sigma) <- rownames(Sigma) <- Tickers
barplotPortfolioRisk(w,Sigma) + scale_y_continuous(labels = scales::percent) + 
  ggtitle("Portfolio weight and risk distribution")

# 6. Most Sharpe Ratio Portfolio
MSRP_portfolio_fun <- function(dataset, risk_free = 0, ...) {
  X <- CalculateReturns(dataset$adjusted, "log")[-1]
  mu <- colMeans(X)
  Sigma <- cov(X)
  w <- Variable(nrow(Sigma))
  objective <- Minimize(quad_form(w, Sigma))
  constraints <- list(t(mu - risk_free) %*% w == 1, w >= 0)
  prob <- Problem(objective, constraints)
  result <- CVXR::solve(prob)
  weights_raw <- as.vector(result$getValue(w))
  weights <- weights_raw / sum(weights_raw)   # normalize to sum = 1
  return(weights)
}

# MSRP's weight and risk distribution
w <- MSRP_portfolio_fun(dataset_list$`dataset 1`)
names(w) <- Tickers
X <- CalculateReturns(dataset_list$`dataset 1`$adjusted,"log")[-1]
Sigma <- cov(X)
colnames(Sigma) <- rownames(Sigma) <- Tickers
barplotPortfolioRisk(w,Sigma) + scale_y_continuous(labels = scales::percent) + 
  ggtitle("Portfolio weight and risk distribution")

# 7. Most Diverse Portfolio
MDP_portfolio_fun <- function(dataset, ...) {
  X <- CalculateReturns(dataset$adjusted, "log")[-1, ]  
  Sigma <- cov(X)
  n <- ncol(Sigma)
  w <- Variable(n)
  objective <- Minimize(quad_form(w, Sigma))
  constraints <- list(w >= 0, sum(w) == 1)
  prob <- Problem(objective, constraints)
  result <- CVXR::solve(prob)
  weights <- as.vector(result$getValue(w))
  weights <- weights / sum(weights)
  return(weights)
}

# MDP's weight and risk distribution
w <- MDP_portfolio_fun(dataset_list$`dataset 1`)
names(w) <- Tickers
X <- CalculateReturns(dataset_list$`dataset 1`$adjusted,"log")[-1]
Sigma <- cov(X)
colnames(Sigma) <- rownames(Sigma) <- Tickers
barplotPortfolioRisk(w, Sigma) + scale_y_continuous(labels=scales::percent) + 
  ggtitle("Most Diversified Portfolio Weight and Risk Distribution")

# 8. Maximum de-correlation portfolio
MDCP_portfolio_fun <- function(dataset, ...) {
  X <- CalculateReturns(dataset$adjusted, "log")[-1]  # compute log returns
  Sigma <- cov(X)
  # Correlation matrix
  C <- diag(1 / sqrt(diag(Sigma))) %*% Sigma %*% diag(1 / sqrt(diag(Sigma)))
  # Define optimization variables
  w <- Variable(nrow(C))
  objective <- Minimize(quad_form(w, C))
  constraints <- list(sum(w) == 1, w >= 0)
  prob <- Problem(objective, constraints)
  # Solve optimization
  result <- try(CVXR::solve(prob), silent = TRUE)
  # Handle solver errors
  if (inherits(result, "try-error") || is.null(result$getValue(w))) {
    warning("MDCP optimization failed, returning equal weights")
    n <- nrow(C)
    return(rep(1/n, n))
  }
  weights <- as.vector(result$getValue(w))
  weights <- weights / sum(weights)  # normalize just in case
  return(weights)
}

# Maximum De-correlation Portfolio's weight and risk distribution
w <- MDCP_portfolio_fun(dataset_list$`dataset 1`)
names(w) <- Tickers
X <- CalculateReturns(dataset_list$`dataset 1`$adjusted,"log")[-1]
Sigma <- cov(X)
C <- diag(1 / sqrt(diag(Sigma))) %*% Sigma %*% diag(1 / sqrt(diag(Sigma)))
colnames(C) <- Tickers
rownames(C) <- Tickers
colnames(Sigma) <- rownames(Sigma) <- Tickers
  # Correlation matrix heatmap (show decorrelation)
corrplot(C, method="color", title = "Correlation Matrix for MDCP", mar=c(0,0,1,0))
  # Weight and Risk barplot
barplotPortfolioRisk(w, Sigma) + scale_y_continuous(labels=scales::percent) + 
  ggtitle("Max Decorrelation Portfolio Weight and Risk Distribution")



# Portfolio list for backtesting
portfolios <- list(
  "EWP"       = EWP_portfolio_fun,
  "GMVP"      = GMVP_portfolio_fun,
  "Markowitz" = Markowitz_portfolio_fun,
  "RPP"       = RPP_portfolio_fun,
  "IVP"       = IVP_portfolio_fun,
  "MSRP"      = MSRP_portfolio_fun,
  "MDP"       = MDP_portfolio_fun,
  "MDCP"      = MDCP_portfolio_fun
)

names(portfolios)
# ----------------------------
# Walk-forward backtests: train and test
# ----------------------------
common_args <- list(
  portfolio_funs      = portfolios,
  benchmark       = c("index"),
  rebalance_every = 21,     
  optimize_every  = 21,     
  lookback        = 126,  
  shortselling    = FALSE,
  show_progress_bar=TRUE,
  cost = list(buy = 0e-4, sell = 0e-4, short = 0e-4, long_leverage = 0e-4)
)

bt_trn <- do.call(portfolioBacktest, c(list(dataset_list = dataset_list_trn), common_args))
bt_tst <- do.call(portfolioBacktest, c(list(dataset_list = dataset_list_tst), common_args))

names(bt_trn)
names(bt_tst)

# ----------------------------
# Reporting: train vs test
# ----------------------------
res_trn <- backtestSummary(bt_trn, summary_fun = median, show_benchmark = TRUE)
res_tst <- backtestSummary(bt_tst, summary_fun = median, show_benchmark = TRUE)

# Interactive tables
summaryTable(res_trn, type = "DT", order_col = "Sharpe ratio", order_dir = "desc")
summaryTable(res_tst, type = "DT", order_col = "Sharpe ratio", order_dir = "desc")

# Barplots and boxplots
summaryBarPlot(res_trn, measures = c("Sharpe ratio", "max drawdown"))
summaryBarPlot(res_tst, measures = c("Sharpe ratio", "max drawdown"))
backtestBoxPlot(bt_trn, measure = "Sharpe ratio", type = c("ggplot2"))
backtestBoxPlot(bt_tst, measure = "Sharpe ratio", type = c("ggplot2"))

# Cumulative return charts against the Benchmark inside each split
backtestChartCumReturn(bt_trn, c("EWP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_tst, c("EWP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_trn, c("GMVP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_tst, c("GMVP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_trn, c("Markowitz", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_tst, c("Markowitz", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_trn, c("RPP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_tst, c("RPP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_trn, c("IVP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_tst, c("IVP", "index"), type = c("ggplot2"))
backtestChartCumReturn(bt_trn, c("MSRP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_tst, c("MSRP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_trn, c("MDP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_tst, c("MDP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_trn, c("MDCP","index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_tst, c("MDCP","index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_trn, c("EWP", "GMVP", "Markowitz", "RPP", "IVP", "MSRP", "MDP", "MDCP","index"), type = c("ggplot2"), dataset_num=1)
backtestChartCumReturn(bt_tst, c("EWP", "GMVP", "Markowitz", "RPP", "IVP", "MSRP", "MDP", "MDCP","index"), type = c("ggplot2"), dataset_num=1)

# Drawdown charts against the Benchmark inside each split
backtestChartDrawdown(bt_trn, c("EWP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdown(bt_tst, c("EWP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdown(bt_trn, c("GMVP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdown(bt_tst, c("GMVP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdown(bt_trn, c("Markowitz", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdown(bt_tst, c("Markowitz", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdownn(bt_trn, c("RPP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdownn(bt_tst, c("RPP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdown(bt_trn, c("IVP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdown(bt_tst, c("IVP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdown(bt_trn, c("MSRP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdown(bt_tst, c("MSRP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdown(bt_trn, c("MDP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdown(bt_tst, c("MDP", "index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdown(bt_trn, c("MDCP","index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdown(bt_tst, c("MDCP","index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdown(bt_trn, c("EWP", "GMVP", "Markowitz", "RPP", "IVP", "MSRP", "MDP", "MDCP","index"), type = c("ggplot2"), dataset_num=1)
backtestChartDrawdown(bt_tst, c("EWP", "GMVP", "Markowitz", "RPP", "IVP", "MSRP", "MDP", "MDCP","index"), type = c("ggplot2"), dataset_num=1)

