#
# This script presents calculation of the option price using Monte Carlo
# and exact formula if dynamics of underlying is ruled by square-root diffusion process
# 

source("required_libraries.R")

#
# Model parameters
#
V0 = 17.5  # initial level of the process
kappa_V = 0.1  
theta_V = 20.0  
sigma_V = 2.0  
zeta_V = 0.0  
r = 0.01  # risk-free interest rate

#
# Option parameters
#
K = 20.0  # strike
T = 1.0  # time horizon, maturity

#
# Simulation and MC parameters
#
M = 50 # number of time points
I= 50000 # number of trajectories
x0 = V0 # starting value of the process

#
# Formula for futures valuation
#
futuresPrice <-function(V0, kappa_V, theta_V, zeta_V, T){
  
  # Futures pricing formula in GL96 model
  # 
  # Arguments:
  #
  #     V0:       current process level
  #     
  #     kappa_V:  mean-reversion factor
  #     
  #     theta_V:  long-run mean of the process
  #     
  #     zeta_V:   volatility risk premium
  #     
  #     T:        time-to-maturity
  #
  #
  # Output:
  #
  #   price of the future 
  
  
  alpha = kappa_V * theta_V
  beta = kappa_V + zeta_V
  price = (alpha / beta * (1 - exp(-beta * T)) + exp(-beta * T) * V0)
  return(price)
  
}

# performance of futuresPrice function
# futuresPrice(V0, kappa_V, theta_V, zeta_V, T)

#
# non-central chi-squared 
#
cx <- function(K, gamma, nu, lambda_V, exact=True){
  
  #  Complementary distribution function of non-central chi-squared density.
  #
  # Args:
  #
  #       K: strike price
  # 
  #       gamma: as defined in the GL96 model
  # 
  #       nu: degrees of freedom
  # 
  #       lambda_V: non-centrality parameter
  # 
  
  out <- (1 - pchisq(gamma * K, nu, lambda_V))
  return(out) 
  
}

#
# exact call price
#
callPriceGL96 <- function(V0, kappa_V, theta_V, sigma_V, zeta_V, T, r, K){
  
  # Call option pricing formula 
  #
  # Args:
  # 
  #     V0: current process level
  # 
  #     kappa_V: mean-reversion factor
  # 
  #     theta_V: long-run mean of the process
  # 
  #     sigma_V: volatility of process
  # 
  #     zeta_V: volatility risk premium
  # 
  #     T: time-to-maturity
  # 
  #     r: risk-free short rate
  # 
  #     K: strike price of the option
  # 
  #
  # Output:
  # 
  #     price of the call option based on the exact formula
  
  D = exp(-r * T)  # discount factor
  
  alpha = kappa_V * theta_V 
  beta = kappa_V + zeta_V
  gamma = 4 * beta / (sigma_V ** 2 * (1 - exp(-beta * T)))
  nu = 4 * alpha / sigma_V ** 2
  lambda_V = gamma * exp(-beta * T) * V0
  
  # the pricing formula
  call = (D * exp(-beta * T) * V0 * cx(K, gamma, nu + 4, lambda_V)
          + D * (alpha / beta) * (1 - exp(-beta * T))
          * cx(K, gamma, nu + 2, lambda_V)
          - D * K * cx(K, gamma, nu, lambda_V))
  return (call)
  
}

# price of the call using closed-form expression
# callPriceGL96(V0, kappa_V, theta_V, sigma_V, zeta_V, T, r, K)

generatePaths <- function(x0, kappa_V, theta_V, sigma_V, zeta_V, T, M, I, r, K){
  
  # 
  # Args:
  # 
  #     x0: starting value of the process to perform simulations
  # 
  #     kappa_V: mean-reversion factor
  # 
  #     theta_V: long-run mean 
  # 
  #     sigma_V: volatility 
  # 
  #     zeta_V: volatility risk premium, used to switch to risk-neutral world 
  # 
  #     T: horizon of simulations, and maturity 
  # 
  #     r: risk-free short rate
  # 
  #     K: strike price of the option
  #
  #     M: number of time points at which values of the process must be calculated
  #
  #     I: number of trajectories of the process to be calculated
  # 
  #
  # Output:
  # 
  #     list object, first element: matrix with simulated trajectories;
  #                  second element: price of the call option calculated using
  #                                  Monte Carlo
  
  # switch to risk-neutral process
  kappa_V <- kappa_V + zeta_V
  
  dt = T / M
  
  x <- matrix(0, M+1, I)
  x[1,] <- x0
  rv <- matrix( rnorm(I*(M+1)), M+1, I) 
  
  d = 4 * kappa_V * theta_V / (sigma_V^2)
  c = (sigma_V^2 * (1 - exp(-kappa_V * dt))) / (4 * kappa_V)
  
  if (d>1){
    
    for (t in 1:M){
      l <- x[t,] * exp(-kappa_V * dt)/c
      chi <- rchisq(I, d-1)
      x[t+1,] <- c*((rv[t+1,] + l^(0.5))^2 + chi)
    }
    
  } else {
    
    for (t in 1:(M)){
      l <- x[t, ] * exp(-kappa_V * dt) / c
      N <-  rpois(I, l/2)
      chi <- rchisq(I, d + 2 * N)
      x[t+1] = c * chi
      
    }
  }
  
  # calculate Monte Carlo price of the call based on generated trajectories
  MCprice <- exp(-r * T) * sum( pmax(as.vector(x[nrow(x),] - K), 
                                     rep(0,ncol(x)) ) ) / I
  # create a list containing generated trajectories, and Monte Carlo price
  # of the option
  out <- list(generatedTrajectories = x,
              MCpriceCall = MCprice) 
  return(out)
} 

# generatePaths(x0, kappa_V, theta_V, sigma_V, zeta_V, T, M, I, r, K)$MCpriceCall


### main part ##################################################################

#
# Typical trajectory of the square-root diffusion process
#

# calculate trajectory
trajectories <- data.frame(generatePaths(x0, kappa_V, theta_V, sigma_V, zeta_V, 
                         T, 10000, 1, r, K)$generatedTrajectories)
colnames(trajectories)<-"process"
timeTraj <- seq(0,1, length.out = 10001)
trajPlot <- cbind(timeTraj, trajectories)

# plotting procedure
ggplot(trajPlot, aes(timeTraj, process)) + 
  geom_line(size=0.3) + 
  xlab("time") + 
  ylab("") + 
  ggtitle("") +
  ggsave("squareRoot.pdf", width = 6, height =4)


#
# Evaluation of option price using MC and exact formula ------------------------
#

# range of strikes to compute option prices
strikes = seq(15,26, length.out = 20)

optionPriceFormula <- rep(0,length(strikes))
optionPriceMC <- rep(0,length(strikes))

for(i in 1:length(strikes)){
  
  optionPriceFormula[i] <- callPriceGL96(V0, kappa_V, theta_V, sigma_V, zeta_V, 
                                         T, r, strikes[i])
  optionPriceMC[i] <- generatePaths(x0, kappa_V, theta_V, sigma_V, zeta_V, 
                                    T, M, I, r, strikes[i])$MCpriceCall
} 

# option prices based on MC and formula
dataToPlot <- data.frame(strikes, optionPriceFormula, optionPriceMC )

# mean absolute error
MAE <- mean(abs(optionPriceFormula-optionPriceMC))

# plotting procedure
ggplot(data=dataToPlot, aes(x=strikes))+
  geom_line(aes(y=optionPriceFormula, colour ="using formula"))+ 
  geom_point(aes(y=optionPriceMC, colour ="using Monte Carlo")) +
  scale_colour_manual("", 
                      breaks = c("using formula", "using Monte Carlo"),
                      values = c("using formula"="green", "using Monte Carlo"="red"))+
  ylab("option values")+
  ggsave("MCformula.pdf", width = 7, height = 4)


