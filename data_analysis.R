##############################################################################
# Evtreme Value Analysis of S&P 500  Data
# Author: Andrew Benson
# Date: 12/14/21
# Description: This file imports a .csv file with historical data
# on returns from the S&P 500. Some graphics are created on the general data
# as well as on the block maximum data. The data are fit to an extreme value
# distribution using the evd package.
##############################################################################

library("evd")
library("tidyverse")
library("lubridate")
library("gridExtra")

#----------------------------------------------------------------------------
# FOR CODE TEST ONLY
#---------------------------------------------------------------------------
#View(portpirie)
#years <- seq(1923, 1987, by = 1)
#plot(years, portpirie, main="Annual max sea levels at Port Pirie, South Australia",
#     xlab="Year", ylab="Series-level (meters)", pch=19)
#pmle = fgev(portpirie)
#print(pmle$estimate)
#confint(mle, level = 0.95)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Data tidying, data analysis, and simple data vizulization
#----------------------------------------------------------------------------
# Import the data
spx <- readr::read_csv("/home/drew/Documents/GradSchool/Stochastic_Processes_I/CourseProject/Data/spx_data.csv")
prices <- spx$`Adj Close`
# Returns will be calculated as ln(price_{t+1} / price_t)
returns <- diff(log(prices), lag=1)
returns <- append(returns, 0, 0)
spx_dates_converted <- lubridate::mdy(spx$Date)
spx["Returns"] <- returns
spx["Year"] <- lubridate::year(spx_dates_converted)

# Filter the data to 1960 - 2020
spx_filtered <- spx %>%
  filter(Year >= 1960, Year <= 2020)

# Time series plot of index level over time
plot(spx$`Adj Close`, type="l", main="S&P 500 Adjusted Closing Price",
     xlab="Year", ylab="Index Level", pch=20, xaxt='n')
axis(1, at=c(1, 2490, 5016, 7544, 10072, 12839, 15103), labels=c("1960", "1970", "1980", "1990", "2000", "2010", "2020"))

# Plot of S&P 500 log returns
plot(spx$Returns, type="l", main="S&P 500 Log Daily Returns",
     xlab="Year", ylab="Return", pch=20, xaxt='n')
axis(1, at=c(1, 2490, 5016, 7544, 10072, 12839, 15103), labels=c("1960", "1970", "1980", "1990", "2000", "2010", "2020"))

# Filter down the data to the years under consideration; 1960-2020
block_max <- spx %>%
  filter(Year >= 1960, Year <= 2020) %>% 
  group_by(Year) %>% 
  summarise(min = min(Returns))

# Plot the minimum value for each year; 1960 - 2020
plot(block_max$Year, block_max$min * 100, main="Annual Minima for S&P 500",
     xlab="Year", ylab="Return in Percent", pch=20)



#----------------------------------------------------------------------------
# This section is more focused on statistical analysis
#----------------------------------------------------------------------------
# Maximum likelihood estimation of the GEV
gev_mle <- fgev(abs(block_max$min))
# MLE parameter estimates
print(gev_mle$estimate)
#loc      scale      shape 
#0.02369486 0.01063665 0.51805068 
gev_loc = gev_mle$estimate[1]
gev_scale = gev_mle$estimate[2]
gev_shape = gev_mle$estimate[3]

# Confidence intervals for MLE estimates
print(confint(gev_mle, level = 0.95))
#  2.5 %     97.5 %
#loc   0.020559498 0.02683023
#scale 0.007878369 0.01339494
#shape 0.220189597 0.81591177

# MLE was run again assuming Frechet distribution
f_mle <- fextreme(abs(block_max$min), list(loc = 0.02369486, scale = 0.01063665, shape = gev_shape), distn = "frechet", mlen = 61, largest = TRUE, std.err = FALSE)
f_loc = f_mle$estimate[1]
f_scale = f_mle$estimate[2] # Turns our the results are identical to the MLE run on GVD

# Histogram of the minimum values for each year; 1960 - 2020
emp_hist <- block_max %>%
  ggplot( aes(x=abs(min))) +
  geom_histogram( binwidth=0.02, fill="#22282e", color="#e9ecef", alpha=0.9) +
  theme_bw() + 
  labs(x = "Return", y = "Count", title = "Histogram of Minimum Returns")

# The frechet distribution with the MLE parameters
X  <-  seq(0.011, 0.22, length = 1000)
y  <- dfrechet(X, loc = gev_loc, scale = gev_scale, gev_shape = gev_mle$estimate[3])
df <- data.frame(x = X,y)
f_dist <- ggplot(df, aes(x = X, y)) +
  geom_line() +
  theme_bw() + 
  labs(x = "Return", y = "Density", title = "Distribution of Minimum Returns")
f_dist
# Plot the histogram and the Frechet distribution side by side
grid.arrange(emp_hist, f_dist, nrow=1, ncol=2)

#----------------------------------------------------------------------------
# Checking the model fit
#----------------------------------------------------------------------------
# Probability plot of the empierical distribution function against
# a Frechet distribution function with MLE estimates
par(mfrow = c(1,2))
u <- sort(abs(block_max$min))
ecdf <- seq_along(u) / length(u)
frechet_estimated <- pfrechet(u, loc = gev_loc, scale = gev_scale, shape = 1)
plot(x = ecdf, y = frechet_estimated, xlim = c(0,1), ylim = c(0,1), 
     xlab = "ECDF",  ylab = "Frechet Fitted", main = "Probability Plot")
abline(0, 1, col = 'red')

u <- sort(abs(block_max$min))[25:61]
ecdf <- seq_along(u) / length(u)
frechet_estimated <- pfrechet(u, loc = gev_loc, scale = gev_scale, shape = 1)
plot(x = ecdf, y = frechet_estimated, xlim = c(0,1), ylim = c(0,1),  
     xlab = "ECDF",  ylab = "Frechet Fitted", main = "Probability Plot")
abline(0, 1, col = 'red')


# Check to see if the Frechet distribution is a good fit
# Generate a 61 random draws from a frechet distribution with the MLE parameters
frechet_sample <- rfrechet(61, loc = gev_loc, shape = gev_scale, scale = gev_shape)
# Run the Kolmogrov Smirnov test  
res <- stats::ks.test(abs(block_max$min), frechet_sample)
pander::pander(res)


m <- mean(spx$Returns)
s <- sd(spx$Returns)
u <- sort(abs(block_max$min))
frechet_estimated <- pfrechet(u, loc = gev_loc, scale = gev_scale, shape = 1)

norm_estimates <- (1 - stats::pnorm(u, mean = m, sd = s))
difference <- abs(norm_estimates - frechet_estimated)
plot(frechet_estimated, norm_estimates)
