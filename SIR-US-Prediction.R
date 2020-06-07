## Title: SIR-US-Prediction
## Author: MFR
## Date: June 2, 2020
## File Name: SIR-US-Prediction.R
## Loading required libraries including deSolve package for SIR Model
# install.packages("deSolve")
library(deSolve)
library(tidyverse)
library(data.table)
library(magrittr)
library(lubridate)
library(gridExtra)
library(kableExtra)
library(plotly)
library(DT)
library(ggplot2)
library(dplyr)
####################################################################################
## Downloading dataset from Kaggle for initial analysis
## https://www.kaggle.com/sudalairajkumar/novel-corona-virus-2019-dataset/discussion/139288
## Loading Hmisc library translating SAS datasets into S
library(Hmisc)
## Reading COVID19_line_list_data.csv from local drive
corona <- read.csv("C:/Users/MfR/COVID19-Capstone/Kaggle/COVID19_line_list_data.csv")
glimpse(corona)

## The number of rows and columns in the dataset
dim(corona)

# if death variable isn't 0, then the patient has died
corona$death_dummy <- as.integer(corona$death != 0)

# death rate of our dataset
sum(corona$death_dummy) / nrow(corona)

# Age Analysis
dead = subset(corona, death_dummy == 1)
alive = subset(corona, death_dummy == 0)
mean(dead$age, na.rm=TRUE)
mean(alive$age, na.rm=TRUE)

t.test(dead$age, alive$age, alternative="two.sided", conf.level = 0.95)

# Gender Analysis
men = subset(corona, gender == "male")
women = subset(corona, gender == "female")
mean(men$death_dummy, na.rm=TRUE)
mean(women$death_dummy, na.rm=TRUE)
t.test(men$death_dummy, women$death_dummy, alternative="two.sided", conf.level = 0.95)

v3 <- rbind(men,women) # combine objects as rows

summary(v3)
dim(v3)
## Loading caret package for Machine Learning
library(caret)

y <- v3$gender
summary(y)
x <- v3$age
summary(x)

## Splitting data into Training and test sets
set.seed(2007)
test_index <- createDataPartition(y, times = 1, p = 0.2, list = FALSE)

test_set <- v3[test_index, ]
train_set <- v3[-test_index, ]

y_hat <- sample(c("male", "female"), length(test_index), replace = TRUE)

y_hat <- sample(c("male", "female"), length(test_index), replace = TRUE) %>%
  factor(levels = levels(test_set$gender))

# The overall accuracy is simply defined as the overall proportion that is predicted correctly:

mean(y_hat == test_set$gender)

table(predicted = y_hat, actual = test_set$gender)

# If we study this table closely, it reveals a problem. If we compute the accuracy separately for each sex, we get:

test_set %>%
  mutate(y_hat = y_hat) %>% group_by(gender) %>%
  summarise(accuracy = mean(y_hat == gender))

# Not surprisingly, our accuracy is about 55%. We are guessing!
prev <- mean(y == "female")
prev

cm <- confusionMatrix(data = y_hat, reference = test_set$gender)
cm$overall["Accuracy"]

cm$byClass[c("Sensitivity","Specificity", "Prevalence")]

####################################################################################
## Loading covid19 Library that updates COVID-19 Data from JHU-CCSE:
## install.packages("devtools")
devtools::install_github("irays/covid19")

library(covid19)
tail(covid19)
## The number of rows and columns in the dataset with data types in different columns
dim(covid19)
str(covid19)
####################################################################################
## Countrywise Tabular Data
covid19 %>%
  filter(Country.Region != "Others") %>%
  group_by(Country.Region, type) %>%
  summarise(total_cases = sum(cases)) %>%
  pivot_wider(names_from = type, values_from = total_cases) %>%
  arrange(- confirmed) %>%
  filter(confirmed >= 25) %>%
  mutate(death_rate = death / confirmed)  %>%
  datatable(rownames = FALSE,
            colnames = c("Confirmed", "Death","Recovered", "Death Rate")) %>%
  formatPercentage("death_rate", 2)
####################################################################################
# Coronavirus - Worldwide Cumulative Distribution as Bar Diagram

`%>%` <- magrittr::`%>%`

# extract the cumulative incidence
df <- covid19 %>%
  dplyr::filter(Country.Region == "US") %>%
  dplyr::group_by(date, type) %>%
  dplyr::summarise(total = sum(cases, na.rm = TRUE)) %>%
  tidyr::pivot_wider(
    names_from = type,
    values_from = total
  ) %>%
  dplyr::arrange(date) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(active = confirmed - death - recovered) %>%
  dplyr::mutate(
    confirmed_cum = cumsum(confirmed),
    death_cum = cumsum(death),
    recovered_cum = cumsum(recovered),
    active_cum = cumsum(active)
  )

total_cases <- covid19 %>%
  group_by(type) %>%
  summarise(cases = sum(cases)) %>%
  mutate(type = factor(type, levels = c("confirmed", "recovered","death")))

total_cases

plot_ly(data = total_cases,
        x = ~ type,
        y = ~cases,
        type = 'bar',
        text = ~ paste(type, cases, sep = ": "),
        hoverinfo = 'text') %>%
  layout(title = "Coronavirus - Cases Distribution",
         yaxis = list(title = "Number of Cases"),
         xaxis = list(title = "Case Type"),
         hovermode = "compare")

###################################################################################
# Table of the ten countries with the highest confirmed cases.
# We will use the `datatable` function from the **DT** package to view the table:

confirmed_country <- covid19 %>%
  filter(type == "confirmed") %>%
  group_by(Country.Region) %>%
  summarise(total_cases = sum(cases)) %>%
  mutate(perc = total_cases / sum(total_cases)) %>%
  arrange(-total_cases)

confirmed_country %>%
  head(10) %>%
  datatable(rownames = FALSE,
            colnames = c("Country", "Cases", "Perc of Total")) %>%
  formatPercentage("perc", 2)

##################################################################################
## The pie plot summarize the distribution of the worldwide confirmed cases:

covid19 %>%
  filter(type == "confirmed") %>%
  group_by(Country.Region) %>%
  summarise(total_cases = sum(cases)) %>%
  arrange(-total_cases) %>%
  mutate(country = factor(Country.Region, levels = Country.Region)) %>%
  ungroup() %>%
  plot_ly(labels = ~ country,
          values = ~ total_cases,
          type = "pie",
          textposition = 'inside',
          textinfo = 'label+percent',
          insidetextfont = list(color = '#FFFFFF'),
          hoverinfo = 'text',
          text = ~ paste(country, "<br />",
                         "Number of confirmed cases: ", total_cases, sep = "")) %>%
  layout(title = "Coronavirus - Confirmed Cases")

####################################################################################

# SIR Model
SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta * I * S / N
    dI <- beta * I * S / N - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

####################################################################
# put the daily cumulative incidence numbers for US from
# February 29 to May 04 into a vector called Infected

sir_start_date <- "2020-02-29"
sir_end_date <- "2020-06-04"
# active = confirmed-death-recovered
# active_cum = cumsum(active)
Infected <- subset(df, date >= ymd(sir_start_date) & date <= ymd(sir_end_date))$active_cum

# Create an incrementing Day vector the same length as our
# cases vector
Day <- 1:(length(Infected))

# now specify initial values for N, S, I and R
# N is the population of New York City (8.399M), Washington (7.615) and Oregon (4.218)
N <- 20232000
# N <- 331002651 # (total US population)
# The initial values of the variables need to be defined in a named vector:
init <- c(
  S = N - Infected[1],
  I = Infected[1],
  R = 0
)


## Define a function to calculate the quantity is known as the
# residual mean squared error (RMSE):
## Passing in parameters beta and gamma that are to be optimised for the best fit to the
## incidence data.

RMSE <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = Day, func = SIR, parms = parameters)
  fit <- out[, 3]
  # sum((Infected - fit)^2)
  sqrt((1/N)*sum((Infected - fit)^2))
}

# Now find the values of beta and gamma that give the smallest RMSE, which represents
# the best fit to the data. Start with values of 0.5 for each, and constrain them to
# the interval 0 to 1.0

Opt <- optim(c(0.50, 0.50),
             RMSE,
             method = "L-BFGS-B",
             lower = c(0, 0),
             upper = c(1, 1)
)

# check for convergence
Opt$message
# Convergence is confirmed. Now we can examine the fitted values for beta and gamma

Opt_par <- setNames(Opt$par, c("beta", "gamma"))
Opt_par

# time in days for predictions
# t <- 1:as.integer(ymd(sir_end_date) - 1 - ymd(sir_start_date)) # odd date
t <- 1:as.integer(ymd(sir_end_date) + 1 - ymd(sir_start_date)) # odd date
# t <- 1:as.integer(ymd(sir_end_date) - ymd(sir_start_date)) # even date

# get the fitted values from our SIR model
fitted_cumulative_incidence <- data.frame(ode(
  y = init, times = t,
  func = SIR, parms = Opt_par
))

# add a Date column and the observed incidence data

fitted_cumulative_incidence <- fitted_cumulative_incidence %>%
  mutate(
    Date = ymd(sir_start_date) + days(t),
    Country = "US",
    cumulative_incident_cases = Infected
  )

# plot the data

fitted_cumulative_incidence %>%
  ggplot(aes(x = Date)) +
  geom_line(aes(y = I), colour = "red") +
  geom_point(aes(y = cumulative_incident_cases), colour = "blue") +
  labs(
    y = "Cumulative Infected",
    title = "COVID-19 SIR Fitted (Red) vs Infected (Blue), US",
    subtitle = "(Red = fitted from SIR model, blue = JHUCCSE)"
  ) +
  theme_minimal()
# The following graph is similar than the previous one, except that the y-axis is
# measured on a log scale.
fitted_cumulative_incidence %>%
  ggplot(aes(x = Date)) +
  geom_line(aes(y = I), colour = "red") +
  geom_point(aes(y = cumulative_incident_cases), colour = "blue") +
  labs(
    y = "Cumulative Infected",
    title = "COVID-19 SIR Fitted (Red) vs Infected (Blue), US",
    subtitle = "(Red = SIR Fitted Curve, blue = JHUCCSE)"
  ) +
  theme_minimal() +
  scale_y_log10(labels = scales::comma)

# Reproduction number

Opt_par

R0 <- as.numeric(Opt_par[1] / Opt_par[2])
R0

# A $R_{0}$ of 1.3 means that, on average in US, 1.3 persons are infected for
# each infected person.

# time in days for predictions
t <- 1:150

# get the fitted values from our SIR model
fitted_cumulative_incidence <- data.frame(ode(
  y = init, times = t,
  func = SIR, parms = Opt_par
))

# add a Date column and join the observed incidence data
fitted_cumulative_incidence <- fitted_cumulative_incidence %>%
  mutate(
    Date = ymd(sir_start_date) + days(t),
    Country = "US",
    cumulative_incident_cases = c(Infected, rep(NA, length(t) - length(Infected)))
  )

# plot the data
fitted_cumulative_incidence %>%
  ggplot(aes(x = Date)) +
  geom_line(aes(y = I), colour = "red") +
  geom_line(aes(y = S), colour = "black") +
  geom_line(aes(y = R), colour = "green") +
  geom_point(aes(y = cumulative_incident_cases),
             colour = "blue"
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(y = "Persons", title = "COVID-19 SIR Fitted (Red) vs Infected (Blue), US") +
  scale_colour_manual(name = "", values = c(
    red = "red", black = "black",
    green = "green", blue = "blue"
  ), labels = c(
    "Susceptible",
    "Recovered", "Observed", "Infectious"
  )) +
  theme_minimal()

# plot the data
fitted_cumulative_incidence %>%
  ggplot(aes(x = Date)) +
  geom_line(aes(y = I, colour = "red")) +
  geom_line(aes(y = S, colour = "black")) +
  geom_line(aes(y = R, colour = "green")) +
  geom_point(aes(y = cumulative_incident_cases, colour = "blue")) +
  scale_y_log10(labels = scales::comma) +
  labs(
    y = "Persons",
    title = "COVID-19 SIR Fitted (Red) vs Infected (Blue), US"
  ) +
  scale_colour_manual(
    name = "",
    values = c(red = "red", black = "black", green = "green", blue = "blue"),
    labels = c("Susceptible", "Observed", "Recovered", "Infectious")
  ) +
  theme_minimal()

# peak of pandemic with date
fit <- fitted_cumulative_incidence

fit[fit$I == max(fit$I), c("Date", "I")]

# RMSE with optimized parameter
RMSE(Opt_par)
####################################################################################
