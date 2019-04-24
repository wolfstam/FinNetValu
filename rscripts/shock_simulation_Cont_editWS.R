################################################################################
## Generic framework to simulate financial contagion due to
## overlapping portfolios following Greenwood
##
## The model is formulated in terms of matrix operations and should
## include leverage targeting (Greenwood), maximum leverage targeting
## (Cont) as well as default contagion (Caccioli)
################################################################################

## Model follows Cont (eq. (21) and (22))
## Changes:
## * Separate bank balance sheet size and portfolio holdings (as in Greenwood)
## * Flexible parametrization to allow for non-linear market impact and leverage targeting rules

library(tidyverse)
library(readxl)
library(lubridate)

state <- function (capital, exposures) {
    ## Exposures corresponds to Pi matrix of Cont

    num_banks = nrow(exposures)
    num_assets = ncol(exposures)

    bank_sizes = rowSums(exposures)
    portfolios = exposures / bank_sizes
    
    list(capital = capital,
         Pi = exposures)
         ## bank_balance_sheet = bank_sizes,
         ## portfolios = portfolios)
}

params <- function (market_impact, leverage_strategy, alpha) {
    list(market_impact = market_impact,
         leverage_strategy = leverage_strategy,
         alpha = alpha)
}

sim_step <- function (state, loss, params) {
    Pi <- state$Pi ## diag(state$bank_balance_sheet) %*% state$portfolios

    Gamma <- params$leverage_strategy(state, loss)
    ## Compute new exposure values (13)
    Psi <- params$market_impact(c(Gamma %*% Pi), colSums(Pi))
    Pi_next <- diag(1 - Gamma) %*% Pi %*% diag(1 - Psi)

    ## Compute losses
    ## M <- (1 - Gamma) * rowSums(Pi - Pi_next)
    ## R <- Gamma * rowSums(Pi - (params$alpha*Pi_next + (1-params$alpha)*Pi))
    M <- (1 - Gamma) * (Pi %*% Psi)
    R <- (1 - (1 - params$alpha) * Gamma) * (Pi %*% Psi)
    L <- c(M + R)

    ## Update state
    state_new = modifyList(state,
                           list(Pi = Pi_next,
                                capital = pmax(state$capital - L, 0)))

    list(state = state_new, loss = L)
}

simulate <- function (state, loss, params, iter=100) {
    banks <- rownames(state$Pi)
    if (is.null(banks))
        banks <- 1:nrow(state$Pi)
    
    df <- data_frame(bank = banks, loss = loss, step = 0,
                     capital = state$capital, assets = rowSums(state$Pi))
    sl <- list(state = state, loss = loss)
    for (i in 1:iter) {
        sl_next <- sim_step(sl$state, sl$loss, params)
        if (all(sl_next$loss == 0))
            break
        df <- bind_rows(df, data_frame(bank = banks, loss = sl_next$loss, step = i,
                                       capital = sl_next$state$capital, assets = rowSums(sl_next$state$Pi)))
        sl <- sl_next
    }
    df
}   

market_impact_linear <- function (D) {
    function (q, S) {
        min(1, q / D)
    }
}

market_impact_exp <- function (D) {
    function (q, S) {
        1 - exp(- q / D)
    }
}        

market_impact_frac <- function (q, S) {
    1 - exp(- 1.0536 * q / (S + .Machine$double.eps))
}

market_impact_cont <- function (D, B, S0) {
    delta <- (1 - B / S0) * D
    function (q, S) {
        (1 - B / S) * (1 - exp(- q / delta))
    }
}
        
leverage_cont <- function (lambda_max, lambda_b) {
    function (state, loss) {
        ass <- rowSums(state$Pi)
        cap <- state$capital
        lambda <- ass / (cap + .Machine$double.eps)
        Gamma <- ((lambda_b - 1) * loss - cap * (lambda_b - lambda)) / (ass + .Machine$double.eps)
        pmin(Gamma * (lambda > lambda_max), 1)
    }
}

leverage_target <- function (leverage_target) {
    function (state, loss) {
    }
}

leverage_default <- function (state, loss) {
    as.numeric(!(state$capital > 0))
}

shock_assets <- function (state, assets, reduction_factor) {
    Pi_shock <- state$Pi
    Pi_shock[,assets] <- reduction_factor * Pi_shock[,assets]
    loss <- rowSums(state$Pi) - rowSums(Pi_shock)
    cap <- pmax(state$capital - loss, 0)
    list(state = modifyList(state, list(capital = cap, Pi = Pi_shock)),
         loss = loss)
}

# #source("EBA_Cont.R")
# source("EBA_Cont_editWS.R")

#--------------------------------------------------------------#
# tests by WS 
# two bank example from Shaanning 2017 (thesis)
Pi = cbind(c(90, 70))
cap = c(4, 4.5)
rownames(Pi) <- c("Bank A","Bank B")
#colnames(cap) <- c("Bank A","Bank B")
colnames(Pi) <- c("Asset 1")
D = 0.4 * (50/0.02) * sqrt(20)
B = 0.5
S0 = 1
alpha = .5
lambda_max = 33.3
lambda_b = 0.95 * lambda_max
#--------------------------------------------------------------#

banks <- rownames(Pi)
N <- nrow(Pi)
M <- ncol(Pi)
state_cont <- state(cap, Pi)

params_default <- params(market_impact_cont(D, B, S0), leverage_cont(lambda_max, lambda_b), alpha)

#--------------------------------------------------------------#

# seems to correspond to epsilon = [0.2, 0.] and Theta from two bank example
sl <- shock_assets(state_cont, colnames(Pi), c(0.9778, 1.))
Pi <- sl$state$Pi ## diag(state$bank_balance_sheet) %*% state$portfolios

Gamma <- params_default$leverage_strategy(sl$state, sl$loss)
## Compute new exposure values (13)
Psi <- params_default$market_impact(c(Gamma %*% Pi), colSums(Pi))
Pi_next <- diag(1 - Gamma) %*% Pi %*% (1 - Psi)

## Compute losses
## M <- (1 - Gamma) * rowSums(Pi - Pi_next)
## R <- Gamma * rowSums(Pi - (params$alpha*Pi_next + (1-params$alpha)*Pi))
M <- (1 - Gamma) * (Pi %*% Psi)
R <- (1 - (1 - params_default$alpha) * Gamma) * (Pi %*% Psi)
L <- c(M + R)
pmax(sl$state$capital - L, 0)

sl_next <- sim_step(sl$state, sl$loss, params_default)

# -> it seems theta is not used, nor is S updated, I don't know how to test the two bank example in this way 
# without being able to specify theta and epsilon -> all computed values in between, e.g. Gamma, Psi, Pi_next, L
# are weird numbers, partly with different value ranges, which I do not know how to interpret
#--------------------------------------------------------------#

# shock_banks <- function (state, banks) {
#     N <- nrow(state$Pi)
# 
#     loss <- rep(0, N)
#     names(loss) <- rownames(state$Pi)
#     loss[banks] <- state$capital[banks]
#     cap <- state$capital
#     cap[banks] <- 0
#     list(state = modifyList(state, list(capital = cap)),
#          loss = loss)
# }
# 
# 
# sim_df <- data_frame(shocked_bank = banks,
#                      loss = 
#                          lapply(banks, function (b) {
#                              sl <- shock_banks(state_cont, b)
#                              simulate(sl$state, sl$loss, params_default) %>%
#                                  ## Compute total loss
#                                  group_by(bank) %>%
#                                  arrange(desc(step)) %>%
#                                  summarize(loss = sum(loss), capital = first(capital))
#                          })) %>%
#     unnest()
# 
# sim_df <- data_frame(shocked_asset = colnames(Pi),
#                      loss = 
#                          lapply(shocked_asset, function (a) {
#                              sl <- shock_assets(state_cont, a, 0.95)
#                              simulate(sl$state, sl$loss, params_default)
#                          })) %>%
#     unnest()
# 
# ## Investigate number of failed banks
# sim_df %>%
#     ## Compute total loss
#     group_by(bank, shocked_asset) %>%
#     arrange(step) %>%
#     summarize(loss = sum(loss), capital = last(capital), assets = last(assets)) %>%
#     ungroup() %>%
#     group_by(shocked_asset) %>%
#     mutate(solvent = capital > 0) %>%
#     summarize(num_failed = sum(1 - solvent)) %>%
#     ggplot(aes(shocked_asset, num_failed)) +
#     geom_point()
# 
# sim_df %>%
#     ## Compute total loss
#     group_by(bank, shocked_asset) %>%
#     arrange(step) %>%
#     summarize(loss = sum(loss), capital = last(capital), assets = last(assets)) %>%
#     ungroup() %>%
#     ## group_by(bank) %>%
#     ## summarize(loss = mean(loss)) %>%
#     full_join(data_frame(bank = names(cap), equity = cap)) %>%
#     mutate(insolvent = !(capital > 0), illiquid = !(assets > 0)) %>%
#     ggplot(aes(bank, loss/equity, color=interaction(insolvent, illiquid))) +
#     geom_point() +
#     ## scale_y_log10() +
#     coord_flip() 
#     ## geom_jitter()
# 
# ## Centrality vs total loss
# Omega <- Pi %*% diag(1e-13, nrow=ncol(Pi)) %*% t(Pi)
# centrality <- eigen(Omega)$vectors[,1]
# sim_df %>%
#     full_join(data_frame(bank = rownames(Omega), centrality = centrality)) %>%
#     group_by(bank, shocked_asset, step > 0) %>%
#     arrange(step) %>%
#     ## filter(step > 0) %>% ## Only fire sale losses
#     summarize(loss = sum(loss), capital = last(capital), assets = last(assets), centrality = unique(centrality)) %>%
#     ungroup() %>%
#     ggplot(aes(centrality, loss, color=`step > 0`)) +
#     geom_point()
