library(tidyverse)
library(igraph)
library(stringr)
library(ggthemes)

draw_network <- function (N, M, num_edges) {
    r.graph <- sample_bipartite(N,M,type = "gnm", m=num_edges, directed = FALSE)
    as_incidence_matrix(r.graph)
}

calc_AM <- function (A, total_ext_assets) {
    total_ext_assets / (rowSums(A) + .Machine$double.eps) * A ## TODO: Fix me
}

bank_states <- function (solvent, total_ext_assets, total_liq_assets, deposits) {
    list(solvent = solvent,
         total_ext_assets = total_ext_assets,
         total_liq_assets = total_liq_assets,
         deposits = deposits,
         capital = total_ext_assets + total_liq_assets - deposits)
}

asset_states <- function (market_depth, prices) {
    list(market_depth = market_depth,
         prices = prices)
}

failing_banks <- function (AM, bank_states, asset_states) {
    ext_asset_value <- AM %*% asset_states$prices
    failing <- bank_states$solvent & (ext_asset_value + bank_states$total_liq_assets - bank_states$deposits < 0)
    print("failed banks: ")
    print(failing)
    ## Return failing banks
    failing
}

sell_assets <- function (AM, failing_banks, asset_states, params) {
    asset_amounts <- colSums(AM) + .Machine$double.eps
    frac_liquidated <- colSums(AM[failing_banks, , drop=FALSE])/asset_amounts
    print("frac_liquidated: ")
    print(frac_liquidated)
    
    ## Return new asset prices
    asset_states$prices * exp(- params$Alpha * frac_liquidated) ## TODO: include market depth
    print("new prices: ")
    print(asset_states$prices * exp(- params$Alpha * frac_liquidated))
    
}

propagate_asset_shock <- function (AM, bank_states, asset_states, params) {
    failed <- failing_banks(AM, bank_states, asset_states)
    if (any(failed)) {
        asset_prices <- sell_assets(AM, failed, asset_states, params)
        propagate_asset_shock(AM,
                              modifyList(bank_states, list(solvent=bank_states$solvent & !failed)),
                              modifyList(asset_states, list(prices=asset_prices)),
                              params)
    } else {
        list(banks = bank_states, assets = asset_states)
    }
}

propagate_bank_shock <- function (AM, bank_states, asset_states, params) {
    failed <- !bank_states$solvent
    asset_prices <- sell_assets(AM, failed, asset_states, params)
    propagate_asset_shock(AM,
                          modifyList(bank_states, list(solvent=bank_states$solvent & !failed)),
                          modifyList(asset_states, list(prices=asset_prices)),
                          params)
}

global_params <-
    list(## Bank parameters
        Liquid.asset.frac = 0.2,
        Ext.asset.frac = 0.8,
        Capital.frac = 0.1,
        ## Asset shock parameters
        Reduction.factor = 0.3,
        Alpha = 1.0536,
        ## Bailout parameters
        No.of.padded.banks = 5,
        New.capital.frac = 1)

init_states <- function (A, params) {
    N <- nrow(A)
    M <- ncol(A)
    Total.asset <- rep(1, N) ## c(rep(5, 20), rep(1, N-20))
    Total.ext.asset <- params$Ext.asset.frac * Total.asset 
    Total.liquid.asset <- params$Liquid.asset.frac * Total.asset
    Capital <- params$Capital.frac * Total.asset
    Deposits <- Total.ext.asset + Total.liquid.asset - Capital

    ## Fix banks with no connections to external assets
    Total.liquid.asset[rowSums(A) == 0] <- 1
    bs <- bank_states(rep(TRUE, N),
                      Total.ext.asset,
                      Total.liquid.asset,
                      Deposits)

    as <- asset_states(rep(1, M), rep(1, M))
    list(banks = bs, assets = as)
}

shock_assets <- function (asset_states, shocked_assets, params) {
    asset_states$prices[shocked_assets] <- params$Reduction.factor*asset_states$prices[shocked_assets]
    print("shocked asset_states$prices")
    print(asset_states$prices)
    asset_states
}

shock_banks <- function (bank_states, shocked_banks) {
    bank_states$solvent[shocked_banks] <- FALSE
    bank_states
}

run_asset_shocks <- function (graphml_files, scenario, params) {
    results <- data_frame(graph = character(0),
                          scenario = character(0),
                          shocked_asset = integer(0),
                          failed_banks = list())
    for (gf in graphml_files) {
        print(paste("Processing", gf))
        g <- read_graph(gf, "graphml")
        A <- as_incidence_matrix(g)

        states <- init_states(A, params)
        AM <- calc_AM(A, states$banks$total_ext_assets)
        
        if (scenario == "none") {
            ## Do nothing
        } else if (scenario == "rnd") {
            ## Randomly bailout banks
            random_padding <- sample(1:nrow(A),params$No.of.padded.banks,replace=F)

            total <- states$banks$total_ext_assets + states$banks$total_liq_assets 
            states$banks$deposits[random_padding] <-
                total[random_padding] - params$New.capital.frac*total[random_padding]
        } else if (scenario == "cont") {
            Omega <- AM %*% diag(1/states$assets$market_depth) %*% t(AM)
            lead.eig <- abs( eigen(Omega)$vectors[,1] )
            cont_padding <- order(lead.eig, decreasing=TRUE)[1:params$No.of.padded.banks]
            
            total <- states$banks$total_ext_assets + states$banks$total_liq_assets 
            states$banks$deposits[cont_padding] <-
                total[cont_padding] - params$New.capital.frac*total[cont_padding]
        } else if (scenario == "my") {
            Omega <- AM %*% diag(1/states$assets$market_depth) %*% t(AM)
            Omega.rel <- diag(params$Capital.frac, nrow=nrow(Omega)) %*% Omega
            cent <- rowSums(Omega.rel) * colSums(Omega.rel)
            my_padding <- order(cent, decreasing=TRUE)[1:params$No.of.padded.banks]

            total <- states$banks$total_ext_assets + states$banks$total_liq_assets 
            states$banks$deposits[my_padding] <-
                total[my_padding] - params$New.capital.frac*total[my_padding]            
        } else if (scenario == "sasi") {
            ## Compute bank centralities
            proj <- bipartite_projection(g)
            banks.one.mode <- proj[[1]]
            assets.one.mode <- proj[[2]]
            N <- nrow(A)
            M <- ncol(A)
            Bank.centrality <- rep(0, N)
            for(i in 1:(M-1)) {
                for(j in (i+1):M) {
                    if(sum(AM[,i]*AM[,j])==0) next
                    Bank.centrality <- Bank.centrality + AM[,i]*AM[,j]/sum(AM[,i]*AM[,j])*
                        edge_betweenness(assets.one.mode,e=E(assets.one.mode)[i%--%j],
                                         directed=FALSE,weights=NULL)
                }
            }
            ## Bank.centrality <- Bank.centrality/(.Machine[[1]] + (rowSums(A) * (rowSums(A) - 1))/2)

            padding <- order(Bank.centrality, decreasing=TRUE)[1:params$No.of.padded.banks]
            
            total <- states$banks$total_ext_assets + states$banks$total_liq_assets 
            states$banks$deposits[padding] <-
                total[padding] - params$New.capital.frac*total[padding]
        } else {
            stop("Unknown scenario ", scenario)
        }

        for (idx in 1:ncol(A)) {
            ## Shock asset i
            as <- shock_assets(states$assets, idx, params)
            
            sim <- propagate_asset_shock(AM, states$banks, as, params)
            ## Store simulation results
            results <- bind_rows(results,
                                 list(graph=gf, scenario=scenario, shocked_asset=idx,
                                      failed_banks=list(which(!sim$banks$solvent))))
        }
    }
    results
}

#--------------------------------------#
# Run simple example step by step using functions of this script

# A = matrix(c(0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1), nrow=4, ncol=7)
# A = matrix(c(0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0), nrow=3, ncol=5)
A = matrix(c(0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0), nrow=3, ncol=5)

params = global_params
states <- init_states(A, params)
# print(states$banks$total_ext_assets)
# print(states$banks$total_liq_assets)
# print(states$banks$deposits)

AM <- calc_AM(A, states$banks$total_ext_assets)

# results data frame
results <- data_frame(shocked_asset = integer(0),
                      failed_banks = list())

for (idx in 1:ncol(A)) {
  print(paste0("IDX: ", idx))
  ## Shock asset i
  as <- shock_assets(states$assets, idx, params)

  sim <- propagate_asset_shock(AM, states$banks, as, params)

  ## Store simulation results
  results <- bind_rows(results,
                     list(shocked_asset=idx,
                          failed_banks=list(which(!sim$banks$solvent))))
}
