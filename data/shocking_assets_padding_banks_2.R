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
    ## Return failing banks
    failing
}

sell_assets <- function (AM, failing_banks, asset_states, params) {
    asset_amounts <- colSums(AM) + .Machine$double.eps
    frac_liquidated <- colSums(AM[failing_banks, , drop=FALSE])/asset_amounts
    ## Return new asset prices
    asset_states$prices * exp(- params$Alpha * frac_liquidated) ## TODO: include market depth
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

test_this <- function (A) {
    Total.asset <- matrix(data = 1, nrow = N, ncol = 1)
    Total.ext.asset <- Ext.asset.frac * Total.asset 
    Total.liquid.asset <- Liquid.asset.frac * Total.asset
    Capital <- Capital.frac * Total.asset
    Deposits <- Total.ext.asset + Total.liquid.asset - Capital

    AM <- calc_AM(A, c(Total.ext.asset))
    Total.liquid.asset[rowSums(A) == 0] <- 1
    bs <- bank_states(rep(TRUE, N),
                      Total.ext.asset,
                      Total.liquid.asset,
                      Deposits)

    tada <- matrix(FALSE, nrow=M, ncol=N)
    for (i in 1:M) {
        ## Shock asset i
        as <- asset_states(NULL, rep(1, M))
        as$prices[i] <- Reduction.factor*as$prices[i]
        
        foo <- propagate_asset_shock(AM, bs, as, list(alpha=Alpha))
        tada[i,] <- !foo$banks$solvent
    }
    tada
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

create_and_save_graphs <- function (N, M, bank_degrees) {
    sapply(1:length(bank_degrees),
           function (id) {
               bd <- bank_degrees[id]
               g_name <- paste0("Banks", N, "Assets", M, "bank_degree", bd, "ID", id, ".graphml")
               E <- floor(N*bd)
               g <- sample_bipartite(N,M,type = "gnm", m=E, directed = FALSE)
               write_graph(g, g_name, "graphml")
           })
    "Done"
}

## create_and_save_graphs(25, 30, rep(c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,8,10,12,14), each=50))

## rnd_bailout <- run_asset_shocks(dir(".", pattern=".graphml"), "rnd", global_params)

str_extract_num_after <- function (str, after_pat) {
    num_regex = "\\d+(.\\d+)?"
    str %>%
        str_extract(paste0(after_pat, num_regex)) %>%
        str_extract(num_regex) %>%
        as.numeric()
}

## tada %>%
##     mutate(deg = str_extract_num_after(graph, "bank_degree")) %>%
##     mutate(num_failed = unlist(Map(length, failed_banks))) %>%
##     group_by(graph, scenario, deg) %>%
##     summarize(failed = mean(num_failed)) %>%
##     ggplot(aes(x=deg, y=failed, fill=scenario)) +
##     geom_boxplot(aes(group=interaction(deg,scenario))) +
##     scale_color_colorblind() +
##     geom_smooth(aes(x=deg, color=scenario), se=FALSE)

## Load Cont data
source("../EBA_Cont.R")

cont_states <- function (Pi, cap) {
    N <- nrow(Pi)
    banks <- rownames(Pi)
    bank_assets <- rowSums(Pi)
    names(bank_assets) <- banks
    bs <- bank_states(rep(TRUE, N),
                      bank_assets[banks],
                      rep(0, N),
                      pmax(bank_assets[banks] - cap[banks], 0))
    M <- ncol(Pi)
    as <- asset_states(rep(1e13, M), rep(1, M))
    list(banks = bs, assets = as)
}

random_padding <- function (N) {
    sample(1:N,N,replace=F)
}

cont_padding <- function (AM, market_depth) {
    Omega <- AM %*% diag(1/market_depth) %*% t(AM)
    lead.eig <- abs( eigen(Omega)$vectors[,1] )
    order(lead.eig, decreasing=TRUE)
}

my_padding <- function (AM, market_depth, cap) {
    Omega <- AM %*% diag(1/market_depth) %*% t(AM)
    banks <- rownames(Omega)
    Omega.rel <- diag(1/cap[banks]) %*% Omega
    cent <- rowSums(Omega.rel) * colSums(Omega.rel)
    order(cent, decreasing=TRUE)
}

sasi_padding <- function (AM) {
    ## Compute bank centralities
    g <- graph_from_incidence_matrix(AM)
    proj <- bipartite_projection(g)
    banks.one.mode <- proj[[1]]
    assets.one.mode <- proj[[2]]
    N <- nrow(AM)
    M <- ncol(AM)
    Bank.centrality <- rep(0, N)
    for(i in 1:(M-1)) {
        print(i)
        for(j in (i+1):M) {
            if(sum(AM[,i]*AM[,j])==0) next
            Bank.centrality <- Bank.centrality + AM[,i]*AM[,j]/sum(AM[,i]*AM[,j])*
                edge_betweenness(assets.one.mode,e=E(assets.one.mode)[i%--%j],
                                 directed=FALSE,weights=NULL)
        }
    }
    ## Bank.centrality <- Bank.centrality/(.Machine[[1]] + (rowSums(A) * (rowSums(A) - 1))/2)
    
    order(Bank.centrality, decreasing=TRUE)
}

paddings <- list(rand = random_padding(nrow(Pi)),
                 cont = cont_padding(Pi, rep(1, ncol(Pi))),
                 my   = my_padding(Pi, rep(1e-13, ncol(Pi)), cap),
                 sasi = sasi_padding(Pi))

run_asset_shocks_cont <- function (params, scenario, padding_order="none") {
    results <- data_frame(scenario = character(0),
                          shocked_asset = integer(0),
                          failed_banks = list(),
                          params = list(),
                          num_padded = integer(0),
                          pad_total_cost = numeric(0))

    states <- cont_states(Pi, cap)
    AM <- Pi

    if (padding_order == "none") {
        ## No padding ... do nothing
        pad_total_cost <- 0
    } else {
        padded_banks<- padding_order[1:params$No.of.padded.banks]
        total <- states$banks$total_ext_assets + states$banks$total_liq_assets 
        states$banks$deposits[padded_banks] <-
            total[padded_banks] - params$New.capital.frac*total[padded_banks]
        pad_total_cost <- sum((params$New.capital.frac - params$Capital.frac)*total[padded_banks])
    }
    
    for (idx in 1:ncol(AM)) {
        ## Shock asset i
        as <- shock_assets(states$assets, idx, params)
        
        sim <- propagate_asset_shock(AM, states$banks, as, params)
        ## Store simulation results
        results <- bind_rows(results,
                             list(scenario=scenario, shocked_asset=idx,
                                  failed_banks=list(which(!sim$banks$solvent)),
                                  pad_total_cost = pad_total_cost,
                                  num_padded = params$No.of.padded.banks,
                                  params = list(params)))
    }
    results
}

global_params <- modifyList(global_params, list(Reduction.factor=0.3))
df <- run_asset_shocks_cont(modifyList(global_params, list(No.of.padded.banks = 0)), "baseline", "none")
for (np in seq(5,50,by=5)) {
    params <- modifyList(global_params, list(No.of.padded.banks = np))
    padding <- modifyList(paddings, list(rand = random_padding(nrow(Pi))))

    for (scenario in names(padding)) {
        print(paste("Padding",np,"in scenario",scenario))
        df <- df %>%
            bind_rows(run_asset_shocks_cont(params, scenario, padding[[scenario]]))
    }
}

df %>%
    mutate(num_failed = unlist(Map(length, failed_banks))) %>%
    group_by(scenario, num_padded) %>%
    summarize(failed = mean(num_failed), cost = unique(pad_total_cost)) %>%
    ggplot(aes(scenario, failed)) +
    geom_point(aes(size=cost)) +
    geom_text(aes(label=num_padded), nudge_x = 0.1)

df %>%
    mutate(num_failed = unlist(Map(length, failed_banks))) %>%
    group_by(scenario, num_padded) %>%
    summarize(failed = max(num_failed), cost = unique(pad_total_cost)) %>%
    ggplot(aes(cost, failed)) +
    geom_point(aes(color=scenario)) +
    geom_text(aes(label=num_padded), nudge_y = 0.2) +
    scale_color_colorblind()

df %>%
    mutate(num_failed = unlist(Map(length, failed_banks))) %>%
    group_by(scenario, num_padded) %>%
    ggplot(aes(num_padded, num_failed)) +
    geom_jitter(aes(color=pad_total_cost)) +
    facet_wrap(~scenario)

## Compute equity lost by failed banks
bc <- data_frame(bank = rownames(Pi),
                 id = seq_along(bank),
                 equity = cap[bank])

df %>%
    mutate(eq_loss =  unlist(Map(function (fb)
        sum(unlist(Map(function (i)
            bc %>% filter(id == i) %>% select(equity),
            fb))),
        failed_banks))) %>%
    group_by(scenario, num_padded) %>%
    summarize(loss = max(eq_loss), cost = unique(pad_total_cost)) %>%
    ggplot(aes(cost, loss)) +
    geom_point(aes(color=scenario)) +
    geom_text(aes(label=num_padded), nudge_y = 0.2) +
    scale_color_colorblind()

df %>%
    mutate(num_failed = unlist(Map(length, failed_banks))) %>%
    group_by(scenario, num_padded) %>%
    ggplot(aes(num_padded, num_failed)) +
    geom_jitter(aes(color=scenario), height=0, width=2) +
    geom_smooth(aes(group=scenario, color=scenario), se=FALSE) +
    scale_color_colorblind()

## Padding assets ... at random first
run_asset_shocks_cont_pad_assets <- function (params, scenario, padding_order="none") {
    results <- data_frame(scenario = character(0),
                          shocked_asset = integer(0),
                          failed_banks = list(),
                          params = list(),
                          num_padded = integer(0),
                          pad_total_cost = numeric(0))

    states <- cont_states(Pi, cap)
    AM <- Pi

    if (padding_order == "none") {
        ## No padding ... do nothing
        pad_total_cost <- 0
    } else {
        padded_assets <- padding_order[1:params$No.of.padded.assets]

        params$Alpha <- rep(params$Alpha, ncol(Pi))
        params$Alpha[padded_assets] <- 0
        pad_total_cost <- sum(colSums(Pi)[padded_assets])
    }
    
    for (idx in 1:ncol(AM)) {
        ## Shock asset i
        as <- shock_assets(states$assets, idx, params)
        
        sim <- propagate_asset_shock(AM, states$banks, as, params)
        ## Store simulation results
        results <- bind_rows(results,
                             list(scenario=scenario, shocked_asset=idx,
                                  failed_banks=list(which(!sim$banks$solvent)),
                                  pad_total_cost = pad_total_cost,
                                  num_padded = params$No.of.padded.assets,
                                  params = list(params)))
    }
    results
}

paddings <- list(rand = random_padding(ncol(Pi)))
                 
global_params <- modifyList(global_params, list(Reduction.factor=0.3))
df <- run_asset_shocks_cont_pad_assets(modifyList(global_params, list(No.of.padded.assets = 0)), "baseline", "none")
for (np in seq(5,100,by=5)) {
    params <- modifyList(global_params, list(No.of.padded.assets = np))
    padding <- modifyList(paddings, list(rand = random_padding(ncol(Pi))))

    for (scenario in names(padding)) {
        print(paste("Padding",np,"in scenario",scenario))
        df <- df %>%
            bind_rows(run_asset_shocks_cont_pad_assets(params, scenario, padding[[scenario]]))
    }
}
