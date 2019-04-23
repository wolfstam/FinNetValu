"""
    FSModel(A, aᵉ, aⁱ, d; α=1.0536)

Implements a simple firesales model. Given the binary [NxM] adjacency matrix of
of a bipartite graph with connections between `N` banks and `M` external assets,
`aᵉ` total external assets per bank, `aⁱ` total liquid assets per bank as well
as deposits of each bank, `d` (= (aᵉ + aⁱ - capital)).
"""
struct FSModel <: FinancialModel

    AM
    aᵉ
    aⁱ
    d
    α
    N
    M
    function FSModel(A::AbstractMatrix, aᵉ::AbstractVector,
                    aⁱ::AbstractVector, d::AbstractVector; α=1.0536)
        @argcheck size(A, 1) == size(aᵉ, 1) == size(aⁱ, 1) == size(d, 1)

        # create adjacency matrix weighted by the holdings of external assets
        AM = calcAM(A, aᵉ)
        new(AM, aᵉ, aⁱ, d, α, size(A, 1), size(A, 2))
    end
end

##############################################
# Implementation of FinancialModel interface #
##############################################

numfirms(net::FSModel) = net.N

function valuation(net::FSModel, x, a)
    failing_ids = calcfailingbanks(net, x)
    Γ = fracliq(net, x, failing_ids)
    return vcat(.!failing_ids .* solvent(net, x),
                compprices(net, x, Γ))
end

function valuation!(y, net::FSModel, x, a)
    failing_ids = calcfailingbanks(net, x)
    Γ = fracliq(net, x, failing_ids)
    y .= vcat(.!failing_ids .* solvent(net, x),
              compprices(net, x, Γ))
end

function init(net::FSModel, a::AbstractVector)
    return vcat(trues(numfirms(net)), a)
end

function solvent(net::FSModel, x)
    return solvview(net, x)
end

##########################
# Model specific methods #
##########################

"""
    compprices(net, x, Γ)

Computes the fire sale prices of the assets at time t given the
prices at time t-1 and the fraction of each asset liquidated at time t, `Γ`.
"""
function compprices(net::FSModel, x, Γ)
    p = priceview(net, x)
    return [p[i] * exp(-net.α * Γ[i]) for i in 1:numextassets(net)]
end

"""
    fracliq(net, x, failing_ids)

Computes the fraction of assets to be liquidated at time t, given the list
of failing banks, `failing_ids`. In other words the fraction of assets of
failing banks.
"""
function fracliq(net::FSModel, x, failing_ids)
    return sum(net.AM[failing_ids,:], dims=1)./(sum(net.AM, dims=1) .+ eps(Float64))
end

"""
    calcfailingbanks(net, x)

Returns a boolean array indicating whether a bank has failed, i.e. become
insolvent, due to a round of fire sales.
"""
function calcfailingbanks(net::FSModel, x)
    return solvview(net, x) .* ((net.AM * priceview(net, x) + net.aⁱ - net.d) .< 0)
end

"""
    calcAM(A, aᵉ)

Creates the adjacency matrix [NxM] weighted with total external assets if a bank
has a share in these, otherwise 0. In other words the amount of asset
j ϵ (1,..,M) held by bank i ϵ (1,...,N).
"""
function calcAM(A::AbstractMatrix, aᵉ::AbstractVector)
    return aᵉ./(sum(A, dims=2) .+ eps(Float64)) .* A
end

"""
    solvview(net, x)

View of the solvency state of each bank of the state variable, `x`.
"""
function solvview(net::FSModel, x)
    return convert.(Bool, view(x, 1:numfirms(net)))
end

"""
    priceview(net, x)

View of the asset prices of the state variable `x`.
"""
function priceview(net::FSModel, x)
    return view(x, numfirms(net)+1:(numfirms(net)+numextassets(net)))
end

numextassets(net::FSModel) = net.M
