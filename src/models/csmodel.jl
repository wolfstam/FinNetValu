"""
Implements the Cont and Shaaning firesales model from Cont & Shaanning, 2017.
"""
# TODO: - clean up redundant code
#       - refactor to incorporate Pi and C (maybe S) in x
#       - how to handle a in this model?

struct CSModel <: FinancialModel

    Π_k
    C_ϵ
    S_k
    I_ϵ
    B
    δ
    λ_max
    λ_target
    α

    """
    """
    function CSModel(Π::AbstractMatrix, C::AbstractVector,
                    Θ::AbstractMatrix, ϵ::AbstractVector,
                    B::AbstractVector, S_0::AbstractVector,
                    ADV::AbstractVector, σ::AbstractVector, c, τ,
                    λ_max; λ_target=0.95*λ_max, α=0.5)
        @argcheck size(Π, 1) == size(Θ, 1) == size(C, 1)
        @argcheck size(B, 1) == size(S_0, 1) == size(ADV, 1) == size(σ, 1) == size(Π, 2)
        @argcheck size(ϵ, 1) == size(Θ, 2)
        @argcheck 1. < λ_max
        @argcheck λ_target <= λ_max
        @argcheck 0. <= α <= 1.

        # loss in the value of illiquid holdings due to stress scenario
        I_ϵ = [sum(Θ[i,:]) - sum(ϵ.*Θ[i,:]) for i in 1:size(Θ,1)]
        # decrease in equity due to stress scenario
        C_ϵ = [max(C[i] - sum(ϵ.*Θ[i,:]), 0) for i in 1:size(C,1)]

        D = marketdepth(ADV, σ, c, τ)
        δ = [(1 - B[i]/S_0[i])*D[i] for i in 1:size(Π, 2)]

        new(copy(Π), C_ϵ, copy(S_0), I_ϵ, B, δ, λ_max, λ_target, α)
    end
end

##############################################
# Implementation of FinancialModel interface #
##############################################

numfirms(net::CSModel) = size(net.C_ϵ, 1)

function valuation!(y, net::CSModel, x, a)
    # compute the deleverage proportion
    Γ = delevprop(net, x)
    # compute market impact of previous sales
    Ψ = marketimpact(net, sum(Γ.*net.Π_k, dims=1)[:])
    # update market prices in asset classes due to previous sales
    updateprices(net, Ψ)

    # # store Π_k before updating it
    # Π_k = copy(net.Π_k)

    # compute total Loss due to fire sale round
    L = compLoss(net, Γ, Ψ)

    # update Π_k incorporating new asset prices
    # (hcat(x...)' required to transform list of arrays to 2D array)
    net.Π_k .= compPi(net, Γ, Ψ)
    # net.Π_k .= hcat([(1. - Γ[i]) .* net.Π_k[i, :] .* (1. .- Ψ) for i in 1:numfirms(net)]...)'

    # update current (equity) state, i.e. C_k
    # y .= [max(x[i] - (1. - (1. - net.α)*Γ[i]) * sum(Π_k[i] .* Ψ), 0.) for i in 1:numfirms(net)]
    y .= [max(x[i] - L[i], 0.) for i in 1:numfirms(net)]
    return y
end

function valuation(net::CSModel, x, a)
    # compute the deleverage proportion
    Γ = delevprop(net, x)
    # compute market impact of previous sales
    Ψ = marketimpact(net, sum(Γ.*net.Π_k, dims=1)[:])
    # update market prices in asset classes due to previous sales
    updateprices(net, Ψ)

    # # store Π_k before updating it
    # Π_k = copy(net.Π_k)

    # compute total Loss due to fire sale round
    L = compLoss(net, Γ, Ψ)

    # update Π_k incorporating new asset prices
    # (hcat(x...)' required to transform list of arrays to 2D array)
    net.Π_k .= compPi(net, Γ, Ψ)
    # net.Π_k .= hcat([(1. - Γ[i]) .* net.Π_k[i, :] .* (1. .- Ψ) for i in 1:numfirms(net)]...)'

    # return [max(x[i] - (1. - (1. - net.α)*Γ[i]) * sum(net.Π_k[i] .* Ψ), 0.) for i in 1:numfirms(net)]
    return [max(x[i] - L[i], 0.) for i in 1:numfirms(net)]
end

function init(net::CSModel, a)
    return copy(net.C_ϵ)
end

function solvent(net::CSModel, x)
    return x .> zero(eltype(net.C_ϵ))
end

function illiquid(net::CSModel)
    return vcat(sum(net.Π_k, dims=2)...) .== zero(eltype(net.C_ϵ))
end

##########################
# Model specific methods #
##########################
"""
Computes the Pi_{k} matrix of the next round, i.e. Pi_{k+1}
"""
function compPi(net::CSModel, Γ, Ψ)
    # update Π_k incorporating new asset prices
    # (hcat(x...)' required to transform list of arrays to 2D array)
    return hcat([(1. - Γ[i]) .* net.Π_k[i, :] .* (1. .- Ψ) for i in 1:numfirms(net)]...)'
end

"""
Computes the total Loss due to a fire sales round.
"""
function compLoss(net::CSModel, Γ, Ψ)
    return [(1. - (1. - net.α)*Γ[i]) * sum(net.Π_k[i] .* Ψ) for i in 1:numfirms(net)]
end

"""
Computes the deleverage proportion.
"""
function delevprop(net::CSModel, x)
    delev = zeros(size(x))
    for i in 1:size(net.Π_k, 1)
        # only if the leverage is higher than the system wide max leverage,
        # the institute is not illiquid,
        # and the institute is not insolvent, does it sell a proportion of its
        # security assets
        if (((sum(net.Π_k[i,:])+net.I_ϵ[i])/x[i]) > net.λ_max) && (sum(net.Π_k[i,:]) > 0.) && (x[i] > 0.)
            delev[i] = min((sum(net.Π_k[i,:]) + net.I_ϵ[i]- net.λ_target*x[i])
                            /sum(net.Π_k[i,:]), 1)
        end
    end
    return delev
end

"""
Computes the market impact for each asset class.
"""
function marketimpact(net::CSModel, q::AbstractVector)
    return [(1-net.B[i]/net.S_k[i])*(1-exp(-q[i]/net.δ[i])) for i in 1:(size(q, 1))]
end

"""
Updates the market prices due to asset sales.
"""
function updateprices(net::CSModel, Ψ)
    net.S_k .= [net.S_k[μ] * (1. - Ψ[μ]) for μ in 1:size(net.Π_k, 2)]
end

"""
Computes the marketdepth for each asset class.
"""
function marketdepth(ADV::AbstractVector, σ::AbstractVector, c, τ)
    return [c*(ADV[i]/σ[i])*sqrt(τ) for i in 1:size(σ, 1)]
end

numsectypes(net::CSModel) = size(net.Π_k, 2)
