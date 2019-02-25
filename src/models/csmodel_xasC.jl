"""
Implements the Cont and Shaaning firesales model from Cont & Shaanning, 2017.
"""
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
    N
    M

    """
    """
    function CSModel(Π_0::AbstractMatrix, C::AbstractVector,
                    Θ::AbstractMatrix, ϵ::AbstractVector,
                    B::AbstractVector, S_0::AbstractVector,
                    ADV::AbstractVector, σ::AbstractVector, c, τ,
                    λ_max; λ_target=0.95*λ_max, α=0.5)
        @argcheck size(Π_0, 1) == size(Θ, 1) == size(C, 1)
        @argcheck (size(B, 1) == size(S_0, 1) == size(ADV, 1) == size(σ, 1)
                    == size(Π_0, 2))
        @argcheck size(ϵ, 1) == size(Θ, 2)
        @argcheck 1. < λ_max
        @argcheck λ_target <= λ_max
        @argcheck 0. <= α <= 1.

        # loss in the value of illiquid holdings due to stress scenario
        I_ϵ = [sum(Θ[i,:]) - sum(ϵ.*Θ[i,:]) for i in 1:size(Θ,1)]
        # decrease in equity due to stress scenario
        C_ϵ = [max(C[i] - sum(ϵ.*Θ[i,:]), 0) for i in 1:size(C,1)]

        D = marketdepth(ADV, σ, c, τ)
        δ = [(1 - B[i]/S_0[i])*D[i] for i in 1:size(Π_0, 2)]

        new(copy(Π_0), C_ϵ, copy(S_0), I_ϵ, B, δ, λ_max, λ_target, α,
            size(C, 1), size(Π_0, 2))
    end
end

##############################################
# Implementation of FinancialModel interface #
##############################################

numfirms(net::CSModel) = net.N

function valuation!(y, net::CSModel, x, a)
    # compute the deleverage proportion
    Γ = delevprop(net, x)
    # compute market impact of previous sales
    Ψ = marketimpact(net, sum(Γ.*net.Π_k, dims=1)[:])
    # update market prices in asset classes due to previous sales
    net.S_k .= compPrices(net, Ψ)

    # compute total Loss due to fire sale round, using Π_{k}
    L = compLoss(net, Γ, Ψ)

    # update Π_{k+1} incorporating new asset prices
    net.Π_k .= compΠ(net, Γ, Ψ)

    # update current (equity) state, i.e. C_k
    y .= compC(net, x, L)
end

function valuation(net::CSModel, x, a)
    # compute the deleverage proportion
    Γ = delevprop(net, x)
    # compute market impact of previous sales
    Ψ = marketimpact(net, sum(Γ.*net.Π_k, dims=1)[:])
    # update market prices in asset classes due to previous sales
    net.S_k .= compPrices(net, Ψ)

    # compute total Loss due to fire sale round, using Π_{k}
    L = compLoss(net, Γ, Ψ)

    # update Π_{k+1} incorporating new asset prices
    net.Π_k .= compΠ(net, Γ, Ψ)

    return compC(net, x, L)
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
    compΠ(net, Γ, Ψ)

Computes the Pi_{k} matrix of the next round, i.e. Pi_{k+1}
"""
function compΠ(net::CSModel, Γ, Ψ)
    # update Π_k incorporating new asset prices
    # (hcat(x...)' required to transform list of arrays to 2D array)
    return hcat([(1. - Γ[i]) .* net.Π_k[i, :] .* (1. .- Ψ)
                for i in 1:numfirms(net)]...)'
end

# TODO: possibly move this calculation into compC?
"""
    compLoss(net, Γ, Ψ)

Computes the total Loss due to a fire sales round.
"""
function compLoss(net::CSModel, Γ, Ψ)
    return [(1. - (1. - net.α)*Γ[i]) * sum(net.Π_k[i] .* Ψ)
            for i in 1:numfirms(net)]
end

"""
    compC(net, x, L)

Computes the equiry of each bank of the next round given the losses due to a
round of fire sales. I.e. returns C_{k+1}
"""
function compC(net::CSModel, x::AbstractVector, L)
    # update current (equity) state, i.e. C_k
    return [max(x[i] - L[i], 0.) for i in 1:numfirms(net)]
end

"""
    compPrices(net, Ψ)

Computes the market prices due to asset sales.
"""
function compPrices(net::CSModel, Ψ)
    return [net.S_k[μ] * (1. - Ψ[μ]) for μ in 1:size(net.Π_k, 2)]
end

"""
    leverageratio(net, x)

Computes the leverage ratio of every bank.
"""
function leverageratio(net::CSModel, x)
    return [(sum(net.Π_k[i,:])+net.I_ϵ[i])/x[i] for i in 1:numfirms(net)]
end

"""
    delevprop(net, x)

Computes the deleverage proportion for each bank.
"""
function delevprop(net::CSModel, x)
    delev = zeros(size(x))
    λ = leverageratio(net, x)
    for i in 1:numfirms(net)
        # only if the leverage is higher than the system wide max leverage,
        # the institute is not illiquid,
        # (and the institute is not insolvent, -> commented out for now)
        # does it sell a proportion of its' security assets
        if ((λ[i] > net.λ_max)
            && (sum(net.Π_k[i,:]) > 0.)) #&& (x[i] > 0.)
            delev[i] = min((sum(net.Π_k[i,:]) + net.I_ϵ[i]- net.λ_target*x[i])
                            /sum(net.Π_k[i,:]), 1)
        end
    end
    return delev
end

"""
    marketimpact(net, x, q)

Computes the market impact for each asset class due to asset sales.
"""
function marketimpact(net::CSModel, q::AbstractVector)
    return [(1-net.B[i]/net.S_k[i])*(1-exp(-q[i]/net.δ[i]))
            for i in 1:(size(q, 1))]
end

"""
    marketdepth(ADV, σ, c, τ)

Computes the marketdepth for each asset class.
"""
function marketdepth(ADV::AbstractVector, σ::AbstractVector, c, τ)
    return [c*(ADV[i]/σ[i])*sqrt(τ) for i in 1:size(σ, 1)]
end

numsectypes(net::CSModel) = net.M
