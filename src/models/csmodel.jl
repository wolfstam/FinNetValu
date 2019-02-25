"""
Implements the Cont and Shaaning firesales model from Cont & Shaanning, 2017.
"""
# TODO:
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
    N
    M

    """
    Π_0 : [N x M]
    C   : [N,]
    Θ   : [N x K]
    ϵ   : [K,]
    B   : [M,]
    S_0 : [M,]
    ADV : [M,]
    σ   : [M,]
    """
    function CSModel(Π_0::AbstractMatrix, C::AbstractVector,
                    Θ::AbstractMatrix, ϵ::AbstractVector,
                    B::AbstractVector, S_0::AbstractVector,
                    ADV::AbstractVector, σ::AbstractVector, c, τ,
                    λ_max; λ_target=0.95*λ_max, α=0.5)
        @argcheck size(Π_0, 1) == size(Θ, 1) == size(C, 1)
        @argcheck size(B, 1) == size(S_0, 1) == size(ADV, 1) == size(σ, 1) == size(Π_0, 2)
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

        new(Π_0, C_ϵ, S_0, I_ϵ, B, δ, λ_max, λ_target, α, size(C, 1), size(Π_0, 2))
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
    Ψ = marketimpact(net, x, sum(Γ .* Πview(net, x), dims=1)[:])

    # compute C_{k+1}, Π_{k+1} and S_{k+1}, where C_{k+1} and S_{k+1} make use
    # of Π_{k}
    y .= vcat(vec(hcat(compC(net, x, compLoss(net, x, Γ, Ψ)),
            compΠ(net, x, Γ, Ψ))),
            compS(net, x, Ψ))
end

function valuation(net::CSModel, x::AbstractVector, a)
    # compute the deleverage proportion
    Γ = delevprop(net, x)
    # compute market impact of previous sales, using S_{k}
    Ψ = marketimpact(net, x, sum(Γ .* Πview(net, x), dims=1)[:])

    # compute C_{k+1}, Π_{k+1} and S_{k+1}, where C_{k+1} and S_{k+1} make use
    # of Π_{k}
    return vcat(vec(hcat(compC(net, x, compLoss(net, x, Γ, Ψ)),
            compΠ(net, x, Γ, Ψ))),
            compS(net, x, Ψ))
end

function init(net::CSModel, a)
    return vcat(vec(hcat(copy(net.C_ϵ), copy(net.Π_k))), copy(net.S_k))
end

function solvent(net::CSModel, x::AbstractVector)
    return Cview(net, x) .> zero(eltype(Cview(net, x)))
end

function illiquid(net::CSModel, x::AbstractVector)
    return vcat(sum(Πview(net, x), dims=2)...) .== zero(eltype(Cview(net, x)))
end

##########################
# Model specific methods #
##########################
"""
    compΠ(net, x, Γ, Ψ)

Computes the Pi_{k} matrix of the next round, i.e. Pi_{k+1}
"""
function compΠ(net::CSModel, x::AbstractVector, Γ, Ψ)
    # update Π_k incorporating new asset prices
    # (hcat(x...)' required to transform list of arrays to 2D array)
    return hcat([(1. - Γ[i]) .* Πview(net, x)[i, :] .* (1. .- Ψ)
                for i in 1:numfirms(net)]...)'
end

"""
    compLoss(net, x, Γ, Ψ)

Computes the total Loss due to a fire sales round.
"""
function compLoss(net::CSModel, x::AbstractVector, Γ, Ψ)
    return [(1. - (1. - net.α)*Γ[i]) * sum(Πview(net, x)[i] .* Ψ)
            for i in 1:numfirms(net)]
end

"""
    compC(net, x, L)

Computes the equiry of each bank of the next round given the losses due to a
round of fire sales. I.e. returns C_{k+1}
"""
function compC(net::CSModel, x::AbstractVector, L)
    # update current (equity) state, i.e. C_k
    return [max(Cview(net, x)[i] - L[i], 0.) for i in 1:numfirms(net)]
end

"""
    compS(net, x, Ψ)

Computes the market prices due to asset sales.
"""
function compS(net::CSModel, x::AbstractVector, Ψ)
    return [Sview(net, x)[μ] * (1. - Ψ[μ]) for μ in 1:numsecassets(net)]
end

"""
    leverageratio(net, x)

Computes the leverage ratio of every bank.
"""
function leverageratio(net::CSModel, x)
    return [(sum(Πview(net, x)[i,:])+net.I_ϵ[i])/Cview(net, x)[i]
            for i in 1:numfirms(net)]
end

"""
    delevprop(net, x)

Computes the deleverage proportion for each bank.
"""
function delevprop(net::CSModel, x)
    delev = zeros(size(Cview(net, x)))
    λ = leverageratio(net, x)
    for i in 1:numfirms(net)
        # only if the leverage is higher than the system wide max leverage,
        # the institute is not illiquid,
        # (and the institute is not insolvent, -> commented out for now)
        # does it sell a proportion of its' security assets
        if ((λ[i] > net.λ_max)
            && (sum(Πview(net, x)[i,:]) > 0.)) #&& (Cview(x)[i] > 0.)
            delev[i] = min((sum(Πview(net, x)[i,:])
                            + net.I_ϵ[i] - net.λ_target * Cview(net, x)[i])
                            /sum(Πview(net, x)[i,:]),
                            1)
        end
    end
    return delev
end

"""
    marketimpact(net, x, q)

Computes the market impact for each asset class due to asset sales.
"""
function marketimpact(net::CSModel, x::AbstractVector, q::AbstractVector)
    return [(1-net.B[i]/Sview(net, x)[i])*(1-exp(-q[i]/net.δ[i]))
            for i in 1:(size(q, 1))]
end

"""
    marketdepth(ADV, σ, c, τ)

Computes the marketdepth for each asset class.
"""
function marketdepth(ADV::AbstractVector, σ::AbstractVector, c, τ)
    return [c*(ADV[i]/σ[i])*sqrt(τ) for i in 1:size(σ, 1)]
end

"""
    Πview(net, x)

View the security assets matrix in the state variable x.
"""
function Πview(net::CSModel, x::AbstractVector)
    N = numfirms(net)
    M = numsecassets(net)
    return reshape(view(x, N+1:(M*N+N)), N, M)
end

"""
    Cview(net, x)

View the equity vector in the state variable x.
"""
function Cview(net::CSModel, x::AbstractVector)
    return view(x, 1:numfirms(net))
end

"""
    Sview(net, x)

View the prices matrix of each asset in the state variable x.
"""
function Sview(net::CSModel, x::AbstractVector)
    N = numfirms(net)
    M = numsecassets(net)
    return view(x, (N+N*M)+1:(N+M*N+M))
end

numsecassets(net::CSModel) = net.M
