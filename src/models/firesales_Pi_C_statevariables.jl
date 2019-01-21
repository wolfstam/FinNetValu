"""
"""
struct FireSalesModel <: FinancialModel

    Π_k
    C_k
    S_k
    I_ϵ
    B
    δ
    λ_max
    λ_target
    α

    """
    """
    function FireSalesModel(Π::AbstractMatrix, C::AbstractVector,
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
        C_k = [max(C[i] - sum(ϵ.*Θ[i,:]), 0) for i in 1:size(C,1)]

        D = marketdepth(ADV, σ, c, τ)
        δ = [(1 - B[i]/S_0[i])*D[i] for i in 1:size(Π, 2)]

        new(Π, C_k, S_0, I_ϵ, B, δ, λ_max, λ_target, α)
    end
end

##############################################
# Implementation of FinancialModel interface #
##############################################

numfirms(net::FireSalesModel) = size(net.C_k, 1)

function valuation!(y, net::FireSalesModel, x, a)
    # update market prices due to sales in asset classes
    updateprices(net)

    # print(y)
    # print("\n \n")

    Γ = delevprop(net)
    # print(Γ)
    # print("\n \n")
    Ψ = marketimpact(net, sum(delevprop(net).*net.Π_k, dims=1)[:])
    # hcat(x...)' required to transform list of arrays to 2D array
    y[1] .= hcat([(1. - Γ[i]) .* net.Π_k[i, :] .* (1. .- Ψ) for i in 1:numfirms(net)]...)'
    y[2] .= [max(net.C_k[i] - (1. - (1. - net.α)*Γ[i]) * sum(net.Π_k[i] .* Ψ), 0.) for i in 1:numfirms(net)]

    # print(y)
    # print("\n \n \n")

    return y
end

function valuation(net::FireSalesModel, x, a)
    # update market prices due to sales in asset classes
    updateprices(net)

    Γ = delevprop(net)
    Ψ = marketimpact(net, sum(delevprop(net).*net.Π_k, dims=1)[:])
    # hcat(x...)' required to transform list of arrays to 2D array
    return Array[
            hcat([(1. - Γ[i]) .* net.Π_k[i, :] .* (1. .- Ψ) for i in 1:numfirms(net)]...)',
            [max(net.C_k[i] - (1. - (1. - net.α)*Γ[i]) * sum(net.Π_k[i] .* Ψ), 0.) for i in 1:numfirms(net)]
                 ]
end

function init(net::FireSalesModel)
    return Array[net.Π_k, net.C_k]
end

function solvent(net::FireSalesModel)
    return net.C_k .> zero(eltype(net.C_k))
end

function illiquid(net::FireSalesModel)
    return vec(sum(net.Π_k, dims=2) .== zero(eltype(net.C_k)))
end

##########################
# Model specific methods #
##########################
function delevprop(net::FireSalesModel)
    delev = zeros(size(net.C_k))
    for i in 1:size(net.Π_k, 1)
        if ((sum(net.Π_k[i,:])+net.I_ϵ[i])/net.C_k[i]) > net.λ_max
            delev[i] = min((sum(net.Π_k[i,:]) + net.I_ϵ[i]- net.λ_target*net.C_k[i])
                            /sum(net.Π_k[i,:]), 1)
        end
    end
    return delev
end

function marketimpact(net::FireSalesModel, q::AbstractVector)
    return [(1-net.B[i]/net.S_k[i])*(1-exp(-q[i]/net.δ[i])) for i in 1:(size(q, 1))]
end

"""
Updates the market prices due to asset sales.
"""
function updateprices(net::FireSalesModel)
    # q = sum(delevprop(net).*net.Π_k, dims=1)[:]
    # return q
    Ψ = marketimpact(net, sum(delevprop(net).*net.Π_k, dims=1)[:])
    net.S_k .= [net.S_k[μ] * (1. - Ψ[μ]) for μ in 1:size(net.Π_k, 2)]
end

function marketdepth(ADV::AbstractVector, σ::AbstractVector, c, τ)
    return [c*(ADV[i]/σ[i])*sqrt(τ) for i in 1:size(σ, 1)]
end

numsectypes(net::FireSalesModel) = size(net.Π_k, 2)
