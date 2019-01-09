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

    """
    """
    function FireSalesModel(Π::AbstractMatrix, C::AbstractVector,
                            Θ::AbstractMatrix, ϵ::AbstractVector,
                            B::AbstractVector, S_0::AbstractVector,
                            ADV::AbstractVector, σ::AbstractVector, c, τ,
                            λ_max, λ_target=0.95*λ_max)
        @argcheck size(Π, 1) == size(Θ, 1) == size(C, 1)
        @argcheck size(B, 1) == size(S_0, 1) == size(ADV, 1) == size(σ, 1) == size(Π, 2)
        @argcheck size(ϵ, 1) == size(Θ, 2)
        @argcheck λ_target <= λ_max

        C_k = [max(C[i] - sum(ϵ.*Θ[i,:]), 0) for i in 1:size(C,1)]
        I_ϵ = [sum(Θ[i,:]) - sum(ϵ.*Θ[i,:]) for i in 1:size(Θ,1)]

        D = marketdepth(ADV, σ, c, τ)
        δ = [(1 - B[i]/S_0[i])*D[i] for i in 1:size(Π, 2)]

        new(Π, C_k, S_0, I_ϵ, B, δ, λ_max, λ_target)
    end
end

##############################################
# Implementation of FinancialModel interface #
##############################################

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

function marketimpact(net::FireSalesModel, q)
    return [(1-net.B[i]/net.S_k[i])*(1-exp(-q[i]/net.δ[i])) for i in 1:(size(q, 1))]
end

function marketdepth(ADV::AbstractVector, σ::AbstractVector, c, τ)
    return [c*(ADV[i]/σ[i])*sqrt(τ) for i in 1:size(σ, 1)]
end
