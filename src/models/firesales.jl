"""
"""
struct FireSalesModel <: FinancialModel

    Π
    C
    I_ϵ
    λ_max
    λ_target

    """
    """
    function FireSalesModel(Π::AbstractMatrix, C::AbstractVector,
                            Θ::AbstractMatrix, ϵ::AbstractVector,
                            λ_max, λ_target = 0.95*λ_max)
        @argcheck size(Π, 1) == size(Θ, 1) == size(C, 1)
        @argcheck size(ϵ, 1) == size(Θ, 2)
        @argcheck λ_target <= λ_max

        I_ϵ = sum(Θ, dims=2) - sum(repeat(ϵ, outer=[1, 2])'.*Θ, dims=2)

        new(Π, C, I_ϵ, λ_target, λ_max)
    end
end

##############################################
# Implementation of FinancialModel interface #
##############################################

##########################
# Model specific methods #
##########################
function delevprop(net::FireSalesModel, Π, C)
end

function marketimpact(nett::FireSalesModel, q, S)
end
