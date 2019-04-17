"""
Implements the firesales model of Sasidevan Vijayakumar & Nils Bertschinger.
"""

"""
"""
struct VBModel <: FinancialModel

    AM
    aᵉ
    aⁱ
    d
    α
    N
    M
    function VBModel(A::AbstractMatrix, aᵉ::AbstractVector,
                    aⁱ::AbstractVector, d::AbstractVector; α=1.0536)
        @argcheck size(A, 1) == size(aᵉ, 1) == size(aⁱ, 1) == size(d, 1)

        AM = calcAM(A, aᵉ)
        new(AM, aᵉ, aⁱ, d, α, size(A, 1), size(A, 2))
    end
end

##############################################
# Implementation of FinancialModel interface #
##############################################

numfirms(net::VBModel) = net.N

function valuation(net::VBModel, x, a)
    failing_ids = calcfailingbanks(net, x)
    Γ = fracliq(net, x, failing_ids)
    # convert boolean array of failing banks to float64 array of solvent banks
    return vcat(convert.(Float64, .!failing_ids),
                compprices(net::VBModel, x, Γ))
end

function valuation!(y, net::VBModel, x, a)
    failing_ids = calcfailingbanks(net, x)
    Γ = fracliq(net, x, failing_ids)
    # convert boolean array of failing banks to float64 array of solvent banks
    y .= vcat(convert.(Float64, .!failing_ids),
                compprices(net::VBModel, x, Γ))
end

function init(net::VBModel, a::AbstractVector)
    return vcat(ones(net.N), a)
end

function solvent(net::VBModel, x)
    return solvview(net, x)
end
##########################
# Model specific methods #
##########################

function compprices(net::VBModel, x, Γ)
    p = priceview(net, x)
    return [p[i] * exp(-net.α * Γ[i]) for i in 1:numextassets(net)]
end

function fracliq(net::VBModel, x, failing_ids)
    return sum(net.AM[failing_ids,:], dims=1)./(sum(net.AM, dims=1) .+ eps(Float64))
end

function calcfailingbanks(net::VBModel, x)
    return convert.(Bool,
            solvview(net, x) .* (net.AM * priceview(net, x) + net.aⁱ - net.d .< 0))
end

function calcAM(A::AbstractMatrix, aᵉ::AbstractVector)
    return aᵉ./(sum(A, dims=2) .+ eps(Float64)) .* A
end

function solvview(net::VBModel, x)
    return view(x, 1:net.N)
end

function priceview(net::VBModel, x)
    return view(x, net.N+1:(net.N+net.M))
end

numextassets(net::VBModel) = net.M
