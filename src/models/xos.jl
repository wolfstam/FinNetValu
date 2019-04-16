import Base.show
using LinearAlgebra

"""
    XOSModel(N, Mˢ, Mᵈ, Mᵉ, d)

Financial network of `N` firms with investment portfolios `Mˢ`, `Mᵈ`
and `Mᵉ` of holding fractions in counterparties equity, debt and
external assets respectively. Nominal debt `d` is due at maturity.
"""
struct XOSModel{T1,T2,T3,U} <: FinancialModel
    N::Int64
    Mˢ::T1
    Mᵈ::T2
    Mᵉ::T3
    d::U

    function XOSModel(Mˢ::T1, Mᵈ::T2, Mᵉ::T3, d::AbstractVector) where {T1,T2,T3}
        @assert isleft_substochastic(Mˢ)
        @assert isleft_substochastic(Mᵈ)
        @assert all(d .>= 0)
        new{T1,T2,T3,typeof(d)}(length(d), Mˢ, Mᵈ, Mᵉ, d)
    end
end

function show(io::IO, net::XOSModel)
    if any(net.Mˢ .> 0) & any(net.Mᵈ .> 0)
        msg = "equity and debt"
    elseif any(net.Mˢ .> 0)
        msg = "equity"
    elseif any(net.Mᵈ .> 0)
        msg = "debt"
    else
        msg = "no"
    end
    print(io, "XOS model of N = ",
          numfirms(net), " firms with ",
          msg, " cross holdings.")
end

##############################################
# Implementation of FinancialModel interface #
##############################################

numfirms(net::XOSModel) = net.N

function valuation!(y, net::XOSModel, x, a)
    tmp = net.Mᵉ * a .+ net.Mˢ * equityview(net, x) .+ net.Mᵈ * debtview(net, x)
    equityview(net, y) .= max.(zero(eltype(x)), tmp .- net.d)
    debtview(net, y)   .= min.(net.d, tmp)
end

function valuation(net::XOSModel, x, a)
    tmp = net.Mᵉ * a .+ net.Mˢ * equityview(net, x) .+ net.Mᵈ * debtview(net, x)
    vcat(max.(zero(eltype(x)), tmp .- net.d),
         min.(net.d, tmp))
end

function fixjacobian(net::XOSModel, a, x = fixvalue(net, a))
    ## Uses analytical formulas for speedup
    ξ = solvent(net, x)
    eins = one(eltype(ξ))
    dVdx = vcat(hcat(Diagonal(ξ) * net.Mˢ, Diagonal(ξ) * net.Mᵈ),
                hcat(Diagonal(eins .- ξ) * net.Mˢ, Diagonal(eins .- ξ) * net.Mᵈ))
    dVda = vcat(Diagonal(ξ), Diagonal(1.0 .- ξ)) * net.Mᵉ
    (I - dVdx) \ Matrix(dVda) ## Note: RHS needs to be dense
end

"""
Returns [dx/dL, dx/dMˢ, dx/dMᵈ, dx/dd]
"""
function fixjacobiannet(net::XOSModel, a, x = fixvalue(net, a))
    ξ = solvent(net, x)
    eins = one(eltype(ξ))
    N = numfirms(net)
    eye = Matrix(I, N, N)
    r = debtview(net, x)
    s = equityview(net, x)

    dVdL = zeros(2*N, N, N)
    for i in 1:N
        for k in 1:N
            θ_k = r[k]/net.d[k]
            for l in 1:N
                dVdL[i, k, l] = -ξ[i]*(i==k) + ξ[i]*((i==l)*θ_k - net.Mᵈ[i, k]*θ_k)
                dVdL[i+N, k, l] = ξ[i]*(i==k) + (1-ξ[i])*((i==l)*θ_k - net.Mᵈ[i, k]*θ_k)
            end
        end
    end

    dVdMᵈ = zeros(2*N, N, N)
    dVdMˢ = zeros(2*N, N, N)
    for l in 1:N
        dVdMᵈ[:, :, l] = vcat(Diagonal(ξ) * (r[l] .* eye), Diagonal(eins .- ξ) * (r[l] .* eye))
        dVdMˢ[:, :, l] = vcat(Diagonal(ξ) * (s[l] .* eye), Diagonal(eins .- ξ) * (s[l] .* eye))
    end

    dVdd = vcat(Diagonal(-ξ) .* eye, Diagonal(ξ) .* eye)

    dVdx = vcat(hcat(Diagonal(ξ) * net.Mˢ, Diagonal(ξ) * net.Mᵈ),
             hcat(Diagonal(eins .- ξ) * net.Mˢ, Diagonal(eins .- ξ) * net.Mᵈ))

    dxdL = reshape((I - dVdx) \ reshape(dVdL, (2*N, N*N)), (2*N, N, N))
    dxdMˢ = reshape((I - dVdx) \ reshape(dVdMˢ, (2*N, N*N)), (2*N, N, N))
    dxdMᵈ = reshape((I - dVdx) \ reshape(dVdMᵈ, (2*N, N*N)), (2*N, N, N))
    dxdd = (I - dVdx) \ dVdd

    return [dxdL, dxdMˢ, dxdMᵈ, dxdd]
end

function solvent(net::XOSModel, x)
    equityview(net, x) .> zero(eltype(x))
end

function init(net::XOSModel, a)
    vcat(max.(a .- net.d, 0), net.d)
end

##########################
# Model specific methods #
##########################

"""
    equityview(net, x)

View the equity part of `x` which can be a state vector of Jacobian
matrix.
"""
function equityview end

"""
    debtview(net, x)

View the debt part of `x` which can be a state vector of Jacobian
matrix.
"""
function debtview end

equityview(net::XOSModel, x::AbstractVector) = view(x, 1:numfirms(net))
equityview(net::XOSModel, x::AbstractMatrix) = view(x, 1:numfirms(net), :)

debtview(net::XOSModel, x::AbstractVector) = begin N = numfirms(net); view(x, (N+1):(2*N)) end
debtview(net::XOSModel, x::AbstractMatrix) = begin N = numfirms(net); view(x, (N+1):(2*N), :) end
