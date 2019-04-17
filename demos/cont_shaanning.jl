using FinNetValu
using LaTeXStrings
using Plots
pyplot()

Π = hcat([90.; 70.])
Θ = [10. 0.; 0. 20.]
C = [4., 4.5]
ϵ = [0.15, 0.0]
λ_max = 33
λ_target = 0.95*λ_max
S = [1.]
B = [0.5*S[1]]
ADV = [50.]
σ = [0.02]
c = 0.4
τ = 20.
α = .5

csmodel = CSModel(Π, C, B, S, ADV, σ, c, τ, λ_max, λ_target=λ_target, α=α)
# intial a, i.e. apply shock to network
a = FinNetValu.init_a(csmodel, Θ, ϵ)
# initial current state
x = FinNetValu.init(csmodel, a)

fp = fixvalue(csmodel, a)
println(fp)
println(solvent(csmodel, fp))
println(illiquid(csmodel, fp))

# by default insolvent banks are allowed to sell the rest of their liquid assets
# to change this, set the boolean to false in the csmodel creation
# this may produce different dynamics
csmodel = CSModel(Π, C, B, S, ADV, σ, c, τ, λ_max, λ_target=λ_target, α=α, insolsell=false)
a = FinNetValu.init_a(csmodel, Θ, ϵ)
x = FinNetValu.init(csmodel, a)

fp = fixvalue(csmodel, a)
println(fp)
println(solvent(csmodel, fp))
println(illiquid(csmodel, fp))
