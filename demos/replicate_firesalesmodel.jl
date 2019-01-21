## Replicate two bank example from Shaanning 2017 thesis

using FinNetValu
using Gadfly
using DataFrames

Π = hcat([90.; 70.])
Θ = [10. 0.; 0. 20.]
C = [4., 4.5]
λ_max = 40. #55. # random, not specified in thesis
λ_target = 0.9*λ_max
ϵ = [0.2, 0.0]
S = [1.] # random, not specified in thesis
B = [0.5*S[1]]
ADV = [50.]
σ = [0.02]
c = 0.4
τ = 20.

#-------------------------------------------------------------------------------
# replicate Figure 1.2.1

"""
    runkrounds(fsmodel, k)
Simulates k rounds of deleveraging given the initial FireSalesModel. Returns the
deleveraging proportions of 1st bank and the losses of second bank over k
rounds.
"""
function runkrounds(fsmodel::FireSalesModel, k::Int)
    Γ = ones(k)
    L = ones(k)
    # intial y
    y = ones(size(C))
    # initial current state
    x = FinNetValu.init(fsmodel, ones(numfirms(fsmodel)))

    for j in 1:k
        # deleveraging proportion before sales
        Γ[j] = delevprop(fsmodel, x)[1]
        tmp = copy(x[2])
        # current state updated after deleveraging round
        x = valuation(fsmodel, x, copy(y))
        L[j] = tmp-x[2]
        print(x)
    end
    return Γ, L
end

# fsmodel = FireSalesModel(Π, C, Θ, [0.3, 0.], B, S, ADV, σ, c, τ, λ_max)
# println(runkrounds(fsmodel, 5))

"""
    varyshockskrounds(; ϵ1_all=collect(0:0.005:0.45), k=5)
Sweeps over a range of initial shock values of asset class 1 and simulates k
rounds of deleveraging for each shock value.
"""
function varyshockskrounds(;ϵ1_all=collect(0:0.005:0.45), k=5)
    Γ = ones(size(ϵ1_all, 1), k)
    L = ones(size(ϵ1_all, 1), k)

    for i in 1:size(ϵ1_all, 1)
        fsmodel = FireSalesModel(Π, C, Θ, [ϵ1_all[i], 0.], B, S, ADV, σ, c, τ,
                                λ_max, λ_target=λ_target)

        Γ[i, :], L[i, :] = runkrounds(fsmodel, k)
    end
    return Γ, L, ϵ1_all
end

# simulate k rounds of deleveraging and fire sales
Γ, L, x = varyshockskrounds(k=5)
row, col = size(Γ)

# plot deleveraging proportions of Bank A over k rounds
data = vcat([DataFrame(x = x*100, M = Γ[:,i], yc = "Γ_$(i)") for i in 1:col]...)
fig121a = plot(data, x = :x, y = :M, color = :yc,
                Geom.line,
                Guide.xlabel("Initial shock ϵ (%)"),
                Guide.ylabel("Γ"),
                Guide.title("Bank A"));

# plot indirect contagion loss suffered by Bank B over k rounds
data = vcat([DataFrame(x = x*100, M = L[:,i], yc = "L_$(i)") for i in 1:col]...)
fig121b = plot(data, x = :x, y = :M, color = :yc,
                Geom.line,
                Guide.xlabel("Initial shock ϵ (%)"),
                Guide.ylabel("Loss"),# in monetary units"),
                Guide.title("Bank B"));

fig121 = vstack(fig121a, fig121b);
draw(SVG("demos/plots/Shaanning_Fig_121.svg", 10cm, 15cm), fig121)

#-------------------------------------------------------------------------------
# replicate Figure 1.2.4

function visualizeinsolvencyilliquidity(;ϵ1_all=collect(0:0.01:0.45),
                                        sf=collect(0.3:0.01:1.))
    map = zeros(size(sf, 1), size(ϵ1_all, 1), 2)
    for i in 1:size(sf, 1)
        for j in 1:size(ϵ1_all, 1)
            fsmodel = FireSalesModel(Π, C, Θ, [ϵ1_all[j], 0.], B, S, ADV, σ,
                                    c*sf[i], τ, λ_max, λ_target=λ_target)
            # run fire sales cascade, store fixed point
            fp = fixvalue(fsmodel, [0., 0.])

            # compute whether banks illiquid or solvent after cascade
            solv = solvent(fsmodel, fp)
            illiq = illiquid(fsmodel)
            # give specific code for illiquid and insolvent
            map[i, j, illiq] .= 1
            map[i, j, .!solv] .= 2
        end
    end
    return map
end

heatmap = visualizeinsolvencyilliquidity()

# WARNING!! couldn't manage to change yticks to lie between 0.3 and 1.
# Also: color coding is: 2 = insolvent, 1 = illiquid, 0 = all righty
fig124a = spy(heatmap[:,:,1],
                Guide.colorkey(labels=["liquid", "illiquid", "insolvent"]),
                Guide.xlabel("Initial shock ϵ (%)"),
                Guide.ylabel("Market Depth Scaling Factor"),# in monetary units"),
                Guide.title("Bank A"));
fig124b = spy(heatmap[:,:,2],
                Guide.colorkey(labels=["liquid", "illiquid", "insolvent"]),
                Guide.xlabel("Initial shock ϵ (%)"),
                Guide.ylabel("Market Depth Scaling Factor"),# in monetary units"),
                Guide.title("Bank B"));

fig124 = vstack(fig124a, fig124b);
draw(SVG("demos/plots/Shaanning_Fig_124.svg", 10cm, 15cm), fig124)
