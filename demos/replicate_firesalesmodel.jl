## Replicate two bank example from Shaanning 2017 thesis

using FinNetValu
using LaTeXStrings
using Plots
pyplot()

Π = hcat([90.; 70.])
Θ = [10. 0.; 0. 20.]
C = [4., 4.5]
λ_max = 33. #55. # not specified in thesis
λ_target = 0.97 * λ_max # also not specified in thesis
ϵ = [0.2, 0.0]
S = [1.] # not specified in thesis
B = [0.5*S[1]]
ADV = [50.]
σ = [0.02]
c = 0.4
τ = 20.

#-------------------------------------------------------------------------------
# replicate Figure 1.2.1
#-------------------------------------------------------------------------------

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
    end
    return Γ, L
end

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

p1 = plot(x*100, Γ, xlabel=L"Initial shock $\epsilon$ (%)", ylabel=L"$\Gamma$",
            title="Bank A", label=hcat(["Γ_$(i)" for i in 1:col]...));
p2 = plot(x*100, L, xlabel=L"Initial shock $\epsilon$ (%)", ylabel="Loss",
            title="Bank B", label=hcat(["L_$(i)" for i in 1:col]...));
p = plot(p1,p2, layout=(2,1), size=(600, 600))

savefig(p, "demos/plots/Shaanning_Fig_121.png")

#-------------------------------------------------------------------------------
# replicate Figure 1.2.4
#-------------------------------------------------------------------------------

"""
    visualizeinsolvencyilliquidity(;ϵ1_all=collect(0:0.01:0.45),
                                        sf=collect(0.3:0.01:1.))
Sweeps over all combinations of ϵ1_all and sf, the market depth scaling factor,
and runs the cascade of fire sales until a fixed point is reached. A heatmap is
created indicating whether a bank is insolvent (=-60), illiquid (=0) or all
right (=60).
"""
function visualizeinsolvencyilliquidity(;ϵ1_all=collect(0:0.01:0.45),
                                        sf=collect(0.3:0.01:1.))
    map = 60*ones(Int64, size(sf, 1), size(ϵ1_all, 1), 2)
    for i in 1:size(sf, 1)
        for j in 1:size(ϵ1_all, 1)
            fsmodel = FireSalesModel(Π, C, Θ, [ϵ1_all[j], 0.], B, S, ADV, σ,
                                    c*sf[i], τ, λ_max, λ_target=λ_target)
            # run fire sales cascade, store fixed point
            fp = fixvalue(fsmodel, [0., 0.])

            # compute whether banks illiquid or solvent after cascade
            # give specific code for illiquid and insolvent
            map[i, j, illiquid(fsmodel)] .= 0.
            map[i, j, .!solvent(fsmodel, fp)] .= -60.
        end
    end
    return map, ϵ1_all, sf
end

insolmap, xs, ys = visualizeinsolvencyilliquidity()

# green = all right, white = illiquid, red = insolvent
p1 = heatmap(xs*100, ys, insolmap[:,:,1], yflip=true,
            c=cgrad([:red,:white,:green]),
            colorbar = :right,
            grid=false,
            border=nothing,
            xlabel=L"Initial shock $\epsilon$ (%)",
            ylabel="Market depth scaling factor",
            title="Bank A");
p2 = heatmap(xs*100, ys, insolmap[:,:,2], yflip=true,
            c=cgrad([:red,:white,:green]),
            colorbar = :right,
            grid=false,
            border=nothing,
            xlabel=L"Initial shock $\epsilon$ (%)",
            ylabel="Market depth scaling factor",
            title="Bank B");
p = plot(p1,p2, layout=(2,1), size=(600, 600), dpi=100)

savefig(p, "demos/plots/Shaanning_Fig_124.png")
