using FinNetValu
import Distributions.Uniform
using Plots
pyplot(legend=false)

## global params
N = 25  # banks
M = 30  # assets
# bank parameters
liq_asset_frac = 0.2
ext_asset_frac = 0.8
cap_frac = 0.1
# Asset shock parameters
loss = 0.3
α = 1.0536
# # Bailout parameters
# no_of_padded_banks = 5
# new_capital_frac = 1

#---------------------------------
for i in 1:50
    p = 1.0/N
    A = erdosrenyi(N, M, p)
    # randomly choose an asset and shock by 30%
    a = ones(M)
    assetID = rand(1:M)
    a[assetID] *= loss

    total_assets = ones(N)
    aᵉ = total_assets * ext_asset_frac
    aⁱ = total_assets * liq_asset_frac
    cap = total_assets * cap_frac
    d = aᵉ + aⁱ - cap
    # Fix banks with no connections to external assets
    aⁱ[vec(sum(A, dims=2) .== 0)] .= 1

    ## Test with fixvalue() function
    # create network model
    fsmodel = FSModel(A, aᵉ, aⁱ, d, α=α)
    # run fixed point calculation
    fp = fixvalue(fsmodel, a, m=0)
    print(N - sum(solvent(fsmodel, fp)), " ")
end

#-------------------------------

total_assets = ones(N)
aᵉ = total_assets * ext_asset_frac
aⁱ = total_assets * liq_asset_frac
cap = total_assets * cap_frac
d = aᵉ + aⁱ - cap
nsim = 500

# create bipartite erdos renyi graphs representing 30 asset classes and 25 banks
# vary the average node degree between 0 and 14, for each degree value create
# 50 graphs
k_array = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 8, 10, 12, 14]
failed_banks = ones(length(k_array), nsim)
for i in 1:length(k_array)
    println(k_array[i])
    p = k_array[i]/N
    for j in 1:nsim
        A = erdosrenyi(N, M, p)
        aⁱ = total_assets * liq_asset_frac
        # Fix banks with no connections to external assets
        aⁱ[vec(sum(A, dims=2) .== 0)] .= 1

        # create network model
        fsmodel = FSModel(A, aᵉ, aⁱ, d, α=α)
        # randomly choose an asset and shock by 30%
        a = ones(M)
        a[rand(1:M, 1)] *= loss
        # run fixed point calculation
        fp = fixvalue(fsmodel, a, m=0)

        failed_banks[i, j] = N - sum(solvent(fsmodel, fp))
    end
end
# println(failed_banks)

p = scatter(repeat(k_array, outer=(1,nsim)) .+
            rand(Uniform(-0.15, 0.15), length(k_array), nsim),
            failed_banks .+
            rand(Uniform(-0.4, 0.4), length(k_array), nsim),
            xlabel = "Average degree <k>",
            ylabel = "Failed Banks",
            ylim=[-1., 25.],
            color="black",
            legend=false,
            tickfont=font(25),
            guidefont=font(25),
            size=(1400, 900));
savefig(p, "demos/plots/firesales_erdosrenyi.png")

# h = histogram(failed_banks[4,:]);
# savefig(h, "demos/plots/firesales_erdosrenyi_hist.png")
#
# println(failed_banks[2,:])
