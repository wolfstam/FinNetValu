## Replicate the figures from Manuel Schnabel's thesis

using FinNetValu
using DataFrames
using Missings: disallowmissing
using LaTeXStrings
using CSV
using Plots
pyplot()

"""
Takes a DataFrame object, convert it to an Array and drop support for missing.
"""
function disablemissing_converttoarray(a::DataFrame)
    return disallowmissing(convert(Array, a))
end

"""
Parses an Array of floats represented in strings to floats.
"""
function parsearraytofloat(a::Array)
    return map(x->tryparse(Float64, x), a)
end

"""
    load_csmatrices(C_csvpath, Π_csvpath, Θ_csvpath)

Receives the file paths to the csv files containing the information for the
C, Π and Θ matrices of the Cont and Shaanning model. The function returns these
matrices as well as string arrays containing the codes for the banks, security
assets and illiquid assets.
"""
function load_csmatrices(C_csvpath, Π_csvpath, Θ_csvpath)

    # load Capital, Pi and Theta data file
    tmpC = CSV.read(C_csvpath; header=false)
    tmpPi = CSV.read(Π_csvpath; header=false)
    tmpTheta = CSV.read(Θ_csvpath; header=false)

    # BANKING_GROUP_CODE from EBA 2011 results
    # enforcing String array representation without Missings as datatype
    bank_ids = disallowmissing(tmpPi[2:end, 1])
    # convert DataFrame object to Array, removing support for Missings
    Π_ids = vec(disablemissing_converttoarray(tmpPi[1:1, 2:end]))
    Θ_ids = vec(disablemissing_converttoarray(tmpTheta[1:1, 2:end]))

    # extract only value matrix, convert to array, remove support for Missing,
    # parse strings to floats (same for Θ)
    Π = parsearraytofloat(disablemissing_converttoarray(tmpPi[2:end, 2:end]))
    Θ = parsearraytofloat(disablemissing_converttoarray(tmpTheta[2:end, 2:end]))
    C = parsearraytofloat(tmpC[2:end, 2])

    return(C, Π, Θ, bank_ids, Π_ids, Θ_ids)
end

# file paths to csv files
C_fp = "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/Cont_C.csv"
Π_fp = "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/Cont_Pi.csv"
Θ_fp = "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/Cont_Theta.csv"

C, Π, Θ, bank_ids, Π_ids, Θ_ids = load_csmatrices(C_fp, Π_fp, Θ_fp)

#------------------------------------------------------------------------------#
## German stress test scenario (Figure 1)

DE_bank_ids = [occursin("DE", x) for x in bank_ids]
DE_Π_ids = [occursin("DE", x) for x in Π_ids]
# commercial and residential mortgage exposures are the only illiquid asset
# class within the German market
DE_Θ_ids = [occursin("DE", x) for x in Θ_ids]

# C_DE = C[DE_bank_ids]
# Π_DE = Π[DE_bank_ids, DE_Π_ids]
# Θ_DE = Θ[DE_bank_ids, DE_Θ_ids]

# sweep over several shock sizes
DE_ϵ = collect(0:0.005:0.1)
ϵ_array = repeat(zeros(size(Θ_ids)), 1, size(DE_ϵ, 1))
# apply shocks only to german illiquid assets
for i in 1:size(DE_ϵ, 1)
    ϵ_array[DE_Θ_ids, i] = DE_ϵ[i]*ones(sum(DE_Θ_ids))
end

# marketdepth was 10^11 across all asset classes therefore:
τ = 1
c = 1
σ = ones(size(Π, 2))
ADV = 10^11*ones(size(Π, 2))

#TODO: double check this (see bottom page 17 of Draft)
# leverage constraints
λ_max = 30
λ_target = 0.95*λ_max

S = ones(size(Π, 2))
B = 0.5*S

#TODO: not specified in draft
α = 0.5

# test leverage of all banks
leverage_ratio = (sum(Π, dims=2) + sum(Θ, dims=2))./C
# linearly scale capital of banks above the threshold of Basel III
for i in 1:size(C,1)
    if leverage_ratio[i] > λ_target
        C[i] *= leverage_ratio[i] * 1/λ_target
    end
end

# num_insolv = ones(size(ϵ_array, 2))
# for j in 1:size(ϵ_array, 2)
#     csmodel = CSModel(Π, C, Θ, ϵ_array[:, j], B, S, ADV, σ,
#                             c, τ, λ_max, λ_target=λ_target, α=α)
#     # run fire sales cascade, store fixed point
#     fp = fixvalue(csmodel, [0., 0.])
#     num_insolv[j] = sum(.!solvent(csmodel, fp))
# end
#
# p1 = plot(DE_ϵ*100, num_insolv, xlabel="Shock size in %", ylabel="Number of failed banks",
#             title="German stress test scenario", size=(400,300))

#------------------------------------------------------------------------------#
# sanity check

"""
    runkrounds(csmodel, k)
Simulates k rounds of valuation given the initial CSModel. Returns the
systems state.
"""
function runkrounds(csmodel::CSModel, k::Int)

    # intial y
    y = ones(size(C))
    # initial current state
    x = FinNetValu.init(csmodel, ones(numfirms(csmodel)))

    for j in 1:k
        # current state updated after deleveraging round
        x = valuation(csmodel, x, copy(y))
    end
    return x
end

# run for only k rounds of deleveraging
num_insolv = ones(size(ϵ_array, 2))
for j in 1:size(ϵ_array, 2)
    csmodel = CSModel(Π, C, Θ, ϵ_array[:, j], B, S, ADV, σ,
                            c, τ, λ_max, λ_target=λ_target, α=α)
    # run fire sales cascade, store fixed point
    x = runkrounds(csmodel, 1)
    num_insolv[j] = sum(x .<= 0.0)
end

p1 = plot(DE_ϵ*100, num_insolv, xlabel="Shock size in %", ylabel="Number of failed banks",
            title="German stress test scenario", size=(400,300))


# run for only k rounds of deleveraging
num_insolv = ones(8)
for j in 1:8
    csmodel = CSModel(Π, C, Θ, ϵ_array[:, 10], B, S, ADV, σ,
                            c, τ, λ_max, λ_target=λ_target, α=α)
    x = runkrounds(csmodel, j)
    num_insolv[j] = sum(x .<= 0.0)
end

csmodel = CSModel(Π, C, Θ, 0.1*ones(size(Θ_ids)), B, S, ADV, σ,
                        c, τ, λ_max, λ_target=λ_target, α=α)
x = runkrounds(csmodel, 2)
print(sum(x .<= 0.0))
