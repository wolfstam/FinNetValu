## Replicate the figures from Manuel Schnabel's thesis

using FinNetValu
using DataFrames
using Missings: disallowmissing
using LaTeXStrings
using CSV
using Plots
pyplot()

"""
Takes a DataFrame object, converts it to an Array and drops support for missing.
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
    # C = parsearraytofloat(tmpC[2:end, 2])
    C = parsearraytofloat(disablemissing_converttoarray(tmpC)[2:end])

    return(C, Π, Θ, bank_ids, Π_ids, Θ_ids)
end

#------------------------------------------------------------------------------#
# The code in this block uses the data from Manuels R workspace sheet,
# "Rworkspace_before_initial_shock.RData"
#------------------------------------------------------------------------------#
# file paths to csv files
C_fp = "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/MSchnabel_preprocessing/Cont_C_MSchnabel.csv"
Π_fp = "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/MSchnabel_preprocessing/Cont_Pi_MSchnabel.csv"
Θ_fp = "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/MSchnabel_preprocessing/Cont_Theta_MSchnabel.csv"

C, Π, Θ, bank_ids, Π_ids, Θ_ids = load_csmatrices(C_fp, Π_fp, Θ_fp)

#------------------------------------------------------------------------------#
# The code in this block uses data from EBA_Cont_editWS.R, i.e. not further
# processed as in Manuel's thesis concernign differences to total exposures
#------------------------------------------------------------------------------#
#
# # file paths to csv files
# C_fp = "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/Cont_C.csv"
# Π_fp = "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/Cont_Pi.csv"
# Θ_fp = "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/Cont_Theta.csv"
#
# C, Π, Θ, bank_ids, Π_ids, Θ_ids = load_csmatrices(C_fp, Π_fp, Θ_fp)
#
#------------------------------------------------------------------------------#
# # German stress test scenario (Figure 1 of Manuel's thesis)
#------------------------------------------------------------------------------#

DE_bank_ids = [occursin("DE", x) for x in bank_ids]
# DE_Π_ids = [occursin("DE", x) for x in Π_ids]

# commercial and residential mortgage exposures are the only illiquid asset
# class within the German market
DE_Θ_ids = [occursin("DE", x) for x in Θ_ids]
# sanity check, stress test all
# DE_Θ_ids[:] .= true

# sweep over several shock sizes
DE_ϵ = collect(0:0.005:.1)
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

λ_max = 33
λ_target = 0.95*λ_max

S = ones(size(Π, 2))
B = 0.5*S

α = 0.5

# test leverage of all banks
leverage_ratio = (sum(Π, dims=2) + sum(Θ, dims=2))./C
print(leverage_ratio)
# # and linearly scale capital of banks above the threshold of Basel III
# for i in 1:size(C,1)
#     if leverage_ratio[i] > λ_target
#         C[i] *= leverage_ratio[i] * 1/λ_target
#     end
# end

num_insolv = ones(size(ϵ_array, 2))
for j in 1:size(ϵ_array, 2)
    csmodel = CSModel(Π, C, B, S, ADV, σ,
                            c, τ, λ_max, λ_target=λ_target, α=α, insolsell=false)
    # run fire sales cascade, store fixed point
    fp = fixvalue(csmodel, init_a(csmodel, Θ, ϵ_array[:, j]))
    num_insolv[j] = sum(.!solvent(csmodel, fp))
end

p = plot(DE_ϵ*100, num_insolv, xlabel="Shock size in %", ylabel="Number of insolvent banks",
            title="German stress test scenario", size=(400,300))
savefig(p, "demos/plots/Schnabel_Fig_1.png")

#------------------------------------------------------------------------------#
# sanity check

"""
    runkrounds(csmodel, k)
Simulates k rounds of valuation given the initial CSModel. Returns the
systems state and insolvency index B.
"""
function runkrounds(csmodel::CSModel, Θ, ϵ, k::Int)

    # intial a
    a = init_a(csmodel, Θ, ϵ)
    # initial current state
    x = FinNetValu.init(csmodel, a)

    for j in 1:k
        # current state updated after deleveraging round
        x = valuation(csmodel, x, a)
    end

    insol_id = .!solvent(csmodel, x)
    tmpΠ = FinNetValu.Πview(csmodel, x)
    ii_b = sum(tmpΠ[insol_id, :])/sum(Π[insol_id, :])

    return x, ii_b
end

# run for only k rounds of deleveraging
num_insolv = ones(size(ϵ_array, 2))
ii_b = ones(size(ϵ_array, 2))

for j in 1:size(ϵ_array, 2)
    csmodel = CSModel(Π, C, B, S, ADV, σ,
                            c, τ, λ_max, λ_target=λ_target, α=α, insolsell=false)
    # run fire sales cascade, store fixed point
    x, ii_b[j] = runkrounds(csmodel, Θ, ϵ_array[:, j], 20)
    num_insolv[j] = sum(.!solvent(csmodel, x))# + sum(illiquid(csmodel, x))
end

p1 = plot(DE_ϵ*100, num_insolv, xlabel="Shock size in %", ylabel="Number of insolvent banks",
            title="German stress test scenario", size=(400,300))
p2 = plot(DE_ϵ*100, ii_b, xlabel="Shock size in %", ylabel="Insolvency index B",
            title="German stress test scenario", size=(400,300))

#------------------------------------------------------------------------------#
# Figure 5

# C = C.*0.95
# leverage_ratio = (sum(Π, dims=2) + sum(Θ, dims=2))./(C.*0.95)
# print(leverage_ratio)

DE_commercial_Θ_ids = [occursin("commercial.real.estate_DE", x) for x in Θ_ids]
ϵ_array_5 = zeros(size(Θ_ids))
ϵ_array_5[DE_commercial_Θ_ids] .= 0.1

csmodel = CSModel(Π, C, B, S, ADV, σ,
                        c, τ, λ_max, λ_target=λ_target, α=α, insolsell=false)

# intial a
a = init_a(csmodel, Θ, ϵ_array_5)
# initial current state
x = FinNetValu.init(csmodel, a)

delev = zeros(10,1)
insol = zeros(10,1)
# for k in 1:10
#     delev[k] = sum(FinNetValu.delevprop(csmodel, x) .> 0.0)
#     x = valuation(csmodel, x, a)
# end

delev[1] = sum(FinNetValu.delevprop(csmodel, x, a) .> 0.0)
insol[1] = sum(.!solvent(csmodel, x))
x = valuation(csmodel, x, a)
delev[2] = sum(FinNetValu.delevprop(csmodel, x, a) .> 0.0)
insol[2] = sum(.!solvent(csmodel, x))
x = valuation(csmodel, x, a)
delev[3] = sum(FinNetValu.delevprop(csmodel, x, a) .> 0.0)
insol[3] = sum(.!solvent(csmodel, x))
tmp = FinNetValu.Cview(csmodel, x)
x = valuation(csmodel, x, a)
delev[4] = sum(FinNetValu.delevprop(csmodel, x, a) .> 0.0)
insol[4] = sum(.!solvent(csmodel, x))
x = valuation(csmodel, x, a)
delev[5] = sum(FinNetValu.delevprop(csmodel, x, a) .> 0.0)
insol[5] = sum(.!solvent(csmodel, x))
x = valuation(csmodel, x, a)
delev[6] = sum(FinNetValu.delevprop(csmodel, x, a) .> 0.0)
insol[6] = sum(.!solvent(csmodel, x))
x = valuation(csmodel, x, a)
delev[7] = sum(FinNetValu.delevprop(csmodel, x, a) .> 0.0)
insol[7] = sum(.!solvent(csmodel, x))
x = valuation(csmodel, x, a)
delev[8] = sum(FinNetValu.delevprop(csmodel, x, a) .> 0.0)
insol[8] = sum(.!solvent(csmodel, x))

print(delev)
print(insol)

#------------------------------------------------------------------------------#
# Figure 3

DE_commercial_Θ_ids = [occursin("commercial.real.estate_DE", x) for x in Θ_ids]

hsbc_id = [occursin("GB089", x) for x in bank_ids]
bnp_id = [occursin("FR013", x) for x in bank_ids]

hmap = ones(5,5)

tmpC = deepcopy(C)

for bailout in 1:10
    tmpC[hsbc_id] .+= bailout*(120000000000/5)

    for bailout2 in 1:10
        tmpC[bnp_id] .+= bailout2*(120000000000/5)
        csmodel = CSModel(Π, tmpC, B, S, ADV, σ,
                            c, τ, λ_max, λ_target=λ_target, α=α, insolsell=false)
        a = init_a(csmodel, Θ, DE_commercial_Θ_ids.*0.1)
        fp = fixvalue(csmodel, a)
        sol = .!solvent(csmodel, fp)
        hmap[bailout, bailout2] = sol[hsbc_id][1] + 2*sol[bnp_id][1]
    end
end

print(hmap)
