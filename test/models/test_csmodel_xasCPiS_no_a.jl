Π = [100 20 30;
    30 50 200]
C = [100, 200]
Θ = [20 30 40 100;
    60 410 100 10]
ϵ = [0.1, 0.0, 0.5, 0.0]
λ_max = 1.2
q = [1., 3., 1.]
S = [1., 3., 2.]
B = [0.2, 0.2, 0.3]
ADV = [100, 100, 90]
σ = [0.02, 0.03, 0.1]

Π2 = hcat([90.; 70.])
Θ2 = [10. 0.; 0. 20.]
C2 = [4., 4.5]
λ_max2 = 40.
ϵ2 = [0.2, 0.0]
S2 = [1.]
B2 = [0.2]
ADV2 = [100.]
σ2 = [0.02]

# rounded values
Π2_2 = hcat([54.47; 69.46])
C2_2 = [1.4415, 3.9603]
S2_2 = [0.99229]

csmodel = FinNetValu.CSModel(Π, C, Θ, ϵ, B, S, ADV, σ, 1., 1., λ_max)
csmodel2 = FinNetValu.CSModel(Π, [1, -1], Θ, ϵ, B, S, ADV, σ, 1, 1, λ_max)
csmodel3 = FinNetValu.CSModel(Π, [100, -1], Θ, ϵ, B, S, ADV, σ, 1, 1, λ_max)
csmodel4 = FinNetValu.CSModel(Π2, C2, Θ2, [0., 0.], B2, S2, ADV2, σ2, 1, 1,
                                λ_max2)
csmodel5 = FinNetValu.CSModel(Π2, C2, Θ2, [0.4, 0.4], B2, S2, ADV2, σ2, 1, 1,
                                λ_max2)
csmodel6 = FinNetValu.CSModel(Π2, C2, Θ2, ϵ2, B2, S2, ADV2, σ2, 1, 1, λ_max2)
csmodel7 = FinNetValu.CSModel(hcat([0; 70]), C2, Θ2, [0.37, 0.], B2, S2, ADV2,
                                σ2, 1, 1, λ_max2)
csmodel8 = FinNetValu.CSModel(Π2, C2, Θ2, [0.2, 0.], [0.5], [1.], [50.], σ2,
                                0.4, 20., 33.3; λ_target=0.95*33.3, α=0.5)
csmodel8_2 = FinNetValu.CSModel(Π2_2, C2, Θ2, [0.2, 0.], [0.5], S2_2, [50.], σ2,
                                0.4, 20., 33.3; λ_target=0.95*33.3, α=0.5)
csmodel8_3 = deepcopy(csmodel8)
csmodel9 = FinNetValu.CSModel(Π, [0, 4.5], Θ, ϵ, B, S, ADV, σ, 1, 1, λ_max)


# TODO all a's so far only dummies, need to refactor csmodel to incorporate Π
# and C as x
a = zeros(size(csmodel.C_ϵ))
x = FinNetValu.init(csmodel, a)
y = deepcopy(x)
a4 = zeros(size(csmodel4.C_ϵ))
x4 = FinNetValu.init(csmodel4, a4)
y4 = deepcopy(x4)
a5 = zeros(size(csmodel5.C_ϵ))
x5 = FinNetValu.init(csmodel5, a5)
a6 = zeros(size(csmodel5.C_ϵ))
x6 = FinNetValu.init(csmodel6, a6)
y6 = deepcopy(x6)
a7 = zeros(size(csmodel7.C_ϵ))
x7 = FinNetValu.init(csmodel7, a7)
y7 = deepcopy(x7)
a8 = zeros(size(csmodel8.C_ϵ))
x8 = FinNetValu.init(csmodel8, a8)
y8 = deepcopy(x8)
a8_2 = zeros(size(csmodel8_2.C_ϵ))
x8_2 = FinNetValu.init(csmodel8_2, a8_2)
y8_2 = deepcopy(x8_2)
a9 = zeros(size(csmodel9.C_ϵ))
x9 = FinNetValu.init(csmodel9, a)
y9 = deepcopy(x9)

# run 1 round of valuation
a8_3 = zeros(size(csmodel8_3.C_ϵ))
x8_3 = FinNetValu.valuation(csmodel8_3, FinNetValu.init(csmodel8_3, a8_3), a8_3)
y8_3 = deepcopy(x8_3)

# run 1 more round of valuation
csmodel8_4 = deepcopy(csmodel8_3)
a8_4 = zeros(size(csmodel8_4.C_ϵ))
x8_4 = FinNetValu.valuation(csmodel8_4, x8_3, a8_4)

@testset "csmodel" begin
    @test @isdefined csmodel
    @test fieldnames(typeof(csmodel)) == (:Π_k, :C_ϵ, :S_k, :I_ϵ, :B, :δ,
                                            :λ_max, :λ_target, :α, :N, :M)

    @test_throws ArgumentError FinNetValu.CSModel(Π', C, Θ, ϵ, B, S, ADV,
                                                            σ, 1, 1, 0.5)
    @test_throws ArgumentError FinNetValu.CSModel(Π', C, Θ, ϵ, B, S, ADV,
                                                        σ, 1, 1, λ_max)
    @test_throws ArgumentError FinNetValu.CSModel(Π, C, Θ', ϵ, B, S, ADV,
                                                        σ, 1, 1, λ_max)
    @test_throws ArgumentError FinNetValu.CSModel(Π, C, Θ, ϵ,
                                                        rand(1:4, size(Π, 2)+1),
                                                        S, ADV, σ, 1, 1, λ_max)
    @test_throws ArgumentError FinNetValu.CSModel(Π, C, Θ, ϵ, B,
                                                        rand(1:4, size(Π, 2)+1),
                                                        ADV, σ, 1, 1, λ_max)
    @test_throws ArgumentError FinNetValu.CSModel(Π, C, Θ,
                                                        rand(1:4, size(Θ, 2)+1),
                                                        B, S, ADV, σ, 1, 1, λ_max)
    @test_throws ArgumentError FinNetValu.CSModel(Π, C, Θ, ϵ, B, S,
                                                        [100], σ, 1, 1, λ_max)
    @test_throws ArgumentError FinNetValu.CSModel(Π, C, Θ, ϵ, B, S, ADV,
                                                        [0.02, 0.02], 1, 1, λ_max)
    @test_throws ArgumentError FinNetValu.CSModel(Π, C, Θ, ϵ, B, S, ADV,
                                                        σ, 1, 1, λ_max, λ_target=λ_max+1)
    @test_throws ArgumentError FinNetValu.CSModel(Π, C, Θ, ϵ, B, S, ADV,
                                                        σ, 1, 1, λ_max, α=2.)

    @test csmodel.I_ϵ' == [168. 524.]
    @test csmodel8.I_ϵ' == [8. 20.]

    @test csmodel2.C_ϵ == [0, 0]
    @test csmodel.C_ϵ == [78, 144]
    @test csmodel3.C_ϵ == [78, 0]
    @test csmodel8.C_ϵ == [2., 4.5]

    @test csmodel.δ == [(1-B[i]/S[i])*(ADV[i]/σ[i]) for i in 1:size(Π, 2)]
    @test round(csmodel8.δ[1], digits=2) == 2236.07
end

@testset "init" begin
    @test FinNetValu.init(csmodel6, a6) != false
    @test FinNetValu.init(csmodel6, a6) == vcat(vec(hcat(csmodel6.C_ϵ, csmodel6.Π_k)), csmodel6.S_k)
end

@testset "Piview" begin
    @test FinNetValu.Πview(csmodel, x) == csmodel.Π_k
    @test FinNetValu.Πview(csmodel6, x6) == csmodel6.Π_k
end

@testset "Cview" begin
    @test FinNetValu.Cview(csmodel, x) == csmodel.C_ϵ
    @test FinNetValu.Cview(csmodel6, x6) == csmodel6.C_ϵ
end

@testset "Sview" begin
    @test FinNetValu.Sview(csmodel, x) == csmodel.S_k
    @test FinNetValu.Sview(csmodel6, x6) == csmodel6.S_k
end

@testset "leverageratio" begin
    @test FinNetValu.leverageratio(csmodel, x) != false
    @test FinNetValu.leverageratio(csmodel8, x8) == [49, 20]
    @test round.(FinNetValu.leverageratio(csmodel8_3, x8_3), digits=2) == [43.63, 22.59]
end

@testset "deleverageproportion" begin
    @test FinNetValu.delevprop(csmodel, x) != false
    @test FinNetValu.delevprop(csmodel4, x4) == zeros(2)
    @test FinNetValu.delevprop(csmodel6, x6) == [22/90, 0.]
    @test FinNetValu.delevprop(csmodel7, x7)[1] == 0.
    @test round.(FinNetValu.delevprop(csmodel8, x8), digits=2) == [0.39, 0.]
    @test round.(FinNetValu.delevprop(csmodel8_3, x8_3), digits=2) == [0.32, 0.]
    @test round.(FinNetValu.delevprop(csmodel8_4, x8_4), digits=2) == [0.15, 0.]
    @test FinNetValu.delevprop(csmodel9, x9)[1] == 0.
end

@testset "marketimpact" begin
    @test FinNetValu.marketimpact(csmodel, x, q) != false
    @test FinNetValu.marketimpact(csmodel, x, zeros(size(Π,2))) == zeros(size(Π,2))
    @test round.(FinNetValu.marketimpact(csmodel8, x8, [34.73]), digits=5) == [0.00771]
end

@testset "marketdepth" begin
    @test FinNetValu.marketdepth(ADV, σ, 1, 1) != false
    @test FinNetValu.marketdepth(ADV, σ, 1, 1) == ADV./σ
    @test round.(FinNetValu.marketdepth([50.], [0.02], 0.4, 20.), digits=2) == [4472.14]
end

@testset "numsecassets" begin
    @test FinNetValu.numsecassets(csmodel) == 3
end

@testset "numfirms" begin
    @test FinNetValu.numfirms(csmodel) != false
    @test FinNetValu.numfirms(csmodel) == 2
end

@testset "solvent" begin
    @test FinNetValu.solvent(csmodel4, x4) != false
    @test FinNetValu.solvent(csmodel4, x4) == [true, true]
    @test FinNetValu.solvent(csmodel5, x5) == [false, false]
end

@testset "illiquid" begin
    @test FinNetValu.illiquid(csmodel6, x6) != false
    @test FinNetValu.illiquid(csmodel7, x7) == [true, false]
end

@testset "compS" begin
    @test FinNetValu.compS(deepcopy(csmodel6), x6, [1.]) != false
    @test FinNetValu.compS(deepcopy(csmodel4), x4, [0.]) == csmodel4.S_k
    @test round.(FinNetValu.compS(deepcopy(csmodel8), x8, FinNetValu.marketimpact(csmodel8, x8, [34.73])), digits=5) == [0.99229]
end

@testset "compΠ" begin
    @test FinNetValu.compΠ(deepcopy(csmodel6), x6, [0. 0.], [1.]) != false
    @test FinNetValu.compΠ(deepcopy(csmodel6), x6, [0. 0.], [0.]) == Π2
    @test round.(FinNetValu.compΠ(csmodel8, x8, [0.39 0.], [0.00771]), digits=2) == hcat([54.48; 69.46])
    @test round.(FinNetValu.compΠ(csmodel8_2, x8_2, [0.3096 0.], [0.0037]), digits=2) == hcat([37.47; 69.20])
end

@testset "compLoss" begin
    @test FinNetValu.compLoss(csmodel6, x6, [0. 0.], [0.]) != false
    @test FinNetValu.compLoss(csmodel4, x4, [0. 0.], [0.]) == [0., 0.]
    @test round.(FinNetValu.compLoss(csmodel8, x8, [0.39 0.], [0.00771]), digits=4) == [0.5586, 0.5397]
    @test round.(FinNetValu.compLoss(csmodel8_2, x8_2, [0.3096 0.], [0.0037]), digits=2) == [0.17, 0.26]
end

@testset "compC" begin
    @test FinNetValu.compC(csmodel6, x6, [0. 0.]) != false
    @test round.(FinNetValu.compC(csmodel8, x8, [0.56 0.54]), digits=2) == [1.44, 3.96]
end

@testset "valuation" begin
    @test FinNetValu.valuation(deepcopy(csmodel6), x6, a6) != false

    @test size(FinNetValu.valuation(deepcopy(csmodel6), x6, a6)) == size(x6)
    @test size(FinNetValu.valuation(deepcopy(csmodel), x, a)) == size(x)

    @test FinNetValu.valuation(deepcopy(csmodel6), x6, a6) != x6
    @test FinNetValu.valuation(deepcopy(csmodel4), x4, a4) == x4

    @test round.(FinNetValu.Cview(csmodel8, FinNetValu.valuation(deepcopy(csmodel8), x8, a8)),
                digits=2) == [1.44, 3.96]
    @test round.(FinNetValu.Cview(csmodel8_3, FinNetValu.valuation(deepcopy(csmodel8_3), x8_3, a8_3)),
                digits=2) == [1.26, 3.7]

    @test round.(FinNetValu.Πview(csmodel8, FinNetValu.valuation(deepcopy(csmodel8), x8, a8)),
                digits=2) == hcat([54.84; 69.46])
    @test round.(FinNetValu.Πview(csmodel8_3, FinNetValu.valuation(deepcopy(csmodel8_3), x8_3, a8_3)),
                digits=2) == hcat([37.42; 69.20])

    @test round.(FinNetValu.Sview(csmodel8, FinNetValu.valuation(deepcopy(csmodel8), x8, a8)),
                digits=2) == [0.99]
    @test round.(FinNetValu.Sview(csmodel8_3, FinNetValu.valuation(deepcopy(csmodel8_3), x8_3, a8_3)),
                digits=3) == [0.989]
end

@testset "valuation!" begin
    @test FinNetValu.valuation!(y6, deepcopy(csmodel6), x6, a6) != false

    FinNetValu.valuation!(y6, deepcopy(csmodel6), x6, a6)
    @test size(y6) == size(x6)
    FinNetValu.valuation!(y, deepcopy(csmodel), x, a)
    @test size(y[1]) == size(x[1])

    FinNetValu.valuation!(y8, deepcopy(csmodel8), x8, a8)
    @test round.(FinNetValu.Cview(csmodel8, y8), digits=2) == [1.44, 3.96]
    @test round.(FinNetValu.Πview(csmodel8, y8), digits=2) == hcat([54.84; 69.46])
    @test round.(FinNetValu.Sview(csmodel8, y8), digits=2) == [0.99]

    FinNetValu.valuation!(y8, deepcopy(csmodel8), y8, a8)
    @test round.(FinNetValu.Cview(csmodel8, y8), digits=2) == [1.26, 3.7]
    @test round.(FinNetValu.Πview(csmodel8, y8), digits=2) == hcat([37.42; 69.20])
    @test round.(FinNetValu.Sview(csmodel8, y8), digits=3) == [0.989]
end
