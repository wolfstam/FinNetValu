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

csmodel = FinNetValu.CSModel(Π, C, Θ, ϵ, B, S, ADV, σ, 1., 1., λ_max)
csmodel2 = FinNetValu.CSModel(Π, [1, -1], Θ, ϵ, B, S, ADV, σ, 1, 1, λ_max)
csmodel3 = FinNetValu.CSModel(Π, [100, -1], Θ, ϵ, B, S, ADV, σ, 1, 1, λ_max)
csmodel4 = FinNetValu.CSModel(Π2, C2, Θ2, [0., 0.], B2, S2, ADV2, σ2, 1, 1, λ_max2)
csmodel5 = FinNetValu.CSModel(Π2, C2, Θ2, [0.4, 0.4], B2, S2, ADV2, σ2, 1, 1, λ_max2)
csmodel6 = FinNetValu.CSModel(Π2, C2, Θ2, ϵ2, B2, S2, ADV2, σ2, 1, 1, λ_max2)
csmodel7 = FinNetValu.CSModel(hcat([0; 70]), C2, Θ2, [0.37, 0.], B2, S2, ADV2, σ2, 1, 1, λ_max2)

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

@testset "csmodel" begin
    @test @isdefined csmodel
    @test fieldnames(typeof(csmodel)) == (:Π_k, :C_ϵ, :S_k, :I_ϵ, :B, :δ,
                                            :λ_max, :λ_target, :α)

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

    @test csmodel2.C_ϵ == [0, 0]
    @test csmodel.C_ϵ == [78, 144]
    @test csmodel3.C_ϵ == [78, 0]

    @test csmodel.δ == [(1-B[i]/S[i])*(ADV[i]/σ[i]) for i in 1:size(Π, 2)]

end

@testset "deleverageproportion" begin
    @test FinNetValu.delevprop(csmodel, x) != false
    @test FinNetValu.delevprop(csmodel4, x4) == zeros(2)
    @test FinNetValu.delevprop(csmodel5, x5) == ones(2)
    @test FinNetValu.delevprop(csmodel6, x6) == [22/90, 0.]
    @test FinNetValu.delevprop(csmodel7, x7)[1] == 0.
end

@testset "marketimpact" begin
    @test FinNetValu.marketimpact(csmodel, q) != false
    @test FinNetValu.marketimpact(csmodel, zeros(size(Π,2))) == zeros(size(Π,2))
end

@testset "marketdepth" begin
    @test FinNetValu.marketdepth(ADV, σ, 1, 1) != false
    @test FinNetValu.marketdepth(ADV, σ, 1, 1) == ADV./σ
end

@testset "numsecurities" begin
    @test FinNetValu.numsectypes(csmodel) == 3
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
    @test FinNetValu.illiquid(csmodel6) != false
    @test FinNetValu.illiquid(csmodel7) == [true, false]
end

@testset "init" begin
    @test FinNetValu.init(csmodel6, a6) != false
    @test FinNetValu.init(csmodel6, a6) == csmodel6.C_ϵ
end

@testset "updateprices" begin
    @test FinNetValu.updateprices(deepcopy(csmodel6), [1.]) != false
    @test FinNetValu.updateprices(deepcopy(csmodel4), [0.]) == csmodel4.S_k
end

@testset "valuation!" begin
    @test FinNetValu.valuation!(y6, deepcopy(csmodel6), x6, a6) != false

    @test size(FinNetValu.valuation!(y6, deepcopy(csmodel6), x6, a6)) == size(y6)
    @test size(FinNetValu.valuation!(y, deepcopy(csmodel), x, a)) == size(y)

    @test FinNetValu.valuation!(y6, deepcopy(csmodel6), x6, a6) != x6
    @test FinNetValu.valuation!(y4, deepcopy(csmodel4), x4, a4) == x4
end

@testset "valuation" begin
    @test FinNetValu.valuation(deepcopy(csmodel6), x6, a6) != false

    @test size(FinNetValu.valuation(deepcopy(csmodel6), x6, a6)) == size(y6)
    @test size(FinNetValu.valuation(deepcopy(csmodel), x, a)) == size(y)

    @test FinNetValu.valuation(deepcopy(csmodel6), x6, a6) != x6
    @test FinNetValu.valuation(deepcopy(csmodel4), x4, a4) == x4
end
