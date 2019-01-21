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

fsmodel = FinNetValu.FireSalesModel(Π, C, Θ, ϵ, B, S, ADV, σ, 1., 1., λ_max)
fsmodel2 = FinNetValu.FireSalesModel(Π, [1, -1], Θ, ϵ, B, S, ADV, σ, 1, 1, λ_max)
fsmodel3 = FinNetValu.FireSalesModel(Π, [100, -1], Θ, ϵ, B, S, ADV, σ, 1, 1, λ_max)
fsmodel4 = FinNetValu.FireSalesModel(Π2, C2, Θ2, [0., 0.], B2, S2, ADV2, σ2, 1, 1, λ_max2)
fsmodel5 = FinNetValu.FireSalesModel(Π2, C2, Θ2, [0.4, 0.4], B2, S2, ADV2, σ2, 1, 1, λ_max2)
fsmodel6 = FinNetValu.FireSalesModel(Π2, C2, Θ2, ϵ2, B2, S2, ADV2, σ2, 1, 1, λ_max2)
fsmodel7 = FinNetValu.FireSalesModel(hcat([0; 70]), C2, Θ2, [0.37, 0.], B2, S2, ADV2, σ2, 1, 1, λ_max2)

x = FinNetValu.init(fsmodel)
y = a = deepcopy(x)
x4 = FinNetValu.init(fsmodel4)
y4 = a4 = deepcopy(x4)
x5 = FinNetValu.init(fsmodel5)
x6 = FinNetValu.init(fsmodel6)
y6 = a6 = deepcopy(x6)
x7 = FinNetValu.init(fsmodel7)
y7 = a7 = deepcopy(x7)

@testset "firesales" begin
    @test @isdefined fsmodel
    @test fieldnames(typeof(fsmodel)) == (:Π_k, :C_ϵ, :S_k, :I_ϵ, :B, :δ,
                                            :λ_max, :λ_target, :α)

    @test_throws ArgumentError FinNetValu.FireSalesModel(Π', C, Θ, ϵ, B, S, ADV,
                                                            σ, 1, 1, 0.5)
    @test_throws ArgumentError FinNetValu.FireSalesModel(Π', C, Θ, ϵ, B, S, ADV,
                                                        σ, 1, 1, λ_max)
    @test_throws ArgumentError FinNetValu.FireSalesModel(Π, C, Θ', ϵ, B, S, ADV,
                                                        σ, 1, 1, λ_max)
    @test_throws ArgumentError FinNetValu.FireSalesModel(Π, C, Θ, ϵ,
                                                        rand(1:4, size(Π, 2)+1),
                                                        S, ADV, σ, 1, 1, λ_max)
    @test_throws ArgumentError FinNetValu.FireSalesModel(Π, C, Θ, ϵ, B,
                                                        rand(1:4, size(Π, 2)+1),
                                                        ADV, σ, 1, 1, λ_max)
    @test_throws ArgumentError FinNetValu.FireSalesModel(Π, C, Θ,
                                                        rand(1:4, size(Θ, 2)+1),
                                                        B, S, ADV, σ, 1, 1, λ_max)
    @test_throws ArgumentError FinNetValu.FireSalesModel(Π, C, Θ, ϵ, B, S,
                                                        [100], σ, 1, 1, λ_max)
    @test_throws ArgumentError FinNetValu.FireSalesModel(Π, C, Θ, ϵ, B, S, ADV,
                                                        [0.02, 0.02], 1, 1, λ_max)
    @test_throws ArgumentError FinNetValu.FireSalesModel(Π, C, Θ, ϵ, B, S, ADV,
                                                        σ, 1, 1, λ_max, λ_target=λ_max+1)
    @test_throws ArgumentError FinNetValu.FireSalesModel(Π, C, Θ, ϵ, B, S, ADV,
                                                        σ, 1, 1, λ_max, α=2.)

    @test fsmodel.I_ϵ' == [168. 524.]

    @test fsmodel2.C_ϵ == [0, 0]
    @test fsmodel.C_ϵ == [78, 144]
    @test fsmodel3.C_ϵ == [78, 0]

    @test fsmodel.δ == [(1-B[i]/S[i])*(ADV[i]/σ[i]) for i in 1:size(Π, 2)]

end

@testset "deleverageproportion" begin
    @test FinNetValu.delevprop(fsmodel, x) != false
    @test FinNetValu.delevprop(fsmodel4, x4) == zeros(2)
    @test FinNetValu.delevprop(fsmodel5, x5) == ones(2)
    @test FinNetValu.delevprop(fsmodel6, x6) == [22/90, 0.]
    @test FinNetValu.delevprop(fsmodel7, x7)[1] == 0.
end

@testset "marketimpact" begin
    @test FinNetValu.marketimpact(fsmodel, q) != false
    @test FinNetValu.marketimpact(fsmodel, zeros(size(Π,2))) == zeros(size(Π,2))
end

@testset "marketdepth" begin
    @test FinNetValu.marketdepth(ADV, σ, 1, 1) != false
    @test FinNetValu.marketdepth(ADV, σ, 1, 1) == ADV./σ
end

@testset "numsecurities" begin
    @test FinNetValu.numsectypes(fsmodel) == 3
end

@testset "numfirms" begin
    @test FinNetValu.numfirms(fsmodel) != false
    @test FinNetValu.numfirms(fsmodel) == 2
end

@testset "solvent" begin
    @test FinNetValu.solvent(fsmodel4, x4) != false
    @test FinNetValu.solvent(fsmodel4, x4) == [true, true]
    @test FinNetValu.solvent(fsmodel5, x5) == [false, false]
end

@testset "illiquid" begin
    @test FinNetValu.illiquid(fsmodel6) != false
    @test FinNetValu.illiquid(fsmodel7) == [true, false]
end

@testset "init" begin
    @test FinNetValu.init(fsmodel6, a6) != false
    @test FinNetValu.init(fsmodel6, a6) == fsmodel6.C_ϵ
end

@testset "updateprices" begin
    @test FinNetValu.updateprices(deepcopy(fsmodel6), [1.]) != false
    @test FinNetValu.updateprices(deepcopy(fsmodel4), [0.]) == fsmodel4.S_k
end

@testset "valuation!" begin
    @test FinNetValu.valuation!(y6, deepcopy(fsmodel6), x6, a6) != false

    @test size(FinNetValu.valuation!(y6, deepcopy(fsmodel6), x6, a6)) == size(y6)
    @test size(FinNetValu.valuation!(y, deepcopy(fsmodel), x, a)) == size(y)

    @test FinNetValu.valuation!(y6, deepcopy(fsmodel6), x6, a6) != x6
    @test FinNetValu.valuation!(y4, deepcopy(fsmodel4), x4, a4) == x4
end

@testset "valuation" begin
    @test FinNetValu.valuation(deepcopy(fsmodel6), x6, a6) != false

    @test size(FinNetValu.valuation(deepcopy(fsmodel6), x6, a6)) == size(y6)
    @test size(FinNetValu.valuation(deepcopy(fsmodel), x, a)) == size(y)

    @test FinNetValu.valuation(deepcopy(fsmodel6), x6, a6) != x6
    @test FinNetValu.valuation(deepcopy(fsmodel4), x4, a4) == x4
end
