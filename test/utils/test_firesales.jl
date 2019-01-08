Π = [100 20 30;
    30 50 200]
C = [100, 200]
Θ = [20 30 40 100;
    60 410 100 10]
ϵ = [0.1, 0.0, 0.5, 0.0]
λ_max = 1.2
q = [1. 3. 1.]
S = [1., 3., 2.]
B = [0.2, 0.2, 0.3]
ADV = [100, 100, 90]
σ = [0.02, 0.03, 0.1]

fsmodel = FinNetValu.FireSalesModel(Π, C, Θ, ϵ, B, S, ADV, σ, 1, 1, λ_max)
fsmodel2 = FinNetValu.FireSalesModel(Π, [1, -1], Θ, ϵ, B, S, ADV, σ, 1, 1, λ_max)
fsmodel3 = FinNetValu.FireSalesModel(Π, [100, -1], Θ, ϵ, B, S, ADV, σ, 1, 1, λ_max)
fsmodel4 = FinNetValu.FireSalesModel(Π, [1000, 1000], Θ, ϵ, B, S, ADV, σ, 1, 1, λ_max)
fsmodel5 = FinNetValu.FireSalesModel(Π, C, Θ, [0.9, 0.9, 0.9, 0.9], B, S, ADV, σ, 1, 1, λ_max)

@testset "firesales" begin
    @test @isdefined fsmodel
    @test fieldnames(typeof(fsmodel)) == (:Π_k, :C_k, :S_k, :I_ϵ, :B, :δ,
                                            :λ_max, :λ_target)

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
                                                        σ, 1, 1, λ_max, λ_max+1)

    @test fsmodel.I_ϵ' == [168. 524.]

    @test fsmodel2.C_k == [0, 0]
    @test fsmodel.C_k == [78, 144]
    @test fsmodel3.C_k == [78, 0]

    @test fsmodel.δ == [(1-B[i]/S[i])*(ADV[i]/σ[i]) for i in 1:size(Π, 2)]

end

@testset "deleverageproportion" begin
    @test FinNetValu.delevprop(fsmodel) != false
    @test FinNetValu.delevprop(fsmodel4) == zeros(2)
    # # requires mock
    # @test FinNetValu.delevprop(fsmodel5) == ones(2)
end

@testset "marketimpact" begin
    @test FinNetValu.marketimpact(fsmodel, q) != false
    @test FinNetValu.marketimpact(fsmodel, zeros(size(Π,2))) == zeros(size(Π,2))
end

@testset "marketdepth" begin
    @test FinNetValu.marketdepth(ADV, σ, 1, 1) != false
    @test FinNetValu.marketdepth(ADV, σ, 1, 1) == ADV./σ
end
