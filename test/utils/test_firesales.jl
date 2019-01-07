# Lᵉ = [9., 4., 2.]
# A = [0 0.5 0;
#      0 0 0.5;
#      0.5 0 0]
# Aᵉ = [10., 5., 3.]
# enmodel = FinNetValu.EisenbergNoeModel(Lᵉ, A')

Π = [100 20 30;
    30 50 200]
C = [100, 200]
Θ = [20 30 40 100;
    60 410 100 10]
ϵ = [0.1, 0.0, 0.5, 0.0]
λ_max = 1.2
I_test = [168. 524.]
q = [1. 3.]
S = [1. 3.]
fsmodel = FinNetValu.FireSalesModel(Π, C, Θ, ϵ, λ_max)

@testset "firesales" begin
    # @test @isdefined enmodel
    # @test enmodel.name == "Eisenberg & Noe"
    # @test fieldnames(typeof(enmodel)) == (:name, :N, :A, :l, :𝕍ᵉ, :𝕍)

    @test @isdefined fsmodel
    @test fieldnames(typeof(fsmodel)) == (:Π, :C, :I_ϵ, :λ_max, :λ_target)

    @test_throws ArgumentError FinNetValu.FireSalesModel(Π', C, Θ, ϵ, λ_max)
    @test_throws ArgumentError FinNetValu.FireSalesModel(Π, C, Θ', ϵ, λ_max)
    @test_throws ArgumentError FinNetValu.FireSalesModel(Π, C, Θ, ϵ[1:end-1],
                                                        λ_max)
    @test_throws ArgumentError FinNetValu.FireSalesModel(Π, C, Θ, ϵ[1:end-1],
                                                        λ_max, λ_max+1)

    @test fsmodel.I_ϵ' == I_test
end

@testset "deleverageproportion" begin
    @test FinNetValu.delevprop(fsmodel, Π, C) != false
end

@testset "marketimpact" begin
    @test FinNetValu.marketimpact(fsmodel, q, S) != false
end
