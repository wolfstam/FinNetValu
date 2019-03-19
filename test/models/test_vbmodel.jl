
A = [0 1 0 0 0;
    0 0 0 1 0;
    1 0 0 0 0]
aᵉ = [0.8, 0.8, 0.8]
aⁱ = [0.2, 0.2, 0.2]
d = [0.1, 0.1, 0.1]
ϵ = 0.05
p = ones(5)
p_ϵ = (1 - ϵ).*p
AM = [0 .8 0 0 0;
    0 0 0 .8 0;
    .8 0 0 0 0]

vbmodel = FinNetValu.VBModel(A, aᵉ, aⁱ, d)
x = FinNetValu.init(vbmodel, p_ϵ)

@testset "vbmodel" begin
    @test @isdefined vbmodel
    @test fieldnames(typeof(vbmodel)) == (:AM, :aᵉ, :aⁱ, :d, :α, :N, :M)

    @test_throws ArgumentError FinNetValu.VBModel(A, [1., 1.], aⁱ, d)
    @test_throws ArgumentError FinNetValu.VBModel(A, aᵉ, [1., 1.], d)
    @test_throws ArgumentError FinNetValu.VBModel(A, aᵉ, aⁱ, [1., 1., 1., 1.])

    @test vbmodel.N == 3
    @test vbmodel.M == 5

    @test round.(vbmodel.AM, digits=4) == round.(AM, digits=4)
end

@testset "numfirms" begin
    @test FinNetValu.numfirms(vbmodel) == 3
end

@testset "numextassets" begin
    @test FinNetValu.numextassets(vbmodel) == 5
end

@testset "init" begin
    @test FinNetValu.init(vbmodel, p_ϵ) != false
    @test FinNetValu.init(vbmodel, p_ϵ)[1:vbmodel.N] == ones(3)
    @test FinNetValu.init(vbmodel, p_ϵ)[vbmodel.N+1:end] == p_ϵ
end

@testset "solvent" begin
    @test FinNetValu.solvent(vbmodel, x) == ones(3)
end

@testset "solvview" begin
    @test FinNetValu.solvview(vbmodel, x) == ones(3)
end

@testset "priceview" begin
    @test FinNetValu.priceview(vbmodel, x) == p_ϵ
end

@testset "calcAM" begin
    @test FinNetValu.calcAM(A, aᵉ) != false
    @test round.(FinNetValu.calcAM(A, aᵉ), digits=4) == round.(AM, digits=4)
end

@testset "calcfailingbanks" begin
    @test FinNetValu.calcfailingbanks(vbmodel, x) != false
    @test FinNetValu.calcfailingbanks(vbmodel,
                                    FinNetValu.init(vbmodel, ones(5))) == falses(3)
    @test FinNetValu.calcfailingbanks(vbmodel,
                                    FinNetValu.init(vbmodel, -0.5*ones(5))) == trues(3)
    @test FinNetValu.calcfailingbanks(vbmodel, vcat([1., 0., 1.], p_ϵ)) == [false, false, false]
end

@testset "fracliq" begin
    @test FinNetValu.fracliq(vbmodel, x, falses(3)) != false
    @test FinNetValu.fracliq(vbmodel,
                            FinNetValu.init(vbmodel, ones(5)),
                            falses(3)) == [0. 0. 0. 0. 0.]
    @test round.(FinNetValu.fracliq(vbmodel,
                                    FinNetValu.init(vbmodel, -0.5*ones(5)),
                                    trues(3)), digits=4) == round.([1. 1. 0. 1. 0.], digits=4)
end

@testset "compprices" begin
    @test FinNetValu.compprices(vbmodel, x, zeros(5)) != false
    @test FinNetValu.compprices(vbmodel, x, zeros(5)) == p_ϵ
end

@testset "valuation" begin
    @test FinNetValu.valuation(vbmodel, x, p_ϵ) != false
    @test FinNetValu.valuation(vbmodel, x, p_ϵ) == vcat(ones(3), p_ϵ)
end

@testset "valuation!" begin
    @test FinNetValu.valuation!(ones(size(x)), vbmodel, x, p_ϵ) != false
    @test FinNetValu.valuation!(ones(size(x)), vbmodel, x, p_ϵ) == vcat(ones(3), p_ϵ)
end
