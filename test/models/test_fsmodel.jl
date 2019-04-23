
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

fsmodel = FinNetValu.FSModel(A, aᵉ, aⁱ, d)
x = FinNetValu.init(fsmodel, p_ϵ)
x2 = deepcopy(x)
x2[1:3] .= [false, false, false]

@testset "fsmodel" begin
    @test @isdefined fsmodel
    @test fieldnames(typeof(fsmodel)) == (:AM, :aᵉ, :aⁱ, :d, :α, :N, :M)

    @test_throws ArgumentError FinNetValu.FSModel(A, [1., 1.], aⁱ, d)
    @test_throws ArgumentError FinNetValu.FSModel(A, aᵉ, [1., 1.], d)
    @test_throws ArgumentError FinNetValu.FSModel(A, aᵉ, aⁱ, [1., 1., 1., 1.])

    @test fsmodel.N == 3
    @test fsmodel.M == 5

    @test round.(fsmodel.AM, digits=4) == round.(AM, digits=4)
end

@testset "numfirms" begin
    @test FinNetValu.numfirms(fsmodel) == 3
end

@testset "numextassets" begin
    @test FinNetValu.numextassets(fsmodel) == 5
end

@testset "init" begin
    @test FinNetValu.init(fsmodel, p_ϵ) != false
    # @test FinNetValu.init(fsmodel, p_ϵ)[1:fsmodel.N] == ones(3)
    @test FinNetValu.init(fsmodel, p_ϵ)[1:fsmodel.N] == trues(3)
    @test FinNetValu.init(fsmodel, p_ϵ)[fsmodel.N+1:end] == p_ϵ
end

@testset "solvent" begin
    # @test FinNetValu.solvent(fsmodel, x) == ones(3)
    @test FinNetValu.solvent(fsmodel, x) == trues(3)
end

@testset "solvview" begin
    # @test FinNetValu.solvview(fsmodel, x) == ones(3)
    @test FinNetValu.solvview(fsmodel, x) == trues(3)
end

@testset "priceview" begin
    @test FinNetValu.priceview(fsmodel, x) == p_ϵ
end

@testset "calcAM" begin
    @test FinNetValu.calcAM(A, aᵉ) != false
    @test round.(FinNetValu.calcAM(A, aᵉ), digits=4) == round.(AM, digits=4)
end

@testset "calcfailingbanks" begin
    @test FinNetValu.calcfailingbanks(fsmodel, x) != false
    @test FinNetValu.calcfailingbanks(fsmodel,
                                    FinNetValu.init(fsmodel, ones(5))) == falses(3)
    @test FinNetValu.calcfailingbanks(fsmodel,
                                    FinNetValu.init(fsmodel, -0.5*ones(5))) == trues(3)
    @test FinNetValu.calcfailingbanks(fsmodel, vcat([1., 0., 1.], p_ϵ)) == [false, false, false]
end

@testset "fracliq" begin
    @test FinNetValu.fracliq(fsmodel, x, falses(3)) != false
    @test FinNetValu.fracliq(fsmodel,
                            FinNetValu.init(fsmodel, ones(5)),
                            falses(3)) == [0. 0. 0. 0. 0.]
    @test round.(FinNetValu.fracliq(fsmodel,
                                    FinNetValu.init(fsmodel, -0.5*ones(5)),
                                    trues(3)), digits=4) == round.([1. 1. 0. 1. 0.], digits=4)
end

@testset "compprices" begin
    @test FinNetValu.compprices(fsmodel, x, zeros(5)) != false
    @test FinNetValu.compprices(fsmodel, x, zeros(5)) == p_ϵ
end

@testset "valuation" begin
    @test FinNetValu.valuation(fsmodel, x, p_ϵ) != false
    @test FinNetValu.valuation(fsmodel, x, p_ϵ) == vcat(trues(3), p_ϵ)
    fp = FinNetValu.valuation(fsmodel, x2, p_ϵ)
    @test FinNetValu.solvent(fsmodel, fp) == [false, false, false]
end

@testset "valuation!" begin
    @test FinNetValu.valuation!(ones(size(x)), fsmodel, x, p_ϵ) != false
    @test FinNetValu.valuation!(ones(size(x)), fsmodel, x, p_ϵ) == vcat(trues(3), p_ϵ)
    fp = ones(size(x2))
    FinNetValu.valuation!(fp, fsmodel, x2, p_ϵ)
    @test FinNetValu.solvent(fsmodel, fp) == [false, false, false]
end
