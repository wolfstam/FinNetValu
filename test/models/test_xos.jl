# using FinNetValu
using LinearAlgebra

xosmodel = FinNetValu.XOSModel([[0.0 0.2 0.3 0.1];
                   [0.2 0.0 0.2 0.1];
                   [0.1 0.1 0.0 0.3];
                   [0.1 0.1 0.1 0.0]],
                  [[0.0 0.0 0.1 0.0];
                   [0.0 0.0 0.0 0.1];
                   [0.1 0.0 0.0 0.1];
                   [0.0 0.1 0.0 0.0]],
                  LinearAlgebra.I,
                  [0.8, 0.8, 0.8, 0.8])

a = [2.0, 0.5, 0.6, 0.6]
# J = FinNetValu.fixjacobiannet(xosmodel, a)

xosmodel_ex6_9 = FinNetValu.XOSModel(
                   [[0.0 0.0 0.3];
                   [0.4 0.0 0.1];
                   [0.0 0.3 0.0]],
                   [[0. 0. 0.];
                    [0. 0. 0.2];
                    [0. 0. 0.]],
                  LinearAlgebra.I,
                  [4., 1., (11*47+16)/47])
a_ex6_9 = [1.0, 3., 11.]
# J_ex6_9 = FinNetValu.fixjacobiannet(xosmodel_ex6_9, a_ex6_9)

A = hcat(vcat(2 .*ones(4,4), 3 .*ones(4,4)),
         vcat(4 .*ones(4,4), 5 .*ones(4,4)),
         vcat(6 .*ones(4,4), 7 .*ones(4,4)))

# @testset "dsdMˢview" begin
#     @test FinNetValu.dsdMˢview(xosmodel, A) != false
#     @test FinNetValu.dsdMˢview(xosmodel, A) == 2 .*ones(4,4)
# end
#
# @testset "drdMˢview" begin
#     @test FinNetValu.drdMˢview(xosmodel, A) != false
#     @test FinNetValu.drdMˢview(xosmodel, A) == 3 .*ones(4,4)
# end
#
# @testset "dsdMᵈview" begin
#     @test FinNetValu.dsdMᵈview(xosmodel, A) != false
#     @test FinNetValu.dsdMᵈview(xosmodel, A) == 4 .*ones(4,4)
# end
#
# @testset "drdMᵈview" begin
#     @test FinNetValu.drdMᵈview(xosmodel, A) != false
#     @test FinNetValu.drdMᵈview(xosmodel, A) == 5 .*ones(4,4)
# end
#
# @testset "dsddview" begin
#     @test FinNetValu.dsddview(xosmodel, A) != false
#     @test FinNetValu.dsddview(xosmodel, A) == 6 .*ones(4,4)
# end
#
# @testset "drddview" begin
#     @test FinNetValu.drddview(xosmodel, A) != false
#     @test FinNetValu.drddview(xosmodel, A) == 7 .*ones(4,4)
# end

@testset "fixjacobiannet" begin
    @test FinNetValu.fixjacobiannet(xosmodel, a) != false
    @test size(FinNetValu.fixjacobiannet(xosmodel, a)[1]) == (8, 4, 4)
    @test size(FinNetValu.fixjacobiannet(xosmodel, a)[2]) == (8, 4, 4)
    @test size(FinNetValu.fixjacobiannet(xosmodel, a)[3]) == (8, 4)
    @test size(FinNetValu.fixjacobiannet(xosmodel_ex6_9, a_ex6_9)[1]) == (6, 3, 3)

    @test FinNetValu.fixjacobiannet(xosmodel_ex6_9, a_ex6_9)[3][2, 3] == 0.1/0.97
    @test FinNetValu.fixjacobiannet(xosmodel_ex6_9, a_ex6_9)[3][3, 3] == -0.94/0.97
    @test FinNetValu.fixjacobiannet(xosmodel_ex6_9, a_ex6_9, zeros(6))[3] == zeros(6, 3)
end
