α = FinNetValu.Spec('α', (6, 1), 1)

α2 = FinNetValu.Spec('α2', (6, 2, 3), 1)

@testset "Spec" begin
    @test @isdefined α
    @test @isdefined α2
end
