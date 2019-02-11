using mlmcTauLeap, Test, Distributions, Random

@testset "mlmcMain" begin
    Random.seed!(123)
    out, conf, timer = mlmc_main(3,2,5,5,true)
    @test abs(out - 3700) <100 #based on value from matlab code
    @test conf>0 #s.d. should be positive
    @test conf<10
end
