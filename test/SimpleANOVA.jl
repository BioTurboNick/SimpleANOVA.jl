@testset "1-way ANOVAs" begin
    expected = [AnovaValue("Total", 4822.4575, 19),
                AnovaResult("A", 4684.9975, 3, 1561.66583, 181.773995, 1.43627554e-12),
                AnovaFactor("Error", 137.46, 16, 8.59125)]

    @testset "Replicate Vectors" begin
        observations = Array{Vector{Float64}, 1}(undef, 4)
        observations[1] = [60.8, 57.0, 65.0, 58.6, 61.7]
        observations[2] = [68.7, 67.7, 74.9, 66.3, 69.8]
        observations[3] = [102.6, 102.1, 100.2, 96.5, 100.4]
        observations[4] = [87.9, 84.2, 83.1, 85.7, 90.3]

        results = anova(observations)

        @test all(expected .≈ results.effects)
    end

    @testset "Replicate First Dimension" begin
        observations = [60.8 68.7 102.6 87.9;
                        57.0 67.7 102.1 84.2;
                        65.0 74.9 100.2 83.1;
                        58.6 66.3 96.5 85.7;
                        61.7 69.8 100.4 90.3]

        results = anova(observations, [replicates])

        @test all(expected .≈ results.effects)
    end

end
