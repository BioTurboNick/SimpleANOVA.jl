@testset "Contrast Tests" begin
    @testset "Single-Factor Contrast" begin
        data = [AnovaValue(     "Total", 4617.60, 39),
                AnovaResult("Treatment", 3497.60,  4, 874.40, 27.33, 0, 0), # p value and effect size unnecessary for this test
                AnovaFactor(    "Error", 1120.00, 35,  32.00)]

        result = AnovaData(data, data[1], 1, [5], 8, [AnovaFactor(r.name, r.ss, r.df, r.ms) for r ∈ data[2:2]], [data[end]], [4.00; 10.00; 11.00; 24.00; 29.00])

        @test all(contrast(result, [1, 1, 1, 2, 2]) .≈ (99.008334, 9.6573624e-12, 3.2114433))
    end
end
