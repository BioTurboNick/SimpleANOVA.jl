@testset "Contrast Tests" begin
    @testset "Single-Factor Contrast" begin
        @testset "Single Contrast" begin
            data = [AnovaValue(     "Total", 4617.60, 39),
                    AnovaResult("Treatment", 3497.60,  4, 874.40, 27.33, 0, 0), # p value and effect size unnecessary for this test
                    AnovaFactor(    "Error", 1120.00, 35,  32.00)]

            result = AnovaData(data, [AnovaFactor(r.name, r.ss, r.df, r.ms) for r ∈ data[2:2]], [5], 8, [4.00; 10.00; 11.00; 24.00; 29.00])

            expected = AnovaContrastResult(2640.2222, 35, 99.008334, 9.6573624e-12, 0.85954768)

            @test all(contrast(result, [1, 1, 1, 2, 2]) .≈ expected)
        end

        @testset "Difference Contrasts" begin
            data = [AnovaValue(     "Total", 54340.4, 39),
                    AnovaResult("Treatment", 26135.9,  3, 8711.97, 11.1199, 0, 0), # p value and effect size unnecessary for this test
                    AnovaFactor(    "Error", 28204.4, 36,  783.457)]
            result = AnovaData(data, [AnovaFactor(r.name, r.ss, r.df, r.ms) for r ∈ data[2:2]], [4], 10, [280.90, 240.75, 208.89, 239.96])

            expected = AnovaContrastResults([AnovaContrastResult(26044.011,   36, 24.931819,   1.53404354e-5, 0.63966842),
                                             AnovaContrastResult( 2665.05625, 36,  2.26777496, 0.140814326,   0.2434352),
                                             AnovaContrastResult( 9653.449,   36,  6.1608033,  0.017862323,   0.38226473)])

            @test all(differencecontrasts(result) ≈ expected)
        end

        @testset "Difference Contrasts (Reverse)" begin
            data = [AnovaValue(     "Total", 54340.4, 39),
                    AnovaResult("Treatment", 26135.9,  3, 8711.97, 11.1199, 0, 0), # p value and effect size unnecessary for this test
                    AnovaFactor(    "Error", 28204.4, 36,  783.457)]

            result = AnovaData(data, [AnovaFactor(r.name, r.ss, r.df, r.ms) for r ∈ data[2:2]], [4], 10, [280.90, 240.75, 208.89, 239.96])

            expected = AnovaContrastResults([AnovaContrastResult(  126.261778, 36,  0.120869854, 0.73011755,   0.057846875),
                                             AnovaContrastResult(26972.44225,  36, 22.9516466,   2.8424295e-5, 0.62396317),
                                             AnovaContrastResult(16120.225,    36, 10.28788115,  0.0028101042, 0.47144314)])

            @test all(differencecontrasts(result, true) ≈ expected)
        end

        @testset "Simple Contrasts (Default)" begin
            data = [AnovaValue(     "Total", 54340.4, 39),
                    AnovaResult("Treatment", 26135.9,  3, 8711.97, 11.1199, 0, 0), # p value and effect size unnecessary for this test
                    AnovaFactor(    "Error", 28204.4, 36,  783.457)]

            result = AnovaData(data, [AnovaFactor(r.name, r.ss, r.df, r.ms) for r ∈ data[2:2]], [4], 10, [280.90, 240.75, 208.89, 239.96])

            expected = AnovaContrastResults([AnovaContrastResult(16120.225, 36, 10.28788115, 0.0028101042, 0.47144314),
                                             AnovaContrastResult(51854.401, 36, 33.093329,   1.4877478e-6, 0.69207342),
                                             AnovaContrastResult(16760.836, 36, 10.6967172,  0.002368843,  0.47861035)])

            @test all(simplecontrasts(result) ≈ expected)
        end

        @testset "Simple Contrasts (Factor Level 3)" begin
            data = [AnovaValue(     "Total", 54340.4, 39),
                    AnovaResult("Treatment", 26135.9,  3, 8711.97, 11.1199, 0, 0), # p value and effect size unnecessary for this test
                    AnovaFactor(    "Error", 28204.4, 36,  783.457)]

            result = AnovaData(data, [AnovaFactor(r.name, r.ss, r.df, r.ms) for r ∈ data[2:2]], [4], 10, [280.90, 240.75, 208.89, 239.96])

            expected = AnovaContrastResults([AnovaContrastResult(51854.401, 36, 33.093329,  1.4877478e-6, 0.69207342),
                                             AnovaContrastResult(10150.596, 36,  6.4780811, 0.0153496883, 0.39051772),
                                             AnovaContrastResult( 9653.449, 36,  6.1608033, 0.017862323,  0.38226473)])

            @test all(simplecontrasts(result, 3) ≈ expected)
        end

        @testset "Repeated Contrasts" begin
            data = [AnovaValue(     "Total", 54340.4, 39),
                    AnovaResult("Treatment", 26135.9,  3, 8711.97, 11.1199, 0, 0), # p value and effect size unnecessary for this test
                    AnovaFactor(    "Error", 28204.4, 36,  783.457)]

            result = AnovaData(data, [AnovaFactor(r.name, r.ss, r.df, r.ms) for r ∈ data[2:2]], [4], 10, [280.90, 240.75, 208.89, 239.96])

            expected = AnovaContrastResults([AnovaContrastResult(16120.225, 36, 10.28788115, 0.0028101042, 0.47144314),
                                             AnovaContrastResult(10150.596, 36,  6.4780811,  0.0153496883, 0.39051772),
                                             AnovaContrastResult( 9653.449, 36,  6.1608033,  0.017862323,   0.38226473)])

            @test all(repeatedcontrasts(result) ≈ expected)
        end
    end
end
