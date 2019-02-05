@testset "Posthoc Tests" begin
    #=
    Factor A
       1       2       3       4       5
    28.2    39.6    46.3    41.0    56.3
    33.2    40.8    42.1    44.1    54.1
    36.4    37.9    43.5    46.4    59.4
    34.6    37.1    48.8    40.2    62.7
    29.1    43.6    43.7    38.6    60.0
    31.0    42.4    40.1    36.3    57.3
    =#

    observations = [28.2 39.6 46.3 41.0 56.3;
                    33.2 40.8 42.1 44.1 54.1;
                    36.4 37.9 43.5 46.4 59.4;
                    34.6 37.1 48.8 40.2 62.7;
                    29.1 43.6 43.7 38.6 60.0;
                    31.0 42.4 40.1 36.3 57.3]

    results = anova(observations)

    @testset "Tukey HSD" begin
        expected = [AnovaPosthocFactor("A",
                       [AnovaPosthocComparison((1,2),  8.15,       25, 1.27574815,  6.3884083,  0.0011293111),
                        AnovaPosthocComparison((1,3), 12.0,        25, 1.27574815,  9.4062453,  5.3421475e-6),
                        AnovaPosthocComparison((1,4),  9.0166667,  25, 1.27574815,  7.0677482,  0.00033392423),
                        AnovaPosthocComparison((1,5), 26.216667,   25, 1.27574815, 20.550033,   1.06215037e-12),
                        AnovaPosthocComparison((2,3),  3.85,       25, 1.27574815,  3.017837,   0.237621746),
                        AnovaPosthocComparison((2,4),  0.86666667, 25, 1.27574815,  0.67933994, 0.98848032),
                        AnovaPosthocComparison((2,5), 18.0666667,  25, 1.27574815, 14.161625,   3.0022479e-9),
                        AnovaPosthocComparison((3,4),  2.9833333,  25, 1.27574815,  2.3384971,  0.479109962),
                        AnovaPosthocComparison((3,5), 14.2166667,  25, 1.27574815, 11.1437878,  2.93742795e-7),
                        AnovaPosthocComparison((4,5), 17.2,        25, 1.27574815, 13.482285,   8.0080987e-9)])]

        posthocresults = tukey(results)
        hsd(results)
        honestlysignificantdifference(results)
        multiplecomparison(results)

        @test all(expected .≈ posthocresults.factorcomparisons)
    end

    @testset "Newman-Keuls" begin
        expected = [AnovaPosthocFactor("A",
                       [AnovaPosthocComparison((1,2),  8.15,       25, 1.27574815,  6.3884083,  0.000130222167),
                        AnovaPosthocComparison((1,3), 12.0,        25, 1.27574815,  9.4062453,  1.67382736e-6),
                        AnovaPosthocComparison((1,4),  9.0166667,  25, 1.27574815,  7.0677482,  0.000207629065),
                        AnovaPosthocComparison((1,5), 26.216667,   25, 1.27574815, 20.550033,   1.06215037e-12),
                        AnovaPosthocComparison((2,3),  3.85,       25, 1.27574815,  3.017837,   0.042841263),
                        AnovaPosthocComparison((2,4),  0.86666667, 25, 1.27574815,  0.67933994, 0.88116408),
                        AnovaPosthocComparison((2,5), 18.0666667,  25, 1.27574815, 14.161625,   1.8229017e-9),
                        AnovaPosthocComparison((3,4),  2.9833333,  25, 1.27574815,  2.3384971,  0.11071968),
                        AnovaPosthocComparison((3,5), 14.2166667,  25, 1.27574815, 11.1437878,  9.1097817e-8),
                        AnovaPosthocComparison((4,5), 17.2,        25, 1.27574815, 13.482285,   8.33539571e-10)])]

        posthocresults = snk(results)
        newmankeuls(results)
        studentnewmankeuls(results)
        newmankeulsmultiplerange(results)

        @test all(expected .≈ posthocresults.factorcomparisons)
    end

    #=@testset "WSD" begin
        error("Not done")
        expected = [AnovaPosthocComparison("1×2",  8.15,       25, 1.27574815,  6.3884083,  0.0011293111),
                    AnovaPosthocComparison("1×3", 12.0,        25, 1.27574815,  9.4062453,  5.3421475e-6),
                    AnovaPosthocComparison("1×4",  9.0166667,  25, 1.27574815,  7.0677482,  0.00033392423),
                    AnovaPosthocComparison("1×5", 26.216667,   25, 1.27574815, 20.550033,   1.06215037e-12),
                    AnovaPosthocComparison("2×3",  3.85,       25, 1.27574815,  3.017837,   0.237621746),
                    AnovaPosthocComparison("2×4",  0.86666667, 25, 1.27574815,  0.67933994, 0.98848032),
                    AnovaPosthocComparison("2×5", 18.0666667,  25, 1.27574815, 14.161625,   3.0022479e-9),
                    AnovaPosthocComparison("3×4",  2.9833333,  25, 1.27574815,  2.3384971,  0.479109962),
                    AnovaPosthocComparison("3×5", 14.2166667,  25, 1.27574815, 11.1437878,  2.93742795e-7),
                    AnovaPosthocComparison("4×5", 17.2,        25, 1.27574815, 13.482285,   8.0080987e-9)]

        posthocresults = wsd(results)
        whollysignificantdifference(results)
        multiplecomparisonandrange(results)

        @test all(expected .≈ posthocresults.factorcomparisons)
    end=#
end
