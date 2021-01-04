@testset "ANOVA Tests" begin
    @testset "Multidimensional Array" begin
        @testset "1-way ANOVA" begin
            expected = [AnovaValue( "Total", 4822.4575, 19),
                        AnovaResult(    "A", 4684.9975,  3, 1561.66583, 181.773995, 1.43627554e-12, 0.9644332),
                        AnovaFactor("Error",  137.46,   16,    8.59125)]

            @testset "Replicate Vectors" begin
                observations = Array{Vector{Float64}, 1}(undef, 4)
                observations[1] = [ 60.8,  57.0,  65.0, 58.6,  61.7]
                observations[2] = [ 68.7,  67.7,  74.9, 66.3,  69.8]
                observations[3] = [102.6, 102.1, 100.2, 96.5, 100.4]
                observations[4] = [ 87.9,  84.2,  83.1, 85.7,  90.3]

                results = anova(observations)

                @test all(expected .≈ results.effects)
            end

            @testset "Replicate First Dimension" begin
                observations = [60.8 68.7 102.6 87.9;
                                57.0 67.7 102.1 84.2;
                                65.0 74.9 100.2 83.1;
                                58.6 66.3  96.5 85.7;
                                61.7 69.8 100.4 90.3]

                results = anova(observations)

                @test all(expected .≈ results.effects)
            end
        end

        @testset "2-way ANOVA" begin
            observations1 = Array{Vector{Float64}, 2}(undef, 2, 2)
            observations1[1,1] = [16.5, 18.4, 12.7, 14.0, 12.8]
            observations1[1,2] = [14.5, 11.0, 10.8, 14.3, 10.0]
            observations1[2,1] = [39.1, 26.2, 21.3, 35.8, 40.2]
            observations1[2,2] = [32.0, 23.8, 28.8, 25.0, 29.3]

            observations2 = cat(hcat([16.5, 18.4, 12.7, 14.0, 12.8], [39.1, 26.2, 21.3, 35.8, 40.2]),
                                hcat([14.5, 11.0, 10.8, 14.3, 10.0], [32.0, 23.8, 28.8, 25.0, 29.3]), dims = 3)

            @testset "Fixed-effects (Model I)" begin
                expected = [AnovaValue( "Total", 1827.6975, 19),
                            AnovaResult(    "A",   70.3125,  1,   70.3125,  3.0706495,  0.098856175,  0.025621074),
                            AnovaResult(    "B", 1386.1125,  1, 1386.1125, 60.533556,   7.9430782e-7, 0.73663535),
                            AnovaResult("A × B",    4.9005,  1,    4.9005,  0.21401199, 0.64987001,  -0.0097253817),
                            AnovaFactor("Error",  366.372,  16,   22.89825)]

                @testset "Replicate Vectors" begin
                    results = anova(observations1)
                    @test all(expected .≈ results.effects)
                end

                @testset "Replicate First Dimension" begin
                    results = anova(observations2)
                    @test all(expected .≈ results.effects)
                end
            end

            @testset "Random-effects (Model II)" begin
                expected = [AnovaValue( "Total", 1827.6975, 19),
                            AnovaResult(    "A",   70.3125,  1,   70.3125,  14.3480257,  0.16431864,  0.039894829),
                            AnovaResult(    "B", 1386.1125,  1, 1386.1125, 282.85124,    0.037808553, 0.842402253),
                            AnovaResult("A × B",    4.9005,  1,    4.9005,   0.21401199, 0.64987001, -0.021953683),
                            AnovaFactor("Error",  366.372,  16,   22.89825)]

                @testset "Replicate Vectors" begin
                    results = anova(observations1, [random, random])
                    @test all(expected .≈ results.effects)
                end

                @testset "Replicate First Dimension" begin
                    results = anova(observations2, [random, random])
                    @test all(expected .≈ results.effects)
                end
            end

            @testset "Mixed-effects (Model III)" begin
                expected = [AnovaValue( "Total", 1827.6975, 19),
                            AnovaResult(    "A",   70.3125,  1,   70.3125, 14.3480257,  0.16431864,   0.0203534123),
                            AnovaResult(    "B", 1386.1125,  1, 1386.1125, 60.533556,   7.9430782e-7, 0.84834776),
                            AnovaResult("A × B",    4.9005,  1,    4.9005,  0.21401199, 0.64987001,  -0.0112002576),
                            AnovaFactor("Error",  366.372,  16,   22.89825)]

                @testset "Replicate Vectors" begin
                    results = anova(observations1, [random])
                    @test all(expected .≈ results.effects)
                end

                @testset "Replicate First Dimension" begin
                    results = anova(observations2, [random])
                    @test all(expected .≈ results.effects)
                end
            end
        end

        @testset "2-way ANOVA without replicates" begin
            observations1 = Array{Vector{Float64}, 2}(undef, 3, 4)
            observations1[1,1] = [123]
            observations1[1,2] = [138]
            observations1[1,3] = [110]
            observations1[1,4] = [151]
            observations1[2,1] = [145]
            observations1[2,2] = [165]
            observations1[2,3] = [140]
            observations1[2,4] = [167]
            observations1[3,1] = [156]
            observations1[3,2] = [176]
            observations1[3,3] = [185]
            observations1[3,4] = [175]

            observations2 = [123 138 110 151; 145 165 140 167; 156 176 185 175]

            @testset "Fixed-effects (Model I)" begin
                expected = [AnovaValue(     "Total", 5594.9167,  11),
                            AnovaResult(        "A", 1116.91667,  3,  372.30556,  2.6316513, 0.144561215,  0.120720546),
                            AnovaResult(        "B", 3629.1667,   2, 1814.58333, 12.8264284, 0.0068110648, 0.58333253),
                            AnovaFactor("Remainder",  848.83333,  6,  141.472222)]

                @testset "Replicate Vectors" begin
                    results = anova(observations1)
                    @test all(expected .≈ results.effects)
                end

                @testset "Replicate First Dimension" begin
                    results = anova(observations2, hasreplicates = false)
                    @test all(expected .≈ results.effects)
                end
            end

            @testset "Random-effects (Model II)" begin
                expected = [AnovaValue(     "Total", 5594.9167,  11),
                            AnovaResult(        "A", 1116.91667,  3,  372.30556,  2.6316513, 0.144561215,  0.120849876),
                            AnovaResult(        "B", 3629.1667,   2, 1814.58333, 12.8264284, 0.0068110648, 0.65695214),
                            AnovaFactor("Remainder",  848.83333,  6,  141.472222)]

                @testset "Replicate Vectors" begin
                    results = anova(observations1, [random, random])
                    @test all(expected .≈ results.effects)
                end

                @testset "Replicate First Dimension" begin
                    results = anova(observations2, [random, random], hasreplicates = false)
                    @test all(expected .≈ results.effects)
                end
            end

            @testset "Mixed-effects (Model III)" begin
                expected = [AnovaValue(     "Total", 5594.9167,  11),
                            AnovaResult(        "A", 1116.91667,  3,  372.30556,  2.6316513, 0.144561215,  0.093461097),
                            AnovaResult(        "B", 3629.1667,   2, 1814.58333, 12.8264284, 0.0068110648, 0.67741863),
                            AnovaFactor("Remainder",  848.83333,  6,  141.472222)]

                @testset "Replicate Vectors" begin
                    results = anova(observations1, [random])
                    @test all(expected .≈ results.effects)
                end

                @testset "Replicate First Dimension" begin
                    results = anova(observations2, [random], hasreplicates = false)
                    @test all(expected .≈ results.effects)
                end
            end
        end

        @testset "1-way ANOVA with 1 nested factor" begin
            expected = [AnovaValue( "Total", 71.666667, 11),
                        AnovaResult(    "A", 61.166667,  2, 30.583333, 61.166667,   0.0037032412, 0.83371824),
                        AnovaResult(    "B",  1.5,       3,  0.5,       0.33333333, 0.80220227,  -0.083140878),
                        AnovaFactor("Error",  9.0,       6,  1.5)]

            @testset "Replicate Vectors" begin
                observations = Array{Vector{Float64}, 2}(undef, 2, 3)
                observations[1,1] = [102, 104]
                observations[1,2] = [108, 110]
                observations[1,3] = [104, 106]
                observations[2,1] = [103, 104]
                observations[2,2] = [109, 108]
                observations[2,3] = [105, 107]

                results = anova(observations, [nested])
                @test all(expected .≈ results.effects)
            end

            @testset "Replicate First Dimension" begin
                observations = cat(hcat([102, 104], [103, 104]),
                                   hcat([108, 110], [109, 108]),
                                   hcat([104, 106], [105, 107]), dims = 3)

                results = anova(observations, [nested])
                @test all(expected .≈ results.effects)
            end
        end

        @testset "1-way ANOVA with 2 nested factors" begin
            expected = [AnovaValue( "Total", 101.58753,   47),
                        AnovaResult(    "A",   7.5108875,  2, 3.7554438, 2.0223607, 0.27789847,   0.036705354),
                        AnovaResult(    "B",   5.5708812,  3, 1.8569604, 0.99348408, 0.45717102, -0.0007064175),
                        AnovaResult(    "C",  11.2148375,  6, 1.8691396, 0.87059412, 0.52587983, -0.032229523),
                        AnovaFactor("Error",  77.290925,  36, 2.14697014)]

            @testset "Replicate Vectors" begin
                observations = Array{Vector{Float64}, 3}(undef, 2, 2, 3)
                observations[1,1,1] = [3.17, 4.41, 1.81, 1.74]
                observations[2,1,1] = [2.81, 4.98, 2.62, 2.53]
                observations[1,2,1] = [3.0,  3.02, 4.73, 1.77]
                observations[2,2,1] = [4.13, 0.71, 3.18, 3.34]
                observations[1,1,2] = [2.42, 1.28, 1.4,  2.56]
                observations[2,1,2] = [1.26, 1.08, 1.42, 0.85]
                observations[1,2,2] = [0.36, 0.35, 2.64, 3.75]
                observations[2,2,2] = [3.86, 2.53, 3.97, 3.03]
                observations[1,1,3] = [5.24, 2.24, 0.18, 3.06]
                observations[2,1,3] = [4.04, 4.14, 0.33, 4.61]
                observations[1,2,3] = [6.06, 1.61, 2.25, 2.44]
                observations[2,2,3] = [0.02, 3.95, 0.87, 2.0]


                results = anova(observations, [nested, nested])
                @test all(expected .≈ results.effects)
            end

            @testset "Replicate First Dimension" begin
                observations = cat(cat(hcat([3.17, 4.41, 1.81, 1.74], [2.81, 4.98, 2.62, 2.53]),
                                       hcat([3.0,  3.02, 4.73, 1.77], [4.13, 0.71, 3.18, 3.34]), dims = 3),
                                   cat(hcat([2.42, 1.28, 1.4,  2.56], [1.26, 1.08, 1.42, 0.85]),
                                       hcat([0.36, 0.35, 2.64, 3.75], [3.86, 2.53, 3.97, 3.03]), dims = 3),
                                   cat(hcat([5.24, 2.24, 0.18, 3.06], [4.04, 4.14, 0.33, 4.61]),
                                       hcat([6.06, 1.61, 2.25, 2.44], [0.02, 3.95, 0.87, 2.0]), dims = 3), dims = 4)

                results = anova(observations, [nested, nested])
                @test all(expected .≈ results.effects)
            end
        end

        @testset "2-way ANOVA with 1 nested factor, no replicates" begin
            observations1 = Array{Vector{Float64}, 3}(undef, 5, 2, 2)
            observations1[1,1,1] = [16.5]
            observations1[2,1,1] = [18.4]
            observations1[3,1,1] = [12.7]
            observations1[4,1,1] = [14.0]
            observations1[5,1,1] = [12.8]
            observations1[1,2,1] = [39.1]
            observations1[2,2,1] = [26.2]
            observations1[3,2,1] = [21.3]
            observations1[4,2,1] = [35.8]
            observations1[5,2,1] = [40.2]
            observations1[1,1,2] = [14.5]
            observations1[2,1,2] = [11.0]
            observations1[3,1,2] = [10.8]
            observations1[4,1,2] = [14.3]
            observations1[5,1,2] = [10.0]
            observations1[1,2,2] = [32.0]
            observations1[2,2,2] = [23.8]
            observations1[3,2,2] = [28.8]
            observations1[4,2,2] = [25.0]
            observations1[5,2,2] = [29.3]

            observations2 = cat(hcat([16.5, 18.4, 12.7, 14.0, 12.8], [39.1, 26.2, 21.3, 35.8, 40.2]),
                                hcat([14.5, 11.0, 10.8, 14.3, 10.0], [32.0, 23.8, 28.8, 25.0, 29.3]), dims = 3)

            @testset "Fixed-effects (Model I)" begin
               expected = [AnovaValue( "Total", 1827.6975, 19),
                           AnovaResult(    "A",   70.3125,  1,   70.3125,  3.0706495,  0.098856175,  0.025621074),
                           AnovaResult(    "B", 1386.1125,  1, 1386.1125, 60.533556,   7.9430782e-7, 0.73663535),
                           AnovaResult("A × B",    4.9005,  1,    4.9005,  0.21401199, 0.64987001,  -0.0097253817),
                           AnovaFactor(    "C",  366.372,  16,   22.89825)]

               @testset "Replicate Vectors" begin
                   results = anova(observations1, [nested])
                   @test all(expected .≈ results.effects)
               end

               @testset "Replicate First Dimension" begin
                   results = anova(observations2, [nested], hasreplicates = false)
                   @test all(expected .≈ results.effects)
               end
            end

            @testset "Random-effects (Model II)" begin
               expected = [AnovaValue( "Total", 1827.6975, 19),
                           AnovaResult(    "A",   70.3125,  1,   70.3125,  14.3480257,  0.16431864,  0.039894829),
                           AnovaResult(    "B", 1386.1125,  1, 1386.1125, 282.85124,    0.037808553, 0.84240225),
                           AnovaResult("A × B",    4.9005,  1,    4.9005,   0.21401199, 0.64987001, -0.021953683),
                           AnovaFactor(    "C",  366.372,  16,   22.89825)]

               @testset "Replicate Vectors" begin
                   results = anova(observations1, [nested, random, random])
                   @test all(expected .≈ results.effects)
               end

               @testset "Replicate First Dimension" begin
                   results = anova(observations2, [nested, random, random], hasreplicates = false)
                   @test all(expected .≈ results.effects)
               end
            end

            @testset "Mixed-effects (Model III)" begin
               expected = [AnovaValue( "Total", 1827.6975, 19),
                           AnovaResult(    "A",   70.3125,  1,   70.3125, 14.3480257,  0.16431864,   0.0203534123),
                           AnovaResult(    "B", 1386.1125,  1, 1386.1125, 60.533556,   7.9430782e-7, 0.84834776),
                           AnovaResult("A × B",    4.9005,  1,    4.9005,  0.21401199, 0.64987001,  -0.0112002576),
                           AnovaFactor(    "C",  366.372,  16,   22.89825)]

               @testset "Replicate Vectors" begin
                   results = anova(observations1, [nested, random])
                   @test all(expected .≈ results.effects)
               end

               @testset "Replicate First Dimension" begin
                   results = anova(observations2, [nested, random], hasreplicates = false)
                   @test all(expected .≈ results.effects)
               end
            end
        end

        @testset "2-way ANOVA with 1 nested factor, replicates" begin
            observations1 = Array{Vector{Float64}, 3}(undef, 5, 2, 2)
            observations1[1,1,1] = [16.5, 16.3]
            observations1[2,1,1] = [18.4, 18.9]
            observations1[3,1,1] = [12.7, 13.4]
            observations1[4,1,1] = [14.0, 13.8]
            observations1[5,1,1] = [12.8, 12.5]
            observations1[1,2,1] = [39.1, 39.0]
            observations1[2,2,1] = [26.2, 25.4]
            observations1[3,2,1] = [21.3, 21.5]
            observations1[4,2,1] = [35.8, 35.0]
            observations1[5,2,1] = [40.2, 40.5]
            observations1[1,1,2] = [14.5, 15.0]
            observations1[2,1,2] = [11.0, 11.3]
            observations1[3,1,2] = [10.8, 10.3]
            observations1[4,1,2] = [14.3, 14.3]
            observations1[5,1,2] = [10.0, 10.4]
            observations1[1,2,2] = [32.0, 32.8]
            observations1[2,2,2] = [23.8, 23.0]
            observations1[3,2,2] = [28.8, 28.3]
            observations1[4,2,2] = [25.0, 24.0]
            observations1[5,2,2] = [29.3, 28.4]

            observations2 = cat(cat(hcat([16.5, 16.3], [18.4, 18.9], [12.7, 13.4], [14.0, 13.8], [12.8, 12.5]),
                                    hcat([39.1, 39.0], [26.2, 25.4], [21.3, 21.5], [35.8, 35.0], [40.2, 40.5]), dims = 3),
                                cat(hcat([14.5, 15.0], [11.0, 11.3], [10.8, 10.3], [14.3, 14.3], [10.0, 10.4]),
                                    hcat([32.0, 32.8], [23.8, 23.0], [28.8, 28.3], [25.0, 24.0], [29.3, 28.4]), dims = 3), dims = 4)

            @testset "Fixed-effects (Model I)" begin
               expected = [AnovaValue( "Total", 3608.231, 39),
                           AnovaResult(    "A",  144.4,    1,  144.4,      3.0540486,   0.099697096,   0.026567684),
                           AnovaResult(    "B", 2692.881,  1, 2692.881,   56.954221,    1.17181416e-6, 0.72372875),
                           AnovaResult("A × B",   11.236,  1,   11.236,    0.237640515, 0.63252737,   -0.0098605873),
                           AnovaResult(    "C",  756.504, 16,   47.2815, 294.588785,    3.4334068e-20, 0.2578079),
                           AnovaFactor("Error",    3.21,  20,    0.1605)]

               @testset "Replicate Vectors" begin
                   results = anova(observations1, [nested])
                   @test all(expected .≈ results.effects)
               end

               @testset "Replicate First Dimension" begin
                   results = anova(observations2, [nested])
                   @test all(expected .≈ results.effects)
               end
            end

            @testset "Random-effects (Model II)" begin
                expected = [AnovaValue( "Total", 3608.231, 39),
                            AnovaResult(    "A",  144.4,    1,   144.4,     12.8515486,   0.17318116,    0.04139207),
                            AnovaResult(    "B", 2692.881,  1,  2692.881,  239.66545,     0.04106525,    0.83354988),
                            AnovaResult("A × B",   11.236,  1,    11.236,    0.237640515, 0.63252737,   -0.0224084264),
                            AnovaResult(    "C",  756.504, 16,    47.2815, 294.588785,    3.4334068e-20, 0.146468694),
                            AnovaFactor("Error",    3.21,  20,     0.1605)]

               @testset "Replicate Vectors" begin
                   results = anova(observations1, [nested, random, random])
                   @test all(expected .≈ results.effects)
               end

               @testset "Replicate First Dimension" begin
                   results = anova(observations2, [nested, random, random])
                   @test all(expected .≈ results.effects)
               end
            end

            @testset "Mixed-effects (Model III)" begin
                expected = [AnovaValue( "Total", 3608.231, 39),
                            AnovaResult(    "A",  144.4,    1,  144.4,     12.8515486,   0.17318116,    0.0211334126),
                            AnovaResult(    "B", 2692.881,  1, 2692.881,   56.954221,    1.17181416e-6, 0.83972464),
                            AnovaResult("A × B",   11.236,  1,   11.236,    0.237640515, 0.63252737,   -0.0114409964),
                            AnovaResult(    "C",  756.504, 16,   47.2815, 294.588785,    3.4334068e-20, 0.1495640769),
                            AnovaFactor("Error",    3.21,  20,    0.1605)]

               @testset "Replicate Vectors" begin
                   results = anova(observations1, [nested, random])
                   @test all(expected .≈ results.effects)
               end

               @testset "Replicate First Dimension" begin
                   results = anova(observations2, [nested, random])
                   @test all(expected .≈ results.effects)
               end
            end
        end

        @testset "3-way ANOVA" begin
            observations1 = Array{Vector{Float64}, 3}(undef, 2, 3, 3)
            observations1[1,1,1] = [1.9, 1.8, 1.6, 1.4]
            observations1[1,1,2] = [2.1, 2.0, 1.8, 2.2]
            observations1[1,1,3] = [1.1, 1.2, 1.0, 1.4]
            observations1[1,2,1] = [2.3, 2.1, 2.0, 2.6]
            observations1[1,2,2] = [2.4, 2.6, 2.7, 2.3]
            observations1[1,2,3] = [2.0, 2.1, 1.9, 2.2]
            observations1[1,3,1] = [2.9, 2.8, 3.4, 3.2]
            observations1[1,3,2] = [3.6, 3.1, 3.4, 3.2]
            observations1[1,3,3] = [2.9, 2.8, 3.0, 3.1]
            observations1[2,1,1] = [1.8, 1.7, 1.4, 1.5]
            observations1[2,1,2] = [2.3, 2.0, 1.9, 1.7]
            observations1[2,1,3] = [1.4, 1.0, 1.3, 1.2]
            observations1[2,2,1] = [2.4, 2.7, 2.4, 2.6]
            observations1[2,2,2] = [2.0, 2.3, 2.1, 2.4]
            observations1[2,2,3] = [2.4, 2.6, 2.3, 2.2]
            observations1[2,3,1] = [3.0, 3.1, 3.0, 2.7]
            observations1[2,3,2] = [3.1, 3.0, 2.8, 3.2]
            observations1[2,3,3] = [3.2, 2.9, 2.8, 2.9]

            observations2 = cat(cat(hcat([1.9, 1.8, 1.6, 1.4], [1.8, 1.7, 1.4, 1.5]),
                                    hcat([2.3, 2.1, 2.0, 2.6], [2.4, 2.7, 2.4, 2.6]),
                                    hcat([2.9, 2.8, 3.4, 3.2], [3.0, 3.1, 3.0, 2.7]), dims = 3),
                                cat(hcat([2.1, 2.0, 1.8, 2.2], [2.3, 2.0, 1.9, 1.7]),
                                    hcat([2.4, 2.6, 2.7, 2.3], [2.0, 2.3, 2.1, 2.4]),
                                    hcat([3.6, 3.1, 3.4, 3.2], [3.1, 3.0, 2.8, 3.2]), dims = 3),
                                cat(hcat([1.1, 1.2, 1.0, 1.4], [1.4, 1.0, 1.3, 1.2]),
                                    hcat([2.0, 2.1, 1.9, 2.2], [2.4, 2.6, 2.3, 2.2]),
                                    hcat([2.9, 2.8, 3.0, 3.1], [3.2, 2.9, 2.8, 2.9]), dims = 3), dims = 4)

            @testset "All Fixed" begin
                expected = [AnovaValue(     "Total", 30.355,       71),
                            AnovaResult(        "A",  1.8175,       2,  0.90875,      24.475062,    2.71459995e-8, 0.057358295),
                            AnovaResult(        "B", 24.655833,     2, 12.3279167,   332.02369,     4.5550981e-31, 0.8088138),
                            AnovaResult(        "C",  0.008888889,  1,  0.008888889,   0.239401496, 0.62662035,   -0.0009292123),
                            AnovaResult(    "A × B",  1.10166667,   4,  0.27541667,    7.4177057,   7.7516681e-5,  0.031361677),
                            AnovaResult(    "A × C",  0.37027778,   2,  0.18513889,    4.9862843,   0.0102990988,  0.009739973),
                            AnovaResult(    "B × C",  0.17527778,   2,  0.08763889,    2.3603491,   0.104056404,   0.0033238381),
                            AnovaResult("A × B × C",  0.220555556,  4,  0.055138889,   1.4850374,   0.21958095,    0.002370253),
                            AnovaFactor(    "Error",  2.005,       54,  0.03712963)]


                @testset "Replicate Vectors" begin
                    results = anova(observations1)
                    @test all(expected .≈ results.effects)
                end

                @testset "Replicate First Dimension" begin
                    results = anova(observations2)
                    @test all(expected .≈ results.effects)
                end
            end

            @testset "All Random" begin
                expected = [AnovaValue(     "Total", 30.355,       71),
                            AnovaResult(        "A",  1.8175,       2,  0.90875,      2.2415211,   0.211958974,  0.035028707),
                            AnovaResult(        "B", 24.655833,     2, 12.3279167,   40.036536,    0.0022084069, 0.83651337),
                            AnovaResult(        "C",  0.008888889,  1,  0.008888889,  0.040842374, 0.85716044,  -0.0096850896),
                            AnovaResult(    "A × B",  1.10166667,   4,  0.27541667,   4.9949622,   0.074190827,  0.045989677),
                            AnovaResult(    "A × C",  0.37027778,   2,  0.18513889,   3.3576826,   0.139349696,  0.018094299),
                            AnovaResult(    "B × C",  0.17527778,   2,  0.08763889,   1.58942065,  0.31046402,   0.0045235748),
                            AnovaResult("A × B × C",  0.220555556,  4,  0.055138889,  1.4850374,   0.21958095,   0.0075199598),
                            AnovaFactor(    "Error",  2.005,       54,  0.03712963)]

                @testset "Replicate Vectors" begin
                    results = anova(observations1, [random, random, random])
                    @test all(expected .≈ results.effects)
                end

                @testset "Replicate First Dimension" begin
                    results = anova(observations2, [random, random, random])
                    @test all(expected .≈ results.effects)
                end
            end

            @testset "1 random, 2 fixed" begin
                expected = [AnovaValue(     "Total", 30.355,       71),
                            AnovaResult(        "A",  1.8175,       2,  0.90875,       4.9084771,   0.16924835,   0.047662611),
                            AnovaResult(        "B", 24.655833,     2, 12.3279167,   140.667195,    0.0070587972, 0.80623914),
                            AnovaResult(        "C",  0.008888889,  1,  0.008888889,   0.239401496, 0.62662035,  -0.0018601531),
                            AnovaResult(    "A × B",  1.10166667,   4,  0.27541667,    4.9949622,   0.074190827,  0.029018388),
                            AnovaResult(    "A × C",  0.37027778,   2,  0.18513889,    4.9862843,   0.0102990988, 0.0194980636),
                            AnovaResult(    "B × C",  0.17527778,   2,  0.08763889,    2.3603491,   0.104056404,  0.006653859),
                            AnovaResult("A × B × C",  0.220555556,  4,  0.055138889,   1.4850374,   0.21958095,   0.0047449151),
                            AnovaFactor(    "Error",  2.005,       54,  0.03712963)]

                @testset "Replicate Vectors" begin
                    results = anova(observations1, [random])
                    @test all(expected .≈ results.effects)
                end

                @testset "Replicate First Dimension" begin
                    results = anova(observations2, [random])
                    @test all(expected .≈ results.effects)
                end
            end

            @testset "2 random, 1 fixed" begin
                expected = [AnovaValue(     "Total", 30.355,       71),
                            AnovaResult(        "A",  1.8175,       2,  0.90875,       2.2415211,   0.211958974, 0.023628362),
                            AnovaResult(        "B", 24.655833,     2, 12.3279167,   140.667195,    0.0070587972, 0.8619071),
                            AnovaResult(        "C",  0.008888889,  1,  0.008888889,   0.101426307, 0.78030599, -0.0036968215),
                            AnovaResult(    "A × B",  1.10166667,   4,  0.27541667,    4.9949622,   0.074190827, 0.031022005),
                            AnovaResult(    "A × C",  0.37027778,   2,  0.18513889,    3.3576826,   0.139349696, 0.012205379),
                            AnovaResult(    "B × C",  0.17527778,   2,  0.08763889,    2.3603491,   0.104056404, 0.0071132844),
                            AnovaResult("A × B × C",  0.220555556,  4,  0.055138889,   1.4850374,   0.21958095,  0.0050725346),
                            AnovaFactor(    "Error",  2.005,       54,  0.03712963)]

                @testset "Replicate Vectors" begin
                    results = anova(observations1, [random, random])
                    @test all(expected .≈ results.effects)
                end

                @testset "Replicate First Dimension" begin
                    results = anova(observations2, [random, random])
                    @test all(expected .≈ results.effects)
                end
            end
        end

        @testset "Repeated measures ANOVA with no among-subject factors" begin
            observations1 = Array{Vector{Float64}, 2}(undef, 2, 2)
            observations1[1,1] = [8]
            observations1[2,1] = [7]
            observations1[3,1] = [1]
            observations1[4,1] = [6]
            observations1[1,2] = [9]
            observations1[2,2] = [5]
            observations1[3,2] = [2]
            observations1[4,2] = [5]
            observations1[1,3] = [6]
            observations1[2,3] = [2]
            observations1[3,3] = [3]
            observations1[4,3] = [8]
            observations1[1,4] = [5]
            observations1[2,4] = [3]
            observations1[3,4] = [1]
            observations1[4,4] = [9]
            observations1[1,5] = [8]
            observations1[2,5] = [4]
            observations1[3,5] = [5]
            observations1[4,5] = [8]
            observations1[1,6] = [7]
            observations1[2,6] = [5]
            observations1[3,6] = [6]
            observations1[4,6] = [7]
            observations1[1,7] = [10]
            observations1[2,7] = [2]
            observations1[3,7] = [7]
            observations1[4,7] = [2]
            observations1[1,8] = [12]
            observations1[2,8] = [6]
            observations1[3,8] = [8]
            observations1[4,8] = [1]

            observations2 = [8  9  6  5  8  7 10 12
                             7  5  2  3  4  5  2  6
                             1  2  3  1  5  6  7  8 
                             6  5  8  9  8  7  2  1]
            
            # NOT CORRECT VALUES YET
            expected = [AnovaValue( "Total", 1827.6975, 19),
                        AnovaResult(    "A",   70.3125,  1,   70.3125,  3.0706495,  0.098856175,  0.025621074),
                        AnovaResult(    "B", 1386.1125,  1, 1386.1125, 60.533556,   7.9430782e-7, 0.73663535),
                        AnovaResult("A × B",    4.9005,  1,    4.9005,  0.21401199, 0.64987001,  -0.0097253817),
                        AnovaFactor("Error",  366.372,  16,   22.89825)]
 
            @testset "Replicate Vectors" begin
                results = anova(observations1, [within, subject])
                @test all(expected .≈ results.effects)
            end

            @testset "Replicate First Dimension" begin
                results = anova(observations2, [within, subject])
                @test all(expected .≈ results.effects)
            end
        end

        observations = cat([  1  38  10;  43  20   9; 15  20   6; 40  28  20;  8 11  27; 17  17   9;  30  15  19;  34  27  12;  34  24 12; 26  23 21; 1  28  33;   7  26  23; 22  34  21; 30  32  17; 40 24  15; 15  29  13;  20  30  16;   9  24  17;  14  34 19; 15  23 29],
        [  6  -5 -14;  30 -12 -10; 15 -15 -16; 30  -4 -10; 12 -2   5; 17  -6  -6;  21  -2 -20;  23  -7 -12;  20 -10 -9; 27 -15 -6; -19 -13  -2; -18 -16 -17; -8 -23 -19; -6 -22 -11; -6 -9 -10; -9 -18 -17; -17 -17  -4; -12 -15  -4; -11 -14 -1; -6 -15 -1],
        [  5   4  -2;   8   4 -13; 12   6   1; 19   0   2;  8  6  -5; 15   6 -13;  21  16   3;  28   7   2;  26  12  4; 27  14  0; -10  13   9;   6  19   5;  4  14   0;  3  21   4;  0 19   2;  4   7   8;   9  12  10;  -5  18   8;   7  20 12; 13  15 10], dims = 3)
        AnovaFactor("B", 21628.677777777782, 2.0, 10814.338888888891)
        AnovaFactor("C", 2092.344444444445, 2.0, 1046.1722222222224)
        AnovaFactor("B × C", 2624.422222222216, 4.0, 656.105555555554)
        AnovaValue("Within Subjects", 40390.88888888888, 160.0)
        AnovaValue("Total", 42310.99444444444, 179.0)
        AnovaFactor("S/B", 3352.8777777777673, 38.0, 88.23362573099388)
        AnovaFactor("S/C", 7785.87777777778, 38.0, 204.8915204678363)
        AnovaValue("S/BC", 2906.6888888888916, 76.0)

        @testset "Repeated measures ANOVA with two within-subject factors and one among-subjects factors" begin
            observations = cat(cat([  1  38  10;  43  20   9; 15  20   6; 40  28  20;  8 11  27; 17  17   9;  30  15  19;  34  27  12;  34  24 12; 26  23 21],
                                   [  6  -5 -14;  30 -12 -10; 15 -15 -16; 30  -4 -10; 12 -2   5; 17  -6  -6;  21  -2 -20;  23  -7 -12;  20 -10 -9; 27 -15 -6],
                                   [  5   4  -2;   8   4 -13; 12   6   1; 19   0   2;  8  6  -5; 15   6 -13;  21  16   3;  28   7   2;  26  12  4; 27  14  0], dims = 3),
                               cat([  1  28  33;   7  26  23; 22  34  21; 30  32  17; 40 24  15; 15  29  13;  20  30  16;   9  24  17;  14  34 19; 15  23 29],
                                   [-19 -13  -2; -18 -16 -17; -8 -23 -19; -6 -22 -11; -6 -9 -10; -9 -18 -17; -17 -17  -4; -12 -15  -4; -11 -14 -1; -6 -15 -1],
                                   [-10  13   9;   6  19   5;  4  14   0;  3  21   4;  0 19   2;  4   7   8;   9  12  10;  -5  18   8;   7  20 12; 13  15 10], dims = 3), dims = 4)
            
            # NOT CORRECT VALUES YET
            expected = [AnovaValue( "Total", 1827.6975, 19),
                        AnovaResult(    "A",   70.3125,  1,   70.3125,  3.0706495,  0.098856175,  0.025621074),
                        AnovaResult(    "B", 1386.1125,  1, 1386.1125, 60.533556,   7.9430782e-7, 0.73663535),
                        AnovaResult("A × B",    4.9005,  1,    4.9005,  0.21401199, 0.64987001,  -0.0097253817),
                        AnovaFactor("Error",  366.372,  16,   22.89825)]
 
            @testset "Replicate First Dimension" begin
                results = anova(observations2, [within, subject])
                @test all(expected .≈ results.effects)
            end
        end
    end

    @testset "Vectors" begin
        observations = [1.9, 1.8, 1.6, 1.4, 2.1, 2.0, 1.8, 2.2, 1.1, 1.2, 1.0, 1.4,
                        2.3, 2.1, 2.0, 2.6, 2.4, 2.6, 2.7, 2.3, 2.0, 2.1, 1.9, 2.2,
                        2.9, 2.8, 3.4, 3.2, 3.6, 3.1, 3.4, 3.2, 2.9, 2.8, 3.0, 3.1,
                        1.8, 1.7, 1.4, 1.5, 2.3, 2.0, 1.9, 1.7, 1.4, 1.0, 1.3, 1.2,
                        2.4, 2.7, 2.4, 2.6, 2.0, 2.3, 2.1, 2.4, 2.4, 2.6, 2.3, 2.2,
                        3.0, 3.1, 3.0, 2.7, 3.1, 3.0, 2.8, 3.2, 3.2, 2.9, 2.8, 2.9]

        factorassignments = [[repeat([1], 36); repeat([2], 36)],
                             repeat([repeat([1], 12); repeat([2], 12); repeat([3], 12)], 2),
                             repeat([repeat([1], 4); repeat([2], 4); repeat([3], 4)], 6)]

        expected = [AnovaValue(     "Total", 30.355,       71),
                    AnovaResult(        "A",  1.8175,       2,  0.90875,      24.475062,    2.71459995e-8, 0.057358295),
                    AnovaResult(        "B", 24.655833,     2, 12.3279167,   332.02369,     4.5550981e-31, 0.8088138),
                    AnovaResult(        "C",  0.008888889,  1,  0.008888889,   0.239401496, 0.62662035,   -0.0009292123),
                    AnovaResult(    "A × B",  1.10166667,   4,  0.27541667,    7.4177057,   7.7516681e-5,  0.031361677),
                    AnovaResult(    "A × C",  0.37027778,   2,  0.18513889,    4.9862843,   0.0102990988,  0.009739973),
                    AnovaResult(    "B × C",  0.17527778,   2,  0.08763889,    2.3603491,   0.104056404,   0.0033238381),
                    AnovaResult("A × B × C",  0.220555556,  4,  0.055138889,   1.4850374,   0.21958095,    0.002370253),
                    AnovaFactor(    "Error",  2.005,       54,  0.03712963)]

        @testset "3-way ANOVA ordered" begin
            results = anova(observations, factorassignments)
            @test all(expected .≈ results.effects)
        end

        @testset "3-way ANOVA shuffled" begin
            randorder = [65,  6, 31, 70, 45,  2, 41, 59,  3, 39, 40, 47, 61, 55, 13, 68, 67, 25,
                         60, 66, 44, 33, 42, 52, 56, 15, 49, 43, 69, 10, 51, 58, 19, 22,  8, 37,
                         71, 23, 46, 24, 54, 53, 16, 18, 11, 32, 29, 35, 63, 34, 64, 62, 28, 72,
                         30, 38, 17, 21,  9,  5, 50, 12,  7, 36, 48, 27, 14, 20,  4,  1, 57, 26]
            shuffledobservations = observations[randorder]
            shuffledfactorassignments = [fa[randorder] for fa ∈ factorassignments]

            results = anova(shuffledobservations, shuffledfactorassignments)

            @test all(expected .≈ results.effects)
        end

        @testset "3-way ANOVA tolerates nonconsecutive factor levels" begin
            nonconsecutivefactorassignments = [fa .* 2 for fa ∈ factorassignments]
            results = anova(observations, nonconsecutivefactorassignments)
            @test all(expected .≈ results.effects)
        end

        @testset "3-way ANOVA tolerates noninteger factor levels" begin
            factorassignmentstrings = [repeat(["1"], 36); repeat(["2"], 36)]
            stringsinfactorassignments = [[factorassignmentstrings]; factorassignments[2:3]]
            results = anova(observations, stringsinfactorassignments)
            @test all(expected .≈ results.effects)
        end

        @testset "2-way ANOVA rejects nested factor with exclusive levels" begin
            factorlevels = factorassignments .|> unique .|> sort
            nfactorlevels = length.(factorlevels)
            nlevels = [4; nfactorlevels]
            nestedfactorassignment = (repeat(1:4, Int(length(observations) / 4)) .+
                                      sum([(factorassignments[i] .- 1) .* prod(nlevels[1:i]) for i ∈ 1:length(factorassignments)])) ./ 4 .|> ceil .|> Int
            nestedfactorassignments = [[nestedfactorassignment]; factorassignments[2:3]]

            @test_throws Exception anova(observations, nestedfactorassignments, [nested])
        end
    end

    @testset "DataFrame" begin
        observations = [1.9, 1.8, 1.6, 1.4, 2.1, 2.0, 1.8, 2.2, 1.1, 1.2, 1.0, 1.4,
                        2.3, 2.1, 2.0, 2.6, 2.4, 2.6, 2.7, 2.3, 2.0, 2.1, 1.9, 2.2,
                        2.9, 2.8, 3.4, 3.2, 3.6, 3.1, 3.4, 3.2, 2.9, 2.8, 3.0, 3.1,
                        1.8, 1.7, 1.4, 1.5, 2.3, 2.0, 1.9, 1.7, 1.4, 1.0, 1.3, 1.2,
                        2.4, 2.7, 2.4, 2.6, 2.0, 2.3, 2.1, 2.4, 2.4, 2.6, 2.3, 2.2,
                        3.0, 3.1, 3.0, 2.7, 3.1, 3.0, 2.8, 3.2, 3.2, 2.9, 2.8, 2.9]

        factorassignments = [[repeat([1], 36); repeat([2], 36)],
                             repeat([repeat([1], 12); repeat([2], 12); repeat([3], 12)], 2),
                             repeat([repeat([1], 4); repeat([2], 4); repeat([3], 4)], 6)]

        expected = [AnovaValue(                       "Total", 30.355,       71),
                    AnovaResult(                    "FactorA",  1.8175,       2,  0.90875,      24.475062,    2.71459995e-8, 0.057358295),
                    AnovaResult(                    "FactorB", 24.655833,     2, 12.3279167,   332.02369,     4.5550981e-31, 0.8088138),
                    AnovaResult(                    "FactorC",  0.008888889,  1,  0.008888889,   0.239401496, 0.62662035,   -0.0009292123),
                    AnovaResult(          "FactorA × FactorB",  1.10166667,   4,  0.27541667,    7.4177057,   7.7516681e-5,  0.031361677),
                    AnovaResult(          "FactorA × FactorC",  0.37027778,   2,  0.18513889,    4.9862843,   0.0102990988,  0.009739973),
                    AnovaResult(          "FactorB × FactorC",  0.17527778,   2,  0.08763889,    2.3603491,   0.104056404,   0.0033238381),
                    AnovaResult("FactorA × FactorB × FactorC",  0.220555556,  4,  0.055138889,   1.4850374,   0.21958095,    0.002370253),
                    AnovaFactor(                      "Error",  2.005,       54,  0.03712963)]

        @testset "3-way ANOVA shuffled" begin
            df = DataFrame([[observations]; factorassignments], [:observations; :FactorC; :FactorB; :FactorA])
            results = anova(df, :observations, [:FactorC; :FactorB; :FactorA])
            @test all(expected .≈ results.effects)
        end
    end
end
