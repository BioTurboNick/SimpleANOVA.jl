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
                            AnovaResult(    "A",   70.3125,  1,   70.3125,  3.0706495,  0.098856175,  0.093819147),
                            AnovaResult(    "B", 1386.1125,  1, 1386.1125, 60.533556,   7.9430782e-7, 0.74853381),
                            AnovaResult("A × B",    4.9005,  1,    4.9005,  0.21401199, 0.64987001,  -0.040907022),
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
                            AnovaResult(    "A",   70.3125,  1,   70.3125,  14.3480257,  0.16431864,  0.034672057),
                            AnovaResult(    "B", 1386.1125,  1, 1386.1125, 282.85124,    0.037808553, 0.73212043),
                            AnovaResult("A × B",    4.9005,  1,    4.9005,   0.21401199, 0.64987001, -0.0095398248),
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
                            AnovaResult(    "A",   70.3125,  1,   70.3125, 14.3480257,  0.16431864,   0.035006009),
                            AnovaResult(    "B", 1386.1125,  1, 1386.1125, 60.533556,   7.9430782e-7, 0.7560050),
                            AnovaResult("A × B",    4.9005,  1,    4.9005,  0.21401199, 0.64987001,  -0.0099811084),
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
                            AnovaResult(        "A", 1116.91667,  3,  372.30556,  2.6316513, 0.144561215,  0.28972875),
                            AnovaResult(        "B", 3629.1667,   2, 1814.58333, 12.8264284, 0.0068110648, 0.66342108),
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
                            AnovaResult(        "A", 1116.91667,  3,  372.30556,  2.6316513, 0.144561215,  0.120720546),
                            AnovaResult(        "B", 3629.1667,   2, 1814.58333, 12.8264284, 0.0068110648, 0.58333253),
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
                            AnovaResult(        "A", 1116.91667,  3,  372.30556,  2.6316513, 0.144561215,  0.120720546),
                            AnovaResult(        "B", 3629.1667,   2, 1814.58333, 12.8264284, 0.0068110648, 0.66342108),
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
                        AnovaResult(    "A", 61.166667,  2, 30.583333, 61.166667,   0.0037032412, 0.80044346),
                        AnovaResult(    "B",  1.5,       3,  0.5,       0.33333333, 0.80220227,  -0.2),
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
                        AnovaResult(    "A",   7.5108875,  2, 3.7554438, 2.0223607, 0.27789847,   0.036110684),
                        AnovaResult(    "B",   5.5708812,  3, 1.8569604, 0.99348408, 0.45717102, -0.00036050442),
                        AnovaResult(    "C",  11.2148375,  6, 1.8691396, 0.87059412, 0.52587983, -0.0164476185),
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
               expected = [AnovaValue(     "Total", 1827.6975, 19),
                           AnovaResult(        "A",   70.3125,  1,   70.3125,  3.0706495,  0.098856175,  0.093819147),
                           AnovaResult(        "B", 1386.1125,  1, 1386.1125, 60.533556,   7.9430782e-7, 0.74853381),
                           AnovaResult(    "A × B",    4.9005,  1,    4.9005,  0.21401199, 0.64987001,  -0.040907022),
                           AnovaFactor("Remainder",  366.372,  16,   22.89825)]

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
               expected = [AnovaValue(     "Total", 1827.6975, 19),
                           AnovaResult(        "A",   70.3125,  1,   70.3125,  14.3480257,  0.16431864,  0.034672057),
                           AnovaResult(        "B", 1386.1125,  1, 1386.1125, 282.85124,    0.037808553, 0.73212043),
                           AnovaResult(    "A × B",    4.9005,  1,    4.9005,   0.21401199, 0.64987001, -0.0095398248),
                           AnovaFactor("Remainder",  366.372,  16,   22.89825)]

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
               expected = [AnovaValue(     "Total", 1827.6975, 19),
                           AnovaResult(         "A",   70.3125,  1,   70.3125, 14.3480257,  0.16431864,   0.035006009),
                           AnovaResult(         "B", 1386.1125,  1, 1386.1125, 60.533556,   7.9430782e-7, 0.7560050),
                           AnovaResult(     "A × B",    4.9005,  1,    4.9005,  0.21401199, 0.64987001,  -0.0099811084),
                           AnovaFactor("Remainder",  366.372,  16,   22.89825)]

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
                           AnovaResult(    "A",  144.4,    1,  144.4,      3.0540486,   0.099697096,   0.113261094),
                           AnovaResult(    "B", 2692.881,  1, 2692.881,   56.954221,    1.17181416e-6, 0.77675692),
                           AnovaResult("A × B",   11.236,  1,   11.236,    0.237640515, 0.63252737,   -0.04976526),
                           AnovaResult(    "C",  756.504, 16,   47.2815, 294.588785,    3.4334068e-20, 0.9915566),
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
                            AnovaResult(    "A",  144.4,    1,   144.4,     12.8515486,   0.17318116,    0.037626308),
                            AnovaResult(    "B", 2692.881,  1,  2692.881,  239.66545,     0.04106525,    0.7577153),
                            AnovaResult("A × B",   11.236,  1,    11.236,    0.237640515, 0.63252737,   -0.010184878),
                            AnovaResult(    "C",  756.504, 16,    47.2815, 294.588785,    3.4334068e-20, 0.21302926),
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
                            AnovaResult(    "A",  144.4,    1,  144.4,     12.8515486,   0.17318116,    0.03801347),
                            AnovaResult(    "B", 2692.881,  1, 2692.881,   56.954221,    1.17181416e-6, 0.78506533),
                            AnovaResult("A × B",   11.236,  1,   11.236,    0.237640515, 0.63252737,   -0.0106962797),
                            AnovaResult(    "C",  756.504, 16,   47.2815, 294.588785,    3.4334068e-20, 0.223725856),
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
                            AnovaResult(        "A",  1.8175,       2,  0.90875,      24.475062,    2.71459995e-8, 0.39470429),
                            AnovaResult(        "B", 24.655833,     2, 12.3279167,   332.02369,     4.5550981e-31, 0.9019137),
                            AnovaResult(        "C",  0.008888889,  1,  0.008888889,   0.239401496, 0.62662035,   -0.010676655),
                            AnovaResult(    "A × B",  1.10166667,   4,  0.27541667,    7.4177057,   7.7516681e-5,  0.262830006),
                            AnovaResult(    "A × C",  0.37027778,   2,  0.18513889,    4.9862843,   0.0102990988,  0.09969129),
                            AnovaResult(    "B × C",  0.17527778,   2,  0.08763889,    2.3603491,   0.104056404,   0.036411574),
                            AnovaResult("A × B × C",  0.220555556,  4,  0.055138889,   1.4850374,   0.21958095,    0.02623946),
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
                            AnovaResult(        "A",  1.8175,       2,  0.90875,      2.2415211,   0.211958974,  0.034966576),
                            AnovaResult(        "B", 24.655833,     2, 12.3279167,   40.036536,    0.0022084069, 0.83502961),
                            AnovaResult(        "C",  0.008888889,  1,  0.008888889,  0.040842374, 0.85716044,  -0.007250933),
                            AnovaResult(    "A × B",  1.10166667,   4,  0.27541667,   4.9949622,   0.074190827,  0.030605402),
                            AnovaResult(    "A × C",  0.37027778,   2,  0.18513889,   3.3576826,   0.139349696,  0.0090311023),
                            AnovaResult(    "B × C",  0.17527778,   2,  0.08763889,   1.58942065,  0.31046402,   0.0022577756),
                            AnovaResult("A × B × C",  0.220555556,  4,  0.055138889,  1.4850374,   0.21958095,   0.0025022071),
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
                            AnovaResult(        "A",  1.8175,       2,  0.90875,       4.9084771,   0.16924835,   0.31727666),
                            AnovaResult(        "B", 24.655833,     2, 12.3279167,   140.667195,    0.0070587972, 0.8871463),
                            AnovaResult(        "C",  0.008888889,  1,  0.008888889,   0.239401496, 0.62662035,  -0.0090684744),
                            AnovaResult(    "A × B",  1.10166667,   4,  0.27541667,    4.9949622,   0.074190827,  0.220538135),
                            AnovaResult(    "A × C",  0.37027778,   2,  0.18513889,    4.9862843,   0.0102990988, 0.095055451),
                            AnovaResult(    "B × C",  0.17527778,   2,  0.08763889,    2.3603491,   0.104056404,  0.032438379),
                            AnovaResult("A × B × C",  0.220555556,  4,  0.055138889,   1.4850374,   0.21958095,   0.023132043),
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
                            AnovaResult(        "A",  1.8175,       2,  0.90875,       2.2415211,   0.211958974, 0.034245054),
                            AnovaResult(        "B", 24.655833,     2, 12.3279167,   140.667195,    0.0070587972, 0.86231611),
                            AnovaResult(        "C",  0.008888889,  1,  0.008888889,   0.101426307, 0.78030599, -0.0027739319),
                            AnovaResult(    "A × B",  1.10166667,   4,  0.27541667,    4.9949622,   0.074190827, 0.0310367264),
                            AnovaResult(    "A × C",  0.37027778,   2,  0.18513889,    3.3576826,   0.139349696, 0.0091583783),
                            AnovaResult(    "B × C",  0.17527778,   2,  0.08763889,    2.3603491,   0.104056404, 0.00355833),
                            AnovaResult("A × B × C",  0.220555556,  4,  0.055138889,   1.4850374,   0.21958095,  0.0025374709),
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

#=
Repeated measures with 1 within 1 between

observations2 = cat(hcat([10, 15, 20], [11, 16, 21], [12, 15, 19]),
                    hcat([8, 13, 19], [7, 14, 18], [10, 13, 16]),
                    hcat([14, 20, 24], [13, 18, 23], [15, 18, 25]), dims = 3)
=#

        @testset "Repeated measures ANOVA with no among-subject factors" begin
            observations1 = Array{Vector{Float64}, 2}(undef, 4, 8)
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
            
            expected = [AnovaValue(     "Total", 253.875, 31),
                        AnovaResult(        "A",  83.125,  3, 27.708333,   3.793806, 0.025570297, 0.23878518),
                        AnovaFactor(        "S",  17.375,  7,  2.48214286),
                        AnovaFactor("Remainder", 153.375, 21,  7.3035714)]
 
            @testset "Replicate Vectors" begin
                results = anova(observations1, [fixed, subject], factornames = ["A", "S"], hasreplicates = false)
                @test all(expected .≈ results.effects)
            end
            
            @testset "Replicate First Dimension" begin
                results = anova(observations2, [fixed, subject], factornames = ["A", "S"], hasreplicates = false)
                @test all(expected .≈ results.effects)
            end
        end

        @testset "Repeated measures ANOVA with two within-subject factors and no among-subjects factors" begin
            observations = cat([  1  38  10;  43  20   9; 15  20   6; 40  28  20;  8 11  27; 17  17   9;  30  15  19;  34  27  12;  34  24 12; 26  23 21; 1  28  33;   7  26  23; 22  34  21; 30  32  17; 40 24  15; 15  29  13;  20  30  16;   9  24  17;  14  34 19; 15  23 29],
                            [  6  -5 -14;  30 -12 -10; 15 -15 -16; 30  -4 -10; 12 -2   5; 17  -6  -6;  21  -2 -20;  23  -7 -12;  20 -10 -9; 27 -15 -6; -19 -13  -2; -18 -16 -17; -8 -23 -19; -6 -22 -11; -6 -9 -10; -9 -18 -17; -17 -17  -4; -12 -15  -4; -11 -14 -1; -6 -15 -1],
                            [  5   4  -2;   8   4 -13; 12   6   1; 19   0   2;  8  6  -5; 15   6 -13;  21  16   3;  28   7   2;  26  12  4; 27  14  0; -10  13   9;   6  19   5;  4  14   0;  3  21   4;  0 19   2;  4   7   8;   9  12  10;  -5  18   8;   7  20 12; 13  15 10], dims = 3)
                                
            observations = permutedims(observations, (2,3,1))

            expected = [AnovaValue(     "Total", 42310.994,   179),
                        AnovaResult(        "A", 21628.678,     2, 10814.339,    122.564825, 2.6801966e-17, 0.56072406),
                        AnovaResult(        "B",  2092.34444,   2,  1046.17222,    5.105981, 0.0108629307,  0.09100631),
                        AnovaResult(    "A × B",  2624.4222,    4,   656.10556,  17.1549224, 4.5890403e-10, 0.12820473),
                        AnovaFactor(    "A × S",  3352.8778,   38,    88.233626),
                        AnovaFactor(    "B × S",  7785.8778,   38,   204.89152),
                        AnovaFactor(        "S",  1920.10556,  19,   101.058187),
                        AnovaFactor("Remainder",  2906.6889,   76,    38.245906)]

            @testset "Replicate First Dimension" begin
                results = anova(observations, [fixed, fixed, subject], factornames = ["B", "A", "S"], hasreplicates = false)
                @test all(expected .≈ results.effects)
            end
        end

        @testset "Repeated measures ANOVA with two within-subject factors and one among-subjects factors" begin
            observations = cat(cat([  1  38  10;  43  20   9; 15  20   6; 40  28  20;  8 11  27],
                                   [  6  -5 -14;  30 -12 -10; 15 -15 -16; 30  -4 -10; 12 -2   5],
                                   [  5   4  -2;   8   4 -13; 12   6   1; 19   0   2;  8  6  -5], dims = 3),
                               cat([  1  28  33;   7  26  23; 22  34  21; 30  32  17; 40 24  15],
                                   [-19 -13  -2; -18 -16 -17; -8 -23 -19; -6 -22 -11; -6 -9 -10],
                                   [-10  13   9;   6  19   5;  4  14   0;  3  21   4;  0 19   2], dims = 3), dims = 4)
            
            observations = permutedims(observations, (2,3,1,4))

            expected = [AnovaValue(     "Total", 23021.6,      89),
                        AnovaResult(        "A",   106.71111,   1,  106.71111,   1.3357441,  0.28113388,     0.0039775976),
                        AnovaResult(        "B", 11800.8667,    2, 5900.4333,  132.826163,   1.08455885e-10, 0.63553779),
                        AnovaResult(        "C",   864.26667,   2,  432.13333,   2.6831321,  0.098885156,    0.074690905),
                        AnovaResult(    "A × B",  1554.8222,    2,  777.4111,   17.5005,     9.3827532e-5,   0.17916062),
                        AnovaResult(    "A × C",  1504.6222,    2,  752.3111,    4.671128,   0.0252462965,   0.14970369),
                        AnovaResult(    "B × C",  1483.86667,   4,  370.96667,   8.206975,   0.000113071502, 0.16248534),
                        AnovaResult("A × B × C",   333.244445,  4,   83.31111,   1.84310954, 0.14485915,     0.0221924946),
                        AnovaFactor(     "S[A]",   639.11111,   8,   79.88889),
                        AnovaFactor( "B × S[A]",   710.75556,  16,   44.422222),
                        AnovaFactor( "C × S[A]",  2576.8889,   16,  161.055556),
                        AnovaFactor("Remainder",  1446.44444,  32,   45.201389)]
 
            @testset "Replicate First Dimension" begin
                results = anova(observations, [fixed, fixed, subject], factornames = ["C", "B", "S", "A"], hasreplicates = false)
                @test all(expected .≈ results.effects)
            end
        end

        @testset "Subject or block as first factor throws error" begin
            observations2 = cat(hcat([16.5, 18.4, 12.7, 14.0, 12.8], [39.1, 26.2, 21.3, 35.8, 40.2]),
                                hcat([14.5, 11.0, 10.8, 14.3, 10.0], [32.0, 23.8, 28.8, 25.0, 29.3]), dims = 3)
            
            @test_throws ErrorException anova(observations2, [subject], factornames = ["S", "A"]) # should change to ArgumentError in code
            @test_throws ErrorException anova(observations2, [block], factornames = ["S", "A"])
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
                    AnovaResult(        "A",  1.8175,       2,  0.90875,      24.475062,    2.71459995e-8, 0.39470429),
                    AnovaResult(        "B", 24.655833,     2, 12.3279167,   332.02369,     4.5550981e-31, 0.9019137),
                    AnovaResult(        "C",  0.008888889,  1,  0.008888889,   0.239401496, 0.62662035,   -0.010676655),
                    AnovaResult(    "A × B",  1.10166667,   4,  0.27541667,    7.4177057,   7.7516681e-5,  0.262830006),
                    AnovaResult(    "A × C",  0.37027778,   2,  0.18513889,    4.9862843,   0.0102990988,  0.09969129),
                    AnovaResult(    "B × C",  0.17527778,   2,  0.08763889,    2.3603491,   0.104056404,   0.036411574),
                    AnovaResult("A × B × C",  0.220555556,  4,  0.055138889,   1.4850374,   0.21958095,    0.02623946),
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

        @testset "Replicates" begin
            factorassignments = [[repeat([1], 36); repeat([2], 36)],
                                 repeat([repeat([1], 12); repeat([2], 12); repeat([3], 12)], 2),
                                 repeat([repeat([1], 4); repeat([2], 4); repeat([3], 4)], 6)]

            expected = [AnovaValue(     "Total", 30.355,       71),
                        AnovaResult(        "A",  1.8175,       2,  0.90875,      24.475062,    2.71459995e-8, 0.39470429),
                        AnovaResult(        "B", 24.655833,     2, 12.3279167,   332.02369,     4.5550981e-31, 0.9019137),
                        AnovaResult(        "C",  0.008888889,  1,  0.008888889,   0.239401496, 0.62662035,   -0.010676655),
                        AnovaResult(    "A × B",  1.10166667,   4,  0.27541667,    7.4177057,   7.7516681e-5,  0.262830006),
                        AnovaResult(    "A × C",  0.37027778,   2,  0.18513889,    4.9862843,   0.0102990988,  0.09969129),
                        AnovaResult(    "B × C",  0.17527778,   2,  0.08763889,    2.3603491,   0.104056404,   0.036411574),
                        AnovaResult("A × B × C",  0.220555556,  4,  0.055138889,   1.4850374,   0.21958095,    0.02623946),
                        AnovaFactor(    "Error",  2.005,       54,  0.03712963)]

            @testset "3-way ANOVA shuffled" begin
                df = DataFrame([[observations]; factorassignments], [:observations; :C; :B; :A])
                results = anova(df, :observations, [:C; :B; :A])
                @test all(expected .≈ results.effects)
            end
        end

        @testset "No replicates" begin
            observations = cat(hcat([16.5, 18.4, 12.7, 14.0, 12.8], [39.1, 26.2, 21.3, 35.8, 40.2]),
                               hcat([14.5, 11.0, 10.8, 14.3, 10.0], [32.0, 23.8, 28.8, 25.0, 29.3]), dims = 3) |> vec

            df = DataFrame(
                 C = repeat(1:5, 4),
                 B = repeat(1:2, inner=5, outer=2),
                 A = repeat(1:2, inner = 10),
                 observations = observations
             )

             expected = [AnovaValue(     "Total", 1827.6975, 19),
                         AnovaResult(        "A",   70.3125,  1,   70.3125,  3.0706495,  0.098856175,  0.093819147),
                         AnovaResult(        "B", 1386.1125,  1, 1386.1125, 60.533556,   7.9430782e-7, 0.74853381),
                         AnovaResult(    "A × B",    4.9005,  1,    4.9005,  0.21401199, 0.64987001,  -0.040907022),
                         AnovaFactor("Remainder",  366.372,  16,   22.89825)]

            @testset "3-way ANOVA shuffled" begin
                results = anova(df, :observations, [:C; :B; :A])
                @test all(expected .≈ results.effects)
            end
        end
    end

    @testset "CategoricalArrays" begin
        df = DataFrame(
                 C = categorical(repeat(1:2; inner=6)),
                 B = categorical(repeat(1:3, inner=2, outer=2)),
                 A = categorical(repeat(1:2, 6)),
                 observations = round.(rand(Float64, 12), digits=1)
              )
        anova(df, :observations, [:C; :B; :A]) # just looking to not fail
    end
end
