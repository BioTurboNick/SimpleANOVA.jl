@testset "Pretest Tests" begin
    @testset "Levene" begin
        @testset "Multidimensional Array" begin
            @testset "1-way ANOVA" begin
                expected = [AnovaValue(  "Total", 48.24912, 19),
                            AnovaResult("Groups",  1.69792,  3, 0.56597333, 0.19452932, 0.89858, -0.137424287),
                            AnovaFactor( "Error", 46.5512,  16, 2.90945)]

                @testset "Replicate Vectors" begin
                    observations = Array{Vector{Float64}, 1}(undef, 4)
                    observations[1] = [ 60.8,  57.0,  65.0, 58.6,  61.7]
                    observations[2] = [ 68.7,  67.7,  74.9, 66.3,  69.8]
                    observations[3] = [102.6, 102.1, 100.2, 96.5, 100.4]
                    observations[4] = [ 87.9,  84.2,  83.1, 85.7,  90.3]

                    results = levene(observations)

                    @test all(expected .≈ results.effects)
                end

                @testset "Replicate First Dimension" begin
                    observations = [60.8 68.7 102.6 87.9;
                                    57.0 67.7 102.1 84.2;
                                    65.0 74.9 100.2 83.1;
                                    58.6 66.3  96.5 85.7;
                                    61.7 69.8 100.4 90.3]

                    results = levene(observations)

                    @test all(expected .≈ results.effects)
                end
            end

            @testset "2-way ANOVA" begin
                expected = [AnovaValue(  "Total", 135.172,   19),
                            AnovaResult("Groups",  89.24992,  3,  29.749973, 10.365375, 0.00049480716, 0.58416608),
                            AnovaFactor( "Error",  45.92208, 16,   2.87013)]

                @testset "Replicate Vectors" begin
                    observations = Array{Vector{Float64}, 2}(undef, 2, 2)
                    observations[1,1] = [16.5, 18.4, 12.7, 14.0, 12.8]
                    observations[1,2] = [14.5, 11.0, 10.8, 14.3, 10.0]
                    observations[2,1] = [39.1, 26.2, 21.3, 35.8, 40.2]
                    observations[2,2] = [32.0, 23.8, 28.8, 25.0, 29.3]

                    results = levene(observations)

                    @test all(expected .≈ results.effects)
                end

                @testset "Replicate First Dimension" begin
                    observations = cat(hcat([16.5, 18.4, 12.7, 14.0, 12.8], [39.1, 26.2, 21.3, 35.8, 40.2]),
                                        hcat([14.5, 11.0, 10.8, 14.3, 10.0], [32.0, 23.8, 28.8, 25.0, 29.3]), dims = 3)

                    results = levene(observations)

                    @test all(expected .≈ results.effects)
                end
            end

            @testset "2-way ANOVA without replicates" begin
                @testset "Fixed-effects (Model I)" begin
                    @testset "Replicate Vectors" begin
                        observations = Array{Vector{Float64}, 2}(undef, 3, 4)
                        observations[1,1] = [123]
                        observations[1,2] = [138]
                        observations[1,3] = [110]
                        observations[1,4] = [151]
                        observations[2,1] = [145]
                        observations[2,2] = [165]
                        observations[2,3] = [140]
                        observations[2,4] = [167]
                        observations[3,1] = [156]
                        observations[3,2] = [176]
                        observations[3,3] = [185]
                        observations[3,4] = [175]

                        @test_throws ErrorException levene(observations)
                    end
                end
            end

            @testset "3-way ANOVA" begin
                expected = [AnovaValue(  "Total", 0.50277778,  71),
                            AnovaResult("Groups", 0.077777778, 17, 0.0045751634, 0.58131488, 0.8912173, -0.109700816),
                            AnovaFactor( "Error", 0.425,       54, 0.0078703704)]

                @testset "Replicate Vectors" begin
                    observations = Array{Vector{Float64}, 3}(undef, 2, 3, 3)
                    observations[1,1,1] = [1.9, 1.8, 1.6, 1.4]
                    observations[1,1,2] = [2.1, 2.0, 1.8, 2.2]
                    observations[1,1,3] = [1.1, 1.2, 1.0, 1.4]
                    observations[1,2,1] = [2.3, 2.1, 2.0, 2.6]
                    observations[1,2,2] = [2.4, 2.6, 2.7, 2.3]
                    observations[1,2,3] = [2.0, 2.1, 1.9, 2.2]
                    observations[1,3,1] = [2.9, 2.8, 3.4, 3.2]
                    observations[1,3,2] = [3.6, 3.1, 3.4, 3.2]
                    observations[1,3,3] = [2.9, 2.8, 3.0, 3.1]
                    observations[2,1,1] = [1.8, 1.7, 1.4, 1.5]
                    observations[2,1,2] = [2.3, 2.0, 1.9, 1.7]
                    observations[2,1,3] = [1.4, 1.0, 1.3, 1.2]
                    observations[2,2,1] = [2.4, 2.7, 2.4, 2.6]
                    observations[2,2,2] = [2.0, 2.3, 2.1, 2.4]
                    observations[2,2,3] = [2.4, 2.6, 2.3, 2.2]
                    observations[2,3,1] = [3.0, 3.1, 3.0, 2.7]
                    observations[2,3,2] = [3.1, 3.0, 2.8, 3.2]
                    observations[2,3,3] = [3.2, 2.9, 2.8, 2.9]

                    results = levene(observations)

                    @test all(expected .≈ results.effects)
                end

                @testset "Replicate First Dimension" begin
                    observations = cat(cat(hcat([1.9, 1.8, 1.6, 1.4], [1.8, 1.7, 1.4, 1.5]),
                                            hcat([2.3, 2.1, 2.0, 2.6], [2.4, 2.7, 2.4, 2.6]),
                                            hcat([2.9, 2.8, 3.4, 3.2], [3.0, 3.1, 3.0, 2.7]), dims = 3),
                                        cat(hcat([2.1, 2.0, 1.8, 2.2], [2.3, 2.0, 1.9, 1.7]),
                                            hcat([2.4, 2.6, 2.7, 2.3], [2.0, 2.3, 2.1, 2.4]),
                                            hcat([3.6, 3.1, 3.4, 3.2], [3.1, 3.0, 2.8, 3.2]), dims = 3),
                                        cat(hcat([1.1, 1.2, 1.0, 1.4], [1.4, 1.0, 1.3, 1.2]),
                                            hcat([2.0, 2.1, 1.9, 2.2], [2.4, 2.6, 2.3, 2.2]),
                                            hcat([2.9, 2.8, 3.0, 3.1], [3.2, 2.9, 2.8, 2.9]), dims = 3), dims = 4)

                    results = levene(observations)

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

            expected = [AnovaValue(  "Total", 0.50277778,  71),
                        AnovaResult("Groups", 0.077777778, 17, 0.0045751634, 0.58131488, 0.8912173, -0.109700816),
                        AnovaFactor( "Error", 0.425,       54, 0.0078703704)]

            @testset "3-way ANOVA ordered" begin
                results = levene(observations, factorassignments)
                @test all(expected .≈ results.effects)
            end

            @testset "3-way ANOVA shuffled" begin
                randorder = [65,  6, 31, 70, 45,  2, 41, 59,  3, 39, 40, 47, 61, 55, 13, 68, 67, 25,
                             60, 66, 44, 33, 42, 52, 56, 15, 49, 43, 69, 10, 51, 58, 19, 22,  8, 37,
                             71, 23, 46, 24, 54, 53, 16, 18, 11, 32, 29, 35, 63, 34, 64, 62, 28, 72,
                             30, 38, 17, 21,  9,  5, 50, 12,  7, 36, 48, 27, 14, 20,  4,  1, 57, 26]
                shuffledobservations = observations[randorder]
                shuffledfactorassignments = [fa[randorder] for fa ∈ factorassignments]

                results = levene(shuffledobservations, shuffledfactorassignments)

                @test all(expected .≈ results.effects)
            end

            @testset "3-way ANOVA tolerates nonconsecutive factor levels" begin
                nonconsecutivefactorassignments = [fa .* 2 for fa ∈ factorassignments]

                results = levene(observations, nonconsecutivefactorassignments)

                @test all(expected .≈ results.effects)
            end

            @testset "3-way ANOVA tolerates noninteger factor levels" begin
                factorassignmentstrings = [repeat(["1"], 36); repeat(["2"], 36)]
                stringsinfactorassignments = [[factorassignmentstrings]; factorassignments[2:3]]

                results = levene(observations, stringsinfactorassignments)

                @test all(expected .≈ results.effects)
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

            expected = [AnovaValue(  "Total", 0.50277778,  71),
                        AnovaResult("Groups", 0.077777778, 17, 0.0045751634, 0.58131488, 0.8912173, -0.109700816),
                        AnovaFactor( "Error", 0.425,       54, 0.0078703704)]

            @testset "3-way ANOVA shuffled" begin
                df = DataFrame([[observations]; factorassignments], [:observations; :FactorC; :FactorB; :FactorA])

                results = levene(df, :observations, [:FactorC; :FactorB; :FactorA])

                @test all(expected .≈ results.effects)
            end
        end
    end
end
