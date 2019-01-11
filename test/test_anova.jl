@testset "1-way ANOVA" begin
    expected = [AnovaValue( "Total", 4822.4575, 19),
                AnovaResult(    "A", 4684.9975,  3, 1561.66583, 181.773995, 1.43627554e-12),
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
                    AnovaResult(    "A",   70.3125,  1,   70.3215,  3.0706495,  0.098856175),
                    AnovaResult(    "B", 1386.1125,  1, 1386.1125, 60.533556,   7.9430782e-7),
                    AnovaResult("A × B",    4.9005,  1,    4.9005,  0.21401199, 0.64987001),
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
                    AnovaResult(    "A",   70.3125,  1,   70.3215,  14.3480257,  0.16431864),
                    AnovaResult(    "B", 1386.1125,  1, 1386.1125, 282.85124,    0.037808553),
                    AnovaResult("A × B",    4.9005,  1,    4.9005,   0.21401199, 0.64987001),
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
                    AnovaResult(    "A",   70.3125,  1,   70.3215, 14.3480257,  0.164318644),
                    AnovaResult(    "B", 1386.1125,  1, 1386.1125, 60.533556,   7.9430782e-7),
                    AnovaResult("A × B",    4.9005,  1,    4.9005,  0.21401199, 0.64987001),
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

    expected = [AnovaValue(     "Total", 5594.9167,  11),
                AnovaResult(        "A", 1116.91667,  3,  372.30556,  2.6316513, 0.144561215),
                AnovaResult(        "B", 3629.1667,   2, 1814.58333, 12.8264284, 0.0068110648),
                AnovaFactor("Remainder",  848.83333,  6,  141.472222)]

    @testset "Fixed-effects (Model I)" begin
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
                    AnovaResult(        "A",  1.8175,       2,  0.90875,      24.475062,    2.71459995e-8),
                    AnovaResult(        "B", 24.655833,     2, 12.3279167,   332.02369,     4.5550981e-31),
                    AnovaResult(        "C",  0.008888889,  1,  0.008888889,   0.239401496, 0.62662035),
                    AnovaResult(    "A × B",  1.10166667,   4,  0.27541667,    7.4177057,   7.7516681e-5),
                    AnovaResult(    "A × C",  0.37027778,   2,  0.18513889,    4.9862843,   0.0102990988),
                    AnovaResult(    "B × C",  0.17527778,   2,  0.08763889,    2.3603491,   0.104056404),
                    AnovaResult("A × B × C",  0.220555556,  4,  0.055138889,   1.4850374,   0.21958095),
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
                    AnovaResult(        "A",  1.8175,       2,  0.90875,      4.1754946,   0.179725175),
                    AnovaResult(        "B", 24.655833,     2, 12.3279167,   40.036536,    0.0022084069),
                    AnovaResult(        "C",  0.008888889,  1,  0.008888889,  0.021925317, 0.88874495),
                    AnovaResult(    "A × B",  1.10166667,   4,  0.27541667,   4.9949622,   0.074190827),
                    AnovaResult(    "A × C",  0.37027778,   2,  0.18513889,   3.3576826,   0.139349696),
                    AnovaResult(    "B × C",  0.17527778,   2,  0.08763889,   1.58942065,  0.31046402),
                    AnovaResult("A × B × C",  0.220555556,  4,  0.055138889,  1.4850374,   0.21958095),
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
                    AnovaResult(        "A",  1.8175,       2,  0.90875,      4.9084771,   0.16924835),
                    AnovaResult(        "B", 24.655833,     2, 12.3279167,   44.760968,    0.001829334),
                    AnovaResult(        "C",  0.008888889,  1,  0.008888889,  0.239401496, 0.62662035),
                    AnovaResult(    "A × B",  1.10166667,   4,  0.27541667,   4.9949622,   0.074190827),
                    AnovaResult(    "A × C",  0.37027778,   2,  0.18513889,   4.9862843,   0.0102990988),
                    AnovaResult(    "B × C",  0.17527778,   2,  0.08763889,   2.3603491,   0.104056404),
                    AnovaResult("A × B × C",  0.220555556,  4,  0.055138889,  1.4850374,   0.21958095),
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
                    AnovaResult(        "A",  1.8175,       2,  0.90875,      4.1754946,    0.179725175),
                    AnovaResult(        "B", 24.655833,     2, 12.3279167,   44.760968,    0.001829334),
                    AnovaResult(        "C",  0.008888889,  1,  0.008888889,  0.032274332, 0.8661604),
                    AnovaResult(    "A × B",  1.10166667,   4,  0.27541667,   4.9949622,   0.074190827),
                    AnovaResult(    "A × C",  0.37027778,   2,  0.18513889,   3.3576826,   0.139349696),
                    AnovaResult(    "B × C",  0.17527778,   2,  0.08763889,   2.3603491,   0.104056404),
                    AnovaResult("A × B × C",  0.220555556,  4,  0.055138889,  1.4850374,   0.21958095),
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