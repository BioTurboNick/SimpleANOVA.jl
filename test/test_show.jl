@testset "Display Tests" begin
    @testset "AnovaData" begin
        data = [AnovaValue( "Total", 1827.6975, 19),
                AnovaResult(    "A",   70.3125,  1,   70.3125,  3.0706495,  0.098856175),
                AnovaResult(    "B", 1386.1125,  1, 1386.1125, 60.533556,   7.9430782e-7),
                AnovaResult("A × B",    4.9005,  1,    4.9005,  0.21401199, 0.64987001),
                AnovaFactor(    "C",  366.372,  16,   22.89825)]
        result = AnovaData(data, AnovaFactor[], AnovaFactor[], Float64[], 1)
        expectedlines = ["",
                         "Analysis of Variance Results",
                         "",
                         "Effect         SS  DF         MS          F           p",
                         "-------------------------------------------------------",
                         " Total  1827.7     19                                  ",
                         "     A    70.3125   1    70.3125   3.07065   0.0988562 ",
                         "     B  1386.11     1  1386.11    60.5336    7.94308e-7",
                         " A × B     4.9005   1     4.9005   0.214012  0.64987   ",
                         "     C   366.372   16    22.8983                       ",
                         ""]

        @test sprint(show, result) == join(expectedlines, "\n")
    end
end
