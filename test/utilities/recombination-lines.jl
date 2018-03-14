@testset "recombination-lines.jl" begin
    # http://www.cv.nrao.edu/course/astr534/Recombination.html
    @test abs(BPJSpec.hydrogen(109, 1) - 5.0089u"GHz") < 0.0001u"GHz"
end

