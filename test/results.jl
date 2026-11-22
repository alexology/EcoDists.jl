using DataFrames: DataFrame, names, nrow, ncol

@testset "Readable results and optional DataFrames" begin
    ids = ["a0","a1","b0"]
    model = HierarchicalBeta(ids,["A","A","B"],["AAAA","AAAT","CCCC"])
    a,b = community(model,["a0","b0"]),community(model,["a1","b0"])
    result = hierarchical_components(model,a,b)
    @test result isa HierarchicalResult
    @test result.species isa SpeciesContributions
    @test result.species.species == ["A","B"]
    @test result.species.contribution == [0.125,0.0]
    @test result.species[1] == (species="A",contribution=0.125)
    @test sum(result.species.contribution) == result.molecular
    @test Matrix(result.species) == Any["A" 0.125;"B" 0.0]
    @test Matrix(result) == [0.125 0.0 0.125 0.125 1.0 2.0 2.0]
    @test occursin("species",sprint(show,MIME"text/plain"(),result.species))
    @test occursin("contribution",sprint(show,MIME"text/plain"(),result.species))
    @test occursin("taxonomic",sprint(show,MIME"text/plain"(),result))
    @test occursin("Shared-species contributions",sprint(show,MIME"text/plain"(),result))
    @test occursin("more rows",sprint(show,MIME"text/plain"(),result.species;context=(:limit=>true,:displaysize=>(5,80))))

    @test Base.get_extension(EcoDists,:EcoDistsDataFramesExt) !== nothing
    summary = DataFrame(result)
    details = DataFrame(result.species)
    @test names(summary) == ["total","taxonomic","molecular","shared_mean","molecular_fraction","nshared","nunion"]
    @test nrow(summary) == 1 && ncol(summary) == 7
    @test summary.total == [result.total]
    @test summary.nshared == [2]
    @test eltype(summary.nshared) == Int
    @test names(details) == ["species","contribution"]
    @test details.species == ["A","B"]
    @test details.contribution == [0.125,0.0]
    details.contribution[1] = 99
    details.species[1] = "changed"
    @test result.species.contribution[1] == 0.125
    @test result.species.species[1] == "A"

    disjoint = hierarchical_components(model,community(model,["a0"]),community(model,["b0"]))
    @test isempty(disjoint.species)
    @test size(Matrix(disjoint.species)) == (0,2)
    frame = DataFrame(disjoint.species)
    @test names(frame) == ["species","contribution"]
    @test eltype(frame.species) == String
    @test eltype(frame.contribution) == Float64
    @test size(frame) == (0,2)
    empty = hierarchical_components(model,community(model,String[]),community(model,String[]))
    @test isnan(DataFrame(empty).shared_mean[1])
    @test isnan(DataFrame(empty).molecular_fraction[1])

    # Rows follow the fitted mapping, independently of input sample order.
    reverse_model = HierarchicalBeta(["b","a"],["B","A"];local_method=:asv_jaccard)
    same = community(reverse_model,["a","b"])
    @test species_contributions(reverse_model,same,same).species == ["B","A"]
end
