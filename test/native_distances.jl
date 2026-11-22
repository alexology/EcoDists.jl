using Statistics: cor, cov

@testset "Julia-native mixed-trait API" begin
    blocks = [TraitBlock([0.0,1.0,4.0],:quantitative;names=["height"]),
              TraitBlock(["annual","perennial","annual"],:nominal;names=["life_cycle"]),
              TraitBlock([12,1,6],:circular;periods=12,names=["month"]),
              TraitBlock([3 1 0;0 2 2;1 0 3],:fuzzy;col_blocks=[3],names=["diet"])]
    for (name,code) in pairs(EcoDists._TRAIT_KINDS)
        kwargs = name == :circular ? (;periods=2) : (;)
        values = name in (:dichotomous,:multichoice) ? [1,0] : [1,2]
        a,b = TraitBlock(values,name;kwargs...),TraitBlock(values,code;kwargs...)
        @test isequal(a.data,b.data)
        @test a.kind == b.kind
    end
    d = Pavoine(blocks;weights=[2,1,1,3])
    x = traitdata(d)
    parts = trait_distances(d)
    @test length(parts) == 4
    for i in 1:4
        metric = traitmetric(d,i)
        @test metric isa TraitDistance
        @test metric isa SemiMetric
        @test pairwise(metric,x;dims=1) == parts[i]
        @test pairwise(traitmetric(d,d.names[i]),x;dims=1) == parts[i]
        @test pairwise(metric,permutedims(x);dims=2) == parts[i]
        @test trait_distances(d,x[1:2,:])[i] == parts[i][1:2,1:2]
        buffer = zeros(3,3)
        @test pairwise!(metric,buffer,x;dims=1) === buffer
        @test buffer == parts[i]
        @test colwise(metric,permutedims(x[1:2,:]),permutedims(x[2:3,:])) ==
            [parts[i][1,2],parts[i][2,3]]
    end
    @test evaluate(traitmetric(d,"height"),x[1,:],x[2,:]) == 0.25
    @test evaluate(traitmetric(d,"month"),x[1,:],x[2,:]) ≈ sqrt(1/6)
    @test parts[4][1,2] ≈ sqrt(1-1/sqrt(20))
    @test_throws ArgumentError traitmetric(d,"unknown")
    @test_throws BoundsError traitmetric(d,5)
    @test_throws DimensionMismatch evaluate(traitmetric(d,1),[1.0],[2.0])
    weights = [2,1,1,3]
    total = pairwise(d,x;dims=1)
    for i in 1:3, j in 1:3
        @test total[i,j] ≈ sqrt(sum(weights[k]*parts[k][i,j]^2 for k in 1:4)/sum(weights))
    end
    lower(m) = [m[2,1],m[3,1],m[3,2]]
    for squared in (true,false)
        correlations = trait_correlations(d;squared=squared)
        values = hcat(lower.(parts)...)
        global_values = lower(total)
        if squared
            values = values.^2
            global_values = global_values.^2
        end
        @test correlations.covariance ≈ cov(values)
        @test correlations.correlation ≈ cor(values)
        @test correlations.global_correlation ≈ vec(cor(values,global_values))
        @test correlations.names == d.names
    end
    @test isequal(trait_distances(blocks),trait_distances(Pavoine(blocks)))
    @test trait_correlations(blocks).names == d.names

    # Hellinger remains resolvable near identical probabilities with no
    # numerical threshold, instead of cancelling in 1 - affinity.
    delta = 1e-12
    fuzzy = TraitBlock([1-delta delta;1-2delta 2delta],:F)
    hellinger = Pavoine(fuzzy;fuzzy=:hellinger,tol=0)
    value = evaluate(hellinger,traitdata(hellinger)[1,:],traitdata(hellinger)[2,:])
    @test isfinite(value) && 0 < value < 1e-6
    @test isapprox(value,(sqrt(2)-1)*sqrt(delta/2);rtol=1e-8)
end
