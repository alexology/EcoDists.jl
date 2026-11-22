using Test
using EcoDists
using Distances: evaluate, pairwise, pairwise!, colwise, SemiMetric
using DelimitedFiles: readdlm

const FIXTURES = joinpath(@__DIR__, "fixtures")
readfixture(name) = readdlm(joinpath(FIXTURES, name * ".tsv"), '\t', Float64)

function input(name; nominal=false)
    a = readdlm(joinpath(FIXTURES,name * ".tsv"), '\t', nominal ? Any : Float64)
    return map(x -> x == "NaN" ? missing : x, a)
end

function mixed_blocks(prefix="input")
    return [TraitBlock(input(prefix*"_Q"), :Q),
            TraitBlock(input(prefix*"_O"), :O),
            TraitBlock(input(prefix*"_N"; nominal=true), :N),
            TraitBlock(input(prefix*"_D"), :D),
            TraitBlock(input(prefix*"_F"),:fuzzy;col_blocks=[3,2]),
            TraitBlock(input(prefix*"_B"),:multichoice;col_blocks=[3,2]),
            TraitBlock(input(prefix*"_C"),:circular;periods=[12,5],origin=[1,1])]
end

function compare_values(actual, expected; atol=2e-12)
    @test size(actual) == size(expected)
    for (x,y) in zip(actual,expected)
        @test (isnan(x) && isnan(y)) || isapprox(x,y; rtol=2e-11, atol=atol)
    end
end

function reference_case(name,d)
    @testset "$name" begin
        compare_values(pairwise(d),readfixture(name*"_global"))
        n = size(traitdata(d),1)
        lower(m) = [m[i,j] for j in 1:n for i in (j+1):n]
        compare_values(hcat(lower.(trait_distances(d))...),readfixture(name*"_traits"))
        for squared in (true,false)
            result = trait_correlations(d; squared=squared)
            tag = name * (squared ? "_sq_" : "_raw_")
            # ade4's internal kdist.cor fuzzy path thresholds after sqrt;
            # ldist.ktab thresholds before sqrt. Raw correlations retain
            # roundoff from nominally identical fuzzy profiles in woangers.
            atol = squared ? 2e-12 : 1e-9
            compare_values(result.covariance, readfixture(tag*"paircov"); atol=atol)
            compare_values(result.correlation, readfixture(tag*"paircor"); atol=atol)
            compare_values(reshape(result.global_correlation,:,1), readfixture(tag*"glocor"); atol=atol)
        end
    end
end

@testset "EcoDists" begin
    @testset "ade4 numerical references" begin
        blocks = mixed_blocks()
        reference_case("mixed", Pavoine(blocks))
        reference_case("missing", Pavoine(mixed_blocks("missing")))
        reference_case("noscale", Pavoine(blocks; scaling=:none))
        reference_case("q_sd", Pavoine(first(blocks); scaling=:sd))
        reference_case("manhattan", Pavoine(blocks; quantitative=:manhattan))
        reference_case("podani", Pavoine(blocks; ordinal=:podani))
        for (name,method) in (("fuzzy1",:manly),("fuzzy3",:rogers),("fuzzy4",:edwards))
            reference_case(name,Pavoine(blocks[5]; fuzzy=method))
        end
        for (i,method) in enumerate(EcoDists._BINARY_METHODS)
            reference_case("binary$i", Pavoine(TraitBlock(input("binary_input"),:multichoice;col_blocks=[4]); binary=method))
        end
        raw = readdlm(joinpath(FIXTURES,"woangers_input.tsv"), '\t', Any)
        raw = map(x -> x == "NaN" ? missing : x, raw)
        wo = [TraitBlock(raw[:,[1,2,7,8]],:N), TraitBlock(raw[:,3],:circular;periods=12,origin=1),
              TraitBlock(raw[:,4:6],:O), TraitBlock(raw[:,9:19],:fuzzy;col_blocks=[3,3,5]),
              TraitBlock(raw[:,20:21],:Q)]
        reference_case("woangers",Pavoine(wo))
    end

    @testset "Distances.jl interface" begin
        d = Pavoine(mixed_blocks())
        x = traitdata(d)
        expected = pairwise(d)
        @test d isa SemiMetric
        @test pairwise(d,x; dims=1) ≈ expected
        @test pairwise(d,permutedims(x); dims=2) ≈ expected
        @test pairwise(d,[copy(x[i,:]) for i in axes(x,1)]) ≈ expected
        output = zeros(size(expected))
        @test pairwise!(d,output,x; dims=1) === output
        @test output ≈ expected
        @test pairwise(d,x[1:2,:],x[3:6,:];dims=1) ≈ expected[1:2,3:6]
        @test evaluate(d,x[1,:],x[3,:]) ≈ expected[1,3]
        @test d(x[1,:],x[3,:]) ≈ expected[1,3]
        @test colwise(d,permutedims(x[1:2,:]),permutedims(x[3:4,:])) ≈
              [expected[1,3],expected[2,4]]
        @test pairwise(Pavoine(mixed_blocks())) ≈ expected
        @test length(trait_distances(mixed_blocks())) == 14
        @test length(trait_correlations(mixed_blocks()).names) == 14
        x[1,1] = -100
        @test traitdata(d)[1,1] != -100
        @test_throws DimensionMismatch d([1.0],[2.0])
    end

    @testset "Paper formulae and weights" begin
        # Equal weights per trait, including one entire fuzzy group.
        blocks = [TraitBlock([0,4,2],:Q), TraitBlock([:a,:b,:a],:N)]
        d = Pavoine(blocks; weights=[3,1])
        @test pairwise(d)[1,2] == 1
        @test pairwise(d)[1,3] ≈ sqrt(3/16)
        @test pairwise(Pavoine(blocks;weights=[3e300,1e300])) ≈ pairwise(d)
        @test pairwise(Pavoine(blocks;weights=[1,0]))[1,3] == 0.5
        q = TraitBlock([0,1],:Q)
        f = TraitBlock([1 0 0 0; 0 1 0 0],:fuzzy;col_blocks=[4])
        @test pairwise(Pavoine([q,f]))[1,2] == 1
        # Gower's original quantitative choice is sqrt(abs difference).
        @test pairwise(Pavoine(TraitBlock([0,1,4],:Q); quantitative=:manhattan))[1,2] == 0.5
        @test pairwise(Pavoine(TraitBlock([1,12,6,7],:circular;periods=12,origin=1)))[1,2] ≈ sqrt(1/6)
        @test pairwise(Pavoine(TraitBlock([1,3,5],:circular;periods=5,origin=1)))[1,2] == 1
        @test pairwise(Pavoine(TraitBlock([1,3,5],:circular;periods=5,origin=1)))[1,3] ≈ sqrt(1/2)
        f2 = TraitBlock([1 0; 1 1; 0 1],:fuzzy;col_blocks=[2])
        @test pairwise(Pavoine(f2))[1,2] ≈ sqrt(1-inv(sqrt(2)))
        @test pairwise(Pavoine(f2;fuzzy=:hellinger))[1,2] ≈ sqrt(1-inv(sqrt(2)))
        # Weighted aggregation is reconstructed independently from trait matrices.
        w = collect(1.0:14.0)
        d = Pavoine(mixed_blocks("missing");weights=w)
        parts = trait_distances(d)
        actual = pairwise(d)
        for j in 1:6, i in (j+1):6
            available = [k for k in 1:14 if !isnan(parts[k][i,j])]
            @test actual[i,j] ≈ sqrt(sum(w[k]*parts[k][i,j]^2 for k in available)/sum(w[available]))
        end
    end

    @testset "Missingness and degenerate traits" begin
        d = Pavoine(TraitBlock([0,0,1],:D))
        @test isnan(pairwise(d)[1,2])
        @test pairwise(d)[1,3] == 1
        @test pairwise(d)[1,1] == 0
        @test isnan(trait_distances(d)[1][1,2])
        @test all(iszero, pairwise(Pavoine(TraitBlock(fill(2,3),:O))))
        @test all(iszero, pairwise(Pavoine(TraitBlock(fill(2,3),:O);ordinal=:podani)))
        @test all(iszero, pairwise(Pavoine(TraitBlock(fill(2,3),:Q);scaling=:sd)))
        @test all(iszero, pairwise(Pavoine(TraitBlock([1],:Q))))
        @test isnan(trait_correlations(TraitBlock([1],:Q)).covariance[1,1])
        @test isnan(trait_correlations(TraitBlock([1,1,1],:Q)).global_correlation[1])
        f = TraitBlock([1.0 2; 0 0; NaN 1],:fuzzy;col_blocks=[2])
        @test f.data[1,:] ≈ [1/3,2/3]
        @test all(isnan,f.data[2:3,:])
        @test isnan(pairwise(Pavoine(f))[1,2])
        @test TraitBlock([2 0; 0 3],:multichoice;col_blocks=[2]).data == [1 0;0 1]
        @test TraitBlock([2 0;0 3],:multichoice;col_blocks=2, names="diet").names == ["diet"]
        @test TraitBlock([1e308 1e308;1 0],:fuzzy;col_blocks=[2]).data[1,:] == [0.5,0.5]
        @test pairwise(Pavoine(TraitBlock([1 1;2 2;1 0],:fuzzy;col_blocks=[2])))[1,2] == 0
        @test pairwise(Pavoine(TraitBlock([1 1;2 2;1 0],:fuzzy;col_blocks=[2]);tol=0))[1,2] == 0
        # A missing nominal level remains unavailable rather than a category.
        @test isnan(TraitBlock(["a",missing,"b"],:N).data[2,1])
    end

    @testset "Input validation" begin
        @test_throws ArgumentError Pavoine(TraitBlock[])
        @test_throws ArgumentError TraitBlock([1,2],:unknown)
        @test_throws ArgumentError TraitBlock([Inf,1],:Q)
        @test_throws ArgumentError TraitBlock([0,2],:D)
        @test_throws ArgumentError TraitBlock(["a","b"],:O)
        @test_throws ArgumentError TraitBlock(zeros(0,1),:Q)
        @test_throws ArgumentError TraitBlock([1,2],:C)
        @test_throws ArgumentError TraitBlock([1,13],:circular;periods=12,origin=1)
        @test_throws ArgumentError TraitBlock([1,2.5],:circular;periods=12,origin=1)
        @test_throws ArgumentError TraitBlock([1,1],:circular;periods=1,origin=1)
        @test_throws DimensionMismatch TraitBlock(ones(3,4),:fuzzy;col_blocks=[2,1])
        @test_throws ArgumentError TraitBlock(ones(3,4),:fuzzy;col_blocks=[0,4])
        @test_throws ArgumentError TraitBlock([-1 1;1 1],:fuzzy;col_blocks=[2])
        @test_throws ArgumentError Pavoine(TraitBlock([missing,missing],:Q))
        @test_throws ArgumentError Pavoine(TraitBlock(zeros(2,2),:fuzzy;col_blocks=[2]))
        @test_throws DimensionMismatch Pavoine([TraitBlock([1,2],:Q),TraitBlock([1],:N)])
        q = TraitBlock([1,2],:Q)
        @test_throws DimensionMismatch Pavoine(q;weights=[1,2])
        for weights in ([0],[-1],[Inf],[NaN])
            @test_throws ArgumentError Pavoine(q;weights=weights)
        end
        for kwargs in ((;scaling=:bad),(;quantitative=:bad),(;ordinal=:bad),
                       (;fuzzy=:bad),(;binary=:bad),(;tol=-1),(;tol=Inf))
            @test_throws ArgumentError Pavoine(q;kwargs...)
        end
    end
end

# Independent of all existing ade4 numerical reference tests.
include("hierarchical.jl")
include("native_distances.jl")
