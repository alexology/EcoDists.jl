using Random: MersenneTwister, rand, randperm
using SparseArrays: sparse, spzeros

# Independent assignment dynamic program on replicated uniform atoms. Used
# only on small supports (lcm(m,n) <= 12), never by the production solver.
function exhaustive_uniform_cost(costs,sites)
    m,n = size(costs)
    units = lcm(m,n)
    rows = repeat(collect(1:m); inner=div(units,m))
    cols = repeat(collect(1:n); inner=div(units,n))
    dp = fill(typemax(Int),1<<units)
    dp[1] = 0
    for mask in 0:((1<<units)-2)
        dp[mask+1] == typemax(Int) && continue
        i = count_ones(mask)+1
        for j in 1:units
            bit = 1<<(j-1)
            mask & bit != 0 && continue
            next = mask | bit
            dp[next+1] = min(dp[next+1], dp[mask+1]+costs[rows[i],cols[j]])
        end
    end
    return dp[end]/(units*sites)
end

@testset "Hierarchical molecular beta diversity" begin
    @testset "Fixed-mask nucleotide p-distance" begin
        d = NucleotidePDistance(["ACGT","aggt","TCGA"])
        @test evaluate(d,"ACGT","AGGT") == 1/4
        @test evaluate(d,"ACGT","TCGA") == 2/4
        @test pairwise(d,["ACGT","AGGT"]) == [0 0.25;0.25 0]
        @test_throws ArgumentError pairwise(d,["ACG?"])
        @test size(pairwise(d,String[])) == (0,0)
        @test alignment_mask(d) == [1,2,3,4]
        mask = alignment_mask(d); empty!(mask)
        @test length(alignment_mask(d)) == 4
        d = NucleotidePDistance(["AAAA","AANA","AAAT"];ambiguities=:mask)
        @test alignment_mask(d) == [1,2,4]
        @test d("AAAA","AAAT") == 1/3  # masked also for this unambiguous pair
        @test_throws ArgumentError d("NAAA","AAAA") # never refit on a pair
        @test NucleotidePDistance(["A-","AT"];gaps=:state)("A-","AT") == 0.5
        @test NucleotidePDistance(["A-","AT"];gaps=:mask)("A-","AT") == 0
        @test NucleotidePDistance(["AN","AT"];mask=[1])("AN","AT") == 0
        @test NucleotidePDistance(["AN","AT"];mask=[true,false])("AN","AT") == 0
        for sequences in (["AN","AT"], ["A-","AT"], ["AU","AT"], ["A?","AT"],
                          ["A.","AT"], [" A","AA"], ["AC",missing], ["",""], String[])
            @test_throws ArgumentError NucleotidePDistance(sequences)
        end
        @test_throws DimensionMismatch NucleotidePDistance(["A","AC"])
        @test_throws DimensionMismatch d("AAA","AAAA")
        @test_throws ArgumentError NucleotidePDistance(["NN","NN"];ambiguities=:mask)
        @test_throws ArgumentError NucleotidePDistance(["AA"];mask=Int[])
        @test_throws ArgumentError NucleotidePDistance(["AA"];mask=[2,1])
        @test_throws ArgumentError NucleotidePDistance(["AA"];mask=[1,1])
        @test_throws ArgumentError NucleotidePDistance(["AA"];mask=[0])
        @test_throws DimensionMismatch NucleotidePDistance(["AA"];mask=[true])
        @test_throws ArgumentError NucleotidePDistance(["A?","AT"];mask=[1])
        @test_throws ArgumentError NucleotidePDistance(["AA"];gaps=:delete)
        @test_throws ArgumentError NucleotidePDistance(["AA"];ambiguities=:ignore)
        @test_throws ArgumentError NucleotidePDistance(["ſA","SA"];ambiguities=:mask)
    end

    @testset "Exact uniform transportation" begin
        rng = MersenneTwister(481)
        for m in 1:4, n in 1:4, trial in 1:6
            c = rand(rng,0:8,m,n)
            actual = EcoDists._uniform_transport_cost(c,8)
            @test actual == exhaustive_uniform_cost(c,8)
            @test actual == EcoDists._uniform_transport_cost(permutedims(c),8)
        end
        # Greedy nearest matching yields 5/8 for these DNA sequences; OT is 3/8.
        d = HierarchicalBeta(["a","b","c","d"],fill("sp",4),["AAAA","AATT","AAAT","TTAA"])
        @test d(community(d,["a","b"]),community(d,["c","d"])) == 3/8
        # Distinct identifiers with duplicate sequences retain uniform ASV mass.
        d = HierarchicalBeta(["a","b","copy"],fill("sp",3),["AAAA","TTTT","AAAA"])
        a,b = community(d,["a","b"]),community(d,["a","b","copy"])
        @test d(a,b) == 1/6
        @test d(community(d,["a"]),community(d,["copy"])) == 0
        @test d(community(d,["a"]),community(d,["a","b"])) == 1/2
    end

    ids = ["a1","a2","a3","b1","b2","c1","d1"]
    assignments = ["a","a","a","b","b","c","d"]
    seqs = ["AAAA","AAAT","TTTT","CCCC","CCCT","AA","GGGGG"]
    for method in (:asv_jaccard,:wasserstein)
        @testset "$method decomposition and interface" begin
            d = method == :wasserstein ? HierarchicalBeta(ids,assignments,seqs) :
                HierarchicalBeta(ids,assignments;local_method=method)
            a = community(d,["a1","b1"])
            b = community(d,["a2","b1"])
            local_d = method == :wasserstein ? 1/4 : 1.0
            parts = hierarchical_components(d,a,b)
            @test parts.total == local_d/2
            @test parts.taxonomic == 0
            @test parts.molecular == local_d/2
            @test parts.total == parts.taxonomic+parts.molecular
            @test parts.molecular == sum(parts.species.contribution)
            @test collect(parts.species) == [(species="a",contribution=local_d/2),(species="b",contribution=0.0)]
            @test parts.shared_mean == local_d/2
            @test parts.molecular_fraction == 1
            @test molecular_component(d,a,b) == parts.molecular
            @test taxonomic_component(d,a,b) == parts.taxonomic
            @test species_contributions(d,a,b) == parts.species
            @test shared_species_mean(d,a,b) == parts.shared_mean
            @test molecular_fraction(d,a,b) == parts.molecular_fraction
            @test evaluate(d,a,b) == d(a,b) == parts.total
            @test d(a,a) == 0
            @test d(a,community(d,["a1","b1","a1"])) == 0
            disjoint = community(d,["c1","d1"])
            @test d(a,disjoint) == taxonomic_component(d,a,disjoint) == 1
            @test molecular_component(d,a,disjoint) == 0
            @test isnan(shared_species_mean(d,a,disjoint))
            empty = community(d,String[])
            @test d(empty,empty) == 0
            @test d(a,empty) == d(empty,a) == 1
            @test isnan(molecular_fraction(d,empty,empty))
            @test isnan(molecular_fraction(d,a,a))
            one = community(d,["a2","c1"])
            shared = hierarchical_components(d,a,one)
            @test shared.taxonomic == 2/3
            @test shared.molecular == local_d/3
            @test shared.shared_mean == local_d
            @test shared.nshared == 1 && shared.nunion == 3
            @test shared.molecular < parts.molecular
            @test shared.total == shared.taxonomic+shared.molecular
            # A shared species with a single ASV has zero local dissimilarity.
            @test d(community(d,["b1","c1"]),community(d,["b1","d1"])) == 2/3
            samples = [a,b,one,disjoint,empty]
            matrix = pairwise(d,samples)
            @test matrix ≈ permutedims(matrix)
            buffer = zeros(5,5)
            @test pairwise!(d,buffer,samples) === buffer
            @test buffer == matrix
            @test pairwise(d,samples[1:2],samples[3:5]) == matrix[1:2,3:5]
            if method == :wasserstein
                @test ground_distances(d,"a") == [0 0.25 1;0.25 0 0.75;1 0.75 0]
                copied = ground_distances(d,"a"); copied[1,2]=99
                @test ground_distances(d,"a")[1,2] == 0.25
                @test_throws ArgumentError ground_distances(d,"unknown")
            else
                @test_throws ArgumentError ground_distances(d,"a")
            end
        end
    end

    @testset "Incidence preparation, sparse matrices, and input validation" begin
        d = HierarchicalBeta(ids,assignments,seqs)
        x = [1 0 1;0 1 0;1 0 0;0 1 1;0 0 0;0 1 0;0 0 0]
        dense = communities(d,x)
        @test pairwise(d,communities(d,sparse(x))) == pairwise(d,dense)
        @test pairwise(d,communities(d,Bool.(x))) == pairwise(d,dense)
        @test d(dense[1],community(d,["a1","a3"])) == 0
        order = [7,2,5,1,6,4,3]
        @test pairwise(d,communities(d,x[order,:];asv_ids=ids[order])) == pairwise(d,dense)
        @test length(communities(d,zeros(Int,7,0))) == 0
        @test d(only(communities(d,zeros(Int,7,1))),community(d,String[])) == 0
        @test_throws ArgumentError communities(d,2x)
        @test_throws ArgumentError communities(d,sparse(2x))
        @test_throws ArgumentError communities(d,fill(NaN,7,1))
        @test_throws ArgumentError communities(d,fill(missing,7,1))
        @test_throws DimensionMismatch communities(d,zeros(2,3))
        @test_throws ArgumentError communities(d,zeros(2,3);asv_ids=["a1","a1"])
        @test_throws ArgumentError community(d,["unknown"])
        @test_throws ArgumentError community(d,"a1")
        other = HierarchicalBeta(ids,assignments,seqs)
        @test_throws ArgumentError d(dense[1],community(other,["a1"]))
        @test_throws ArgumentError taxonomic_component(d,dense[1],community(other,["a1"]))
        @test_throws ArgumentError HierarchicalBeta(["a","a"],["sp","sp"],["A","A"])
        @test_throws ArgumentError HierarchicalBeta(["a"],[missing],["A"])
        @test_throws ArgumentError HierarchicalBeta(["a"],[nothing],["A"])
        @test_throws ArgumentError HierarchicalBeta(["a"],[" "],["A"])
        @test_throws DimensionMismatch HierarchicalBeta(["a"],String[],["A"])
        @test_throws ArgumentError HierarchicalBeta(["a"],["sp"])
        @test_throws ArgumentError HierarchicalBeta(["a"],["sp"],[missing])
        @test_throws DimensionMismatch HierarchicalBeta(["a","b"],["sp","sp"],["A","AA"])
        @test_throws ArgumentError HierarchicalBeta(["a"],["sp"],["A"];masks=Dict("unknown"=>[1]))
        @test_throws ArgumentError HierarchicalBeta(["a"],["sp"];local_method=:bad)
        @test_throws ArgumentError HierarchicalBeta(["a"],["sp"],["A"];unassigned=:bad)
        drop = HierarchicalBeta(["a","u"],["sp",missing],["A",missing];unassigned=:drop)
        @test drop.dropped_asvs == ["u"]
        @test drop(community(drop,["a","u"]),community(drop,["a"])) == 0
        @test drop(community(drop,["u"]),community(drop,String[])) == 0
        @test length(communities(drop,sparse([1 0;1 1]))) == 2
        empty_model = HierarchicalBeta(String[],String[],String[])
        @test empty_model(community(empty_model,String[]),community(empty_model,String[])) == 0
        masked = HierarchicalBeta(["a","b","c"],fill("sp",3),["AAAA","AANA","AAAT"];ambiguities=:mask)
        @test ground_distances(masked,"sp")[1,3] == 1/3
        explicit = HierarchicalBeta(["a","b"],fill("sp",2),["AN","AT"];masks=Dict("sp"=>[1]))
        @test explicit(community(explicit,["a"]),community(explicit,["b"])) == 0
        # Large, sparse incidence stays sparse during preparation.
        many = HierarchicalBeta(["id$i" for i in 1:20_000],fill("sp",20_000);local_method=:asv_jaccard)
        sparse_input = sparse([1,20_000],[1,2],[1,1],20_000,2)
        prepared = communities(many,sparse_input)
        @test many(prepared[1],prepared[2]) == 1
        @test sum(length(c.species) for c in prepared) == 2
        @test sum(length.(vcat((c.repertoires for c in prepared)...))) == 2
    end

    @testset "Exhaustive small communities and randomized metric properties" begin
        rng = MersenneTwister(22092018)
        for method in (:asv_jaccard,:wasserstein)
            all_ids = ["v$i" for i in 1:5]
            labels = ["a","a","b","b","c"]
            sequences = ["AAAA","ATTT","CCCC","CCCT","GGGG"]
            d = method == :wasserstein ? HierarchicalBeta(all_ids,labels,sequences) :
                HierarchicalBeta(all_ids,labels;local_method=method)
            samples = [community(d,[all_ids[i] for i in 1:5 if mask & (1<<(i-1)) != 0]) for mask in 0:31]
            matrix = pairwise(d,samples)
            @test all(x -> 0 <= x <= 1, matrix)
            @test matrix == permutedims(matrix)
            @test all(matrix[i,i] == 0 for i in 1:32)
            for i in 1:32, j in 1:32
                # Check both calls directly: pairwise's symmetry optimization
                # alone would hide an accidentally asymmetric implementation.
                @test d(samples[i],samples[j]) == d(samples[j],samples[i])
                parts = hierarchical_components(d,samples[i],samples[j])
                @test parts.total == parts.taxonomic+parts.molecular
                @test parts.molecular == sum(parts.species.contribution;init=0.0)
                @test matrix[i,j] ≈ parts.total
                sa = Set(labels[k] for k in 1:5 if (i-1) & (1<<(k-1)) != 0)
                sb = Set(labels[k] for k in 1:5 if (j-1) & (1<<(k-1)) != 0)
                u = union(sa,sb)
                expected = isempty(u) ? 0.0 : length(symdiff(sa,sb))/length(u)
                @test parts.taxonomic == expected
            end
            for i in 1:32, j in 1:32, k in 1:32
                @test matrix[i,k] <= matrix[i,j]+matrix[j,k]+1e-14
            end
            # Independent random aligned DNA, unequal supports, and shared ASVs.
            ids2 = ["r$i" for i in 1:15]
            species2 = repeat(["a","b","c"];inner=5)
            seq2 = [String(rand(rng,collect("ACGT"),20)) for _ in 1:15]
            random_model = method == :wasserstein ? HierarchicalBeta(ids2,species2,seq2) :
                HierarchicalBeta(ids2,species2;local_method=method)
            random_samples = [community(random_model,ids2[rand(rng,Bool,15)]) for _ in 1:45]
            distances = pairwise(random_model,random_samples)
            for _ in 1:1500
                i,j,k = rand(rng,1:45,3)
                @test distances[i,k] <= distances[i,j]+distances[j,k]+1e-13
                @test distances[i,j] == distances[j,i]
            end
        end
    end

    @testset "Upper bound under accumulated floating-point roundoff" begin
        ids = ["$(prefix)$i" for i in 1:8 for prefix in ("a","b")]
        species = ["sp$i" for i in 1:8 for _ in 1:2]
        push!(ids,"extra");push!(species,"extra")
        d = HierarchicalBeta(ids,species;local_method=:asv_jaccard)
        a = community(d,vcat(["a$i" for i in 1:8],["extra"]))
        b = community(d,["b$i" for i in 1:8])
        parts = hierarchical_components(d,a,b)
        @test d(a,b) == 1
        @test 0 <= parts.total <= 1
        @test parts.total == parts.taxonomic+parts.molecular
        @test parts.molecular == sum(parts.species.contribution)
        @test all(v -> isapprox(v,1/9;atol=2e-15),parts.species.contribution)
    end
end
