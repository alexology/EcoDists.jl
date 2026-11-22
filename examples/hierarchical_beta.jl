# Synthetic incidence data: taxonomic and molecular turnover are controlled
# separately. No reads or individual haplotype-frequency estimates are used.
using EcoDists
using Distances: pairwise
using SparseArrays: sparse

function hierarchical_example()
    ids = ["a0","a1","a2","b0","b1","b2","c0","d0"]
    species = ["A","A","A","B","B","B","C","D"]
    sequences = ["AAAA","AAAT","TTTT","CCCC","CCCT","GGGG","ACAC","TGTG"]
    names = ["reference","identical","small_molecular","large_molecular","taxonomic_turnover","taxonomic_only"]
    detected = [["a0","b0"],["a0","b0"],["a1","b1"],
                ["a2","b2"],["a1","c0"],["a0","c0"]]
    rows,cols = Int[],Int[]
    for (j,sample) in enumerate(detected), id in sample
        push!(rows,something(findfirst(==(id),ids)))
        push!(cols,j)
    end
    incidence = sparse(rows,cols,fill(true,length(rows)),length(ids),length(names))
    molecular = HierarchicalBeta(ids,species,sequences)
    repertoires = HierarchicalBeta(ids,species;local_method=:asv_jaccard)
    samples = communities(molecular,incidence)
    asv_samples = communities(repertoires,incidence)
    println("sample\ttaxonomic_Jaccard\thierarchical_ASV_Jaccard\thierarchical_W1\tmolecular_component")
    for j in 2:length(samples)
        r = hierarchical_components(molecular,samples[1],samples[j])
        println(names[j], '\t', r.taxonomic, '\t', repertoires(asv_samples[1],asv_samples[j]),
                '\t', r.total, '\t', r.molecular)
        @assert r.total == r.taxonomic+r.molecular
        @assert sum(r.species.contribution) == r.molecular
    end
    control = hierarchical_components(molecular,samples[1],samples[6])
    @assert control.molecular == 0
    @assert control.total == control.taxonomic  # taxonomic Jaccard recovered
    small = hierarchical_components(molecular,samples[1],samples[3])
    large = hierarchical_components(molecular,samples[1],samples[4])
    turnover = hierarchical_components(molecular,samples[1],samples[5])
    @assert small.taxonomic == large.taxonomic == 0
    @assert 0 < small.total < large.total
    @assert small.shared_mean == turnover.shared_mean == 0.25
    @assert turnover.molecular < small.molecular  # fewer shared species / union
    @assert isapprox(turnover.molecular,1/12)
    println("\nReadable decomposition for taxonomic_turnover:")
    show(stdout,MIME"text/plain"(),turnover)
    println()
    # Optional, after loading DataFrames at the top level:
    # DataFrame(turnover)          # one-row summary
    # DataFrame(turnover.species)  # species contributions
    println("\nPairwise hierarchical Wasserstein p-distances:")
    display(pairwise(molecular,samples))
    return (; molecular,repertoires,samples,incidence)
end

hierarchical_example()
