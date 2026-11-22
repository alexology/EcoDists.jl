const _SpeciesContributionRow = NamedTuple{(:species,:contribution),Tuple{String,Float64}}

"""
    SpeciesContributions <: AbstractVector

Shared-species contributions as a two-column table in fitted species order.
Access columns with `.species` and `.contribution`, and rows with `table[i]`.
`sum(table.contribution)` is the molecular component. Zero contributions are
retained. `Matrix(table)` includes species labels and contributions in two
columns (element type Any). With DataFrames loaded, `DataFrame(table)` gives
typed `species::String` and `contribution::Float64` columns, including for empties.
"""
struct SpeciesContributions <: AbstractVector{_SpeciesContributionRow}
    species::Vector{String}
    contribution::Vector{Float64}
end
SpeciesContributions() = SpeciesContributions(String[],Float64[])
Base.size(t::SpeciesContributions) = (length(t.species),)
Base.IndexStyle(::Type{SpeciesContributions}) = IndexLinear()
Base.getindex(t::SpeciesContributions,i::Int) =
    (;species=t.species[i],contribution=t.contribution[i])
Base.Matrix(t::SpeciesContributions) = hcat(t.species,t.contribution)

Base.show(io::IO,t::SpeciesContributions) = print(io,"SpeciesContributions(",length(t)," rows)")
function Base.show(io::IO,::MIME"text/plain",t::SpeciesContributions)
    n = length(t)
    print(io,n,"×2 SpeciesContributions")
    visible = get(io,:limit,false) ? min(n,max(0,displaysize(io)[1]-5)) : n
    labels = [repr(t.species[i]) for i in 1:visible]
    width = max(7,maximum(textwidth.(labels);init=0))
    print(io,'\n',rpad("species",width),"  contribution")
    for i in 1:visible
        print(io,'\n',rpad(labels[i],width),"  ",t.contribution[i])
    end
    visible < n && print(io,"\n… ",n-visible," more rows")
end

"""
    HierarchicalResult

Readable result of `hierarchical_components`. Fields are `total`, `taxonomic`,
`molecular`, `species::SpeciesContributions`, `shared_mean`, `molecular_fraction`,
`nshared`, and `nunion`. Property access preserves full-precision numerical values.
`Matrix(result)` is a 1×7 numerical summary in the field order
`total, taxonomic, molecular, shared_mean, molecular_fraction, nshared, nunion`.
With DataFrames loaded, `DataFrame(result)` is a named one-row summary and
`DataFrame(result.species)` is the separate species table. DataFrames is optional.
"""
struct HierarchicalResult
    total::Float64
    taxonomic::Float64
    molecular::Float64
    species::SpeciesContributions
    shared_mean::Float64
    molecular_fraction::Float64
    nshared::Int
    nunion::Int
end

Base.Matrix(r::HierarchicalResult) = reshape(Float64[r.total,r.taxonomic,r.molecular,
    r.shared_mean,r.molecular_fraction,r.nshared,r.nunion],1,7)
Base.show(io::IO,r::HierarchicalResult) =
    print(io,"HierarchicalResult(total=",r.total,", taxonomic=",r.taxonomic,
          ", molecular=",r.molecular,")")
function Base.show(io::IO,::MIME"text/plain",r::HierarchicalResult)
    print(io,"Hierarchical beta diversity")
    for name in (:total,:taxonomic,:molecular,:shared_mean,:molecular_fraction,:nshared,:nunion)
        print(io,'\n',rpad(string(name),20),"  ",getproperty(r,name))
    end
    print(io,"\n\nShared-species contributions\n")
    show(io,MIME"text/plain"(),r.species)
end
