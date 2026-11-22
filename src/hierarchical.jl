"""
    MolecularCommunity

A prepared incidence community built by `community` or `communities`. Stores
only present species and distinct ASV indices. Reuse it across comparisons with
its fitted `HierarchicalBeta`; prepared communities from other models are rejected.
Treat its fields as read-only.
"""
struct MolecularCommunity
    owner::Base.RefValue{Nothing}
    species::Vector{Int}
    repertoires::Vector{Vector{Int}}
end

"""
    HierarchicalBeta(asv_ids, species, sequences=nothing;
                     local_method=:wasserstein, masks=Dict(), gaps=:error,
                     ambiguities=:error, unassigned=:error)

Fit an incidence-based hierarchical community distance, independently of `Pavoine`.
ASV identifiers must be unique strings/symbols and have one species assignment.
Species labels are strings/symbols; `missing`, `nothing`, and empty labels are
unassigned. They raise an error unless `unassigned=:drop`, which explicitly
excludes those ASVs from all subsequent communities. Unknown ASV IDs always error.

`local_method=:asv_jaccard` compares distinct ASV identifier sets and needs no
sequences. `:wasserstein` uses exact uniform optimal transport with nucleotide
p-distance; it requires aligned sequences in ASV order. Each species can have a
different alignment length. `masks` optionally maps species names to fixed masks;
`gaps` and `ambiguities` follow `NucleotidePDistance`. All species-specific ground
matrices are computed once and cached as exact mismatch counts/common denominators.
No read counts or individual haplotype frequencies enter the calculation.

Build prepared samples with `community` or `communities`, then use
`Distances.evaluate(d,a,b)` or `Distances.pairwise(d,samples)`. This symmetric
`Distances.SemiMetric` has values in [0,1]; identical sequences or masked columns
can produce pseudometric behavior. Fitted fields should be treated as read-only.
"""
struct HierarchicalBeta{M} <: Distances.SemiMetric
    asv_ids::Vector{String}
    species::Vector{String}
    lookup::Dict{String,Tuple{Int,Int}}
    alignments::Vector{NucleotidePDistance}
    ground::Vector{_PDistanceMatrix}
    dropped_asvs::Vector{String}
    owner::Base.RefValue{Nothing}
end

function _beta_label(x, name)
    x isa Union{AbstractString,Symbol} || throw(ArgumentError("$name must be a string or symbol"))
    label = String(x)
    isempty(strip(label)) && throw(ArgumentError("$name cannot be empty"))
    return label
end
_unassigned_species(x) = isnothing(x) || ismissing(x) ||
    (x isa Union{AbstractString,Symbol} && isempty(strip(String(x))))

function HierarchicalBeta(asv_ids, species, sequences=nothing;
                          local_method=:wasserstein, masks=Dict(), gaps=:error,
                          ambiguities=:error, unassigned=:error)
    local_method in (:asv_jaccard,:wasserstein) || throw(ArgumentError("unknown local molecular method"))
    unassigned in (:error,:drop) || throw(ArgumentError("unassigned must be :error or :drop"))
    gaps in (:error,:mask,:state) || throw(ArgumentError("unknown gap policy"))
    ambiguities in (:error,:mask) || throw(ArgumentError("unknown ambiguity policy"))
    asv_ids isa AbstractString && throw(ArgumentError("provide a collection of ASV IDs"))
    ids = [_beta_label(x,"ASV ID") for x in asv_ids]
    allunique(ids) || throw(ArgumentError("ASV identifiers must be unique"))
    assignments = collect(species)
    length(ids) == length(assignments) || throw(DimensionMismatch("one species assignment is required per ASV"))
    seqs = isnothing(sequences) ? nothing : collect(sequences)
    if local_method == :wasserstein
        isnothing(seqs) && throw(ArgumentError("Wasserstein requires aligned nucleotide sequences"))
        length(seqs) == length(ids) || throw(DimensionMismatch("one sequence is required per ASV"))
    elseif !isempty(masks) || !isnothing(seqs)
        throw(ArgumentError("sequences and alignment masks are only used by :wasserstein"))
    end
    groups = Vector{Int}[]
    names = String[]
    name_index = Dict{String,Int}()
    lookup = Dict{String,Tuple{Int,Int}}()
    dropped = String[]
    for i in eachindex(ids)
        label = assignments[i]
        if _unassigned_species(label)
            unassigned == :error && throw(ArgumentError("ASV $(ids[i]) has no species assignment"))
            push!(dropped,ids[i]); lookup[ids[i]] = (0,0)
            continue
        end
        name = _beta_label(label,"species label")
        k = get(name_index,name,0)
        if k == 0
            push!(names,name); push!(groups,Int[])
            k = length(names); name_index[name] = k
        end
        push!(groups[k],i)
        lookup[ids[i]] = (k,length(groups[k]))
    end
    mask_map = Dict{String,Any}()
    for (key,value) in pairs(masks)
        name = _beta_label(key,"mask species")
        haskey(name_index,name) || throw(ArgumentError("mask supplied for unknown species $name"))
        haskey(mask_map,name) && throw(ArgumentError("duplicate mask species $name"))
        mask_map[name] = value
    end
    alignments = NucleotidePDistance[]
    ground = _PDistanceMatrix[]
    if local_method == :wasserstein
        sequence_values = something(seqs)
        for (k,rows) in enumerate(groups)
            selected = sequence_values[rows]
            metric = NucleotidePDistance(selected; mask=get(mask_map,names[k],nothing),
                                         gaps=gaps, ambiguities=ambiguities)
            push!(alignments,metric)
            push!(ground,_p_distance_matrix(metric,selected))
        end
    end
    return HierarchicalBeta{local_method}(ids,names,lookup,alignments,ground,dropped,Ref(nothing))
end

function _beta_lookup(d::HierarchicalBeta,id)
    key = _beta_label(id,"ASV ID")
    haskey(d.lookup,key) || throw(ArgumentError("unknown ASV ID: $key"))
    return d.lookup[key]
end

function _prepare_community(d::HierarchicalBeta, entries)
    by_species = Dict{Int,Vector{Int}}()
    for (s,asv) in entries
        iszero(s) && continue
        push!(get!(by_species,s,Int[]),asv)
    end
    present = sort!(collect(keys(by_species)))
    repertoires = [sort!(unique!(by_species[s])) for s in present]
    return MolecularCommunity(d.owner,present,repertoires)
end

"""
    community(d::HierarchicalBeta, asv_ids)

Prepare one community from detected ASV IDs. Repeated identifiers count once;
known ASVs excluded by `unassigned=:drop` are ignored, and unknown IDs error.
An empty collection explicitly represents an empty community.
"""
function community(d::HierarchicalBeta, asv_ids)
    asv_ids isa AbstractString && throw(ArgumentError("provide a collection of ASV IDs"))
    return _prepare_community(d,(_beta_lookup(d,id) for id in asv_ids))
end

function _incidence_rows(d, incidence, asv_ids)
    Base.require_one_based_indexing(incidence)
    ids = collect(asv_ids)
    length(ids) == size(incidence,1) || throw(DimensionMismatch("ASV IDs must match incidence rows"))
    allunique(_beta_label.(ids,"ASV ID")) || throw(ArgumentError("incidence row ASV IDs must be unique"))
    return [_beta_lookup(d,id) for id in ids]
end
function _check_incidence(x)
    x isa Real && isfinite(x) && (x == 0 || x == 1) ||
        throw(ArgumentError("incidence must contain only finite 0/1 values; convert read counts explicitly"))
    return x == 1
end

"""
    communities(d::HierarchicalBeta, incidence; asv_ids=d.asv_ids)

Prepare samples from an ASV-by-sample 0/1 or Bool matrix. `asv_ids` labels rows;
it can select or reorder fitted ASVs. Counts, negatives and missing values are
rejected. CSC sparse inputs are traversed through stored entries without a dense
ASV-by-sample copy. Explicitly stored zeros are ignored.
"""
function communities(d::HierarchicalBeta, incidence::AbstractMatrix; asv_ids=d.asv_ids)
    rows = _incidence_rows(d,incidence,asv_ids)
    result = MolecularCommunity[]
    for j in axes(incidence,2)
        entries = Tuple{Int,Int}[]
        for i in axes(incidence,1)
            _check_incidence(incidence[i,j]) && push!(entries,rows[i])
        end
        push!(result,_prepare_community(d,entries))
    end
    return result
end
function communities(d::HierarchicalBeta, incidence::SparseMatrixCSC; asv_ids=d.asv_ids)
    rows = _incidence_rows(d,incidence,asv_ids)
    row_indices = rowvals(incidence)
    values = nonzeros(incidence)
    result = MolecularCommunity[]
    for j in axes(incidence,2)
        entries = Tuple{Int,Int}[]
        for p in nzrange(incidence,j)
            _check_incidence(values[p]) && push!(entries,rows[row_indices[p]])
        end
        push!(result,_prepare_community(d,entries))
    end
    return result
end

"""
    ground_distances(d::HierarchicalBeta, species)

Return a copy of a species' cached p-distance matrix, in its ASV input order.
Only available for `local_method=:wasserstein`. No sequence re-encoding occurs.
"""
function ground_distances(d::HierarchicalBeta{:wasserstein},species)
    name = _beta_label(species,"species label")
    k = findfirst(==(name),d.species)
    isnothing(k) && throw(ArgumentError("unknown species: $name"))
    return Matrix(d.ground[k])
end
ground_distances(::HierarchicalBeta{:asv_jaccard},species) =
    throw(ArgumentError("ASV Jaccard does not use molecular ground distances"))

function _shared_count(a,b)
    i,j,shared = 1,1,0
    while i <= length(a) && j <= length(b)
        if a[i] == b[j]
            shared += 1; i += 1; j += 1
        elseif a[i] < b[j]
            i += 1
        else
            j += 1
        end
    end
    return shared
end
function _local_beta(::HierarchicalBeta{:asv_jaccard},s,a,b)
    shared = _shared_count(a,b)
    union_size = length(a)+length(b)-shared
    return (union_size-shared)/union_size
end
function _local_beta(d::HierarchicalBeta{:wasserstein},s,a,b)
    a == b && return 0.0
    cached = d.ground[s]
    return _uniform_transport_cost(view(cached.mismatches,a,b),cached.sites)
end

function _check_communities(d,a,b)
    a.owner === d.owner && b.owner === d.owner ||
        throw(ArgumentError("communities must be prepared with this fitted HierarchicalBeta"))
    return nothing
end

function _beta_parts(d,a,b; contributions=nothing)
    _check_communities(d,a,b)
    shared = _shared_count(a.species,b.species)
    union_size = length(a.species)+length(b.species)-shared
    taxonomic = union_size == 0 ? 0.0 : (union_size-shared)/union_size
    i,j = 1,1
    molecular = 0.0
    local_sum = 0.0
    while i <= length(a.species) && j <= length(b.species)
        sa,sb = a.species[i],b.species[j]
        if sa == sb
            local_distance = _local_beta(d,sa,a.repertoires[i],b.repertoires[j])
            contribution = local_distance/union_size
            local_sum += local_distance
            molecular += contribution
            !isnothing(contributions) && (contributions[d.species[sa]] = contribution)
            i += 1; j += 1
        elseif sa < sb
            i += 1
        else
            j += 1
        end
    end
    return (; taxonomic, molecular, local_sum, shared, union_size)
end

function (d::HierarchicalBeta)(a::MolecularCommunity,b::MolecularCommunity)
    parts = _beta_parts(d,a,b)
    molecular = parts.union_size == 0 ? 0.0 : parts.local_sum/parts.union_size
    return min(1.0,parts.taxonomic + molecular)
end
Distances.result_type(::HierarchicalBeta, ::Type, ::Type) = Float64

"""
    hierarchical_components(d::HierarchicalBeta, a, b)

Return `total`, `taxonomic`, `molecular`, `species` (a Dict of shared-species
contributions), `shared_mean`, `molecular_fraction`, `nshared`, and `nunion`.
`total = taxonomic + molecular`; `taxonomic` is incidence species Jaccard and
`molecular = sum(values(species))`. Each shared species contributes `d_s/nunion`.
Compute this once when several diagnostics are needed, avoiding repeated transport.

Empty/empty has total zero, empty/nonempty has total one. `shared_mean` is NaN
without shared species; `molecular_fraction` is NaN when total dissimilarity is zero.
The fraction otherwise equals `molecular/total` and is not a frequency estimate.
"""
function hierarchical_components(d::HierarchicalBeta,a::MolecularCommunity,b::MolecularCommunity)
    species = Dict{String,Float64}()
    parts = _beta_parts(d,a,b;contributions=species)
    molecular = sum(values(species); init=0.0)
    taxonomic = parts.taxonomic
    total = taxonomic+molecular
    # Summing normalized positive contributions can exceed the mathematical
    # upper bound by roundoff. Correct only that last-bit excess in one term,
    # preserving both exact floating-point addition identities in the result.
    if total > 1.0
        key = first(keys(species))
        for (name,value) in species
            value > species[key] && (key = name)
        end
        while total > 1.0
            value = species[key]
            species[key] = max(0.0,min(prevfloat(value),value-(total-1.0)))
            molecular = sum(values(species);init=0.0)
            total = taxonomic+molecular
        end
    end
    shared_mean = parts.shared == 0 ? NaN : parts.local_sum/parts.shared
    molecular_fraction = total == 0 ? NaN : molecular/total
    return (; total,taxonomic,molecular,species,shared_mean,molecular_fraction,
             nshared=parts.shared,nunion=parts.union_size)
end

"""
    taxonomic_component(d::HierarchicalBeta, a, b)

Species-incidence Jaccard, with empty/empty defined as zero. Does not solve OT.
"""
function taxonomic_component(d::HierarchicalBeta,a::MolecularCommunity,b::MolecularCommunity)
    _check_communities(d,a,b)
    shared = _shared_count(a.species,b.species)
    union_size = length(a.species)+length(b.species)-shared
    return union_size == 0 ? 0.0 : (union_size-shared)/union_size
end

"""
    molecular_component(d::HierarchicalBeta, a, b)

Sum of within-shared-species dissimilarities divided by the species union size.
"""
molecular_component(d::HierarchicalBeta,a::MolecularCommunity,b::MolecularCommunity) =
    hierarchical_components(d,a,b).molecular

"""
    species_contributions(d::HierarchicalBeta, a, b)

Dictionary mapping each shared species to its contribution to the molecular
component, including zero contributions. Species absent from either site are omitted.
"""
species_contributions(d::HierarchicalBeta,a::MolecularCommunity,b::MolecularCommunity) =
    hierarchical_components(d,a,b).species

"""
    shared_species_mean(d::HierarchicalBeta, a, b)

Conditional mean of local molecular dissimilarities among shared species;
NaN if there are no shared species. Equals `molecular*nunion/nshared` when defined.
"""
shared_species_mean(d::HierarchicalBeta,a::MolecularCommunity,b::MolecularCommunity) =
    hierarchical_components(d,a,b).shared_mean

"""
    molecular_fraction(d::HierarchicalBeta, a, b)

Fraction of hierarchical dissimilarity attributable to the molecular component;
NaN when total dissimilarity is zero, otherwise `molecular/total`.
"""
molecular_fraction(d::HierarchicalBeta,a::MolecularCommunity,b::MolecularCommunity) =
    hierarchical_components(d,a,b).molecular_fraction
