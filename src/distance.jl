struct TraitSpec
    kind::Symbol
    columns::UnitRange{Int}
    period::Int
    scale::Float64
    ties::Dict{Float64,Int}
end

"""
    Pavoine(blocks; weights=nothing, scaling=:range, quantitative=:euclidean,
            ordinal=:rank, fuzzy=:chord, binary=:jaccard, tol=1e-8)

Fit Pavoine et al.'s mixed-variables coefficient to a vector or tuple of
`TraitBlock`s. Each conceptual trait has equal weight by default, including
each fuzzy/multichoice group. Weights must be finite, nonnegative, and not all zero.

`scaling` is `:range`, `:sd` (sample standard deviation), or `:none` for Q;
O always uses sample midranks, with `ordinal=:rank` or `:podani`.
`quantitative=:euclidean` uses squared differences; `:manhattan` uses absolute
differences as squared trait distances (the square root is taken at aggregation).
Fuzzy choices are `:chord`, `:manly`, `:rogers`, `:hellinger`, and `:edwards`.
Binary choices are documented in the methods guide.

The result implements `Distances.SemiMetric` on the **prepared numerical rows**
returned by `traitdata`, not on raw trait rows. Use
`Distances.pairwise(d, traitdata(d); dims=1)` for the fitted distance matrix.
`Distances.pairwise(d)` is the convenience form using all fitted observations.
Undefined pairs return `NaN`; all self distances are zero.
Missing traits and joint absences for D are omitted from each pair's denominator.
"""
struct Pavoine <: Distances.SemiMetric
    data::Matrix{Float64}
    specs::Vector{TraitSpec}
    names::Vector{String}
    weights::Vector{Float64}
    quantitative::Symbol
    ordinal::Symbol
    fuzzy::Symbol
    binary::Symbol
    tol::Float64
end

const _BINARY_METHODS = (:jaccard, :simplematching, :sokalsneath, :rogerstanimoto,
                         :sorensen, :gowerlegendre, :ochiai, :sokalsneath13,
                         :phi, :s2)

function _midranks!(v)
    positions = findall(!isnan, v)
    order = sort(positions; by=i -> v[i])
    ties = Dict{Float64,Int}()
    k = 1
    while k <= length(order)
        stop = k
        while stop < length(order) && v[order[stop+1]] == v[order[k]]
            stop += 1
        end
        rank = (k + stop) / 2
        ties[rank] = stop - k + 1
        for l in k:stop
            v[order[l]] = rank
        end
        k = stop + 1
    end
    return ties
end

function Pavoine(blocks::Union{AbstractVector{TraitBlock},Tuple{Vararg{TraitBlock}}};
                 weights=nothing, scaling=:range, quantitative=:euclidean,
                 ordinal=:rank, fuzzy=:chord, binary=:jaccard, tol=1e-8)
    isempty(blocks) && throw(ArgumentError("at least one trait block is required"))
    scaling in (:range, :sd, :none) || throw(ArgumentError("unknown scaling"))
    quantitative in (:euclidean, :manhattan) || throw(ArgumentError("unknown quantitative method"))
    ordinal in (:rank, :podani) || throw(ArgumentError("unknown ordinal method"))
    fuzzy in (:chord, :manly, :rogers, :hellinger, :edwards) || throw(ArgumentError("unknown fuzzy method"))
    binary in _BINARY_METHODS || throw(ArgumentError("unknown binary method"))
    tol isa Real && isfinite(tol) && tol >= 0 || throw(ArgumentError("tol must be finite and nonnegative"))
    n = size(first(blocks).data,1)
    all(b -> size(b.data,1) == n, blocks) || throw(DimensionMismatch("blocks must have the same number of rows"))
    if scaling == :sd && length(unique(b.kind for b in blocks)) > 1
        @warn "standard deviation scaling is not recommended for mixed trait types"
    end
    data = hcat((b.data for b in blocks)...)
    specs = TraitSpec[]
    names = String[]
    offset = 0
    for b in blocks
        for (j,g) in enumerate(b.groups)
            columns = (first(g)+offset):(last(g)+offset)
            v = view(data,:,first(columns))
            observed = filter(!isnan, vec(data[:,columns]))
            isempty(observed) && throw(ArgumentError("trait $(b.names[j]) is entirely unavailable"))
            ties = Dict{Float64,Int}()
            scale = 1.0
            if b.kind == :O
                ties = _midranks!(v)
                r = filter(!isnan, v)
                low, high = extrema(r)
                if ordinal == :rank
                    scale = high - low
                    scale = scale < tol || iszero(scale) ? 1.0 : scale
                    v .= (v .- low) ./ scale
                    scale = 1.0
                else
                    scale = high - low - (ties[high]-1)/2 - (ties[low]-1)/2
                    iszero(scale) && (scale = 1.0)
                end
            elseif b.kind == :Q
                if scaling == :range
                    low, high = extrema(observed)
                    scale = high - low
                    isfinite(scale) || throw(ArgumentError("quantitative range overflows Float64"))
                    scale = scale < tol || iszero(scale) ? 1.0 : scale
                    v .= (v .- low) ./ scale
                elseif scaling == :sd
                    scale = length(observed) > 1 ? std(observed) : 0.0
                    scale = scale < tol || iszero(scale) ? 1.0 : scale
                    isfinite(scale) || throw(ArgumentError("quantitative standard deviation overflows Float64"))
                    v .= (v .- mean(observed)) ./ scale
                end
                scale = 1.0
            end
            period = b.kind == :C ? b.periods[j] : 0
            push!(specs, TraitSpec(b.kind, columns, period, scale, ties))
            push!(names, b.names[j])
        end
        offset += size(b.data,2)
    end
    # Names are qualified only when different blocks reuse their default labels.
    if !allunique(names)
        names = [string(i, ":", name) for (i,name) in enumerate(names)]
    end
    w = isnothing(weights) ? ones(length(specs)) : Float64.(collect(weights))
    length(w) == length(specs) || throw(DimensionMismatch("one weight is required per conceptual trait"))
    all(x -> isfinite(x) && x >= 0, w) && any(>(0), w) ||
        throw(ArgumentError("weights must be finite, nonnegative, and not all zero"))
    w ./= maximum(w)
    return Pavoine(data, specs, names, w, quantitative, ordinal, fuzzy, binary, Float64(tol))
end
Pavoine(block::TraitBlock; kwargs...) = Pavoine((block,); kwargs...)

"""
    traitdata(d::Pavoine)

Return a copy of the fitted numerical observations (rows). Nominal levels are
encoded, Q is scaled, O is ranked, F is normalized, and C is zero-based.
Use these rows with `Distances.evaluate`, `pairwise`, `pairwise!`, or `colwise`.
Missing entries are encoded as `NaN`. Reuse a fitted model for subsets rather
than refitting, which would change sample ranges and ranks.
"""
traitdata(d::Pavoine) = copy(d.data)

Distances.result_type(::Pavoine, ::Type, ::Type) = Float64
Distances.pairwise(d::Pavoine) = Distances.pairwise(d,d.data;dims=1)

function _binary_sq(method, a, b, c, z)
    mismatch = b + c
    total = a + mismatch + z
    s = if method == :sokalsneath
        a / (a + 2mismatch)
    elseif method == :sokalsneath13
        a*z / sqrt((a+b)*(a+c)*(z+b)*(z+c))
    elseif method == :phi
        (a*z-b*c) / sqrt((a+b)*(a+c)*(z+b)*(z+c))
    else # S2 has unit self-similarity, handled by the caller.
        a / total
    end
    return max(0.0, 1-s)
end

function _sqtrait(d::Pavoine, s::TraitSpec, a, b)
    g = s.columns
    for j in g
        (isnan(a[j]) || isnan(b[j])) && return NaN
    end
    k = first(g)
    x, y = a[k], b[k]
    q = if s.kind == :Q || (s.kind == :O && d.ordinal == :rank)
        metric = d.quantitative == :euclidean ? Distances.SqEuclidean() : Distances.Cityblock()
        Distances.evaluate(metric,x,y)
    elseif s.kind == :O
        x == y ? 0.0 : (abs(x-y) - (s.ties[x]-1)/2 - (s.ties[y]-1)/2) / s.scale
    elseif s.kind == :N
        Float64(Distances.evaluate(Distances.Hamming(),x,y))
    elseif s.kind == :D
        x == 0 && y == 0 ? NaN : Float64(Distances.evaluate(Distances.Hamming(),x,y))
    elseif s.kind == :C
        Distances.evaluate(Distances.PeriodicEuclidean(s.period),x,y) / fld(s.period,2)
    elseif s.kind == :F
        all(j -> a[j] == b[j], g) && return 0.0
        av,bv = view(a,g),view(b,g)
        if d.fuzzy == :chord
            Distances.evaluate(Distances.CosineDist(),av,bv)
        elseif d.fuzzy == :manly
            Distances.evaluate(Distances.TotalVariation(),av,bv)
        elseif d.fuzzy == :rogers
            Distances.evaluate(Distances.SqEuclidean(),av,bv) / 2
        else
            # Stable Hellinger formula: avoids cancellation near identical
            # distributions and keeps nonnegative squared distances at tol=0.
            z = Distances.evaluate(Distances.SqEuclidean(),
                                   (sqrt(x) for x in av),(sqrt(y) for y in bv)) / 2
            d.fuzzy == :edwards ? sqrt(z) : z
        end
    else
        av,bv = view(a,g),view(b,g)
        if d.binary == :jaccard
            return _trait_threshold(d,s,Distances.evaluate(Distances.Jaccard(),av,bv))
        elseif d.binary == :rogerstanimoto
            value = Distances.evaluate(Distances.RogersTanimoto(),
                                       (x != 0 for x in av),(y != 0 for y in bv))
            return _trait_threshold(d,s,value)
        elseif d.binary == :sorensen
            return _trait_threshold(d,s,Distances.evaluate(Distances.BrayCurtis(),av,bv))
        elseif d.binary in (:simplematching,:gowerlegendre)
            value = Distances.evaluate(Distances.Hamming(),av,bv) / length(g)
            return _trait_threshold(d,s,d.binary == :gowerlegendre ? 2value : value)
        elseif d.binary == :ochiai
            return _trait_threshold(d,s,Distances.evaluate(Distances.CosineDist(),av,bv))
        end
        both, onlya, onlyb, neither = 0.0, 0.0, 0.0, 0.0
        for j in g
            if a[j] > 0
                b[j] > 0 ? (both += 1) : (onlya += 1)
            else
                b[j] > 0 ? (onlyb += 1) : (neither += 1)
            end
        end
        _binary_sq(d.binary, both, onlya, onlyb, neither)
    end
    return _trait_threshold(d,s,q)
end

function _trait_threshold(d::Pavoine,s::TraitSpec,value)
    q = max(0.0,value)
    # Numerical zero convention; applied per conceptual trait, independently
    # of how a caller groups columns into blocks.
    return (s.kind == :F ? q < d.tol : sqrt(q) < d.tol) ? 0.0 : q
end

function (d::Pavoine)(a::AbstractVector, b::AbstractVector)
    Base.require_one_based_indexing(a,b)
    length(a) == length(b) == size(d.data,2) ||
        throw(DimensionMismatch("vectors must have the fitted number of columns; use traitdata(d)"))
    a === b && return 0.0
    numerator = 0.0
    denominator = 0.0
    for i in eachindex(d.specs)
        w = d.weights[i]
        iszero(w) && continue
        q = _sqtrait(d, d.specs[i], a, b)
        if !isnan(q)
            numerator += w*q
            denominator += w
        end
    end
    return iszero(denominator) ? NaN : sqrt(numerator / denominator)
end
