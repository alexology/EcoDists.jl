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
returned by `traitdata`, not on raw trait rows. Use `dist_ktab(d)` for the fitted
distance matrix. Undefined pairs return `NaN`; all self distances are zero.
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

function _binary_sq(method, a, b, c, z)
    mismatch = b + c
    total = a + mismatch + z
    s = if method == :jaccard
        a / (a + mismatch)
    elseif method == :simplematching
        (a + z) / total
    elseif method == :sokalsneath
        a / (a + 2mismatch)
    elseif method == :rogerstanimoto
        (a + z) / (a + 2mismatch + z)
    elseif method == :sorensen
        2a / (2a + mismatch)
    elseif method == :gowerlegendre
        (a - mismatch + z) / total
    elseif method == :ochiai
        a / sqrt((a+b)*(a+c))
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
        z = abs(x-y)
        d.quantitative == :euclidean ? z*z : z
    elseif s.kind == :O
        x == y ? 0.0 : (abs(x-y) - (s.ties[x]-1)/2 - (s.ties[y]-1)/2) / s.scale
    elseif s.kind == :N
        Float64(x != y)
    elseif s.kind == :D
        x == 0 && y == 0 ? NaN : Float64(x != y)
    elseif s.kind == :C
        delta = abs(x-y)
        min(delta, s.period-delta) / fld(s.period,2)
    elseif s.kind == :F
        all(j -> a[j] == b[j], g) && return 0.0
        if d.fuzzy == :chord
            dotxy = 0.0
            xx = 0.0
            yy = 0.0
            for j in g
                dotxy += a[j]*b[j]
                xx += a[j]^2
                yy += b[j]^2
            end
            1 - dotxy / sqrt(xx) / sqrt(yy)
        elseif d.fuzzy == :manly
            sum(abs(a[j]-b[j]) for j in g) / 2
        elseif d.fuzzy == :rogers
            sum((a[j]-b[j])^2 for j in g) / 2
        else
            z = max(0.0, 1 - sum(sqrt(a[j]*b[j]) for j in g))
            d.fuzzy == :edwards ? sqrt(z) : z
        end
    else
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
    q = max(0.0, q)
    # ade4 thresholds fuzzy squared distances before taking their square root.
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

"""
    dist_ktab(d::Pavoine)
    dist_ktab(blocks; kwargs...)

Return the symmetric matrix of mixed trait distances. Observations are ordered
as the input rows. The diagonal is zero; pairs without shared available traits
have distance `NaN`. Keywords are forwarded to `Pavoine`.
"""
dist_ktab(d::Pavoine) = Distances.pairwise(d, d.data; dims=1)
dist_ktab(blocks; kwargs...) = dist_ktab(Pavoine(blocks; kwargs...))
