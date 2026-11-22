_ismissing(x) = ismissing(x) || (x isa AbstractFloat && isnan(x))
_matrix(x::AbstractVector) = reshape(x, :, 1)
_matrix(x::AbstractMatrix) = x

function _numeric(x)
    a = _matrix(x)
    Base.require_one_based_indexing(a)
    all(>(0), size(a)) || throw(ArgumentError("trait data must have rows and columns"))
    out = Matrix{Float64}(undef, size(a))
    for i in eachindex(a)
        v = a[i]
        if _ismissing(v)
            out[i] = NaN
        else
            v isa Real || throw(ArgumentError("numeric traits require real values"))
            z = Float64(v)
            isfinite(z) || throw(ArgumentError("trait values must be finite or missing"))
            out[i] = z
        end
    end
    return out
end

function _names(names, kind, n)
    out = isnothing(names) ? [string(kind, i) for i in 1:n] :
          names isa AbstractString ? [String(names)] : String.(collect(names))
    length(out) == n || throw(DimensionMismatch("one name is required per trait"))
    allunique(out) || throw(ArgumentError("trait names must be unique within a block"))
    return out
end

function _groups(widths, p)
    w = widths isa Integer ? [widths] : collect(widths)
    !isempty(w) && all(x -> x isa Integer && !(x isa Bool) && x > 0, w) ||
        throw(ArgumentError("col_blocks must contain positive integers"))
    sum(w) == p || throw(DimensionMismatch("col_blocks must sum to the number of columns"))
    ends = cumsum(Int.(w))
    return [(ends[k] - w[k] + 1):ends[k] for k in eachindex(w)]
end

"""
    TraitBlock(data, kind; names=nothing, col_blocks=nothing, periods=nothing, origin=1)

Describe observations in rows and traits in columns. `kind` is one of `:Q`
(quantitative), `:O` (numeric ordered levels), `:N` (nominal), `:D` (asymmetric
presence/absence), `:F` (fuzzy), `:B` (multichoice), or `:C` (circular).
Vectors describe a single trait. Fuzzy and multichoice traits use contiguous
groups of modality columns specified by `col_blocks` (default: one group).
`names` labels conceptual traits, not modalities. Circular traits require the
number of levels in `periods` and use levels `origin:origin+periods-1`.

`missing` and numeric `NaN` denote unavailable values. Fuzzy rows are normalized
within each group; multichoice values are converted to presence/absence.
Partly missing or all-zero fuzzy/multichoice groups become entirely unavailable.
Nominal values are encoded separately within each column; their numerical codes
have no order. Input arrays are copied.
"""
struct TraitBlock
    data::Matrix{Float64}
    kind::Symbol
    groups::Vector{UnitRange{Int}}
    names::Vector{String}
    periods::Vector{Int}
end

function TraitBlock(data::Union{AbstractVector,AbstractMatrix}, kind::Symbol;
                    names=nothing, col_blocks=nothing, periods=nothing, origin=1)
    kind in (:Q, :O, :N, :D, :F, :B, :C) || throw(ArgumentError("unknown trait kind: $kind"))
    raw = _matrix(data)
    if kind == :N
        Base.require_one_based_indexing(raw)
        all(>(0), size(raw)) || throw(ArgumentError("trait data must have rows and columns"))
        a = fill(NaN, size(raw))
        for j in axes(raw, 2)
            codes = Dict{Any,Int}()
            for i in axes(raw, 1)
                v = raw[i,j]
                if !_ismissing(v)
                    v isa Real && !isfinite(v) && throw(ArgumentError("nominal values cannot be infinite"))
                    a[i,j] = get!(codes, v, length(codes) + 1)
                end
            end
        end
    else
        a = _numeric(raw)
    end
    p = size(a, 2)
    grouped = kind in (:F, :B)
    !grouped && !isnothing(col_blocks) && throw(ArgumentError("col_blocks is only valid for F and B"))
    groups = grouped ? _groups(isnothing(col_blocks) ? [p] : col_blocks, p) : [j:j for j in 1:p]
    labels = _names(names, kind, length(groups))
    if kind == :D
        all(x -> isnan(x) || x == 0 || x == 1, a) ||
            throw(ArgumentError("dichotomous values must be 0 or 1"))
    elseif grouped
        all(x -> isnan(x) || x >= 0, a) || throw(ArgumentError("F and B values must be nonnegative"))
        for g in groups, i in axes(a,1)
            row = view(a,i,g)
            if any(isnan, row) || all(iszero, row)
                row .= NaN
            elseif kind == :F
                # Scaling first avoids overflow when large affinities are added.
                row ./= maximum(row)
                row ./= sum(row)
            else
                row .= row .> 0
            end
        end
    end
    cycles = Int[]
    if kind == :C
        isnothing(periods) && throw(ArgumentError("circular traits require periods"))
        per = _expand(periods, p, "periods")
        ori = _expand(origin, p, "origin")
        all(x -> x isa Real && isfinite(x) && isinteger(x) && x >= 2, per) ||
            throw(ArgumentError("circular periods must be integers of at least 2"))
        all(x -> x isa Real && isfinite(x) && isinteger(x), ori) ||
            throw(ArgumentError("circular origins must be finite integers"))
        cycles = Int.(per)
        for j in 1:p, i in axes(a,1)
            v = a[i,j]
            isnan(v) && continue
            isinteger(v) && ori[j] <= v < ori[j] + per[j] ||
                throw(ArgumentError("circular values must be integer levels within the declared cycle"))
            a[i,j] -= ori[j]
        end
    elseif !isnothing(periods)
        throw(ArgumentError("periods is only valid for C"))
    end
    return TraitBlock(a, kind, groups, labels, cycles)
end

function _expand(x, n, name)
    out = x isa Real ? fill(x, n) : collect(x)
    length(out) == n || throw(DimensionMismatch("$name must have one entry per column"))
    return out
end

"""
    prep_fuzzy(data, col_blocks; names=nothing)

Prepare fuzzy affinity groups for `Pavoine`, like ade4's `prep.fuzzy`.
Each group is one trait regardless of its number of modalities. Row weights
used by ade4 to store auxiliary column frequencies do not affect these distances.
"""
prep_fuzzy(data, col_blocks; names=nothing) =
    TraitBlock(data, :F; col_blocks=col_blocks, names=names)

"""
    prep_binary(data, col_blocks; names=nothing)

Prepare multichoice groups for `Pavoine`, like ade4's `prep.binary`.
Positive values become one; zero remains zero. An all-zero group is unavailable.
"""
prep_binary(data, col_blocks; names=nothing) =
    TraitBlock(data, :B; col_blocks=col_blocks, names=names)

"""
    prep_circular(data; rangemin=nothing, rangemax=nothing, names=nothing)
    prep_circular(data, rangemin, rangemax; names=nothing)

Prepare discrete circular traits with inclusive endpoints, like ade4's
`prep.circular`. Endpoints may be scalars or vectors. Omitted endpoints use
observed extrema; specify them for known cycles, e.g. 1 and 12 for months.
"""
function prep_circular(data; rangemin=nothing, rangemax=nothing, names=nothing)
    a = _numeric(data)
    observed = [filter(!isnan, a[:,j]) for j in axes(a,2)]
    any(isempty, observed) && throw(ArgumentError("circular columns cannot be entirely missing"))
    lo = isnothing(rangemin) ? minimum.(observed) : _expand(rangemin, size(a,2), "rangemin")
    hi = isnothing(rangemax) ? maximum.(observed) : _expand(rangemax, size(a,2), "rangemax")
    return TraitBlock(a, :C; periods=hi .- lo .+ 1, origin=lo, names=names)
end
prep_circular(data, rangemin, rangemax; names=nothing) =
    prep_circular(data; rangemin=rangemin, rangemax=rangemax, names=names)
