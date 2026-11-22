"""
    TraitDistance <: Distances.SemiMetric

A single conceptual trait of a fitted `Pavoine` model. Construct it with
`traitmetric(d,index)` or `traitmetric(d,name)` and evaluate it on prepared rows
from `traitdata(d)`. Uses the same scaling, ranks and local metric as its parent
model, independently of the trait's global weight.
"""
struct TraitDistance <: Distances.SemiMetric
    parent::Pavoine
    spec::TraitSpec
end

"""
    traitmetric(d::Pavoine, index_or_name)

Return a `TraitDistance` that works with `Distances.evaluate`, `pairwise`,
`pairwise!`, and `colwise`. Input vectors use the full fitted representation,
including when this metric selects only a subset of its columns.
"""
traitmetric(d::Pavoine,index::Integer) = TraitDistance(d,d.specs[index])
function traitmetric(d::Pavoine,name::AbstractString)
    index = findfirst(==(name),d.names)
    isnothing(index) && throw(ArgumentError("unknown trait name: $name"))
    return traitmetric(d,index)
end

function (metric::TraitDistance)(a::AbstractVector,b::AbstractVector)
    Base.require_one_based_indexing(a,b)
    length(a) == length(b) == size(metric.parent.data,2) ||
        throw(DimensionMismatch("vectors must have the fitted number of columns; use traitdata(d)"))
    a === b && return 0.0
    return sqrt(_sqtrait(metric.parent,metric.spec,a,b))
end
Distances.result_type(::TraitDistance,::Type,::Type) = Float64

"""
    trait_distances(d::Pavoine, data=d.data; dims=1)
    trait_distances(blocks; kwargs...)

Return one matrix per conceptual trait, in `d.names` order, using
`Distances.pairwise` on `TraitDistance` objects. Observations are rows by default.
A supplied matrix must use the prepared encoding from `traitdata(d)`; subsets
reuse fitted scaling and ranks. Joint D absences and unavailable pairs are NaN
off the diagonal. Weights affect only the global distance. With blocks, keywords
are forwarded to `Pavoine`.
"""
trait_distances(d::Pavoine,data::AbstractMatrix=d.data;dims::Int=1) =
    [Distances.pairwise(traitmetric(d,i),data;dims=dims) for i in eachindex(d.specs)]
trait_distances(blocks;kwargs...) = trait_distances(Pavoine(blocks;kwargs...))

function _pair_stats(x,y)
    indices = [i for i in eachindex(x,y) if isfinite(x[i]) && isfinite(y[i])]
    length(indices) >= 2 || return (NaN,NaN)
    a,b = x[indices], y[indices]
    covariance = cov(a,b)
    correlation = all(==(first(a)),a) || all(==(first(b)),b) ? NaN : cor(a,b)
    return covariance, correlation
end

"""
    trait_correlations(d::Pavoine; squared=true)
    trait_correlations(blocks; squared=true, kwargs...)

Analyze trait contributions using Pavoine et al. (2009), Eq. 9. Return a named
tuple with `covariance`, `correlation`, `global_correlation`, and trait `names`.
Statistics use unordered distinct observation pairs, squared distances by
default, and pairwise complete observations. Constant vectors or fewer than two
available pairs give NaN correlations. Covariances use the sample denominator.
Distances are evaluated directly with Distances.jl objects; no dense distance
matrices are constructed just to extract their lower triangles.
"""
function trait_correlations(d::Pavoine;squared::Bool=true)
    n = size(d.data,1)
    pairs = [(i,j) for j in 1:n for i in (j+1):n]
    transform(x) = squared ? x*x : x
    vectors = [[transform(Distances.evaluate(traitmetric(d,k),view(d.data,i,:),view(d.data,j,:)))
                for (i,j) in pairs] for k in eachindex(d.specs)]
    global_vector = [transform(Distances.evaluate(d,view(d.data,i,:),view(d.data,j,:)))
                     for (i,j) in pairs]
    p = length(vectors)
    covariance = Matrix{Float64}(undef,p,p)
    correlation = similar(covariance)
    for j in 1:p, i in j:p
        c,r = _pair_stats(vectors[i],vectors[j])
        covariance[i,j] = covariance[j,i] = c
        correlation[i,j] = correlation[j,i] = r
    end
    global_correlation = [last(_pair_stats(v,global_vector)) for v in vectors]
    return (;covariance,correlation,global_correlation,names=copy(d.names))
end
trait_correlations(blocks;squared::Bool=true,kwargs...) =
    trait_correlations(Pavoine(blocks;kwargs...);squared=squared)
