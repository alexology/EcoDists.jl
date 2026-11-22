"""
    ldist_ktab(d::Pavoine)
    ldist_ktab(blocks; kwargs...)

Return one symmetric distance matrix per conceptual trait, in model order.
Weights affect the mixed distance only. Joint D absences and unavailable
observations produce `NaN` off the diagonal. Trait names are in `d.names`.
"""
function ldist_ktab(d::Pavoine)
    n = size(d.data,1)
    out = Matrix{Float64}[]
    for s in d.specs
        m = zeros(n,n)
        for j in 1:n, i in (j+1):n
            q = _sqtrait(d, s, view(d.data,i,:), view(d.data,j,:))
            m[i,j] = m[j,i] = sqrt(q)
        end
        push!(out,m)
    end
    return out
end
ldist_ktab(blocks; kwargs...) = ldist_ktab(Pavoine(blocks; kwargs...))

function _pair_stats(x,y)
    indices = [i for i in eachindex(x,y) if isfinite(x[i]) && isfinite(y[i])]
    length(indices) >= 2 || return (NaN,NaN)
    a,b = x[indices], y[indices]
    covariance = cov(a,b)
    correlation = all(==(first(a)),a) || all(==(first(b)),b) ? NaN : cor(a,b)
    return covariance, correlation
end

"""
    kdist_cor(d::Pavoine; squared=true)
    kdist_cor(blocks; squared=true, kwargs...)

Compute trait covariance (`paircov`), trait correlation (`paircor`), and
trait-to-global correlation (`glocor`) over unordered pairs of distinct
observations, as in Pavoine et al. Eq. 9 and ade4's `kdist.cor`.
The result also contains `names`. By default distances are squared first.
Pairwise complete observations are used; fewer than two available pairs or a
constant vector gives `NaN` correlation. Covariances use the sample denominator.
"""
function kdist_cor(d::Pavoine; squared::Bool=true)
    matrices = ldist_ktab(d)
    global_matrix = dist_ktab(d)
    n = size(d.data,1)
    pairs = [(i,j) for j in 1:n for i in (j+1):n]
    transform(x) = squared ? x*x : x
    vectors = [[transform(m[i,j]) for (i,j) in pairs] for m in matrices]
    global_vector = [transform(global_matrix[i,j]) for (i,j) in pairs]
    p = length(vectors)
    paircov = Matrix{Float64}(undef,p,p)
    paircor = similar(paircov)
    for j in 1:p, i in j:p
        covariance, correlation = _pair_stats(vectors[i], vectors[j])
        paircov[i,j] = paircov[j,i] = covariance
        paircor[i,j] = paircor[j,i] = correlation
    end
    glocor = [last(_pair_stats(v,global_vector)) for v in vectors]
    return (; paircov, paircor, glocor, names=copy(d.names))
end
kdist_cor(blocks; squared::Bool=true, kwargs...) =
    kdist_cor(Pavoine(blocks; kwargs...); squared=squared)
