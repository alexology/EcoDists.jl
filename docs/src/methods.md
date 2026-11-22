# Methods

## Aggregation and preparation

The global coefficient is the weighted root mean square of available trait
distances (Pavoine et al. 2009, Eq. 5). No additional square root or Euclidean
correction is applied after that aggregation. The defaults reproduce ade4's
noninteractive choices: Euclidean quantitative/rank distances, normalized chord
fuzzy distance, and square-root Jaccard multichoice distance.

`TraitBlock` copies input arrays and preserves conceptual groups. All numerical
work uses Float64. Nominal encodings only test equality and have no ordering.
Ordered traits must be numeric; average ranks are computed over available
observations, including ties. A group's partially missing or all-zero affinity
vector represents an unavailable trait, matching the preparation functions in ade4.

`Pavoine` stores fitted numerical data and trait metadata. `traitdata` returns
a copy suitable for the standard Distances.jl API. Model fields and prepared
data should be treated as read-only. This version supports comparisons of fitted
observations and their subsets; it does not provide an encoder for new raw rows.

## Quantitative and ordinal choices

For quantitative traits, `scaling=:range` divides by the observed range,
`:sd` centers and divides by the sample standard deviation, and `:none` leaves
values unchanged. Use `:none` with previously normalized population ranges.
Standard deviation scaling can produce distances above one and warns when used
with mixed types. Constant traits contribute zero distance and retain their weight.

For differences `z` after normalization:

| `quantitative` | Squared trait distance |
|---|---|
| `:euclidean` | `abs(z)^2` |
| `:manhattan` | `abs(z)` |

The Manhattan choice reproduces Gower's original quantitative coefficient within
the paper's root-mean-square framework. It yields a square root of the normalized
absolute difference for a single trait.

`ordinal=:rank` applies the quantitative rule to range-normalized midranks.
`ordinal=:podani` instead uses Podani's tie-corrected dissimilarity as the squared
trait distance. For distinct ranks `r` and `s`, with tie counts `t_r` and `t_s`,
the numerator is `abs(r-s) - (t_r-1)/2 - (t_s-1)/2`. The denominator uses the
extreme ranks and their tie counts. Equal ranks have zero distance. The Podani
option is symmetric but need not obey the triangle inequality.

## Circular traits

For a cycle of `m` equally spaced integer levels and absolute level difference
`h`, the squared trait distance is
`min(h, m-h) / floor(m/2)`. This is Eq. 8 with the correction for odd `m`.
The number of levels is `rangemax-rangemin+1`, including both endpoints.
The maximum theoretical distance is one for both even and odd cycles.
Continuous angles and unequally spaced circular levels are outside this interface.

## Fuzzy choices

Let `p` and `q` be nonnegative affinity vectors normalized to sum to one,
`c = dot(p,q)/(norm(p)*norm(q))`, and `A = sum(sqrt.(p .* q))`.

| `fuzzy` | Squared trait distance | Interpretation |
|---|---|---|
| `:chord` | `1-c` | Normalized Orloci chord; paper Eqs. 11–12, ade4 default |
| `:manly` | `sum(abs.(p-q))/2` | Square root of total variation after aggregation |
| `:rogers` | `sum((p-q).^2)/2` | Normalized Euclidean/Rogers |
| `:hellinger` | `1-A` | Standard normalized Hellinger |
| `:edwards` | `sqrt(1-A)` | ade4 interactive Edwards convention |

The Edwards choice intentionally follows `dist.ktab`/`ldist.ktab`: ade4 takes
another square root of `dist.prop(method=5)`. Thus its trait distance is the
fourth root of `1-A`. Use `:hellinger` for the standard normalized Hellinger
distance. Neither option is the paper's default chord distance.

## Multichoice binary choices

Write `a` for shared presences, `b` and `c` for disagreements in the two directions,
and `z` for shared absences over a modality group. The trait distance is
`sqrt(1-s)`, using these similarities:

| `binary` | Similarity `s` | ade4 method |
|---|---|---|
| `:jaccard` | `a/(a+b+c)` | 1 |
| `:simplematching` | `(a+z)/(a+b+c+z)` | 2 |
| `:sokalsneath` | `a/(a+2(b+c))` | 3 |
| `:rogerstanimoto` | `(a+z)/(a+2(b+c)+z)` | 4 |
| `:sorensen` | `2a/(2a+b+c)` | 5 |
| `:gowerlegendre` | `(a-b-c+z)/(a+b+c+z)` | 6 |
| `:ochiai` | `a/sqrt((a+b)(a+c))` | 7 |
| `:sokalsneath13` | `az/sqrt((a+b)(a+c)(z+b)(z+c))` | 8 |
| `:phi` | `(az-bc)/sqrt((a+b)(a+c)(z+b)(z+c))` | 9 |
| `:s2` | `a/(a+b+c+z)`, unit self-similarity | 10 |

Methods 6 and 9 can exceed distance one. Methods 8 and 9 can be undefined for
constant modality profiles; unavailable results are excluded from aggregation.
S2 can give positive distance between distinct observations with identical
profiles. These options follow ade4's coefficients and do not promise a metric
or Euclidean geometry. The `SemiMetric` interface exposes symmetry; it does not
assert the triangle inequality. The default Jaccard choice avoids those issues
for available, nonempty profiles.

## Contributions and missing observations

`ldist_ktab` produces one distance matrix per conceptual trait. `kdist_cor` uses
the lower triangle, excluding the diagonal, to compute sample covariances,
correlations between traits, and correlations with the global distance. The
default `squared=true` is the paper's Eq. 9. `squared=false` uses distances directly.
Pairwise complete observations are selected independently for each statistic.
An insufficient sample or zero variance yields `NaN` correlation. Constant
vectors may still have zero covariance. No significance test is implied.

## Numerical conventions and compatibility

Numerical reference tests use ade4 1.7-23 and R 4.5.0. They cover complete and
missing synthetic data, both quantitative choices, Podani ordinal distances,
all ade4 fuzzy choices, all ten binary coefficients, and `woangers` (80 species,
13 traits). They compare global matrices, every trait-distance vector, and
covariances/correlations of both raw and squared distances.

Deliberate boundary conventions:

- Entirely unavailable conceptual traits raise `ArgumentError`; they are not
  silently removed, so trait names and user-supplied weight positions remain fixed.
- Constant ordinal and quantitative traits have distance zero. Some ade4 paths
  divide by zero for constant ordinal traits; this is handled explicitly here.
- Invalid levels, infinite values, negative affinities, invalid group sizes,
  and invalid weights are rejected before pairwise evaluation.
- `tol` defaults to `1e-8`. For fuzzy traits the squared trait distance is
  thresholded before its square root, as in ade4. Other individual trait distances
  below `tol` become zero. Identical fuzzy profiles have exact zero distance.
  ade4 also applies tolerances to block sums in some paths, so results extremely
  close to a threshold can differ when traits are grouped differently. EcoDists
  applies thresholds per trait consistently, independent of block layout.
- There is no interactive `scann` prompt: all methods are keyword choices.
- ade4's internal `kdist.cor` fuzzy path thresholds after the square root,
  whereas `ldist.ktab` thresholds before it. On `woangers`, retaining numerical
  noise for identical fuzzy profiles changes raw-distance covariances/correlations
  by less than `5e-10`. EcoDists consistently derives contributions from its
  per-trait matrices. Reference tests allow `1e-9` absolute error for raw
  contributions and `2e-12` for squared contributions and distance matrices.
- Weights are an explicit extension to the ade4 interface, already present in
  the paper's formula. Finite weights are internally rescaled without changing
  their relative values, avoiding overflow in the weighted mean.
- Fuzzy auxiliary column frequencies and `row.w` metadata used elsewhere in ade4
  are not stored; they do not enter this trait-distance coefficient.

Ordination (PCoA/DPCoA), clustering, and derived diversity indices in the case
study are applications of these distances, not implementations in this package.
