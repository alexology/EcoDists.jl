# EcoDists.jl

Mixed ecological trait distances following Pavoine, Vallet, Dufour, Gachet and
Daniel (2009), built on [Distances.jl](https://github.com/JuliaStats/Distances.jl).
This package implements the mixed-variables coefficient, distances per trait,
and the covariance/correlation analysis of trait contributions implemented by
the `dist.ktab` family in [ade4](https://adeverse.github.io/ade4/reference/dist.ktab.html).

This is part of an experimental effort using AI-assisted development to bring
ecological and bioinformatics methods into a common Julia environment.
High performance is a development objective. A limited solver comparison for the
hierarchical family is documented in its methods guide; it is not a whole-pipeline
performance claim.

The independent `HierarchicalBeta` family combines species-incidence Jaccard
with within-shared-species ASV Jaccard or exact Wasserstein p-distance.
It preserves the existing Pavoine API and returns an additive taxonomic/molecular
decomposition, with per-species contributions. See the
[hierarchical methods guide](docs/src/hierarchical.md) and runnable
[synthetic metabarcoding example](examples/hierarchical_beta.jl).

## Installation

Julia 1.10 or later is required. Install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/alexology/EcoDists.jl")
```

EcoDists is available from GitHub and is not registered in Julia General.
Its runtime dependencies are Distances, Statistics and the SparseArrays standard
library. R is needed only to
regenerate reference test data, not to use or test the package.

## Example

```julia
using EcoDists
using Distances: pairwise, evaluate

blocks = [
    TraitBlock([10.0, 20.0, 15.0], :Q; names=["height"]),
    TraitBlock(["annual", "perennial", "annual"], :N; names=["life_cycle"]),
    TraitBlock([1, 3, 2], :O; names=["ordered_class"]),
    prep_circular([12, 1, 6], 1, 12; names=["flowering_month"]),
    prep_fuzzy([3 1 0; 0 2 2; 1 0 3], [3]; names=["diet"]),
]

d = Pavoine(blocks)                  # fit ranges, ranks and encodings once
D = dist_ktab(d)                     # observations × observations
per_trait = ldist_ktab(d)             # one matrix per conceptual trait
contributions = kdist_cor(d)         # paircov, paircor, glocor, names

X = traitdata(d)                     # prepared numeric rows, copied
@assert pairwise(d, X; dims=1) ≈ D
@assert evaluate(d, X[1,:], X[2,:]) ≈ D[1,2]

weighted = dist_ktab(blocks; weights=[2, 1, 1, 1, 3])
gower = dist_ktab(blocks; quantitative=:manhattan)
```

`Distances.evaluate` and `pairwise` take the **prepared rows from `traitdata`**.
Raw values, independently encoded nominal values, and newly observed ordinal
levels cannot be passed directly to a fitted distance. To inspect a subset,
subset `X` and reuse `d`; refitting would change sample ranges and midranks.

## Trait types

| Code | Meaning | Default trait distance |
|---|---|---|
| `:Q` | Quantitative | Absolute difference divided by observed range |
| `:O` | Ordinal, supplied as numeric ordered levels | Range-scaled midrank difference |
| `:N` | Nominal, including symmetric binary traits | 0 for agreement, 1 otherwise |
| `:D` | Asymmetric presence/absence | 0 for shared presence, 1 for disagreement; joint absence excluded |
| `:F` | Fuzzy affinities or proportions | Normalized Orloci chord distance |
| `:B` | Multichoice nominal, encoded by modality columns | Square root of Jaccard dissimilarity |
| `:C` | Discrete evenly spaced circular levels | Paper Eq. 8, including odd-cycle correction |

Observations occupy rows. Each scalar column counts as one trait. Each fuzzy
or multichoice group counts as **one trait**, regardless of modality count.
For example, `prep_fuzzy(A, [3, 5])` describes two traits using eight columns.
Use `:N` when shared absences should count as agreement and `:D` when they should
be excluded. Known circular limits must be specified even if some levels are absent.

## Missing values and weights

For each pair, the global distance is

```math
D_{ij} = \sqrt{\frac{\sum_k w_k\,\delta_{ijk}\,d_{ijk}^2}
                         {\sum_k w_k\,\delta_{ijk}}}.
```

`missing` and numeric `NaN` mark unavailable observations. A partly missing
or all-zero fuzzy/multichoice group becomes unavailable as a whole. A joint
absence of a dichotomous trait is also excluded for that pair. The weights of
the remaining traits are renormalized per pair. Pairs with no available traits
return `NaN`; self distances are zero. Entirely unavailable traits are rejected
explicitly. Weights apply to conceptual traits and may be zero.

Equal input weights do not imply equal contributions. `kdist_cor(d).glocor`
is the correlation of each squared trait-distance vector with the squared global
distance vector, following Eq. 9. These correlations are not additive percentages.

## Methods and validation

See the [methods guide](docs/src/methods.md) for alternative formulas,
normalization choices, and differences from ade4.
All exported objects have Julia docstrings. A runnable example is in
[examples/mixed_traits.jl](examples/mixed_traits.jl).

```julia
using Pkg
Pkg.test("EcoDists")
```

Tests include independently generated ade4 1.7-23 reference results for seven
trait types, missing observations, alternative scaling and metrics, all ten
binary coefficients, and the 80-species `woangers` example from the paper.
Reference regeneration: `Rscript tools/make_reference.R` from the package root.
The numerical reference version is recorded in `test/fixtures/VERSION`.

Validated on Julia 1.13.3: 136,761 test assertions passed (54,634 original mixed-trait
checks and 82,127 hierarchical checks); Aqua 8/8 passed,
JET reported no possible errors, ExplicitImports checks passed, and every export
has a docstring. Julia 1.10 compatibility is declared but has not been executed.
The methods guide documents tolerance differences in raw-distance correlations.

The additional hierarchical tests cover exact transport, fixed alignment masks,
sparse incidence, empty communities, taxonomic/molecular decomposition, duplicate
sequences, and exhaustive and randomized triangle inequalities. An optional
HiGHS LP oracle and reproducible solver benchmark are in `benchmark/`; these
optimization packages are not runtime or unit-test dependencies.

The Pavoine implementation covers trait-distance construction and trait-contribution
analysis. HierarchicalBeta adds incidence community comparisons with molecular
decomposition. Ordination, clustering, and diversity analyses discussed as
applications in the Pavoine paper are left to other packages. Missing data can make
the resulting distance matrix non-Euclidean; no Euclidean correction is applied.

## Reference

Pavoine S., Vallet J., Dufour A.-B., Gachet S. and Daniel H. (2009).
On the challenge of treating various types of variables: application for
improving the measurement of functional diversity. *Oikos* **118**, 391–402.
[doi:10.1111/j.1600-0706.2008.16668.x](https://doi.org/10.1111/j.1600-0706.2008.16668.x).

EcoDists is released under the [MIT license](LICENSE).
