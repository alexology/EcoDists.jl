# Hierarchical molecular beta diversity

## Scope

`HierarchicalBeta` is an incidence community dissimilarity with species identity
at one level and molecular repertoires within shared species at the next. It is
a separate family from the existing `Pavoine` mixed-trait coefficient; its formula
does not modify `Pavoine`, Gower scaling, or the ade4 reference methods.
Combining different sources of diversity is related to the general motivation of
[Pavoine et al. (2009)](https://doi.org/10.1111/j.1600-0706.2008.16668.x).
The hierarchical coefficient here is the definition supplied for EcoDists, rather
than a coefficient attributed to that paper.

## Definition and exact decomposition

Let `S_A` and `S_B` be the species with at least one detected, assigned ASV in two
communities, `U = |S_A union S_B|`, and `C = S_A intersection S_B`. For each shared
species, a local dissimilarity `d_s` lies in [0,1]. For `U > 0`,

```math
D_H(A,B) = \frac{|S_A\triangle S_B|+\sum_{s\in C}d_s(A,B)}{U}
          = D_{\mathrm{inter}}+D_{\mathrm{intra}},
```

```math
D_{\mathrm{inter}}=\frac{|S_A\triangle S_B|}{U},\qquad
D_{\mathrm{intra}}=\sum_{s\in C}c_s,\qquad c_s=\frac{d_s}{U}.
```

The interspecific component is exactly species-incidence Jaccard. A species that
is not shared contributes to the taxonomic component, never the molecular one.
The molecular component is zero whenever all shared-species local distances are
zero, recovering taxonomic Jaccard. It is bounded by `|C|/U`; the total is bounded
by one. Conditional molecular turnover is `mean_shared = sum(d_s)/|C|`, so
`D_intra = (|C|/U)*mean_shared`. **At fixed conditional turnover**, reducing the
shared-species proportion reduces the molecular component. This is not a general
monotonicity guarantee if the identities or local distances of shared species change.

`hierarchical_components(d,a,b)` returns:

| Field | Meaning |
|---|---|
| `total` | Hierarchical dissimilarity |
| `taxonomic` | Species-incidence Jaccard |
| `molecular` | Molecular component |
| `species` | Dictionary of contributions `d_s/U`, including zeros for shared species |
| `shared_mean` | Conditional mean local dissimilarity among shared species |
| `molecular_fraction` | `molecular/total` |
| `nshared`, `nunion` | Species counts for the comparison |

Both empty communities have total, taxonomic and molecular components zero.
Exactly one empty community gives total and taxonomic components one, molecular
component zero. `shared_mean` is `NaN` without shared species. The molecular
fraction is `NaN` when total dissimilarity is zero and otherwise lies in [0,1].

The implementation returns `total == taxonomic + molecular` and
`molecular == sum(values(species))` for the decomposition object. Ordinary
floating-point summation can overshoot the theoretical upper bound by a few ulps;
in that boundary case one species contribution is adjusted down by the rounding
excess. The scalar distance and decomposition can differ in their last bit from
summation order. No biological threshold or approximate OT regularization is used.

## Fitted mapping and incidence interface

```julia
using EcoDists
using Distances: evaluate, pairwise
using SparseArrays: sparse

ids = ["a1","a2","b1"]
species = ["species_A","species_A","species_B"]
aligned = ["ACGT","ACGA","TTTT"]
d = HierarchicalBeta(ids,species,aligned)  # default local_method=:wasserstein
a = community(d,["a1","b1"])
b = community(d,["a2","b1"])
evaluate(d,a,b)                          # 0.125
hierarchical_components(d,a,b)

X = sparse([1 0; 0 1; 1 1])             # ASV rows, sample columns
samples = communities(d,X)
pairwise(d,samples)
```

Identifiers and species labels are strings or symbols. Fitting rejects duplicate
ASV IDs and inconsistent input lengths. An ASV has one fixed species assignment.
Communities deduplicate repeated detections of an identifier. Missing species
assignments (`missing`, `nothing`, or empty labels) error by default;
`unassigned=:drop` explicitly excludes such ASVs, recorded in `d.dropped_asvs`.
They are not pooled as an artificial species. Encode other unassignment sentinels
as `missing` before fitting; text such as `"unknown"` is otherwise an ordinary
species label. Unknown ASV IDs error even under the drop policy. Missing sequences
are not allowed for assigned species in Wasserstein mode. Sequences belonging to
explicitly dropped ASVs are not fitted.

Incidence entries must be finite 0/1 or Bool. The interface rejects read counts;
convert them to incidence explicitly, e.g. `counts .> 0`, after any detection
filtering required by your protocol. `communities(d,X; asv_ids=row_labels)` accepts
subsets or reordered rows. Sparse CSC inputs are traversed by stored entries,
without densifying the ASV-by-sample matrix. Prepared communities store only
present species and ASVs, and are bound to their fitted model. Prepare them once;
`pairwise(d, samples)`, `pairwise!`, and `evaluate` then reuse the mapping and ground
matrices. Raw incidence matrices are not passed directly to `pairwise(d,X)`.

Independent accessors are `taxonomic_component`, `molecular_component`,
`species_contributions`, `shared_species_mean`, and `molecular_fraction`.
The scalar call `d(a,b)` or `evaluate(d,a,b)` gives total dissimilarity.
`taxonomic_component` skips OT entirely. If several diagnostics are required,
call `hierarchical_components` once to avoid repeated transport solves.

## Alignment and nucleotide p-distance

Sequences must already be homologous and aligned within each species; this package
does not infer alignments, reading frames, homology, or orthology. All sequences
assigned to one species must have the same positive alignment length. Different
species can have different lengths. Ground distances are mismatch counts divided
by the number of retained alignment columns, with no evolutionary correction:

```math
p(x,y)=\frac{1}{L}\sum_{k\in M} I(x_k\ne y_k),\qquad L=|M|.
```

`NucleotidePDistance(sequences; ...)` fits a single fixed mask `M` using all input
sequences. `alignment_mask(metric)` returns its retained column indices. The
metric also implements `Distances.evaluate` and `pairwise` on aligned strings.
Wasserstein fitting encodes the aligned strings and caches the species matrices
once. Subsequent community comparisons do not encode or compare strings again.

| Input or policy | Behavior |
|---|---|
| A/C/G/T, either case | Canonical nucleotide states |
| `mask` | Boolean vector or sorted unique retained column indices |
| Default `gaps=:error` | Reject `-` in any retained column |
| `gaps=:state` | Treat `-` as a fifth distinct state; gap-gap agrees |
| `gaps=:mask` | Exclude every column with a gap in any fitted sequence of that species |
| Default `ambiguities=:error` | Reject any retained IUPAC ambiguity |
| `ambiguities=:mask` | Exclude every column with R/Y/S/W/K/M/B/D/H/V/N in any fitted sequence |
| `?`, `.`, U, whitespace, non-ASCII, other invalid symbols | Always reject, including outside the selected mask |

For `HierarchicalBeta`, `masks=Dict("species_A" => mask_A, ...)` supplies masks
per species. A valid ambiguity or gap outside an explicit mask is allowed.
The denominator always uses the final retained length. Entirely unavailable
alignments or masks retaining no columns error. There is **no pairwise deletion**.
A new string passed to a fitted nucleotide metric cannot remove additional
columns: any disallowed state at a retained column errors. Adding sequences to
the species pool and refitting with a masking policy can change the mask and
therefore previously computed distances; preserve the fitted pool and masks for
comparisons across analyses. Masking is an explicit complete-species alignment
operation, not an imputation of unknown bases.

## Local choices and interpretation

`local_method=:asv_jaccard` uses ordinary Jaccard dissimilarity of the nonempty
ASV identifier sets within a shared species. It needs no sequences and treats
different identifiers as distinct even if their nucleotide sequences coincide.

`local_method=:wasserstein` places mass `1/m` on each of the `m` distinct ASV
identifiers detected for a species at site A, and mass `1/n` on each of the `n`
at site B. It minimizes `sum(gamma[i,j]*p(i,j))` subject to these uniform
marginals and nonnegative transport. Shared ASVs are included normally; their
ground cost is zero, but the two repertoires can assign different masses to them.
All mass must be transported, including when `m != n`. No nearest-neighbor or
entropically regularized substitute is used.

ASVs with identical retained sequences have zero ground cost. They are retained
as distinct identifier atoms, so two identical sequences among three IDs have
total mass 2/3; they are not silently collapsed to one equally weighted atom.
For example `{AAAA,TTTT}` versus `{AAAA,AAAA,TTTT}` has W1 = 1/6, whereas
two singleton repertoires containing different IDs with the same sequence have
distance zero. This is pseudometric behavior on ID-labeled repertoires. On the
induced distributions over retained nucleotide states, W1 is a metric.

The hierarchical type uses `Distances.SemiMetric` conservatively to expose
symmetry without promising identity of indiscernibles for sequence-equivalent
repertoires. Tests check triangle inequalities using metric local Jaccard and
fixed-mask p-distance/W1, including unequal supports and shared ASVs. Randomized
tests are evidence of implementation behavior, not a formal proof for arbitrary
alternative local coefficients (which are not accepted by this API).

## Exact transport and dependencies

The built-in solver is a successive shortest augmenting path min-cost-flow solver
specialized to uniform distributions and integer nucleotide mismatch counts.
Set `F = lcm(m,n)`, give each source `F/m` integer units and each target `F/n`
units. Residual arcs allow previous assignments to be reversed; shortest paths
use integer reduced costs and dual potentials. The objective is accumulated in
Int128 and divided by `F*L` only at the end. Flow capacities and costs are integral,
so there is no optimization tolerance, entropy parameter, or rounding of a
floating LP plan. The final returned scalar is Float64. The cached ground matrix
stores integer numerators and its common denominator; `ground_distances(d,species)`
returns a numerical copy without re-encoding sequences.

The only new runtime dependency is the Julia standard library `SparseArrays`.
`Random` is a standard-library test dependency. No external solver is required
to use the package or execute the unit/property tests.

An independent LP implementation using
[HiGHS.jl's MathOptInterface API](https://github.com/jump-dev/HiGHS.jl)
is provided in `benchmark/transport.jl`.
[ExactOptimalTransport.jl](https://juliaoptimaltransport.github.io/ExactOptimalTransport.jl/)
was also inspected: its discrete `emd` API uses a supplied LP optimizer.
Using HiGHS directly avoids requiring that wrapper and its broader distribution
dependencies for the benchmark. HiGHS and MathOptInterface are optional benchmark
dependencies only. The internal solver was selected after the comparison below;
the benchmark also serves as an independent numerical oracle.

### Reproduce the benchmark

```text
julia benchmark/setup.jl
julia --project=benchmark benchmark/transport.jl
```

On Windows, Julia 1.13.3, Intel i7-9750H, HiGHS.jl 1.27.0 / HiGHS 1.16.0,
MathOptInterface 1.54.0, HiGHS using one thread, 2026-11-22: medians of seven
warmed runs with model construction included. Costs come from independent
257-column aligned DNA sequences. Julia allocations do not include allocations
inside the native HiGHS solver, and are not process peak memory.

| ASVs A × B | Integer flow (ms) | HiGHS LP (ms) | Flow Julia bytes | LP Julia bytes |
|---|---:|---:|---:|---:|
| 3 × 4 | 0.0047 | 0.7267 | 7,040 | 10,928 |
| 8 × 13 | 0.0474 | 1.0882 | 25,813 | 40,423 |
| 32 × 32 | 0.5675 | 5.3039 | 159,327 | 306,788 |
| 32 × 47 | 2.0009 | 7.3785 | 362,564 | 494,098 |
| 64 × 97 | 18.6566 | 32.9997 | 880,521 | 1,734,855 |
| 128 × 128 | 31.4308 | 125.3514 | 3,031,631 | 5,003,265 |

The independently formulated HiGHS LP and the integer solver agreed on 250
seeded random rectangular transport problems, maximum absolute difference
`1.11e-16`. Unit tests additionally compare to a replicated-assignment dynamic
program for small supports. These timings cover the measured supports only;
they do not establish performance for much larger species repertoires or whole
metabarcoding pipelines. The existing Pavoine implementation was not benchmarked.

## Computational complexity

For species `s` with `K_s` fitted ASVs and retained alignment length `L_s`, fitting
requires `O(K_s^2 L_s)` comparisons and `O(K_s^2)` cached ground storage, plus
temporary sequence encodings. Jaccard mode does not allocate ground matrices.
Sparse incidence preparation takes `O(number_of_rows + number_of_samples + nnz)`
lookup/traversal work, plus grouping and sorting the detected IDs. Stored zeros
are checked but ignored. Prepared communities use space proportional to detected
ASVs and species, not total ASVs times total samples.

Within a shared species, Jaccard is linear in its two sorted repertoire sizes.
For W1, write `V=m+n+2`, `E=mn+m+n`, and `F=lcm(m,n)`. The implemented min-cost-flow
method needs `O(E+V)` temporary memory. Each shortest-path step is `O(V^2+E)`;
at least one integer flow unit is sent per augmentation, giving the conservative
bound `O(F*(V^2+E))`. One-element supports have a direct exact mean-cost path.
Very large repertoires may justify a network-simplex backend in future; none is
claimed tested here. `pairwise` also allocates the requested sample-by-sample
output; use `evaluate` or batches if that output is too large. No unbounded cache
of site-pair OT solutions is maintained, and comparisons do not mutate fitted data.

## Detection and taxonomic uncertainty

This coefficient compares **observed molecular repertoires**. Uniform ASV masses
are an incidence convention, not estimates of individual counts, allele
frequencies, or population haplotype frequencies. Replicate-read abundance is not
a biological frequency proxy in this implementation. Primer bias, unequal biomass,
limited sequencing depth, errors, chimeras and denoising/filtering choices can
remove real variants or introduce apparent variants. These limitations are central
to [Elbrecht et al. (2018)](https://doi.org/10.7717/peerj.4644). Assess detection
consistency and uncertainty before interpreting molecular turnover biologically.

Assignments are assumed fixed and correct. Splitting or merging species changes
both hierarchy levels; contamination, paralogy and nonhomologous alignments can
mislead interpretation. Dropping unassigned ASVs explicitly changes the observed
community. The package does not propagate taxonomic assignment probabilities,
infer populations, or estimate genetic differentiation statistics such as F_ST.
The related framework of
[Hardy and Senterre (2007)](https://doi.org/10.1111/j.1365-2745.2007.01222.x)
connects genetic and community diversity through additive partitioning; the
incidence Jaccard/W1 coefficient here is a different construction.

## Validation and example

`test/hierarchical.jl` covers fixed masks, gaps/ambiguities, invalid alignments,
missing assignments, duplicate sequences, unequal repertoire sizes, shared ASVs,
empty communities, dense/sparse incidence, and cross-model rejection. It checks
boundedness, symmetry, additive decomposition, taxonomic Jaccard recovery, every
triple of the 32 communities from a five-ASV pool, and seeded random triples on
independently generated DNA. The existing ade4 reference tests execute unchanged.
On Julia 1.13.3, the complete suite passed 136,761 assertions: 54,634 original
mixed-trait checks and 82,127 hierarchical checks. Aqua passed 8/8, JET found no
possible errors, ExplicitImports checks passed, and every export has a docstring.
Julia 1.10 compatibility remains declared but was not executed in this validation.

Run `julia --project=. examples/hierarchical_beta.jl` for a synthetic comparison
of taxonomic Jaccard, hierarchical ASV Jaccard and hierarchical W1. Controlled
scenarios change molecular divergence while holding species composition fixed,
then change shared-species proportions at fixed conditional molecular turnover.
Assertions verify the decomposition and Jaccard recovery.

## References

- Pavoine S., Vallet J., Dufour A.-B., Gachet S. and Daniel H. (2009).
  On the challenge of treating various types of variables: application for
  improving the measurement of functional diversity. *Oikos* 118, 391–402.
  [DOI](https://doi.org/10.1111/j.1600-0706.2008.16668.x).
- Hardy O.J. and Senterre B. (2007). Characterizing the phylogenetic structure
  of communities by an additive partitioning of phylogenetic diversity.
  *Journal of Ecology* 95, 493–506.
  [DOI](https://doi.org/10.1111/j.1365-2745.2007.01222.x).
- Elbrecht V., Vamos E.E., Steinke D. and Leese F. (2018). Estimating
  intraspecific genetic diversity from community DNA metabarcoding data.
  *PeerJ* 6, e4644. [DOI](https://doi.org/10.7717/peerj.4644).
