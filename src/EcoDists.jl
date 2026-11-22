"""
    EcoDists

Mixed trait distances and trait contributions following Pavoine et al. (2009),
with the Distances.jl interface and counterparts of ade4's `dist.ktab` family.
The independent `HierarchicalBeta` family compares species incidence and
within-species molecular repertoires using ASV Jaccard or exact uniform transport.
"""
module EcoDists

import Distances
using Statistics: mean, std, cov, cor
using SparseArrays: SparseMatrixCSC, nzrange, rowvals, nonzeros

export TraitBlock, prep_fuzzy, prep_binary, prep_circular,
       Pavoine, traitdata, dist_ktab, ldist_ktab, kdist_cor
export NucleotidePDistance, alignment_mask, HierarchicalBeta, MolecularCommunity,
       community, communities, ground_distances, hierarchical_components,
       taxonomic_component, molecular_component, species_contributions,
       shared_species_mean, molecular_fraction

include("traits.jl")
include("distance.jl")
include("contributions.jl")
include("nucleotide.jl")
include("transport.jl")
include("hierarchical.jl")

end
