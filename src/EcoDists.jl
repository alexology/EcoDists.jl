"""
    EcoDists

Mixed trait distances and trait contributions following Pavoine et al. (2009),
composed from Distances.jl primitives. R/ade4 is used only as an independent
numerical reference in optional fixture-generation tools.
The independent `HierarchicalBeta` family compares species incidence and
within-species molecular repertoires using ASV Jaccard or exact uniform transport.
"""
module EcoDists

import Distances
using Statistics: mean, std, cov, cor
using SparseArrays: SparseMatrixCSC, nzrange, rowvals, nonzeros

export TraitBlock, Pavoine, traitdata, TraitDistance, traitmetric,
       trait_distances, trait_correlations
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
