"""
    EcoDists

Mixed trait distances and trait contributions following Pavoine et al. (2009),
with the Distances.jl interface and counterparts of ade4's `dist.ktab` family.
"""
module EcoDists

import Distances
using Statistics: mean, std, cov, cor

export TraitBlock, prep_fuzzy, prep_binary, prep_circular,
       Pavoine, traitdata, dist_ktab, ldist_ktab, kdist_cor

include("traits.jl")
include("distance.jl")
include("contributions.jl")

end
