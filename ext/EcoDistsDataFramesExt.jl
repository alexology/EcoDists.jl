module EcoDistsDataFramesExt

import EcoDists
import DataFrames

function DataFrames.DataFrame(t::EcoDists.SpeciesContributions;copycols::Bool=true)
    return DataFrames.DataFrame((;species=t.species,contribution=t.contribution);copycols=copycols)
end

function DataFrames.DataFrame(r::EcoDists.HierarchicalResult;copycols::Bool=true)
    return DataFrames.DataFrame((;total=[r.total],taxonomic=[r.taxonomic],
        molecular=[r.molecular],shared_mean=[r.shared_mean],
        molecular_fraction=[r.molecular_fraction],nshared=[r.nshared],nunion=[r.nunion]);
        copycols=copycols)
end

end
