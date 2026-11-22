using EcoDists
using Distances: pairwise

blocks = [
    TraitBlock([10.0,20.0,15.0], :quantitative; names=["height"]),
    TraitBlock(["annual","perennial","annual"], :nominal; names=["life_cycle"]),
    TraitBlock([12,1,6], :circular; periods=12,names=["flowering_month"]),
    TraitBlock([3 1 0;0 2 2;1 0 3], :fuzzy; col_blocks=[3],names=["diet"]),
]
d = Pavoine(blocks)
x = traitdata(d)
display(pairwise(d,x; dims=1))
@assert trait_distances(d)[1] == pairwise(traitmetric(d,"height"),x; dims=1)
println("Trait contributions (correlations of squared distances):")
for (name, contribution) in zip(d.names,trait_correlations(d).global_correlation)
    println(name, ": ", contribution)
end
