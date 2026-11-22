using EcoDists
using Distances: pairwise

blocks = [
    TraitBlock([10.0,20.0,15.0], :Q; names=["height"]),
    TraitBlock(["annual","perennial","annual"], :N; names=["life_cycle"]),
    prep_circular([12,1,6], 1,12; names=["flowering_month"]),
    prep_fuzzy([3 1 0;0 2 2;1 0 3], [3]; names=["diet"]),
]
d = Pavoine(blocks)
x = traitdata(d)
display(dist_ktab(d))
@assert pairwise(d,x; dims=1) ≈ dist_ktab(d)
println("Trait contributions (correlations of squared distances):")
for (name, contribution) in zip(d.names,kdist_cor(d).glocor)
    println(name, ": ", contribution)
end
