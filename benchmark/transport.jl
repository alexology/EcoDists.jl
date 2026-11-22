# Optional independent LP oracle and benchmark. HiGHS and MOI are benchmark-only.
using EcoDists
using Random: MersenneTwister, rand
using Statistics: median
using Printf: @printf
import HiGHS
import MathOptInterface as MOI

function highs_transport(counts,sites)
    m,n = size(counts)
    model = HiGHS.Optimizer()
    MOI.set(model,MOI.Silent(),true)
    MOI.set(model,MOI.RawOptimizerAttribute("threads"),1)
    MOI.set(model,MOI.RawOptimizerAttribute("solver"),"simplex")
    MOI.set(model,MOI.RawOptimizerAttribute("primal_feasibility_tolerance"),1e-9)
    MOI.set(model,MOI.RawOptimizerAttribute("dual_feasibility_tolerance"),1e-9)
    variables,_ = MOI.add_constrained_variables(model,fill(MOI.GreaterThan(0.0),m*n))
    x = reshape(variables,m,n)
    for i in 1:m
        terms = [MOI.ScalarAffineTerm(1.0,x[i,j]) for j in 1:n]
        MOI.add_constraint(model,MOI.ScalarAffineFunction(terms,0.0),MOI.EqualTo(1/m))
    end
    for j in 1:n
        terms = [MOI.ScalarAffineTerm(1.0,x[i,j]) for i in 1:m]
        MOI.add_constraint(model,MOI.ScalarAffineFunction(terms,0.0),MOI.EqualTo(1/n))
    end
    terms = [MOI.ScalarAffineTerm(counts[i,j]/sites,x[i,j]) for j in 1:n for i in 1:m]
    MOI.set(model,MOI.ObjectiveSense(),MOI.MIN_SENSE)
    MOI.set(model,MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),MOI.ScalarAffineFunction(terms,0.0))
    MOI.optimize!(model)
    MOI.get(model,MOI.TerminationStatus()) == MOI.OPTIMAL || error("HiGHS did not certify an optimal solution")
    return MOI.get(model,MOI.ObjectiveValue())
end

function validate_transport()
    rng = MersenneTwister(20261122)
    worst = 0.0
    comparisons = 0
    for trial in 1:250
        m,n = rand(rng,1:16,2)
        sites = 257
        costs = rand(rng,0:sites,m,n)
        reference = highs_transport(costs,sites)
        actual = EcoDists._uniform_transport_cost(costs,sites)
        error = abs(actual-reference)
        @assert error <= 1e-9
        worst = max(worst,error)
        comparisons += 1
    end
    println("Independent HiGHS validation: ",comparisons," cases; max absolute difference ",worst)
end

function benchmark_transport()
    println("Julia ",VERSION,"; HiGHS.jl ",pkgversion(HiGHS),"; MOI ",pkgversion(MOI))
    println("CPU: ",first(Sys.cpu_info()).model,"; HiGHS threads=1")
    println("Costs: 257-column aligned-DNA mismatch counts; model construction included; 7 warmed repetitions")
    println("m\tn\tinteger_flow_ms\thighs_lp_ms\tinteger_flow_bytes\thighs_julia_bytes")
    rng = MersenneTwister(4701)
    for (m,n) in ((3,4),(8,13),(32,32),(32,47),(64,97),(128,128))
        sites = 257
        sequences = [String(rand(rng,collect("ACGT"),sites)) for _ in 1:(m+n)]
        metric = NucleotidePDistance(sequences)
        ground = EcoDists._p_distance_matrix(metric,sequences)
        counts = ground.mismatches[1:m,(m+1):(m+n)]
        f = () -> EcoDists._uniform_transport_cost(counts,sites)
        h = () -> highs_transport(counts,sites)
        @assert isapprox(f(),h();atol=1e-9,rtol=1e-9)
        timings = [(@timed f()) for _ in 1:7]
        oracle = [(@timed h()) for _ in 1:7]
        @printf("%d\t%d\t%.4f\t%.4f\t%d\t%d\n",m,n,
                1000median(t.time for t in timings),1000median(t.time for t in oracle),
                median(t.bytes for t in timings),median(t.bytes for t in oracle))
    end
end

validate_transport()
benchmark_transport()
