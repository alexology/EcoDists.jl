# Exact successive shortest augmenting paths for uniform discrete transport.
# Supplies/demands are scaled to integer units using lcm(m,n); ground costs
# are integer nucleotide mismatch counts. All residual arithmetic is integral.
function _flow_arc!(adjacency, targets, capacities, costs, from, to, cap, cost)
    push!(targets,to); push!(capacities,cap); push!(costs,cost)
    push!(adjacency[from],length(targets))
    push!(targets,from); push!(capacities,0); push!(costs,-cost)
    push!(adjacency[to],length(targets))
end

function _uniform_transport_cost(counts::AbstractMatrix{Int}, sites::Int)
    m,n = size(counts)
    m > 0 && n > 0 || throw(ArgumentError("transport supports must be nonempty"))
    sites > 0 || throw(ArgumentError("transport denominator must be positive"))
    all(x -> 0 <= x <= sites, counts) || throw(ArgumentError("ground mismatch counts must be in 0:sites"))
    if m == 1 || n == 1
        return Float64(sum(Int128(x) for x in counts)) / (Int128(m)*n*sites)
    end
    divisor = gcd(m,n)
    units = Base.Checked.checked_mul(div(m,divisor),n)
    source, sink = m+n+1, m+n+2
    vertices = sink
    adjacency = [Int[] for _ in 1:vertices]
    targets, capacities, costs = Int[], Int[], Int[]
    for i in 1:m
        _flow_arc!(adjacency,targets,capacities,costs,source,i,div(n,divisor),0)
    end
    for j in 1:n, i in 1:m
        _flow_arc!(adjacency,targets,capacities,costs,i,m+j,units,counts[i,j])
    end
    for j in 1:n
        _flow_arc!(adjacency,targets,capacities,costs,m+j,sink,div(m,divisor),0)
    end
    potentials = zeros(Int,vertices)
    distances = fill(typemax(Int),vertices)
    previous = zeros(Int,vertices)
    visited = falses(vertices)
    remaining = units
    objective = Int128(0)
    while remaining > 0
        fill!(distances,typemax(Int)); fill!(visited,false); fill!(previous,0)
        distances[source] = 0
        for _ in 1:vertices
            u = 0
            best = typemax(Int)
            for v in 1:vertices
                if !visited[v] && distances[v] < best
                    u, best = v, distances[v]
                end
            end
            iszero(u) && break
            visited[u] = true
            # Finish all reachable nodes before updating dual potentials.
            for e in adjacency[u]
                capacities[e] == 0 && continue
                v = targets[e]
                candidate = best + costs[e] + potentials[u] - potentials[v]
                if candidate < distances[v]
                    distances[v] = candidate
                    previous[v] = e
                end
            end
        end
        previous[sink] != 0 || error("internal transport error: no residual augmenting path")
        for v in 1:vertices
            distances[v] != typemax(Int) && (potentials[v] += distances[v])
        end
        delta = remaining
        v = sink
        while v != source
            e = previous[v]
            delta = min(delta,capacities[e])
            v = targets[isodd(e) ? e+1 : e-1]
        end
        v = sink
        while v != source
            e = previous[v]
            reverse = isodd(e) ? e+1 : e-1
            capacities[e] -= delta
            capacities[reverse] += delta
            objective += Int128(delta)*costs[e]
            v = targets[reverse]
        end
        remaining -= delta
    end
    return Float64(objective) / (Int128(units)*sites)
end
