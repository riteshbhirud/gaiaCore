module EpidemicSim

using DataStructures: PriorityQueue, enqueue!, dequeue_pair!
using Statistics: mean, median, quantile
using StatsBase: wsample
using NearestNeighbors: KDTree, knn
using Random
using DataFrames

mutable struct Agent
    id::Int
    home_node::Int
    work_node::Int
    pos::Int
    direction::Symbol
    infected::Bool
    infected_in::Symbol
    transit_dependency::Float64
    income_quantile::Int
    route::Vector{Int}
    route_idx::Int
    on_vehicle::Int
    met_infected_on_vehicle::Int
    n_trips::Int
end

mutable struct Vehicle
    id::Int
    route::Vector{Int}
    pos_idx::Int
    passengers::Vector{Int}
    max_load::Int
    run_every::Float64
    line_id::Int
end

struct SimResult
    n_agents::Int
    frequency_min::Float64
    infection_curve::Vector{Float64}
    pct_transit_users::Float64
    final_infected_frac::Float64
    infection_by_income::Dict{Int,Vector{Float64}}
    infection_by_transit::Dict{Symbol,Vector{Float64}}
    tract_infections::Dict{Int,Int}
end

struct CityGraph
    n_nodes::Int
    coords::Matrix{Float64}
    neighbors::Vector{Vector{Int}}
    edge_weight::Dict{Tuple{Int,Int},Float64}
    transit_lines::Vector{Vector{Int}}
    transit_stops::Set{Int}
    node_to_location_id::Vector{Int}
end

function build_city_graph(tracts::DataFrame; k_neighbors::Int = 6)
    n = nrow(tracts)
    lats = Float64.(tracts.latitude)
    lons = Float64.(tracts.longitude)
    coords = hcat(lats, lons)'

    tree = KDTree(coords)

    neighbors = [Int[] for _ in 1:n]
    edge_weight = Dict{Tuple{Int,Int},Float64}()

    for i in 1:n
        idxs, _ = knn(tree, coords[:, i], k_neighbors + 1)
        for idx in idxs
            idx == i && continue
            dlat = (lats[idx] - lats[i]) * 111_000
            dlon = (lons[idx] - lons[i]) * 111_000 * cosd(lats[i])
            dist_m = sqrt(dlat^2 + dlon^2)
            walk_time = dist_m / 1.25

            if idx ∉ neighbors[i]
                push!(neighbors[i], idx)
                edge_weight[(i, idx)] = walk_time
            end
            if i ∉ neighbors[idx]
                push!(neighbors[idx], i)
                edge_weight[(idx, i)] = walk_time
            end
        end
    end

    pop_dens = coalesce.(tracts.pop_density, 0.0)
    sorted_idx = sortperm(pop_dens, rev=true)
    top_tracts = sorted_idx[1:min(60, n)]

    transit_lines = Vector{Int}[]
    transit_stops = Set{Int}()

    n_lines = 4
    stops_per_line = min(12, div(length(top_tracts), n_lines))

    for line_idx in 1:n_lines
        start = (line_idx - 1) * stops_per_line + 1
        stop = min(line_idx * stops_per_line, length(top_tracts))
        start > length(top_tracts) && break

        line_nodes = top_tracts[start:stop]

        if line_idx % 2 == 1
            sort!(line_nodes, by = i -> lats[i])
        else
            sort!(line_nodes, by = i -> lons[i])
        end

        push!(transit_lines, line_nodes)
        union!(transit_stops, line_nodes)

        for j in 1:(length(line_nodes)-1)
            a, b = line_nodes[j], line_nodes[j+1]
            dlat = (lats[b] - lats[a]) * 111_000
            dlon = (lons[b] - lons[a]) * 111_000 * cosd(lats[a])
            dist_m = sqrt(dlat^2 + dlon^2)
            transit_time = dist_m / 8.33
            for (x, y) in [(a, b), (b, a)]
                existing = get(edge_weight, (x, y), Inf)
                edge_weight[(x, y)] = min(existing, transit_time)
                if y ∉ neighbors[x]
                    push!(neighbors[x], y)
                end
            end
        end
    end

    return CityGraph(n, coords, neighbors, edge_weight, transit_lines,
                     transit_stops, tracts.location_id)
end

function dijkstra_path(g::CityGraph, src::Int, dst::Int)
    dist = fill(Inf, g.n_nodes)
    prev = fill(0, g.n_nodes)
    dist[src] = 0.0
    pq = PriorityQueue{Int,Float64}()
    enqueue!(pq, src, 0.0)

    while !isempty(pq)
        u, d = dequeue_pair!(pq)
        u == dst && break
        d > dist[u] && continue
        for v in g.neighbors[u]
            w = get(g.edge_weight, (u, v), Inf)
            nd = dist[u] + w
            if nd < dist[v]
                dist[v] = nd
                prev[v] = u
                if haskey(pq, v)
                    pq[v] = nd
                else
                    enqueue!(pq, v, nd)
                end
            end
        end
    end

    path = Int[]
    u = dst
    while u != 0
        pushfirst!(path, u)
        u = prev[u]
    end
    return isempty(path) || path[1] != src ? [src, dst] : path
end

function find_route_with_transit(g::CityGraph, src::Int, dst::Int, freq_sec::Float64)
    dist = fill(Inf, g.n_nodes)
    prev = fill(0, g.n_nodes)
    dist[src] = 0.0
    pq = PriorityQueue{Int,Float64}()
    enqueue!(pq, src, 0.0)

    while !isempty(pq)
        u, d = dequeue_pair!(pq)
        u == dst && break
        d > dist[u] && continue
        for v in g.neighbors[u]
            w = get(g.edge_weight, (u, v), Inf)
            wait = 0.0
            if v ∈ g.transit_stops && u ∉ g.transit_stops
                wait = freq_sec / 2.0
            end
            nd = dist[u] + w + wait
            if nd < dist[v]
                dist[v] = nd
                prev[v] = u
                if haskey(pq, v)
                    pq[v] = nd
                else
                    enqueue!(pq, v, nd)
                end
            end
        end
    end

    path = Int[]
    u = dst
    while u != 0
        pushfirst!(path, u)
        u = prev[u]
    end
    return isempty(path) || path[1] != src ? [src, dst] : path
end

function find_walk_route(g::CityGraph, src::Int, dst::Int)
    dist = fill(Inf, g.n_nodes)
    prev = fill(0, g.n_nodes)
    dist[src] = 0.0
    pq = PriorityQueue{Int,Float64}()
    enqueue!(pq, src, 0.0)

    while !isempty(pq)
        u, d = dequeue_pair!(pq)
        u == dst && break
        d > dist[u] && continue
        for v in g.neighbors[u]
            dlat = (g.coords[1, v] - g.coords[1, u]) * 111_000
            dlon = (g.coords[2, v] - g.coords[2, u]) * 111_000 * cosd(g.coords[1, u])
            dist_m = sqrt(dlat^2 + dlon^2)
            w = dist_m / 1.25
            nd = dist[u] + w
            if nd < dist[v]
                dist[v] = nd
                prev[v] = u
                if haskey(pq, v)
                    pq[v] = nd
                else
                    enqueue!(pq, v, nd)
                end
            end
        end
    end

    path = Int[]
    u = dst
    while u != 0
        pushfirst!(path, u)
        u = prev[u]
    end
    return isempty(path) || path[1] != src ? [src, dst] : path
end

function route_time(g::CityGraph, path::Vector{Int})
    t = 0.0
    for i in 1:(length(path)-1)
        t += get(g.edge_weight, (path[i], path[i+1]), 1000.0)
    end
    return t
end

function route_uses_transit(g::CityGraph, path::Vector{Int})
    for node in path
        node ∈ g.transit_stops && return true
    end
    return false
end

function choose_route(g::CityGraph, src::Int, dst::Int, freq_sec::Float64)
    walk_path = find_walk_route(g, src, dst)

    t_walk = 0.0
    for i in 1:(length(walk_path)-1)
        u, v = walk_path[i], walk_path[i+1]
        dlat = (g.coords[1, v] - g.coords[1, u]) * 111_000
        dlon = (g.coords[2, v] - g.coords[2, u]) * 111_000 * cosd(g.coords[1, u])
        t_walk += sqrt(dlat^2 + dlon^2) / 1.25
    end

    best_transit_time = Inf
    transit_path = walk_path

    for line in g.transit_lines
        src_stop_idx = 0
        dst_stop_idx = 0
        src_stop_dist = Inf
        dst_stop_dist = Inf

        for (li, stop) in enumerate(line)
            d_src = _haversine_dist(g, src, stop)
            d_dst = _haversine_dist(g, dst, stop)
            if d_src < src_stop_dist
                src_stop_dist = d_src
                src_stop_idx = li
            end
            if d_dst < dst_stop_dist
                dst_stop_dist = d_dst
                dst_stop_idx = li
            end
        end

        (src_stop_idx == 0 || dst_stop_idx == 0) && continue
        src_stop_idx == dst_stop_idx && continue

        walk_to = src_stop_dist / 1.25
        wait = freq_sec / 2.0
        lo, hi = minmax(src_stop_idx, dst_stop_idx)
        ride = 0.0
        for j in lo:(hi-1)
            ride += _haversine_dist(g, line[j], line[j+1]) / 8.33
            ride += 30.0
        end
        walk_from = dst_stop_dist / 1.25

        total = walk_to + wait + ride + walk_from

        if total < best_transit_time
            best_transit_time = total
            if src_stop_idx < dst_stop_idx
                transit_path = vcat([src], line[src_stop_idx:dst_stop_idx], [dst])
            else
                transit_path = vcat([src], reverse(line[dst_stop_idx:src_stop_idx]), [dst])
            end
        end
    end

    t_transit = best_transit_time

    σ_min = max(min(t_walk, t_transit), 1.0)
    α_walk = exp(-t_walk / σ_min)
    α_transit = exp(-t_transit / σ_min)
    total_α = α_walk + α_transit

    p_transit = α_transit / total_α

    if rand() < p_transit && isfinite(t_transit)
        return transit_path, true
    else
        return walk_path, false
    end
end

function _haversine_dist(g::CityGraph, i::Int, j::Int)
    dlat = (g.coords[1, j] - g.coords[1, i]) * 111_000
    dlon = (g.coords[2, j] - g.coords[2, i]) * 111_000 * cosd(g.coords[1, i])
    return sqrt(dlat^2 + dlon^2)
end

function run_simulation(g::CityGraph, agents::Vector{Agent};
                        freq_min::Float64,
                        n_commutes::Int = 200,
                        p0::Float64 = 0.001,
                        max_interactions_street::Int = 3,
                        max_interactions_transit::Int = 10,
                        passenger_limit::Int = 10,
                        stop_at_frac::Float64 = 0.95)
    freq_sec = freq_min * 60.0
    n_agents = length(agents)

    for a in agents
        a.pos = a.home_node
        a.direction = :to_work
        a.route = Int[]
        a.route_idx = 0
        a.on_vehicle = 0
        a.met_infected_on_vehicle = 0
        a.n_trips = 0
    end

    n_patient_zero = max(1, round(Int, 0.01 * n_agents))
    patient_zeros = randperm(n_agents)[1:n_patient_zero]
    for idx in patient_zeros
        agents[idx].infected = true
        agents[idx].infected_in = :street
    end

    used_transit = falses(n_agents)
    tract_infections = Dict{Int,Int}()
    infection_curve = Float64[]
    infection_by_income = Dict{Int,Vector{Float64}}()
    for q in 1:5
        infection_by_income[q] = Float64[]
    end
    infection_by_transit = Dict{Symbol,Vector{Float64}}()
    infection_by_transit[:high] = Float64[]
    infection_by_transit[:low] = Float64[]

    td_sorted = sort([a.transit_dependency for a in agents])
    td_threshold_high = td_sorted[max(1, round(Int, 0.75 * n_agents))]
    td_threshold_low = td_sorted[max(1, round(Int, 0.25 * n_agents))]

    node_agents = Dict{Int,Vector{Int}}()

    for commute in 1:n_commutes
        for (i, a) in enumerate(agents)
            src = a.pos
            dst = a.direction == :to_work ? a.work_node : a.home_node
            if src == dst
                a.route = [src]
                a.route_idx = 1
                continue
            end
            path, uses_t = choose_route(g, src, dst, freq_sec)
            a.route = path
            a.route_idx = 1
            if uses_t
                used_transit[i] = true
            end
            a.n_trips += 1
        end

        max_steps = maximum(length(a.route) for a in agents)

        for step in 1:max_steps
            empty!(node_agents)

            for (i, a) in enumerate(agents)
                if a.route_idx < length(a.route)
                    a.route_idx += 1
                    a.pos = a.route[a.route_idx]
                end
                if a.pos ∈ g.transit_stops && a.route_idx < length(a.route) &&
                   a.route[a.route_idx + 1] ∈ g.transit_stops
                    a.on_vehicle = 1
                else
                    if a.on_vehicle > 0 && a.met_infected_on_vehicle > 0
                        if !a.infected
                            p_ttc = 19.0 * p0
                            n_met = min(a.met_infected_on_vehicle, max_interactions_transit)
                            p_infect = 1.0 - (1.0 - p_ttc)^n_met
                            if rand() < p_infect
                                a.infected = true
                                a.infected_in = :transit
                                tract_infections[a.home_node] = get(tract_infections, a.home_node, 0) + 1
                            end
                        end
                        a.met_infected_on_vehicle = 0
                    end
                    a.on_vehicle = 0
                end

                if !haskey(node_agents, a.pos)
                    node_agents[a.pos] = Int[]
                end
                push!(node_agents[a.pos], i)
            end

            for (node, agent_indices) in node_agents
                n_at_node = length(agent_indices)
                n_at_node <= 1 && continue

                n_infected = count(i -> agents[i].infected, agent_indices)
                n_infected == 0 && continue
                n_infected == n_at_node && continue

                for i in agent_indices
                    a = agents[i]
                    if a.on_vehicle > 0
                        if !a.infected
                            a.met_infected_on_vehicle += min(n_infected, max_interactions_transit)
                        end
                    else
                        if !a.infected
                            inf_contacts = min(n_infected, max_interactions_street)
                            p_infect = 1.0 - (1.0 - p0)^inf_contacts
                            if rand() < p_infect
                                a.infected = true
                                a.infected_in = :street
                                tract_infections[a.home_node] = get(tract_infections, a.home_node, 0) + 1
                            end
                        end
                    end
                end
            end
        end

        for a in agents
            if a.on_vehicle > 0 && a.met_infected_on_vehicle > 0 && !a.infected
                p_ttc = 19.0 * p0
                n_met = min(a.met_infected_on_vehicle, max_interactions_transit)
                p_infect = 1.0 - (1.0 - p_ttc)^n_met
                if rand() < p_infect
                    a.infected = true
                    a.infected_in = :transit
                    tract_infections[a.home_node] = get(tract_infections, a.home_node, 0) + 1
                end
            end
            a.met_infected_on_vehicle = 0
            a.on_vehicle = 0
        end

        for a in agents
            dst = a.direction == :to_work ? a.work_node : a.home_node
            a.pos = dst
            a.direction = a.direction == :to_work ? :to_home : :to_work
        end

        frac_infected = count(a -> a.infected, agents) / n_agents
        push!(infection_curve, frac_infected)

        for q in 1:5
            q_agents = filter(a -> a.income_quantile == q, agents)
            isempty(q_agents) && (push!(infection_by_income[q], 0.0); continue)
            push!(infection_by_income[q], count(a -> a.infected, q_agents) / length(q_agents))
        end

        high_td = filter(a -> a.transit_dependency >= td_threshold_high, agents)
        low_td = filter(a -> a.transit_dependency <= td_threshold_low, agents)
        push!(infection_by_transit[:high],
              isempty(high_td) ? 0.0 : count(a -> a.infected, high_td) / length(high_td))
        push!(infection_by_transit[:low],
              isempty(low_td) ? 0.0 : count(a -> a.infected, low_td) / length(low_td))

        frac_infected >= stop_at_frac && break
    end

    pct_transit = 100.0 * count(used_transit) / n_agents
    final_frac = count(a -> a.infected, agents) / n_agents

    return SimResult(n_agents, freq_min, infection_curve, pct_transit,
                     final_frac, infection_by_income, infection_by_transit,
                     tract_infections)
end

function run_multi(g, agents_template, freq_min; n_runs::Int = 20, kwargs...)
    results = SimResult[]
    for run in 1:n_runs
        agents = [Agent(a.id, a.home_node, a.work_node, a.pos, a.direction,
                        false, :not, a.transit_dependency, a.income_quantile,
                        Int[], 0, 0, 0, 0) for a in agents_template]
        push!(results, run_simulation(g, agents; freq_min, kwargs...))
    end
    return results
end

function average_curves(results::Vector{SimResult})
    max_len = maximum(length(r.infection_curve) for r in results)
    curves = zeros(max_len, length(results))
    for (j, r) in enumerate(results)
        for i in 1:max_len
            curves[i, j] = i <= length(r.infection_curve) ? r.infection_curve[i] : r.infection_curve[end]
        end
    end
    return vec(mean(curves, dims=2))
end

function average_stratified(results::Vector{SimResult}, key)
    all_keys = keys(first(results).infection_by_income)
    out = Dict()
    for k in all_keys
        curves_k = [r.infection_by_income[k] for r in results]
        max_len = maximum(length.(curves_k))
        mat = zeros(max_len, length(curves_k))
        for (j, c) in enumerate(curves_k)
            for i in 1:max_len
                mat[i, j] = i <= length(c) ? c[i] : c[end]
            end
        end
        out[k] = vec(mean(mat, dims=2))
    end
    return out
end

function average_transit_stratified(results::Vector{SimResult})
    out = Dict{Symbol,Vector{Float64}}()
    for k in [:high, :low]
        curves_k = [r.infection_by_transit[k] for r in results]
        max_len = maximum(length.(curves_k))
        mat = zeros(max_len, length(curves_k))
        for (j, c) in enumerate(curves_k)
            for i in 1:max_len
                mat[i, j] = i <= length(c) ? c[i] : c[end]
            end
        end
        out[k] = vec(mean(mat, dims=2))
    end
    return out
end

end
