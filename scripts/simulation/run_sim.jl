using Pkg
Pkg.activate(@__DIR__)

using GaiaDB
using DataFrames
using Statistics
using StatsBase
using Random
using Printf

const BASE_URL = "http://localhost:3000"
const OUTPUT_DIR = joinpath(@__DIR__, "output")
mkpath(OUTPUT_DIR)

const N_AGENTS = 2000
const σ_WALK = 15.0
const B1 = 5
const B2 = 5
const LAMBDA1 = 1.0 / 100.0
const LAMBDA2 = 2.0 / 100.0
const CS = 0.01
const N_COMMUTES = 350
const STOP_FRAC = 0.99

println("Loading SDoH data and initializing agents...")
flush(stdout)

Random.seed!(42)

tracts = GaiaDB._load_tract_sdoh(BASE_URL)
tracts = dropmissing(tracts, :pop_density)
tracts = filter(row -> row.pop_density > 0, tracts)

weights = Float64.(tracts.pop_density)
weights ./= sum(weights)

valid_incomes = collect(skipmissing(tracts.median_income))
quintile_fn = GaiaDB._income_to_quintile(valid_incomes)
pct_nv = coalesce.(tracts.pct_no_vehicle, 0.0)
med_inc = coalesce.(tracts.median_income, median(valid_incomes))
n_tracts = nrow(tracts)

home_idx = wsample(1:n_tracts, weights, N_AGENTS)

mutable struct Agent
    id::Int
    income_quantile::Int
    transit_dependency::Float64
    infected::Bool
    infected_in::Symbol
end

function make_agents()
    agents = Agent[]
    for i in 1:N_AGENTS
        hi = home_idx[i]
        td = clamp(pct_nv[hi] / 100.0 + randn() * 0.05, 0.0, 1.0)
        iq = quintile_fn(med_inc[hi])
        push!(agents, Agent(i, iq, td, false, :not))
    end
    n_p0 = max(1, round(Int, CS * N_AGENTS))
    for idx in randperm(N_AGENTS)[1:n_p0]
        agents[idx].infected = true
    end
    return agents
end

println("  $(N_AGENTS) agents initialized from $(n_tracts) Cook County tracts")
println("  Transit dep mean: $(round(mean(pct_nv[home_idx])/100, digits=3))")
flush(stdout)

function run_one(σ_transit::Float64; seed::Int=0)
    seed > 0 && Random.seed!(seed)
    agents = make_agents()

    σ_min = min(σ_WALK, σ_transit)

    denom = exp(-σ_WALK / σ_min) + exp(-σ_transit / σ_min)
    α_walk_base = exp(-σ_WALK / σ_min) / denom
    α_transit_base = exp(-σ_transit / σ_min) / denom

    n = length(agents)
    infection_curve = Float64[]
    q_curves = Dict(q => Float64[] for q in 1:5)
    td_curves = Dict(k => Float64[] for k in [:high, :low])

    tds = sort([a.transit_dependency for a in agents])
    td_high_thresh = tds[max(1, round(Int, 0.75 * n))]
    td_low_thresh = tds[max(1, round(Int, 0.25 * n))]

    total_transit_choices = 0
    total_choices = 0

    for t in 1:N_COMMUTES
        x_t = count(a -> a.infected, agents) / n

        for a in agents
            td_bonus = a.transit_dependency * 0.4
            α_transit_agent = clamp(α_transit_base + td_bonus, 0.01, 0.99)
            α_walk_agent = 1.0 - α_transit_agent

            chose_transit = rand() < α_transit_agent

            if a.infected
                total_transit_choices += chose_transit ? 1 : 0
                total_choices += 1
                continue
            end

            if chose_transit
                total_transit_choices += 1
                p_inf = 1.0 - exp(-B2 * LAMBDA2 * α_transit_agent * x_t)
                if rand() < p_inf
                    a.infected = true
                    a.infected_in = :transit
                end
            else
                p_inf = 1.0 - exp(-B1 * LAMBDA1 * α_walk_agent * x_t)
                if rand() < p_inf
                    a.infected = true
                    a.infected_in = :street
                end
            end
            total_choices += 1
        end

        frac = count(a -> a.infected, agents) / n
        push!(infection_curve, frac)

        for q in 1:5
            qa = filter(a -> a.income_quantile == q, agents)
            push!(q_curves[q], isempty(qa) ? 0.0 : count(a -> a.infected, qa) / length(qa))
        end

        high_td = filter(a -> a.transit_dependency >= td_high_thresh, agents)
        low_td = filter(a -> a.transit_dependency <= td_low_thresh, agents)
        push!(td_curves[:high], isempty(high_td) ? 0.0 : count(a -> a.infected, high_td) / length(high_td))
        push!(td_curves[:low], isempty(low_td) ? 0.0 : count(a -> a.infected, low_td) / length(low_td))

        frac >= STOP_FRAC && break
    end

    pct_transit = total_choices > 0 ? 100.0 * total_transit_choices / total_choices : 0.0

    return (infection_curve=infection_curve, pct_transit=pct_transit,
            final_frac=count(a -> a.infected, agents) / n,
            α_transit=α_transit_base, α_walk=α_walk_base,
            q_curves=q_curves, td_curves=td_curves)
end

freq_to_sigma(freq_min) = 2.0 * freq_min

frequencies = Float64.(1:25)
n_runs = 50

summary_rows = DataFrame(
    freq=Float64[], sigma_transit=Float64[], run=Int[],
    pct_transit=Float64[], final_infected=Float64[],
    q1_infected=Float64[], q5_infected=Float64[],
    high_td_infected=Float64[], low_td_infected=Float64[],
    alpha_transit=Float64[],
    time_to_50=Int[], time_to_90=Int[], time_to_99=Int[],
    q1_time_to_50=Int[], q5_time_to_50=Int[],
    high_td_time_to_50=Int[], low_td_time_to_50=Int[])

curve_rows = DataFrame(freq=Float64[], commute=Int[],
    frac_infected=Float64[], q1=Float64[], q5=Float64[],
    high_td=Float64[], low_td=Float64[])

function time_to(curve, thresh)
    idx = findfirst(>=(thresh), curve)
    return isnothing(idx) ? length(curve) : idx
end

for freq in frequencies
    σ_t = freq_to_sigma(freq)
    @printf("Freq %2.0f min (σ_t=%4.1f): ", freq, σ_t)
    flush(stdout)

    all_curves = Vector{Vector{Float64}}()
    all_q1 = Vector{Vector{Float64}}()
    all_q5 = Vector{Vector{Float64}}()
    all_ht = Vector{Vector{Float64}}()
    all_lt = Vector{Vector{Float64}}()

    for run in 1:n_runs
        r = run_one(σ_t; seed=42000 + run*100 + round(Int, freq))
        push!(all_curves, r.infection_curve)
        push!(all_q1, r.q_curves[1])
        push!(all_q5, r.q_curves[5])
        push!(all_ht, r.td_curves[:high])
        push!(all_lt, r.td_curves[:low])

        t100 = min(100, length(r.infection_curve))

        push!(summary_rows, (
            freq=freq, sigma_transit=σ_t, run=run,
            pct_transit=r.pct_transit,
            final_infected=r.final_frac,
            q1_infected=r.q_curves[1][t100],
            q5_infected=r.q_curves[5][t100],
            high_td_infected=r.td_curves[:high][t100],
            low_td_infected=r.td_curves[:low][t100],
            alpha_transit=r.α_transit,
            time_to_50=time_to(r.infection_curve, 0.5),
            time_to_90=time_to(r.infection_curve, 0.9),
            time_to_99=time_to(r.infection_curve, 0.99),
            q1_time_to_50=time_to(r.q_curves[1], 0.5),
            q5_time_to_50=time_to(r.q_curves[5], 0.5),
            high_td_time_to_50=time_to(r.td_curves[:high], 0.5),
            low_td_time_to_50=time_to(r.td_curves[:low], 0.5)))
    end

    function pad_avg(curves)
        max_len = maximum(length.(curves))
        mat = zeros(max_len, length(curves))
        for (j, c) in enumerate(curves)
            for i in eachindex(1:max_len)
                mat[i, j] = i <= length(c) ? c[i] : c[end]
            end
        end
        vec(mean(mat, dims=2))
    end

    avg = pad_avg(all_curves)
    aq1 = pad_avg(all_q1)
    aq5 = pad_avg(all_q5)
    aht = pad_avg(all_ht)
    alt = pad_avg(all_lt)

    for i in eachindex(avg)
        push!(curve_rows, (freq=freq, commute=i,
            frac_infected=avg[i], q1=aq1[i], q5=aq5[i],
            high_td=aht[i], low_td=alt[i]))
    end

    fd = filter(r -> r.freq == freq, summary_rows)
    @printf("α_t=%.3f  transit=%.0f%%  t99=%.0f  HTD_t50=%.0f  LTD_t50=%.0f  gap=%.1f\n",
        mean(fd.alpha_transit), mean(fd.pct_transit),
        mean(fd.time_to_99),
        mean(fd.high_td_time_to_50), mean(fd.low_td_time_to_50),
        mean(fd.low_td_time_to_50) - mean(fd.high_td_time_to_50))
    flush(stdout)
end

println("\n--- Analytical verification (MacRury Eq 3 approximation) ---")
println("σ_transit | α_walk | α_transit |    A    | tf≈1/A*ln(...)")
for freq in frequencies
    σ_t = freq_to_sigma(freq)
    σ_min = min(σ_WALK, σ_t)
    denom = exp(-σ_WALK / σ_min) + exp(-σ_t / σ_min)
    α1 = exp(-σ_WALK / σ_min) / denom
    α2 = exp(-σ_t / σ_min) / denom
    A = B1 * LAMBDA1 * α1^2 + B2 * LAMBDA2 * α2^2
    tf = (1/A) * log(0.99/0.01 * 0.99/0.01)
    @printf("  %5.1f   | %.3f  |  %.3f    | %.5f | %.0f\n", σ_t, α1, α2, A, tf)
end
flush(stdout)

function write_csv(path, df)
    open(path, "w") do io
        println(io, join(names(df), ","))
        for row in eachrow(df)
            println(io, join([row[c] for c in names(df)], ","))
        end
    end
end

write_csv(joinpath(OUTPUT_DIR, "summary.csv"), summary_rows)
write_csv(joinpath(OUTPUT_DIR, "curves.csv"), curve_rows)

tract_rows = DataFrame(
    location_id=tracts.location_id,
    lat=tracts.latitude, lon=tracts.longitude,
    pop_density=coalesce.(tracts.pop_density, 0.0),
    pct_no_vehicle=pct_nv,
    median_income=med_inc)
write_csv(joinpath(OUTPUT_DIR, "tracts_sdoh.csv"), tract_rows)

println("\n" * "="^70)
println("Done! Results in: $OUTPUT_DIR")
println("Run: julia --project=. plot_results.jl")
flush(stdout)
