using Pkg
Pkg.activate(@__DIR__)

using CairoMakie
using DataFrames
using Statistics

const OUTPUT_DIR = joinpath(@__DIR__, "output")

function read_csv(path)
    lines = readlines(path)
    header = Symbol.(split(lines[1], ","))
    data = [split(l, ",") for l in lines[2:end]]
    isempty(data) && return DataFrame()
    df = DataFrame()
    for (i, col) in enumerate(header)
        vals = [tryparse(Float64, row[i]) for row in data]
        if all(!isnothing, vals)
            df[!, col] = Float64.(vals)
        else
            df[!, col] = [row[i] for row in data]
        end
    end
    return df
end

summary_df = read_csv(joinpath(OUTPUT_DIR, "summary.csv"))
curves_df = read_csv(joinpath(OUTPUT_DIR, "curves.csv"))
tracts_df = read_csv(joinpath(OUTPUT_DIR, "tracts_sdoh.csv"))

frequencies = sort(unique(summary_df.freq))
println("Loaded data for $(length(frequencies)) frequencies: $(Int.(frequencies[1])):$(Int.(frequencies[end]))")

freq_stats = DataFrame(
    freq = frequencies,
    sigma_transit = [mean(filter(r -> r.freq == f, summary_df).sigma_transit) for f in frequencies],
    pct_transit = [mean(filter(r -> r.freq == f, summary_df).pct_transit) for f in frequencies],
    alpha_transit = [mean(filter(r -> r.freq == f, summary_df).alpha_transit) for f in frequencies],
    time_to_99 = [mean(filter(r -> r.freq == f, summary_df).time_to_99) for f in frequencies],
    time_to_50 = [mean(filter(r -> r.freq == f, summary_df).time_to_50) for f in frequencies],
    inf_at_100 = [100*mean(filter(r -> r.freq == f, summary_df).final_infected) for f in frequencies],
    q1_t50 = [mean(filter(r -> r.freq == f, summary_df).q1_time_to_50) for f in frequencies],
    q5_t50 = [mean(filter(r -> r.freq == f, summary_df).q5_time_to_50) for f in frequencies],
    htd_t50 = [mean(filter(r -> r.freq == f, summary_df).high_td_time_to_50) for f in frequencies],
    ltd_t50 = [mean(filter(r -> r.freq == f, summary_df).low_td_time_to_50) for f in frequencies],
)

println("Generating Figure 1: Overall infection dynamics...")

fig1 = Figure(size=(1200, 500), fontsize=14)

ax1 = Axis(fig1[1, 1],
    xlabel = "Commute iteration",
    ylabel = "Fraction infected",
    title = "Infection Dynamics by Transit Frequency")

sel_freqs = [1.0, 5.0, 8.0, 15.0, 25.0]
colors = cgrad(:RdYlBu_5, 5, categorical=true)
for (i, freq) in enumerate(sel_freqs)
    fc = filter(r -> r.freq == freq, curves_df)
    lines!(ax1, fc.commute, fc.frac_infected,
           label = "$(Int(freq)) min (σ₂=$(Int(freq*2)))",
           color = colors[i], linewidth = 2.5)
end
axislegend(ax1, position = :rb, framevisible=true)

ax2a = Axis(fig1[1, 2],
    xlabel = "Transit headway (min)",
    ylabel = "Iterations to 99% infected",
    title = "MacRury Figure 4 Reproduction",
    yticklabelcolor = :red)

ax2b = Axis(fig1[1, 2],
    ylabel = "% infected after 100 iterations",
    yticklabelcolor = :blue,
    yaxisposition = :right)
hidespines!(ax2b)
hidexdecorations!(ax2b)

inf_at_100_pct = Float64[]
for freq in frequencies
    fc = filter(r -> r.freq == freq, curves_df)
    t100 = min(100, nrow(fc))
    push!(inf_at_100_pct, 100.0 * fc.frac_infected[t100])
end

scatterlines!(ax2a, frequencies, freq_stats.time_to_99,
    color = :red, linewidth = 2.5, markersize = 8)
scatterlines!(ax2b, frequencies, inf_at_100_pct,
    color = :blue, linewidth = 2.5, markersize = 8, marker = :diamond)

save(joinpath(OUTPUT_DIR, "figure1_overall_infection.pdf"), fig1)
println("  → figure1_overall_infection.pdf")

println("Generating Figure 2: Stratified infection curves...")

disparity_gaps = freq_stats.ltd_t50 .- freq_stats.htd_t50
peak_disparity_idx = argmax(disparity_gaps)
peak_freq = frequencies[peak_disparity_idx]
println("  Max disparity frequency: $(Int(peak_freq)) min (gap=$(round(disparity_gaps[peak_disparity_idx], digits=1)) iters)")

fig2 = Figure(size=(1200, 500), fontsize=14)

fc_peak = filter(r -> r.freq == peak_freq, curves_df)

ax3 = Axis(fig2[1, 1],
    xlabel = "Commute iteration",
    ylabel = "Fraction infected",
    title = "Infection by Income Quintile ($(Int(peak_freq)) min headway)")

lines!(ax3, fc_peak.commute, fc_peak.q1,
       label = "Q1 — lowest income", color = :darkred, linewidth = 2.5)
lines!(ax3, fc_peak.commute, fc_peak.q5,
       label = "Q5 — highest income", color = :darkgreen, linewidth = 2.5)
lines!(ax3, fc_peak.commute, fc_peak.frac_infected,
       label = "Overall", color = :gray40, linewidth = 1.5, linestyle = :dash)
axislegend(ax3, position = :rb, framevisible=true)

ax4 = Axis(fig2[1, 2],
    xlabel = "Commute iteration",
    ylabel = "Fraction infected",
    title = "Infection by Transit Dependency ($(Int(peak_freq)) min headway)")

lines!(ax4, fc_peak.commute, fc_peak.high_td,
       label = "High transit dep. (top 25%)", color = :red, linewidth = 2.5)
lines!(ax4, fc_peak.commute, fc_peak.low_td,
       label = "Low transit dep. (bottom 25%)", color = :blue, linewidth = 2.5)
lines!(ax4, fc_peak.commute, fc_peak.frac_infected,
       label = "Overall", color = :gray40, linewidth = 1.5, linestyle = :dash)
axislegend(ax4, position = :rb, framevisible=true)

save(joinpath(OUTPUT_DIR, "figure2_stratified_infection.pdf"), fig2)
println("  → figure2_stratified_infection.pdf")

println("Generating Figure 3: Disparity across frequencies...")

fig3 = Figure(size=(1200, 500), fontsize=14)

ax5 = Axis(fig3[1, 1],
    xlabel = "Transit headway (min)",
    ylabel = "Iterations to 50% infected",
    title = "Transit Dependency Disparity")

scatterlines!(ax5, frequencies, freq_stats.htd_t50,
    color = :red, linewidth = 2.5, markersize = 10,
    label = "High transit dep. (top 25%)")
scatterlines!(ax5, frequencies, freq_stats.ltd_t50,
    color = :blue, linewidth = 2.5, markersize = 10,
    label = "Low transit dep. (bottom 25%)")
band!(ax5, frequencies, freq_stats.htd_t50, freq_stats.ltd_t50,
      color = (:red, 0.15))
axislegend(ax5, position = :rb, framevisible=true)

ax6 = Axis(fig3[1, 2],
    xlabel = "Transit headway (min)",
    ylabel = "Iterations to 50% infected",
    title = "Income Disparity")

scatterlines!(ax6, frequencies, freq_stats.q1_t50,
    color = :darkred, linewidth = 2.5, markersize = 10,
    label = "Q1 — lowest income")
scatterlines!(ax6, frequencies, freq_stats.q5_t50,
    color = :darkgreen, linewidth = 2.5, markersize = 10,
    label = "Q5 — highest income")
band!(ax6, frequencies, freq_stats.q1_t50, freq_stats.q5_t50,
      color = (:orange, 0.15))
axislegend(ax6, position = :rb, framevisible=true)

save(joinpath(OUTPUT_DIR, "figure3_disparity.pdf"), fig3)
println("  → figure3_disparity.pdf")

println("Generating Figure 4: Tract heatmap...")

fig4 = Figure(size=(1000, 750), fontsize=14)
ax7 = Axis(fig4[1, 1],
    xlabel = "Longitude", ylabel = "Latitude",
    title = "Transit Dependency by Census Tract — Cook County, IL",
    aspect = DataAspect(),
    backgroundcolor = :gray95)

pop_dens = Float64.(tracts_df.pop_density)
max_dens = maximum(pop_dens)
dot_sizes = 3.0 .+ 12.0 .* (pop_dens ./ max_dens)

sc = scatter!(ax7, tracts_df.lon, tracts_df.lat,
    color = tracts_df.pct_no_vehicle,
    colormap = :inferno,
    colorrange = (0, 60),
    markersize = dot_sizes,
    strokewidth = 0.3,
    strokecolor = :gray50)

Colorbar(fig4[1, 2], sc, label = "% Households Without Vehicle")

save(joinpath(OUTPUT_DIR, "figure4_tract_heatmap.pdf"), fig4)
println("  → figure4_tract_heatmap.pdf")

println("\n" * "="^70)
println("All figures saved to: $OUTPUT_DIR")
println("="^70)

peak_f = frequencies[argmax(inf_at_100_pct)]
sweet_f = frequencies[argmax(freq_stats.time_to_99)]
println("\nKey results:")
println("  Peak infection at iteration 100: $(Int(peak_f)) min headway ($(round(maximum(inf_at_100_pct), digits=1))%)")
println("  Slowest infection spread (sweet spot): $(Int(sweet_f)) min headway ($(round(maximum(freq_stats.time_to_99), digits=0)) iterations to 99%)")
println("  Max TD disparity at $(Int(frequencies[peak_disparity_idx])) min: high-TD reaches 50% $(round(freq_stats.htd_t50[peak_disparity_idx], digits=0)) iters vs low-TD $(round(freq_stats.ltd_t50[peak_disparity_idx], digits=0)) iters (gap=$(round(disparity_gaps[peak_disparity_idx], digits=1)))")
