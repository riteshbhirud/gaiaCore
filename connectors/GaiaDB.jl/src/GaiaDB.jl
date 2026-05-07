module GaiaDB

using HTTP
using JSON3
using DataFrames
using Statistics: median
using StatsBase: wsample

export list_locations, get_data_sources, get_exposures,
       initialize_agents_from_sdoh, EpiAgent

function _get(base_url::AbstractString, endpoint::AbstractString;
              schema::AbstractString = "backbone",
              params::Dict{String,String} = Dict{String,String}())
    url = rstrip(base_url, '/') * "/" * endpoint
    headers = Pair{String,String}[]
    if schema == "working"
        push!(headers, "Accept-Profile" => "working")
    end
    resp = HTTP.get(url, headers; query = params, status_exception = true)
    JSON3.read(resp.body)
end

function _to_dataframe(rows)
    isempty(rows) && return DataFrame()
    cols = keys(first(rows))
    DataFrame([col => [_unwrap(getproperty(r, col)) for r in rows] for col in cols])
end

_unwrap(v::JSON3.Object) = Dict(pairs(v))
_unwrap(v::JSON3.Array)  = collect(v)
_unwrap(v) = v

function list_locations(base_url::AbstractString;
                        city::Union{AbstractString,Nothing} = nothing,
                        state::Union{AbstractString,Nothing} = nothing,
                        county::Union{AbstractString,Nothing} = nothing,
                        limit::Int = 1000)
    params = Dict{String,String}("limit" => string(limit))
    !isnothing(city)   && (params["city"]   = "eq.$city")
    !isnothing(state)  && (params["state"]  = "eq.$state")
    !isnothing(county) && (params["county"] = "eq.$county")
    _to_dataframe(_get(base_url, "location"; schema = "working", params))
end

function get_data_sources(base_url::AbstractString; kwargs...)
    params = Dict{String,String}(string(k) => "eq.$v" for (k, v) in kwargs)
    _to_dataframe(_get(base_url, "data_source"; params))
end

function get_exposures(base_url::AbstractString;
                       person_id::Union{Int,Nothing} = nothing,
                       location_id::Union{Int,Nothing} = nothing,
                       limit::Int = 1000)
    params = Dict{String,String}("limit" => string(limit))
    !isnothing(person_id)   && (params["person_id"]   = "eq.$person_id")
    !isnothing(location_id) && (params["location_id"] = "eq.$location_id")
    _to_dataframe(_get(base_url, "external_exposure"; schema = "working", params))
end

Base.@kwdef mutable struct EpiAgent
    id::Int
    home_location_id::Int
    work_location_id::Int
    home_lat::Float64
    home_lon::Float64
    work_lat::Float64
    work_lon::Float64
    infection_status::Symbol = :susceptible
    transit_dependency::Float64
    income_quantile::Int
end

function _get_all(base_url::AbstractString, endpoint::AbstractString;
                  schema::AbstractString = "backbone",
                  params::Dict{String,String} = Dict{String,String}(),
                  page_size::Int = 999)
    all_rows = Any[]
    offset = 0
    while true
        p = copy(params)
        p["limit"]  = string(page_size)
        p["offset"] = string(offset)
        batch = _get(base_url, endpoint; schema, params = p)
        append!(all_rows, batch)
        length(batch) < page_size && break
        offset += page_size
    end
    return all_rows
end

function _load_tract_sdoh(base_url::AbstractString; county::AbstractString = "Cook")
    loc_rows = _get_all(base_url, "location";
                        schema = "working",
                        params = Dict{String,String}(
                            "county" => "eq.$county",
                            "select" => "location_id,latitude,longitude"))
    isempty(loc_rows) && error("No locations found for county=$county")
    locs = _to_dataframe(loc_rows)

    exp_rows = _get_all(base_url, "external_exposure";
                        schema = "working",
                        params = Dict{String,String}(
                            "select" => "location_id,exposure_source_value,value_as_number"))
    isempty(exp_rows) && error("No exposure records found")
    exps = _to_dataframe(exp_rows)

    loc_ids = Set(locs.location_id)
    exps = filter(row -> row.location_id in loc_ids, exps)

    label_to_col = Dict(
        "Population Density (per sq km)"         => :pop_density,
        "Median Household Income (USD)"          => :median_income,
        "Percent Transit-Dependent (No Vehicle)" => :pct_no_vehicle,
        "Percent Uninsured"                      => :pct_uninsured,
    )

    n = nrow(locs)
    for col in values(label_to_col)
        locs[!, col] = Vector{Union{Float64,Missing}}(missing, n)
    end
    lid_to_row = Dict(locs.location_id[i] => i for i in 1:n)

    for r in eachrow(exps)
        col = get(label_to_col, r.exposure_source_value, nothing)
        isnothing(col) && continue
        idx = get(lid_to_row, r.location_id, nothing)
        isnothing(idx) && continue
        locs[idx, col] = Float64(r.value_as_number)
    end

    return locs
end

function _income_to_quintile(incomes::AbstractVector{<:Real})
    sorted = sort(collect(skipmissing(incomes)))
    isempty(sorted) && return _ -> 3
    breaks = [sorted[max(1, ceil(Int, length(sorted) * q))] for q in (0.2, 0.4, 0.6, 0.8)]
    function classify(val)
        for (i, b) in enumerate(breaks)
            val <= b && return i
        end
        return 5
    end
    return classify
end

function initialize_agents_from_sdoh(base_url::AbstractString, n_agents::Int;
                                     county::AbstractString = "Cook")
    tracts = _load_tract_sdoh(base_url; county)

    tracts = dropmissing(tracts, :pop_density)
    tracts = filter(row -> row.pop_density > 0, tracts)
    nrow(tracts) == 0 && error("No valid tracts with population density data")

    weights = Float64.(tracts.pop_density)
    weights ./= sum(weights)

    valid_incomes = collect(skipmissing(tracts.median_income))
    quintile_fn = _income_to_quintile(valid_incomes)

    n_tracts = nrow(tracts)
    loc_ids  = tracts.location_id
    lats     = tracts.latitude
    lons     = tracts.longitude
    pct_novehicle = coalesce.(tracts.pct_no_vehicle, 0.0)
    med_incomes   = coalesce.(tracts.median_income, median(valid_incomes))

    agents = Vector{EpiAgent}(undef, n_agents)

    home_indices = wsample(1:n_tracts, weights, n_agents)
    work_indices = wsample(1:n_tracts, weights, n_agents)

    for i in 1:n_agents
        while work_indices[i] == home_indices[i]
            work_indices[i] = wsample(1:n_tracts, weights)
        end
    end

    for i in 1:n_agents
        hi = home_indices[i]
        wi = work_indices[i]

        tract_p = clamp(pct_novehicle[hi] / 100.0, 0.001, 0.999)
        κ = 20.0
        α = tract_p * κ
        β_param = (1.0 - tract_p) * κ
        td = clamp(α / (α + β_param) + randn() * sqrt(α * β_param / ((α + β_param)^2 * (α + β_param + 1))), 0.0, 1.0)

        agents[i] = EpiAgent(
            id                = i,
            home_location_id  = loc_ids[hi],
            work_location_id  = loc_ids[wi],
            home_lat          = lats[hi],
            home_lon          = lons[hi],
            work_lat          = lats[wi],
            work_lon          = lons[wi],
            transit_dependency = td,
            income_quantile   = quintile_fn(med_incomes[hi]),
        )
    end

    return agents
end

end
