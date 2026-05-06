# GaiaDB.jl

Julia client for the gaiaCore PostgREST API. Returns query results as DataFrames.

## Setup

```julia
using Pkg
Pkg.activate("connectors/GaiaDB.jl")
Pkg.instantiate()
```

## Usage

```julia
using GaiaDB

# List all locations (returns a DataFrame)
df = list_locations("http://localhost:3000")

# Filter by city
df = list_locations("http://localhost:3000"; city="FRESNO", limit=5)

# Get data sources
ds = get_data_sources("http://localhost:3000")

# Get exposure records
ex = get_exposures("http://localhost:3000"; person_id=1)
```

## Running tests

Requires a running gaiaCore instance at `http://localhost:3000` (override with `GAIADB_URL`).

```bash
cd connectors/GaiaDB.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```
