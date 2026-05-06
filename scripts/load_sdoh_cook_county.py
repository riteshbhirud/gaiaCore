import requests
import psycopg2
import sys
import time

DB_HOST = "localhost"
DB_PORT = 5433
DB_NAME = "gaiacore"
DB_USER = "postgres"
DB_PASS = "postgres"

STATE_FIPS = "17"
COUNTY_FIPS = "031"
COUNTY_NAME = "Cook"
STATE_NAME = "IL"

TIGERWEB_URL = (
    "https://tigerweb.geo.census.gov/arcgis/rest/services/"
    "TIGERweb/tigerWMS_ACS2022/MapServer/6/query"
)
ACS_DETAIL_URL = "https://api.census.gov/data/2022/acs/acs5"
ACS_PROFILE_URL = "https://api.census.gov/data/2022/acs/acs5/profile"


def fetch_tract_centroids():
    print("Fetching tract centroids from TIGERweb...")
    tracts = []
    offset = 0
    batch_size = 1000

    while True:
        params = {
            "where": f"STATE='{STATE_FIPS}' AND COUNTY='{COUNTY_FIPS}'",
            "outFields": "GEOID,NAME,CENTLAT,CENTLON,AREALAND",
            "returnGeometry": "false",
            "f": "json",
            "resultRecordCount": batch_size,
            "resultOffset": offset,
        }
        resp = requests.get(TIGERWEB_URL, params=params, timeout=60)
        resp.raise_for_status()
        data = resp.json()

        features = data.get("features", [])
        if not features:
            break

        for f in features:
            a = f["attributes"]
            tracts.append({
                "geoid": a["GEOID"],
                "name": a["NAME"],
                "lat": float(a["CENTLAT"]),
                "lon": float(a["CENTLON"]),
                "arealand_m2": int(a["AREALAND"]),
            })

        if len(features) < batch_size:
            break
        offset += batch_size
        time.sleep(0.5)

    print(f"  Got {len(tracts)} tracts")
    return {t["geoid"]: t for t in tracts}


def fetch_acs_detail(tracts):
    print("Fetching ACS detailed tables (income, vehicles, population)...")
    params = {
        "get": "NAME,B19013_001E,B08141_001E,B08141_002E,B01003_001E",
        "for": "tract:*",
        "in": f"state:{STATE_FIPS} county:{COUNTY_FIPS}",
    }
    resp = requests.get(ACS_DETAIL_URL, params=params, timeout=60)
    resp.raise_for_status()
    rows = resp.json()
    header = rows[0]
    print(f"  Got {len(rows) - 1} rows")

    for row in rows[1:]:
        rec = dict(zip(header, row))
        geoid = rec["state"] + rec["county"] + rec["tract"]
        if geoid not in tracts:
            continue
        t = tracts[geoid]
        t["median_income"] = _to_float(rec["B19013_001E"])
        total_workers = _to_float(rec["B08141_001E"])
        no_vehicle = _to_float(rec["B08141_002E"])
        t["pct_no_vehicle"] = (
            round(100.0 * no_vehicle / total_workers, 2)
            if total_workers and total_workers > 0
            else None
        )
        pop = _to_float(rec["B01003_001E"])
        t["population"] = pop
        area_km2 = t["arealand_m2"] / 1_000_000 if t["arealand_m2"] > 0 else None
        t["pop_density"] = (
            round(pop / area_km2, 2) if pop and area_km2 else None
        )


def fetch_acs_profile(tracts):
    print("Fetching ACS data profiles (percent uninsured)...")
    params = {
        "get": "DP03_0099PE",
        "for": "tract:*",
        "in": f"state:{STATE_FIPS} county:{COUNTY_FIPS}",
    }
    resp = requests.get(ACS_PROFILE_URL, params=params, timeout=60)
    resp.raise_for_status()
    rows = resp.json()
    header = rows[0]
    print(f"  Got {len(rows) - 1} rows")

    for row in rows[1:]:
        rec = dict(zip(header, row))
        geoid = rec["state"] + rec["county"] + rec["tract"]
        if geoid not in tracts:
            continue
        tracts[geoid]["pct_uninsured"] = _to_float(rec["DP03_0099PE"])


def _to_float(val):
    if val is None or val == "" or val == "null":
        return None
    try:
        v = float(val)
        return v if v >= 0 else None
    except (ValueError, TypeError):
        return None


def load_into_gaiacore(tracts):
    conn = psycopg2.connect(
        host=DB_HOST, port=DB_PORT, dbname=DB_NAME,
        user=DB_USER, password=DB_PASS,
    )
    conn.autocommit = False
    cur = conn.cursor()

    valid = [
        t for t in tracts.values()
        if t.get("lat") and t.get("lon") and t.get("population") is not None
    ]
    print(f"\nLoading {len(valid)} tracts into gaiaCore...")

    location_ids = {}
    for t in valid:
        cur.execute(
            """
            INSERT INTO working.location
                (address_1, city, state, zip, county,
                 location_source_value, country_source_value,
                 latitude, longitude, geom)
            VALUES
                (%s, %s, %s, %s, %s,
                 %s, %s,
                 %s, %s, ST_SetSRID(ST_MakePoint(%s, %s), 4326))
            RETURNING location_id
            """,
            (
                t["name"],
                "Chicago",
                STATE_NAME,
                None,
                COUNTY_NAME,
                t["geoid"],
                "UNITED STATES OF AMERICA",
                t["lat"],
                t["lon"],
                t["lon"],
                t["lat"],
            ),
        )
        location_ids[t["geoid"]] = cur.fetchone()[0]

    print(f"  Inserted {len(location_ids)} location records")

    sdoh_measures = [
        ("median_income", "Median Household Income (USD)", "USD"),
        ("pct_no_vehicle", "Percent Transit-Dependent (No Vehicle)", "percent"),
        ("pct_uninsured", "Percent Uninsured", "percent"),
        ("pop_density", "Population Density (per sq km)", "per_km2"),
    ]

    exposure_count = 0
    for t in valid:
        loc_id = location_ids.get(t["geoid"])
        if loc_id is None:
            continue
        for field, label, unit in sdoh_measures:
            value = t.get(field)
            if value is None:
                continue
            cur.execute(
                """
                INSERT INTO working.external_exposure
                    (location_id, exposure_source_value,
                     dose_unit_source_value, value_as_number,
                     modifier_source_value)
                VALUES (%s, %s, %s, %s, %s)
                """,
                (loc_id, label, unit, value, "ACS 5-Year 2022"),
            )
            exposure_count += 1

    conn.commit()
    print(f"  Inserted {exposure_count} external_exposure records")
    cur.close()
    conn.close()


def main():
    print("=" * 60)
    print("Cook County IL - Census Tract SDoH Loader")
    print("=" * 60)
    print()

    tracts = fetch_tract_centroids()
    fetch_acs_detail(tracts)
    fetch_acs_profile(tracts)

    has_income = sum(1 for t in tracts.values() if t.get("median_income") is not None)
    has_vehicle = sum(1 for t in tracts.values() if t.get("pct_no_vehicle") is not None)
    has_unins = sum(1 for t in tracts.values() if t.get("pct_uninsured") is not None)
    has_density = sum(1 for t in tracts.values() if t.get("pop_density") is not None)
    print(f"\nSDoH coverage:")
    print(f"  Median income:      {has_income}/{len(tracts)}")
    print(f"  No vehicle %:       {has_vehicle}/{len(tracts)}")
    print(f"  Uninsured %:        {has_unins}/{len(tracts)}")
    print(f"  Population density: {has_density}/{len(tracts)}")

    load_into_gaiacore(tracts)

    print("\nDone! Verify with:")
    print(f"  curl -H 'Accept-Profile: working' 'http://localhost:3000/location?limit=5&county=eq.{COUNTY_NAME}'")
    print(f"  curl -H 'Accept-Profile: working' 'http://localhost:3000/external_exposure?limit=10'")


if __name__ == "__main__":
    main()
