# JSON Occultation Record - Complete Specification

## Overview

This document describes the complete JSON structure for an asteroid occultation event as computed by IOccultCalc.

**File Example:** `occultation_complete_example.json`

---

## Root Structure

```json
{
  "id": <number>,
  "event_type": <string>,
  "asteroid": {...},
  "star": {...},
  "event": {...},
  "location": {...},
  "prediction_quality": {...},
  "observing_recommendations": {...},
  "scientific_value": {...},
  "computation_metadata": {...},
  "visualization_data": {...},
  "external_links": {...}
}
```

---

## 1. Root Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | number | Unique event identifier | `1` |
| `event_type` | string | Type of event | `"asteroid_occultation"` |

---

## 2. Asteroid Object

Complete information about the occulting asteroid.

```json
"asteroid": {
  "id": <number>,
  "number": <number>,
  "name": <string>,
  "designation": <string>,
  "magnitude": {...},
  "diameter_km": <number>,
  "albedo": <number>,
  "orbital_elements": {...},
  "physical_properties": {...}
}
```

### 2.1 Asteroid Fields

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `id` | number | - | Asteroid number (MPC) |
| `number` | number | - | Same as id |
| `name` | string | - | Asteroid name |
| `designation` | string | - | Full designation with number |
| `diameter_km` | number | km | Asteroid diameter |
| `albedo` | number | - | Geometric albedo (0-1) |

### 2.2 Magnitude Object

```json
"magnitude": {
  "absolute_h": <number>,
  "apparent": <number>,
  "v_band": <number>
}
```

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `absolute_h` | number | mag | Absolute magnitude H |
| `apparent` | number | mag | Apparent magnitude at event time |
| `v_band` | number | mag | V-band magnitude |

### 2.3 Orbital Elements Object

```json
"orbital_elements": {
  "semi_major_axis_au": <number>,
  "eccentricity": <number>,
  "inclination_deg": <number>,
  "longitude_ascending_node_deg": <number>,
  "argument_perihelion_deg": <number>,
  "mean_anomaly_deg": <number>,
  "epoch_jd": <number>,
  "elements_source": <string>,
  "osculating": <boolean>
}
```

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `semi_major_axis_au` | number | AU | Semi-major axis |
| `eccentricity` | number | - | Orbital eccentricity |
| `inclination_deg` | number | deg | Orbital inclination |
| `longitude_ascending_node_deg` | number | deg | Longitude of ascending node (Ω) |
| `argument_perihelion_deg` | number | deg | Argument of perihelion (ω) |
| `mean_anomaly_deg` | number | deg | Mean anomaly (M) |
| `epoch_jd` | number | JD | Epoch of elements |
| `elements_source` | string | - | Source of orbital elements |
| `osculating` | boolean | - | True if osculating elements |

### 2.4 Physical Properties Object

```json
"physical_properties": {
  "rotation_period_hours": <number>,
  "spectral_type": <string>,
  "taxonomy": <string>,
  "shape_model_available": <boolean>
}
```

---

## 3. Star Object

Complete information about the occulted star.

```json
"star": {
  "catalog": <string>,
  "catalog_id": <string>,
  "name": <string>,
  "alternative_designations": [<string>, ...],
  "coordinates": {...},
  "photometry": {...},
  "spectral_info": {...}
}
```

### 3.1 Star Fields

| Field | Type | Description |
|-------|------|-------------|
| `catalog` | string | Star catalog name (e.g., "Gaia DR3") |
| `catalog_id` | string | Full catalog identifier |
| `name` | string | Common name or designation |
| `alternative_designations` | array | Other catalog names |

### 3.2 Coordinates Object

```json
"coordinates": {
  "ra_deg": <number>,
  "dec_deg": <number>,
  "ra_hms": <string>,
  "dec_dms": <string>,
  "epoch": <string>,
  "proper_motion": {...},
  "parallax": {...}
}
```

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `ra_deg` | number | deg | Right Ascension in degrees |
| `dec_deg` | number | deg | Declination in degrees |
| `ra_hms` | string | - | RA in hours:minutes:seconds format |
| `dec_dms` | string | - | Dec in degrees:arcmin:arcsec format |
| `epoch` | string | - | Coordinate epoch (e.g., "J2016.0") |

#### 3.2.1 Proper Motion Object

```json
"proper_motion": {
  "pmra_mas_yr": <number>,
  "pmdec_mas_yr": <number>,
  "has_proper_motion": <boolean>
}
```

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `pmra_mas_yr` | number | mas/yr | Proper motion in RA |
| `pmdec_mas_yr` | number | mas/yr | Proper motion in Dec |
| `has_proper_motion` | boolean | - | True if proper motion available |

#### 3.2.2 Parallax Object

```json
"parallax": {
  "value_mas": <number>,
  "error_mas": <number>,
  "distance_pc": <number>
}
```

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `value_mas` | number | mas | Parallax value |
| `error_mas` | number | mas | Parallax error |
| `distance_pc` | number | pc | Computed distance |

### 3.3 Photometry Object

```json
"photometry": {
  "magnitude_g": <number>,
  "magnitude_v": <number>,
  "magnitude_r": <number>,
  "magnitude_i": <number>,
  "magnitude_j": <number>,
  "magnitude_h": <number>,
  "magnitude_k": <number>,
  "color_index_bv": <number>,
  "color_index_vi": <number>
}
```

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `magnitude_g` | number | mag | Gaia G-band magnitude |
| `magnitude_v` | number | mag | Johnson V magnitude |
| `magnitude_r` | number | mag | R-band magnitude |
| `magnitude_i` | number | mag | I-band magnitude |
| `magnitude_j` | number | mag | 2MASS J magnitude |
| `magnitude_h` | number | mag | 2MASS H magnitude |
| `magnitude_k` | number | mag | 2MASS K magnitude |
| `color_index_bv` | number | mag | B-V color index |
| `color_index_vi` | number | mag | V-I color index |

---

## 4. Event Object

Complete information about the occultation event.

```json
"event": {
  "date_time": {...},
  "circumstances": {...},
  "geometry": {...},
  "path_information": {...},
  "observing_circumstances": {...}
}
```

### 4.1 Date Time Object

```json
"date_time": {
  "jd": <number>,
  "mjd": <number>,
  "iso8601": <string>,
  "calendar": {...},
  "gregorian": <string>
}
```

| Field | Type | Format | Description |
|-------|------|--------|-------------|
| `jd` | number | JD | Julian Date |
| `mjd` | number | MJD | Modified Julian Date |
| `iso8601` | string | ISO 8601 | ISO date-time string |
| `gregorian` | string | - | Human-readable date-time |

#### 4.1.1 Calendar Object

```json
"calendar": {
  "year": <number>,
  "month": <number>,
  "day": <number>,
  "hour": <number>,
  "minute": <number>,
  "second": <number>
}
```

### 4.2 Circumstances Object

```json
"circumstances": {
  "duration_seconds": <number>,
  "duration_uncertainty_seconds": <number>,
  "magnitude_drop": <number>,
  "magnitude_drop_uncertainty": <number>,
  "combined_magnitude": <number>,
  "snr_estimate": <number>,
  "observability_rating": <string>
}
```

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `duration_seconds` | number | s | Event duration |
| `duration_uncertainty_seconds` | number | s | Duration uncertainty (1σ) |
| `magnitude_drop` | number | mag | Magnitude drop during event |
| `magnitude_drop_uncertainty` | number | mag | Drop uncertainty |
| `combined_magnitude` | number | mag | Asteroid+star combined magnitude |
| `snr_estimate` | number | - | Signal-to-noise ratio estimate |
| `observability_rating` | string | - | Rating: excellent/good/moderate/poor |

### 4.3 Geometry Object

```json
"geometry": {
  "closest_approach": {...},
  "shadow_properties": {...},
  "sky_motion": {...}
}
```

#### 4.3.1 Closest Approach Object

```json
"closest_approach": {
  "distance_km": <number>,
  "uncertainty_1sigma_km": <number>,
  "uncertainty_cross_track_km": <number>,
  "uncertainty_along_track_km": <number>
}
```

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `distance_km` | number | km | Closest approach distance (usually 0) |
| `uncertainty_1sigma_km` | number | km | Total uncertainty (1σ) |
| `uncertainty_cross_track_km` | number | km | Cross-track uncertainty |
| `uncertainty_along_track_km` | number | km | Along-track uncertainty |

#### 4.3.2 Shadow Properties Object

```json
"shadow_properties": {
  "shadow_velocity_km_s": <number>,
  "shadow_width_km": <number>,
  "position_angle_deg": <number>,
  "position_angle_uncertainty_deg": <number>
}
```

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `shadow_velocity_km_s` | number | km/s | Shadow velocity on Earth |
| `shadow_width_km` | number | km | Shadow width (= asteroid diameter) |
| `position_angle_deg` | number | deg | Position angle of shadow motion |
| `position_angle_uncertainty_deg` | number | deg | PA uncertainty |

### 4.4 Path Information Object

**CRITICAL FOR MAP VISUALIZATION**

```json
"path_information": {
  "path_center": {...},
  "path_width_km": <number>,
  "path_uncertainty_1sigma_km": <number>,
  "northern_limit": {...},
  "southern_limit": {...},
  "path_coordinates": [...]
}
```

| Field | Type | Description |
|-------|------|-------------|
| `path_width_km` | number | Width of occultation path |
| `path_uncertainty_1sigma_km` | number | Path uncertainty (1σ) |

#### 4.4.1 Path Center Object

```json
"path_center": {
  "latitude_deg": <number>,
  "longitude_deg": <number>,
  "location_name": <string>
}
```

#### 4.4.2 Path Coordinates Array

**For drawing the path on a map:**

```json
"path_coordinates": [
  {"lat": <number>, "lon": <number>},
  {"lat": <number>, "lon": <number>},
  ...
]
```

This array contains points along the centerline of the occultation path.

---

## 5. Location Object

Observer location information.

```json
"location": {
  "name": <string>,
  "code": <string>,
  "country": <string>,
  "region": <string>,
  "coordinates": {...},
  "timezone": <string>,
  "local_time": <string>
}
```

### 5.1 Coordinates Object

```json
"coordinates": {
  "latitude_deg": <number>,
  "longitude_deg": <number>,
  "elevation_m": <number>,
  "geodetic_datum": <string>
}
```

---

## 6. Visualization Data Object

**ESSENTIAL FOR API/UI VISUALIZATION**

```json
"visualization_data": {
  "map": {...},
  "sky_chart": {...}
}
```

### 6.1 Map Object

**For geographic map visualization:**

```json
"map": {
  "center": {"lat": <number>, "lon": <number>},
  "zoom_level": <number>,
  "path_geometry": <string>,
  "shadow_polygon": [[lat, lon], ...]
}
```

| Field | Type | Description |
|-------|------|-------------|
| `center` | object | Map center coordinates |
| `zoom_level` | number | Suggested zoom level (1-18) |
| `path_geometry` | string | GeoJSON geometry type |
| `shadow_polygon` | array | Polygon coordinates for shadow area |

### 6.2 Sky Chart Object

**For celestial map visualization:**

```json
"sky_chart": {
  "ra_center_deg": <number>,
  "dec_center_deg": <number>,
  "field_of_view_arcmin": <number>,
  "magnitude_limit": <number>,
  "show_trajectory": <boolean>,
  "trajectory_points": [...]
}
```

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `ra_center_deg` | number | deg | RA of chart center |
| `dec_center_deg` | number | deg | Dec of chart center |
| `field_of_view_arcmin` | number | arcmin | Field of view |
| `magnitude_limit` | number | mag | Faintest stars to show |
| `show_trajectory` | boolean | - | Show asteroid path |

#### 6.2.1 Trajectory Points Array

```json
"trajectory_points": [
  {"ra": <number>, "dec": <number>, "time": <string>},
  ...
]
```

---

## 7. Additional Metadata

### 7.1 Prediction Quality

```json
"prediction_quality": {
  "confidence_level": <number>,
  "prediction_uncertainty_class": <string>,
  "orbit_quality_code": <string>,
  "last_observation_days_ago": <number>,
  "number_of_observations": <number>,
  "observation_arc_years": <number>,
  "elements_age_days": <number>,
  "star_position_uncertainty_mas": <number>,
  "ephemeris_uncertainty_arcsec": <number>
}
```

### 7.2 Observing Recommendations

```json
"observing_recommendations": {
  "priority": <string>,
  "difficulty": <string>,
  "recommended_equipment": {...},
  "preparation": {...},
  "notes": [<string>, ...]
}
```

### 7.3 Computation Metadata

```json
"computation_metadata": {
  "computed_at": <string>,
  "computation_time_ms": <number>,
  "integration_steps": <number>,
  "convergence_achieved": <boolean>,
  "propagator": <string>,
  "perturbations_included": [<string>, ...],
  "reference_frame": <string>
}
```

---

## API Usage Examples

### Example 1: Drawing Geographic Map

```javascript
// Extract path data
const pathCoords = occultation.event.path_information.path_coordinates;
const center = occultation.visualization_data.map.center;
const shadowPolygon = occultation.visualization_data.map.shadow_polygon;

// Draw on Leaflet/OpenLayers/Google Maps
map.setCenter([center.lat, center.lon]);
map.addPolyline(pathCoords, {color: 'red', width: 2});
map.addPolygon(shadowPolygon, {fillColor: 'rgba(255,0,0,0.2)'});
```

### Example 2: Drawing Sky Chart

```javascript
// Extract star position and trajectory
const star = {
  ra: occultation.star.coordinates.ra_deg,
  dec: occultation.star.coordinates.dec_deg,
  mag: occultation.star.photometry.magnitude_v
};

const trajectory = occultation.visualization_data.sky_chart.trajectory_points;

// Plot on celestial chart
skyChart.plotStar(star.ra, star.dec, star.mag);
skyChart.plotTrajectory(trajectory);
```

### Example 3: Event Time Display

```javascript
// Multiple time formats available
const jd = occultation.event.date_time.jd;
const iso = occultation.event.date_time.iso8601;
const readable = occultation.event.date_time.gregorian;
const local = occultation.location.local_time;

console.log(`Event at: ${readable} (${local} local time)`);
```

---

## Data Type Reference

| Type | JSON Type | Example | Description |
|------|-----------|---------|-------------|
| `<number>` | number | `12.34` | Numeric value |
| `<string>` | string | `"Eros"` | Text string |
| `<boolean>` | boolean | `true` | True or false |
| `<array>` | array | `[1, 2, 3]` | Array of values |
| `<object>` | object | `{"key": "value"}` | Key-value pairs |

---

## Units Reference

| Unit | Description | Example |
|------|-------------|---------|
| `deg` | Degrees | 45.0 |
| `arcmin` | Arc minutes | 30.0 |
| `arcsec` | Arc seconds | 15.0 |
| `mas` | Milliarcseconds | 2.34 |
| `mas/yr` | Milliarcseconds per year | -3.45 |
| `km` | Kilometers | 16.84 |
| `km/s` | Kilometers per second | 18.42 |
| `m` | Meters | 1165 |
| `AU` | Astronomical Units | 1.458 |
| `pc` | Parsecs | 427.4 |
| `mag` | Magnitudes | 14.23 |
| `s` | Seconds | 8.34 |
| `JD` | Julian Date | 2460315.789456 |
| `MJD` | Modified Julian Date | 60315.289456 |

---

## Schema Validation

For JSON schema validation, see: `occultation_schema.json`

---

## Version History

- **v1.0** (2025-12-02): Initial complete specification
- Added visualization_data for API support
- Added path_coordinates for map rendering
- Added trajectory_points for sky chart rendering

---

## Related Files

- `occultation_complete_example.json` - Complete example
- `occultation_schema.json` - JSON Schema for validation
- `OUTPUT_FORMATTER_DOCUMENTATION.md` - Output formatter documentation

---

## Contact

For questions or suggestions about this specification:
- GitHub: https://github.com/manvalan/IOoccultCalc
- Issues: https://github.com/manvalan/IOoccultCalc/issues
