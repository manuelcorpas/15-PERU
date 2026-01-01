#!/usr/bin/env python3
"""
Figure 1: Map of 28 Peruvian Populations with K=5 ADMIXTURE Ancestry
Uses Natural Earth geographic data with accurate population coordinates.

ADMIXTURE K=5 ancestry components (aligned with Figure 4):
- Altiplano Andean (Southern Andean/Aymara-associated)
- Amazonian Native
- European/Coastal
- Northern Andean (Other highlands)
- African

Population data: 1,142 individuals across 28 populations
(3 populations excluded due to insufficient sample size for ADMIXTURE)

Addresses reviewer comments:
1. "Peru" label prominently displayed
2. All sampling sites visible (overlapping points distinguished)
3. Regional delineation (Coastal, Highland, Amazonian)
4. Inset showing Peru's location in South America
5. K=5 ancestry pie charts per population
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Wedge
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import warnings

# Suppress the deprecation warning for gpd.datasets
warnings.filterwarnings('ignore', category=FutureWarning)

# ============================================================================
# CONFIGURATION
# ============================================================================
USE_PIE_CHARTS = True  # Set to False to use shaped markers (original v3 behavior)

# ============================================================================
# LOAD GEOGRAPHIC DATA
# ============================================================================

# Load Natural Earth data - with fallbacks for different geopandas versions
def load_world_data():
    """Load world country boundaries with fallbacks for different geopandas versions."""
    
    # Method 1: Try the classic geopandas.datasets (geopandas < 1.0)
    try:
        return gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    except (AttributeError, Exception):
        pass
    
    # Method 2: Try geodatasets package (geopandas >= 1.0)
    # Note: Use 'naturalearth.cities' parent or download countries directly
    try:
        import geodatasets
        # naturalearth_lowres equivalent in geodatasets
        return gpd.read_file(geodatasets.get_path('naturalearth.countries_110m'))
    except (ImportError, ValueError, Exception):
        pass
    
    # Method 3: Download Natural Earth countries directly from URL
    try:
        url = "https://naciscdn.org/naturalearth/110m/cultural/ne_110m_admin_0_countries.zip"
        gdf = gpd.read_file(url)
        # Rename columns to match expected format
        if 'NAME' in gdf.columns and 'name' not in gdf.columns:
            gdf = gdf.rename(columns={'NAME': 'name', 'CONTINENT': 'continent'})
        elif 'ADMIN' in gdf.columns and 'name' not in gdf.columns:
            gdf = gdf.rename(columns={'ADMIN': 'name', 'CONTINENT': 'continent'})
        return gdf
    except Exception as e:
        print(f"URL download failed: {e}")
        pass
    
    # Method 4: Try GitHub raw GeoJSON
    try:
        url = "https://raw.githubusercontent.com/datasets/geo-countries/master/data/countries.geojson"
        gdf = gpd.read_file(url)
        if 'ADMIN' in gdf.columns:
            gdf = gdf.rename(columns={'ADMIN': 'name'})
        return gdf
    except Exception:
        pass
    
    raise RuntimeError(
        "Could not load Natural Earth countries data.\n"
        "Please try one of these solutions:\n"
        "  1. pip install geodatasets\n"
        "  2. Download manually from: https://www.naturalearthdata.com/downloads/110m-cultural-vectors/\n"
        "     and update the script path accordingly"
    )

world = load_world_data()

# Debug: print available columns to help troubleshoot
print(f"Loaded world data with {len(world)} features")
print(f"Available columns: {list(world.columns)}")

# Handle different column naming conventions
if 'name' not in world.columns:
    if 'NAME' in world.columns:
        world = world.rename(columns={'NAME': 'name'})
    elif 'ADMIN' in world.columns:
        world = world.rename(columns={'ADMIN': 'name'})
    elif 'NAME_EN' in world.columns:
        world = world.rename(columns={'NAME_EN': 'name'})

if 'continent' not in world.columns:
    if 'CONTINENT' in world.columns:
        world = world.rename(columns={'CONTINENT': 'continent'})

# Get Peru and neighboring countries
peru = world[world['name'] == 'Peru']
south_america = world[world['continent'] == 'South America']
neighbors = world[world['name'].isin(['Ecuador', 'Colombia', 'Brazil', 'Bolivia', 'Chile'])]

# ============================================================================
# ADMIXTURE K=5 DATA (from admixture_population_means_K5_ordered.csv)
# Aligned with Figure 4 legend:
#   Ancestry_1 = Altiplano Andean (Southern Andean/Aymara-associated)
#   Ancestry_2 = Amazonian Native
#   Ancestry_3 = European/Coastal
#   Ancestry_4 = Northern Andean (Other highlands)
#   Ancestry_5 = African
# ============================================================================

admixture_data = {
    'LAMAS': {'altiplano': 0.000010, 'amazonian': 0.987690, 'european': 0.012212, 'northern_andean': 0.000078, 'african': 0.000010},
    'NAHUA': {'altiplano': 0.000010, 'amazonian': 0.999960, 'european': 0.000010, 'northern_andean': 0.000010, 'african': 0.000010},
    'MACHIGUENGA': {'altiplano': 0.000010, 'amazonian': 0.999960, 'european': 0.000010, 'northern_andean': 0.000010, 'african': 0.000010},
    'ASHANINKA': {'altiplano': 0.000010, 'amazonian': 0.999960, 'european': 0.000010, 'northern_andean': 0.000010, 'african': 0.000010},
    'MATSES': {'altiplano': 0.000010, 'amazonian': 0.999553, 'european': 0.000010, 'northern_andean': 0.000417, 'african': 0.000010},
    'SHIPIBO': {'altiplano': 0.008162, 'amazonian': 0.969681, 'european': 0.003973, 'northern_andean': 0.018009, 'african': 0.000176},
    'CANDOSHI': {'altiplano': 0.019273, 'amazonian': 0.903211, 'european': 0.000010, 'northern_andean': 0.077496, 'african': 0.000010},
    'AWAJUN': {'altiplano': 0.000010, 'amazonian': 0.894368, 'european': 0.007811, 'northern_andean': 0.097801, 'african': 0.000010},
    'IQUITOS': {'altiplano': 0.033317, 'amazonian': 0.581600, 'european': 0.258251, 'northern_andean': 0.111944, 'african': 0.014887},
    'TUMBES': {'altiplano': 0.041726, 'amazonian': 0.191223, 'european': 0.375806, 'northern_andean': 0.270387, 'african': 0.120858},
    'TRUJILLO': {'altiplano': 0.168178, 'amazonian': 0.105176, 'european': 0.213626, 'northern_andean': 0.370167, 'african': 0.142853},
    'CHACHAPOYAS': {'altiplano': 0.137705, 'amazonian': 0.258828, 'european': 0.183545, 'northern_andean': 0.416652, 'african': 0.003270},
    'LAMBAYEQUE': {'altiplano': 0.098971, 'amazonian': 0.238287, 'european': 0.138987, 'northern_andean': 0.473948, 'african': 0.049806},
    'TALLAN': {'altiplano': 0.051878, 'amazonian': 0.317262, 'european': 0.013198, 'northern_andean': 0.614349, 'african': 0.003313},
    'MOCHES': {'altiplano': 0.014051, 'amazonian': 0.027841, 'european': 0.044305, 'northern_andean': 0.908668, 'african': 0.005135},
    'TACNA': {'altiplano': 0.762635, 'amazonian': 0.072269, 'european': 0.137561, 'northern_andean': 0.013699, 'african': 0.013836},
    'PUNO': {'altiplano': 0.865068, 'amazonian': 0.082078, 'european': 0.046868, 'northern_andean': 0.004637, 'african': 0.001349},
    'QUEROS': {'altiplano': 0.945745, 'amazonian': 0.053950, 'european': 0.000010, 'northern_andean': 0.000285, 'african': 0.000010},
    'UROS': {'altiplano': 0.965716, 'amazonian': 0.021837, 'european': 0.006582, 'northern_andean': 0.005855, 'african': 0.000010},
    'CHOPCCAS': {'altiplano': 0.823188, 'amazonian': 0.000855, 'european': 0.008692, 'northern_andean': 0.167084, 'african': 0.000182},
    'JAQARUS': {'altiplano': 0.718978, 'amazonian': 0.002548, 'european': 0.129049, 'northern_andean': 0.144246, 'african': 0.005178},
    'AYACUCHO': {'altiplano': 0.580537, 'amazonian': 0.057747, 'european': 0.191684, 'northern_andean': 0.167697, 'african': 0.002334},
    'CUSCO': {'altiplano': 0.635210, 'amazonian': 0.043500, 'european': 0.234112, 'northern_andean': 0.076084, 'african': 0.011094},
    'MOQUEGUA': {'altiplano': 0.624672, 'amazonian': 0.052500, 'european': 0.282173, 'northern_andean': 0.013650, 'african': 0.027005},
    'AREQUIPA': {'altiplano': 0.437996, 'amazonian': 0.026133, 'european': 0.506247, 'northern_andean': 0.016235, 'african': 0.013389},
    'LIMA': {'altiplano': 0.352586, 'amazonian': 0.088634, 'european': 0.322483, 'northern_andean': 0.213226, 'african': 0.023071},
    'HUARAZ': {'altiplano': 0.392404, 'amazonian': 0.093215, 'european': 0.182882, 'northern_andean': 0.328893, 'african': 0.002607},
    'AFRODESCENDIENTES': {'altiplano': 0.129592, 'amazonian': 0.036752, 'european': 0.079208, 'northern_andean': 0.086718, 'african': 0.667730},
}

# Ancestry colors - matching Figure 4 ADMIXTURE plot (K=5)
ANCESTRY_COLORS = {
    'altiplano':       '#E07B73',  # Coral/salmon red (Ancestry 1 - Altiplano Andean)
    'amazonian':       '#5B9B9B',  # Teal/turquoise (Ancestry 2 - Amazonian Native)
    'european':        '#6B98C4',  # Blue (Ancestry 3 - European/Coastal)
    'northern_andean': '#E5A83B',  # Orange/yellow (Ancestry 4 - Northern Andean)
    'african':         '#9B7BB5',  # Purple/violet (Ancestry 5 - African)
}

# Map population names to ADMIXTURE keys (only populations with ADMIXTURE data)
ADMIX_KEY_MAP = {
    'Matses': 'MATSES', 'Awajun': 'AWAJUN', 'Shipibo': 'SHIPIBO',
    'Ashaninka': 'ASHANINKA', 'Machiguenga': 'MACHIGUENGA', 'Nahua': 'NAHUA',
    'Candoshi': 'CANDOSHI', 'Lamas': 'LAMAS', 'Iquitos': 'IQUITOS',
    'Chopccas': 'CHOPCCAS', 'Queros': 'QUEROS', 'Chachapoyas': 'CHACHAPOYAS',
    'Huaraz': 'HUARAZ', 'Ayacucho': 'AYACUCHO', 'Uros': 'UROS',
    'Puno': 'PUNO', 'Jaqarus': 'JAQARUS', 'Cusco': 'CUSCO',
    'Arequipa': 'AREQUIPA', 'Moches': 'MOCHES', 'Tallan': 'TALLAN',
    'Trujillo': 'TRUJILLO', 'Lima': 'LIMA', 'Lambayeque': 'LAMBAYEQUE',
    'Tumbes': 'TUMBES', 'Tacna': 'TACNA', 'Moquegua': 'MOQUEGUA',
    'Afrodescendientes': 'AFRODESCENDIENTES',
}

def draw_pie_chart(ax, lon, lat, proportions, colors, radius=0.35, edgecolor='black', linewidth=0.5):
    """Draw a pie chart at the specified location."""
    theta1 = 90  # Start at top
    wedges = []
    for prop, color in zip(proportions, colors):
        if prop > 0.001:  # Skip tiny wedges
            theta2 = theta1 - prop * 360
            wedge = Wedge((lon, lat), radius, theta2, theta1, facecolor=color, 
                          edgecolor=edgecolor, linewidth=linewidth)
            wedges.append(wedge)
            theta1 = theta2
    return wedges

# ============================================================================
# POPULATION DATA WITH ACCURATE COORDINATES AND SAMPLE COUNTS
# 31 populations based on actual sample data
# Coordinates sourced from Peruvian geographic references
# ============================================================================

populations = pd.DataFrame([
    # ==========================================================================
    # AMAZONIAN REGION - Indigenous
    # ==========================================================================
    {"name": "Matses", "lat": -5.30, "lon": -73.85, "region": "Amazonian", 
     "type": "Indigenous", "n": 17,
     "notes": "Loreto, near Brazilian border"},
    
    {"name": "Awajun", "lat": -4.50, "lon": -77.80, "region": "Amazonian", 
     "type": "Indigenous", "n": 25,
     "notes": "Amazonas region, Marañón basin (includes AWARUNA)"},
    
    {"name": "Shipibo", "lat": -8.10, "lon": -74.35, "region": "Amazonian", 
     "type": "Indigenous", "n": 14,
     "notes": "Ucayali region, Yarinacocha area"},
    
    {"name": "Ashaninka", "lat": -10.80, "lon": -74.60, "region": "Amazonian", 
     "type": "Indigenous", "n": 1,
     "notes": "Junín/Ucayali border region"},
    
    {"name": "Machiguenga", "lat": -11.90, "lon": -72.50, "region": "Amazonian", 
     "type": "Indigenous", "n": 3,
     "notes": "Cusco region, Urubamba basin"},
    
    {"name": "Nahua", "lat": -11.35, "lon": -71.80, "region": "Amazonian", 
     "type": "Indigenous", "n": 2,
     "notes": "Madre de Dios region"},
    
    {"name": "Candoshi", "lat": -4.30, "lon": -76.80, "region": "Amazonian", 
     "type": "Indigenous", "n": 17,
     "notes": "Loreto region, Pastaza basin"},
    
    {"name": "Lamas", "lat": -6.42, "lon": -76.52, "region": "Amazonian", 
     "type": "Indigenous", "n": 29,
     "notes": "San Martín region, Quechua-Lamista"},

    # ==========================================================================
    # AMAZONIAN REGION - Mestizo
    # ==========================================================================
    {"name": "Iquitos", "lat": -3.75, "lon": -73.25, "region": "Amazonian", 
     "type": "Mestizo", "n": 68,
     "notes": "Capital of Loreto region"},

    # ==========================================================================
    # HIGHLAND REGION - Indigenous Quechua-speaking
    # ==========================================================================
    {"name": "Chopccas", "lat": -12.60, "lon": -75.00, "region": "Highland", 
     "type": "Indigenous", "n": 56, "language": "Quechua",
     "notes": "Huancavelica region"},
    
    {"name": "Queros", "lat": -13.35, "lon": -70.95, "region": "Highland", 
     "type": "Indigenous", "n": 12, "language": "Quechua",
     "notes": "Cusco region, high altitude"},
    
    {"name": "Chachapoyas", "lat": -6.23, "lon": -77.87, "region": "Amazonian", 
     "type": "Indigenous", "n": 47, "language": "Quechua",
     "notes": "Amazonas region, cloud forest transition zone"},
    
    {"name": "Huaraz", "lat": -9.53, "lon": -77.53, "region": "Highland", 
     "type": "Indigenous", "n": 80, "language": "Quechua",
     "notes": "Ancash region capital, Cordillera Blanca"},
    
    {"name": "Ayacucho", "lat": -13.35, "lon": -74.22, "region": "Highland", 
     "type": "Indigenous", "n": 63, "language": "Quechua",
     "notes": "Ayacucho region"},

    # ==========================================================================
    # HIGHLAND REGION - Indigenous Aymara-speaking
    # ==========================================================================
    {"name": "Uros", "lat": -15.35, "lon": -69.50, "region": "Highland", 
     "type": "Indigenous", "n": 44, "language": "Aymara",
     "notes": "Lake Titicaca floating islands"},
    
    {"name": "Puno", "lat": -16.10, "lon": -70.30, "region": "Highland", 
     "type": "Indigenous", "n": 104, "language": "Aymara",
     "notes": "Puno region, Lake Titicaca shore"},
    
    {"name": "Jaqarus", "lat": -12.50, "lon": -75.92, "region": "Highland", 
     "type": "Indigenous", "n": 20, "language": "Aymara",
     "notes": "Lima region, Yauyos province - isolated Aymara speakers"},

    # ==========================================================================
    # HIGHLAND REGION - Mestizo
    # ==========================================================================
    {"name": "Cusco", "lat": -13.70, "lon": -72.10, "region": "Highland", 
     "type": "Mestizo", "n": 82,
     "notes": "Former Inca capital"},
    
    {"name": "Arequipa", "lat": -16.40, "lon": -71.54, "region": "Highland", 
     "type": "Mestizo", "n": 34,
     "notes": "Second largest city"},

    # ==========================================================================
    # COASTAL REGION - Indigenous
    # ==========================================================================
    {"name": "Moches", "lat": -7.45, "lon": -79.30, "region": "Coastal", 
     "type": "Indigenous", "n": 75,
     "notes": "La Libertad region, descendants of Moche culture"},
    
    {"name": "Tallan", "lat": -5.20, "lon": -80.63, "region": "Coastal", 
     "type": "Indigenous", "n": 42,
     "notes": "Piura region, descendants of Tallán culture"},

    # ==========================================================================
    # COASTAL REGION - Mestizo
    # ==========================================================================
    {"name": "Trujillo", "lat": -8.50, "lon": -79.05, "region": "Coastal", 
     "type": "Mestizo", "n": 69,
     "notes": "Third largest city"},
    
    {"name": "Lima", "lat": -12.05, "lon": -77.04, "region": "Coastal", 
     "type": "Mestizo", "n": 29,
     "notes": "Capital city"},
    
    {"name": "Lambayeque", "lat": -6.77, "lon": -79.84, "region": "Coastal", 
     "type": "Mestizo", "n": 12,
     "notes": "Lambayeque region"},
    
    {"name": "Tumbes", "lat": -3.57, "lon": -80.45, "region": "Coastal", 
     "type": "Mestizo", "n": 33,
     "notes": "Far northern coast, Ecuador border"},
    
    {"name": "Tacna", "lat": -18.01, "lon": -70.25, "region": "Coastal", 
     "type": "Mestizo", "n": 30,
     "notes": "Southern coast, Chile border"},
    
    {"name": "Moquegua", "lat": -17.19, "lon": -70.94, "region": "Coastal", 
     "type": "Mestizo", "n": 31,
     "notes": "Southern coast"},

    # ==========================================================================
    # AFRO-PERUVIAN
    # ==========================================================================
    {"name": "Afrodescendientes", "lat": -13.45, "lon": -76.15, "region": "Coastal", 
     "type": "Afroperuano", "n": 103,
     "notes": "Chincha area, Ica region"},
])

# Fill missing language column
if 'language' not in populations.columns:
    populations['language'] = ''
populations['language'] = populations['language'].fillna('')

# Convert to GeoDataFrame
geometry = [Point(xy) for xy in zip(populations['lon'], populations['lat'])]
pop_gdf = gpd.GeoDataFrame(populations, geometry=geometry, crs="EPSG:4326")

# ============================================================================
# CREATE FIGURE
# ============================================================================

fig = plt.figure(figsize=(16, 14))

# Main map - Peru detail
ax_main = fig.add_axes([0.05, 0.05, 0.70, 0.90])

# South America inset (top left)
ax_sa = fig.add_axes([0.02, 0.72, 0.20, 0.25])

# Legend area (right side)
ax_legend = fig.add_axes([0.76, 0.15, 0.22, 0.70])
ax_legend.axis('off')

# ============================================================================
# PANEL B: MAIN MAP OF PERU
# ============================================================================

# Plot neighboring countries (light grey)
neighbors.plot(ax=ax_main, color='#E8E8E8', edgecolor='#AAAAAA', linewidth=0.5)

# Plot Peru base with slight color
peru.plot(ax=ax_main, color='#FAFAFA', edgecolor='black', linewidth=1.5)

# Add region shading - carefully designed to cover all of Peru
# and place each population in its correct ecological region
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.path import Path

# Get Peru boundary as a matplotlib path for clipping
peru_geom = peru.geometry.values[0]
if peru_geom.geom_type == 'MultiPolygon':
    peru_poly = max(peru_geom.geoms, key=lambda x: x.area)
else:
    peru_poly = peru_geom
peru_coords = list(peru_poly.exterior.coords)

# Define region boundaries based on population locations:
# Coastal: Piura(-80.63), Chiclayo(-79.84), Moches(-79.05), Trujillo(-79.03), 
#          Lima(-77.04), Afroperuanos(-76.15), Ica(-75.73)
# Highland: Cajamarca(-78.51), Chachapoyas(-77.87), Ancash(-77.81), Huaraz(-77.53),
#           Huánuco(-76.24), Jaqaru(-75.92), Huancayo(-75.21), Chopccas(-74.87),
#           Ayacucho(-74.22), Apurímac(-72.88), Cusco(-71.97), Arequipa(-71.54),
#           Q'eros(-70.95), Puno(-70.02), Uros(-69.87)
# Amazonian: Awajún(-77.80), Pucallpa(-74.55), Asháninka(-74.40), Shipibo-Konibo(-74.35),
#            Matzes(-73.85), Iquitos(-73.25), Machiguenga(-72.85), Nahua(-72.25)

# REGION BOUNDARIES - PRECISELY TRACED FROM REFERENCE MAP
# "Regiones geográficas del Perú" - AUTHORITATIVE SOURCE
#
# CRITICAL OBSERVATIONS FROM REFERENCE:
# 1. Costa (yellow) = EXTREMELY thin - barely visible strip, ~30-50km wide
# 2. Sierra (brown) = Narrow Andes corridor, narrower north, wider south (altiplano)
# 3. Selva (green) = DOMINANT - ~60%+ of Peru, extends far west
#
# The Sierra-Selva boundary is traced point-by-point from reference

# ============================================================================
# COASTAL REGION (Costa) - Extended to cover all coastal populations
# Coastal pops: Piura(-80.63), Chiclayo(-79.84), Moches(-79.05), Trujillo(-79.03), 
#               Lima(-77.04), Afroperuanos(-76.15), Ica(-75.73)
#               Tumbes(-80.45, -3.57), Tacna(-70.25, -18.01), Moquegua(-70.94, -17.19)
# Must cover: Northern tip including Tumbes, and southern coast including Tacna/Moquegua
# ============================================================================
coastal_verts = [
    (-82, 2),        # NW ocean
    (-80.0, -3.0),   # Include Tumbes at (-80.45, -3.57)
    (-80.0, -4.0),   # Tumbes area
    (-80.3, -4.8),   # Include Tallan/Piura
    (-79.5, -5.6),   
    (-79.3, -6.4),   # Include Lambayeque
    (-78.8, -7.2),   
    (-78.6, -8.0),   # Include Moches/Trujillo
    (-78.3, -8.8),   
    (-77.8, -9.6),   
    (-77.4, -10.4),  
    (-77.0, -11.2),  
    (-76.6, -12.0),  # Include Lima
    (-76.2, -12.8),  
    (-75.8, -13.2),  # Include Afrodescendientes
    (-75.4, -13.8),  
    (-75.2, -14.4),  
    (-74.6, -15.0),  
    (-74.0, -15.6),  
    (-73.4, -16.2),  
    (-72.8, -16.8),  
    (-72.2, -17.2),  # Include Moquegua at (-70.94, -17.19)
    (-70.6, -17.4),  # Extended east for Moquegua
    (-70.0, -17.8),  # Include Tacna at (-70.25, -18.01)
    (-69.8, -18.35), # Chile border - extended east
    (-68.5, -18.35), # Eastern extent at Chile border
    (-82, -18.5),    # SW ocean
    (-82, 2),
]

# ============================================================================
# HIGHLAND REGION (Sierra) - Extended to cover all highland populations
# Highland pops: Cajamarca(-78.51), Chachapoyas(-77.87, -6.23), Ancash(-77.81), Huaraz(-77.53),
#                Huánuco(-76.24), Jaqaru(-75.92), Huancayo(-75.21), Chopccas(-74.87),
#                Ayacucho(-74.22), Apurímac(-72.88), Cusco(-71.97), Arequipa(-71.54),
#                Q'eros(-70.95), Puno(-70.02, -15.84), Uros(-69.87, -15.82)
# NOTE: Awajún(-77.80, -4.50), Asháninka(-74.40) and Machiguenga(-72.85) are AMAZONIAN
# Highland must be NARROW in the north - no highland at the very tip (that's coastal)
# Titicaca area (Uros, Puno) must be in Highland
# Southern boundary must stay WEST of Tacna/Moquegua (they are Coastal)
# ============================================================================
highland_verts = [
    # Western boundary - starts south of the northern coastal tip
    (-80.0, -3.0),   # Start - north of this is coastal only
    (-80.0, -4.0),   
    (-80.3, -4.8),
    (-79.5, -5.6),
    (-79.3, -6.4),
    (-78.8, -7.2),
    (-78.6, -8.0),
    (-78.3, -8.8),
    (-77.8, -9.6),
    (-77.4, -10.4),
    (-77.0, -11.2),
    (-76.6, -12.0),
    (-76.2, -12.8),
    (-75.8, -13.2),
    (-75.4, -13.8),
    (-75.2, -14.4),
    (-74.6, -15.0),
    (-74.0, -15.6),
    (-73.4, -16.2),
    (-72.8, -16.8),
    (-72.2, -17.2),  # Match coastal boundary
    (-70.6, -17.4),  # Match coastal - Moquegua is coastal
    (-70.0, -17.8),  # Match coastal - Tacna is coastal
    (-69.8, -18.35),
    # Southern border - extends east to include Altiplano/Titicaca
    (-68.5, -18.35),
    # =========================================================================
    # EASTERN BOUNDARY - Must include Uros(-69.50) and Puno(-70.02) in Highland
    # Going NORTH from southern border
    # =========================================================================
    (-68.5, -17.8),  # Extended EAST to include Titicaca
    (-68.5, -17.0),  # Still east to keep Titicaca in Highland
    (-68.8, -16.2),  # Widen Highland here
    (-69.0, -15.6),  # Keep boundary EAST of Uros at -69.50
    (-69.0, -15.0),  # Wide Highland continues
    (-69.5, -14.4),  # Now start narrowing
    (-70.5, -14.0),  # Cusco region
    (-70.6, -13.6),  # Include Q'eros (-70.95)
    (-70.6, -13.0),  
    (-71.6, -12.6),  # Apurímac - curve back west
    (-73.2, -12.2),  # WEST of Machiguenga (-72.85)
    (-74.0, -11.6),  # Ayacucho
    (-74.6, -11.2),  # WEST of Asháninka (-74.40)
    (-75.0, -10.6),  # Chopccas/Huancayo area
    (-75.4, -10.0),  # Jaqaru
    (-75.8, -9.4),   # Huánuco
    (-76.4, -8.8),   # Huaraz/Ancash
    (-76.8, -8.0),   
    (-77.2, -7.2),   # Cajamarca
    (-77.5, -6.8),   # Highland narrows significantly
    (-78.0, -6.2),   # Highland ends WEST of Chachapoyas (-77.87) - Chachapoyas in Amazonian
    (-78.2, -5.6),   # Narrowing
    (-78.4, -5.0),   # WEST of Awajún (-77.80)
    (-78.6, -4.4),   # Narrow highland
    (-78.6, -3.8),   # Very narrow at top
    (-79.2, -3.2),   # Narrow strip
    (-80.0, -3.0),   # Back to start - highland ends here
]

# ============================================================================
# AMAZONIAN REGION (Selva) - Large green region
# Western boundary = Sierra eastern boundary (reversed)
# Must include: Awajún(-77.80, -4.50), Asháninka(-74.40), Machiguenga(-72.85), Nahua(-72.25)
# Must EXCLUDE: Uros, Puno (they are Highland)
# ============================================================================
amazon_verts = [
    # Western boundary (= Sierra eastern, going south) - MUST MATCH highland_verts eastern edge
    (-80.0, -3.0),   # Start at highland northern tip
    (-79.2, -3.2),
    (-78.6, -3.8),
    (-78.6, -4.4),   # Narrow highland - matches highland_verts
    (-78.4, -5.0),   # WEST of Awajún (-77.80)
    (-78.2, -5.6),   # Narrowing
    (-78.0, -6.2),   # Highland ends WEST of Chachapoyas - matches highland_verts
    (-77.5, -6.8),   # Highland narrows significantly - matches highland_verts
    (-77.2, -7.2),
    (-76.8, -8.0),
    (-76.4, -8.8),
    (-75.8, -9.4),
    (-75.4, -10.0),
    (-75.0, -10.6),
    (-74.6, -11.2),  # East of this: Asháninka at -74.40
    (-74.0, -11.6),
    (-73.2, -12.2),  # East of this: Machiguenga at -72.85
    (-71.6, -12.6),
    (-70.6, -13.0),
    (-70.6, -13.6),
    (-70.5, -14.0),
    (-69.5, -14.4),
    (-69.0, -15.0),
    (-69.0, -15.6),  # Boundary stays EAST - Uros/Puno in Highland
    (-68.8, -16.2),
    (-68.5, -17.0),
    (-68.5, -17.8),
    (-68.5, -18.35),
    # Eastern/Northern Peru borders
    (-68.5, -18.35), # Bolivia border
    (-68.5, -14.0),  
    (-68.7, -11.0),  # Brazil border
    (-69.0, -8.0),   
    (-69.5, -5.0),   
    (-70.0, -2.5),   # Colombia border
    (-71.0, -1.0),   
    (-73.0, 0),      # Ecuador border
    (-75.5, -0.2),   
    (-78.0, -0.5),   
    (-79.0, -1.5),   # Connect back to highland start
    (-80.0, -3.0),   # Back to start
]

# Create and add patches with clipping to Peru boundary
clip_path = MplPolygon(peru_coords, closed=True, transform=ax_main.transData)

coastal_patch = MplPolygon(coastal_verts, closed=True, facecolor='#6BAED6', 
                           alpha=0.35, edgecolor='none', zorder=1)
coastal_patch.set_clip_path(clip_path)
ax_main.add_patch(coastal_patch)

highland_patch = MplPolygon(highland_verts, closed=True, facecolor='#D9A84E',
                            alpha=0.35, edgecolor='none', zorder=1)
highland_patch.set_clip_path(clip_path)
ax_main.add_patch(highland_patch)

amazon_patch = MplPolygon(amazon_verts, closed=True, facecolor='#74C476',
                          alpha=0.35, edgecolor='none', zorder=1)
amazon_patch.set_clip_path(clip_path)
ax_main.add_patch(amazon_patch)

# Re-plot Peru outline on top
peru.boundary.plot(ax=ax_main, color='black', linewidth=1.5, zorder=3)

# ============================================================================
# PLOT POPULATIONS WITH SIZE SCALED BY SAMPLE COUNT
# ============================================================================

# Define colors and markers
region_colors = {
    'Coastal': '#1f77b4',      # Blue
    'Highland': '#d62728',      # Red
    'Amazonian': '#2ca02c'      # Green
}

type_markers = {
    'Indigenous': 'o',          # Circle
    'Mestizo': 's',             # Square
    'Afroperuano': '^'        # Triangle
}

# Calculate marker sizes based on sample count
# Scale: min size 40, max size 400
min_n = populations['n'].min()
max_n = populations['n'].max()
def scale_size(n):
    # Linear scaling from 40 to 400
    return 40 + (n - min_n) / (max_n - min_n) * 360

def get_pie_radius(n):
    """Scale pie chart radius from 0.18 to 0.42 based on sample size."""
    return 0.18 + (n - min_n) / (max_n - min_n) * 0.24

# Plot each population
for idx, pop in populations.iterrows():
    if USE_PIE_CHARTS:
        # PIE CHART MODE
        admix_key = ADMIX_KEY_MAP.get(pop['name'])
        radius = get_pie_radius(pop['n'])
        
        if admix_key and admix_key in admixture_data:
            # Draw pie chart for populations with ADMIXTURE data (K=5)
            data = admixture_data[admix_key]
            # Order: Altiplano, Amazonian, European, Northern Andean, African
            proportions = [data['altiplano'], data['amazonian'], data['european'], 
                          data['northern_andean'], data['african']]
            colors = [ANCESTRY_COLORS['altiplano'], ANCESTRY_COLORS['amazonian'], 
                     ANCESTRY_COLORS['european'], ANCESTRY_COLORS['northern_andean'],
                     ANCESTRY_COLORS['african']]
            
            wedges = draw_pie_chart(ax_main, pop['lon'], pop['lat'], proportions, colors, 
                                   radius=radius, edgecolor='black', linewidth=0.5)
            for w in wedges:
                w.set_zorder(10)
                ax_main.add_patch(w)
        # Skip populations without ADMIXTURE data (they shouldn't exist now)
    else:
        # SHAPED MARKER MODE (original v3 behavior)
        color = region_colors[pop['region']]
        marker = type_markers[pop['type']]
        size = scale_size(pop['n'])
        
        ax_main.scatter(pop['lon'], pop['lat'], 
                       c=color, marker=marker, s=size,
                       edgecolors='white', linewidths=1.5,
                       zorder=10, alpha=0.85)

# ============================================================================
# ADD LABELS
# ============================================================================

# Label offsets for each population (adjusted for clarity)
label_config = {
    # Amazonian Indigenous - labels to the right
    'Matses': {'offset': (1, 0.02), 'ha': 'left'},
    'Awajun': {'offset': (1, 0.02), 'ha': 'left'},
    'Shipibo': {'offset': (1, 0.02), 'ha': 'left'},
    'Ashaninka': {'offset': (1, 0.02), 'ha': 'left'},
    'Machiguenga': {'offset': (1, 0.02), 'ha': 'left'},
    'Nahua': {'offset': (1, 0.02), 'ha': 'left'},
    'Candoshi': {'offset': (1, 0.02), 'ha': 'left'},
    'Lamas': {'offset': (1, 0.02), 'ha': 'left'},
    'Chachapoyas': {'offset': (1, 0.02), 'ha': 'left'},  # Amazonas region (cloud forest)
    
    # Amazonian Mestizo
    'Iquitos': {'offset': (1, 0.02), 'ha': 'left'},
    
    # Highland Quechua
    'Chopccas': {'offset': (1, 0.02), 'ha': 'left'},
    'Queros': {'offset': (1, 0.02), 'ha': 'left'},
    'Huaraz': {'offset': (1, 0.02), 'ha': 'left'},
    'Ayacucho': {'offset': (1, 0.02), 'ha': 'left'},
    
    # Highland Aymara
    'Uros': {'offset': (1, 0.02), 'ha': 'left'},
    'Puno': {'offset': (1, 0.02), 'ha': 'left'},
    'Jaqarus': {'offset': (-1, 0.02), 'ha': 'right'},  # Moved to left to avoid overlap
    
    # Highland Mestizo
    'Cusco': {'offset': (1, 0.02), 'ha': 'left'},
    'Arequipa': {'offset': (1, 0.02), 'ha': 'left'},
    
    # Coastal Indigenous - labels to the left
    'Moches': {'offset': (-1, 0.02), 'ha': 'right'},
    'Tallan': {'offset': (-1, 0.02), 'ha': 'right'},
    
    # Coastal Mestizo
    'Trujillo': {'offset': (-1, 0.02), 'ha': 'right'},
    'Lima': {'offset': (-1, 0.02), 'ha': 'right'},
    'Lambayeque': {'offset': (-1, 0.02), 'ha': 'right'},
    'Tumbes': {'offset': (-1, 0.02), 'ha': 'right'},
    'Tacna': {'offset': (1, 0.02), 'ha': 'left'},
    'Moquegua': {'offset': (1, 0.02), 'ha': 'left'},
    
    # Afroperuano
    'Afrodescendientes': {'offset': (-1, 0.02), 'ha': 'right'},
}

for idx, pop in populations.iterrows():
    config = label_config.get(pop['name'], {'offset': (1, 0.05), 'ha': 'left'})
    
    # Get pie radius for this population
    pie_radius = get_pie_radius(pop['n'])
    # Label offset = pie radius + small buffer (0.05 degrees)
    label_distance = pie_radius + 0.05
    
    # Apply direction from config
    base_x, base_y = config['offset']
    if base_x > 0:
        scaled_x = label_distance
    elif base_x < 0:
        scaled_x = -label_distance
    else:
        scaled_x = 0
    scaled_y = base_y
    
    # Consistent typography for all labels - Nature style
    ax_main.annotate(
        pop['name'],
        xy=(pop['lon'], pop['lat']),
        xytext=(pop['lon'] + scaled_x, pop['lat'] + scaled_y),
        fontsize=7,
        fontweight='normal',
        fontfamily='sans-serif',
        ha=config['ha'],
        va='center',
        color='#333333',
        zorder=15
    )

# ============================================================================
# ADD MAP ANNOTATIONS
# ============================================================================

# "PERU" label - subtle watermark
ax_main.text(-74.5, -10.0, 'PERU', fontsize=28, fontweight='bold',
            ha='center', va='center', alpha=0.12, style='italic',
            color='#444444', zorder=1, fontfamily='sans-serif')

# Pacific Ocean - muted
ax_main.text(-82.0, -10.0, 'Pacific\nOcean', fontsize=9, style='italic',
            ha='center', va='center', color='#5B7B8C', alpha=0.6)

# Neighboring countries - subtle, consistent
ax_main.text(-79.5, -1.5, 'ECUADOR', fontsize=7, ha='center', color='#888888')
ax_main.text(-72.0, -1.5, 'COLOMBIA', fontsize=7, ha='center', color='#888888')
ax_main.text(-68.0, -8.0, 'BRAZIL', fontsize=7, ha='center', color='#888888')
ax_main.text(-68.5, -15.5, 'BOLIVIA', fontsize=7, ha='center', color='#888888')
ax_main.text(-70.0, -19.5, 'CHILE', fontsize=7, ha='center', color='#888888')

# Lake Titicaca - muted
from matplotlib.patches import Ellipse
lake = Ellipse((-69.5, -15.8), 1.0, 0.6, color='#5B7B8C', alpha=0.4, zorder=5)
ax_main.add_patch(lake)
ax_main.text(-68.5, -15.3, 'Lake\nTiticaca', fontsize=6, style='italic',
            ha='center', color='#5B7B8C', alpha=0.7)

# Region labels - subtle, no boxes
ax_main.text(-80.2, -9.5, 'COASTAL', fontsize=9, fontweight='medium',
            ha='center', va='center', color='#5B7B8C', alpha=0.6, rotation=75,
            fontfamily='sans-serif')
ax_main.text(-75.5, -11.5, 'HIGHLAND', fontsize=9, fontweight='medium',
            ha='center', va='center', color='#8B7355', alpha=0.6, rotation=45,
            fontfamily='sans-serif')
ax_main.text(-76.0, -5.5, 'AMAZONIAN', fontsize=9, fontweight='medium',
            ha='center', va='center', color='#5B7B5B', alpha=0.6,
            fontfamily='sans-serif')

# ============================================================================
# MAP ELEMENTS
# ============================================================================

# Set extent
ax_main.set_xlim(-83, -67.5)
ax_main.set_ylim(-19, 0.5)

# Refined axis styling
for spine in ax_main.spines.values():
    spine.set_linewidth(0.5)
    spine.set_color('#444444')
ax_main.tick_params(axis='both', which='major', labelsize=8, width=0.5, 
                   length=3, color='#444444', labelcolor='#444444')

# Labels
ax_main.set_xlabel('Longitude', fontsize=9, color='#444444')
ax_main.set_ylabel('Latitude', fontsize=9, color='#444444')

# Grid - subtle
ax_main.grid(True, linestyle='-', alpha=0.15, color='#CCCCCC', linewidth=0.3)

# North arrow - refined
ax_main.annotate('', xy=(-82, 0), xytext=(-82, -1.5),
                arrowprops=dict(arrowstyle='->', lw=1, color='#444444'))
ax_main.text(-82, 0.3, 'N', fontsize=10, fontweight='bold', ha='center', color='#444444')

# Scale bar - refined
scale_lon, scale_lat = -81.5, -17.5
ax_main.plot([scale_lon, scale_lon + 2], [scale_lat, scale_lat], 
            color='#444444', linewidth=2, solid_capstyle='butt')
ax_main.text(scale_lon + 1, scale_lat - 0.5, '200 km', ha='center', fontsize=7, color='#444444')

# Panel label
ax_main.text(-67.0, 0.3, 'B', fontsize=14, fontweight='bold', color='#333333')

# ============================================================================
# SOUTH AMERICA INSET (matching v3 style)
# ============================================================================

# Title above inset
ax_sa.set_title('South America', fontsize=9, fontweight='normal', pad=2, color='#333333')

south_america.plot(ax=ax_sa, color='#D3D3D3', edgecolor='#999999', linewidth=0.3)
peru.plot(ax=ax_sa, color='#CD5C5C', edgecolor='#8B0000', linewidth=0.8)  # Red like v3

ax_sa.set_xlim(-85, -30)
ax_sa.set_ylim(-60, 15)
ax_sa.set_aspect('equal')
ax_sa.axis('off')

# Add "Peru" label on the country
peru_centroid = peru.geometry.centroid.iloc[0]
ax_sa.text(peru_centroid.x - 2, peru_centroid.y + 2, 'Peru', fontsize=7, fontweight='bold',
          color='#8B0000', ha='center', va='center')

# Add red box showing main map extent
from matplotlib.patches import Rectangle
rect = Rectangle((-83, -19), 15.5, 19.5, fill=False, edgecolor='#CD5C5C', linewidth=1.2)
ax_sa.add_patch(rect)

# ============================================================================
# LEGEND
# ============================================================================

# Title
ax_legend.text(0.5, 0.98, 'Legend', fontsize=11, fontweight='bold', 
              ha='center', va='top', transform=ax_legend.transAxes, color='#333333')

# Count populations by region
coastal_n = len(populations[populations['region'] == 'Coastal'])
highland_n = len(populations[populations['region'] == 'Highland'])
amazon_n = len(populations[populations['region'] == 'Amazonian'])

# Count by type
indigenous_n = len(populations[populations['type'] == 'Indigenous'])
mestizo_n = len(populations[populations['type'] == 'Mestizo'])
afro_n = len(populations[populations['type'] == 'Afroperuano'])

if USE_PIE_CHARTS:
    # =========================================================================
    # PIE CHART MODE LEGEND - clean, Nature style
    # =========================================================================
    
    # Ecological Regions (shading)
    y = 0.92
    ax_legend.text(0.05, y, 'Ecological Regions', fontsize=9, fontweight='bold',
                  transform=ax_legend.transAxes, color='#333333')
    
    regions = [
        (f'Coastal (n={coastal_n})', '#6BAED6'),
        (f'Highland (n={highland_n})', '#D9A84E'),
        (f'Amazonian (n={amazon_n})', '#74C476')
    ]
    
    for i, (label, color) in enumerate(regions):
        y = 0.87 - 0.045 * i
        rect = plt.Rectangle((0.05, y - 0.012), 0.06, 0.024,
                              facecolor=color, edgecolor='#666666', linewidth=0.3,
                              alpha=0.5, transform=ax_legend.transAxes, clip_on=False)
        ax_legend.add_patch(rect)
        ax_legend.text(0.14, y, label, fontsize=8, va='center', transform=ax_legend.transAxes, color='#333333')
    
    # Ancestry Components (K=5)
    y = 0.70
    ax_legend.text(0.05, y, 'Ancestry (K=5)', fontsize=9, fontweight='bold',
                  transform=ax_legend.transAxes, color='#333333')
    
    # Ancestry color key - simple small squares (5 ancestries)
    anc_labels = [
        ('Altiplano Andean', ANCESTRY_COLORS['altiplano']),
        ('Amazonian Native', ANCESTRY_COLORS['amazonian']),
        ('European/Coastal', ANCESTRY_COLORS['european']),
        ('Northern Andean', ANCESTRY_COLORS['northern_andean']),
        ('African', ANCESTRY_COLORS['african']),
    ]
    for i, (label, color) in enumerate(anc_labels):
        y = 0.65 - 0.040 * i
        rect = plt.Rectangle((0.05, y - 0.010), 0.04, 0.020,
                              facecolor=color, edgecolor='#666666', linewidth=0.3,
                              transform=ax_legend.transAxes, clip_on=False)
        ax_legend.add_patch(rect)
        ax_legend.text(0.11, y, label, fontsize=7, va='center', transform=ax_legend.transAxes, color='#333333')
    
    # Sample size indicator
    y = 0.42
    ax_legend.text(0.05, y, 'Sample Size', fontsize=9, fontweight='bold',
                  transform=ax_legend.transAxes, color='#333333')
    
    # Use scatter for proper circle sizing
    size_examples = [(1, 'n = 1'), (50, 'n = 50'), (104, 'n = 104')]
    for i, (n, label) in enumerate(size_examples):
        y = 0.37 - 0.050 * i
        # Scale marker size to match pie chart appearance
        marker_size = 25 + (n - 1) / (104 - 1) * 120
        ax_legend.scatter([0.08], [y], s=marker_size, c='#888888', edgecolors='#666666',
                         linewidths=0.3, transform=ax_legend.transAxes, clip_on=False)
        ax_legend.text(0.14, y, label, fontsize=8, va='center', transform=ax_legend.transAxes, color='#333333')
    
    # Total samples
    total_n = populations['n'].sum()
    ax_legend.text(0.05, 0.14, f'Total: {total_n:,} individuals', fontsize=8, 
                  fontweight='bold', transform=ax_legend.transAxes, color='#333333')
    ax_legend.text(0.05, 0.08, f'{len(populations)} populations', fontsize=8,
                  transform=ax_legend.transAxes, color='#333333')

else:
    # =========================================================================
    # SHAPED MARKER MODE LEGEND - original v3 style
    # =========================================================================
    
    # Region colors (shading)
    y_pos = 0.92
    ax_legend.text(0.1, y_pos, 'Ecological Regions:', fontsize=11, fontweight='bold',
                  transform=ax_legend.transAxes)
    
    regions = [
        (f'Coastal (n={coastal_n})', '#6BAED6'),
        (f'Highland (n={highland_n})', '#D9A84E'),
        (f'Amazonian (n={amazon_n})', '#74C476')
    ]
    
    for i, (label, color) in enumerate(regions):
        y = y_pos - 0.05 * (i + 1)
        rect = FancyBboxPatch((0.1, y - 0.012), 0.08, 0.025, 
                              boxstyle="round,pad=0.01", facecolor=color, 
                              edgecolor='black', alpha=0.5,
                              transform=ax_legend.transAxes)
        ax_legend.add_patch(rect)
        ax_legend.text(0.22, y, label, fontsize=9, va='center', transform=ax_legend.transAxes)
    
    # Population types (marker shapes)
    y_pos = 0.70
    ax_legend.text(0.1, y_pos, 'Population Type:', fontsize=11, fontweight='bold',
                  transform=ax_legend.transAxes)
    
    types = [
        (f'Indigenous (n={indigenous_n})', 'o', 'gray'),
        (f'Mestizo (n={mestizo_n})', 's', 'gray'),
        (f'Afroperuano (n={afro_n})', '^', 'gray')
    ]
    
    for i, (label, marker, color) in enumerate(types):
        y = y_pos - 0.05 * (i + 1)
        ax_legend.scatter([0.14], [y], marker=marker, s=80, c=color, 
                         edgecolors='black', linewidths=1,
                         transform=ax_legend.transAxes, clip_on=False)
        ax_legend.text(0.22, y, label, fontsize=9, va='center', transform=ax_legend.transAxes)
    
    # Region color coding (marker colors)
    y_pos = 0.48
    ax_legend.text(0.1, y_pos, 'Region Color Coding:', fontsize=11, fontweight='bold',
                  transform=ax_legend.transAxes)
    
    region_markers = [
        ('Coastal', 'o', '#1f77b4'),
        ('Highland', 'o', '#d62728'),
        ('Amazonian', 'o', '#2ca02c')
    ]
    
    for i, (label, marker, color) in enumerate(region_markers):
        y = y_pos - 0.05 * (i + 1)
        ax_legend.scatter([0.14], [y], marker=marker, s=80, c=color,
                         edgecolors='black', linewidths=1,
                         transform=ax_legend.transAxes, clip_on=False)
        ax_legend.text(0.22, y, label, fontsize=9, va='center', transform=ax_legend.transAxes)
    
    # Sample size indicator
    y_pos = 0.26
    ax_legend.text(0.1, y_pos, 'Sample Size:', fontsize=11, fontweight='bold',
                  transform=ax_legend.transAxes)
    
    size_examples = [(1, 'n=1'), (30, 'n=30'), (104, 'n=104')]
    for i, (n, label) in enumerate(size_examples):
        y = y_pos - 0.05 * (i + 1)
        size = 40 + (n - min_n) / (max_n - min_n) * 360
        ax_legend.scatter([0.14], [y], marker='o', s=size * 0.5, c='gray',
                         edgecolors='black', linewidths=1,
                         transform=ax_legend.transAxes, clip_on=False)
        ax_legend.text(0.22, y, label, fontsize=9, va='center', transform=ax_legend.transAxes)
    
    # Total samples
    total_n = populations['n'].sum()
    ax_legend.text(0.1, 0.06, f'Total: {total_n} individuals', fontsize=10, 
                  fontweight='bold', transform=ax_legend.transAxes)
    ax_legend.text(0.1, 0.01, f'{len(populations)} populations', fontsize=10,
                  transform=ax_legend.transAxes)

# ============================================================================
# SAVE FIGURE
# ============================================================================

# Use current directory for output
import os
output_dir = "ANALYSIS/00-17-FIGURE1"
os.makedirs(output_dir, exist_ok=True)
version = "v5_K5" if USE_PIE_CHARTS else "v3"
png_path = os.path.join(output_dir, f'Figure1_Peru_Map_{version}.png')
pdf_path = os.path.join(output_dir, f'Figure1_Peru_Map_{version}.pdf')

plt.savefig(png_path, dpi=600, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig(pdf_path, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

# ============================================================================
# PRINT SUMMARY
# ============================================================================

print("Figure 1 saved successfully!")
print(f"  PNG: {png_path}")
print(f"  PDF: {pdf_path}")
print(f"  Mode: {'K=5 ADMIXTURE pie charts' if USE_PIE_CHARTS else 'Shaped markers (v3)'}")
print(f"\n{'='*60}")
print("POPULATION SUMMARY")
print(f"{'='*60}")
print(f"Total populations: {len(populations)}")
print(f"Total individuals: {populations['n'].sum()}")
print(f"\nBy type:")
print(f"  Indigenous: {len(populations[populations['type'] == 'Indigenous'])}")
print(f"  Mestizo: {len(populations[populations['type'] == 'Mestizo'])}")
print(f"  Afroperuano: {len(populations[populations['type'] == 'Afroperuano'])}")

print(f"\nBy region:")
for region in ['Coastal', 'Highland', 'Amazonian']:
    region_pops = populations[populations['region'] == region]
    n_pops = len(region_pops)
    n_individuals = region_pops['n'].sum()
    indigenous = len(region_pops[region_pops['type'] == 'Indigenous'])
    mestizo = len(region_pops[region_pops['type'] == 'Mestizo'])
    print(f"  {region}: {n_pops} populations, {n_individuals} individuals (Indigenous: {indigenous}, Mestizo: {mestizo})")

print(f"\nTop 10 populations by sample size:")
for _, pop in populations.nlargest(10, 'n').iterrows():
    print(f"  - {pop['name']}: n={pop['n']} ({pop['region']}, {pop['type']})")

if USE_PIE_CHARTS:
    no_admix = [name for name, key in ADMIX_KEY_MAP.items() if key is None]
    print(f"\nPopulations without ADMIXTURE data (shown as colored circles):")
    for name in no_admix:
        pop_row = populations[populations['name'] == name]
        if len(pop_row) > 0:
            print(f"  - {name}: n={pop_row.iloc[0]['n']}")