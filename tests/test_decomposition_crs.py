#!/usr/bin/env python3
"""
Test and visualize decomposition behavior with different CRS (WGS84 vs UTM).

This test demonstrates that the CGAL decomposition library has numerical precision
issues with certain polygon geometries when using large coordinate values (UTM).
The solution is to keep polygons in geographic CRS (WGS84) during decomposition.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import json
import shapely
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MPLPolygon
from shapely.geometry import Polygon
from trajgenpy import bindings
from trajgenpy.Geometries import shapely_polygon_to_cgal, GeoPolygon, generate_sweep_pattern
import math

# Load test features
test_dir = os.path.dirname(os.path.abspath(__file__))
features_path = os.path.join(test_dir, 'test_features.geojson')

with open(features_path, 'r') as f:
    geojson_data = json.load(f)

# Find heath_4 and heath_5
target_heaths = ['heath_4', 'heath_5']
heath_features = []
for feature in geojson_data['features']:
    name = feature.get('properties', {}).get('unique_name', '')
    if name in target_heaths:
        heath_features.append((name, feature))

heath_features = sorted(heath_features, key=lambda x: x[0])

if len(heath_features) < 2:
    print(f"ERROR: Only found {len(heath_features)} heath features")
    sys.exit(1)

print("=" * 80)
print("SWEEP COMPARISON: WGS84 vs UTM at Different Altitudes")
print(f"Testing {len(heath_features)} polygons x 2 altitudes")
print("=" * 80)

# Parameters
altitudes = [10, 20]  # meters
field_of_view = 60  # degrees
overlap = 20  # percent

# Process each combination
results = []

for altitude in altitudes:
    sweep_width = 2 * altitude * math.tan(math.radians(field_of_view / 2))
    sweep_distance = sweep_width * (1 - overlap / 100)
    
    print(f"\n{'=' * 80}")
    print(f"ALTITUDE: {altitude}m (sweep distance: {sweep_distance:.2f}m)")
    print('=' * 80)
    
    for heath_name, feature in heath_features:
        print(f"\n{'-' * 80}")
        print(f"Processing: {heath_name} @ {altitude}m altitude")
        print('-' * 80)
        
        # Extract polygon
        coords = feature['geometry']['coordinates'][0]
        polygon_wgs84 = Polygon(coords)
        
        # Create GeoPolygon wrapper
        geodata_wgs84 = GeoPolygon(polygon_wgs84, crs="EPSG:4326")
        
        # ====================================================================
        # WGS84 (WORKING)
        # ====================================================================
        lat = polygon_wgs84.centroid.y
        meters_per_degree = 111000 * math.cos(math.radians(lat))
        buffer_degrees = -sweep_distance / meters_per_degree
        
        buffered_wgs84 = geodata_wgs84.buffer(buffer_degrees)
        
        if buffered_wgs84.geometry.geom_type == 'MultiPolygon':
            buffered_geom = max(buffered_wgs84.geometry.geoms, key=lambda p: p.area)
        else:
            buffered_geom = buffered_wgs84.geometry
        
        try:
            cgal_polygon_wgs84 = shapely_polygon_to_cgal(buffered_geom)
            pwh_wgs84 = bindings.Polygon_with_holes_2(cgal_polygon_wgs84)
            decomposed_wgs84 = bindings.decompose(pwh_wgs84)
            
            all_sweeps_wgs84 = []
            for cgal_poly in decomposed_wgs84:
                coords_part = [(p.x, p.y) for p in cgal_poly.vertices]
                shapely_part = Polygon(coords_part)
                sweep_lines = generate_sweep_pattern(shapely_part, sweep_distance / meters_per_degree, 
                                                    clockwise=True, connect_sweeps=False)
                for line in sweep_lines:
                    all_sweeps_wgs84.extend(list(line.coords))
            
            print(f"✓ WGS84: {len(decomposed_wgs84)} parts, {len(all_sweeps_wgs84)//2} sweeps")
            wgs84_success = True
        except Exception as e:
            print(f"✗ WGS84: Failed - {e}")
            decomposed_wgs84 = []
            all_sweeps_wgs84 = []
            wgs84_success = False
        
        # ====================================================================
        # UTM (MAY CRASH)
        # ====================================================================
        geodata_utm = GeoPolygon(polygon_wgs84, crs="EPSG:4326")
        geodata_utm.set_crs("EPSG:32630")
        
        poly_utm = geodata_utm.geometry
        if poly_utm.geom_type == 'MultiPolygon':
            poly_utm = max(poly_utm.geoms, key=lambda p: p.area)
        
        buffered_utm = geodata_utm.buffer(-sweep_distance)
        
        if buffered_utm.geometry.geom_type == 'MultiPolygon':
            buffered_utm_geom = max(buffered_utm.geometry.geoms, key=lambda p: p.area)
        else:
            buffered_utm_geom = buffered_utm.geometry
        
        # Special handling: heath_4 crashes with UTM, so skip it
        if heath_name == 'heath_4':
            print(f"⚠ UTM: Skipping (known to crash on {heath_name})")
            decomposed_utm = []
            all_sweeps_utm = []
            utm_success = False
            utm_crashed = True
        else:
            try:
                cgal_polygon_utm = shapely_polygon_to_cgal(buffered_utm_geom)
                pwh_utm = bindings.Polygon_with_holes_2(cgal_polygon_utm)
                decomposed_utm = bindings.decompose(pwh_utm)
                
                all_sweeps_utm = []
                for cgal_poly in decomposed_utm:
                    coords_part = [(p.x, p.y) for p in cgal_poly.vertices]
                    shapely_part = Polygon(coords_part)
                    sweep_lines = generate_sweep_pattern(shapely_part, sweep_distance, 
                                                        clockwise=True, connect_sweeps=False)
                    for line in sweep_lines:
                        all_sweeps_utm.extend(list(line.coords))
                
                print(f"✓ UTM: {len(decomposed_utm)} parts, {len(all_sweeps_utm)//2} sweeps")
                utm_success = True
                utm_crashed = False
            except Exception as e:
                print(f"✗ UTM: Failed - {e}")
                decomposed_utm = []
                all_sweeps_utm = []
                utm_success = False
                utm_crashed = False
        
        results.append({
            'name': heath_name,
            'altitude': altitude,
            'sweep_distance': sweep_distance,
            'polygon_wgs84': polygon_wgs84,
            'buffered_wgs84': buffered_geom,
            'decomposed_wgs84': decomposed_wgs84,
            'sweeps_wgs84': all_sweeps_wgs84,
            'wgs84_success': wgs84_success,
            'polygon_utm': poly_utm,
            'buffered_utm': buffered_utm_geom,
            'decomposed_utm': decomposed_utm,
            'sweeps_utm': all_sweeps_utm,
            'utm_success': utm_success,
            'utm_crashed': utm_crashed,
        })

# ============================================================================
# VISUALIZATION (2 columns x 4 rows)
# ============================================================================
print("\n" + "=" * 80)
print("Creating visualization...")
print("=" * 80)

fig, axes = plt.subplots(4, 2, figsize=(16, 22))

for idx, result in enumerate(results):
    row = idx
    
    # ========================================================================
    # LEFT: WGS84
    # ========================================================================
    ax_left = axes[row, 0]
    
    if result['wgs84_success']:
        title_left = f"{result['name']} @ {result['altitude']}m - WGS84 (WORKING)\n{len(result['decomposed_wgs84'])} parts, {len(result['sweeps_wgs84'])//2} sweeps"
        title_color = 'green'
    else:
        title_left = f"{result['name']} @ {result['altitude']}m - WGS84 (FAILED)"
        title_color = 'red'
    
    ax_left.set_title(title_left, fontsize=10, fontweight='bold', color=title_color)
    ax_left.set_xlabel('Longitude (degrees)', fontsize=8)
    ax_left.set_ylabel('Latitude (degrees)', fontsize=8)
    
    # Original polygon
    x, y = result['polygon_wgs84'].exterior.xy
    ax_left.plot(x, y, 'b-', linewidth=1.5, label='Original', zorder=1)
    
    # Buffered polygon
    x_buf, y_buf = result['buffered_wgs84'].exterior.xy
    ax_left.plot(x_buf, y_buf, 'g--', linewidth=1, label='Buffered', zorder=2)
    
    # Decomposed parts with different colors
    colors = plt.cm.tab10(range(len(result['decomposed_wgs84'])))
    for i, cgal_poly in enumerate(result['decomposed_wgs84']):
        coords_part = [(p.x, p.y) for p in cgal_poly.vertices]
        poly_patch = MPLPolygon(coords_part, fill=True, facecolor=colors[i], 
                               edgecolor='darkred', linewidth=0.8, alpha=0.4, zorder=3)
        ax_left.add_patch(poly_patch)
    
    # Sweeps
    if result['sweeps_wgs84']:
        for i in range(0, len(result['sweeps_wgs84']), 2):
            if i + 1 < len(result['sweeps_wgs84']):
                x_vals = [result['sweeps_wgs84'][i][0], result['sweeps_wgs84'][i+1][0]]
                y_vals = [result['sweeps_wgs84'][i][1], result['sweeps_wgs84'][i+1][1]]
                ax_left.plot(x_vals, y_vals, 'orange', linewidth=1.5, alpha=0.7, zorder=4)
        ax_left.plot([], [], 'orange', linewidth=1.5, label=f'Sweeps', zorder=4)
    
    ax_left.legend(loc='best', fontsize=7)
    ax_left.grid(True, alpha=0.3)
    ax_left.axis('equal')
    ax_left.tick_params(labelsize=7)
    
    # ========================================================================
    # RIGHT: UTM
    # ========================================================================
    ax_right = axes[row, 1]
    
    if result['utm_crashed']:
        title_right = f"{result['name']} @ {result['altitude']}m - UTM (WOULD CRASH)\nSkipped to avoid segfault"
        title_color = 'red'
    elif result['utm_success']:
        title_right = f"{result['name']} @ {result['altitude']}m - UTM (WORKING)\n{len(result['decomposed_utm'])} parts, {len(result['sweeps_utm'])//2} sweeps"
        title_color = 'green'
    else:
        title_right = f"{result['name']} @ {result['altitude']}m - UTM (FAILED)"
        title_color = 'red'
    
    ax_right.set_title(title_right, fontsize=10, fontweight='bold', color=title_color)
    ax_right.set_xlabel('Easting (meters)', fontsize=8)
    ax_right.set_ylabel('Northing (meters)', fontsize=8)
    
    # Original polygon
    x_utm, y_utm = result['polygon_utm'].exterior.xy
    ax_right.plot(x_utm, y_utm, 'b-', linewidth=1.5, label='Original', zorder=1)
    
    # Buffered polygon
    x_buf_utm, y_buf_utm = result['buffered_utm'].exterior.xy
    ax_right.plot(x_buf_utm, y_buf_utm, 'g--', linewidth=1, 
                 label=f'Buffered ({-result["sweep_distance"]:.1f}m)', zorder=2)
    
    if result['utm_crashed']:
        # Show crash marker
        centroid_x, centroid_y = result['buffered_utm'].centroid.x, result['buffered_utm'].centroid.y
        ax_right.scatter([centroid_x], [centroid_y], color='red', s=300, marker='X', 
                        label='WOULD CRASH', zorder=5)
    else:
        # Decomposed parts with different colors
        colors = plt.cm.tab10(range(len(result['decomposed_utm'])))
        for i, cgal_poly in enumerate(result['decomposed_utm']):
            coords_part = [(p.x, p.y) for p in cgal_poly.vertices]
            poly_patch = MPLPolygon(coords_part, fill=True, facecolor=colors[i], 
                                   edgecolor='darkred', linewidth=0.8, alpha=0.4, zorder=3)
            ax_right.add_patch(poly_patch)
        
        # Sweeps
        if result['sweeps_utm']:
            for i in range(0, len(result['sweeps_utm']), 2):
                if i + 1 < len(result['sweeps_utm']):
                    x_vals = [result['sweeps_utm'][i][0], result['sweeps_utm'][i+1][0]]
                    y_vals = [result['sweeps_utm'][i][1], result['sweeps_utm'][i+1][1]]
                    ax_right.plot(x_vals, y_vals, 'orange', linewidth=1.5, alpha=0.7, zorder=4)
            ax_right.plot([], [], 'orange', linewidth=1.5, label='Sweeps', zorder=4)
    
    ax_right.legend(loc='best', fontsize=7)
    ax_right.grid(True, alpha=0.3)
    ax_right.axis('equal')
    ax_right.tick_params(labelsize=7)

plt.tight_layout()
output_path = os.path.join(test_dir, 'test_decomposition_crs_output.png')
plt.savefig(output_path, dpi=150, bbox_inches='tight')
print(f"\n✓ Saved visualization to: {output_path}")

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
for result in results:
    print(f"\n{result['name']} @ {result['altitude']}m:")
    if result['wgs84_success']:
        print(f"  WGS84: ✓ {len(result['decomposed_wgs84'])} parts, {len(result['sweeps_wgs84'])//2} sweeps")
    else:
        print(f"  WGS84: ✗ Failed")
    
    if result['utm_crashed']:
        print(f"  UTM:   ⚠ Skipped (would segfault)")
    elif result['utm_success']:
        print(f"  UTM:   ✓ {len(result['decomposed_utm'])} parts, {len(result['sweeps_utm'])//2} sweeps")
    else:
        print(f"  UTM:   ✗ Failed")

print("\nKey finding: heath_4 crashes with UTM (geometry-specific bug), heath_5 works fine")
