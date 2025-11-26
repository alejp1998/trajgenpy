#!/usr/bin/env python3
"""
Test decomposition of heath_4 with different simplification tolerances in UTM.

This test explores whether polygon simplification can avoid the segmentation fault
that occurs with heath_4 in UTM coordinates.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import json
from shapely.geometry import Polygon
from trajgenpy import bindings
import matplotlib.pyplot as plt
import numpy as np
from trajgenpy.Geometries import shapely_polygon_to_cgal, GeoPolygon
import math

# Load test features
test_dir = os.path.dirname(os.path.abspath(__file__))
features_path = os.path.join(test_dir, 'test_features.geojson')

with open(features_path, 'r') as f:
    geojson_data = json.load(f)

# Find heath_4 feature
heath_4 = None
for feature in geojson_data['features']:
    if feature.get('properties', {}).get('unique_name') == 'heath_4':
        heath_4 = feature
        break

if not heath_4:
    print("ERROR: heath_4 not found")
    sys.exit(1)

# Extract polygon
coords = heath_4['geometry']['coordinates'][0]
polygon_wgs84 = Polygon(coords)

print("=" * 80)
print("TEST: Heath_4 Simplification Tolerance Analysis (UTM)")
print("=" * 80)

# Parameters
altitude = 20
field_of_view = 60
overlap = 20
sweep_width = 2 * altitude * math.tan(math.radians(field_of_view / 2))
sweep_distance = sweep_width * (1 - overlap / 100)

print(f"\nOriginal WGS84 polygon:")
print(f"  - Vertices: {len(coords)}")
print(f"  - Area: {polygon_wgs84.area:.10f} sq degrees")

# Convert to UTM
geodata_utm = GeoPolygon(polygon_wgs84, crs="EPSG:4326")
geodata_utm.set_crs("EPSG:32630")

poly_utm = geodata_utm.geometry
if poly_utm.geom_type == 'MultiPolygon':
    poly_utm = max(poly_utm.geoms, key=lambda p: p.area)

print(f"\nConverted to UTM Zone 30N:")
print(f"  - Vertices: {len(list(poly_utm.exterior.coords))}")
print(f"  - Area: {poly_utm.area:.2f} m²")

# Buffer in UTM
buffered_utm = geodata_utm.buffer(-sweep_distance)
if buffered_utm.geometry.geom_type == 'MultiPolygon':
    buffered_utm_geom = max(buffered_utm.geometry.geoms, key=lambda p: p.area)
else:
    buffered_utm_geom = buffered_utm.geometry

original_vertices = len(list(buffered_utm_geom.exterior.coords))
original_area = buffered_utm_geom.area

print(f"\nBuffered polygon:")
print(f"  - Vertices: {original_vertices}")
print(f"  - Area: {original_area:.2f} m²")

# Buffer in WGS84 (for comparison) - convert buffered UTM back to WGS84
from pyproj import Transformer
transformer_to_utm = Transformer.from_crs("EPSG:4326", "EPSG:32630", always_xy=True)
transformer_to_wgs = Transformer.from_crs("EPSG:32630", "EPSG:4326", always_xy=True)

buffered_wgs_coords = [transformer_to_wgs.transform(x, y) for x, y in buffered_utm_geom.exterior.coords]
wgs_geom = Polygon(buffered_wgs_coords)

# Test different tolerance values
# Calculate characteristic length for tolerance scaling
characteristic_length = (original_area ** 0.5)

# Test tolerances from very aggressive to extreme
# Testing if extreme simplification can avoid the crash
tolerance_percentages = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0]

print(f"\nCharacteristic length: {characteristic_length:.2f} m")
print(f"\n{'=' * 80}")
print("Testing different simplification tolerances:")
print('=' * 80)

results = []

for pct in tolerance_percentages:
    tolerance = characteristic_length * (pct / 100)
    
    print(f"\n{'-' * 80}")
    print(f"Tolerance: {pct}% ({tolerance:.4f} m)")
    print('-' * 80)
    
    # Simplify the buffered polygon
    simplified = buffered_utm_geom.simplify(tolerance, preserve_topology=True)
    
    if simplified.is_empty or not simplified.is_valid:
        print(f"  ✗ Simplification produced invalid geometry")
        results.append({
            'tolerance_pct': pct,
            'tolerance_m': tolerance,
            'vertices': 0,
            'area_retained': 0,
            'status': 'invalid',
            'decomposed_parts': 0
        })
        continue
    
    vertices = len(list(simplified.exterior.coords))
    area_retained = (simplified.area / original_area) * 100
    
    print(f"  Simplified: {original_vertices} → {vertices} vertices ({area_retained:.1f}% area)")
    
    # Attempt decomposition
    try:
        cgal_polygon = shapely_polygon_to_cgal(simplified)
        pwh = bindings.Polygon_with_holes_2(cgal_polygon)
        
        print(f"  CGAL objects created, attempting decomposition...")
        decomposed = bindings.decompose(pwh)
        
        print(f"  ✓ SUCCESS: Decomposed into {len(decomposed)} parts")
        results.append({
            'tolerance_pct': pct,
            'tolerance_m': tolerance,
            'vertices': vertices,
            'area_retained': area_retained,
            'status': 'success',
            'decomposed_parts': len(decomposed)
        })
        
    except Exception as e:
        print(f"  ✗ Python exception: {e}")
        results.append({
            'tolerance_pct': pct,
            'tolerance_m': tolerance,
            'vertices': vertices,
            'area_retained': area_retained,
            'status': 'exception',
            'decomposed_parts': 0
        })

print(f"\n{'=' * 80}")
print("RESULTS SUMMARY")
print('=' * 80)
print(f"\n{'Tolerance %':<12} {'Tolerance (m)':<15} {'Vertices':<10} {'Area %':<10} {'Parts':<8} {'Status'}")
print('-' * 80)

for r in results:
    status_symbol = '✓' if r['status'] == 'success' else '✗'
    print(f"{r['tolerance_pct']:<12.1f} {r['tolerance_m']:<15.4f} {r['vertices']:<10} "
          f"{r['area_retained']:<10.1f} {r['decomposed_parts']:<8} {status_symbol} {r['status']}")

# Find minimum tolerance that works
successful = [r for r in results if r['status'] == 'success']
if successful:
    min_working = min(successful, key=lambda x: x['tolerance_pct'])
    print(f"\n✓ Minimum working tolerance: {min_working['tolerance_pct']}% ({min_working['tolerance_m']:.4f} m)")
    print(f"  - Reduces vertices: {original_vertices} → {min_working['vertices']}")
    print(f"  - Retains area: {min_working['area_retained']:.1f}%")
    print(f"  - Produces {min_working['decomposed_parts']} decomposed parts")
else:
    print(f"\n✗ No tolerance level avoided the crash")
    print(f"   Conclusion: Simplification cannot fix this geometry-specific bug")

print(f"\n{'=' * 80}")
print("BINARY SEARCH FOR MINIMUM WORKING TOLERANCE")
print('=' * 80)

if successful:
    # Binary search between 5% (crashes) and min_working (works)
    min_working_pct = min_working['tolerance_pct']
    low, high = 5.0, min_working_pct
    precision = 0.5  # Search with 0.5% precision
    
    print(f"\nSearching between {low}% (crashes) and {high}% (works)...")
    
    best_working = min_working
    
    while high - low > precision:
        mid = (low + high) / 2
        tolerance = characteristic_length * (mid / 100)
        
        print(f"\nTesting {mid:.1f}%... ", end='', flush=True)
        
        simplified = buffered_utm_geom.simplify(tolerance, preserve_topology=True)
        
        if simplified.is_empty or not simplified.is_valid:
            print(f"invalid geometry")
            low = mid
            continue
        
        vertices = len(list(simplified.exterior.coords))
        
        # Write test script to subprocess to handle segfaults
        test_script = f"""
import sys
sys.path.insert(0, '/home/alejp/dev/trajgenpy')
from trajgenpy import bindings
from shapely.geometry import Polygon

simplified_coords = {list(simplified.exterior.coords)}
polygon = Polygon(simplified_coords)

def shapely_polygon_to_cgal(polygon):
    coords = list(polygon.exterior.coords)[:-1]
    cgal_polygon = bindings.Polygon_2()
    for x, y in coords:
        cgal_polygon.push_back(bindings.Point_2(x, y))
    return cgal_polygon

cgal_polygon = shapely_polygon_to_cgal(polygon)
pwh = bindings.Polygon_with_holes_2(cgal_polygon)
decomposed = bindings.decompose(pwh)
print(f"{{len(decomposed)}}")
"""
        
        import subprocess
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
            f.write(test_script)
            test_file = f.name
        
        try:
            result = subprocess.run(
                ['python3', test_file],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            import os
            os.unlink(test_file)
            
            if result.returncode == 0:
                num_parts = int(result.stdout.strip())
                print(f"✓ works ({vertices} vertices, {num_parts} parts)")
                high = mid
                best_working = {
                    'tolerance_pct': mid,
                    'tolerance_m': tolerance,
                    'vertices': vertices,
                    'area_retained': (simplified.area / original_area) * 100,
                    'decomposed_parts': num_parts
                }
            else:
                print(f"✗ crashes ({vertices} vertices)")
                low = mid
        except:
            print(f"✗ crashes ({vertices} vertices)")
            low = mid
            import os
            if os.path.exists(test_file):
                os.unlink(test_file)
    
    print(f"\n{'=' * 80}")
    print("REFINED RESULT")
    print('=' * 80)
    print(f"✓ Minimum working tolerance: {best_working['tolerance_pct']:.1f}% ({best_working['tolerance_m']:.4f} m)")
    print(f"  - Reduces vertices: {original_vertices} → {best_working['vertices']}")
    print(f"  - Retains area: {best_working['area_retained']:.1f}%")
    print(f"  - Produces {best_working['decomposed_parts']} decomposed parts")

print(f"\n{'=' * 80}")
print("CONCLUSION")
print('=' * 80)
if successful:
    print(f"✓ Simplification CAN avoid the crash with heath_4 in UTM coordinates!")
    print(f"  Minimum tolerance: ~{best_working['tolerance_pct']:.1f}% of characteristic length")
    print(f"  For this polygon: {best_working['tolerance_m']:.2f} meters")
    print(f"\n  However, this requires significant geometry modification and may:")
    print(f"  - Lose fine geometric details")
    print(f"  - Reduce coverage accuracy")
    print(f"  - Not work for all polygon shapes")
    print(f"\n  Recommended approach: Use WGS84 coordinates for decomposition")
    print(f"  (no simplification needed, preserves all geometry)")
else:
    print(f"✗ Even aggressive simplification cannot avoid the CGAL crash.")
    print(f"   Recommendation: Use WGS84 coordinates for decomposition instead.")

# ============================================================================
# VISUALIZATION
# ============================================================================
print(f"\n{'=' * 80}")
print("GENERATING VISUALIZATION")
print('=' * 80)

# Select tolerances to visualize (including some that crash)
viz_tolerances = [0, 5.0, 9.1, 9.4, 10.0, 20.0, 30.0, 50.0]
n_cols = 4
n_rows = (len(viz_tolerances) + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 4.5 * n_rows))
fig.subplots_adjust(hspace=0.35, wspace=0.15)
axes = axes.flatten() if n_rows * n_cols > 1 else [axes]

# Helper function to decompose in subprocess
def try_decompose_subprocess(polygon_coords):
    test_script = f"""
import sys
sys.path.insert(0, '/home/alejp/dev/trajgenpy')
from trajgenpy import bindings
from shapely.geometry import Polygon
import json

polygon = Polygon({polygon_coords})

def shapely_polygon_to_cgal(polygon):
    coords = list(polygon.exterior.coords)[:-1]
    cgal_polygon = bindings.Polygon_2()
    for x, y in coords:
        cgal_polygon.push_back(bindings.Point_2(x, y))
    return cgal_polygon

def cgal_polygon_to_shapely(cgal_poly):
    coords = [(cgal_poly.vertex(i).x(), cgal_poly.vertex(i).y()) 
              for i in range(cgal_poly.size())]
    coords.append(coords[0])
    return coords

cgal_polygon = shapely_polygon_to_cgal(polygon)
pwh = bindings.Polygon_with_holes_2(cgal_polygon)
decomposed = bindings.decompose(pwh)

result = []
for poly in decomposed:
    coords = cgal_polygon_to_shapely(poly)
    result.append(coords)

print(json.dumps(result))
"""
    
    import subprocess
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
        f.write(test_script)
        test_file = f.name
    
    try:
        result = subprocess.run(
            ['python3', test_file],
            capture_output=True,
            text=True,
            timeout=5
        )
        
        os.unlink(test_file)
        
        if result.returncode == 0:
            decomposed_coords = json.loads(result.stdout.strip())
            return [Polygon(coords) for coords in decomposed_coords], True
        else:
            return None, False
    except:
        if os.path.exists(test_file):
            os.unlink(test_file)
        return None, False

for idx, tol_pct in enumerate(viz_tolerances):
    ax = axes[idx]
    
    # Always plot the original shape first as reference (light gray)
    x_orig, y_orig = buffered_utm_geom.exterior.xy
    ax.fill(x_orig, y_orig, color='lightgray', alpha=0.3, linewidth=0, label='Original UTM')
    ax.plot(x_orig, y_orig, 'gray', linewidth=1, linestyle='--', alpha=0.5)
    
    if tol_pct == 0:
        # Original buffered polygon
        simplified_utm = buffered_utm_geom
        tolerance_m = 0
        vertices_utm = len(list(simplified_utm.exterior.coords))
        
        # No simplification for WGS84 either
        simplified_wgs = wgs_geom
        vertices_wgs = len(list(simplified_wgs.exterior.coords))
    else:
        tolerance_m = characteristic_length * (tol_pct / 100)
        simplified_utm = buffered_utm_geom.simplify(tolerance_m, preserve_topology=True)
        vertices_utm = len(list(simplified_utm.exterior.coords))
        
        # Also simplify in WGS84 (using degree-equivalent tolerance)
        # Convert tolerance from meters to approximate degrees
        tolerance_deg = tolerance_m / 111000  # Rough conversion at this latitude
        simplified_wgs = wgs_geom.simplify(tolerance_deg, preserve_topology=True)
        vertices_wgs = len(list(simplified_wgs.exterior.coords))
    
    # Try to decompose UTM version
    decomposed_utm, utm_success = try_decompose_subprocess(list(simplified_utm.exterior.coords))
    
    # Try to decompose WGS84 version
    decomposed_wgs, wgs_success = try_decompose_subprocess(list(simplified_wgs.exterior.coords))
    
    # Track number of parts for status (only count if success is True)
    utm_parts = len(decomposed_utm) if (utm_success and decomposed_utm) else 0
    wgs_parts = len(decomposed_wgs) if (wgs_success and decomposed_wgs) else 0
    
    # Plot UTM decomposed parts with cool colors (blues/greens)
    if utm_success and decomposed_utm:
        colors_utm = plt.cm.Blues(np.linspace(0.4, 0.9, max(len(decomposed_utm), 3)))
        for i, part in enumerate(decomposed_utm):
            x_part, y_part = part.exterior.xy
            ax.fill(x_part, y_part, color=colors_utm[i % len(colors_utm)], alpha=0.4, 
                   edgecolor='blue', linewidth=1.5, linestyle='-')
    elif not utm_success:
        # Show crash with red X pattern
        x_utm, y_utm = simplified_utm.exterior.xy
        ax.fill(x_utm, y_utm, color='red', alpha=0.15)
    
    # Plot WGS84 decomposed parts with warm colors (reds/oranges) - overlaid
    if wgs_success and decomposed_wgs and tol_pct > 0:
        # Transform WGS84 back to UTM for plotting
        from pyproj import Transformer
        transformer_to_utm = Transformer.from_crs("EPSG:4326", "EPSG:32630", always_xy=True)
        
        colors_wgs = plt.cm.Oranges(np.linspace(0.4, 0.9, max(len(decomposed_wgs), 3)))
        for i, part in enumerate(decomposed_wgs):
            # Transform to UTM for plotting
            wgs_coords = list(part.exterior.coords)
            utm_coords = [transformer_to_utm.transform(x, y) for x, y in wgs_coords]
            utm_x, utm_y = zip(*utm_coords)
            ax.fill(utm_x, utm_y, color=colors_wgs[i % len(colors_wgs)], alpha=0.3, 
                   edgecolor='darkorange', linewidth=1.5, linestyle='--')
    
    # Plot the simplified UTM polygon outline
    x_utm, y_utm = simplified_utm.exterior.xy
    line_style = '-' if utm_success else ':'
    ax.plot(x_utm, y_utm, 'b', linewidth=2, linestyle=line_style, alpha=0.8)
    
    # Plot WGS84 simplified outline (transformed to UTM)
    if tol_pct > 0:
        wgs_coords = list(simplified_wgs.exterior.coords)
        utm_coords_wgs = [transformer_to_utm.transform(x, y) for x, y in wgs_coords]
        utm_x_wgs, utm_y_wgs = zip(*utm_coords_wgs)
        line_style_wgs = '--' if wgs_success else ':'
        ax.plot(utm_x_wgs, utm_y_wgs, 'darkorange', linewidth=2, linestyle=line_style_wgs, alpha=0.8)
    
    # Set title and labels
    if tol_pct == 0:
        ax.set_title(f"Original Buffered\nUTM:{vertices_utm}v WGS:{vertices_wgs}v", 
                    fontsize=12, fontweight='bold')
    else:
        # Build status string
        utm_status = f"✓{utm_parts}p" if utm_success else "✗CRASH"
        wgs_status = f"✓{wgs_parts}p" if wgs_success else "✗CRASH"
        ax.set_title(f"Tolerance: {tol_pct}% ({tolerance_m:.2f}m)\nUTM:{vertices_utm}v WGS:{vertices_wgs}v | UTM:{utm_status} WGS:{wgs_status}", 
                    fontsize=11, fontweight='bold')
    
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('Easting (m)', fontsize=10)
    ax.set_ylabel('Northing (m)', fontsize=10)
    ax.ticklabel_format(style='plain', useOffset=False)
    
    # Add legend only to first plot
    if idx == 0:
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        legend_elements = [
            Patch(facecolor='lightgray', alpha=0.3, edgecolor='gray', linestyle='--', 
                  label='Original (44v)'),
            Line2D([0], [0], color='blue', linewidth=2, linestyle='-', label='UTM boundary'),
            Patch(facecolor='blue', alpha=0.4, edgecolor='blue', label='UTM decomposed'),
            Line2D([0], [0], color='darkorange', linewidth=2, linestyle='--', label='WGS84 boundary'),
            Patch(facecolor='orange', alpha=0.3, edgecolor='darkorange', linestyle='--', 
                  label='WGS84 decomposed'),
            Patch(facecolor='red', alpha=0.15, label='Crash indicator'),
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=8, framealpha=0.9)

# Hide unused subplots
for idx in range(len(viz_tolerances), len(axes)):
    axes[idx].axis('off')

plt.suptitle('Heath_4 Simplification Tolerance Effects on Decomposition (UTM Zone 30N)', 
            fontsize=16, fontweight='bold')

output_file = os.path.join(os.path.dirname(__file__), 'test_simplification_tolerance_viz.png')
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"\n✓ Visualization saved to: {output_file}")
print(f"  File size: {os.path.getsize(output_file) / 1024:.1f} KB")
