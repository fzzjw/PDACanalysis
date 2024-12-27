import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import hdbscan
import matplotlib.patches as patches
from matplotlib.patches import Ellipse
from scipy.stats import chi2
from descartes import PolygonPatch
import alphashape
from shapely.geometry import Polygon
from shapely.plotting import plot_polygon
from shapely.geometry import Point, MultiPoint
from concave_hull import concave_hull, concave_hull_indexes
from scipy.spatial import ConvexHull
from sklearn.covariance import EllipticEnvelope
import openslide

def read_coords_from_txt(filename):
    """
    Read coordinates from a TXT file.
    :param filename: Path to the TXT file.
    :return: List of coordinates as tuples (x, y).
    """
    coords = []
    with open(filename, 'r') as file:
        for line in file:
            x, y = line.strip().split('\t')
            coords.append((float(x), float(y)))
    return coords


def calculate_numbers(centroid, normal_tissue, invasive_front, cancer_center, cancer_island):
    """
    Calculate the number of centroids in each region.
    :param centroid: Coordinates of the centroid (x, y).
    :param normal_tissue: Coordinates of the normal tissue region.
    :param invasive_front: Coordinates of the invasive front region.
    :param cancer_center: Coordinates of the cancer center region.
    :param cancer_island: Coordinates of the cancer island region.
    :return: Tuple with counts of centroids in each region.
    """
    def count_centroids_in_region(centroid, region_coords):
        """
        Check if the centroid is within the given region.
        :param centroid: Coordinates of the centroid (x, y).
        :param region_coords: Coordinates defining the region.
        :return: 1 if the centroid is within the region, otherwise 0.
        """
        region_polygon = Polygon(region_coords)
        return 1 if region_polygon.contains(Point(centroid)) else 0

    # Count centroids in each region
    num_in_normal = count_centroids_in_region(centroid, normal_tissue)
    num_in_invasive_front = count_centroids_in_region(centroid, invasive_front)
    num_in_cancer_center = count_centroids_in_region(centroid, cancer_center)
    num_in_cancer_island = count_centroids_in_region(centroid, cancer_island)

    return num_in_normal, num_in_invasive_front, num_in_cancer_center, num_in_cancer_island



file_name = '' #filename
image_path = '../CD20/' + file_name + '/' + file_name + '.svs'
min_cluster_size=30
min_samples=15
pixel_to_um = 0.205

global_coord_file = '../CD20/' + file_name + '/points/global_coordinates.txt'

output_file = '../CD20/' + file_name + '/output_file/cluster_info.txt'
statistic_file = '../CD20/' + file_name + '/output_file/cluster_statistic.xlsx'

filename_if = '../CD20/' + file_name + '/contours/' + 'Invasive_front_coords.txt'
filename_ct = '../CD20/' + file_name + '/contours/' + 'CT_coords.txt'
filename_n = '../CD20/' + file_name + '/contours/' + 'Normal_coords.txt'
filename_ci = '../CD20/' + file_name + '/contours/' + 'Cancer_island_coords.txt'
filename_whole = '../CD20/' + file_name + '/contours/' + 'Whole_tissue_coords.txt'

# Read region coordinates
invasive_front = read_coords_from_txt(filename_if)
cancer_center = read_coords_from_txt(filename_ct)
normal_tissue = read_coords_from_txt(filename_n)
cancer_island = read_coords_from_txt(filename_ci)
whole_tissue = read_coords_from_txt(filename_whole)

# Calculate pixel area for regions
if_area = Polygon(invasive_front)
ct_area = Polygon(cancer_center)
n_area = Polygon(normal_tissue)
ci_area = Polygon(cancer_island)
w_area = Polygon(whole_tissue)

# Convert pixel area to mm²
if_area_mm = if_area.area * pixel_to_um**2 * 1e-6
ct_area_mm = ct_area.area * pixel_to_um**2 * 1e-6
ci_area_mm = ci_area.area * pixel_to_um**2 * 1e-6
n_area_mm = n_area.area * pixel_to_um**2 * 1e-6
w_area_mm = w_area.area * pixel_to_um**2 * 1e-6

# Load global coordinates
coordinates = np.loadtxt(global_coord_file, delimiter='\t')
coordinates_um = coordinates

# Perform clustering using HDBSCAN
clusterer = hdbscan.HDBSCAN(min_cluster_size, min_samples)
cluster_labels = clusterer.fit_predict(coordinates_um)

# Count unique labels and noise points
unique_labels = set(cluster_labels)
num_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)
num_noise_points = list(cluster_labels).count(-1)

print(f'Number of clusters: {num_clusters}')
print(f'Number of noise points: {num_noise_points}')

# Initialize result storage
cluster_info = []
cluster_statistic = []

# Initialize counters for clusters and areas within regions
if_cl_num = ct_cl_num = ci_cl_num = n_cl_num = 0
if_cl_area = ct_cl_area = ci_cl_area = n_cl_area = 0

# Calculate clusters per unit area
if_cl_num_per = ct_cl_num_per = ci_cl_num_per = n_cl_num_per = 0
if_cl_area_per = ct_cl_area_per = ci_cl_area_per = n_cl_area_per = 0

# Process each cluster
for label in unique_labels:
    if label == -1:
        continue  # Skip noise points

    # Get cell coordinates for the current cluster
    class_member_mask = (cluster_labels == label)
    cluster_coords = coordinates_um[class_member_mask]
    num_cells = len(cluster_coords)

    # Compute concave hull and area
    idxes = concave_hull_indexes(cluster_coords[:, :2], length_threshold=50)
    hull_coords = cluster_coords[idxes]
    polygon = Polygon(hull_coords)
    cluster_area = polygon.area
    cluster_area_mm = cluster_area * pixel_to_um**2 * 1e-6

    # Compute density and centroid
    density = num_cells / cluster_area_mm
    centroid = np.mean(cluster_coords, axis=0)

    # Count the centroid in different regions
    num_in_normal, num_in_invasive_front, num_in_cancer_center, num_in_cancer_island = calculate_numbers(
        centroid, normal_tissue, invasive_front, cancer_center, cancer_island)

    # Update region-specific cluster counts and areas
    if num_in_normal:
        n_cl_num += 1
        n_cl_area += cluster_area_mm
    if num_in_invasive_front:
        if_cl_num += 1
        if_cl_area += cluster_area_mm
    if num_in_cancer_center:
        ct_cl_num += 1
        ct_cl_area += cluster_area_mm
    if num_in_cancer_island:
        ci_cl_num += 1
        ci_cl_area += cluster_area_mm

    # Append cluster info
    cluster_info.append([label, num_cells, cluster_area_mm, density, centroid[0], centroid[1],
                         num_in_normal, num_in_invasive_front, num_in_cancer_center, num_in_cancer_island])

# Compute per unit area metrics
if_cl_num_per = if_cl_num / if_area_mm
ct_cl_num_per = ct_cl_num / ct_area_mm
ci_cl_num_per = ci_cl_num / ci_area_mm
n_cl_num_per = n_cl_num / n_area_mm

if_cl_area_per = if_cl_area / if_area_mm
ct_cl_area_per = ct_cl_area / ct_area_mm
ci_cl_area_per = ci_cl_area / ci_area_mm
n_cl_area_per = n_cl_area / n_area_mm

# Append statistics
cluster_statistic.append([n_cl_num_per, n_cl_area_per,
                          if_cl_num_per, if_cl_area_per,
                          ct_cl_num_per, ct_cl_area_per,
                          ci_cl_num_per, ci_cl_area_per])

# Save basic information to a TXT file
with open(output_file, 'w') as f:
    f.write("ClusterID CellCount ClusterArea(mm^2) CellDensity(/mm^2) ClusterCentroid_X ClusterCentroid_Y NormalTissueCount InvasiveFrontCount CancerCenterCount CancerIslandCount\n")
    for info in cluster_info:
        f.write(f"{info[0]} {info[1]} {info[2]:.2f} {info[3]:.6f} {info[4]:.2f} {info[5]:.2f} {info[6]} {info[7]} {info[8]} {info[9]}\n")

# Save summarized statistics to an Excel file
df = pd.DataFrame(cluster_statistic, columns=["N_cluster_num_per", "N_cluster_area_per",
                                              "IF_cluster_num_per", "IF_cluster_area_per",
                                              "CT_cluster_num_per", "CT_cluster_area_per",
                                              "CI_cluster_num_per", "CI_cluster_area_per"])

df.to_excel(statistic_file, index=False)

print(f"IF has {if_cl_num} clusters")
print(f"CT has {ct_cl_num} clusters")
print(f"CI has {ci_cl_num} clusters")
print(f"N has {n_cl_num} clusters")
