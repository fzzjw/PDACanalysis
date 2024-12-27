import numpy as np
import xml.etree.ElementTree as ET
from shapely.geometry import Polygon, LineString, Point
import openslide
import cv2
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

def read_xml_and_get_coordinates(xml_file):
    """
    Read XML file and extract contour data and coordinates.
    :param xml_file: Path to the XML file.
    :return: Two lists - outer points and line points.
    """
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
    except ET.ParseError as e:
        print(f"Error parsing XML: {e}")
        return [], []

    outer_points = []
    if_line = []

    for annotation in root.findall('.//Annotation'):
        annotation_type = annotation.get('Type')
        for coordinate in annotation.findall('.//Coordinate'):
            x = float(coordinate.get('X'))
            y = float(coordinate.get('Y'))
            if annotation_type == "Spline":
                outer_points.append((x, y))
            elif annotation_type in ["PointSet", "Dot"]:
                if_line.append((x, y))

    return outer_points, if_line


def create_mask_if(image_path, annotations, level=4):
    """
    Create a mask for the invasive front (IF).
    :param image_path: Path to the SVS file.
    :param annotations: Annotation coordinates.
    :param level: Downsample level (default is 4x).
    """
    slide = openslide.OpenSlide(image_path)
    downsample = slide.level_downsamples[level]
    dims = slide.level_dimensions[level]
    mask_if = np.zeros(dims[::-1], dtype=np.uint8)  # Reverse dims for compatibility with numpy
    annotations = [annotations]

    for points in annotations:
        scaled_points = np.array([(x / downsample, y / downsample) for x, y in points], dtype=np.float32)
        cv2.polylines(mask_if, [scaled_points.astype(np.int32)], isClosed=False, color=255, thickness=1)

    return mask_if, downsample


def create_mask_whole_contour(image_path, annotations, level=4):
    """
    Create a mask for the whole tissue.
    :param image_path: Path to the SVS file.
    :param annotations: Annotation coordinates.
    :param level: Downsample level (default is 4x).
    """
    slide = openslide.OpenSlide(image_path)
    downsample = slide.level_downsamples[level]
    dims = slide.level_dimensions[level]
    mask_whole_tissue = np.zeros(dims[::-1], dtype=np.uint8)  # Reverse dims for compatibility with numpy
    annotations = [annotations]

    for points in annotations:
        scaled_points = np.array([(x / downsample, y / downsample) for x, y in points], dtype=np.float32)
        cv2.fillPoly(mask_whole_tissue, [scaled_points.astype(np.int32)], 255)

    return mask_whole_tissue, downsample


def extract_coordinates_from_mask(mask, downsample):
    """
    Extract coordinates from an IF mask.
    :param mask: IF mask as a binary image.
    :param downsample: Downsample factor.
    :return: List of coordinates.
    """
    mask_array = np.array(mask, dtype=np.uint8)
    contours, _ = cv2.findContours(mask_array, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    coords = []
    for contour in contours:
        for point in contour:
            x, y = point[0]
            coords.append((x * downsample + downsample / 2, y * downsample + downsample / 2))

    return coords


def create_if_buffer(if_line_coords, pixel_to_micron_ratio, offset_mm):
    """
    Generate invasive front buffer region.
    :param pixel_to_micron_ratio: Microns per pixel (e.g., 0.205 µm/pixel).
    :param offset_mm: Buffer distance in millimeters.
    """
    offset_pixels = offset_mm / (pixel_to_micron_ratio / 1000)
    line = LineString(if_line_coords)
    outer_polygon = line.buffer(offset_pixels)

    if_buffer = []
    if outer_polygon.is_valid:
        if outer_polygon.geom_type == 'Polygon':
            if_buffer = list(outer_polygon.exterior.coords)
        elif outer_polygon.geom_type == 'MultiPolygon':
            for poly in outer_polygon:
                if_buffer.extend(list(poly.exterior.coords))

    return if_buffer


def create_CT_and_N_with_invasive_front(whole_tissue_contour_coords, if_buffer):
    """
    Generate normal and tumor center regions.
    :param whole_tissue_contour_coords: Whole tissue contour coordinates.
    :param if_buffer: IF region coordinates.
    """
    whole_tissue_contour = Polygon(whole_tissue_contour_coords)
    invasive_front = Polygon(if_buffer)
    split_result = whole_tissue_contour.difference(invasive_front)

    split_coords = list(split_result.geoms)
    up_coords = split_coords[0]
    down_coords = split_coords[1]

    return list(up_coords.exterior.coords), list(down_coords.exterior.coords)


def create_overlay_plot(image_path, red_coords, yellow_coords, green_coords, level=4):
    """
    Visualize regions with an overlay.
    :param image_path: Path to the image file.
    :param red_coords: IF region coordinates.
    :param yellow_coords: Tumor center region coordinates.
    :param green_coords: Normal region coordinates.
    :param level: Downsample level.
    """
    slide = openslide.OpenSlide(image_path)
    downsample = slide.level_downsamples[level]
    dims = slide.level_dimensions[level]

    pil_image = slide.read_region((0, 0), level, dims)
    image_rgb = np.array(pil_image.convert('RGB'))

    fig, ax = plt.subplots(figsize=(dims[0] / 100, dims[1] / 100), dpi=100)
    ax.imshow(image_rgb)

    colors = {
        'bright_yellow': (1, 1, 0, 0.3),
        'peach_orange': (1.0, 0.8, 0.6, 0.3),
        'moss_green': (0.13, 0.37, 0.31, 0.3)
    }

    def add_overlay(coords, color):
        for points in coords:
            scaled_points = [(int(x / downsample), int(y / downsample)) for x, y in points]
            polygon = patches.Polygon(scaled_points, closed=True, color=colors[color])
            ax.add_patch(polygon)

    add_overlay([red_coords], 'bright_yellow')
    add_overlay([yellow_coords], 'peach_orange')
    add_overlay([green_coords], 'moss_green')

    x_pixel_ticks = np.linspace(0, dims[0], num=10)
    y_pixel_ticks = np.linspace(0, dims[1], num=10)
    x_um_ticks = x_pixel_ticks * downsample * 0.205
    y_um_ticks = y_pixel_ticks * downsample * 0.205

    plt.xticks(x_pixel_ticks, labels=[f"{tick:.0f}" for tick in x_um_ticks], fontsize=40)
    plt.yticks(y_pixel_ticks, labels=[f"{tick:.0f}" for tick in y_um_ticks], fontsize=40)

    ax.set_xlabel('X, µm', fontsize=40)
    ax.set_ylabel('Y, µm', fontsize=40)
    plt.show()


def create_invasive_front(whole_tissue_contour_coords, if_buffer):
    """
    Generate the invasive front region.
    :param whole_tissue_contour_coords: Whole tissue contour coordinates.
    :param if_buffer: IF region coordinates.
    """
    whole_tissue_polygon = Polygon(whole_tissue_contour_coords)
    if_polygon = Polygon(if_buffer)
    intersection = whole_tissue_polygon.intersection(if_polygon)

    invasive_front_coords = []
    if intersection.is_valid:
        if intersection.geom_type == 'Polygon':
            invasive_front_coords = list(intersection.exterior.coords)
        elif intersection.geom_type == 'MultiPolygon':
            for poly in intersection:
                invasive_front_coords.extend(list(poly.exterior.coords))

    return invasive_front_coords


def create_cancer_island(invasive_front_coords, center_tumor_coords):
    """
    Merge invasive front and tumor center regions into a single polygon.
    :param invasive_front_coords: IF region coordinates.
    :param center_tumor_coords: Tumor center region coordinates.
    """
    invasive_front_polygon = Polygon(invasive_front_coords)
    center_tumor_polygon = Polygon(center_tumor_coords)
    merged_polygon = invasive_front_polygon.union(center_tumor_polygon)

    cancer_island = []
    if merged_polygon.is_valid:
        if merged_polygon.geom_type == 'Polygon':
            cancer_island = list(merged_polygon.exterior.coords)
        elif merged_polygon.geom_type == 'MultiPolygon':
            for poly in merged_polygon:
                cancer_island.extend(list(poly.exterior.coords))

    return cancer_island


def save_coords_to_txt(coords, filename):
    """
    Save coordinates to a TXT file.
    :param coords: List of coordinates.
    :param filename: File name.
    """
    with open(filename, 'w') as file:
        for coord in coords:
            file.write(f"{coord[0]}\t{coord[1]}\n")


def read_coords_from_txt(filename):
    """
    Read coordinates from a TXT file.
    :param filename: File name.
    """
    coords = []
    with open(filename, 'r') as file:
        for line in file:
            x, y = line.strip().split('\t')
            coords.append((float(x), float(y)))
    return coords

file_name = ''  #file_name

xml_file ='../CD20/' + file_name + '/' + file_name + '.xml'

img = '../CD20/' + file_name + '/' + file_name + '.svs'

filename_if = '../CD20/' + file_name + '/contours/' + 'Invasive_front_coords.txt'
filename_ct = '../CD20/' + file_name + '/contours/' + 'CT_coords.txt'
filename_n = '../CD20/'  + file_name + '/contours/' + 'Normal_coords.txt'
filename_ci = '../CD20/' + file_name + '/contours/' + 'Cancer_island_coords.txt'
filename_whole = '../CD20/' + file_name + '/contours/' + 'Whole_tissue_coords.txt'




pixel_to_micron_ratio, offset_mm = 0.205, 0.5
outer_points, if_line = read_xml_and_get_coordinates(xml_file)

mask_whole_tissue, downsample = create_mask_whole_contour(img, outer_points)
mask_if, downsample = create_mask_if(img, if_line)

if_line_coords = extract_coordinates_from_mask(mask_if, downsample)
whole_tissue_contour_coords = extract_coordinates_from_mask(mask_whole_tissue, downsample)

if_buffer = create_if_buffer(if_line_coords, pixel_to_micron_ratio, offset_mm)
up_contour_coords, down_contour_coords = create_CT_and_N_with_invasive_front(whole_tissue_contour_coords, if_buffer)
invasive_front_coords = create_invasive_front(whole_tissue_contour_coords, if_buffer)

'''
if_buffer: IF buffer region (expanded area around the invasive front)
up_contour_coords, down_contour_coords: Normal tissue or tumor center regions, the exact type depends on visualization
invasive_front_coords: IF region along the boundary line
cancer_island: Combined region of IF and CT (center tumor area)
'''
cancer_island = create_cancer_island(invasive_front_coords, down_contour_coords) # if + ct

create_overlay_plot(img, invasive_front_coords, up_contour_coords, down_contour_coords)

save_coords_to_txt(invasive_front_coords, filename_if)

save_coords_to_txt(down_contour_coords, filename_ct) 

save_coords_to_txt(up_contour_coords, filename_n) 

save_coords_to_txt(cancer_island, filename_ci)
save_coords_to_txt(whole_tissue_contour_coords, filename_whole)

invasive_front = read_coords_from_txt(filename_if)
cancer_center = read_coords_from_txt(filename_ct)
normal_tissue = read_coords_from_txt(filename_n)
cancer_island = read_coords_from_txt(filename_ci)
whole_tissue = read_coords_from_txt(filename_whole)