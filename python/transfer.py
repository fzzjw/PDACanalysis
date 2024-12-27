import openslide
import numpy as np
import os
import re


def write_txt(image, output_file, tile_resoluntion=400, downsample=4):
    """
    Read WSI image and generate a file with tile starting coordinates.
    :param image: Path to the WSI image.
    :param output_file: Path to the output file for tile coordinates.
    :param tile_resoluntion: Resolution of each tile in pixels (default is 400).
    :param downsample: Downsampling factor (default is 4).
    """
    try:
        # Open the WSI image; origin is top-left with x-axis to the right and y-axis downward.
        slide = openslide.open_slide(image)
        column0, row0 = slide.level_dimensions[0]
        step = tile_resoluntion * downsample  # Step size for tiles
        number_row = row0 // step  # Number of rows (Y-axis)
        number_column = column0 // step  # Number of columns (X-axis)

        # Write tile coordinates to the output file
        with open(output_file, 'w') as output_file:
            output_file.write('x_coordinates\ty_coordinates\n')
            output_file.write('---------------------------------------\n')
            for j in range(number_column):
                for i in range(number_row):
                    x = j * step
                    y = i * step
                    output_file.write(f'{x}\t{y}\n')

        print(f"Coordinates have been written to {output_file.name} successfully.")
        print(f"Number of columns: {number_column}")
        print(f"Number of rows: {number_row}")

    except Exception as e:
        print(f"An error occurred: {e}")


def rename_txt_files(sub_downsample_file, local_coord_path):
    """
    Rename local coordinate files based on the starting coordinates.
    :param sub_downsample_file: Path to the file containing sub-region starting coordinates.
    :param local_coord_path: Directory containing local coordinate files for sub-regions.
    """
    # Load the sub-region starting coordinates
    tile_data = np.loadtxt(sub_downsample_file, skiprows=2)

    # Combine x and y coordinates into tuples
    coordinates = list(map(tuple, tile_data))
    
    # Get all .txt files in the directory
    txt_files = [f for f in os.listdir(local_coord_path) if f.endswith('.txt')]
    
    # Template for renaming files
    file_template = "file_{}_x-{}_y-{}.txt"

    # Regex to extract coordinates from filenames
    pattern = re.compile(r'x-(\d+)_y-(\d+)')
    
    # Map coordinates to filenames
    coordinate_to_file = {}
    for filename in txt_files:
        match = pattern.search(filename)
        if match:
            x = float(match.group(1))
            y = float(match.group(2))
            coordinate_to_file[(x, y)] = filename

    # Rename files based on coordinates
    for idx, (x, y) in enumerate(coordinates, start=1):
        if (x, y) in coordinate_to_file:
            old_filename = coordinate_to_file[(x, y)]
            new_filename = file_template.format(idx, int(x), int(y))
            os.rename(os.path.join(local_coord_path, old_filename), os.path.join(local_coord_path, new_filename))
            print(f"Renamed {old_filename} to {new_filename}")
        else:
            print(f"No file found for coordinates ({x}, {y})")
