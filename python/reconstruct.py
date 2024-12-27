import numpy as np
import os
import matplotlib.pyplot as plt
import transfer
import openslide


def load_sub_origin(file_path):
    """
    Load starting coordinates for each tile.
    :param file_path: Path to the txt file.
    :return: An Nx2 numpy array with starting coordinates (x, y) for each tile.
    """
    return np.loadtxt(file_path, delimiter='\t', skiprows=2)


def convert_to_global_coordinates(local_coords, pixel_size_um):
    """
    Convert local coordinates to global coordinates (downsampled).
    :param local_coords: Local coordinate points (Nx2 numpy array).
    :param pixel_size_um: Pixel resolution of the tile in µm/pixel.
    :return: Global coordinate points (Nx2 numpy array, downsampled by 4x).
    """
    # Convert local coordinates from µm to pixels
    global_coords_pixels_downsample = local_coords / pixel_size_um
    return global_coords_pixels_downsample


def reconstruct_coordinates(sub_downsample_file, coords_folder, pixel_size_um):
    """
    Read and convert sub-region coordinates to global coordinates.
    :param sub_downsample_file: Path to the file containing sub-region starting coordinates.
    :param coords_folder: Folder containing sub-region point coordinate files.
    :param pixel_size_um: Pixel resolution of the original WSI in µm/pixel.
    :return: All converted coordinate points (Nx2 numpy array).
    """
    tile_coordinate_data = load_sub_origin(sub_downsample_file)
    coordinates = list(map(tuple, tile_coordinate_data))
    point_all = []

    for idx, (x, y) in enumerate(coordinates, start=1):  # Start from 1
        filename = os.path.join(coords_folder, f'file_{idx}_x-{int(x)}_y-{int(y)}.txt')

        # Check if the file exists and is not empty
        if os.path.exists(filename):
            with open(filename, 'r') as file:
                lines = file.readlines()
                if len(lines) > 1:
                    try:
                        tile_data = np.loadtxt(filename, usecols=(1, 2), skiprows=1)
                        if tile_data.ndim == 1:  # Ensure tile_data is 2D
                            tile_data = np.expand_dims(tile_data, axis=0)
                        
                        global_downsample = convert_to_global_coordinates(tile_data, pixel_size_um)
                        global_downsample[:, 0] += tile_coordinate_data[idx - 1, 0]  # Adjust x coordinates
                        global_downsample[:, 1] += tile_coordinate_data[idx - 1, 1]  # Adjust y coordinates
                        global_sample_coords = global_downsample

                        point_all.extend(global_sample_coords)
                        print(f'Processed file: file_{idx}_x-{int(x)}_y-{int(y)}.txt')
                    except Exception as e:
                        print(f'Error reading {filename}: {e}')
                else:
                    print(f'File is empty or only contains a header: {filename}')
        else:
            print(f'File does not exist: {filename}')

    return np.array(point_all)


def save_global_coordinates(file_path, global_sample_coords):
    """
    Save global coordinate points to a txt file.
    :param file_path: Path to the txt file.
    :param global_sample_coords: Global coordinate points (Nx2 numpy array).
    """
    np.savetxt(file_path, global_sample_coords, delimiter='\t', comments='')


def plot_global_coordinates(file_path, img):
    """
    Read the global_coordinates.txt file and plot the coordinates.
    :param file_path: Path to the coordinate file.
    :param img: Path to the image file.
    """
    slide = openslide.OpenSlide(img)
    ori_width, ori_height = slide.dimensions

    # Calculate image dimensions in µm
    width_um = ori_width * 0.205
    height_um = ori_height * 0.205

    # Load coordinate data
    coordinates = np.loadtxt(file_path)
    x_coords = coordinates[:, 0]
    y_coords = coordinates[:, 1]

    # Create the plot
    plt.figure(figsize=(ori_width / 1000, ori_height / 1000))
    plt.scatter(x_coords, y_coords)

    plt.xlabel('X (µm)', fontsize=80)
    plt.ylabel('Y (µm)', fontsize=80)

    # Set custom ticks (in µm)
    num_ticks = 10
    x_pixel_ticks = np.linspace(0, ori_width, num=num_ticks)
    y_pixel_ticks = np.linspace(0, ori_height, num=num_ticks)
    x_um_ticks = x_pixel_ticks * 0.205
    y_um_ticks = y_pixel_ticks * 0.205

    plt.xticks(x_pixel_ticks, labels=[f"{tick:.0f}" for tick in x_um_ticks], fontsize=80)
    plt.yticks(y_pixel_ticks, labels=[f"{tick:.0f}" for tick in y_um_ticks], fontsize=80)

    # Invert y-axis for proper visualization
    plt.gca().invert_yaxis()

    # Adjust border thickness
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(5)

    # Display the plot
    plt.show()


file_name = '' #file_name

input_path = '../CD20/' + file_name + '/local_coord' 

coords_folder = input_path 

img = '../CD20/' + file_name + '/' + file_name + '.svs' 

sub_downsample_file = output_file = '../CD20/' + file_name + '/input/sub_downsample_CD20.txt' 

plt.rcParams['font.family'] = 'Arial'

transfer.write_txt(img, output_file)

transfer.rename_txt_files(sub_downsample_file, coords_folder)


pixel_size_um = 0.205
downsample = 4



global_coord_file = '../CD20/' + file_name + '/points/global_coordinates.txt' 

global_sample_coords = reconstruct_coordinates(sub_downsample_file, coords_folder, pixel_size_um)

save_global_coordinates(global_coord_file, global_sample_coords)
print("success")