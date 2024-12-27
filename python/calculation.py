from get_area import *
from shapely.geometry import Point, Polygon
import numpy as np
import openslide
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed


def calculate_area_occupancy_rate(image_path, txt_file, min_um_size, middle_um_size, max_um_size, step_size_small,
                                  step_size_large, area_coords, pixel_to_um):
    """
    Calculate the occupancy rate of a given area in an image.
    :param image_path: Path to the SVS image.
    :param txt_file: Path to the TXT file with coordinate points.
    :param min_um_size: Minimum grid size in microns.
    :param middle_um_size: Middle grid size in microns.
    :param max_um_size: Maximum grid size in microns.
    :param step_size_small: Step size for small scales in microns.
    :param step_size_large: Step size for large scales in microns.
    :param area_coords: Coordinates defining the area of interest.
    :param pixel_to_um: Conversion factor from pixels to microns.
    :return: Lists of occupancy rates, occupied counts, and total grid counts.
    """
    slide = openslide.OpenSlide(image_path)
    width, height = slide.level_dimensions[0]

    # Convert grid sizes and step sizes from microns to pixels
    min_grid_size = min_um_size / pixel_to_um
    middle_grid_size = middle_um_size / pixel_to_um
    max_grid_size = max_um_size / pixel_to_um
    step_size_small = step_size_small / pixel_to_um
    step_size_large = step_size_large / pixel_to_um

    # Load coordinate points
    points = np.loadtxt(txt_file)

    # Create a polygon for the irregular area
    polygon = Polygon(area_coords)
    minx, miny, maxx, maxy = polygon.bounds

    results = []

    def process_grid_size(grid_size, step_size):
        """
        Process a single grid size to calculate the occupancy rate.
        :param grid_size: Size of the grid in pixels.
        :param step_size: Step size for traversing the grid.
        """
        num_occupied = 0
        num_grids = 0

        x_ranges = np.arange(minx, maxx, grid_size)
        y_ranges = np.arange(miny, maxy, grid_size)

        x_points = points[:, 0]
        y_points = points[:, 1]

        for y in y_ranges:
            for x in x_ranges:
                x_end = min(x + grid_size, width)
                y_end = min(y + grid_size, height)

                grid_polygon = Polygon([(x, y), (x_end, y), (x_end, y_end), (x, y_end)])
                if polygon.intersects(grid_polygon):
                    in_x_range = (x_points >= x) & (x_points < x_end)
                    in_y_range = (y_points >= y) & (y_points < y_end)
                    points_in_grid = points[in_x_range & in_y_range]

                    if points_in_grid.size > 0:
                        num_occupied += 1
                    num_grids += 1

        occupancy_rate = num_occupied / num_grids if num_grids > 0 else 0
        return grid_size, occupancy_rate, num_occupied, num_grids

    # Process grids in parallel
    with ThreadPoolExecutor() as executor:
        futures = []
        for grid_size in np.arange(min_grid_size, middle_grid_size + step_size_small, step_size_small):
            futures.append(executor.submit(process_grid_size, grid_size, step_size_small))
        for grid_size in np.arange(middle_grid_size + step_size_large, max_grid_size + step_size_large, step_size_large):
            futures.append(executor.submit(process_grid_size, grid_size, step_size_large))

        for future in as_completed(futures):
            try:
                result = future.result()
                results.append(result)
                grid_size, occupancy_rate, _, _ = result
                print(f"Grid size: {grid_size} processed. Occupancy Rate: {occupancy_rate}")
            except Exception as e:
                print(f"Error processing grid size: {e}")

    # Sort results by grid size
    results.sort(key=lambda x: x[0])

    # Extract and return sorted results
    occupancy_rates = [(grid_size, occupancy_rate) for grid_size, occupancy_rate, _, _ in results]
    num_occupied_count = [num_occupied for _, _, num_occupied, _ in results]
    num_grids_count = [num_grids for _, _, _, num_grids in results]

    return occupancy_rates, num_occupied_count, num_grids_count


def save_occupancy_rates_to_txt(grid_sizes, rates, num_occupied_count, num_grids_count, output_file):
    """
    Save occupancy rate data to a TXT file.
    :param grid_sizes: List of grid sizes.
    :param rates: List of occupancy rates.
    :param num_occupied_count: List of occupied grid counts.
    :param num_grids_count: List of total grid counts.
    :param output_file: Path to the output TXT file.
    """
    with open(output_file, 'w') as f:
        f.write("Square_size\tOccupancy_rates\tOccupied_count\tSquare_counts\n")
        f.write("---------------------------------------\n")
        for size, rate, occupied_count, grids_count in zip(grid_sizes, rates, num_occupied_count, num_grids_count):
            f.write(f"{size}\t{rate}\t{occupied_count}\t{grids_count}\n")


def save_FD_to_txt(overall_fd, small_fd, large_fd, output_file):
    """
    Save fractal dimension (FD) data to a TXT file.
    :param overall_fd: Overall FD value.
    :param small_fd: Small-scale FD value.
    :param large_fd: Large-scale FD value.
    :param output_file: Path to the output TXT file.
    """
    with open(output_file, 'w') as f:
        f.write("Overall_FD\tSmall_scale_FD\tLarge_scale_FD\tFD_difference\n")
        f.write("---------------------------------------------------------------------------------------------------------------------\n")
        f.write(f"{overall_fd}\t{small_fd}\t{large_fd}\t{large_fd-small_fd}\n")


def load_occupancy_rates_from_txt(input_file):
    """
    Load occupancy rate data from a TXT file.
    :param input_file: Path to the input TXT file.
    :return: Lists of grid sizes, rates, occupied counts, and total grid counts.
    """
    grid_sizes = []
    rates = []
    num_occupied = []
    num_grids = []
    with open(input_file, 'r') as f:
        next(f)
        next(f)
        for line in f:
            size, rate, o_count, g_count = map(float, line.strip().split())
            grid_sizes.append(size)
            rates.append(rate)
            num_occupied.append(o_count)
            num_grids.append(g_count)
    return grid_sizes, rates, num_occupied, num_grids


def calculate_fd(grid_sizes, num_occupied, small_scales, large_scales, pixel_to_um):
    """
    Calculate fractal dimensions (FD).
    :param grid_sizes: Grid sizes in pixels.
    :param num_occupied: Counts of occupied grids.
    :param small_scales: Range of small scales (in microns).
    :param large_scales: Range of large scales (in microns).
    :param pixel_to_um: Conversion factor from pixels to microns.
    :return: FD values for overall, small scales, and large scales.
    """
    A = 1  # Normalization constant
    grid_sizes_um = np.array([i * pixel_to_um for i in grid_sizes])
    log_a_l = np.log(A / grid_sizes_um)
    log_counts = np.log(np.array(num_occupied))

    # Overall regression
    overall_slope, overall_intercept = np.polyfit(log_a_l, log_counts, 1)
    overall_fd = overall_slope

    # Small-scale regression
    small_mask = (grid_sizes_um >= small_scales[0]) & (grid_sizes_um <= small_scales[1])
    small_log_a_l = log_a_l[small_mask]
    small_log_counts = log_counts[small_mask]
    small_slope, small_intercept = np.polyfit(small_log_a_l, small_log_counts, 1)
    small_fd = small_slope

    # Large-scale regression
    large_mask = (grid_sizes_um >= large_scales[0]) & (grid_sizes_um <= large_scales[1])
    large_log_a_l = log_a_l[large_mask]
    large_log_counts = log_counts[large_mask]
    large_slope, large_intercept = np.polyfit(large_log_a_l, large_log_counts, 1)
    large_fd = large_slope

    return overall_fd, small_fd, large_fd, log_a_l, log_counts, overall_slope, overall_intercept, small_slope, small_intercept, large_slope, large_intercept, grid_sizes_um


def plot_occupancy(grid_sizes, rates, pixel_to_um):
    """
    Plot occupancy rate vs. grid size.
    :param grid_sizes: List of grid sizes in pixels.
    :param rates: List of occupancy rates.
    :param pixel_to_um: Conversion factor from pixels to microns.
    """
    plt.figure(figsize=(10, 10))
    um_sizes = [size * pixel_to_um for size in grid_sizes]
    plt.plot(um_sizes, rates, marker='o')
    plt.xlabel('Grid Size / microns')
    plt.ylabel('Occupancy Rate')
    plt.title('Occupancy Rate vs. Grid Size')
    plt.show()


def plot_fd(log_a_l, log_counts, overall_slope, overall_intercept, small_slope, small_intercept, large_slope,
            large_intercept, small_scales, large_scales, grid_sizes_um, errors=None):
    """
    Plot fractal dimension (FD) regression lines and data points.
    :param log_a_l: Logarithm of normalized grid sizes.
    :param log_counts: Logarithm of occupied grid counts.
    :param overall_slope: Overall regression slope.
    :param overall_intercept: Overall regression intercept.
    :param small_slope: Small-scale regression slope.
    :param small_intercept: Small-scale regression intercept.
    :param large_slope: Large-scale regression slope.
    :param large_intercept: Large-scale regression intercept.
    :param small_scales: Small-scale range (in microns).
    :param large_scales: Large-scale range (in microns).
    :param grid_sizes_um: Grid sizes in microns.
    :param errors: Optional error values for data points.
    """
    plt.figure(figsize=(20, 12))
    if errors is not None:
        plt.errorbar(log_a_l, log_counts, yerr=errors, fmt='o', color='blue', capsize=5, label='Data Points')
    else:
        plt.scatter(log_a_l, log_counts, color='blue', marker='o', label='Data Points')

    plt.plot(log_a_l, overall_slope * log_a_l + overall_intercept, color='red',
             label=f'Overall Fit Line (D={overall_slope:.2f})')

    small_mask = (grid_sizes_um >= small_scales[0]) & (grid_sizes_um <= small_scales[1])
    small_log_a_l = log_a_l[small_mask]
    plt.plot(small_log_a_l, small_slope * small_log_a_l + small_intercept, color='green',
             label=f'Small Scale Fit Line (D={small_slope:.2f})')

    large_mask = (grid_sizes_um >= large_scales[0]) & (grid_sizes_um <= large_scales[1])
    large_log_a_l = log_a_l[large_mask]
    plt.plot(large_log_a_l, large_slope * large_log_a_l + large_intercept, color='purple',
             label=f'Large Scale Fit Line (D={large_slope:.2f})')

    plt.title('Fractal Dimension (FD)')
    plt.xlabel('ln(1 / Square Size (µm))')
    plt.ylabel('ln(Square Counts)')
    plt.legend()
    plt.grid(True)

    grid_sizes_ticks = np.array([10, 20, 50, 100, 200, 400, 600])
    log_grid_sizes_ticks = np.log(1 / grid_sizes_ticks)
    plt.xticks(log_grid_sizes_ticks, grid_sizes_ticks)

    plt.ylim(0, 15)
    plt.yticks(np.arange(0, 16, 1))
    plt.show()


image_path = img

file_name = file_name

txt_file = '../CD20/' + file_name + '/points/global_coordinates.txt' #点坐标

min_um_size = 10
middle_um_size = 40
max_um_size = 600
step_size_small = 1
step_size_large = 20
small_scales = (10, 40)  
large_scales = (200, 600)  
pixel_to_um = 0.205 
#errors = [0.1] * len(num_occupied)

area_name = 'cancer_island' #： invasive_front |    cancer_center   |  cancer_island

if area_name == 'invasive_front':
    area_coords = invasive_front
    occupancy_rates_filepath = '../CD20/' + file_name + '/output_file/' + area_name + '_occupancy_rates.txt'
    FD_txt_filepath = '../CD20/' + file_name + '/output_file/' + area_name + '_FD.txt'
    
    
elif area_name == 'cancer_center':
    area_coords = cancer_center
    occupancy_rates_filepath = '../CD20/' + file_name + '/output_file/' + area_name + '_occupancy_rates.txt'
    FD_txt_filepath = '../CD20/' + file_name + '/output_file/' + area_name + '_FD.txt'
    
elif area_name == 'cancer_island':
    area_coords = cancer_island
    occupancy_rates_filepath = '../CD20/' + file_name + '/output_file/' + area_name + '_occupancy_rates.txt'
    FD_txt_filepath = '../CD20/' + file_name + '/output_file/' + area_name + '_FD.txt'
    
else:
    print('error, check the name')

occupancy_rates,num_occupied_count, num_grids_count = calculate_area_occupancy_rate(image_path, txt_file, min_um_size, middle_um_size, max_um_size, step_size_small, step_size_large, area_coords, pixel_to_um)

grid_sizes, rates = zip(*occupancy_rates)

save_occupancy_rates_to_txt(grid_sizes, rates, num_occupied_count, num_grids_count, occupancy_rates_filepath)

grid_sizes, rates, num_occupied, num_grids = load_occupancy_rates_from_txt(occupancy_rates_filepath)


for grid_size, rate in occupancy_rates:
    print(f"Grid Size: {grid_size}, Occupancy Rate: {rate:.4f}")
#grid_sizes, rates = zip(*occupancy_rates)
overall_fd, small_fd, large_fd, log_a_l, log_counts, overall_slope, overall_intercept, small_slope, small_intercept, large_slope, large_intercept, grid_sizes_um = calculate_fd(grid_sizes, num_occupied, small_scales, large_scales, pixel_to_um)
save_FD_to_txt(overall_fd, small_fd, large_fd, FD_txt_filepath)

plot_occupancy(grid_sizes, rates, pixel_to_um)
plot_fd(log_a_l, log_counts, overall_slope, overall_intercept, small_slope, small_intercept, large_slope, large_intercept, small_scales, large_scales, grid_sizes_um)