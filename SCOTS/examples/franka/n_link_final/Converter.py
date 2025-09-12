import csv
import math
import os

def parse_scs_grid(lines, start_index):
    """
    Parses the grid properties (DIM, ETA, LOWER_LEFT, UPPER_RIGHT) from the SCS file.

    Args:
        lines (list): A list of lines from the SCS file.
        start_index (int): The starting index to look for grid properties.

    Returns:
        tuple: A tuple containing (dim, eta_list, lower_left, upper_right)
               or (None, None, None, None) if parsing fails.
    """
    dim, eta_list, lower_left, upper_right = None, [], [], []
    current_index = start_index

    # Find and parse DIM
    while current_index < len(lines) and "#MEMBER:DIM" not in lines[current_index]:
        current_index += 1
    if current_index >= len(lines):
        return None, None, None, None
    dim = int(lines[current_index + 1].strip())
    current_index += 2

    # Find and parse ETA vector
    while current_index < len(lines) and "#VECTOR:ETA" not in lines[current_index]:
        current_index += 1
    if current_index >= len(lines):
        return None, None, None, None
    current_index += 2 # Skip #BEGIN:n
    for _ in range(dim):
        eta_list.append(float(lines[current_index].strip()))
        current_index += 1
    current_index += 1 # Skip #END

    # Find and parse LOWER_LEFT vector
    while current_index < len(lines) and "#VECTOR:LOWER_LEFT" not in lines[current_index]:
        current_index += 1
    if current_index >= len(lines):
        return None, None, None, None
    current_index += 2 # Skip #BEGIN:n
    for _ in range(dim):
        lower_left.append(float(lines[current_index].strip()))
        current_index += 1
    current_index += 1 # Skip #END

    # Find and parse UPPER_RIGHT vector
    while current_index < len(lines) and "#VECTOR:UPPER_RIGHT" not in lines[current_index]:
        current_index += 1
    if current_index >= len(lines):
        return None, None, None, None
    current_index += 2 # Skip #BEGIN:n
    for _ in range(dim):
        upper_right.append(float(lines[current_index].strip()))
        current_index += 1
    current_index += 1 # Skip #END

    return dim, eta_list, lower_left, upper_right

def grid_index_to_coords(index, dim, eta_list, lower_left, upper_right):
    """
    Converts a single grid index to n-dimensional coordinates.

    Args:
        index (int): The grid point index.
        dim (int): The dimension of the space.
        eta_list (list): A list of grid cell sizes for each dimension.
        lower_left (list): A list of lower-left corner coordinates.
        upper_right (list): A list of upper-right corner coordinates.

    Returns:
        list: The n-dimensional coordinates.
    """
    num_points_per_dim = [
        round((upper_right[i] - lower_left[i]) / eta_list[i]) + 1
        for i in range(dim)
    ]
    
    multi_index = [0] * dim
    
    temp_index = index
    for i in range(dim):
        multi_index[i] = temp_index % num_points_per_dim[i]
        temp_index = temp_index // num_points_per_dim[i]
    
    # Corrected formula: The previous version added eta_list[i] / 2, which
    # calculated the center of the grid cell. The correct method for this
    # dataset is to use the lower-left corner of the grid cell.
    coords = [
        lower_left[i] + multi_index[i] * eta_list[i]
        for i in range(dim)
    ]
    
    return coords


def convert_target_scs(input_file, output_file):
    """
    Converts the target.scs file to a CSV with n-dimensional coordinates.
    """
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
        return

    # Parse grid properties
    dim, eta_list, lower_left, upper_right = parse_scs_grid(lines, 0)
    if not dim:
        print(f"Error: Failed to parse grid properties from '{input_file}'.")
        return

    # Find and parse grid points
    gridpoints = []
    try:
        start_points_index = lines.index("#VECTOR:GRIDPOINTS\n") + 2
        for line in lines[start_points_index:]:
            if line.strip() == "#END":
                break
            gridpoints.append(int(line.strip()))
    except ValueError:
        print(f"Error: '#VECTOR:GRIDPOINTS' not found in '{input_file}'.")
        return

    # Convert grid points to coordinates and write to CSV
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        header = [f"theta_{i+1}" for i in range(dim)]
        writer.writerow(header)
        
        for point_index in gridpoints:
            coords = grid_index_to_coords(point_index, dim, eta_list, lower_left, upper_right)
            writer.writerow(coords)

    print(f"Successfully converted '{input_file}' to '{output_file}'.")


def convert_controller_scs(input_file, output_file):
    """
    Converts the controller.scs file to a CSV mapping state coordinates to input coordinates.
    """
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
        return

    # Parse state space grid
    state_grid_start = lines.index("#SCOTS:STATE_SPACE\n")
    state_dim, state_eta, state_lower, state_upper = parse_scs_grid(lines, state_grid_start)
    if not state_dim:
        print(f"Error: Failed to parse state space grid from '{input_file}'.")
        return

    # Parse input space grid
    input_grid_start = lines.index("#SCOTS:INPUT_SPACE\n")
    input_dim, input_eta, input_lower, input_upper = parse_scs_grid(lines, input_grid_start)
    if not input_dim:
        print(f"Error: Failed to parse input space grid from '{input_file}'.")
        return

    # Find and parse the data matrix
    data_start = lines.index("#MATRIX:DATA\n") + 2 # Skip the line with #BEGIN
    
    # Write to CSV
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        header = [f"state_theta_{i+1}" for i in range(state_dim)] + [f"input_u_{i+1}" for i in range(input_dim)]
        writer.writerow(header)

        for line in lines[data_start:]:
            if line.strip() == "#END":
                break
            parts = line.strip().split()
            state_index = int(parts[0])
            input_indices = [int(p) for p in parts[1:]]

            state_coords = grid_index_to_coords(state_index, state_dim, state_eta, state_lower, state_upper)
            
            for input_index in input_indices:
                input_coords = grid_index_to_coords(input_index, input_dim, input_eta, input_lower, input_upper)
                row = state_coords + input_coords
                writer.writerow(row)

    print(f"Successfully converted '{input_file}' to '{output_file}'.")

if __name__ == "__main__":
    # Define file names
    target_scs_file = '/home/focaslab/Downloads/SCOTS_ros2_v2/SCOTS/examples/franka/n_link_final/n=7/target.scs'
    target_csv_file = '/home/focaslab/Downloads/SCOTS_ros2_v2/SCOTS/examples/franka/n_link_final/n=7/target.csv'
    controller_scs_file = '/home/focaslab/Downloads/SCOTS_ros2_v2/SCOTS/examples/franka/n_link_final/n=7/controller.scs'
    controller_csv_file = '/home/focaslab/Downloads/SCOTS_ros2_v2/SCOTS/examples/franka/n_link_final/n=7/controller.csv'

    # Convert the files
    convert_target_scs(target_scs_file, target_csv_file)
    convert_controller_scs(controller_scs_file, controller_csv_file)
