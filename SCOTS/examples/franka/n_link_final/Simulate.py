import csv
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.spatial import KDTree

def load_data(controller_file, target_file):
    """
    Loads data from controller.csv and target.csv files.
    
    Determines the number of links (n) from the controller data.
    
    Args:
        controller_file (str): Path to the controller CSV file.
        target_file (str): Path to the target CSV file.
    
    Returns:
        tuple: A tuple containing (n, controller_data, target_states_tree)
               where n is the number of links, controller_data is a dictionary
               mapping states to inputs, and target_states_tree is a KDTree
               of target states for efficient lookup.
    """
    try:
        # Load controller data
        with open(controller_file, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            
            # Determine n from the header
            num_cols = len(header)
            if num_cols % 2 != 0:
                print("Error: The number of columns in the controller file is not even.")
                return None, None, None
            n = num_cols // 2
            
            controller_data = {}
            for row in reader:
                state_coords = tuple(float(x) for x in row[:n])
                input_coords = tuple(float(x) for x in row[n:])
                controller_data[state_coords] = input_coords
                
        # Load target data and create KDTree for fast lookup
        with open(target_file, 'r') as f:
            reader = csv.reader(f)
            next(reader) # Skip header
            target_states = [tuple(float(x) for x in row) for row in reader]
            target_states_tree = KDTree(target_states)

    except FileNotFoundError as e:
        print(f"Error: {e}. Please ensure the CSV files are in the same directory.")
        return None, None, None
    except ValueError:
        print("Error: Invalid data format in the CSV files. Make sure all values are numbers.")
        return None, None, None

    return n, controller_data, target_states_tree


def forward_kinematics(theta, link_lengths):
    """
    Calculates the (x, y) coordinates of the robotic arm joints.

    Args:
        theta (list): A list of joint angles in radians.
        link_lengths (list): A list of link lengths.

    Returns:
        tuple: A tuple of lists (x_coords, y_coords) for each joint.
    """
    x_coords = [0.0]
    y_coords = [0.0]
    current_angle = 0.0
    for i in range(len(theta)):
        current_angle += theta[i]
        x_coords.append(x_coords[-1] + link_lengths[i] * math.cos(current_angle))
        y_coords.append(y_coords[-1] + link_lengths[i] * math.sin(current_angle))
    return x_coords, y_coords


def animate_robot(theta_history, link_lengths, target_states):
    """
    Animates the robotic arm's movement using matplotlib.
    """
    n = len(link_lengths)
    fig, ax = plt.subplots()
    
    # Set plot limits
    max_reach = sum(link_lengths)
    ax.set_xlim(-max_reach * 1.1, max_reach * 1.1)
    ax.set_ylim(-max_reach * 1.1, max_reach * 1.1)
    ax.set_aspect('equal', 'box')
    ax.set_title("N-Link Robotic Arm Simulation")
    
    # Plot the target set
    # Corrected conditional check to avoid ValueError
    if target_states.size > 0:
        # For visualization, just plot the end effector positions of target states
        target_coords = [forward_kinematics(state, link_lengths)[0:2] for state in target_states]
        end_effector_x = [coords[0][-1] for coords in target_coords]
        end_effector_y = [coords[1][-1] for coords in target_coords]
        ax.plot(end_effector_x, end_effector_y, 'ro', markersize=2, label="Target Set")
        ax.legend()
        
    line, = ax.plot([], [], 'o-', lw=2, markersize=8)
    
    def init():
        line.set_data([], [])
        return line,

    def update(frame):
        theta = theta_history[frame]
        x_coords, y_coords = forward_kinematics(theta, link_lengths)
        line.set_data(x_coords, y_coords)
        return line,

    anim = animation.FuncAnimation(fig, update, frames=len(theta_history), 
                                   init_func=init, blit=True, interval=50)
    plt.show()
    print("Simulation complete. The plot window has been displayed.")


def main():
    """
    Main function to run the simulation and plot.
    """
    # Define file paths
    controller_file = '/home/focaslab/Downloads/SCOTS_ros2_v2/SCOTS/examples/franka/n_link_final/n=7/controller.csv'
    target_file = '/home/focaslab/Downloads/SCOTS_ros2_v2/SCOTS/examples/franka/n_link_final/n=7/target.csv'
    
    # Load data from CSVs
    n, controller_data, target_tree = load_data(controller_file, target_file)
    if n is None:
        return

    # Extract states for KDTree and map to inputs
    controller_states = list(controller_data.keys())
    controller_inputs = list(controller_data.values())
    controller_tree = KDTree(controller_states)
    
    # Simulation parameters
    dt = 0.05  # Time step
    max_sim_steps = 1000 # Max steps to avoid infinite loops
    
    # Initial state of the arm (all joint angles at 0)
    current_theta = np.zeros(n)
    
    theta_history = [current_theta.copy()]
    
    print(f"Starting simulation for an {n}-link robotic arm...")

    # Simulation loop
    for step in range(max_sim_steps):
        # Check if the current state is in the target set
        # We find the closest target state and check its distance
        dist, _ = target_tree.query(current_theta, k=1)
        if dist < 1e-6: # A small tolerance for floating point comparison
            print(f"Target reached at step {step}!")
            break

        # Find the closest state in the controller's winning domain
        _, idx = controller_tree.query(current_theta, k=1)
        
        # Get the corresponding input
        control_input = controller_inputs[idx]
        
        # Update the state using single integrator dynamics (Euler integration)
        current_theta += np.array(control_input) * dt
        
        # Store the state for animation
        theta_history.append(current_theta.copy())

    if step == max_sim_steps - 1:
        print("Simulation stopped after reaching max simulation steps.")
    
    # Assume all link lengths are 1 for plotting
    link_lengths = [1.0] * n
    
    # Plot the simulation
    animate_robot(theta_history, link_lengths, target_tree.data)


if __name__ == "__main__":
    main()
