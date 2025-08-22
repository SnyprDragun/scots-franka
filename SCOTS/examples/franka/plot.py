# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation


# # =====================================================
# # Part 1: SCS to CSV Converter
# # =====================================================

# class scs_csv_converter():
#     def __init__(self):
#         """Class for converting scs file to csv"""
#         pass
    
#     def read_csv(self, file_name, file_type):
#         """Reads the scs file header and extracts metadata"""
#         self.s_dim = None
#         self.u_dim = None
        
#         with open(file_name, "r") as scs_file:
#             self.text = scs_file.readlines()
    
#         if file_type == "scots": 
#             line = 5
#         else: 
#             line = 3

#         self.s_dim = int(self.text[line][:-1])   # number of states
#         line += 3
        
#         self.s_eta = self.extract_values(line, self.s_dim)
#         line += self.s_dim + 3
        
#         self.s_lb = self.extract_values(line, self.s_dim)
#         line += self.s_dim + 3
        
#         self.s_ub = self.extract_values(line, self.s_dim)
#         line += self.s_dim + 4

#         if file_type == "scots":
#             self.u_dim = int(self.text[line][:-1])  # number of inputs
#             line += 3

#             self.u_eta = self.extract_values(line, self.u_dim)
#             line += self.u_dim + 3
            
#             self.u_lb = self.extract_values(line, self.u_dim)
#             line += self.u_dim + 3
            
#             self.u_ub = self.extract_values(line, self.u_dim)
#             line += self.u_dim + 5
            
#             self.u_shape = (((self.u_ub - self.u_lb) / self.u_eta) + 1).astype(int)
        
#         self.line = line
#         self.s_shape = (((self.s_ub - self.s_lb) / self.s_eta) + 1).astype(int)

    
#     def extract_values(self, start_index, num_indices): 
#         """Extract values from text into numpy arrays"""
#         return np.array([float(self.text[start_index+i][:-1]) for i in range(num_indices)], dtype=float)
    
    
#     def index_to_state(self, index):
#         """Convert state index to actual state values"""
#         coords = np.zeros(self.s_dim)
#         for dim in range(self.s_dim):
#             coords[dim] = index % self.s_shape[dim]
#             index //= self.s_shape[dim]
#         return self.s_lb + coords * self.s_eta
    
    
#     def index_to_input(self, index):
#         """Convert input index to actual input values"""
#         coords = np.zeros(self.u_dim)
#         for dim in range(self.u_dim):
#             coords[dim] = index % self.u_shape[dim]
#             index //= self.u_shape[dim]
#         return self.u_lb + coords * self.u_eta

    
#     def get_state_input_idx(self, curr_line):
#         """Extract state and input indices from one line of .scs"""
#         line_list = curr_line[:-1].split(' ')
#         state_idx = int(line_list[0])
#         if self.u_dim is not None: 
#             input_idx = np.max(np.array([int(i) for i in line_list[1:]], dtype=int))
#         else: 
#             input_idx = None
#         return state_idx, input_idx
            
            
#     def write_csv(self, filename):
#         """Write extracted state/input values to CSV"""
#         csv_file = open(filename, "w+")
#         text = ','.join([f"state{i}" for i in range(self.s_dim)])
#         if self.u_dim is not None: 
#             text += ',' + ','.join([f"input{i}" for i in range(self.u_dim)])
#         csv_file.write(text + '\n')

#         for curr_line in self.text[self.line:-1]:
#             state_idx, input_idx = self.get_state_input_idx(curr_line)
#             state = np.round(self.index_to_state(state_idx), 5)
#             text = ','.join([str(i) for i in state])

#             if self.u_dim is not None:
#                 ctrl_input = np.round(self.index_to_input(input_idx), 5)
#                 text += ',' + ','.join([str(i) for i in ctrl_input])

#             csv_file.write(text + '\n')
#         csv_file.close()


# # =====================================================
# # Part 2: Arm Simulation (Single-Integrator Dynamics)
# # =====================================================

# def arm_controller(states, states_data, inputs_data, states_quant):
#     """
#     Match current state to nearest quantized state from SCOTS and return input.
#     States are just joint angles (no velocities).
#     """
#     diffs = np.abs(states_data - states)
#     mask = np.all(diffs <= (states_quant / 2), axis=1)
#     idx = np.where(mask)[0]
#     if len(idx) > 0:
#         return inputs_data[idx[0], :]
#     else:
#         return np.zeros(states.shape[0])

# def update_states(states, u, dt):
#     """
#     Single-integrator dynamics:
#         q_dot = u
#         q_next = q + u*dt
#     """
#     q = states + u * dt
#     # wrap to [-pi, pi]
#     q = (q + np.pi) % (2 * np.pi) - np.pi
#     return q

# def forward_kinematics(q, L):
#     """Compute joint positions in 2D from joint angles q and link lengths L"""
#     n = len(q)
#     joint_positions = np.zeros((n+1, 2))  # include base at (0,0)
#     total_angle = 0
#     pos = np.array([0.0, 0.0])
#     joint_positions[0, :] = pos
#     for i in range(n):
#         total_angle += q[i]
#         pos = pos + L[i] * np.array([np.cos(total_angle), np.sin(total_angle)])
#         joint_positions[i+1, :] = pos
#     return joint_positions


# def simulate(n=3, dt=0.01, duration=1000, csv_filename="franka_nlink_states.csv"):
#     states = np.loadtxt(csv_filename, delimiter=",", skiprows=1)

#     if states.ndim == 1:
#         states = states.reshape(1, -1)  # ensure 2D

#     fig, ax = plt.subplots()
#     ax.set_xlim(-n, n)
#     ax.set_ylim(-n, n)
#     ax.set_aspect("equal")

#     # Line for arm
#     arm_line, = ax.plot([], [], "o-", lw=2)

#     # One point per joint (including base)
#     joint_points = [ax.plot([], [], "ro")[0] for _ in range(n+1)]

#     def init():
#         arm_line.set_data([], [])
#         for jp in joint_points:
#             jp.set_data([], [])
#         return [arm_line] + joint_points

#     def animate(i):
#         angles = states[i, :n]   # first n entries = joint angles
#         x, y = [0], [0]
#         theta_sum = 0
#         for theta in angles:
#             theta_sum += theta
#             x.append(x[-1] + np.cos(theta_sum))
#             y.append(y[-1] + np.sin(theta_sum))

#         arm_line.set_data(x, y)
#         for j, jp in enumerate(joint_points):
#             jp.set_data([x[j]], [y[j]])
#         return [arm_line] + joint_points

#     ani = animation.FuncAnimation(
#         fig, animate, frames=min(duration, len(states)),
#         init_func=init, blit=True, interval=dt*1000
#     )
#     plt.show()


# # =====================================================
# # Run everything
# # =====================================================
# if __name__ == "__main__":
#     import os

#     scs_file = os.path.expanduser("~/Downloads/SCOTS_ros2_v2/SCOTS/examples/franka/franka_nlink_single_integrator.scs")
#     controller_csv = "franka_nlink_controller.csv"

#     # Step 1: Convert .scs -> controller CSV
#     conv = scs_csv_converter()
#     conv.read_csv(scs_file, "scots")
#     conv.write_csv(controller_csv)

#     # Step 2: Run simulation using SCOTS controller
#     simulate(n=3, dt=0.01, duration=1000, csv_filename=controller_csv)




import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ------------------------------
# Convert SCS output to CSV
# ------------------------------
class SCS2CSV:
    def __init__(self, scs_filename, csv_filename, n):
        self.scs_filename = scs_filename
        self.csv_filename = csv_filename
        self.n = n

    def convert(self):
        states = []
        with open(self.scs_filename, "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split(",")
                if len(parts) < self.n:
                    continue
                try:
                    state = [float(x) for x in parts[:self.n]]
                    states.append(state)
                except ValueError:
                    continue
        states = np.array(states)
        np.savetxt(self.csv_filename, states, delimiter=",", header=",".join([f"theta{i}" for i in range(self.n)]), comments="")
        print(f"Converted {len(states)} states to {self.csv_filename}")


# ------------------------------
# Forward kinematics
# ------------------------------
def forward_kinematics(angles, link_length=1.0):
    """
    Compute the (x,y) positions of each joint and end-effector.
    angles: list of joint angles [theta1,...,theta_n]
    returns: array of joint positions shape (n+1, 2)
    """
    n = len(angles)
    joints = np.zeros((n + 1, 2))  # joint0 at origin
    x, y, theta = 0.0, 0.0, 0.0
    for i in range(n):
        theta += angles[i]
        x += link_length * np.cos(theta)
        y += link_length * np.sin(theta)
        joints[i + 1] = [x, y]
    return joints


# ------------------------------
# Animate function
# ------------------------------
def simulate(n=2, dt=0.01, duration=1000, csv_filename="/home/focaslab/Downloads/SCOTS_ros2_v2/SCOTS/examples/franka/controller.csv"):
    states = np.loadtxt(csv_filename, delimiter=",", skiprows=1)  # skip header

    # Ensure states is 2D (T x n)
    if states.ndim == 1:
        states = states.reshape(1, -1)

    fig, ax = plt.subplots()
    ax.set_xlim(-n - 0.5, n + 0.5)
    ax.set_ylim(-n - 0.5, n + 0.5)
    ax.set_aspect("equal")

    # Lines for arm and scatter for joints
    line, = ax.plot([], [], "o-", lw=2)
    joint_points = [ax.plot([], [], "ro")[0] for _ in range(n + 1)]

    def init():
        line.set_data([], [])
        for jp in joint_points:
            jp.set_data([], [])
        return [line] + joint_points

    def animate(i):
        angles = states[i, :n]  # only first n entries are joint angles
        joints = forward_kinematics(angles)

        line.set_data(joints[:, 0], joints[:, 1])
        for j in range(n + 1):
            joint_points[j].set_data([joints[j, 0]], [joints[j, 1]])
        return [line] + joint_points

    ani = animation.FuncAnimation(fig, animate, frames=len(states),
                                  init_func=init, blit=True, interval=dt * 1000, repeat=False)

    plt.show()


# ------------------------------
# Main
# ------------------------------
if __name__ == "__main__":
    # Example usage: Convert and plot for 3-link arm
    converter = SCS2CSV("/home/focaslab/Downloads/SCOTS_ros2_v2/SCOTS/examples/franka/franka_nlink_single_integrator.scs", "/home/focaslab/Downloads/SCOTS_ros2_v2/SCOTS/examples/franka/controller.csv", n=3)
    converter.convert()

    simulate(n=3, dt=0.05, duration=1000, csv_filename="/home/focaslab/Downloads/SCOTS_ros2_v2/SCOTS/examples/franka/controller.csv")
