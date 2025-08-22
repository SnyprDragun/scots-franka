import rclpy
from rclpy.node import Node
from geometry_msgs.msg import Twist
from turtlesim.msg import Pose
import numpy as np
import pandas as pd
import time
from functools import reduce
from math import atan2

class TM_Pub_Sub(Node):
    def __init__(self):
        super().__init__('tm_pub_sub')
        
        # Initializing variables
        self.tm_vel = Twist()
        self.count = 0
        self.tm_vel.linear.x = 0.0  # Explicitly setting it as a float
        self.tm_vel.linear.y = 0.0
        self.tm_vel.linear.z = 0.0
        self.tm_vel.angular.x = 0.0
        self.tm_vel.angular.y = 0.0
        self.tm_vel.angular.z = 0.0
        self.data_pos = []
        self.ini_x = 0.5
        self.ini_y = 0.5
        self.ini_theta = 0
        self.data_pos.append([self.ini_x, self.ini_y, self.ini_theta])
        
        # Publisher and Subscriber
        self.pub = self.create_publisher(Twist, '/turtle1/cmd_vel', 10)
        self.pose_subscriber = self.create_subscription(Pose, '/turtle1/pose', self.tm_subscribing_function, 10)
        
        # Data initialization
        self.data = pd.read_csv("scots_lab1.csv")
        self.data_in_np = np.array(self.data)
        self.wheel_radius = 0.038
        self.length = 0.23
        self.states = self.data_in_np[:, 0:3]
        self.input = self.data_in_np[:, 3:len(self.data_in_np)]
        self.x_quantization = 0.05
        self.y_quantization = 0.05
        self.theta_quantization = 0.2
        self.count = 0
        self.ini_time = time.time()

    def tm_subscribing_function(self, msg):
        self.tm_pos_x = msg.x
        self.tm_pos_y = msg.y
        self.tm_quart = msg.theta

    def tm_controller(self):
        scots_input = [0, 0]
        while rclpy.ok():
            rclpy.spin_once(self)  # Process incoming messages

            # If within certain regions, stop the robot
            if ((self.tm_pos_x <= 2.6 and self.tm_pos_y <= 2.6 and self.tm_pos_y >= 2.0 and self.tm_pos_x >= 2.0) or
                (self.tm_pos_x >= 1.3 and self.tm_pos_y <= 1.6 and self.tm_pos_y >= 1.3 and self.tm_pos_x <= 1.6)):
                self.tm_vel.linear.x = 0.0
                self.tm_vel.angular.z = 0.0
                self.pub.publish(self.tm_vel)
                if self.tm_pos_x <= 2.6 and self.tm_pos_y <= 2.6 and self.tm_pos_y >= 2 and self.tm_pos_x >= 2:
                    self.get_logger().info("Reached target")
                else:
                    self.get_logger().info("Collision detected")
                break
            else:
                # Find closest matching state
                stat1 = np.where(np.abs(self.tm_pos_x - self.states[:, 0]) <= self.x_quantization / 2)
                stat2 = np.where(np.abs(self.tm_pos_y - self.states[:, 1]) <= self.y_quantization / 2)
                stat3 = np.where(np.abs(self.tm_quart - self.states[:, 2]) <= self.theta_quantization / 2)
                com_state = reduce(np.intersect1d, (stat1, stat2, stat3))
                
                if com_state.size != 0:
                    scots_input = self.input[com_state[0], :]
                else:
                    scots_input = scots_input
                
                # Set velocities
                self.tm_vel.linear.x = scots_input[0]
                self.tm_vel.angular.z = scots_input[1]
                self.pub.publish(self.tm_vel)
                
                # Normalize yaw value
                self.yaw_val = ((2 * np.pi + self.tm_quart) * (self.tm_quart < 0) + self.tm_quart * (self.tm_quart > 0))
                
                # Append data
                self.data = [self.tm_pos_x, self.tm_pos_y, self.yaw_val, scots_input[0], scots_input[1]]
                self.data_pos.append(self.data)
                self.get_logger().info(f"Data: {self.data}")

def main(args=None):
    rclpy.init(args=args)
    tm_obj = TM_Pub_Sub()
    tm_obj.tm_controller()
    rclpy.shutdown()

if __name__ == "__main__":
    main()

