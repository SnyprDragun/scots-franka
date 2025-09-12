clc; 
clear; 
clf;

filename = 'controller_6.csv';

% =========================
% Read CSV (first row skip)
% =========================
data = readmatrix(filename);
data(1,:) = []; % skip header
[nSteps_csv, cols] = size(data);

n = cols/2; % number of joints
states_data = data(:, 1:n);      % states (angles)
input_data  = data(:, n+1:end);  % inputs (velocities)

% =========================
% Simulation parameters
% =========================
T  = 10.0;
dt = 1e-3;
num_steps = round(T/dt);
x = 0.0*(1:n);   % initial state [0.1,0.2,0.3,...]
xi_VCZ = zeros(num_steps, n);
u_VCZ  = zeros(num_steps, n);
if n == 7
    EE_VCZ = zeros(num_steps, 3);
    EE     = zeros(num_steps, 3);
else
    EE_VCZ = zeros(num_steps, 2);
    EE     = zeros(num_steps, 2);
end

time_arr = 0:dt:T-dt;

% =========================
% VCZ trajectory generation
% =========================
for step = 1:num_steps
    diffs = vecnorm(states_data - x, 2, 2);
    [~, closest_idx] = min(diffs);
    u = input_data(closest_idx, :);

    xi_VCZ(step,:) = x;
    u_VCZ(step,:)  = u;

    if n == 7
        % Franka EE
        T_fk = franka_fk(x);
        EE_VCZ(step,:) = T_fk(1:3,4)'; % EE position
    else
        % Generic n-link toy arm
        EE_VCZ(step,:) = VCZ_kinematics(x, ones(1,n));
    end

    x = VCZ_dynamics(x, u, dt);
end

% =========================
% VCZ control law
% =========================
th0    = 0.1*(1:n);
th     = zeros(num_steps, n);
om     = zeros(num_steps, n);
om_d   = zeros(num_steps, n);
torque = zeros(num_steps, n);

if n == 7
    EE_VCZ = zeros(num_steps, 3);
    EE     = zeros(num_steps, 3);
else
    EE_VCZ = zeros(num_steps, 2);
    EE     = zeros(num_steps, 2);
end

th(1,:) = th0;

if n == 7
    T_fk = franka_fk(th0);
    EE(1,:) = T_fk(1:3,4)';
else
    EE(1,:) = VCZ_kinematics(th0, ones(1,n));
end


lam1 = 0.02;
lam2 = 0.1;
taub = linspace(10, 1, n);   % base joint strongest, EE weakest
omb  = 1.0;

for step = 2:num_steps
    q  = [th(step-1,:)'; om(step-1,:)'];
    dq = RR_ode_n(q, torque(step-1,:), n);

    th(step,:) = th(step-1,:) + dt*dq(1:n)';
    om(step,:) = om(step-1,:) + dt*dq(n+1:end)';

    et = norm(th(step,:) - xi_VCZ(step,:))/lam1;
    if norm(th(step,:) - xi_VCZ(step,:)) > 1e-8
        om_d(step,:) = -omb*psi(et).*(th(step,:) - xi_VCZ(step,:)) / norm(th(step,:) - xi_VCZ(step,:));
    end

    eo = norm(om(step,:) - om_d(step,:))/lam2;
    if norm(om(step,:) - om_d(step,:)) > 1e-8
        torque(step,:) = -psi(eo) * ( taub .* (om(step,:) - om_d(step,:)) ) / norm(om(step,:) - om_d(step,:));
    end

    if n == 7
        T_fk = franka_fk(th(step,:));
        EE(step,:) = T_fk(1:3,4)';
    else
        EE(step,:) = VCZ_kinematics(th(step,:), ones(1,n));
    end
end



% =========================
% Plot results
% =========================
figure(1); clf;

for i = 1:n
    subplot(3,n,i) % angles
    hold on;
    plot(time_arr, xi_VCZ(:,i), 'k:', 'LineWidth', 2, ...
        'DisplayName', sprintf('$\\xi_%d(t)$ - VCZ', i));
    plot(time_arr, th(:,i), 'b', 'LineWidth', 2, ...
        'DisplayName', sprintf('$\\theta_%d(t)$ - VCZ control', i));
    xlabel('$$t(s)$$','Interpreter','latex');
    ylabel(sprintf('$$\\xi_%d,\\theta_%d(rad)$$', i, i),'Interpreter','latex');
    legend('Interpreter','latex'); grid on; box on;

    subplot(3,n,n+i) % velocities
    hold on;
    plot(time_arr, om(:,i), 'b', 'LineWidth', 2, ...
        'DisplayName', sprintf('$\\omega_%d(t)$', i));
    xlabel('$$t(s)$$','Interpreter','latex');
    ylabel(sprintf('$$\\omega_%d(rad/s)$$', i),'Interpreter','latex');
    legend('Interpreter','latex'); grid on; box on;

    subplot(3,n,2*n+i) % torques
    hold on;
    plot(time_arr, torque(:,i), 'b', 'LineWidth', 2, ...
        'DisplayName', sprintf('$\\tau_%d(t)$', i));
    xlabel('$$t(s)$$','Interpreter','latex');
    ylabel(sprintf('$$\\tau_%d(Nm)$$', i),'Interpreter','latex');
    legend('Interpreter','latex'); grid on; box on;
end

% =========================
% End-effector trajectory
% =========================
if n == 7
    figure(2); clf; hold on;
    plot3(EE_VCZ(:,1), EE_VCZ(:,2), EE_VCZ(:,3), 'k:', 'LineWidth', 2, 'DisplayName', 'VCZ reference EE');
    plot3(EE(:,1), EE(:,2), EE(:,3), 'b', 'LineWidth', 2, 'DisplayName', 'Controlled EE');
    scatter3(EE_VCZ(1,1), EE_VCZ(1,2), EE_VCZ(1,3), 60, 'r', 'filled', 'DisplayName','Start');
    scatter3(EE_VCZ(end,1), EE_VCZ(end,2), EE_VCZ(end,3), 60, 'g', 'filled', 'DisplayName','End');
    xlabel('$$x (m)$$','Interpreter','latex');
    ylabel('$$y (m)$$','Interpreter','latex');
    zlabel('$$z (m)$$','Interpreter','latex');
    title('End-Effector Trajectory (3D)','Interpreter','latex');
    legend('Interpreter','latex','Location','best');
    grid on; box on; axis equal;
    view([15,15])
else
    % keep your old 2D plot
    figure(2); clf; hold on;
    plot(EE_VCZ(:,1), EE_VCZ(:,2), 'k:', 'LineWidth', 2, 'DisplayName', 'VCZ reference EE');
    plot(EE(:,1), EE(:,2), 'b', 'LineWidth', 2, 'DisplayName', 'Controlled EE');
    scatter(EE_VCZ(1,1), EE_VCZ(1,2), 60, 'r', 'filled', 'DisplayName','Start');
    scatter(EE_VCZ(end,1), EE_VCZ(end,2), 60, 'g', 'filled', 'DisplayName','End');
    xlabel('$$x (m)$$','Interpreter','latex');
    ylabel('$$y (m)$$','Interpreter','latex');
    title('End-Effector Trajectory (2D)','Interpreter','latex');
    legend('Interpreter','latex','Location','best');
    grid on; box on;
end


disp('Simulation complete.');
disp('Script completed successfully!');

%% ========== FUNCTIONS ==========

function p = psi(e)
    a = 5;
    p = tanh(a*e).^3;
end

% === Forward dynamics (integrator form) ===
function x_next = VCZ_dynamics(x,u,dt)
    x_next = x + u*dt; % simple integrator model
end

% === Forward kinematics for n links ===
function pos = VCZ_kinematics(theta, lengths)
    n = length(theta);
    x = 0; y = 0;
    ang = 0;
    for i=1:n
        ang = ang + theta(i);
        x = x + lengths(i)*cos(ang);
        y = y + lengths(i)*sin(ang);
    end
    pos = [x, y];
end

% === n-link robot ODE ===
function dqdt = RR_ode_n(x,u,n)
    % Very simplified n-link decoupled dynamics
    % State: [theta(1..n), omega(1..n)]
    % u: torque (1..n)
    m = 0.1; l = 1; g = 0; % toy params
    th  = x(1:n);
    dth = x(n+1:end);
    ddth = (u' - 0.1*dth); % toy dynamics (damped integrator)
    dqdt = [dth; ddth];
end

function T = franka_fk(theta)
    % Forward kinematics for Franka Emika Panda
    % Input: theta (7x1 or 1x7 joint angles in radians)
    % Output: 4x4 homogeneous transform (EE pose w.r.t base)

    if size(theta,1) == 1
        theta = theta(:); % ensure column vector
    end

    % DH parameters: [a, alpha, d, theta_offset]
    DH = [ 0      -pi/2   0.333   0;
           0       pi/2   0       0;
           0.0825  pi/2   0.316   0;
          -0.0825 -pi/2   0       0;
           0       pi/2   0.384   0;
           0.088   pi/2   0       0;
           0       0      0.107   0 ];
    
    % Initialize transform
    T = eye(4);
    
    for i = 1:7
        a = DH(i,1); 
        alpha = DH(i,2); 
        d = DH(i,3);
        th = theta(i) + DH(i,4);
        
        % DH homogeneous transform
        A = [cos(th), -sin(th)*cos(alpha),  sin(th)*sin(alpha), a*cos(th);
             sin(th),  cos(th)*cos(alpha), -cos(th)*sin(alpha), a*sin(th);
             0,        sin(alpha),          cos(alpha),          d;
             0,        0,                   0,                   1];
         
        T = T * A;
    end
end
