clc; 
clear; 
clf;

filename = 'controller.csv';

% Read CSV file, skipping first row, no headers
data = readmatrix(filename);
data(1,:) = []; % Skip first row

states_data = data(:, 1:2);
input_data = data(:, 3:4);

T = 10.0;
dt = 1e-3;
num_steps = round(T / dt);
x = [0.1, 0.2];
xi_VCZ = zeros(num_steps, 2);
u_VCZ = zeros(num_steps, 2);
EE_VCZ = zeros(num_steps, 2);
time_arr = 0:dt:T-dt;

for step = 1:num_steps
    diffs = vecnorm(states_data - x, 2, 2);
    [~, closest_idx] = min(diffs);
    u = input_data(closest_idx, :);

    xi_VCZ(step, :) = x;
    u_VCZ(step, :) = u;
    EE_VCZ(step, :) = VCZ_kinematics(x(1), x(2));

    x = VCZ_dynamics(x, u, dt);
end

%% VCZ control
% xi = xi;
th0 = [0.1, 0.2];
th = zeros(length(xi_VCZ),2);
om = zeros(length(xi_VCZ),2);
om_d = zeros(length(xi_VCZ),2);
torque = zeros(length(xi_VCZ),2);
EE = zeros(length(xi_VCZ),2);

th(1,:) = th0;
EE(1, :) = VCZ_kinematics(th(1,1), th(1,2));

lam1 = 0.02;
lam2 = 0.1;
% taub = 10;
% omb = 0.2;
taub = 10;
omb = 0.25;

for step=2:num_steps
    q = [th(step-1,:)'; om(step-1,:)'];
    dq = RR_ode(q, torque(step-1,:));
    
    th(step,:) = th(step-1,:) + dt*dq(1:2)';
    om(step,:) = om(step-1,:) + dt*dq(3:4)';
    
    et = norm(th(step,:)-xi_VCZ(step,:))/lam1;
    om_d(step,:) = -omb*psi(et).*(th(step,:)-xi_VCZ(step,:))/norm(th(step,:)-xi_VCZ(step,:));

    eo = norm(om(step,:)-om_d(step,:))/lam2;
    torque(step,:) = -taub*psi(eo).*(om(step,:)-om_d(step,:))/norm(om(step,:)-om_d(step,:));

    EE(step, :) = VCZ_kinematics(th(step,1), th(step,2));
end

% %% SCOTS
% filename = 'scots.csv';
% % Read CSV file and skip first row
% data = readmatrix(filename);
% data(1,:) = []; % remove first row
% 
% states_data = data(:, 1:4);
% input_data  = data(:, 5:6);
% 
% % Simulation parameters
% T = 10.0;
% % dt = 0.01;
% num_steps = round(T / dt);
% x = [0.1, 0.2, 0.0, 0.0]; % theta1, theta2, omega1, omega2
% 
% SCOTS_th     = zeros(num_steps, 4);
% SCOTS_torque   = zeros(num_steps, 2);
% SCOTS_EE  = zeros(num_steps, 2);
% SCOTS_time_arr   = 0:dt:(T - dt);
% 
% for step = 1:num_steps
%     diffs = vecnorm(states_data - x, 2, 2);
%     [~, closest_idx] = min(diffs);
%     u = input_data(closest_idx, :);
% 
%     SCOTS_th(step, :)    = x;
%     SCOTS_torque(step, :)  = u;
%     SCOTS_EE(step, :) = forward_kinematics(x(1), x(2));
% 
%     x = forward_dynamics(x, u, dt);
% end

%% Plots
figure(1);
clf;
% Joint Angle
subplot(2,3,1)
hold on;
rectangle('Position',[0,0.7,10,0.1],'FaceColor',[0 1 0], 'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off')
plot(time_arr, xi_VCZ(:,1), 'DisplayName', '$\xi_1(t)$ - VCZ Center', 'Color', 'k', 'LineStyle', ':', 'LineWidth', 3); 
plot(time_arr, th(:,1), 'DisplayName', '$\theta_1(t)$ - Proposed VCZ method', 'Color', 'b', 'LineWidth', 3); 
% plot(SCOTS_time_arr, SCOTS_th(:,1), 'DisplayName', '$\theta_1(t)$ - SCOTS Controller', 'Color', 'r', 'LineStyle', '--', 'LineWidth', 3); 
set(gca,'FontSize', 18, 'FontName', 'Latin Modern Roman', 'TickLabelInterpreter', 'latex')
xlabel('$$t(s)$$', 'FontSize', 18,Interpreter ='latex');
ylabel('$$\xi_1, \theta_1(rad)$$', 'FontSize', 18,Interpreter ='latex');
legend('FontName', 'Latin Modern Roman','Interpreter', 'latex','Location','best')
grid on; box on;

subplot(2,3,4)
hold on;
rectangle('Position',[0,-0.8,10,0.1],'FaceColor',[0 1 0], 'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
plot(time_arr, xi_VCZ(:,2), 'DisplayName', '$\xi_2(t) - VCZ Center$', 'Color', 'k', 'LineStyle', ':', 'LineWidth', 3);
plot(time_arr, th(:,2), 'DisplayName', '$\theta_2(t)$ - Proposed VCZ method', 'Color', 'b', 'LineWidth', 3); 
% plot(SCOTS_time_arr, SCOTS_th(:,2), 'DisplayName', '$\theta_2(t)$ - SCOTS Controller', 'Color', 'r', 'LineStyle', '--', 'LineWidth', 3);
set(gca,'FontSize', 18, 'FontName', 'Latin Modern Roman', 'TickLabelInterpreter', 'latex')
xlabel('$$t(s)$$', 'FontSize', 18,Interpreter ='latex');
ylabel('$$\xi_2, \theta_2(rad)$$', 'FontSize', 18,Interpreter ='latex');
legend('FontName', 'Latin Modern Roman','Interpreter', 'latex','Location','best')
grid on; box on;

% Joint Velocity
subplot(2,3,2);
hold on;
% plot(time_arr, u_VCZ(:,1), 'DisplayName', '$u_1$', 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
plot(time_arr, om(:,1), 'DisplayName', '$\omega_1(t)$ - Proposed VCZ method', 'Color', 'b', 'LineWidth', 2); 
% plot(SCOTS_time_arr, SCOTS_th(:,3), 'DisplayName', '$\omega_1$ - SCOTS Controller', 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
set(gca,'FontSize',18, 'FontName', 'Latin Modern Roman', 'TickLabelInterpreter', 'latex')
xlabel('$$t(s)$$', 'FontSize', 18,Interpreter ='latex');
ylabel('$$u_1, \omega_1(rad)$$', 'FontSize', 18,Interpreter ='latex');
legend('FontName', 'Latin Modern Roman','Interpreter', 'latex','Location','best')
grid on; box on;

subplot(2,3,5);
hold on;
% plot(time_arr, u_VCZ(:,2), 'DisplayName', '$u_2$', 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
plot(time_arr, om(:,2), 'DisplayName', '$\omega_2(t)$ - Proposed VCZ method', 'Color', 'b', 'LineWidth', 2); 
% plot(SCOTS_time_arr, SCOTS_th(:,4), 'DisplayName', '$\omega_2$ - SCOTS Controller', 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
set(gca,'FontSize',18, 'FontName', 'Latin Modern Roman', 'TickLabelInterpreter', 'latex')
xlabel('$$t(s)$$', 'FontSize', 18,Interpreter ='latex');
ylabel('$$u_2, \omega_2(rad)$$', 'FontSize', 18,Interpreter ='latex');
legend('FontName', 'Latin Modern Roman','Interpreter', 'latex','Location','best')
grid on; box on;

% Control Torque
% subplot(2,3,3);
% hold on;
% plot(time_arr, torque(:,1), 'DisplayName', '$\tau_1$ - Proposed VCZ method', 'LineWidth', 2); 
% plot(time_arr, torque(:,2), 'DisplayName', '$\tau_2$ - Proposed VCZ method', 'LineWidth', 2);
% plot(SCOTS_time_arr, SCOTS_torque(:,1), 'DisplayName', '$\tau_1$ - SCOTS Controller', 'LineWidth', 2); 
% plot(SCOTS_time_arr, SCOTS_torque(:,2), 'DisplayName', '$\tau_2$ - SCOTS Controller', 'LineWidth', 2);
% set(gca,'FontSize',18, 'FontName', 'Latin Modern Roman', 'TickLabelInterpreter', 'latex')
% xlabel('$$t(s)$$', 'FontSize', 18,Interpreter ='latex');
% ylabel('$$\tau_1, \tau_2(Nm)$$', 'FontSize', 18,Interpreter ='latex');
% legend('FontName', 'Latin Modern Roman','Interpreter', 'latex','Location','best')
% grid on; box on;

subplot(2,3,3);
hold on;
plot(time_arr, torque(:,1), 'DisplayName', '$\tau_1$ - Proposed VCZ method', 'Color', 'b', 'LineWidth', 2); 
% plot(SCOTS_time_arr, SCOTS_torque(:,1), 'DisplayName', '$\tau_1$ - SCOTS Controller', 'Color', 'r', 'LineStyle', '--',  'LineWidth', 2); 
set(gca,'FontSize',18, 'FontName', 'Latin Modern Roman', 'TickLabelInterpreter', 'latex')
xlabel('$$t(s)$$', 'FontSize', 18,Interpreter ='latex');
ylabel('$$\tau_1(Nm)$$', 'FontSize', 18,Interpreter ='latex');
legend('FontName', 'Latin Modern Roman','Interpreter', 'latex','Location','best')
grid on; box on;

subplot(2,3,6);
hold on; 
plot(time_arr, torque(:,2), 'DisplayName', '$\tau_2$ - Proposed VCZ method', 'Color', 'b', 'LineWidth', 2);
% plot(SCOTS_time_arr, SCOTS_torque(:,2), 'DisplayName', '$\tau_2$ - SCOTS Controller', 'Color', 'r', 'LineStyle', '--',  'LineWidth', 2);
set(gca,'FontSize',18, 'FontName', 'Latin Modern Roman', 'TickLabelInterpreter', 'latex')
xlabel('$$t(s)$$', 'FontSize', 18,Interpreter ='latex');
ylabel('$$\tau_2(Nm)$$', 'FontSize', 18,Interpreter ='latex');
legend('FontName', 'Latin Modern Roman','Interpreter', 'latex','Location','best')
grid on; box on;

% % EE Trajectory
% subplot(2,3,6);
% hold on;
% % plot(EE_VCZ(:,1), EE_VCZ(:,2), 'DisplayName', 'End-Effector Path for VCZ', 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); hold on;
% plot(EE(:,1), EE(:,2), 'DisplayName', 'Proposed VCZ method', 'Color', 'b', 'LineWidth', 2);
% plot(SCOTS_EE(:,1), SCOTS_EE(:,2), 'DisplayName', 'SCOTS Controller', 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
% set(gca,'FontSize',18, 'FontName', 'Latin Modern Roman', 'TickLabelInterpreter', 'latex')
% % scatter(EE_VCZ(1,1), EE_VCZ(1,2), 50, 'red', 'filled', 'DisplayName', 'Start');
% % scatter(EE_VCZ(end,1), EE_VCZ(end,2), 50, 'green', 'filled', 'DisplayName', 'End');
% xlabel('$x (m)$','FontName', 'Latin Modern Roman','Interpreter', 'latex');
% ylabel('$y (m)$','FontName', 'Latin Modern Roman','Interpreter', 'latex');
% legend('FontName', 'Latin Modern Roman','Interpreter', 'latex','Location','best')
% grid on; box on;

disp('Simulation complete.');
disp('Script completed successfully!');


%%
% figure(2)
% err = zeros(length(xi_VCZ),1);
% 
% for step=1:num_steps
%     err(step) = norm(th(step,:) - xi_VCZ(step,:)); 
% end
% plot(time_arr, err)
%%

function p = psi(e)
    a = 5;
    % a = 3;
    % p = tanh(a*e);
    % p = tanh(a*e).*(1-exp(-(a*e).^2));
    p = tanh(a*e).^3;
end

% === Forward Dynamics ===
function x_next = VCZ_dynamics(x, u, dt)
    theta1 = x(1);
    theta2 = x(2);
    x_next = [theta1 + u(1)*dt, theta2 + u(2)*dt];
end

% === Forward Kinematics ===
function pos = VCZ_kinematics(theta1, theta2, l1, l2)
    if nargin < 3
        l1 = 1.0;
        l2 = 1.0;
    end
    x = l1 * cos(theta1) + l2 * cos(theta1 + theta2);
    y = l1 * sin(theta1) + l2 * sin(theta1 + theta2);
    pos = [x, y];
end

%%
function dxdt = RR_ode(x,u)
    m = 0.1; % Mass of each link
    l = 10; % Length of each link
    g = 0; %9.8;
    u = u';
    th1 = x(1);
    th2 = x(2);
    dth1 = x(3);
    dth2 = x(4);
    M = m*l^2*[5/3 + cos(th2), 1/3 + 1/2*cos(th2); ...
         1/3 + 1/2*cos(th2), 1/3];
    C = m*l^2*sin(th2)*[-0.5*dth2^2 - dth1*dth2; 0.5*dth1^2];
    D = m*g*l*[1.5*cos(th1) + 0.5*cos(th1+th2); 0.5*cos(th1+th2)];
        
    d = 0*ones(2,1);

    dth = [dth1; dth2];
    ddth = M\(- C - D + u + d);
    dxdt = [dth; ddth];
end

%%
%% 
function x_next = forward_dynamics(x, u, dt)
    theta1 = x(1);
    theta2 = x(2);
    omega1 = x(3);
    omega2 = x(4);

    M11 = (5.0 / 3.0) + cos(theta2);
    M12 = (1.0 / 3.0) + 0.5 * cos(theta2);
    M21 = M12;
    M22 = 1.0 / 3.0;

    C1 = -sin(theta2) * (omega1 * omega2 + 0.5 * omega2^2);
    C2 = 0.5 * sin(theta2) * omega1^2;

    detM = M11 * M22 - M12 * M21;
    if abs(detM) < 1e-10
        detM = 1e-10;
    end

    invM11 = M22 / detM;
    invM12 = -M12 / detM;
    invM21 = -M21 / detM;
    invM22 = M11 / detM;

    a1 = invM11 * (u(1) - C1) + invM12 * (u(2) - C2);
    a2 = invM21 * (u(1) - C1) + invM22 * (u(2) - C2);

    x_next = [...
        theta1 + omega1 * dt, ...
        theta2 + omega2 * dt, ...
        omega1 + a1 * dt, ...
        omega2 + a2 * dt ...
    ];
end

function pos = forward_kinematics(theta1, theta2, l1, l2)
    if nargin < 3
        l1 = 1.0;
        l2 = 1.0;
    end
    x = l1 * cos(theta1) + l2 * cos(theta1 + theta2);
    y = l1 * sin(theta1) + l2 * sin(theta1 + theta2);
    pos = [x, y];
end