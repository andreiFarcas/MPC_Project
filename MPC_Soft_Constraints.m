%% 
clc;
clear variables; close all;

%--------------------------------------------------------
% Parameters Definition
%--------------------------------------------------------
a1 = 1.2272; %[cm2] Area of outlet pipe 1
a2 = 1.2272;
a3 = 1.2272;
a4 = 1.2272;

A1 = 380.1327; %[cm2] Cross sectional area of tank 1
A2 = 380.1327; %[cm2] Cross sectional area of tank 2
A3 = 380.1327; %[cm2] Cross sectional area of tank 3
A4 = 380.1327; %[cm2] Cross sectional area of tank 4

gamma1 = 0.6; % Flow distribution constant. Valve 1
gamma2 = 0.7; % Flow distribution constant. Valve 2
g = 981;       %[cm/s2] The acceleration of gravity
rho = 1.00;    %[g/cm3] Density of water

p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];

F1 = 300; F2 = 300;
F3 = 100; F4 = 100;

m10 = 5000; m20 = 5000; m30 = 5000; m40 = 5000;

t0 = 0.0;   % [s] Initial time
tend = 20*60; % [s] Final time
Ts = 1; % [s] sampling time
t = t0:Ts:tend;               % [s] Sample instants


%% Problem5_deterministic
xs0 = [m10; m20; m30; m40]; % [g] Initial guess on xs
us = [F1; F2]; % [cm3/s] Flow rates
d = [F3;F4]; % model disturbances
v=[0;0;0;0]; % measurement noise

[~,~,~,A,B,Bv,C,D]  = StateSpaceModeling(p,xs0,us,d,v);

[~, ~,Ts,Ff,Gg,Gv,Cc,Dd] = DiscreteModeling(A,B,Bv,C,D,F1,F2,[]);

num_steps = 100; % simulation steps
nx = size(Ff, 1); % nr. of states
ny = size(Cc, 1); % nr of outputs
nu = size(Gg, 2); % nr of inputs

N = 2; % prediction horizon in time steps

% Penalisation matrice Qz and Hs
Qz = 5000000*blkdiag(eye(nu), zeros(nu*(N-1)));

central_diag = 2 * ones(1, nu*N); 
central_diag(N) = 1;
inf_diag = -1 * ones(1, nu*N-1);       
sup_diag = -1 * ones(1, nu*N-1);        
Hs = diag(central_diag) + diag(inf_diag, -1) + diag(sup_diag, 1);
Hs(nu*N,nu*N)=1; % shape (nuN x nuN) 
Hs = Hs*50;

% Soft constraint penalties
rho2 = 5e7; % L2 penalty weight
rho1 = 5e6;  % L1 penalty weight

% Noise covariances
Q = eye(size(Gv, 2)); % Process noise 
R = eye(ny); % Measurement noise 

% Initialize parameters
x_true = zeros(nx, 1); % True state
x_hat = zeros(nx, 1); % Estimated state (using Kalman filter)
u_prev = zeros(nu, 1); % Previous control input
dist = zeros(size(Gv, 2) * N, 1); % Disturbance (zero for now)

% Reference will have 3 different setpoints
z_ref_full = zeros(ny, num_steps);
for k = 1:num_steps
    if k <= 30
        z_ref_full(:, k) = [10; 0];
    elseif k <= 60
        z_ref_full(:, k) = [0; 10];
    else
        z_ref_full(:, k) = [5; 5];
    end
end

P_dynamic = eye(nx);

% Input constraints
u_min = [0; 0];     % Minimum input constraints
u_max = [500; 500]; % Maximum input constraints
delta_u_min = [-200; -200]; 
delta_u_max = [200; 200];
z_min = [0;0];
z_max = [6;8];

% Simulation loop
for k = 1:num_steps
    % True system dynamics
    x_true = Ff * x_true + Gg * u_prev;
    y = Cc * x_true; % Measurement
    
    % Reference for this step
    z_ref = repmat(z_ref_full(:, k), N, 1); % Extend reference to prediction horizon
    
    % Dynamic Kalman Filter
    [x_hat, P_dynamic] = dynamicKalmanFilter(Ff, Gg, Gv, Cc, Q, R, y, u_prev, x_hat, P_dynamic);
    
    % Use the regulator function to compute the optimal control
    % Input-Soft Output
    u0 = Regulator_Soft_Constraints(Ff, Gg, Gv, Cc, Dd, N, x_hat, u_prev, z_ref, dist, Qz, Hs, ...
        u_min, u_max, delta_u_min, delta_u_max, z_min, z_max, rho2, rho1);

    % Apply control input to the system
    u_prev = u0; % Save for the next iteration
    
    % Store results for plotting
    x_true_hist(:, k) = x_true;
    x_hat_hist(:, k) = x_hat;
    y_hist(:, k) = y;
    u_hist(:, k) = u0;
end

%% Comparison Input-Output vs Input-Soft Output
% Initialize parameters
x_true = zeros(nx, 1); % True state
x_hat = zeros(nx, 1); % Estimated state (using Kalman filter)
u_prev = zeros(nu, 1); % Previous control input
dist = zeros(size(Gv, 2) * N, 1); % Disturbance (zero for now)

% Simulation loop
for k = 1:num_steps
    % True system dynamics
    x_true = Ff * x_true + Gg * u_prev;
    y = Cc * x_true; % Measurement
    
    % Reference for this step
    z_ref = repmat(z_ref_full(:, k), N, 1); % Extend reference to prediction horizon
    
    % Dynamic Kalman Filter
    [x_hat, P_dynamic] = dynamicKalmanFilter(Ff, Gg, Gv, Cc, Q, R, y, u_prev, x_hat, P_dynamic);
    
    % Use the regulator function to compute the optimal control
    u0 = Regulator_IO_Constraints(Ff, Gg, Gv, Cc, Dd, N, x_hat, u_prev, z_ref, dist, Qz, Hs, u_min, u_max,delta_u_min,delta_u_max,z_min,z_max);
    
    % Apply control input to the system
    u_prev = u0; % Save for the next iteration
    
    % Store results for plotting
    x_true_hist2(:, k) = x_true;
    x_hat_hist2(:, k) = x_hat;
    y_hist2(:, k) = y;
    u_hist2(:, k) = u0;
end

% Plot results
time = 1:num_steps;

% Plot state trajectories
figure;
for i = 1:nx
    subplot(nx, 1, i);
    plot(time, x_true_hist(i, :), 'b', 'LineWidth', 1.5); hold on;
    plot(time, x_hat_hist(i, :), 'r--', 'LineWidth', 1.5);
    title(['State x_', num2str(i)]);
    legend('True State', 'Estimated State');
    xlabel('Time Step');
    ylabel(['x_', num2str(i)]);
    grid on;
end

% Plot Outputs vs Reference
figure;
for i = 1:ny
    subplot(ny, 1, i);
    plot(time, z_ref_full(i, :), 'r--', 'LineWidth', 1.5); hold on;
    plot(time, y_hist(i, :), 'b', 'LineWidth', 1.5); hold on;
    plot(time, y_hist2(i, :), 'g', 'LineWidth', 1.5);hold on;
    plot(time, z_max(i) * ones(1, length(time)), 'k--', 'LineWidth', 1);
    title(['Output y_', num2str(i)]);
    legend('Reference', 'Input-Soft Output Constraints', 'Input-Output Constraints');
    xlabel('Time Step');
    ylabel(['y_', num2str(i)]);
    grid on;
end

% Plot Control Inputs
figure;
for i = 1:nu
    subplot(nu, 1, i);
    plot(time, u_hist(i, :), 'b', 'LineWidth', 1.5); hold on
    plot(time, u_hist2(i, :),  'g','LineWidth', 1.5);
    title(['Control Input u_', num2str(i)]);
    legend('Input-Soft Output Constraints', 'Input Output Constraints');
    xlabel('Time Step');
    ylabel(['u_', num2str(i)]);
    grid on;
end
