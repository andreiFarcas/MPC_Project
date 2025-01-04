%% 
clc;
clear variables; close all

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
Ts = 10; % [s] sampling time
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

N = 10; % prediction horizon in time steps

% Penalisation matrices Qz and Hs
Qz = 3e5*blkdiag(eye(nu), zeros(nu*(N-1)));

central_diag = 2 * ones(1, nu*N); 
central_diag(N) = 1;
inf_diag = -1 * ones(1, nu*N-1);       
sup_diag = -1 * ones(1, nu*N-1);        
Hs = diag(central_diag) + diag(inf_diag, -1) + diag(sup_diag, 1);
Hs(nu*N,nu*N)=1; % shape (nuN x nuN) 
Hs = Hs*1;

% Noise covariances
Q = 0.01*eye(size(Gv, 2)); % Process noise 
R = 0.0001*eye(ny); % Measurement noise 

% Initialize parameters
x_true = zeros(nx, 1); % True state
x_hat = zeros(nx, 1); % Estimated state (using Kalman filter)
u_prev = zeros(nu, 1); % Previous control input
dist = zeros(size(Gv, 2) * N, 1); % Disturbance (zero for now)

% Reference will have 3 different setpoints
z_ref_full = zeros(ny, num_steps);
for k = 1:num_steps+N
    if k <= 30
        z_ref_full(:, k) = [3; 0];
    elseif k <= 60
        z_ref_full(:, k) = [0; 1];
    else
        z_ref_full(:, k) = [-0.5; 0.5];
    end
end

P_dynamic = eye(nx);

% Simulation loop
for k = 1:num_steps
    % True system dynamics
    x_true = Ff * x_true + Gg * u_prev;
    y = Cc * x_true; % Measurement
    
    % Reference for this step
    z_ref = z_ref_full(:,k);
    for i = 1:N-1
        z_ref = [z_ref;z_ref_full(:, k+i)]; % Extend reference to prediction horizon
    end

    % Dynamic Kalman Filter
    [x_hat, P_dynamic] = dynamicKalmanFilter(Ff, Gg, Gv, Cc, Q, R, y, u_prev, x_hat, P_dynamic);
    
    % Use the regulator function to compute the optimal control
    u0 = Regulator(Ff, Gg, Gv, Cc, Dd, N, x_hat, u_prev, z_ref, dist, Qz, Hs);
    
    % Apply control input to the system
    u_prev = u0; % Save for the next iteration
    
    % Store results for plotting
    x_true_hist(:, k) = x_true;
    x_hat_hist(:, k) = x_hat;
    y_hist(:, k) = y;
    u_hist(:, k) = u0;
end

% Plot results
time = 1:num_steps;

% Outputs vs Reference
figure;
for i = 1:ny
    subplot(ny, 1, i);
    plot(time, z_ref_full(i, 1:num_steps), 'r--', 'LineWidth', 1.5); hold on;
    plot(time, y_hist(i, :), 'b', 'LineWidth', 1.5);
    title(['Output y_', num2str(i)]);
    legend('Reference', 'System Output');
    xlabel('Time Step');
    ylabel(['y_', num2str(i)]);
    % axis([0,num_steps,-3,3])
    grid on;
end


figure;
for i = 1:nx
    subplot(nx, 1, i);
    plot(time, x_true_hist(i, :), 'b', 'LineWidth', 1.5); hold on;
    plot(time, x_hat_hist(i, :), 'r--', 'LineWidth', 1.5);
    title(['State x_', num2str(i)]);
    legend('True State', 'Estimated State');
    xlabel('Time Step');
    ylabel(['x_', num2str(i)]);
    % axis([0,num_steps,-1,10])
    grid on;
end

figure;
for i = 1:nu
    subplot(nu, 1, i);
    plot(time, u_hist(i, :), 'k', 'LineWidth', 1.5);
    title(['Control Input u_', num2str(i)]);
    xlabel('Time Step');
    ylabel(['u_', num2str(i)]);
    grid on;
end

%% Problem5_stochastic
% System Matrices - taken from problem5_stochastic
Ff = [0.9500, 0, 0.0927, 0;
      0, 0.9570, 0, 0.0803;
      0, 0, 0.9048, 0;
      0, 0, 0, 0.9179]; % State transition matrix
Gg = [2.8509, 0.0695;
      0.0800, 3.3382;
      0, 1.3916;
      1.8686, 0]; % Input matrix
Gv = [0.2318, 0;
      0, 0.1999;
      4.6386, 0;
      0, 4.6715]; % Process noise matrix
Cc = [0.0026, 0, 0, 0;
      0, 0.0026, 0, 0]; % Measurement matrix
Dd = 0; % No direct feedthrough

N = 10; % prediction horizon in time steps

% Penalisation matrices Qz and Hs
Qz = 3e5*blkdiag(eye(nu), zeros(nu*(N-1)));

central_diag = 2 * ones(1, nu*N); 
central_diag(N) = 1;
inf_diag = -1 * ones(1, nu*N-1);       
sup_diag = -1 * ones(1, nu*N-1);        
Hs = diag(central_diag) + diag(inf_diag, -1) + diag(sup_diag, 1);
Hs(nu*N,nu*N)=1; % shape (nuN x nuN) 
Hs = Hs*1;

% Noise Covariances
Q = 0.001 * eye(size(Gv, 2)); % Process noise 
R = 0.0001 * eye(ny); % Measurement noise 

% Initialize parameters
x_true = zeros(nx, 1); % True state
x_hat = zeros(nx, 1); % Estimated state (using Kalman filter)
u_prev = zeros(nu, 1); % Previous control input
dist = zeros(size(Gv, 2) * N, 1); % Disturbance (zero for now)

% Reference will have 3 different setpoints
z_ref_full = zeros(ny, num_steps);
for k = 1:num_steps+N
    if k <= 30
        z_ref_full(:, k) = [0; 2];
    elseif k <= 60
        z_ref_full(:, k) = [2; 1];
    else
        z_ref_full(:, k) = [0; 0.5];
    end
end

P_dynamic = eye(nx);

% Simulation loop
for k = 1:num_steps
    % True system dynamics
    w_k = sqrt(Q) * randn(size(Gv, 2), 1); % Process noise
    v_k = sqrt(R) * randn(ny, 1); % Measurement noise
    x_true = Ff * x_true + Gg * u_prev + Gv * w_k;
    y = Cc * x_true + v_k; % Measurement
    
    % Reference for this step
    z_ref = z_ref_full(:,k);
    for i = 1:N-1
        z_ref = [z_ref;z_ref_full(:, k+i)]; % We need to extend reference to prediction horizon
    end
    % Dynamic Kalman Filter
    [x_hat, P_dynamic] = dynamicKalmanFilter(Ff, Gg, Gv, Cc, Q, R, y, u_prev, x_hat, P_dynamic);
    
    % Use the regulator function to compute the optimal control
    u0 = Regulator(Ff, Gg, Gv, Cc, Dd, N, x_hat, u_prev, z_ref, dist, Qz, Hs);
    
    % Apply control input to the system
    u_prev = u0; % Save for the next iteration
    
    % Store results for plotting
    x_true_hist(:, k) = x_true;
    x_hat_hist(:, k) = x_hat;
    y_hist(:, k) = y;
    u_hist(:, k) = u0;
end

% Plot results
time = 1:num_steps;

% Outputs vs Reference
figure;
for i = 1:ny
    subplot(ny, 1, i);
    plot(time, z_ref_full(i, 1:num_steps), 'r--', 'LineWidth', 1.5); hold on;
    plot(time, y_hist(i, :), 'b', 'LineWidth', 1.5);
    title(['Output y_', num2str(i)]);
    legend('Reference', 'System Output');
    xlabel('Time Step');
    ylabel(['y_', num2str(i)]);
    grid on;
end


figure;
for i = 1:nx
    subplot(nx, 1, i);
    plot(time, x_true_hist(i, :), 'b', 'LineWidth', 1.5); hold on;
    plot(time, x_hat_hist(i, :), 'r--', 'LineWidth', 1.5);
    title(['State x_', num2str(i)]);
    legend('True State', 'Estimated State');
    xlabel('Time Step');
    ylabel(['x_', num2str(i)]);
    % axis([0,num_steps,-1,10])
    grid on;
end

figure;
for i = 1:nu
    subplot(nu, 1, i);
    plot(time, u_hist(i, :), 'k', 'LineWidth', 1.5);
    title(['Control Input u_', num2str(i)]);
    xlabel('Time Step');
    ylabel(['u_', num2str(i)]);
    grid on;
end