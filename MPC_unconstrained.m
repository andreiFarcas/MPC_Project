clc
% System Matrices - taken from problem5_stochastic
% Discrete time state space model of the modified 4 tank system
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

% Simulation Parameters
num_steps = 100; % simulation steps
nx = size(Ff, 1); % nr. of states
ny = size(Cc, 1); % nr of outputs
nu = size(Gg, 2); % nr of inputs

% Noise Covariances
Q = 0.001 * eye(size(Gv, 2)); % Process noise 
R = 0.0001 * eye(ny); % Measurement noise 

% Initialize parameters
x_true = zeros(nx, 1); % True state
x_hat = zeros(nx, 1); % Estimated state (using Kalman filter)
u_prev = zeros(nu, 1); % Previous control input
N = 20; % prediction horizon in time steps
dist = zeros(size(Gv, 2) * N, 1); % Disturbance (zero for now)

%z_ref = zeros(ny * N, 1); % Reference trajectory (setpoint tracking)
% Reference will have 3 different setpoints
z_ref_full = zeros(ny, num_steps);
for k = 1:num_steps
    if k <= 30
        z_ref_full(:, k) = [1; 0];
    elseif k <= 60
        z_ref_full(:, k) = [0; 1];
    else
        z_ref_full(:, k) = [0.5; 0.5];
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
    z_ref = repmat(z_ref_full(:, k), N, 1); % Extend reference to prediction horizon
    
    % Dynamic Kalman Filter
    [x_hat, P_dynamic] = dynamicKalmanFilter(Ff, Gg, Gv, Cc, Q, R, y, u_prev, x_hat, P_dynamic);
    
    % Use the regulator function to compute the optimal control
    u0 = Regulator(Ff, Gg, Gv, Cc, Dd, N, x_hat, u_prev, z_ref, dist);
    
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
    plot(time, z_ref_full(i, :), 'r--', 'LineWidth', 1.5); hold on;
    plot(time, y_hist(i, :), 'b', 'LineWidth', 1.5);
    title(['Output y_', num2str(i)]);
    legend('Reference', 'System Output');
    xlabel('Time Step');
    ylabel(['y_', num2str(i)]);
    axis([0,num_steps,-3,3])
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
    axis([0,num_steps,-1,10])
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