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
Q = 0.01 * eye(size(Gv, 2)); % Process noise covariance
R = 0.00001 * eye(ny); % Measurement noise covariance (obs smaller values improve the static model)

% Initial Conditions
x_true = zeros(nx, 1); % will store the real states
x_hat_static = zeros(nx, 1); % will store the estimated states (Static Kalman Filter)
x_hat_dynamic = zeros(nx, 1); % will store the estimated states (Dynamic Kalman Filter)
P_static = eye(nx); % State covariance (Static Kalman Filter)
P_dynamic = eye(nx); % State covariance (Dynamic Kalman Filter)

% Input Signal
u = 300*randn(nu, num_steps); % Random inputs

v = sqrt(R) * randn(ny, num_steps); % Measurement noide

w = sqrt(Q) * randn(size(Gv, 2), num_steps); % Process noise

% Storage for Results
x_true_hist = zeros(nx, num_steps);
x_hat_static_hist = zeros(nx, num_steps);
x_hat_dynamic_hist = zeros(nx, num_steps);
y_hist = zeros(ny, num_steps);

% Precompute Static Kalman Gain for static Kalman Filter
P_static = care(Ff', Cc', Gv * Q * Gv', R); 
% Compute steady-state Kalman gain
K_static = P_static * Cc' / (Cc * P_static * Cc' + R);

% Disturbance Matrix
Gd = [0.1, 0;
      0, 0.1;
      0.05, 0.05;
      0.05, 0.1]; % Disturbance input matrix

% Disturbance Covariance
D_cov = 0.001 * eye(size(Gd, 2)); % Disturbance covariance matrix

d = zeros(size(Gd, 2), num_steps); % Initialize disturbance matrix

% Generate disturbance signal with a step change after time step 50
for k = 1:num_steps
    if k > 50
        % After time = 50 we introduce a step in the disturbance
        d(:, k) = 50 * randn(size(Gd, 2), 1) + 300;
    else
        % Before time step 50, the disturbance is small
        d(:, k) = 50 * randn(size(Gd, 2), 1);
    end
end

% Simulation Loop
for k = 1:num_steps
    % True System Dynamics with Disturbances
    x_true = Ff * x_true + Gg * u(:, k) + Gv * w(:, k) + Gd * d(:, k);
    y = Cc * x_true + v(:, k); % Measurement
    
    % Static Kalman Filter (Using Precomputed Kalman Gain)
    x_hat_static = Ff * x_hat_static + Gg * u(:, k); % Prediction
    x_hat_static = x_hat_static + K_static * (y - Cc * x_hat_static); % Update
    
    % Dynamic Kalman Filter
    % Prediction
    x_hat_dynamic = Ff * x_hat_dynamic + Gg * u(:, k);
    P_dynamic = Ff * P_dynamic * Ff' + Gv * Q * Gv';
    
    % Kalman Gain Computation
    K_dynamic = P_dynamic * Cc' / (Cc * P_dynamic * Cc' + R);
    
    % Update
    x_hat_dynamic = x_hat_dynamic + K_dynamic * (y - Cc * x_hat_dynamic);
    P_dynamic = (eye(nx) - K_dynamic * Cc) * P_dynamic;
    
    % Store Results
    x_true_hist(:, k) = x_true;
    x_hat_static_hist(:, k) = x_hat_static;
    x_hat_dynamic_hist(:, k) = x_hat_dynamic;
    y_hist(:, k) = y;
end

% Plot Results
time = 1:num_steps;

figure;
for i = 1:nx
    subplot(nx, 1, i);
    plot(time, x_true_hist(i, :), 'b', 'LineWidth', 1.5); hold on;
    plot(time, x_hat_static_hist(i, :), 'g--', 'LineWidth', 1.5); % Static Kalman Filter
    plot(time, x_hat_dynamic_hist(i, :), 'r--', 'LineWidth', 1.5); % Dynamic Kalman Filter
    title(['State x_', num2str(i)]);
    legend('True State', 'Static KF', 'Dynamic KF');
    xlabel('Time Step');
    ylabel(['x_', num2str(i)]);
    grid on;
end

% Plot Disturbances
figure;
time = 1:num_steps;
for i = 1:size(d, 1)
    subplot(size(d, 1), 1, i);
    plot(time, d(i, :), 'LineWidth', 1.5);
    title(['Disturbance d_', num2str(i)]);
    xlabel('Time Step');
    ylabel(['d_', num2str(i)]);
    grid on;
end