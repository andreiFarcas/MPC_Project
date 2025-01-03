%% Functions
function [x_est_next, P_next] = EKF_ModifiedFourTankSystem(x_hat, u, d, p, Q, R, Ts, P)
    % x_hat : current estimated state (4x1)
    % u : control input vector (mx1)
    % d : disturbances (q x 1)
    % p : system parameters
    % Q : process noise covariance
    % R : measurement noise covariance
    % Ts : sampling time

    Gv = [0.2325, 0; 0, 0.2021; 4.6522, 0; 0, 4.6835];

    % State Prediction
    Ak = Jacobian_F(x_hat, u, d, p); % Jacobian of f(x,u)
    x_pred = x_hat + Ts * ModifiedFourTankSystem(0, x_hat, u, d, p); % Predicted state
    P_pred = Ak * P * Ak' + Gv * Q * Gv'; % Predicted covariance

    % Calculate true water levels
    h_true = x_pred ./ (p(12) * p(5:8)); % Predicted liquid levels
    y_meas = h_true(1:2) + randn(2, 1) .* sqrt(diag(R)); % Noisy measurements

    % State estimation using EKF
    H = diag(1 ./ (p(12) * p(5:8))); % Measurement model Jacobian
    H = H(1:2,:);
    K = P_pred * H' / (H * P_pred * H' + R); % Kalman gain
    y_pred = x_pred ./ (p(12) * p(5:8)); % Predicted measurement
    y_pred = y_pred(1:2);
    x_est_next = x_pred + K * (y_meas - y_pred); % Estimated state for the next time step
    P_next = (eye(4) - K * H) * P_pred; % Update covariance
end

function F = Jacobian_F(x, u, d, p)
    % Compute Jacobian matrix of the system dynamics
    m = x; % State vector
    g = p(11); % Gravity
    rho = p(12); % Density
    A = p(5:8); % Tank cross-sectional areas
    a = p(1:4); % Pipe cross-sectional areas
    h = m ./ (rho * A); % Liquid levels

    % Partial derivatives of outflow terms
    dqout_dm = a .* sqrt(2 * g ./ h) .* (-1 ./ (2 * rho * A));

    % Construct Jacobian matrix
    F = diag([dqout_dm(1) + dqout_dm(3), dqout_dm(2) + dqout_dm(4), -dqout_dm(3), -dqout_dm(4)]);
end


%%
p = [0.1; 0.1; 0.1; 0.1;  % Pipe cross-sectional areas [cm^2]
               50; 50; 50; 50;      % Tank cross-sectional areas [cm^2]
               0.5; 0.5;            % Valve positions [-]
               981;                 % Gravity [cm/s^2]
               1];                  % Density of water [g/cm^3]

t0 = 0.0;   % [s] Initial time
tf = 20*60; % [s] Final time
Ts = 0.1;            % Sampling time [s]
time = 0:Ts:tf;   % Time vector
num_steps = 100;
%% Estimator
x_hat = [10;10;10;10];
u = [200;200];
d = [100;100];
Q = eye(size(4, 2)); % Process noise 
R = eye(2); % Measurement noise 
Ts = 0.1;
P = eye(4);

x_true = x_hat;
x_sim = zeros(4, length(time));
x_est = zeros(4, length(time));
x_sim(:, 1) = x_hat;
x_est(:, 1) = x_hat;
P_next = P;

for k = 1:length(100)
    x_true = x_true + Ts * ModifiedFourTankSystem(time(k), x_true, u, d, p);
    x_sim(:, k+1) = x_true;
    [x_est_next, P_next] = EKF_ModifiedFourTankSystem(x_hat, u, d, p, Q, R, Ts, P_next);
    x_est(:, k+1) = x_est_next;
end

% Plot results
figure;
for i = 1:4
    subplot(4, 1, i);
    plot(time,x_sim(i, :), 'b-', 'LineWidth', 1.5); hold on;
    plot(time,x_est(i, :), 'r--', 'LineWidth', 1.5);
    title(['Tank ', num2str(i), ' State']);
    legend('Estimated State');
    xlabel('Time [s]');
    ylabel('Mass [g]');
end

%% Cost function
function cost = computeCost(U, N, xinit, zref, p, Qz, Hs)
    x = xinit;  % Initial state
    cost = 0;   % Initialize cost

    % Loop over the prediction horizon
    for k = 1:N
        % Extract the control inputs for this step
        u = U((k-1)*2 + 1:k*2); % Assuming 2 control inputs

        % Simulate the system for one time step using ModifiedFourTankSystem
        d = zeros(2, 1); % Disturbances
        Ts = 1; % Time step size
        xdot = ModifiedFourTankSystem(0, x, u, d, p); % Call system dynamics
        x = x + Ts * xdot; % Update state

        % Compute output (liquid levels in tanks)
        y = x ./ (p(12) * p(5:8)); % Compute the outputs

        % Compute cost for deviation from reference using the weighted L2 norm
        cost = cost + (y - zref(:,k))' * Qz * (y - zref(:, k)); 
    end

    % Add a quadratic term for the control inputs to the cost
    U = U(:);
    cost = cost + 0.5 * U' * Hs * U; % Ensure dimensions match
end

% Initial conditions
N = 10;              
xinit = [1000; 1000; 1000; 1000];   % Initial mass in tanks [g]
uprev = [0; 0];                       % Previous control inputs [cm^3/s]
zref = [100; 100; 100; 100];          % Desired mass in tanks [g]
zref = repmat(zref,1,num_steps);

% Disturbance and weighting matrices
dist = zeros(2, 1);                   % Disturbance [cm^3/s]
Qz = diag([10, 10, 10, 10]);          % Output weighting matrix (penalizes errors in outputs)
weights = repmat([1; 1], N, 1);
Hs = diag(weights);

% Limites des contrôles
u_min = 0;                            % Minimum control input [cm^3/s]
u_max = 100;                          % Maximum control input [cm^3/s]

U = [100;100];
U = repmat([100,100],num_steps,1);

cost = computeCost(U(1:N,:), N, xinit, zref, p, Qz, Hs)

%% Optimisation of cost function
N = 10;
umin = zeros(N * 2, 1);         % [cm^3/s]
umax = 500 * ones(N * 2, 1);    % [cm^3/s]

function u_opt = Optimisation(U0, N, xinit, zref, p, Qz, Hs, umin, umax)

    % Define the cost function for fmincon
    cost_function = @(U) computeCost(U, N, xinit, zref, p, Qz, Hs);

    % Set options for fmincon
    options = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'interior-point');

    % Solve the optimization problem using fmincon
    u_opt = fmincon(cost_function, U0, [], [], [], [], umin, umax, [], options);

end

% Initial guess for controls (N steps, 2 controls per step)
U0 = repmat([10; 10], N, 1);

% Call the optimization function
u_opt = Optimisation(U0, N, xinit, zref, p, Qz, Hs, umin, umax);

% Display the result
disp('Optimal controls:');
disp(u_opt(1:2));

%% Estimator + Optimisation

% Reference will have 3 different setpoints
z_ref_full = zeros(2, num_steps);
for k = 1:num_steps
    if k <= 30
        z_ref_full(:, k) = [1; 0];
    elseif k <= 60
        z_ref_full(:, k) = [0; 1];
    else
        z_ref_full(:, k) = [0.5; 0.5];
    end
end

x_hat = [10;10;10;10];
u = [200;200];
d = [100;100];
Q = eye(size(4, 2)); % Process noise 
R = eye(2); % Measurement noise 
Ts = 0.1;
P = eye(4);

x_true = x_hat;
x_sim = zeros(4, length(time));
x_est = zeros(4, length(time));
x_sim(:, 1) = x_hat;
x_est(:, 1) = x_hat;
P_next = P;
U0 = repmat([10; 10], N, 1);

for k = 1:length(time)
    % True systems dynamics
    x_true = x_true + Ts * ModifiedFourTankSystem(time(k), x_true, u, d, p);
    x_sim(:, k+1) = x_true;

    % Extended Kalman Filter
    [x_est_next, P_next] = EKF_ModifiedFourTankSystem(x_hat, u, d, p, Q, R, Ts, P_next);
    x_est(:, k+1) = x_est_next;

    % Use the regulator function to compute the optimal control
    try
        u_opt = Optimisation(U0, N, x_est_next, zref, p, Qz, Hs, umin, umax);
    catch ME
            disp('Error in NMPC:');
            disp(ME.message);
            return; % Exit the simulation or handle the error as needed
    end
    U0 = u_opt;

    % Store results for plotting
    x_true_hist(:, k) = x_true;
    x_hat_hist(:, k) = x_est_next;
    % y_hist(:, k) = y;
    u_hist(:, k) = U0(1:2);
end

%% Plots
nx = 4;
nu = 2;
% Plot results
time = 1:num_steps;



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

% % Outputs vs Reference
% figure;
% for i = 1:ny
%     subplot(ny, 1, i);
%     plot(time, z_ref_full(i, :), 'r--', 'LineWidth', 1.5); hold on;
%     plot(time, y_hist(i, :), 'b', 'LineWidth', 1.5);
%     title(['Output y_', num2str(i)]);
%     legend('Reference', 'System Output');
%     xlabel('Time Step');
%     ylabel(['y_', num2str(i)]);
%     % axis([0,num_steps,-3,3])
%     grid on;
% end

