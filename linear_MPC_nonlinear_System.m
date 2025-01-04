clear variables
clc
close all

%% Parameters Definition
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

%% Obtaining Steady-State Values
% Define the inputs (control variables)
F1 = 300;   % [cm^3/s] Flow rate from pump 1
F2 = 300;   % [cm^3/s] Flow rate from pump 2
F3 = 0;     % [cm^3/s] Flow rate disturbance to tank 3
F4 = 0;     % [cm^3/s] Flow rate disturbance to tank 4
u = [F1; F2];

% Initial guess for the steady-state state variables (liquid masses)
x0_guess = [5000; 5000; 5000; 5000];

% Wrap the function for fsolve
steadyStateFunction = @(x) ModifiedFourTankSystem(0, x, u, [F3; F4], p);

% Solve for steady state
options = optimoptions('fsolve', 'Display', 'iter', 'FunctionTolerance', 1e-8);
x_steady = fsolve(steadyStateFunction, x0_guess, options);

% Compute steady-state heights
A = p(5:8, 1); % Cross-sectional areas of tanks
h_steady = x_steady ./ (rho * A); % [cm] Tank levels at steady state

disp('Steady-state liquid masses (g):');
disp(x_steady);

disp('Steady-state liquid levels (cm):');
disp(h_steady);

%% Closed-Loop Simulations on the nonlinear model

% Simulation parameters
t_f = 20 * 60; % Final simulation time [s]
dt = 1;      % Time step [s]
T = 0:dt:t_f;  % Time vector
X = zeros(length(T)+10, length(x_steady));
X(1, :) = x_steady; % Initial conditions are considered the steady state values
u_prev = zeros(2, 1); % Initial input (zero)

% Reference tank levels [cm]
h1_ref = h_steady(1) * ones(1, length(T));  % Initialize with steady-state value
h2_ref = h_steady(2) * ones(1, length(T));  

% Apply 25% increase after 5 minutes (300 seconds)
h1_ref(T >= 300) = h_steady(1) + 0.25 * h_steady(1);  % 25% increase for h1
h2_ref(T >= 300) = h_steady(2) + 0.25 * h_steady(2);  % 25% increase for h2

% Parameters for system linearization for Kalman
xs0 = x_steady; % [g] Initial guess on xs
us = [F1; F2]; % [cm3/s] Flow rates
d = [F3;F4]; % model disturbances
v=[0;0;0;0]; % measurement noise
x_hat = x_steady; % Initial estimated state
P_dynamic = eye(4);

% Parameters for Regulator
nx = 4; % nr. of states
ny = 2; % nr of outputs
nu = 2; % nr of inputs
N = 10; % prediction horizon in time steps
dist = zeros(2 * N, 1); % Disturbance (zero for now)
Qz = 5000000*blkdiag(eye(nu), zeros(nu*(N-1)));
central_diag = 2 * ones(1, nu*N); 
central_diag(N) = 1;
inf_diag = -1 * ones(1, nu*N-1);       
sup_diag = -1 * ones(1, nu*N-1);        
Hs = diag(central_diag) + diag(inf_diag, -1) + diag(sup_diag, 1);
Hs(nu*N,nu*N)=1; % shape (nuN x nuN) 
Hs = Hs*50;

% Storage of results for plotting
h1_result = zeros(length(T));
h2_result = h1_result;
u1_result = h1_result;
u2_result = h1_result;
x_HAT = zeros(4, length(T));

% Define ODE function for ode15s
ode_fun = @(t, x) ModifiedFourTankSystem(t, x, u_prev, d, p);

for i = 2:length(T)-11
    t = T(i);

    % Compute real system states dx for this step
    
    %dx = ModifiedFourTankSystem(t, X(i, :), u_prev, d, p);
    %dx = dx(:); % ensures dx is a column vector
    %X(i + 1, :) = X(i, :) + dx' * dt; % Update the states

    t_span = [T(i), T(i+1)];  % Current time step for ode15s
    [t_out, x_out] = ode15s(ode_fun, t_span, X(i, :));  % Solve ODE over this time step
    X(i + 1, :) = x_out(end, :);  % Update the state with the result from ode15s

    areas = p(5:8);
    y = X(i, :) ./ (rho * areas);  % basically heights in the tanks
    y = y(1:2);

    % Reference heights for this step
    h1 = h1_ref(i:i + N-1);
    h2 = h2_ref(i:i + N-1);

    % Linearize at this operating point (For Kalman)
    [~,~,~,A,B,Bv,C,D]  = StateSpaceModeling(p,xs0,us,d,v);
    [~, ~,Ts,Ff,Gg,Gv,Cc,Dd] = DiscreteModeling(A,B,Bv,C,D,F1,F2,[]);
    
    Q = eye(size(Gv, 2)); % Process noise
    R = eye(2); % Measurement noise 

    % Use the Kalman to estimate the future state
    [x_hat, P_dynamic] = dynamicKalmanFilter(Ff, Gg, Gv, Cc, Q, R, y, u_prev, x_hat, P_dynamic);

    % Add Regulator
    % Use the regulator function to compute the optimal control
    u0 = Regulator(Ff, Gg, Gv, Cc, Dd, N, x_hat(:,1), u_prev, [h1 h2]', dist, Qz, Hs);

    % Apply control input to the system
    u_prev = u0; % Save for the next iteration

    % Save the data for plotting
    h1_result(i) = y(1);
    h2_result(i) = y(2);
    u1_result(i) = u0(1);
    u2_result(i) = u0(2);
    x_HAT(:, i) = x_hat(:, 1);
end

%% Plotting Tank Heights
h_result = [h1_result(:,1) h2_result(:,1)];
h_ref = [h1_ref; h2_ref]';
figure;
for i = 1:2
    subplot(2, 1, i);
    plot(T, h_ref(:, i), 'b', 'LineWidth', 1.5); hold on
    plot(T, h_result(:, i), 'r--', 'LineWidth', 1.5); 
    title(['Output y_', num2str(i)]);
    legend('Reference', 'System Output');
    xlabel('Time Step');
    ylabel(['y_', num2str(i)]);
    grid on;
end

%% Plotting States X vs Estimated States X_HAT
figure;

for i = 1:4
    subplot(2, 2, i);
    plot(T, X(1:end-10, i), 'b', 'LineWidth', 1.5); hold on; % Plot true states (X)
    plot(T, x_HAT(i, :), 'r--', 'LineWidth', 1.5);          % Plot estimated states (X_HAT)
    title(['State ', num2str(i)]);
    legend('True State X', 'Estimated State x_{hat}');
    xlabel('Time [s]');
    ylabel(['State ', num2str(i)]);
    grid on;
end