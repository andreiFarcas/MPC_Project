%% Initialization and Parameter Definition
clear variables

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

%% Obtaining the steady state values 

% This code is same as problem 3 but also includes F3 and F4 = 100 cm3/s
t0 = 0.0;   % [s] Initial time
t_f = 20*90; % [s] Final time
time = linspace(t0, t_f, 1000);

% Define a step of 10% after 60 seconds
F1 = 300*ones(size(time));
F2 = F1;
F3 = 100*ones(size(time));
F4 = F3;

u = [F1', F2'];

% Simulating
d = [F3', F4'];

% Define time vector
time = linspace(t0, t_f, length(u));

% Interpolate u and d based on current time (to be used in ode simulation)
u_interp = @(t) interp1(time, u, t, 'linear', 'extrap')';
d_interp = @(t) interp1(time, d, t, 'linear', 'extrap')';

% Initial guess for the steady-state state variables (liquid masses)
x0_guess = [2000; 2000; 500; 500];

% Wrap the function for fsolve
steadyStateFunction = @(x) ModifiedFourTankSystem(0, x, [300; 300], [100; 100], p);

% Solve for steady state
options = optimoptions('fsolve', 'Display', 'iter', 'FunctionTolerance', 1e-8);
x_steady = fsolve(steadyStateFunction, x0_guess, options);

% Compute steady-state heights
A = p(5:8, 1); % Cross-sectional areas of tanks
h_steady = x_steady ./ (rho * A); % [cm] Tank levels at steady state

% Display results
disp('Steady-state liquid masses (g):');
disp(x_steady);

disp('Steady-state liquid levels (cm):');
disp(h_steady);

%% Simulation with steps of 10%, 25%, 50% and variable noise level

% Generate noise
noise_level = 0;    % 0 to 10

% Simulating a step response for F1 = F2 = 300 cm3/s

% Define a step of 10% after 60 seconds
F1 = 300*ones(size(time));
F1(time >= 5*60) = 330;
F2 = F1;
F3 = 100*ones(size(time));
F4 = F3;

u = [F1', F2'];
x0 = x_steady;

% Define time vector
time = linspace(t0, t_f, length(u));
n = length(time);

% Measurement noise
if(noise_level > 0)
    R = noise_level * 0.1 * eye(4);
    Lr = chol(R,'lower');
    v = Lr*randn(4,n);  
    v = v';
else
    v = zeros(1000, 4);
end

% Process Noise
k = 6;                      % Number of samples per disturbance value

F3_disturbance = randi([-5*noise_level, 5*noise_level], ceil((n+1)/k), 1); % Disturbance for F3
F4_disturbance = randi([-5*noise_level, 5*noise_level], ceil((n+1)/k), 1); % Disturbance for F4
F3_disturbance = repelem(F3_disturbance, k, 1);
F4_disturbance = repelem(F4_disturbance, k, 1);

% Truncate to match time steps
F3_disturbance = F3_disturbance(1:n);
F4_disturbance = F4_disturbance(1:n);

% Combine baseline and disturbances
d = [F3 + F3_disturbance'; F4 + F4_disturbance'];  % Disturbances for F3 and F4
d = max(0, d)';                                   % Ensure non-negative flows

% Interpolate u and d based on current time (to be used in ode simulation)
u_interp = @(t) interp1(time, u, t, 'linear', 'extrap')';
d_interp = @(t) interp1(time, d, t, 'linear', 'extrap')';

% Solving the simulation 
ode_fun = @(t, x) ModifiedFourTankSystem(t, x, u_interp(t), d_interp(t), p);
[T, X] = ode15s(ode_fun, [t0 t_f], x0);

% Extracting parameters
A = p(5:8,1)';

% Compute the measured variables
H = X ./ (rho * A);
H_noisy = H + v(1:length(H), :);

% Normalize the step response
H_noisy_norm = (H_noisy - h_steady') / (10/100 * 300);

figure("Name", "Noisy Simulations")

subplot(331)
u = [u(:, 1), u(:, 2), d(:, 1), d(:, 2)];
plot(time, u, "LineWidth",1)
xlabel('Time (s)');
ylabel('Input Flow (cm^3/s)');
axis([-50, 1250, -50, 350])
title('Inputs with 10% increase on F1 and F2');
legend('F1', 'F2', 'F3', 'F4');
grid on;
axis([0, length(time), 0, 500])

subplot(332)
plot(T,H_noisy, "LineWidth",1)
xlabel('Time (s)');
ylabel('Height (cm)');
title('Height evolution of tanks level');
legend('h1', 'h2', 'h3', 'h4');
grid on;

subplot(333)
plot(T,H_noisy_norm, "LineWidth",1)
xlabel('Time (s)');
ylabel('Height (cm)');
title('Normalized height evolution of tanks level');
legend('h1', 'h2', 'h3', 'h4');
grid on;


% Define a step of 25% after 60 seconds
F1 = 300*ones(size(time));
F1(time >= 5*60) = 300 + 25/100 * 300;
F2 = F1;

u = [F1', F2'];

% Interpolate u based on current time (to be used in ode simulation)
u_interp = @(t) interp1(time, u, t, 'linear', 'extrap')';

% Solving the simulation 
ode_fun = @(t, x) ModifiedFourTankSystem(t, x, u_interp(t), d_interp(t), p);
[T, X] = ode15s(ode_fun, [t0 t_f], x0);

% Compute the measured variables
A = p(5:8,1)';
H = X ./ (rho * A);

if(noise_level > 0)
    R = noise_level * 0.25 * eye(4);
    Lr = chol(R,'lower');
    v = Lr*randn(4,n);  
    v = v';
end
H_noisy = H + v(1:length(H), :);

% Normalize the step response
H_noisy_norm = (H_noisy - h_steady') / (25/100 * 300);

subplot(334)
u = [u(:, 1), u(:, 2), d(:, 1), d(:, 2)];
plot(time, u, "LineWidth",1)
xlabel('Time (s)');
ylabel('Input Flow (cm^3/s)');
title('Inputs with 25% increase on F1 and F2');
legend('F1', 'F2', 'F3', 'F4');
grid on;
axis([0, length(time), 0, 500])

subplot(335)
plot(T, H_noisy, "LineWidth",1)
xlabel('Time (s)');
ylabel('Height (cm)');
title('Height evolution of tanks level');
legend('h1', 'h2', 'h3', 'h4');
grid on;

subplot(336)
plot(T,H_noisy_norm, "LineWidth",1)
xlabel('Time (s)');
ylabel('Height (cm)');
title('Normalized height evolution of tanks level');
legend('h1', 'h2', 'h3', 'h4');
grid on;

% Define a step of 50% after 60 seconds
F1 = 300*ones(size(time));
F1(time >= 5*60) = 300 + 50/100 * 300;
F2 = F1;

u = [F1', F2'];

% Interpolate u based on current time (to be used in ode simulation)
u_interp = @(t) interp1(time, u, t, 'linear', 'extrap')';

% Solving the simulation 
ode_fun = @(t, x) ModifiedFourTankSystem(t, x, u_interp(t), d_interp(t), p);
[T, X] = ode15s(ode_fun, [t0 t_f], x0);

% Compute the measured variables
H = X ./ (rho * A);

if(noise_level > 0)
    R = noise_level * 0.5 * eye(4);
    Lr = chol(R,'lower');
    v = Lr*randn(4,n);  
    v = v';
end
H_noisy = H + v(1:length(H), :);

% Normalize the step response
H_noisy_norm = (H_noisy - h_steady') / (50/100 * 300);

subplot(337)
u = [u(:, 1), u(:, 2), d(:, 1), d(:, 2)];
plot(time, u, "LineWidth",1)
xlabel('Time (s)');
ylabel('Input Flow (cm^3/s)');
title('Inputs with 50% increase on F1 and F2');
legend('F1', 'F2', 'F3', 'F4');
grid on;
axis([0, length(time), 0, 500])

subplot(338)
plot(T, H_noisy, "LineWidth",1)
xlabel('Time (s)');
ylabel('Height (cm)');
title('Height evolution of tanks level');
legend('h1', 'h2', 'h3', 'h4');
grid on;

subplot(339)
plot(T,H_noisy_norm, "LineWidth",1)
xlabel('Time (s)');
ylabel('Height (cm)');
title('Normalized height evolution of tanks level');
legend('h1', 'h2', 'h3', 'h4');
grid on;

%% System Identification for G11 and G21

% We identify G11 and G21 which is tf from F1 to H1 and from F1 to H2

% Define a step of 50% after 60 seconds on F1
F1 = 300*ones(size(time));
F2 = 300*ones(size(time));
F1(time >= 5*60) = 450;
F3 = 100*ones(size(time));
F4 = F3;

u = [F1', F2'];
x0 = x_steady;

% Define time vector
time = linspace(t0, t_f, length(u));
n = length(time);

% Interpolate u based on current time (to be used in ode simulation)
u_interp = @(t) interp1(time, u, t, 'linear', 'extrap')';

% Solving the simulation 
ode_fun = @(t, x) ModifiedFourTankSystem(t, x, u_interp(t), d_interp(t), p);
[T, X] = ode15s(ode_fun, [t0 t_f], x0);

% Compute the measured variables
A = p(5:8,1)';
H = X ./ (rho * A);

H_norm = (H - h_steady') / (50/100 * 300);

figure
subplot(121)
u = [u(:, 1), u(:, 2), d(:, 1), d(:, 2)];

plot(time, u, "LineWidth",1)
xlabel('Time (s)');
ylabel('Input Flow (cm^3/s)');
title('Inputs with 50% increase on F1 only');
legend('F1', 'F2', 'F3', 'F4');
grid on;
axis([0, length(time), 0, 500])

subplot(122)
plot(T, H_norm, "LineWidth",1)
xlabel('Time (s)');
ylabel('Height (cm)');
title('Identification results');
grid on;

% Identification of G11
yss = H(length(H), 1);
y0 = H(2, 1);

delta_y = yss - y0;
y_at_time_constant = 0.63*delta_y + y0;

T = 410.5 - 300; % From Graph
K = delta_y/150;

G11 = tf(K, [T 1])

% Identification of G21
yss = H(length(H), 2);
y0 = H(2, 2);

delta_y = yss - y0;
y_at_time_constant = 0.63*delta_y + y0;

T = 480 - 300; % From Graph
K = delta_y/150;

G21 = tf(K, [T 1])

% Check identification results

% Convert transfer functions to state-space
%[G11_A, G11_B, G11_C, G11_D] = tf2ss(G11.num{1}, G11.den{1});
%[G21_A, G21_B, G21_C, G21_D] = tf2ss(G21.num{1}, G21.den{1});

%x0_1 = h_steady(1); % Initial state for Tank 1
%x0_2 = h_steady(2); % Initial state for Tank 2

% Simulate the responses with initial conditions
%H1_simulated = lsim(ss(G11_A, G11_B, G11_C, G11_D), F1, time, x0_1); % Tank 1 response
%H2_simulated = lsim(ss(G21_A, G21_B, G21_C, G21_D), F1, time, x0_2); % Tank 2 response

% Simulate step responses
[G11_simulated] = step(G11, 0:t_f-300);
[G21_simulated] = step(G21, 0:t_f-300);

% Define step response start time
step_start_time = 30;

% Pad responses with initial steady-state value
G11_shifted = [zeros(length(0:0.1:step_start_time - 0.1), 1); G11_simulated];
G21_shifted = [zeros(length(0:0.1:step_start_time - 0.1), 1); G21_simulated];
% Plot results
subplot(122)
hold on
plot(G11_shifted,'-g','LineWidth', 1)
plot(G21_shifted,'-k','LineWidth', 1)
legend('h1', 'h2', 'h3', 'h4', 'Identified G11', 'Identified G21');

%% System Identification for G12 and G22

% We identify G12 and G22 which is tf from F2 to H1 and from F2 to H2

% Define a step of 50% after 60 seconds on F1
F1 = 300*ones(size(time));
F2 = 300*ones(size(time));
F2(time >= 5*60) = 450;
F3 = 100*ones(size(time));
F4 = F3;

u = [F1', F2'];
x0 = x_steady;

% Define time vector
time = linspace(t0, t_f, length(u));
n = length(time);

% Interpolate u based on current time (to be used in ode simulation)
u_interp = @(t) interp1(time, u, t, 'linear', 'extrap')';

% Solving the simulation 
ode_fun = @(t, x) ModifiedFourTankSystem(t, x, u_interp(t), d_interp(t), p);
[T, X] = ode15s(ode_fun, [t0 t_f], x0);

% Compute the measured variables
A = p(5:8,1)';
H = X ./ (rho * A);

H_norm = (H - h_steady') / (50/100 * 300);

figure

subplot(121)
u = [u(:, 1), u(:, 2), d(:, 1), d(:, 2)];

plot(time, u, "LineWidth",1)
xlabel('Time (s)');
ylabel('Input Flow (cm^3/s)');
title('Inputs with 50% increase on F2 only');
legend('F1', 'F2', 'F3', 'F4');
grid on;
axis([0, length(time), 0, 500])

subplot(122)
plot(T, H_norm, "LineWidth",1)
xlabel('Time (s)');
ylabel('Height (cm)');
title('Identification results');
grid on;

% Identification of G11
yss = H(length(H), 1);
y0 = H(2, 1);

delta_y = yss - y0;
y_at_time_constant = 0.63*delta_y + y0;

T = 467 - 300; % From Graph
K = delta_y/150;

G12 = tf(K, [T 1])

% Identification of G22
yss = H(length(H), 2);
y0 = H(2, 2);

delta_y = yss - y0;
y_at_time_constant = 0.63*delta_y + y0;

T = 430 - 300; % From Graph
K = delta_y/150;

G22 = tf(K, [T 1])

% Check identification results
% Define step response start time
step_start_time = 30;

[G12_simulated] = step(G12, 0:t_f-300);
[G22_simulated] = step(G22, 0:t_f-300);

% Pad responses with initial steady-state value
G12_shifted = [zeros(length(0:0.1:step_start_time - 0.1), 1); G12_simulated];
G22_shifted = [zeros(length(0:0.1:step_start_time - 0.1), 1); G22_simulated];

% Plot results
subplot(122)
hold on
plot(G12_shifted,'-g','LineWidth', 1)
plot(G22_shifted,'-k','LineWidth', 1)
legend('h1', 'h2', 'h3', 'h4', 'Identified G12', 'Identified G22');

%% Markov Parameters

% Sampling time
Ts = 5; 

% Discretize the transfer functions
G11_d = c2d(G11, Ts, 'zoh'); % Zero-order hold discretization
G21_d = c2d(G21, Ts, 'zoh');
G12_d = c2d(G12, Ts, 'zoh');
G22_d = c2d(G22, Ts, 'zoh');

% Compute impulse response coefficients (Markov parameters)
n_impulse = 500; % Number of Markov parameters to compute
[impulse_G11, t_G11] = impulse(G11_d, n_impulse);
[impulse_G21, t_G21] = impulse(G21_d, n_impulse);
[impulse_G12, t_G12] = impulse(G12_d, n_impulse);
[impulse_G22, t_G22] = impulse(G22_d, n_impulse);

% Plot the impulse response coefficients
figure;

% G11
subplot(2, 2, 1);
stem(t_G11, squeeze(impulse_G11), 'filled', 'LineWidth', 1.5);
xlabel('Time (samples)');
ylabel('Amplitude');
title('Impulse Response Coefficients of G_{11}');
grid on;

% G21
subplot(2, 2, 2);
stem(t_G21, squeeze(impulse_G21), 'filled', 'LineWidth', 1.5);
xlabel('Time (samples)');
ylabel('Amplitude');
title('Impulse Response Coefficients of G_{21}');
grid on;

% G12
subplot(2, 2, 3);
stem(t_G12, squeeze(impulse_G12), 'filled', 'LineWidth', 1.5);
xlabel('Time (samples)');
ylabel('Amplitude');
title('Impulse Response Coefficients of G_{12}');
grid on;

% G22
subplot(2, 2, 4);
stem(t_G22, squeeze(impulse_G22), 'filled', 'LineWidth', 1.5);
xlabel('Time (samples)');
ylabel('Amplitude');
title('Impulse Response Coefficients of G_{22}');
grid on;

%%
% Discretize the Transfer Functions
Ts = 25; % Sampling time
num_samples = 100; % Number of samples

% Discretize G11, G21, G12, G22
G11_d = c2d(G11, Ts, 'impulse');
G21_d = c2d(G21, Ts, 'impulse');
G12_d = c2d(G12, Ts, 'impulse');
G22_d = c2d(G22, Ts, 'impulse');

% State space for G11, G21, G12, G22
[G11_A, G11_B, G11_C, G11_D] = tf2ss(G11.num{1}, G11.den{1});
[G21_A, G21_B, G21_C, G21_D] = tf2ss(G21.num{1}, G21.den{1});
[G12_A, G12_B, G12_C, G12_D] = tf2ss(G12.num{1}, G12.den{1});
[G22_A, G22_B, G22_C, G22_D] = tf2ss(G22.num{1}, G22.den{1});

% Compute Markov parameters (impulse response coefficients for discrete tf)
impulse_G11 = impulse(G11_d, 0:Ts:(num_samples-1)*Ts);
impulse_G21 = impulse(G21_d, 0:Ts:(num_samples-1)*Ts);
impulse_G12 = impulse(G12_d, 0:Ts:(num_samples-1)*Ts);
impulse_G22 = impulse(G22_d, 0:Ts:(num_samples-1)*Ts);

% Impulse response for the state-space models
[impulse_ss_G11, ~] = impulse(ss(G11_A, G11_B, G11_C, G11_D), 0:Ts:(num_samples-1)*Ts);
[impulse_ss_G21, ~] = impulse(ss(G21_A, G21_B, G21_C, G21_D), 0:Ts:(num_samples-1)*Ts);
[impulse_ss_G12, ~] = impulse(ss(G12_A, G12_B, G12_C, G12_D), 0:Ts:(num_samples-1)*Ts);
[impulse_ss_G22, ~] = impulse(ss(G22_A, G22_B, G22_C, G22_D), 0:Ts:(num_samples-1)*Ts);

figure
subplot(2, 2, 1);
stem(0:num_samples-1, impulse_G11, 'b', 'filled'); hold on;
plot(0:num_samples-1, impulse_ss_G11, 'r', 'LineWidth', 1.5);
title('Impulse Response Coefficients of G_{11}');
xlabel('Time (samples)'); ylabel('Amplitude');
legend('TF', 'State Space', 'Location', 'Best');
grid on;


subplot(2, 2, 2);
stem(0:num_samples-1, impulse_G21, 'b', 'filled'); hold on;
plot(0:num_samples-1, impulse_ss_G21, 'r', 'LineWidth', 1.5);
title('Impulse Response Coefficients of G_{21}');
xlabel('Time (samples)'); ylabel('Amplitude');
legend('TF', 'State Space', 'Location', 'Best');
grid on;

subplot(2, 2, 3);
stem(0:num_samples-1, impulse_G12, 'b', 'filled'); hold on;
plot(0:num_samples-1, impulse_ss_G12, 'r', 'LineWidth', 1.5);
title('Impulse Response Coefficients of G_{12}');
xlabel('Time (samples)'); ylabel('Amplitude');
legend('TF', 'State Space', 'Location', 'Best');
grid on;

subplot(2, 2, 4);
stem(0:num_samples-1, impulse_G22, 'b', 'filled'); hold on;
plot(0:num_samples-1, impulse_ss_G22, 'r', 'LineWidth', 1.5);
title('Impulse Response Coefficients of G_{22}');
xlabel('Time (samples)'); ylabel('Amplitude');
legend('TF Markov', 'SS Impulse', 'Location', 'Best');
grid on;