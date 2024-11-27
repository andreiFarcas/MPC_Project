%% 
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

%% 2.1 Deterministic Nonlinear Scenario

% Scenario Definition

t0 = 0.0;   % [s] Initial time
tf = 20*60; % [s] Final time
m10 = 0.0;  % [g] Liquid mass in tank 1 at time t0
m20 = 0.0;
m30 = 0.0;
m40 = 0.0;
F1 = 200;   % [cm3/s] Flow rate from pump 1
F2 = 200;

F3 = 0;
F4 = 0;

x0 = [m10; m20; m30; m40];
u = [F1; F2; F3; F4];

% Simulating
d = [];

% Run simulation using ode15s
[T, X] = ode15s(@(t, x) ModifiedFourTankSystem(t, x, u, d, p), [t0 tf], x0);

% Extracting parameters and state size
[nT,nX] = size(X);
a = p(1:4,1)';
A = p(5:8,1)';

% Compute the measured variables
H = zeros(nT,nX);
for i=1:nT
    H(i,:) = X(i,:)./(rho*A);
end  

% Compute the flows out of each tank
Qout = zeros(nT,nX);
for i=1:nT
    Qout(i,:) = a.*sqrt(2*g*H(i,:));
end

figure("Name","Deterministic Nonlinear Model")

subplot(221)
u = [F1 * ones(nT, 1), F2 * ones(nT, 1), zeros(nT, 1), zeros(nT, 1)];
plot(T, u, "LineWidth",1)
xlabel('Time (s)');
ylabel('Input Flow (cm^3/s)');
axis([-50, 1250, -50, 250])
title('Inputs');
legend('F1', 'F2', 'F3', 'F4');
grid on;

subplot(222)
plot(T, H, "LineWidth",1)
xlabel('Time (s)');
ylabel('Height (cm)');
title('Height evolution of tanks level');
legend('h1', 'h2', 'h3', 'h4');
grid on;

subplot(223)
plot(T, X/1000, "LineWidth",1)
xlabel('Time (s)');
ylabel('Mass (kg)');
title('Mass evolution of tanks level');
legend('x1', 'x2', 'x3', 'x4');
grid on;

subplot(224)
plot(T, Qout, "LineWidth",1)
xlabel('Time (s)');
ylabel('Flow (cm^3/s)');
title('Flow out of tanks level');
legend('q1', 'q2', 'q3', 'q4');
grid on;

%% 2.2 Stochastic Nonlinear Model Scenario

%-------------------------------------------------------
% Simulation scenario
%-------------------------------------------------------
t0 = 0.0;   % [s] Initial time
tf = 20*60; % [s] Final time
m10 = 0.0;  % [g] Liquid mass in tank 1 at time t0
m20 = 0.0;
m30 = 0.0;
m40 = 0.0;
F1 = 300;   % [cm3/s] Flow rate from pump 1
F2 = 250;
F3 = 100;   % [cm3/s] Baseline flow rate disturbance
F4 = 100;

x0 = [m10; m20; m30; m40];

%-------------------------------------------------------
% Deterministic disturbances (F3 and F4 as piecewise constant)
%-------------------------------------------------------
Ts = 10;                    % [s] Sample Time
t = t0:Ts:tf;               % [s] Sample instants
n = length(t);              % Number of time steps
k = 6;                      % Number of samples per disturbance value

% Control inputs (F1 and F2 only)
% Control inputs (F1 and F2 only)
n1 = 60; % Number of time steps for the first half

% Create each column separately
u1 = [0.5 * F1 * ones(n1, 1); F1 * ones(n - n1, 1)];
u2 = [0.5 * F2 * ones(n1, 1); F2 * ones(n - n1, 1)];

% Concatenate the columns to form the input matrix
u = [u1, u2];

% Generate piecewise constant disturbances for F3 and F4
F3_disturbance = randi([-50, 50], ceil((n+1)/k), 1); % Disturbance for F3
F4_disturbance = randi([-50, 50], ceil((n+1)/k), 1); % Disturbance for F4
F3_disturbance = repelem(F3_disturbance, k, 1);
F4_disturbance = repelem(F4_disturbance, k, 1);

% Truncate to match time steps
F3_disturbance = F3_disturbance(1:n);
F4_disturbance = F4_disturbance(1:n);

% Combine baseline and disturbances
d = [F3 + F3_disturbance'; F4 + F4_disturbance'];  % Disturbances for F3 and F4
d = max(0, d);                                   % Ensure non-negative flows

%-------------------------------------------------------
% Plotting the disturbances
%-------------------------------------------------------
figure()
stairs(t, d')
ylim([0, 400])
legend("F3", "F4")
title("Piecewise constant stochastic disturbances on F3 and F4")
xlabel("Time [s]")
ylabel("Flow rate [cm^3/s]")

%-------------------------------------------------------
% Performing simulations
%-------------------------------------------------------
d = d';

% Define time vector
time = linspace(t0, tf, length(u));

% Interpolate u and d based on current time (to be used in ode simulation)
u_interp = @(t) interp1(time, u, t, 'linear', 'extrap');
d_interp = @(t) interp1(time, d, t, 'linear', 'extrap');

% Solving the simulation 
ode_fun = @(t, x) ModifiedFourTankSystem(t, x, u_interp(t), d_interp(t), p);
[T, X] = ode15s(ode_fun, [t0 tf], x0);

% Extracting parameters and state size
[nT, nX] = size(X);
a = p(1:4,1)';
A = p(5:8,1)';

% Compute the measured variables
H = X ./ (rho * A);

% Add the measurement noise
R = eye(4);
Lr = chol(R,'lower');
v = Lr*randn(4,nT);
v = v';
figure()
plot(v)
title("Measurement noise")
legend("Sensor for h1", "Sensor for h2", "Sensor for h3", "Sensor for h4")

H_noisy = H + v(1:nT, :);

% Compute the flows out of each tank
Qout = a .* sqrt(2 * g * H);

figure("Name","Stochastic Nonlinear Model")

subplot(221)
plot(t, [u, d], "LineWidth",1)
xlabel('Time (s)');
ylabel('Input Flow (cm^3/s)');
axis([-50, 1250, -50, 250])
title('Inputs');
legend('F1', 'F2', 'F3', 'F4');
grid on;

subplot(222)
plot(T, H_noisy)
xlabel('Time (s)');
ylabel('Height (cm)');
title('Height evolution of tanks level, with noise applied');
legend('h1', 'h2', 'h3', 'h4');
grid on;

subplot(223)
plot(T, X/1000)
xlabel('Time (s)');
ylabel('Mass (kg)');
title('Mass evolution of tanks level');
legend('x1', 'x2', 'x3', 'x4');
grid on;

subplot(224)
plot(T, Qout)
xlabel('Time (s)');
ylabel('Flow (cm^3/s)');
title('Flow out of tanks level');
legend('q1', 'q2', 'q3', 'q4');
grid on;

%% 2.3. Stochastic Nonlinear Model (SDE)

%-------------------------------------------------------
% Simulation scenario
%-------------------------------------------------------
t0 = 0.0;   % [s] Initial time
tf = 20*60; % [s] Final time
m10 = 0.0;  % [g] Liquid mass in tank 1 at time t0
m20 = 0.0;
m30 = 0.0;
m40 = 0.0;
F1 = 250;   % [cm3/s] Flow rate from pump 1
F2 = 200;
F3 = 100;   % [cm3/s] Baseline flow rate disturbance
F4 = 100;

x0 = [m10; m20; m30; m40];

% Control inputs (F1 and F2 only)
n1 = 60; % Number of time steps for the first half

% Create each column separately
u1 = [0.5 * F1 * ones(n1, 1); F1 * ones(n - n1, 1)];
u2 = [0.5 * F2 * ones(n1, 1); F2 * ones(n - n1, 1)];

% Concatenate the columns to form the input matrix
u = [u1, u2];

%-------------------------------------------------------
% Process Noise (Brownian Motion for F3 and F4)
%-------------------------------------------------------
Q = [20^2 0; 0 40^2];       % Covariance of the Brownian motion (units^2/s)
Lq = chol(Q, 'lower');       % Cholesky decomposition for sampling
dt = Ts;                     % Time step size
dw = Lq * randn(2, n);       % Generate Brownian increments (scaled noise)
dw = dw / 10;              % Scale down for numerical stability
w = cumsum(dw, 2);           % Cumulative sum to compute Wiener process
w = w';                      % Transpose for time along rows, disturbances in columns

% Add baseline flow rates
F3_stochastic = F3 + w(:, 1); 
F4_stochastic = F4 + w(:, 2); 

% Ensure non-negative flows
F3_stochastic = max(0, F3_stochastic);
F4_stochastic = max(0, F4_stochastic);

% Combine into disturbances matrix
d = [F3_stochastic, F4_stochastic];

% Plotting the disturbances
figure()
plot(w(:,1), 'DisplayName', 'dF3')
hold on
plot(w(:,2), 'DisplayName', 'dF4')
legend("dF3", "dF4")
title("Brownian Motion for Stochastic Process Noise")
xlabel("Time Steps")
ylabel("Magnitude")

% Simulation with stochastic disturbances


% Define time vector
time = linspace(t0, tf, length(u));

% Interpolate u and d based on current time (to be used in ode simulation)
u_interp = @(t) interp1(time, u, t, 'linear', 'extrap');
d_interp = @(t) interp1(time, d, t, 'linear', 'extrap');

% Solving the simulation 
ode_fun = @(t, x) ModifiedFourTankSystem(t, x, u_interp(t), d_interp(t), p);
[T, X] = ode15s(ode_fun, [t0 tf], x0);

[nT,nX] = size(X);
a = p(1:4,1)';
A = p(5:8,1)';

% Compute the measured variables and add measurement noise
H = zeros(nT,nX);
for i=1:nT
    H(i,:) = X(i,:)./(rho*A);
end  

R = eye(4);
Lr = chol(R,'lower');
v = Lr*randn(4,nT);
v = v';

H_noisy = H + v(1:length(T), :);

% Compute the flows out of each tank
Qout = zeros(nT,nX);
for i=1:nT
    Qout(i,:) = a.*sqrt(2*g*H(i,:));
end

figure()
subplot(221)
plot(t, [u, d], "LineWidth",1)
xlabel('Time (s)');
ylabel('Input Flow (cm^3/s)');
axis([-50, 1250, -50, 250])
title('Inputs');
legend('F1', 'F2', 'F3', 'F4');
grid on;

subplot(222)
plot(T, H_noisy)
xlabel('Time (s)');
ylabel('Height (cm)');
title('Height evolution of tanks level, with noise applied');
legend('h1', 'h2', 'h3', 'h4');
grid on;

subplot(223)
plot(T, X/1000)
xlabel('Time (s)');
ylabel('Mass (kg)');
title('Mass evolution of tanks level');
legend('x1', 'x2', 'x3', 'x4');
grid on;

subplot(224)
plot(T, Qout)
xlabel('Time (s)');
ylabel('Flow (cm^3/s)');
title('Flow out of tanks level');
legend('q1', 'q2', 'q3', 'q4');
grid on;
