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
t_f = 20*60; % [s] Final time
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

%% Simulating Steps on Deterministic Model

% Simulating a step response for F1 = F2 = 300 cm3/s

% Define a step of 10% after 60 seconds
F1 = 300*ones(size(time));
F1(time >= 5*60) = 330;
F2 = F1;
F3 = 100*ones(size(time));
F4 = F3;

u = [F1', F2'];
x0 = x_steady;

% Simulating
d = [F3', F4'];

% Define time vector
time = linspace(t0, t_f, length(u));

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

figure("Name", "Deterministic Simulations")

subplot(321)
u = [u(:, 1), u(:, 2), d(:, 1), d(:, 2)];
plot(time, u, "LineWidth",1)
xlabel('Time (s)');
ylabel('Input Flow (cm^3/s)');
axis([-50, 1250, -50, 350])
title('Inputs with 10% increase on F1 and F2');
legend('F1', 'F2', 'F3', 'F4');
grid on;

subplot(322)
plot(T, H, "LineWidth",1)
xlabel('Time (s)');
ylabel('Height (cm)');
title('Height evolution of tanks level');
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

subplot(323)
u = [u(:, 1), u(:, 2), d(:, 1), d(:, 2)];
plot(time, u, "LineWidth",1)
xlabel('Time (s)');
ylabel('Input Flow (cm^3/s)');
title('Inputs with 25% increase on F1 and F2');
legend('F1', 'F2', 'F3', 'F4');
grid on;

subplot(324)
plot(T, H, "LineWidth",1)
xlabel('Time (s)');
ylabel('Height (cm)');
title('Height evolution of tanks level');
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

subplot(325)
u = [u(:, 1), u(:, 2), d(:, 1), d(:, 2)];
plot(time, u, "LineWidth",1)
xlabel('Time (s)');
ylabel('Input Flow (cm^3/s)');
title('Inputs with 50% increase on F1 and F2');
legend('F1', 'F2', 'F3', 'F4');
grid on;

subplot(326)
plot(T, H, "LineWidth",1)
xlabel('Time (s)');
ylabel('Height (cm)');
title('Height evolution of tanks level');
legend('h1', 'h2', 'h3', 'h4');
grid on;