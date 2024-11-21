%% Simulation scenario
clear variables
clc

% Parameters
a1 = 1.2272;     %[cm2] Area of outlet pipe 1
a2 = 1.2272;     %[cm2] Area of outlet pipe 2
a3 = 1.2272;     %[cm2] Area of outlet pipe 3
a4 = 1.2272;     %[cm2] Area of outlet pipe 4

A1 = 380.1327;   %[cm2] Cross sectional area of tank 1
A2 = 380.1327;   %[cm2] Cross sectional area of tank 2
A3 = 380.1327;   %[cm2] Cross sectional area of tank 3
A4 = 380.1327;   %[cm2] Cross sectional area of tank 4

g = 981;         %[cm/s2] The acceleration of gravity

gamma1 = 0.45;  % Flow distribution constant. Valve 1
gamma2 = 0.40;  % Flow distribution constant. Valve 2

rho = 1; % Density [g/cm3]

p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];

h1_0 = 0.0; % [cm] Liquid level in tank 1 at time t0
h2_0 = 0.0; % [cm] Liquid level in tank 2 at time t0
h3_0 = 0.0; % [cm] Liquid level in tank 3 at time t0
h4_0 = 0.0; % [cm] Liquid level in tank 4 at time t0

% Computing the masses at time 0:
m1_0 = rho * A1 * h1_0;
m2_0 = rho * A2 * h2_0;
m3_0 = rho * A3 * h3_0;
m4_0 = rho * A4 * h4_0;

% Initial state vector
x0 = [m1_0; m2_0; m3_0; m4_0];

% Simulation time
t0 = 0.0; % [s] Initial time
tf = 20 * 60; % [s] Final time
t = [t0 tf];

% Pump flow rates
F1 = 300 * ones(1, tf); % [cm3/s] Flow rate from pump 1
F2 = 300 * ones(1, tf); % [cm3/s] Flow rate from pump 2
u = [F1; F2];


% Run simulation using ode15s
[t, x] = ode15s(@(t, x) QuadrupleTankProcess(t, x, u, p, tf), t, x0);

% We now have the time at which it solved the equations t and
% the solutions for each time point in each row of x

% Extract masses from the result
m1 = x(:, 1);
m2 = x(:, 2);
m3 = x(:, 3);
m4 = x(:, 4);

% Compute liquid levels (h) from masses (m)
h1 = m1 / (rho * A1);
h2 = m2 / (rho * A2);
h3 = m3 / (rho * A3);
h4 = m4 / (rho * A4);

% Plot results
figure;
plot(t / 60, h1, 'r-', 'LineWidth', 2); hold on;
plot(t / 60, h2, 'g-', 'LineWidth', 2);
plot(t / 60, h3, 'b-', 'LineWidth', 2);
plot(t / 60, h4, 'k-', 'LineWidth', 2);

xlabel('Time [min]');
ylabel('Liquid height [cm]');

legend('Tank 1', 'Tank 2', 'Tank 3', 'Tank 4');
title('Quadruple Tank System Simulation');
grid on; hold off;

%% Now simulate for 5, 10, 25% step increase in F1

% Simulation time
t0 = 0.0; % [s] Initial time
tf = 60 * 60; % [s] Final time 
t = t0:1:tf;

% Pump flow rates
F1 = 300 * ones(1, tf+1); % [cm3/s] Flow rate from pump 1
F2 = 300 * ones(1, tf+1); % [cm3/s] Flow rate from pump 2

% Custom F1
for i = 1:(20*60)
    F1(i) = F1(i) + 0.05 * F1(i);
end
for i = (20*60+1):(40*60)
    F1(i) = F1(i) + 0.1 * F1(i);
end
for i = (40*60+1):(60*60)
    F1(i) = F1(i) + 0.25 * F1(i);
end

u = [F1; F2];

% Run simulation using ode15s
[t, x] = ode15s(@(t, x) QuadrupleTankProcess(t, x, u, p, tf), t, x0);

% We now have the time at which it solved the equations t and
% the solutions for each time point in each row of x

% Extract masses from the result
m1 = x(:, 1);
m2 = x(:, 2);
m3 = x(:, 3);
m4 = x(:, 4);

% Compute liquid levels (h) from masses (m)
h1 = m1 / (rho * A1);
h2 = m2 / (rho * A2);
h3 = m3 / (rho * A3);
h4 = m4 / (rho * A4);

% Plot results
figure;
subplot(121)

plot(t/60, F1, 'r-', 'LineWidth', 2); hold on
plot(t/60, F2, 'g-', 'LineWidth', 2); hold on

xlabel('Time [min]');
ylabel('Input flows [cm3/s]');
legend('F1', 'F2');
axis([0, tf/60, 250, 400]); grid;
subplot(122)

plot(t / 60, h1, 'r-', 'LineWidth', 2); hold on;
plot(t / 60, h2, 'g-', 'LineWidth', 2);
plot(t / 60, h3, 'b-', 'LineWidth', 2);
plot(t / 60, h4, 'k-', 'LineWidth', 2);

xlabel('Time [min]');
ylabel('Liquid height [cm]');

legend('Tank 1', 'Tank 2', 'Tank 3', 'Tank 4');
title('Quadruple Tank System Simulation with variable F1');
grid on; hold off;