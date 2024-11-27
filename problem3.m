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

%% Obtaining the steady state values

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

% Display results
disp('Steady-state liquid masses (g):');
disp(x_steady);

disp('Steady-state liquid levels (cm):');
disp(h_steady);

%% Computing the transfer functions

% Simulating a step response for F1 = F2 = 300 cm3/s

t0 = 0.0;   % [s] Initial time
t_f = 20*60; % [s] Final time

F1 = 330;   % [cm3/s] Flow rate from pump 1
F2 = 330;

F3 = 0;
F4 = 0;

x0 = x_steady;

u = [F1; F2];

% Simulating
d = [F3; F4];

% Run simulation using ode15s
[T, X] = ode15s(@(t, x) ModifiedFourTankSystem(t, x, u, d, p), [t0 t_f], x0);

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

figure("Name", "System Identification")

subplot(131)
u = [F1 * ones(nT, 1), F2 * ones(nT, 1), zeros(nT, 1), zeros(nT, 1)];
plot(T, u, "LineWidth",1)
xlabel('Time (s)');
ylabel('Input Flow (cm^3/s)');
axis([-50, 1250, -50, 350])
title('Inputs');
legend('F1', 'F2', 'F3', 'F4');
grid on;

subplot(132)
plot(T, H, "LineWidth",1)
xlabel('Time (s)');
ylabel('Height (cm)');
title('Height evolution of tanks level');
axis([-50, 600, 20, 50])
legend('h1', 'h2');
grid on;

% We read from the graph:
h1_ss = 24.67;
h1_f = 29.85;
time_constant_value_1 = 0.632 * (h1_f - h1_ss) + h1_ss;

h2_ss = 36.85;
h2_f = 44.59;
time_constant_value_2 = 0.632 * (h2_f - h2_ss) + h2_ss;

T_1 = 80; % (s) read from graph at time_constant_value_1
K_1 = (h1_f - h1_ss)/(330 - 300)-0.082;
H1 = K_1 * tf(1, [T_1 1]);

T_2 = 100; % (s) read from graph at time_constant_value_2
K_2 = (h2_f - h2_ss)/(330 - 300)-0.123;
H2 = K_2 * tf(1, [T_2 1]);

% Check identification results

% Define the input signal (step change from 300 to 330)
t = 0 : 20*60; % Time vector
F1 = 330 * ones(size(t));
F2 = F1;

% Convert transfer functions to state-space
[H1_A, H1_B, H1_C, H1_D] = tf2ss(H1.num{1}, H1.den{1});
[H2_A, H2_B, H2_C, H2_D] = tf2ss(H2.num{1}, H2.den{1});

% Steady-state heights
h0 = [24.6, 36.8]; % Initial heights for Tank 1 and Tank 2

% Calculate initial states corresponding to h0
x0_1 = (h0(1) - H1_D * F1(1)) / H1_C; % Initial state for Tank 1
x0_2 = (h0(2) - H2_D * F2(1)) / H2_C; % Initial state for Tank 2

% Simulate the responses with initial conditions
H1_simulated = lsim(ss(H1_A, H1_B, H1_C, H1_D), F1, t, x0_1); % Tank 1 response
H2_simulated = lsim(ss(H2_A, H2_B, H2_C, H2_D), F2, t, x0_2); % Tank 2 response

% Plot results
subplot(133)

% Plot Nonlinear system responses (h1 and h2)
plot(T, H(:, 1), 'b', 'LineWidth', 1); 
hold on;
plot(T, H(:, 2), 'r', 'LineWidth', 1);  

% Plot Identified system responses (h1 and h2)
plot(t, H1_simulated, 'g','LineWidth', 1);
plot(t, H2_simulated, 'k','LineWidth', 1);

% Add labels, title, legend, and grid
xlabel('Time (s)');
ylabel('Tank Height (cm)');
title('Identification results');
legend('h1 nonlinear system', 'h2 nonlinear system', 'h1 identified system', 'h2 identified system', 'Location', 'Best');
grid on;

%%

% pidTuner was used to obtain the following:

%%%%% For H2:

% P controller
C_P_2 = 79.8383;

% PI controller
Kp_PI_2 = 130.1929; 
Ki_2 = 3.1435;  
C_PI_2 = Kp_PI_2 + tf(Ki_2, [1 0]);  % PI controller transfer function

% PID controller
Kp_PID_2 = 74.415; 
Ki_PID_2 = 6.0065; 
Kd_2 = 62.6315; 
C_PID_2 = Kp_PID_2 + tf(Ki_PID_2, [1 0]) + tf([Kd_2 0], 1);  % PID controller transfer function

% Closed-loop systems for P, PI, and PID controllers
H2_P = feedback(C_P_2 * H2, 1);  
H2_PI = feedback(C_PI_2 * H2, 1);  
H2_PID = feedback(C_PID_2 * H2, 1);  

figure;
subplot(2,3,4);
step(H2_P);
title('H2 Step Response with P Controller');

subplot(2,3,5);
step(H2_PI);
title('H2 Step Response with PI Controller');

subplot(2,3,6);
step(H2_PID);
title('H2 Step Response with PID Controller');

%%%%% For H1:
% Proportional controller (P)
C_P = 124.4322;

% PI controller
Kp_PI = 73.9979; 
Ki = 2.3444;  
C_PI = Kp_PI + tf(Ki, [1 0]);  % PI controller transfer function

% PID controller
Kp_PID = 73.9979; 
Ki_PID = 2.3444; 
Kd = 0; 
C_PID = Kp_PID + tf(Ki_PID, [1 0]) + tf([Kd 0], 1);  % PID controller transfer function

% Closed-loop systems for P, PI, and PID controllers
H1_P = feedback(C_P * H1, 1);  
H1_PI = feedback(C_PI * H1, 1);  
H1_PID = feedback(C_PID * H1, 1);  

% Step response for each system

subplot(2,3,1);
step(H1_P);
title('H1 Step Response with P Controller');

subplot(2,3,2);
step(H1_PI);
title('H1 Step Response with PI Controller');

subplot(2,3,3);
step(H1_PID);
title('H1 Step Response with PID Controller');

%% Closed-Loop Simulations on the nonlinear model

% Reference tank levels [cm]
h1_ref = h_steady(1) + 25/100 * h_steady(1); % 5% increase from steady state
h2_ref = h_steady(2) + 25/100 * h_steady(2);

% Basically these heigths are steady state ones

% PID Controller Gains
Kp_PID = [Kp_PID, Kp_PID_2];  
Ki_PID = [Ki_PID, Ki_PID_2]; 
Kd_PID = [Kd, Kd_2];

% P Controller:
Kp_p = [C_P, C_P_2]; 
Ki_p = [0, 0];
Kd_p = [0, 0];

% PI Controller
Kp_PI = [Kp_PI, Kp_PI_2];  
Ki_PI = [Ki, Ki_2]; 
Kd_PI = [0, 0];

%% Simulate the Closed-Loop System
t_f = 20 * 60; % Final simulation time [s]
dt = 0.1;      % Time step [s]
T = 0:dt:t_f;  % Time vector
X = zeros(length(T), length(x0));
X(1, :) = x0; % Initial conditions are considered the steady state values

% Initialize variables
e_prev = [0; 0]; % Initial error
e_int = [0; 0];  % Initial integral term
e_prev_time = 0; % Initial time for error computation

% Time loop
for i = 2:length(T)
    t = T(i);

    % Call the dynamics function with PID Controllers Parameters
    [dx, F1, F2, e_prev, e_int, e_prev_time] = closedLoopDynamics(t, X(i-1, :)', u, p, [h1_ref; h2_ref], Kp_PID, Ki_PID, Kd_PID, e_prev, e_int, e_prev_time);

    % Update state using Euler's method
    X(i, :) = X(i-1, :) + (dx' * dt);
end

% Extract and Plot Results
H = X ./ (rho * A);

figure("Name", "Closed-Loop Simulation")

subplot(311)
plot(T, H(:,1:2), "LineWidth", 1);
xlabel('Time (s)');
ylabel('Height (cm)');
title('Height evolution of tank levels with PID Controller');
legend('h1', 'h2');
grid on;

% Now simulate only for PI Controller
X = zeros(length(T), length(x0));
X(1, :) = x0; % Initial conditions are considered the steady state values

% Initialize variables
e_prev = [0; 0]; % Initial error
e_int = [0; 0];  % Initial integral term
e_prev_time = 0; % Initial time for error computation

% Time loop
for i = 2:length(T)
    t = T(i);

    % Call the dynamics function with PI Controllers Parameters
    [dx, F1, F2, e_prev, e_int, e_prev_time] = closedLoopDynamics(t, X(i-1, :)', u, p, [h1_ref; h2_ref], Kp_PI, Ki_PI, Kd_PI, e_prev, e_int, e_prev_time);

    % Update state using Euler's method
    X(i, :) = X(i-1, :) + (dx' * dt);
end

% Extract and Plot Results
H = X ./ (rho * A);

subplot(312)
plot(T, H(:,1:2), "LineWidth", 1);
xlabel('Time (s)');
ylabel('Height (cm)');
title('Height evolution of tank levels with PI Controller');
axis([0, 1200, 20, 50])
legend('h1', 'h2');
grid on;

% Now Simulate with only P controllers

X = zeros(length(T), length(x0));
X(1, :) = x0; % Initial conditions are considered the steady state values

% Initialize variables
e_prev = [0; 0]; % Initial error
e_int = [0; 0];  % Initial integral term
e_prev_time = 0; % Initial time for error computation

% Time loop
for i = 2:length(T)
    t = T(i);

    % Call the dynamics function with P Controllers Parameters
    [dx, F1, F2, e_prev, e_int, e_prev_time] = closedLoopDynamics(t, X(i-1, :)', u, p, [h1_ref; h2_ref], Kp_p, Ki_p, Kd_p, e_prev, e_int, e_prev_time);

    % Update state using Euler's method
    X(i, :) = X(i-1, :) + (dx' * dt);
end

% Extract and Plot Results
H = X ./ (rho * A);

subplot(313)
plot(T, H(:,1:2), "LineWidth", 1);
xlabel('Time (s)');
ylabel('Height (cm)');
title('Height evolution of tank levels with P Controller');
axis([0, 1200, 20, 50])
legend('h1', 'h2');
grid on;
