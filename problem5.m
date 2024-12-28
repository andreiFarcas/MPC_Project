%% 
clc;
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

F1 = 300; F2 = 300;
F3 = 100; F4 = 100;

m10 = 5000;
m20 = 5000;
m30 = 5000;
m40 = 5000;

t0 = 0.0;   % [s] Initial time
tend = 20*60; % [s] Final time
Ts = 10;                    % [s] Sample Time
t = t0:Ts:tend;               % [s] Sample instants
n = length(t);           % Number of time steps
k = 6;
n_markov = 100;

%% Deterministic model
xs0 = [m10; m20; m30; m40]; % [g] Initial guess on xs
us = [F1; F2]; % [cm3/s] Flow rates
d = [F3;F4]; % model disturbances
v=[0;0;0;0]; % measurement noise

[xs0,ys0,zs0,A,B,Bv,C,D]  = StateSpaceModeling(p,xs0,us,d,v);
xs0
ys0
[taus0, gains0,Ts0,Ff,Gg,Gv,Cc,Dd] = DiscreteModeling(A,B,Bv,C,D,F1,F2,'Deterministic Model');
taus0
gains0
Ts0

markov_params0 = MarkovParam(Ff,Gg,Cc,Dd,n_markov,'Deterministic Model');

%% Piecewise stochastic disturbances

% Generate piecewise constant disturbances for F3 and F4
F3_disturbance = randi([-50, 50], ceil((n+1)/k), 1); % Disturbance for F3
F4_disturbance = randi([-50, 50], ceil((n+1)/k), 1); % Disturbance for F4
F3_disturbance = repelem(F3_disturbance, k, 1);
F4_disturbance = repelem(F4_disturbance, k, 1);
F3_disturbance = F3_disturbance(1:n); 
F4_disturbance = F4_disturbance(1:n);

% Combine baseline and disturbances
d = [F3 + F3_disturbance'; F4 + F4_disturbance'];  % Disturbances for F3 and F4
d = max(0, d);                                   % Ensure non-negative flows

% NOISE MEASUREMENTS
R = eye(4);
Lr = chol(R,'lower');
v = Lr*randn(4,n);
v = v';

xs0 = [m10; m20; m30; m40]; % [g] Initial guess on xs
us = [F1; F2]; % [cm3/s] Flow rates
ds =  mean(d,2); % Average model disturbances noise
vs = mean(v,1)'; % Measurement noise

[xs1,ys1,zs,A,B,Bv,C,D]  = StateSpaceModeling(p,xs0,us,ds,vs);
xs1
ys1
[taus1, gains1,Ts1,Ff,Gg,Gv,Cc,Dd] = DiscreteModeling(A,B,Bv,C,D,F1,F2,'Stochastic1 Model');
taus1
gains1
Ts1


markov_params1 = MarkovParam(Ff,Gg,Cc,Dd,n_markov,'Stochastic1 Model');

%% Brownian motion

% F3,F4 DISTURBANCES
Q = [2^2 0; 0 4^2];          % Covariance of the Brownian motion (units^2/s)
Lq = chol(Q, 'lower');       % Cholesky decomposition for sampling
dt = Ts;                     % Time step size
dw = Lq * randn(2, n);       % Generate Brownian increments (scaled noise)
w = cumsum(dw, 2);           % Cumulative sum to compute Wiener process
w = w';                      % Transpose for time along rows, disturbances in columns

F3_disturbance= max(0, F3 +  w(:, 1)); 
F4_disturbance = max(0, F4 + w(:, 2));
d = [F3_disturbance, F4_disturbance];

% NOISE MEASUREMENTS
R = eye(4);
Lr = chol(R,'lower');
v = Lr*randn(4,n);
v = v';

% STEADY STATE DEFINITION
xs0 = [m10; m20; m30; m40]; % [g] Initial guess on xs
us = [F1; F2]; % [cm3/s] Flow rates
ds =  mean(d,1)'; % Average model disturbances noise
vs = mean(v,1)'; % Measurement noise

[xs2,ys2,zs2,A,B,Bv,C,D]  = StateSpaceModeling(p,xs0,us,ds,vs);
xs2
ys2

[taus2, gains2,Ts2,Ff,Gg,Gv,Cc,Dd] = DiscreteModeling(A,B,Bv,C,D,F1,F2,'Stochastic2 Model');
taus2
gains2
Ts2

markov_params2 = MarkovParam(Ff,Gg,Cc,Dd,n_markov,'Stochastic2 Model');