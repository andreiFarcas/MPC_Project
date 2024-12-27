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
%% 5.1 Continuous-time linearized models
% FOR LINEAR DETERMINISTIC MODEL 


% Steady State Definition
xs0 = [5000; 5000; 5000; 5000]; % [g] Initial guess on xs
us = [F1; F2]; % [cm3/s] Flow rates
d = [F3;F4]; % model disturbances
v=[0;0;0;0]; % measurement noise
steadyStateFunction = @(x) ModifiedFourTankSystem(0, x, us, d, p); % Define the function handle for steady-state
xs = fsolve(steadyStateFunction, xs0) % Solve using fsolve
ys = FourTankSystemSensor(xs,p,v) % [cm]
zs = FourTankSystemOutput(xs(1:2),p); % [cm]
h1s = ys(1); h2s = ys(2); h3s = ys(3); h4s = ys(4);

% Linearization
A = zeros(4,4);
A(1,1) = -(a1/A1) * sqrt(g*A1*rho / (2 * xs(1)));
A(1,3) = (a3/A3) * sqrt(g*A3*rho / (2 * xs(3)));
A(2,2) = -(a2/A2) * sqrt(g*A2*rho / (2 * xs(2)));
A(2,4) = (a4/A4) * sqrt(g*A4*rho / (2 * xs(4)));
A(3,3) = -(a3/A3) * sqrt(g*A3*rho / (2 * xs(3)));
A(4,4) = -(a4/A4) * sqrt(g*A4*rho / (2 * xs(4)));

B=[rho*gamma1 0;0 rho*gamma2; 0 rho*(1-gamma2); rho*(1-gamma1) 0];

Bv = [0 0; 0 0; rho 0; 0 rho];
C = [1/(rho*A1) 0 0 0; 0 1/(rho*A2) 0 0];

D = 0;

%% 5.2 Continuous-time transfer functions 
% FOR LINEAR DETERMINISTIC MODEL 
sys = ss(A,[B, Bv],C,D); 
G = tf(sys)

% Tracer les diagrammes de Bode
figure;
step(0.1*F1*G(1,1), 0.1*F2*G(1,2), 0.1*F1*G(2,1), 0.1*F2*G(2,2));
legend('G11', 'G12', 'G21', 'G22');
title('Step response of 0.1*F1 and 0.1*F2');
%% 5.3 Poles and gains
gains = dcgain(G) % Steady-state gain (DC gain)
poles = pole(G);
taus = -1 ./ real(poles) % Time constants

%% 5.4 Discrete-time state space models
tau_min = min(taus);
Ts = tau_min/10; %[s]
[Ff,Gg] = c2d(A,B,Ts);
[~,Gv] = c2d(A,Bv,Ts);
Cc = C;
Dd=D; 

%% 5.5 Markov parameters
n_markov = 100;
markov_params = zeros(n_markov, size(Cc, 1), size(Gg, 2));

markov_params(1, :, :) = Dd;  
for k = 2:n_markov
    markov_params(k, :, :) = Cc * (Ff^(k-1)) * Gg;  
end

figure;
for i = 1:size(Cc, 1)
    for j = 1:size(Gg, 2)
        subplot(size(Cc, 1), size(Gg, 2), (i-1)*size(Gg, 2) + j);
        stem(0:n_markov-1, squeeze(markov_params(:, i, j)));
        title(sprintf('Markov Parameter G_{%d%d}', i, j));
        xlabel('Time step (k)');
        ylabel('Markov Parameter');
    end
end