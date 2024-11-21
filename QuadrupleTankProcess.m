function x_dot = QuadrupleTankProcess(t, x, u, p, tf)
%%%
%
% Takes as input the masses in the 4 tanks and the flow rates F1 and F2
% Basically returns the differential equations that describe the process
% 
%%%

% Extract variables from the state vector
m1 = x(1); % mass in tank 1
m2 = x(2); % mass in tank 2
m3 = x(3); % mass in tank 3
m4 = x(4); % mass in tank 4

% Extract flows from input (u)
% F1 = u(1, t); F2 = u(2, t); - Doesn't work because Matlab doesn't allow
% "continuous" time to index a matrix

% Extract the time vector
t_data = linspace(0, tf, length(u)); % p(13) is the final time

% Interpolate flows from input (u) at the current time `t`
F1 = interp1(t_data, u(1, :), t, 'linear'); 
F2 = interp1(t_data, u(2, :), t, 'linear'); 

% Extract the rest of parameters
a1 = p(1);     %[cm2] Area of outlet pipe 1
a2 = p(2);     %[cm2] Area of outlet pipe 2
a3 = p(3);     %[cm2] Area of outlet pipe 3
a4 = p(4);     %[cm2] Area of outlet pipe 4

A1 = p(5);   %[cm2] Cross sectional area of tank 1
A2 = p(6);   %[cm2] Cross sectional area of tank 2
A3 = p(7);   %[cm2] Cross sectional area of tank 3
A4 = p(8);   %[cm2] Cross sectional area of tank 4

g = p(11);         %[cm/s2] The acceleration of gravity

gamma1 = p(9);  % Flow distribution constant. Valve 1
gamma2 = p(10);  % Flow distribution constant. Valve 2

rho = p(12); % Density [g/cm3]

% Compute the liquid heights from masses
h1 = m1 / (rho * A1);
h2 = m2 / (rho * A2);
h3 = m3 / (rho * A3);
h4 = m4 / (rho * A4);

% Computing the input flows for each tank
q_in1 = gamma1 * F1;
q_in2 = gamma2 * F2;
q_in3 = (1 - gamma2) * F2;
q_in4 = (1 - gamma1) * F1;

% Computing the output flows of each tank
q1 = a1 * sqrt(2 * g * h1);
q2 = a2 * sqrt(2 * g * h2);
q3 = a3 * sqrt(2 * g * h3);
q4 = a4 * sqrt(2 * g * h4);

% Defining the differential equations of the model
x_dot = zeros(4, 1);
x_dot(1) = rho * q_in1 + rho * q3 - rho * q1;
x_dot(2) = rho * q_in2 + rho * q4 - rho * q2;
x_dot(3) = rho * q_in3 - rho * q3;
x_dot(4) = rho * q_in4 - rho * q4;

end