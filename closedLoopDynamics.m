% Define the dynamic control system
function [dx, F1, F2, e_prev, e_int, e_prev_time] = closedLoopDynamics(t, x, u, p, h_ref, Kp, Ki, Kd, e_prev, e_int, e_prev_time)
    % Constants
    rho = p(12);
    A = p(5:8);
    g = p(11);

    % Compute tank heights
    h = x ./ (rho * A); % Convert liquid masses to heights

    % Compute errors
    e = h_ref - h(1:2); % Error for h1 and h2

    % Handle the first step differently because we have no derivative yet
    if t == 0
        e_dot = [0; 0];
    else
        e_dot = (e - e_prev) / (t - e_prev_time); % Derivative error
    end

    % Integral error update
    e_int = e_int + e * (t - e_prev_time); % Accumulate error over time

    % PID Control Law
    F1 = Kp(1) * e(1) + Ki(1) * e_int(1) + Kd(1) * e_dot(1);
    F2 = Kp(2) * e(2) + Ki(2) * e_int(2) + Kd(2) * e_dot(2);

    % Ensure control flows will never be -
    F1 = max(0, F1);
    F2 = max(0, F2);

    u = [F1; F2; u(3); u(4)];

    % Compute system states dx
    dx = ModifiedFourTankSystem(t, x, u, [], p);
    
    % We try to ensure always dx is a column vector
    dx = dx(:);

    % Update previous errors and time with current values 
    e_prev = e;
    e_prev_time = t;
end
